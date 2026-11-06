#!/usr/bin/env Rscript
# DNA coral virome heatmap: one TSV source-data matrix, one PDF.
# Run: Rscript run_heatmap.R inputs results
# In RStudio: setwd() to this package directory, then source this file.
# Inputs are read-only. Calculation, category ordering, labels, colors and PDF
# geometry follow the supplied Project-aware analysis.
# Source data contain the exact plotted log1p matrix and full group labels.
# Missing values are written as text NA; numeric zeros are retained.
# Source data are not rounded. The heatmap uses ln(1 + arithmetic mean TPM).

required_packages <- c("tidyverse", "pheatmap", "RColorBrewer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "),
       ". Use install.packages(c('tidyverse', 'pheatmap', 'RColorBrewer')).")
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
}

validate_key <- function(x, label) {
  if (anyNA(x) || any(trimws(as.character(x)) == "") || anyDuplicated(x)) {
    stop(label, " must contain unique, non-missing, non-empty identifiers.")
  }
}

# Export the exact final heatmap matrix with full group IDs as sample columns.
write_heatmap_tsv <- function(x, row_info, column_labels, path) {
  if (!is.matrix(x) || nrow(x) != nrow(row_info) || ncol(x) != length(column_labels)) {
    stop("Matrix, row metadata and column labels do not align.")
  }
  if (anyDuplicated(column_labels) || anyNA(column_labels) || any(!nzchar(column_labels))) {
    stop("Heatmap sample labels must be unique and non-empty.")
  }
  text_values <- c(column_labels, row_info$Row_ID, row_info$Phylum, row_info$Label)
  # R converts these string escapes to TSV control characters before matching.
  if (any(grepl("[\t\r\n]", text_values))) {
    stop("Labels contain a tab or newline and cannot be represented safely in TSV.")
  }
  values <- as.data.frame(x, check.names = FALSE)
  names(values) <- column_labels
  out <- data.frame(
    Row_Order = seq_len(nrow(x)),
    Row_ID = row_info$Row_ID,
    Phylum = row_info$Phylum,
    Label = row_info$Label,
    values,
    check.names = FALSE
  )
  for (j in seq.int(5L, ncol(out))) {
    out[[j]] <- sprintf("%.17g", x[, j - 4L])
  }
  write.table(out, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA",
              fileEncoding = "UTF-8")
  invisible(path)
}

main <- function(input_dir = "inputs", output_dir = "results") {
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
input_files <- c("dna_tpm.tsv", "dna_virus_taxonomy.tsv",
                 "sample_metadata.tsv")
input_paths <- file.path(input_dir, input_files)
if (!all(file.exists(input_paths))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_paths)], collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
if (identical(input_dir, output_dir)) stop("Use a separate output directory.")
fig_dir <- output_dir

# ==============================================================================
# DNA coral virome heatmap calculation
# Project-aware DNA heatmap with hierarchical annotations and scale 0-18.
# ==============================================================================



# ==============================================================================
# 1. Read data
# ==============================================================================
message("[1/7] Reading and validating inputs")
tpm_df <- read_tsv(file.path(input_dir, "dna_tpm.tsv"), show_col_types = FALSE)
tax_df <- read_tsv(file.path(input_dir, "dna_virus_taxonomy.tsv"), show_col_types = FALSE)
meta_df <- read_tsv(file.path(input_dir, "sample_metadata.tsv"), show_col_types = FALSE)

if(!"Sample" %in% colnames(meta_df) && "SampleID" %in% colnames(meta_df)) {
  meta_df <- meta_df %>% rename(Sample = SampleID)
}


# Added validation: prevent accidental join multiplication or invalid TPM input.
# No values are imputed, rounded, normalised or silently discarded here.
require_columns(tpm_df, "Genome", "TPM table")
require_columns(tax_df, c("Genome", "Phylum", "Family"), "Taxonomy table")
require_columns(meta_df, c("Sample", "Project", "Coral_Host_Species", "Ocean_Basin",
                           "Sample_Type", "Status"), "Metadata table")
for (x in list(tpm_df, tax_df, meta_df)) {
  if (nrow(readr::problems(x))) stop("Input parsing problems: inspect the source TSV files.")
}
validate_key(tpm_df$Genome, "TPM Genome")
validate_key(tax_df$Genome, "Taxonomy Genome")
validate_key(str_trim(as.character(meta_df$Sample)), "Metadata Sample after trimming")
if (anyNA(meta_df$Project) || any(str_trim(as.character(meta_df$Project)) == "")) {
  stop("Project must be provided for every metadata sample.")
}
sample_columns <- setdiff(names(tpm_df), "Genome")
if (!length(sample_columns) || !nrow(tpm_df)) stop("The TPM table is empty.")
if (!all(vapply(tpm_df[sample_columns], is.numeric, logical(1)))) {
  stop("Every TPM sample column must be numeric.")
}
if (any(!vapply(tpm_df[sample_columns],
                function(x) all(is.finite(x) & x >= 0), logical(1)))) {
  stop("TPM must be finite and non-negative; missing values are not assumed to be zero.")
}
if (length(setdiff(tpm_df$Genome, tax_df$Genome))) {
  stop("Some TPM Genome IDs have no taxonomy record; resolve them before analysis.")
}
# The || delimiter is part of the original compound identifiers.
id_fields <- c(as.character(tax_df$Phylum), as.character(tax_df$Family),
               as.character(meta_df$Project), as.character(meta_df$Coral_Host_Species),
               as.character(meta_df$Ocean_Basin),
               as.character(meta_df$Sample_Type), as.character(meta_df$Status))
if (any(grepl("||", id_fields, fixed = TRUE), na.rm = TRUE)) {
  stop("An input label contains the reserved compound-ID delimiter ||.")
}
# ==============================================================================
# 2. Clean taxonomy
# ==============================================================================
message("[2/7] Assigning DNA taxonomy categories")
is_invalid <- function(x) { is.na(x) | x == "" | tolower(x) %in% c("unclassified", "unknown", "nan", "viruses") }

tax_basic <- tax_df %>%
  rowwise() %>%
  mutate(
    Row_Label_Raw = case_when(
      !is_invalid(Family) ~ Family,
      !is_invalid(Phylum) ~ paste0("unclass_", Phylum),
      TRUE ~ "Unclassified DNA Virus"
    ),
    Phylum_Block = ifelse(is_invalid(Phylum), "Unclassified", Phylum)
  ) %>% ungroup()

# Rank Uroviricota categories by total input TPM
genome_sums <- rowSums(tpm_df %>% select(-Genome))
uro_stats <- data.frame(Genome = tpm_df$Genome, Total_TPM = genome_sums) %>%
  left_join(tax_basic, by = "Genome") %>%
  filter(Phylum_Block == "Uroviricota") %>% 
  group_by(Row_Label_Raw) %>% 
  summarise(Sum_TPM = sum(Total_TPM)) %>% 
  arrange(desc(Sum_TPM))

top20_families <- uro_stats$Row_Label_Raw[1:min(20, nrow(uro_stats))]

tax_clean <- tax_basic %>%
  mutate(
    Row_Label = case_when(
      Phylum_Block == "Uroviricota" & Row_Label_Raw %in% top20_families ~ Row_Label_Raw,
      Phylum_Block == "Uroviricota" & !(Row_Label_Raw %in% top20_families) ~ "Others_Uroviricota",
      TRUE ~ Row_Label_Raw
    )
  ) %>% select(Genome, Phylum_Block, Row_Label)

# ==============================================================================
# 3. Metadata and hierarchical groups
# ==============================================================================
message("[3/7] Cleaning metadata and applying Project groups")

tissue_levels <- c("Coral tissue", "Coral mucus", "Coral skeleton", "Coral fragment", "Other") 
detail_levels <- c(
  "Healthy/Normal", "Bleached", 
  "Disease (Bleaching lesion)", "Disease(SCTLD)", "Disease(Dark Spot Syndrome)", "Disease (Pathogen-induced disease)", 
  "Stressed(Heat stress)", "Stressed(Acidification stress)", "Stressed(Combined)", 
  "Other", "Unkonwn" 
)
broad_levels <- c("Healthy", "Bleached", "Disease", "Stressed", "Other", "Unknown")

# Health-status recoding
get_detail_status <- function(s) {
  s_norm <- str_trim(s)
  if (grepl("SCTLD", s_norm, ignore.case=T)) return("Disease(SCTLD)")
  if (grepl("Dark Spot", s_norm, ignore.case=T)) return("Disease(Dark Spot Syndrome)")
  if (grepl("Bleaching lesion", s_norm, ignore.case=T)) return("Disease (Bleaching lesion)")
  if (grepl("Pathogen", s_norm, ignore.case=T)) return("Disease (Pathogen-induced disease)")
  if (grepl("Bleached", s_norm, ignore.case=T)) return("Bleached")
  has_heat <- grepl("Heat", s_norm, ignore.case=T); has_acid <- grepl("Acidification", s_norm, ignore.case=T)
  if (has_heat & has_acid) return("Stressed(Combined)")
  if (has_heat) return("Stressed(Heat stress)")
  if (has_acid) return("Stressed(Acidification stress)")
  if (grepl("Healthy|Normal", s_norm, ignore.case=T)) return("Healthy/Normal")
  return("Other")
}

get_broad_status <- function(detail) {
  if (detail == "Healthy/Normal") return("Healthy")
  if (detail == "Bleached") return("Bleached")
  if (grepl("Disease", detail)) return("Disease")
  if (grepl("Stressed", detail)) return("Stressed")
  return("Other")
}

meta_pre <- meta_df %>%
  mutate(
    Sample = str_trim(as.character(Sample)),
    Project = str_trim(as.character(Project)),
    Species = ifelse(is.na(Coral_Host_Species)|Coral_Host_Species=="", "Unknown Species", str_trim(Coral_Host_Species)),
    Ocean_Basin = ifelse(is.na(Ocean_Basin)|Ocean_Basin=="", "Unknown Ocean", str_trim(Ocean_Basin)),
    Sample_Type_Raw = str_to_sentence(ifelse(is.na(Sample_Type), "Other", str_trim(Sample_Type))),
    Sample_Type = factor(ifelse(Sample_Type_Raw %in% tissue_levels, Sample_Type_Raw, "Other"), levels = tissue_levels),
    Detail_Status = sapply(Status, get_detail_status)
  )

# Rank species by geographical and health-status combinations
species_ranking <- meta_pre %>%
  group_by(Species) %>%
  summarise(
    Unique_Combos = n_distinct(paste(Ocean_Basin, Sample_Type, Detail_Status)),
    Has_Disease = any(grepl("Disease", Detail_Status)),
    .groups = "drop"
  ) %>%
  arrange(desc(Unique_Combos), desc(Has_Disease), Species)

ordered_species_names <- species_ranking$Species

meta_clean <- meta_pre %>%
  mutate(
    Species = factor(Species, levels = ordered_species_names),
    Broad_Status = factor(sapply(Detail_Status, get_broad_status), levels = broad_levels),
    Detail_Status = factor(Detail_Status, levels = detail_levels)
  ) %>%
  arrange(Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, Project)

# ==============================================================================
# 4. Aggregate data (sum, then mean)
# ==============================================================================
message("[4/7] Summing TPM within samples, then averaging within Project-aware groups")

sample_family_totals <- tpm_df %>%
  pivot_longer(-Genome, names_to = "Sample", values_to = "TPM") %>%
  filter(Sample %in% meta_clean$Sample) %>%
  left_join(tax_clean, by = "Genome") %>%
  inner_join(meta_clean, by = "Sample") %>%
  group_by(Sample, Project, Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, Phylum_Block, Row_Label) %>%
  summarise(Sample_Family_Total = sum(TPM, na.rm = TRUE), .groups = "drop")

plot_data <- sample_family_totals %>%
  group_by(Project, Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, Phylum_Block, Row_Label) %>%
  summarise(Mean_TPM = mean(Sample_Family_Total, na.rm = TRUE), .groups = "drop")

heatmap_data <- plot_data %>%
  mutate(Col_ID = paste(Project, Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, sep = "||"),
         Row_ID = paste(Phylum_Block, Row_Label, sep = "||"))

heatmap_matrix_raw <- heatmap_data %>%
  select(Row_ID, Col_ID, Mean_TPM) %>%
  pivot_wider(names_from = Col_ID, values_from = Mean_TPM, values_fill = 0) %>%
  column_to_rownames("Row_ID")

# ==============================================================================
# 5. Heatmap settings and annotation order
# ==============================================================================
message("[5/7] Configuring annotations and the shared heatmap scale")

col_info <- data.frame(Col_ID = colnames(heatmap_matrix_raw)) %>%
  separate(Col_ID, into = c("Project", "Species", "Location", "Tissue", "Broad", "Detail"), sep = "\\|\\|", remove = FALSE) %>%
  mutate(
    Species = factor(Species, levels = ordered_species_names),
    Tissue = factor(Tissue, levels = tissue_levels),
    Broad = factor(Broad, levels = broad_levels),
    Detail = factor(Detail, levels = detail_levels)
  ) %>%
  arrange(Species, Location, Tissue, Broad, Detail, Project) %>%
  mutate(Project_Species = paste(Project, Species, sep = "_"))

row_info <- data.frame(Row_ID = rownames(heatmap_matrix_raw)) %>%
  separate(Row_ID, into = c("Phylum", "Label"), sep = "\\|\\|", remove = FALSE) %>%
  arrange(Phylum, Label)

mat <- heatmap_matrix_raw[row_info$Row_ID, col_info$Col_ID, drop=FALSE]
rownames(mat) <- row_info$Label
mat_log <- log1p(as.matrix(mat))

# pheatmap annotation display order
# The final annotation column is displayed at the top
ann_col <- data.frame(
  "Detail Status" = col_info$Detail,
  "Status"        = col_info$Broad,
  "Tissue"        = col_info$Tissue,
  "Ocean Basin"   = col_info$Location,
  check.names = FALSE
)
rownames(ann_col) <- col_info$Col_ID

# Annotation colors
unique_oceans <- unique(col_info$Location)
ocean_colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(unique_oceans)), unique_oceans)

detail_colors <- c(
  "Healthy/Normal"="#7EC87E", "Bleached"="#FCC085", 
  "Disease (Bleaching lesion)"="#FC8D3B", "Disease(SCTLD)"="#800025", 
  "Disease(Dark Spot Syndrome)"="#BC0025", "Disease (Pathogen-induced disease)"="#E26E7C", 
  "Stressed(Heat stress)"="#756BB0", "Stressed(Acidification stress)"="#9E58F5", "Stressed(Combined)"="black",
  "Other"="#D9D9D9"
)

ann_colors <- list(
  "Ocean Basin"   = ocean_colors,
  "Tissue"        = c("Coral tissue"="#D95F02", "Coral mucus"="#1B9E77", "Coral skeleton"="#8D40CA", "Coral fragment"="#E72989", "Other"="#999999"),
  "Status"        = c("Healthy"="#7FC87E", "Bleached"="#FCC085", "Disease"="#E21A1C", "Stressed"="#3C15D3", "Other"="#999999"),
  "Detail Status" = detail_colors,
  "Phylum"        = setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(unique(row_info$Phylum))), unique(row_info$Phylum))
)



  # ============================================================================
  # 6. Export the exact plotted matrix as one plain TSV
message("[6/7] Writing Heatmap_Log1p source TSV")
heat_colors <- colorRampPalette(c("white", "orange", "red", "black"))(100)
heat_breaks <- seq(0, 18, length.out = 100)
stopifnot(
  identical(rownames(ann_col), colnames(mat_log)),
  identical(rownames(mat_log), as.character(row_info$Label)),
  !anyDuplicated(row_info$Row_ID),
  !anyDuplicated(col_info$Col_ID),
  all(is.finite(mat_log)),
  all(mat_log >= 0),
  identical(unname(mat_log), unname(log1p(as.matrix(mat))))
)
source_labels <- as.character(col_info$Col_ID)
source_data_path <- file.path(output_dir, "DNA_Heatmap_Log1p_SourceData.tsv")
write_heatmap_tsv(mat_log, row_info, source_labels, source_data_path)

# ==============================================================================
# 7. Draw PDF
# ==============================================================================
message("[7/7] Drawing PDF")
dynamic_width <- max(16, ncol(mat) * 0.45) 
dynamic_height <- max(12, nrow(mat) * 0.35)

pdf(file.path(fig_dir, "DNA_Virome_Global_Ocean_Hierarchy.pdf"),
    width = dynamic_width, height = dynamic_height)

pheatmap(mat_log,
         color = colorRampPalette(c("white", "orange", "red", "black"))(100),
         breaks = seq(0, 18, length.out = 100),
         border_color = NA,
         cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_col, 
         annotation_row = data.frame(Phylum=row_info$Phylum, row.names=rownames(mat)),
         annotation_colors = ann_colors,
         gaps_col = which(col_info$Species[-1] != col_info$Species[-nrow(col_info)]),
         gaps_row = which(row_info$Phylum[-1] != row_info$Phylum[-nrow(row_info)]),
         labels_col = col_info$Project_Species,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize = 10, fontsize_row = 9, 
         angle_col = 45,
         cellwidth = 12, cellheight = 10,
         main = "Global DNA Virome: Project_Species | Basin > Tissue > Status (Scale 0-18)"
)

dev.off()
message("PDF: DNA_Virome_Global_Ocean_Hierarchy.pdf")

  message("Heatmap_Log1p source TSV: ", source_data_path)
  invisible(list(mean_TPM = as.matrix(mat), log1p_mean_TPM = mat_log,
    row_info = row_info, col_info = col_info, labels_col = col_info$Project_Species,
    column_annotation = ann_col, annotation_colors = ann_colors,
    top20_categories = top20_families, species_ranking = species_ranking,
    heatmap_colors = heat_colors, heatmap_breaks = heat_breaks,
    pdf_width_inches = dynamic_width, pdf_height_inches = dynamic_height,
    source_data_path = source_data_path))
}

# Run from this package directory, or pass input/output directories explicitly.
args <- commandArgs(trailingOnly = TRUE)
if (interactive()) args <- character()
if (length(args) > 2L) stop("Usage: Rscript run_heatmap.R [input_dir] [output_dir]")
analysis <- main(if (length(args) >= 1L) args[[1L]] else "inputs",
                 if (length(args) >= 2L) args[[2L]] else "results")
