#!/usr/bin/env Rscript
# RNA coral virome heatmap: one TSV source-data matrix, one PDF.
# Run: Rscript run_heatmap.R inputs results
# In RStudio: setwd() to this package directory, then source this file.
# Inputs are read-only. Calculation, category ordering, labels, colors and PDF
# geometry follow the supplied Project-aware analysis.
# Source data contain the exact plotted log1p matrix and full group labels.
# Missing values are written as text NA; numeric zeros are retained.
# Source data are not rounded. The heatmap uses ln(1 + arithmetic mean TPM).

required_packages <- c("tidyverse", "pheatmap", "RColorBrewer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "),
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
make_col_id <- function(x) {
  with(x, paste(Species, Ocean_Basin, Sample_Type, Broad_Status,
                Detail_Status, Group_Project, sep = "||"))
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
  input_files <- c("rna_virus_tpm.tsv", "rna_virus_taxonomy.tsv",
                   "sample_metadata.tsv")
  input_paths <- file.path(input_dir, input_files)
  if (!all(file.exists(input_paths))) {
    stop("Missing inputs: ", paste(input_files[!file.exists(input_paths)], collapse = ", "))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (identical(input_dir, output_dir)) stop("Use a separate output directory.")
  fig_dir <- output_dir

  message("[1/7] Reading and validating inputs")
  tpm_df <- read_tsv(input_paths[1], show_col_types = FALSE)
  tax_df <- read_tsv(input_paths[2], show_col_types = FALSE)
  meta_df <- read_tsv(input_paths[3], show_col_types = FALSE)
  require_columns(tpm_df, "Genome", "TPM table")
  require_columns(tax_df, c("Genome", "Phylum", "Order", "Family"), "Taxonomy")
  require_columns(meta_df, c("Sample", "Project", "Coral_Host_Species", "Ocean_Basin",
                             "Sample_Type", "Status"), "Metadata")
  for (x in list(tpm_df, tax_df, meta_df)) {
    if (nrow(readr::problems(x))) stop("Parsing problems in an input table.")
  }
  validate_key(tpm_df$Genome, "TPM Genome")
  validate_key(tax_df$Genome, "Taxonomy Genome")
  validate_key(str_trim(as.character(meta_df$Sample)), "Metadata Sample after trimming")
  sample_columns <- setdiff(names(tpm_df), "Genome")
  if (!length(sample_columns) || !nrow(tpm_df)) stop("Empty TPM input.")
  if (!all(vapply(tpm_df[sample_columns], is.numeric, logical(1)))) {
    stop("TPM sample columns must be numeric.")
  }
  if (any(!vapply(tpm_df[sample_columns], function(x) all(is.finite(x) & x >= 0), logical(1)))) {
    stop("TPM must be finite and non-negative; missing values are not assumed to be zero.")
  }
  if (length(setdiff(tpm_df$Genome, tax_df$Genome))) stop("A TPM Genome lacks taxonomy.")
  id_fields <- c(unlist(tax_df[vapply(tax_df, is.character, logical(1))]),
                 unlist(meta_df[vapply(meta_df, is.character, logical(1))]))
  if (any(grepl("||", id_fields, fixed = TRUE), na.rm = TRUE)) {
    stop("An input label contains the reserved compound-ID delimiter ||.")
  }
  # RNA-specific taxonomy hierarchy; no Top20 filtering or category pooling.
  message("[2/7] Assigning RNA taxonomy categories")
  is_invalid <- function(x) {
    is.na(x) | x == "" | tolower(x) %in% c("unclassified", "unknown", "nan", "viruses")
  }
  tax_clean <- tax_df %>% rowwise() %>% mutate(
    Row_Label = case_when(
      !is_invalid(Family) ~ Family,
      !is_invalid(Order) ~ paste0("unclass_", Order),
      !is_invalid(Phylum) ~ paste0("unclass_", Phylum),
      TRUE ~ "Unclassified RNA Virus"),
    Phylum_Block = ifelse(is_invalid(Phylum), "Unclassified", Phylum)
  ) %>% ungroup() %>% select(Genome, Phylum_Block, Row_Label)

  message("[3/7] Cleaning metadata and applying explicit Project rules")
  tissue_levels <- c("Coral tissue", "Coral mucus", "Coral skeleton", "Coral fragment", "Other")
  detail_levels <- c(
    "Healthy/Normal", "Bleached", "Diseased(black band disease)", "Diseased(SCTLD)",
    "Diseased (white plague disease)", "Stressed(Heat stress)",
    "Stressed(Acidification stress)", "Stressed(Cold stress)", "Stressed(Combined)",
    "Other", "Unknown")
  broad_levels <- c("Healthy", "Bleached", "Diseased", "Stressed", "Other", "Unknown")
  # Recognition order matches the supplied RNA script, including Mix -> Other.
  get_detail_status <- function(s) {
    if (is.na(s) || s == "") return("Other")
    s_norm <- str_trim(as.character(s))
    if (grepl("Mix", s_norm, ignore.case = TRUE)) return("Other")
    if (grepl("Healthy|Normal", s_norm, ignore.case = TRUE)) return("Healthy/Normal")
    if (grepl("Bleached", s_norm, ignore.case = TRUE)) return("Bleached")
    if (grepl("black band", s_norm, ignore.case = TRUE)) return("Diseased(black band disease)")
    if (grepl("SCTLD", s_norm, ignore.case = TRUE)) return("Diseased(SCTLD)")
    if (grepl("white plague", s_norm, ignore.case = TRUE)) return("Diseased (white plague disease)")
    has_heat <- grepl("Heat", s_norm, ignore.case = TRUE)
    has_acid <- grepl("Acidification", s_norm, ignore.case = TRUE)
    has_cold <- grepl("Cold", s_norm, ignore.case = TRUE)
    if (has_heat & has_acid) return("Stressed(Combined)")
    if (has_heat) return("Stressed(Heat stress)")
    if (has_acid) return("Stressed(Acidification stress)")
    if (has_cold) return("Stressed(Cold stress)")
    if (grepl("Disease|Diseased", s_norm, ignore.case = TRUE)) return("Other")
    "Other"
  }
  get_broad_status <- function(detail) {
    if (detail == "Healthy/Normal") return("Healthy")
    if (detail == "Bleached") return("Bleached")
    if (grepl("Diseased", detail)) return("Diseased")
    if (grepl("Stressed", detail)) return("Stressed")
    "Other"
  }
  # Exact supplied accession list, NOT a broad PRJNA405* prefix match.
  avila_projects <- c(
    "PRJNA405605", "PRJNA405606", "PRJNA405607", "PRJNA405608", "PRJNA405609", "PRJNA405610",
    "PRJNA405611", "PRJNA405612", "PRJNA405613", "PRJNA405614", "PRJNA405615", "PRJNA405616",
    "PRJNA405617", "PRJNA405618", "PRJNA405952", "PRJNA405953", "PRJNA405954", "PRJNA405955",
    "PRJNA405956", "PRJNA405957", "PRJNA405958", "PRJNA405995", "PRJNA405996", "PRJNA405997")
  meta_pre <- meta_df %>% mutate(
    Sample = str_trim(as.character(Sample)),
    Project = ifelse(is.na(Project) | Project == "", "Unknown Project", str_trim(as.character(Project))),
    Analysis_Project = ifelse(Project %in% avila_projects, "Avila2021", Project),
    Species = ifelse(is.na(Coral_Host_Species) | Coral_Host_Species == "",
                     "Unknown Species", str_trim(as.character(Coral_Host_Species))),
    Ocean_Basin = ifelse(is.na(Ocean_Basin) | Ocean_Basin == "",
                         "Unknown Ocean", str_trim(as.character(Ocean_Basin))),
    Sample_Type = factor(ifelse(Sample_Type %in% tissue_levels, Sample_Type, "Other"),
                          levels = tissue_levels),
    Detail_Status_Raw = sapply(Status, get_detail_status)
  ) %>% mutate(
    Detail_Status = factor(Detail_Status_Raw, levels = detail_levels),
    Broad_Status = factor(sapply(Detail_Status_Raw, get_broad_status), levels = broad_levels),
    Group_Project = Analysis_Project)
  # Species ranking remains independent of Project and follows the RNA script.
  species_ranking <- meta_pre %>% group_by(Species) %>% summarise(
    Unique_Combos = n_distinct(paste(Ocean_Basin, Sample_Type, Detail_Status)),
    .groups = "drop") %>% arrange(desc(Unique_Combos), Species)
  meta_clean <- meta_pre %>% mutate(
    Species = factor(Species, levels = species_ranking$Species)
  ) %>% arrange(Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, Group_Project)
  if (anyNA(meta_clean[c("Species", "Sample_Type", "Broad_Status", "Detail_Status")])) {
    stop("An unexpected metadata level produced NA after recoding.")
  }
  used_samples <- intersect(sample_columns, meta_clean$Sample)
  if (!length(used_samples)) stop("No TPM samples match metadata.")
  meta_used <- meta_clean %>% filter(Sample %in% used_samples)
  meta_used$Col_ID <- make_col_id(meta_used)
  group_counts <- meta_used %>% group_by(Col_ID) %>% summarise(
    n_samples = n(), n_original_projects = n_distinct(Project),
    Original_Projects = paste(sort(unique(Project)), collapse = ";"), .groups = "drop")

  message("[4/7] Summing TPM within samples, then averaging within Project-aware groups")
  sample_family_totals <- tpm_df %>%
    pivot_longer(-Genome, names_to = "Sample", values_to = "TPM") %>%
    filter(Sample %in% meta_clean$Sample) %>%
    left_join(tax_clean, by = "Genome") %>% inner_join(meta_clean, by = "Sample") %>%
    group_by(Sample, Group_Project, Species, Ocean_Basin, Sample_Type,
             Broad_Status, Detail_Status, Phylum_Block, Row_Label) %>%
    summarise(Sample_Family_Total = sum(TPM, na.rm = TRUE), .groups = "drop")
  plot_data <- sample_family_totals %>%
    group_by(Group_Project, Species, Ocean_Basin, Sample_Type, Broad_Status,
             Detail_Status, Phylum_Block, Row_Label) %>%
    summarise(Mean_TPM = mean(Sample_Family_Total, na.rm = TRUE), .groups = "drop")
  heatmap_data <- plot_data
  heatmap_data$Col_ID <- make_col_id(heatmap_data)
  heatmap_data$Row_ID <- with(heatmap_data, paste(Phylum_Block, Row_Label, sep = "||"))
  heatmap_matrix_raw <- heatmap_data %>% select(Row_ID, Col_ID, Mean_TPM) %>%
    pivot_wider(names_from = Col_ID, values_from = Mean_TPM, values_fill = 0) %>%
    column_to_rownames("Row_ID")
  col_info <- data.frame(Col_ID = colnames(heatmap_matrix_raw)) %>%
    separate(Col_ID, into = c("Species", "Location", "Tissue", "Broad", "Detail", "Group_Project"),
             sep = "\\|\\|", remove = FALSE) %>%
    mutate(Species = factor(Species, levels = species_ranking$Species)) %>%
    # Preserve supplied RNA final-column ordering: these four fields are character.
    arrange(Species, Location, Tissue, Broad, Detail, Group_Project) %>%
    mutate(Analysis_Project = Group_Project,
           Display_Label = paste(Group_Project, as.character(Species), sep = "_"))
  row_info <- data.frame(Row_ID = rownames(heatmap_matrix_raw)) %>%
    separate(Row_ID, into = c("Phylum", "Label"), sep = "\\|\\|", remove = FALSE) %>%
    arrange(Phylum, Label)
  mat <- heatmap_matrix_raw[row_info$Row_ID, col_info$Col_ID, drop = FALSE]
  rownames(mat) <- row_info$Label
  mat_log <- log1p(as.matrix(mat))  # natural logarithm; no percent or z-score transform

  message("[5/7] Configuring original RNA annotations and the shared DNA heatmap scale")
  ann_col <- data.frame(
    "Detail Status" = col_info$Detail, "Status" = col_info$Broad,
    "Tissue" = col_info$Tissue, "Ocean Basin" = col_info$Location, check.names = FALSE)
  rownames(ann_col) <- col_info$Col_ID
  ocean_colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(
    length(unique(col_info$Location))), unique(col_info$Location))
  detail_colors <- c(
    "Healthy/Normal" = "#A6DBA0", "Bleached" = "#FDDBC7",
    "Diseased(black band disease)" = "#800026", "Diseased(SCTLD)" = "#BD0026",
    "Diseased (white plague disease)" = "#6A3D9A", "Stressed(Heat stress)" = "#756BB1",
    "Stressed(Acidification stress)" = "#9E9AC8", "Stressed(Cold stress)" = "#9ECAE1",
    "Stressed(Combined)" = "black", "Other" = "#D9D9D9", "Unknown" = "#252525")
  ann_colors <- list(
    "Ocean Basin" = ocean_colors,
    "Tissue" = c("Coral tissue" = "#D95F02", "Coral mucus" = "#1B9E77",
                 "Coral skeleton" = "#8D40CA", "Coral fragment" = "#E72989", "Other" = "#999999"),
    "Status" = c("Healthy" = "#7FC87E", "Bleached" = "#FCC085", "Diseased" = "#E21A1C",
                 "Stressed" = "#3C15D3", "Other" = "#999999", "Unknown" = "#252525"),
    "Detail Status" = detail_colors,
    "Phylum" = setNames(colorRampPalette(brewer.pal(8, "Set2"))(
      length(unique(row_info$Phylum))), unique(row_info$Phylum)))
  # Exact DNA vectors (100 colors and 100 break values). Kept separate in exports.
  # This intentionally supersedes the 101 break values in the supplied RNA script.
  heat_colors <- colorRampPalette(c("white", "orange", "red", "black"))(100)
  heat_breaks <- seq(0, 18, length.out = 100)
  gaps_col <- which(col_info$Species[-1] != col_info$Species[-nrow(col_info)])
  gaps_row <- which(row_info$Phylum[-1] != row_info$Phylum[-nrow(row_info)])


  # ============================================================================
  # 6. Export the exact plotted matrix as one plain TSV
  message("[6/7] Writing Heatmap_Log1p source TSV")
  stopifnot(
    identical(rownames(ann_col), colnames(mat_log)),
    identical(rownames(mat_log), as.character(row_info$Label)),
    !anyDuplicated(row_info$Row_ID),
    !anyDuplicated(col_info$Col_ID),
    all(is.finite(mat_log)),
    all(mat_log >= 0),
    identical(unname(mat_log), unname(log1p(as.matrix(mat))))
  )
  source_labels <- paste(
    col_info$Analysis_Project, as.character(col_info$Species),
    col_info$Location, col_info$Tissue, col_info$Broad, col_info$Detail,
    sep = "||"
  )
  source_data_path <- file.path(output_dir, "RNA_Heatmap_Log1p_SourceData.tsv")
  write_heatmap_tsv(mat_log, row_info, source_labels, source_data_path)

  # Large page geometry and RNA font settings are deliberately retained.
  message("[7/7] Drawing PDF")
  dynamic_width <- max(16, ncol(mat) * 0.48)
  dynamic_height <- max(12, nrow(mat) * 0.35)
  plot_title <- "Global RNA Virome: Project_Species | Basin > Tissue > Status (Scale 0-18)"
  ht <- pheatmap(mat_log, color = heat_colors, breaks = heat_breaks, border_color = NA,
    cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = ann_col,
    annotation_row = data.frame(Phylum = row_info$Phylum, row.names = rownames(mat)),
    annotation_colors = ann_colors, gaps_col = gaps_col, gaps_row = gaps_row,
    labels_col = col_info$Display_Label, show_rownames = TRUE, show_colnames = TRUE,
    fontsize = 10, fontsize_col = 8, angle_col = 45, cellwidth = 12, cellheight = 10,
    main = plot_title, silent = TRUE)
  pdf(file.path(fig_dir, "RNA_Virome_Global_ProjectAware.pdf"),
      width = dynamic_width, height = dynamic_height, family = "Helvetica")
  tryCatch({grid.newpage(); grid.draw(ht$gtable)}, finally = dev.off())
  message("Heatmap_Log1p source TSV: ", source_data_path)
  invisible(list(mean_TPM = as.matrix(mat), log1p_mean_TPM = mat_log,
    row_info = row_info, col_info = col_info, labels_col = col_info$Display_Label,
    column_annotation = ann_col, annotation_colors = ann_colors,
    species_ranking = species_ranking, avila_projects = avila_projects,
    heatmap_colors = heat_colors, heatmap_breaks = heat_breaks,
    gaps_col = gaps_col, gaps_row = gaps_row,
    pdf_width_inches = dynamic_width, pdf_height_inches = dynamic_height,
    plot_title = plot_title, plot_object = ht, source_data_path = source_data_path))
}

# Run from this package directory, or pass input/output directories explicitly.
args <- commandArgs(trailingOnly = TRUE)
if (interactive()) args <- character()
if (length(args) > 2L) stop("Usage: Rscript run_heatmap.R [input_dir] [output_dir]")
analysis <- main(if (length(args) >= 1L) args[[1L]] else "inputs",
                 if (length(args) >= 2L) args[[2L]] else "results")
