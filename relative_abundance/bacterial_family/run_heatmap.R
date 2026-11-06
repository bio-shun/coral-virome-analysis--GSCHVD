#!/usr/bin/env Rscript
# Coral MAG family heatmap and a single plain TSV of its log1p matrix.
# Run from this directory:
#   Rscript run_heatmap.R inputs results
# Defaults: inputs/ and results/. Inputs are never modified.
# Group ID: Project||Species||Ocean_Basin||Sample_Type||Broad_Status||Detail_Status
# Retained: 51 specified families + Others, original recoding, mean over all
# MAG x sample TPM cells within a family/group, >=1 maximum mean TPM filter,
# log1p, original row/column hierarchy, annotations, palette and fixed 0-15 scale.
# No zero-total-sample filter is added. There is no family-summing step.

required_packages <- c("tidyverse", "pheatmap", "RColorBrewer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "))
}
validate_key <- function(x, label) {
  if (anyNA(x) || any(!nzchar(trimws(as.character(x)))) || anyDuplicated(x)) {
    stop(label, " must contain unique, non-missing, non-empty IDs.")
  }
}
check_label_fields <- function(x) {
  # R converts these string escapes to TSV control characters before matching.
  if (anyNA(x) || any(!nzchar(x)) || any(grepl("||", x, fixed = TRUE)) ||
      any(grepl("[\t\r\n]", x))) {
    stop("Labels cannot be empty or contain ||, tabs or newlines.")
  }
}
write_heatmap_tsv <- function(x, row_info, path) {
  stopifnot(identical(rownames(x), row_info$Family))
  values <- as.data.frame(x, check.names = FALSE)
  for (j in seq_len(ncol(values))) values[[j]] <- sprintf("%.17g", x[, j])
  out <- data.frame(Phylum = row_info$Phylum, Family = row_info$Family,
                    values, check.names = FALSE)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
}

main <- function(input_dir = "inputs", output_dir = "results") {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  input_files <- c("bacterial_mag_tpm.tsv", "bacterial_mag_taxonomy.tsv", "sample_metadata.tsv")
  if (!all(file.exists(file.path(input_dir, input_files)))) stop("Missing input files.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (identical(input_dir, output_dir)) stop("Use a separate output directory.")
  message("[1/7] Reading and validating the three input tables")
  tpm_df <- read_tsv(file.path(input_dir, input_files[1]), show_col_types = FALSE)
  tax_df <- read_tsv(file.path(input_dir, input_files[2]), show_col_types = FALSE)
  meta_df <- read_tsv(file.path(input_dir, input_files[3]), show_col_types = FALSE)
  require_columns(tpm_df, "Genome", "TPM")
  require_columns(tax_df, c("Genome", "Phylum", "Family"), "Taxonomy")
  require_columns(meta_df, c("Sample", "Project", "Coral_Host_Species", "Ocean_Basin",
                            "Sample_Type", "Status"), "Metadata")
  for (x in list(tpm_df, tax_df, meta_df)) {
    if (nrow(readr::problems(x))) stop("Input parsing problems; inspect source tables.")
  }
  validate_key(tpm_df$Genome, "TPM Genome")
  validate_key(tax_df$Genome, "Taxonomy Genome")
  validate_key(meta_df$Sample, "Metadata Sample")
  if (any(as.character(meta_df$Sample) != str_trim(as.character(meta_df$Sample)))) {
    stop("Whitespace in sample IDs; resolve before joining.")
  }
  if (!nrow(tpm_df) || ncol(tpm_df) < 2L) stop("Empty TPM input.")
  tpm_values <- tpm_df %>% select(-Genome)
  if (!all(vapply(tpm_values, is.numeric, logical(1)))) stop("TPM columns must be numeric.")
  if (any(vapply(tpm_values, function(x) any(is.infinite(x) | x < 0, na.rm = TRUE), logical(1)))) {
    stop("TPM cannot contain negative or infinite values.")
  }
  if (length(setdiff(tpm_df$Genome, tax_df$Genome))) stop("Some MAG IDs lack taxonomy.")
  message("[2/7] Applying the original 51-family definition")

  target_families <- c(
    "Propionibacteriaceae", "Mycobacteriaceae", "Rhodobacteraceae", "JAQWUY01", 
    "Chlorobiaceae", "Alteromonadaceae", "Flavobacteriaceae", "Endozoicomonadaceae", 
    "Vibrionaceae", "WYAZ01", "UBA2968", "JAKSBE01", "Burkholderiaceae", 
    "Pyrinomonadaceae", "Flammeovirgaceae", "Halieaceae", "Amoebophilaceae", 
    "Rhodanobacteraceae", "UBA8231", "Nitrosopumilaceae", "Caminicellaceae", 
    "Ichthyobacteriaceae", "Oleiphilaceae", "Phormidesmiaceae", "Puniceispirillaceae", 
    "Woeseiaceae", "HK1", "Aldehydirespiratoraceae", "Pseudohongiellaceae", 
    "Sphingomonadaceae", "Nocardioidaceae", "Pseudomonadaceae", "UBA5704", 
    "JAKSDN01", "Desulfobacteraceae", "Methyloligellaceae", "JAAXHH01", 
    "Ketobacteraceae", "Cohaesibacteraceae", "Desulfolunaceae", "Arcobacteraceae", 
    "UBA6930", "Poseidoniaceae", "Alcanivoracaceae", 
    "Stappiaceae", "Cellvibrionaceae", "Micavibrionaceae", "Beggiatoaceae", 
    "Desulfovibrionaceae", "Magnetovibrionaceae", "Denitrovibrionaceae"
  )

  tax_clean <- tax_df %>%
    mutate(
      Display_Family = if_else(Family %in% target_families, Family, "Others (incl. Unknown)"),
      Phylum_Block = if_else(Display_Family == "Others (incl. Unknown)", "Z_Others", Phylum)
    ) %>%
    select(Genome, Phylum_Block, Display_Family)

  # ==============================================================================
  # 3. 元数据清洗与层级定义 (状态归并逻辑)
  # ==============================================================================
  message("[3/7] Applying the original metadata and status rules")

  ocean_levels <- c("Mediterranean Sea", "Pacific Ocean", "Red Sea", "Arabian Gulf", 
                    "Atlantic Ocean", "Caribbean Sea", "Indian Ocean")
  tissue_levels <- c("Coral tissue", "Coral mucus", "Coral skeleton", "Coral fragment", "Other")
  status_levels <- c("Healthy", "Bleached", "Disease", "Stressed", "Other")
  detail_levels <- c(
    "Healthy/Normal", "Bleached", 
    "Disease (Bleaching lesion)", "Disease(SCTLD)", "Disease(Dark Spot Syndrome)", "Disease (Pathogen-induced disease)", 
    "Stressed(Heat stress)", "Stressed(Acidification stress)", "Stressed(Combined)", "Other"
  )

  map_detail_status <- function(s) {
    s_norm <- str_trim(as.character(s))
    if (grepl("SCTLD", s_norm, ignore.case=T)) return("Disease(SCTLD)")
    if (grepl("Dark Spot", s_norm, ignore.case=T)) return("Disease(Dark Spot Syndrome)")
    if (grepl("Bleaching lesion", s_norm, ignore.case=T)) return("Disease (Bleaching lesion)")
    if (grepl("Pathogen", s_norm, ignore.case=T)) return("Disease (Pathogen-induced disease)")
    if (grepl("Bleached", s_norm, ignore.case=T)) return("Bleached")
    h <- grepl("Heat", s_norm, ignore.case=T); a <- grepl("Acidification", s_norm, ignore.case=T)
    if (h & a) return("Stressed(Combined)")
    if (h) return("Stressed(Heat stress)")
    if (a) return("Stressed(Acidification stress)")
    if (grepl("Healthy|Normal", s_norm, ignore.case=T)) return("Healthy/Normal")
    return("Other")
  }

  map_broad_status <- function(detail) {
    if (detail == "Healthy/Normal") return("Healthy")
    if (detail == "Bleached") return("Bleached")
    if (grepl("Disease", detail)) return("Disease")
    if (grepl("Stressed", detail)) return("Stressed")
    return("Other")
  }

  meta_clean <- meta_df %>%
    mutate(
      Project = str_trim(as.character(Project)),
      Ocean_Basin = case_when(
        grepl("Mediterranean", Ocean_Basin, ignore.case=T) ~ "Mediterranean Sea",
        grepl("Pacific", Ocean_Basin, ignore.case=T) ~ "Pacific Ocean",
        grepl("Red Sea", Ocean_Basin, ignore.case=T) ~ "Red Sea",
        grepl("Arabian", Ocean_Basin, ignore.case=T) ~ "Arabian Gulf",
        grepl("Atlantic", Ocean_Basin, ignore.case=T) ~ "Atlantic Ocean",
        grepl("Caribbean", Ocean_Basin, ignore.case=T) ~ "Caribbean Sea",
        grepl("Indian", Ocean_Basin, ignore.case=T) ~ "Indian Ocean",
        TRUE ~ "Pacific Ocean"
      ),
      Sample_Type_Fixed = case_when(
        grepl("tissue", Sample_Type, ignore.case=T) ~ "Coral tissue",
        grepl("mucus", Sample_Type, ignore.case=T) ~ "Coral mucus",
        grepl("skeleton", Sample_Type, ignore.case=T) ~ "Coral skeleton",
        grepl("fragment", Sample_Type, ignore.case=T) ~ "Coral fragment",
        TRUE ~ "Other"
      ),
      Detail_Status = sapply(Status, map_detail_status),
      Broad_Status = sapply(Detail_Status, map_broad_status),
      Species = ifelse(is.na(Coral_Host_Species)|Coral_Host_Species=="", "Unknown", str_trim(Coral_Host_Species))
    ) %>%
    mutate(
      Ocean_Basin = factor(Ocean_Basin, levels = ocean_levels),
      Sample_Type = factor(Sample_Type_Fixed, levels = tissue_levels),
      Broad_Status = factor(Broad_Status, levels = status_levels),
      Detail_Status = factor(Detail_Status, levels = detail_levels)
    )

  # ==============================================================================
  # 4. 数据聚合 (分组求均值与过滤)
  # ==============================================================================
  message("[4/7] Computing mean TPM by Project and the original group fields")

  check_label_fields(c(meta_clean$Project, meta_clean$Species,
                       as.character(meta_clean$Ocean_Basin), as.character(meta_clean$Sample_Type),
                       as.character(meta_clean$Broad_Status), as.character(meta_clean$Detail_Status),
                       tax_clean$Phylum_Block, tax_clean$Display_Family))

  plot_data_raw <- tpm_df %>%
    pivot_longer(-Genome, names_to = "Sample", values_to = "TPM") %>%
    inner_join(meta_clean, by = "Sample") %>%
    left_join(tax_clean, by = "Genome") %>%
    group_by(Project, Ocean_Basin, Species, Sample_Type, Broad_Status, Detail_Status, Phylum_Block, Display_Family) %>%
    summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

  # Preserve the original max() rule, but stop rather than silently discard
  # a whole family if future inputs produce undefined grouped means.
  if (any(!is.finite(plot_data_raw$Mean_TPM))) {
    stop("Non-finite grouped means: inspect missing-only groups before filtering.")
  }

  # 过滤低丰度行 (Mean TPM 最高值 < 1 且不是 Others 则剔除)
  active_families <- plot_data_raw %>%
    group_by(Display_Family) %>%
    summarise(max_tpm = max(Mean_TPM), .groups = "drop") %>%
    filter(max_tpm >= 1 | Display_Family == "Others (incl. Unknown)") %>%
    pull(Display_Family)

  plot_data <- plot_data_raw %>% filter(Display_Family %in% active_families)

  # 构造矩阵与 ID
  heatmap_matrix_raw <- plot_data %>%
    mutate(Col_ID = paste(Project, Species, Ocean_Basin, Sample_Type, Broad_Status, Detail_Status, sep = "||"),
           Row_ID = paste(Phylum_Block, Display_Family, sep = "||")) %>%
    select(Row_ID, Col_ID, Mean_TPM) %>%
    pivot_wider(names_from = Col_ID, values_from = Mean_TPM, values_fill = 0) %>%
    column_to_rownames("Row_ID")

  # 排序行与列
  row_info <- data.frame(Row_ID = rownames(heatmap_matrix_raw)) %>%
    separate(Row_ID, into = c("Phylum", "Family"), sep = "\\|\\|", remove = FALSE) %>%
    arrange(Phylum, Family)

  col_info <- data.frame(Col_ID = colnames(heatmap_matrix_raw)) %>%
    separate(Col_ID, into = c("Project", "Species", "Ocean", "Tissue", "Broad", "Detail"), sep = "\\|\\|", remove = FALSE) %>%
    mutate(Ocean = factor(Ocean, levels = ocean_levels), Tissue = factor(Tissue, levels = tissue_levels),
           Broad = factor(Broad, levels = status_levels), Detail = factor(Detail, levels = detail_levels)) %>%
    arrange(Ocean, Species, Tissue, Broad, Detail, Project)

  mat <- heatmap_matrix_raw[row_info$Row_ID, col_info$Col_ID, drop=FALSE]
  mat_ln <- log1p(as.matrix(mat))
  rownames(mat_ln) <- row_info$Family
  stopifnot(!anyDuplicated(rownames(mat_ln)), !anyDuplicated(colnames(mat_ln)),
            all(is.finite(mat_ln)), all(mat_ln >= 0),
            all(lengths(strsplit(colnames(mat_ln), "||", fixed = TRUE)) == 6L))

  # ==============================================================================
  # 5. 配置颜色与注释
  # ==============================================================================
  message("[5/7] Retaining the original annotations and fixed 0-15 scale")

  ann_col <- col_info %>% select(Detail, Broad, Tissue, Ocean) %>%
    rename(`Detail Status` = Detail, Status = Broad, Tissue = Tissue, `Ocean Basin` = Ocean)
  rownames(ann_col) <- col_info$Col_ID

  ann_colors <- list(
    "Ocean Basin"   = c("Mediterranean Sea"="#E41A1C", "Pacific Ocean"="#377EB8", "Red Sea"="#4DAF4A", 
                        "Arabian Gulf"="#984EA3", "Atlantic Ocean"="#FF7F00", "Caribbean Sea"="#FFFF33", "Indian Ocean"="#A65628"),
    "Tissue"        = c("Coral tissue"="#D95F02", "Coral mucus"="#1B9E77", "Coral skeleton"="#8D40CA", "Coral fragment"="#E72989", "Other"="#999999"),
    "Status"        = c("Healthy"="#7FC87E", "Bleached"="#FCC085", "Disease"="#E21A1C", "Stressed"="#3C15D3", "Other"="#999999"),
    "Detail Status" = c("Healthy/Normal"="#7EC87E", "Bleached"="#FCC085", "Disease (Bleaching lesion)"="#FC8D3B", 
                        "Disease(SCTLD)"="#800025", "Disease(Dark Spot Syndrome)"="#BC0025", 
                        "Disease (Pathogen-induced disease)"="#E26E7C", "Stressed(Heat stress)"="#756BB0", 
                        "Stressed(Acidification stress)"="#9E58F5", "Stressed(Combined)"="black", "Other"="#D9D9D9"),
    "Phylum"        = setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(unique(row_info$Phylum))), unique(row_info$Phylum))
  )

  # ==============================================================================
  # 6. 绘图输出
  # ==============================================================================
  message("[6/7] Drawing Project_Species labels")

  pdf_path <- file.path(output_dir, "Bacteria_Family_Project_Ocean_Heatmap.pdf")
  pdf(pdf_path, width = 24, height = 16, family = "Helvetica")
  heatmap <- pheatmap(mat_ln,
           color = colorRampPalette(c("white", "orange", "red", "black"))(100),
           breaks = seq(0, 15, length.out = 100),
           border_color = NA,
           cluster_cols = FALSE, cluster_rows = FALSE,
           annotation_col = ann_col, 
           annotation_row = data.frame(Phylum=row_info$Phylum, row.names=rownames(mat_ln)),
           annotation_colors = ann_colors,
           gaps_col = which(col_info$Ocean[-1] != col_info$Ocean[-nrow(col_info)] | 
                              col_info$Species[-1] != col_info$Species[-nrow(col_info)]),
           gaps_row = which(row_info$Phylum[-1] != row_info$Phylum[-nrow(row_info)]),
           labels_col = paste(col_info$Project, col_info$Species, sep = "_"),
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize = 10, fontsize_row = 9, angle_col = "45",
           cellwidth = 12, cellheight = 11,
           main = "Global Coral Bacterial MAGs: Ocean > Species > Tissue > Status",
           silent = TRUE
  )
  grid.newpage()
  grid.draw(heatmap$gtable)
  dev.off()


  message("[7/7] Exporting the single log1p source-data table")
  source_data_path <- file.path(output_dir, "Bacteria_Family_Heatmap_Log1p_SourceData.tsv")
  write_heatmap_tsv(mat_ln, row_info, source_data_path)
  roundtrip <- read_tsv(source_data_path, show_col_types = FALSE, name_repair = "minimal")
  roundtrip_matrix <- as.matrix(roundtrip[, -(1:2)])
  stopifnot(identical(names(roundtrip)[-(1:2)], colnames(mat_ln)),
            identical(roundtrip$Family, rownames(mat_ln)),
            isTRUE(all.equal(unname(roundtrip_matrix), unname(mat_ln), tolerance = 1e-14)))
  matched_samples <- intersect(names(tpm_df)[-1L], meta_clean$Sample)
  message("Complete: ", nrow(mat_ln), " rows x ", ncol(mat_ln), " groups; ",
          length(matched_samples), " metadata-matched samples")
  message("TPM samples without metadata (excluded as in the original): ",
          length(setdiff(names(tpm_df)[-1L], meta_clean$Sample)))
  message("Source data: ", source_data_path)
  message("PDF: ", pdf_path)
  invisible(list(source_matrix = mat_ln, mean_TPM = as.matrix(mat),
    row_info = row_info, col_info = col_info, plot_data_raw = plot_data_raw,
    active_families = active_families, target_families = target_families,
    annotation_col = ann_col, annotation_colors = ann_colors,
    heatmap_colors = colorRampPalette(c("white", "orange", "red", "black"))(100),
    heatmap_breaks = seq(0, 15, length.out = 100), final_grob = heatmap$gtable,
    matched_samples = matched_samples, pdf_width_inches = 24, pdf_height_inches = 16,
    input_md5 = tools::md5sum(file.path(input_dir, input_files)),
    source_data_path = source_data_path, pdf_path = pdf_path))
}

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) args <- character()
if (length(args) > 2L) {
  stop("Usage: Rscript run_heatmap.R [input_dir] [output_dir]")
}
analysis <- main(if (length(args) >= 1L) args[[1L]] else "inputs",
                 if (length(args) >= 2L) args[[2L]] else "results")
