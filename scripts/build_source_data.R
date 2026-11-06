#!/usr/bin/env Rscript
# Assemble all figure values and statistics without recomputing or rounding them.
# Usage: Rscript scripts/build_source_data.R [analysis_root] [output_directory]
# analysis_root is the repository or an output tree with the same module layout.
# master_source_data.tsv stores one original table cell per row. raw_value is
# authoritative; numeric_value is a convenience column. Source row/column names
# permit a lossless connection to each copied, figure-specific source table.
args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", ".UTF-8")
flag <- grep("^--file=", commandArgs(), value = TRUE)
release <- if (length(flag)) dirname(dirname(normalizePath(sub("^--file=", "", flag[1]), winslash = "/"))) else getwd()
root <- normalizePath(if (length(args)) args[1] else release, winslash = "/", mustWork = TRUE)
output <- if (length(args) >= 2L) args[2] else file.path(root, "source_data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
spec <- data.frame(module = character(), panel = character(), file = character(), role = character(), measurement = character())
add <- function(module, panel, file, role, measurement) {
  spec <<- rbind(spec, data.frame(module, panel, file, role, measurement))
}
add("alpha_diversity", "DNA_RNA", "alpha_diversity/results/alpha_diversity_source.tsv", "figure_values", "Shannon index; Chao1 richness; count-based rarefaction")
add("alpha_diversity", "DNA_RNA", "alpha_diversity/results/alpha_diversity_statistics.tsv", "statistics", "Wilcoxon P and within-stratum BH adjusted P")
add("amg_profiles", "heatmap", "amg_profiles/results/amg_heatmap_log1p.tsv", "figure_values", "ln(1 + mean positive-record TPM)")
for (kind in c("gene", "group"))
  add("amg_profiles", "heatmap", paste0("amg_profiles/results/", kind, "_annotations.csv"), "annotation", "Displayed annotations and order")
for (workflow in c("nested", "factor_evaluation")) {
  base <- file.path("beta_diversity", workflow, "Source_Data")
  if (!dir.exists(file.path(root, base))) base <- file.path("beta_diversity", workflow, "source_data")
  files <- list.files(file.path(root, base), pattern = "\\.tsv$", recursive = TRUE)
  files <- files[!grepl("dictionary|manifest|checksum|permutation_execution|figure_index", files, ignore.case = TRUE)]
  if (!length(files)) stop("No Beta figure source tables in ", base)
  for (f in files) add("beta_diversity", paste(workflow, tools::file_path_sans_ext(f), sep = ":"),
                      gsub("\\\\", "/", file.path(base, f)),
                      if (grepl("points", f)) "figure_values" else "statistics",
                      "sqrt(TPM) Bray-Curtis PCoA; PERMANOVA/dbRDA/PERMDISP; see module dictionary")
}
add("lefse_analysis", "DNA_RNA", "lefse_analysis/results/source_data.tsv", "figure_values_and_statistics", "CPM; signed LDA score; KW P; BH within stratum; explicit selection/display flags")
add("lyso_analysis", "lifestyle", "lyso_analysis/results/phage_lifestyle_stats.csv", "figure_values", "Lifestyle composition; see original column names")
for (assay in c("dna", "rna", "bacterial_family")) {
  filename <- switch(assay, dna = "DNA_Heatmap_Log1p_SourceData.tsv", rna = "RNA_Heatmap_Log1p_SourceData.tsv", bacterial_family = "Bacteria_Family_Heatmap_Log1p_SourceData.tsv")
  add("relative_abundance", assay, paste0("relative_abundance/", assay, "/results/", filename), "figure_values", "ln(1 + mean TPM); see module for averaging denominator")
}
add("viral_network", "a_f", "viral_network/results/network_analysis_master_results.tsv", "figure_values_and_statistics", "Metric-specific units, sample counts, intervals and P values in original columns")
for (f in list.files(file.path(root, "viral_network", "figure_source_data"), pattern = "\\.csv$"))
  add("viral_network", tools::file_path_sans_ext(f), paste0("viral_network/figure_source_data/", f), "panel_source", "Network nodes, edges, topology, resampling and display audit")
add("virus_host_ratio", "VHR", "virus_host_ratio/results/source_data.tsv", "figure_values_and_statistics", "Linkage-informed viral/host abundance signal ratio; log10 ratio; within-context and global BH P")
read_source <- function(path) {
  x <- read.table(path, header = TRUE, sep = if (grepl("\\.csv$", path)) "," else "\t",
             quote = "\"", comment.char = "", colClasses = "character", check.names = FALSE,
             na.strings = character(), fill = FALSE, fileEncoding = "UTF-8")
  x
}
# Double embedded quotes so standard TSV readers and spreadsheets can recover
# JSON-like status annotations without a special backslash-escape setting.
write_tsv <- function(x, name) write.table(x, file.path(output, name), sep = "\t", quote = TRUE,
                                         qmethod = "double", row.names = FALSE, na = "")
pick <- function(x, candidates) {
  z <- intersect(candidates, names(x))
  if (!length(z)) rep("", nrow(x)) else x[[z[1]]]
}
master <- list(); statistics <- list(); dictionaries <- list(); indices <- list()
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]; path <- file.path(root, s$file)
  if (!file.exists(path)) stop("Missing required figure source: ", s$file)
  x <- read_source(path)
  record <- pick(x, c("Record_ID", "record_id", "SampleID", "Sample", "Comparison_ID", "Family", "Gene", "Taxon", "ID"))
  record[record == ""] <- paste0("row:", which(record == ""))
  cells <- lapply(seq_along(x), function(j) {
    value <- x[[j]]
    numeric <- grepl("^[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$", value, perl = TRUE)
    # A valid header-only source table contributes no cells, while its fields,
    # dimensions and checksum remain in the dictionary and table index.
    data.frame(module = rep(s$module, nrow(x)), panel = rep(s$panel, nrow(x)),
               source_file = rep(s$file, nrow(x)),
               source_row = seq_len(nrow(x)), record_id = record, field = rep(names(x)[j], nrow(x)),
               raw_value = value, numeric_value = ifelse(numeric, value, ""), stringsAsFactors = FALSE)
  })
  master[[i]] <- do.call(rbind, cells)
  dictionaries[[i]] <- data.frame(module = s$module, panel = s$panel, source_file = s$file,
                                  field = names(x), measurement_context = s$measurement)
  # Include prefixed fields such as Dispersion_P and Dispersion_Q as well as
  # P/Q at the start of a field. Resolution fields describe the experiment,
  # rather than an observed test P value, and stay only in the full source.
  pcols <- names(x)[grepl("(^[pq]$|(^|[._])[pq]([._]|$)|p.?value|p.?adjust|^KW_P|Pr\\(>F\\))", names(x), ignore.case = TRUE)]
  pcols <- pcols[!grepl("format|signif|label|resolution|minimum", pcols, ignore.case = TRUE)]
  for (p in pcols) {
    values <- suppressWarnings(as.numeric(x[[p]]))
    keep <- is.finite(values) & values >= 0 & values <= 1
    if (!any(keep)) next
    statistics[[length(statistics) + 1L]] <- data.frame(module = s$module, panel = s$panel,
      source_file = s$file, source_row = which(keep), record_id = record[keep],
      assay = pick(x, c("Assay", "Virome", "Virus"))[keep],
      comparison = pick(x, c("comparison", "Comparison", "Condition", "Context", "Comparison_ID"))[keep],
      entity = pick(x, c("Family", "entity_1", "Factor", "Index", "Metric", "metric"))[keep],
      p_value_field = p, p_value = x[[p]][keep],
      effect_field = paste(intersect(c("R2", "R2_adj", "F", "LDA_Signed", "estimate", "Delta_Median"), names(x)), collapse = ";"),
      stringsAsFactors = FALSE)
  }
  copied <- file.path(output, "tables", s$file)
  dir.create(dirname(copied), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(path, copied, overwrite = TRUE)) stop("Could not copy ", s$file)
  indices[[i]] <- data.frame(s, rows = nrow(x), columns = ncol(x),
                             sha256 = digest::digest(file = path, algo = "sha256"),
                             source_data_copy = paste0("tables/", s$file))
}
write_tsv(do.call(rbind, master), "master_source_data.tsv")
statistics_table <- if (length(statistics)) do.call(rbind, statistics) else
  data.frame(module = character(), panel = character(), source_file = character(),
             source_row = integer(), record_id = character(), assay = character(),
             comparison = character(), entity = character(), p_value_field = character(),
             p_value = character(), effect_field = character())
write_tsv(statistics_table, "statistical_summary.tsv")
write_tsv(do.call(rbind, indices), "source_table_index.tsv")
write_tsv(do.call(rbind, dictionaries), "field_dictionary.tsv")
writeLines(c(
  "# Figure Source Data", "",
  "master_source_data.tsv is the combined result table: one original source-table cell per row. raw_value preserves the source text (including numerical precision); numeric_value is populated only for numerical literals. Missing values retain the original source notation. source_row is 1-based and excludes the header.", "",
  "Combined TSV files use standard doubled-quote escaping. A standard CSV/TSV reader recovers embedded quotation marks in status annotations. Import identifier and raw_value columns as text to preserve their original spelling and precision in spreadsheet software.", "",
  "source_table_index.tsv maps each module/panel to its original table, dimensions, units/measurement context and SHA-256. tables/ contains complete figure-specific source tables. field_dictionary.tsv lists all original fields; module READMEs and the Beta dictionaries define statistical methods and exclusions.", "",
  "statistical_summary.tsv gathers reported P-value fields with module/panel/row provenance. It is an index of reported tests, not a new family of hypothesis tests; repeated values across figure summaries and detailed panel tables are intentional. Adjustments apply only to the families defined by each analysis. Selection flags in the original LEfSe/VHR/network tables must be respected.", "",
  "No across-module mean, extra correction or unit conversion is applied. Figure values, statistics and annotations keep their original scientific meanings. Manuscript figure numbers should be assigned after the final figure order is set.", "",
  "Regenerate with: Rscript scripts/build_source_data.R <analysis_root> <source_data_directory>"
), file.path(output, "README.md"))
cat("Source Data complete:", nrow(spec), "tables;", sum(vapply(master, nrow, integer(1))), "cells.\n")
