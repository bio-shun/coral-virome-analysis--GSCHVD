# Coral virome alpha diversity: count-based DNA and RNA analysis.
# Usage: Rscript scripts/run_alpha_diversity.R [input_directory] [output_directory]
# Default directories are resolved relative to this module, independent of getwd().
#
# Analysis: integer counts; seed 123; sample total >= 1,000; rarefaction to the
# minimum retained depth; Shannon and Chao1; Healthy-vs-status Wilcoxon tests.
# Eligible groups contain >= 3 samples per status within project, basin, tissue
# and coral species. BH correction is applied within each stratum and index.
# DNA uses seed 123. The recorded RNA random state reproduces the original RNA
# rarefaction independently of the DNA matrix; see README for its provenance.
# Use plot_alpha_diversity.R to render both calculated assays.
args <- commandArgs(trailingOnly = TRUE)
script_flag <- grep('^--file=', commandArgs(), value = TRUE)
module_root <- if (length(script_flag)) dirname(dirname(normalizePath(sub('^--file=', '', script_flag[[1]]), winslash = '/', mustWork = TRUE))) else getwd()
input_dir <- normalizePath(if (length(args) >= 1) args[[1]] else file.path(module_root, 'inputs'), winslash = '/', mustWork = TRUE)
output_root <- if (length(args) >= 2) args[[2]] else file.path(module_root, 'results')
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = '/', mustWork = TRUE)
new_count_file <- file.path(input_dir, 'dna_counts.tsv')
metadata_file <- file.path(input_dir, 'sample_metadata.tsv')
rna_count_file <- file.path(input_dir, 'rna_counts.tsv')
rna_metadata_file <- file.path(input_dir, 'rna_sample_metadata.tsv')
rna_state_file <- file.path(input_dir, 'rna_rarefaction_rng_state.tsv')
setwd(output_root)
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(RColorBrewer)
  library(stringr)
})
set.seed(123)

# Harmonize metadata labels before selecting eligible ecological strata.
global_colors <- c(
  "Healthy"                          = "#33a02c",
  "Bleached"                         = "#e31a1c",
  "Diseased(SCTLD)"                  = "#542788",
  "Diseased(black band disease)"     = "#40004b",
  "Diseased(white plague disease)"   = "#9970ab",
  "Diseased(Dark Spot Syndrome)"     = "#762a83",
  "Diseased(Bleaching lesion)"       = "#FD8D3C",
  "Diseased(Pathogen-induced)"       = "#8073ac",
  "Stressed(Heat stress)"            = "#ff7f00",
  "Stressed(Acidification stress)"   = "#fdbf6f",
  "Stressed(Combined)"               = "#b15928",
  "Stressed(Cold stress)"            = "#1f78b4",
  "Other"                            = "#999999"
)
get_refined_status <- function(s) {
  s_norm <- str_trim(as.character(s))
  if (is.na(s_norm) || s_norm == "") {
    return("Other")
  }
  if (grepl("Mix", s_norm, ignore.case = TRUE)) {
    return("Other")
  }
  if (grepl("Healthy|Normal", s_norm, ignore.case = TRUE)) {
    return("Healthy")
  }
  if (grepl("Bleached", s_norm, ignore.case = TRUE)) {
    return("Bleached")
  }
  if (grepl("SCTLD", s_norm, ignore.case = TRUE)) {
    return("Diseased(SCTLD)")
  }
  if (grepl("black band", s_norm, ignore.case = TRUE)) {
    return("Diseased(black band disease)")
  }
  if (grepl("white plague", s_norm, ignore.case = TRUE)) {
    return("Diseased(white plague disease)")
  }
  if (grepl("Dark Spot", s_norm, ignore.case = TRUE)) {
    return("Diseased(Dark Spot Syndrome)")
  }
  if (grepl("Bleaching lesion", s_norm, ignore.case = TRUE)) {
    return("Diseased(Bleaching lesion)")
  }
  if (grepl("Pathogen", s_norm, ignore.case = TRUE)) {
    return("Diseased(Pathogen-induced)")
  }
  has_heat <- grepl(
    "Heat",
    s_norm,
    ignore.case = TRUE
  )
  has_acid <- grepl(
    "Acidification",
    s_norm,
    ignore.case = TRUE
  )
  has_cold <- grepl(
    "Cold",
    s_norm,
    ignore.case = TRUE
  )
  if (has_heat && has_acid) {
    return("Stressed(Combined)")
  }
  if (has_heat) {
    return("Stressed(Heat stress)")
  }
  if (has_acid) {
    return("Stressed(Acidification stress)")
  }
  if (has_cold) {
    return("Stressed(Cold stress)")
  }
  if (grepl("Disease|Diseased", s_norm, ignore.case = TRUE)) {
    return("Other")
  }
  return("Other")
}
avila_projects <- c(
  "PRJNA405605",
  "PRJNA405606",
  "PRJNA405607",
  "PRJNA405608",
  "PRJNA405609",
  "PRJNA405610",
  "PRJNA405611",
  "PRJNA405612",
  "PRJNA405613",
  "PRJNA405614",
  "PRJNA405615",
  "PRJNA405616",
  "PRJNA405617",
  "PRJNA405618",
  "PRJNA405952",
  "PRJNA405953",
  "PRJNA405954",
  "PRJNA405955",
  "PRJNA405956",
  "PRJNA405957",
  "PRJNA405958",
  "PRJNA405995",
  "PRJNA405996",
  "PRJNA405997"
)
run_alpha_analysis <- function(
    analysis_type,
    count_file,
    meta_file
) {
  cat("\n\n")
  cat("============================================================\n")
  cat("Starting analysis:", analysis_type, "\n")
  cat("============================================================\n")
  if (!file.exists(count_file)) {
    warning(
      paste0(
        analysis_type,
        ": count file not found: ",
        count_file
      )
    )
    return(NULL)
  }
  if (!file.exists(meta_file)) {
    warning(
      paste0(
        analysis_type,
        ": metadata file not found: ",
        meta_file
      )
    )
    return(NULL)
  }
  counts <- read_tsv(
    count_file,
    show_col_types = FALSE
  )
  counts <- counts %>%
    column_to_rownames(
      var = colnames(counts)[1]
    )
  meta <- read.table(
    meta_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(meta)[1] <- "SampleID"
  required_cols <- c(
    "SampleID",
    "Project",
    "Coral_Host_Species",
    "Ocean_Basin",
    "Sample_Type",
    "Status"
  )
  missing_cols <- setdiff(
    required_cols,
    colnames(meta)
  )
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        analysis_type,
        ": metadata is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  clean_meta <- meta %>%
    mutate(
      SampleID = str_trim(
        as.character(SampleID)
      ),
      Project = ifelse(
        is.na(Project) |
          str_trim(as.character(Project)) == "",
        "Unknown_Project",
        str_trim(
          as.character(Project)
        )
      ),
      Analysis_Project = ifelse(
        Project %in% avila_projects,
        "Avila2021",
        Project
      ),
      Species = ifelse(
        is.na(Coral_Host_Species) |
          str_trim(as.character(Coral_Host_Species)) == "",
        "Unknown",
        str_trim(
          as.character(Coral_Host_Species)
        )
      ),
      Ocean_Basin = ifelse(
        is.na(Ocean_Basin) |
          str_trim(as.character(Ocean_Basin)) == "",
        "Unknown Ocean",
        str_trim(
          as.character(Ocean_Basin)
        )
      ),
      Clean_Type = str_to_sentence(
        gsub(
          "^Coral\\s+",
          "",
          ifelse(
            is.na(Sample_Type),
            "Other",
            str_trim(
              as.character(Sample_Type)
            )
          )
        )
      ),
      Status_Group = sapply(
        Status,
        get_refined_status
      )
    ) %>%
  filter(
    !Status %in% c(
      "Treated/Recovery",
      "Unkonwn",
      "Unknown"
    )
  ) %>%
  filter(
    Status_Group != "Other"
  )
  project_mapping <- clean_meta %>%
    distinct(
      Project,
      Analysis_Project
    ) %>%
    arrange(
      Analysis_Project,
      Project
    )
  write.csv(
    project_mapping,
    paste0(
      "QC_",
      analysis_type,
      "_Project_to_AnalysisProject.csv"
    ),
    row.names = FALSE
  )
  valid_samples <- intersect(
    colnames(counts),
    clean_meta$SampleID
  )
  cat(
    analysis_type,
    "- matched samples:",
    length(valid_samples),
    "\n"
  )
  if (length(valid_samples) < 3) {
    warning(
      paste0(
        analysis_type,
        ": fewer than three matched samples."
      )
    )
    return(NULL)
  }
  count_subset <- counts[
    ,
    valid_samples,
    drop = FALSE
  ]
  count_subset <- count_subset[
    rowSums(count_subset) > 0,
    ,
    drop = FALSE
  ]
  if (nrow(count_subset) == 0) {
    warning(
      paste0(
        analysis_type,
        ": no non-zero viral features."
      )
    )
    return(NULL)
  }
  sample_sums <- colSums(
    count_subset
  )
  valid_depth_samples <- names(
    sample_sums[
      sample_sums >= 1000
    ]
  )
  cat(
    analysis_type,
    "- samples with total count >= 1000:",
    length(valid_depth_samples),
    "\n"
  )
  if (length(valid_depth_samples) < 3) {
    warning(
      paste0(
        analysis_type,
        ": fewer than three samples after depth filtering."
      )
    )
    return(NULL)
  }
  rarefy_depth <- min(
    sample_sums[
      valid_depth_samples
    ]
  )
  cat(
    analysis_type,
    "- rarefaction depth:",
    rarefy_depth,
    "\n"
  )
  rarefied_mat <- rrarefy(
    t(
      count_subset[
        ,
        valid_depth_samples,
        drop = FALSE
      ]
    ),
    sample = rarefy_depth
  )
  alpha_data <- data.frame(
    SampleID = rownames(
      rarefied_mat
    ),
    Shannon = diversity(
      rarefied_mat,
      index = "shannon"
    ),
    Chao1 = estimateR(
      rarefied_mat
    )[
      "S.chao1",
    ],
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      clean_meta,
      by = "SampleID"
    )
  status_counts <- alpha_data %>%
    group_by(
      Analysis_Project,
      Ocean_Basin,
      Clean_Type,
      Species,
      Status_Group
    ) %>%
    summarise(
      n = n(),
      .groups = "drop"
    )
  write.csv(
    status_counts,
    paste0(
      "QC_",
      analysis_type,
      "_Status_Sample_Counts_ByProject.csv"
    ),
    row.names = FALSE
  )
  valid_statuses <- status_counts %>%
    filter(
      n >= 3
    )
  valid_strata <- valid_statuses %>%
    group_by(
      Analysis_Project,
      Ocean_Basin,
      Clean_Type,
      Species
    ) %>%
    summarise(
      Has_Healthy = any(
        Status_Group == "Healthy"
      ),
      N_Status = n_distinct(
        Status_Group
      ),
      .groups = "drop"
    ) %>%
    filter(
      Has_Healthy,
      N_Status >= 2
    )
  write.csv(
    valid_strata,
    paste0(
      "QC_",
      analysis_type,
      "_Eligible_Comparison_Strata.csv"
    ),
    row.names = FALSE
  )
  final_data <- alpha_data %>%
    semi_join(
      valid_statuses,
      by = c(
        "Analysis_Project",
        "Ocean_Basin",
        "Clean_Type",
        "Species",
        "Status_Group"
      )
    ) %>%
    semi_join(
      valid_strata,
      by = c(
        "Analysis_Project",
        "Ocean_Basin",
        "Clean_Type",
        "Species"
      )
    )
  if (nrow(final_data) == 0) {
    warning(
      paste0(
        analysis_type,
        ": no eligible Healthy-vs-status comparison."
      )
    )
    return(NULL)
  }
  project_overlap <- final_data %>%
    distinct(
      Analysis_Project,
      Ocean_Basin,
      Clean_Type,
      Species
    ) %>%
    group_by(
      Ocean_Basin,
      Clean_Type,
      Species
    ) %>%
    summarise(
      N_Project = n_distinct(
        Analysis_Project
      ),
      .groups = "drop"
    )
  final_data <- final_data %>%
    left_join(
      project_overlap,
      by = c(
        "Ocean_Basin",
        "Clean_Type",
        "Species"
      )
    ) %>%
    mutate(
      Display_Species = ifelse(
        N_Project > 1,
        paste0(
          Analysis_Project,
          "_",
          Species
        ),
        Species
      )
    )
  repeated_project_check <- final_data %>%
    filter(
      N_Project > 1
    ) %>%
    distinct(
      Analysis_Project,
      Ocean_Basin,
      Clean_Type,
      Species,
      Display_Species
    ) %>%
    arrange(
      Ocean_Basin,
      Clean_Type,
      Species,
      Analysis_Project
    )
  write.csv(
    repeated_project_check,
    paste0(
      "QC_",
      analysis_type,
      "_Species_Split_ByProject.csv"
    ),
    row.names = FALSE
  )
  write.csv(
    final_data,
    paste0(
      "SourceData_",
      analysis_type,
      "_Alpha_Diversity_ProjectAware.csv"
    ),
    row.names = FALSE
  )
  plot_long <- final_data %>%
    pivot_longer(
      cols = c(
        Shannon,
        Chao1
      ),
      names_to = "Index",
      values_to = "Value"
    )
  stats_results <- plot_long %>%
    group_by(
      Analysis_Project,
      Ocean_Basin,
      Clean_Type,
      Species,
      Index
    ) %>%
    group_modify(
      ~ {
        sub_df <- .x
        if (
          "Healthy" %in% sub_df$Status_Group &&
          n_distinct(sub_df$Status_Group) > 1
        ) {
          res <- compare_means(
            Value ~ Status_Group,
            data = sub_df,
            method = "wilcox.test",
            ref.group = "Healthy"
          )
          res <- res %>%
            mutate(
              p.adj.BH = p.adjust(
                p,
                method = "BH"
              ),
              p.adj.signif = case_when(
                p.adj.BH < 0.001 ~ "***",
                p.adj.BH < 0.01 ~ "**",
                p.adj.BH < 0.05 ~ "*",
                p.adj.BH < 0.1 ~ "ns (trend)",
                TRUE ~ "ns"
              )
            )
          return(res)
        }
        return(
          tibble()
        )
      }
    ) %>%
    ungroup()
  write.csv(
    stats_results,
    paste0(
      "SourceData_",
      analysis_type,
      "_Alpha_Pvalues_ProjectAware_BH.csv"
    ),
    row.names = FALSE
  )
  display_order <- final_data %>%
    distinct(
      Ocean_Basin,
      Clean_Type,
      Species,
      Analysis_Project,
      Display_Species
    ) %>%
    arrange(
      Ocean_Basin,
      Clean_Type,
      Species,
      Analysis_Project
    ) %>%
    pull(
      Display_Species
    ) %>%
    unique()
  plot_long <- plot_long %>%
    mutate(
      Display_Species = factor(
        Display_Species,
        levels = display_order
      )
    )
  dodge_val <- 0.8
  box_width <- 0.6
  p <- ggplot(
    plot_long,
    aes(
      x = Display_Species,
      y = Value,
      fill = Status_Group
    )
  ) +
    geom_boxplot(
      width = box_width,
      outlier.shape = NA,
      alpha = 0.8,
      lwd = 0.5,
      position = position_dodge(
        width = dodge_val,
        preserve = "single"
      )
    ) +
    geom_point(
      aes(
        group = Status_Group
      ),
      position = position_jitterdodge(
        jitter.width = 0.1,
        dodge.width = dodge_val
      ),
      size = 2.5,
      alpha = 0.5,
      shape = 21,
      color = "black",
      stroke = 0.3
    ) +
    facet_grid(
      Index ~ Ocean_Basin + Clean_Type,
      scales = "free",
      space = "free_x"
    ) +
    scale_fill_manual(
      values = global_colors,
      drop = FALSE
    ) +
    theme_bw() +
    labs(
      x = "Coral taxon",
      y = paste0(
        analysis_type,
        " viral alpha diversity"
      ),
      fill = "Status"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "italic",
        size = 10,
        color = "black"
      ),
      strip.text = element_text(
        size = 9,
        face = "bold"
      ),
      strip.background = element_rect(
        fill = "gray98"
      ),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
  n_display <- n_distinct(
    plot_long$Display_Species
  )
  plot_width <- n_display * 0.8 + 6
  cat("\n")
  cat(
    analysis_type,
    "- eligible samples:",
    nrow(final_data),
    "\n"
  )
  cat(
    analysis_type,
    "- eligible Project/ecological strata:",
    nrow(valid_strata),
    "\n"
  )
  cat(
    analysis_type,
    "- repeated species requiring Project_Species labels:",
    n_distinct(repeated_project_check$Species),
    "\n"
  )
  cat(
    analysis_type,
    "- statistical comparisons:",
    nrow(stats_results),
    "\n"
  )
  return(
    list(
      plot = p,
      width = plot_width,
      final_data = final_data,
      stats = stats_results,
      valid_strata = valid_strata,
      rarefy_depth = rarefy_depth,
      matched_samples = length(valid_samples),
      samples_after_depth = length(valid_depth_samples),
      split_species = repeated_project_check
    )
  )
}

# Validate both count matrices without rounding, scaling or substituting TPM.
validate_counts <- function(path, expected_features = NULL) {
  counts <- read_tsv(path, show_col_types = FALSE, name_repair = 'minimal')
  values <- as.matrix(counts[, -1])
  stopifnot(ncol(counts) > 1L, nrow(counts) > 0L,
            is.null(expected_features) || nrow(counts) == expected_features,
            !anyNA(counts[[1]]), !anyDuplicated(counts[[1]]),
            !anyDuplicated(names(counts)), is.numeric(values),
            all(is.finite(values)), all(values >= 0),
            all(values == floor(values)))
}
validate_counts(new_count_file, 24083L)
validate_counts(rna_count_file)

# Fix the R random generators and seed immediately before DNA rarefaction.
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
set.seed(123)
res_dna <- run_alpha_analysis('DNA', new_count_file, metadata_file)
stopifnot(!is.null(res_dna))

# Restore the RNA stream captured after the original seed-123 DNA analysis.
# This keeps historical RNA rarefaction independent of the revised DNA counts.
rna_rng_state <- scan(rna_state_file, what = integer(), quiet = TRUE)
stopifnot(length(rna_rng_state) == 626L, rna_rng_state[[1]] == 10403L,
          !anyNA(rna_rng_state))
assign('.Random.seed', rna_rng_state, envir = .GlobalEnv)
res_rna <- run_alpha_analysis('RNA', rna_count_file, rna_metadata_file)
stopifnot(!is.null(res_rna))

source_columns <- c('SampleID', 'Project', 'Analysis_Project', 'Species', 'Ocean_Basin', 'Clean_Type', 'Status', 'Status_Group', 'Display_Species', 'N_Project', 'Shannon', 'Chao1')
source_all <- bind_rows(mutate(select(res_dna$final_data, all_of(source_columns)), Assay = 'DNA', .before = 1), mutate(select(res_rna$final_data, all_of(source_columns)), Assay = 'RNA', .before = 1))
write.table(source_all, file.path(output_root, 'alpha_diversity_source.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA', fileEncoding = 'UTF-8')
stats_rna <- res_rna$stats
stats_rna$p.format <- as.character(stats_rna$p.format)
stats_all <- bind_rows(mutate(res_dna$stats, Assay = 'DNA', .before = 1), mutate(stats_rna, Assay = 'RNA', .before = 1))
write.table(stats_all, file.path(output_root, 'alpha_diversity_statistics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA', fileEncoding = 'UTF-8')
run_summary <- data.frame(
  Dataset = c('DNA', 'RNA'),
  Matched_samples = c(res_dna$matched_samples, res_rna$matched_samples),
  Samples_after_depth_filter = c(res_dna$samples_after_depth, res_rna$samples_after_depth),
  Rarefaction_depth = c(res_dna$rarefy_depth, res_rna$rarefy_depth),
  Plotted_samples = c(nrow(res_dna$final_data), nrow(res_rna$final_data)),
  Eligible_strata = c(nrow(res_dna$valid_strata), nrow(res_rna$valid_strata)),
  Statistical_comparisons = c(nrow(res_dna$stats), nrow(res_rna$stats)))
write.table(run_summary, file.path(output_root, 'run_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_root, 'sessionInfo.txt'))
print(run_summary)
message('Analysis completed: ', output_root)
