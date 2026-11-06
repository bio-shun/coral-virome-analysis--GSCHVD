#!/usr/bin/env Rscript

# Complete, standalone workflow for reproducing the fixed manuscript dataset.
# This script is not a general-purpose network-analysis package.
# Statistical and figure bodies are generated from 01_network_analysis.R and
# 02_network_figure.R; the fixed-result assertions are intentional.
#
# Usage:
#   Rscript run_network.R plot-only <analysis_package_dir> <output_dir> <lefse_source_data.tsv>
#   Rscript run_network.R full <input_dir> <output_dir> <lefse_source_data.tsv>
# LEfSe is a separate upstream calculation. Its current output supplies only
# the post hoc node annotation; it never enters network construction.

derive_lda_annotation <- function(lefse_source, snapshot_dir) {
  # Match the complete biological context rather than a potentially reordered
  # comparison identifier. Fail explicitly if the current analysis has changed.
  lefse <- read.delim(lefse_source, check.names = FALSE, stringsAsFactors = FALSE,
                      na.strings = "")
  required <- c("Virome", "Condition", "Family", "Comparison_ID", "Analysis_Project",
                "Ocean_Basin", "Coral_Host_Species", "Sample_Type", "Locations",
                "Enriched_In", "LDA_Signed", "LDA_Abs", "KW_P", "KW_P_BH_within_stratum",
                "Healthy_mean_CPM", "Condition_mean_CPM", "Healthy_n", "Condition_n",
                "Retained", "Display_Flag", "Direction_Matches_Group_Mean",
                "Healthy_Samples", "Condition_Samples")
  if (!all(required %in% names(lefse)))
    stop("LEfSe source data lack required fields: ", paste(setdiff(required, names(lefse)), collapse = ", "))
  matched <- which(lefse$Virome == "DNA" & lefse$Condition == "SCTLD" &
                     lefse$Family == "Microviridae" & lefse$Analysis_Project == "PRJNA925892" &
                     lefse$Ocean_Basin == "Atlantic Ocean" &
                     lefse$Coral_Host_Species == "Orbicella faveolata" &
                     lefse$Sample_Type == "Coral tissue")
  if (length(matched) != 1L) stop("Expected exactly one current LEfSe Microviridae annotation row; found ", length(matched))
  row <- lefse[matched, , drop = FALSE]
  numeric_fields <- c("LDA_Signed", "LDA_Abs", "KW_P", "KW_P_BH_within_stratum",
                      "Healthy_mean_CPM", "Condition_mean_CPM", "Healthy_n", "Condition_n")
  if (!all(is.finite(as.numeric(unlist(row[1, numeric_fields], use.names = FALSE)))) ||
      !isTRUE(row$Retained[[1]]) || !isTRUE(row$Display_Flag[[1]]) ||
      !isTRUE(row$Direction_Matches_Group_Mean[[1]]) ||
      row$LDA_Abs[[1]] < 2.5 || row$Enriched_In[[1]] != "SCTLD" ||
      row$Healthy_n[[1]] != 20L || row$Condition_n[[1]] != 69L)
    stop("Current LEfSe annotation does not satisfy the fixed manuscript context and display rule.")
  columns <- c("Family", "Comparison_ID", "Analysis_Project", "Ocean_Basin",
               "Coral_Host_Species", "Sample_Type", "Locations", "Enriched_In",
               "LDA_Signed", "LDA_Abs", "KW_P", "KW_P_BH_within_stratum",
               "Healthy_mean_CPM", "Condition_mean_CPM", "Healthy_n", "Condition_n")
  annotation <- row[, columns, drop = FALSE]
  names(annotation)[1] <- "LDA_Family"
  annotation$LDA_Display_Threshold <- 2.5
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
  source_snapshot <- file.path(snapshot_dir, "lefse_source_data.tsv")
  if (!identical(tolower(normalizePath(lefse_source, winslash = "/", mustWork = TRUE)),
                 tolower(normalizePath(source_snapshot, winslash = "/", mustWork = FALSE)))) {
    if (!file.copy(lefse_source, source_snapshot, overwrite = TRUE)) stop("Failed to snapshot current LEfSe output.")
  }
  annotation_path <- file.path(snapshot_dir, "LDA_Microviridae_annotation.csv")
  write.csv(annotation, annotation_path, row.names = FALSE, na = "")
  write.table(row, file.path(snapshot_dir, "lefse_annotation_source_row.tsv"),
              sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  # Keep portable provenance: the complete source table is copied alongside
  # the derived annotation, so verification does not require a personal path.
  provenance <- data.frame(source_file = "lefse_source_data.tsv",
                           source_md5 = unname(tools::md5sum(source_snapshot)),
                           comparison_id = row$Comparison_ID, virome = row$Virome,
                           condition = row$Condition, family = row$Family,
                           project = row$Analysis_Project, ocean_basin = row$Ocean_Basin,
                           coral_host_species = row$Coral_Host_Species,
                           sample_type = row$Sample_Type, display_threshold = 2.5,
                           role = "Post hoc node annotation; excluded from network construction")
  write.table(provenance, file.path(snapshot_dir, "lefse_annotation_provenance.tsv"),
              sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  annotation_path
}

run_network_analysis <- function(stage_args) {
  
  # Fixed manuscript-dataset reproduction module; not a general analysis API.
  
  # ==============================================================================
  # 0. PATHS / TARGET CONTEXT
  # ==============================================================================
  
  args <- stage_args
  usage <- paste(
    "Usage:",
    "  Rscript run_network.R full <input_dir> <output_package_dir>",
    "",
    "This module reproduces the fixed manuscript analysis from four network input files.",
    sep = "\n"
  )
  if (length(args) == 1L && args[[1]] %in% c("-h", "--help")) {
    cat(usage, "\n")
    quit(status = 0L)
  }
  if (length(args) != 2L) {
    stop(usage, call. = FALSE)
  }
  
  INPUT_DIR <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
  PACKAGE_DIR <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
  OUTPUT_DIR <- file.path(PACKAGE_DIR, "process_data")
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(igraph)
    library(ggraph)
    library(patchwork)
  })
  
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  TARGET_PROJECT <- "PRJNA925892"
  TARGET_OCEAN <- "Atlantic Ocean"
  TARGET_SPECIES <- "Orbicella faveolata"
  TARGET_SAMPLE_TYPE <- "Coral tissue"
  
  N_RESAMPLES <- 500L
  N_PERMUTATIONS <- 999L
  
  set.seed(20260813)
  
  cat("\n============================================================\n")
  cat("PROJECT-AWARE VIRAL ASSOCIATION NETWORK\n")
  cat("============================================================\n")
  cat("Project: ", TARGET_PROJECT, "\n", sep = "")
  cat("Sampling region: ", TARGET_OCEAN, "\n", sep = "")
  cat("Coral taxon: ", TARGET_SPECIES, "\n", sep = "")
  cat("Compartment: ", TARGET_SAMPLE_TYPE, "\n", sep = "")
  cat("Input: ", INPUT_DIR, "\n", sep = "")
  cat("Output: ", OUTPUT_DIR, "\n", sep = "")
  cat("============================================================\n\n")
  
  
  # ==============================================================================
  # 1. HELPERS
  # ==============================================================================
  
  log_msg <- function(...) {
    cat(
      format(Sys.time(), "%H:%M:%S"),
      paste0(...),
      "\n"
    )
  }
  
  
  fix_taxonomy_levels <- function(tax_df) {
  
    cols <- c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family"
    )
  
    for (col in cols) {
  
      if (!col %in% names(tax_df)) {
        tax_df[[col]] <- "Unclassified"
      }
    }
  
    tax_df %>%
      mutate(
        across(
          all_of(cols),
          ~ifelse(
            is.na(.) | . == "",
            "Unclassified",
            .
          )
        )
      ) %>%
      mutate(
        Family_Smart = if_else(
          !grepl(
            "Unclassified",
            Family,
            ignore.case = TRUE
          ),
          Family,
          NA_character_
        )
      )
  }
  
  
  clr_transform <- function(
      x,
      pseudo = 1e-6
  ) {
  
    x <- as.matrix(x)
  
    rel <- sweep(
      x,
      2,
      pmax(
        colSums(x),
        .Machine$double.eps
      ),
      "/"
    )
  
    lx <- log(
      rel + pseudo
    )
  
    sweep(
      lx,
      2,
      colMeans(lx),
      "-"
    )
  }
  
  
  bh_upper <- function(pmat) {
  
    p_adj <- matrix(
      1,
      nrow(pmat),
      ncol(pmat),
      dimnames = dimnames(pmat)
    )
  
    idx <- which(
      upper.tri(pmat)
    )
  
    p_adj[idx] <- p.adjust(
      pmat[idx],
      method = "BH"
    )
  
    p_adj[
      lower.tri(p_adj)
    ] <- t(p_adj)[
      lower.tri(p_adj)
    ]
  
    diag(p_adj) <- 0
  
    p_adj
  }
  
  
  correlation_matrices <- function(x) {
  
    n <- ncol(x)
  
    r <- suppressWarnings(
      cor(
        t(x),
        method = "spearman"
      )
    )
  
    r[
      is.na(r)
    ] <- 0
  
    p <- matrix(
      1,
      nrow(x),
      nrow(x),
      dimnames = list(
        rownames(x),
        rownames(x)
      )
    )
  
    if (nrow(x) > 1) {
  
      for (
        i in seq_len(
          nrow(x) - 1
        )
      ) {
  
        for (
          j in (i + 1):nrow(x)
        ) {
  
          tt <- suppressWarnings(
            r[i, j] *
              sqrt(
                (n - 2) /
                  pmax(
                    1 - r[i, j]^2,
                    .Machine$double.eps
                  )
              )
          )
  
          p[i, j] <- 2 *
            pt(
              -abs(tt),
              df = max(
                n - 2,
                1
              )
            )
        }
      }
    }
  
    p[
      lower.tri(p)
    ] <- t(p)[
      lower.tri(p)
    ]
  
    list(
      r = r,
      p = p,
      p_adj = bh_upper(p)
    )
  }
  
  
  make_network <- function(
      x,
      group
  ) {
  
    if (nrow(x) < 3) {
      stop(
        "At least three taxonomic groups are required."
      )
    }
  
    cm <- correlation_matrices(x)
  
    keep <- upper.tri(cm$r) &
      abs(cm$r) >= 0.7 &
      cm$p_adj <= 0.05
  
    ij <- which(
      keep,
      arr.ind = TRUE
    )
  
    if (nrow(ij) > 0) {
  
      edges <- tibble(
        From = rownames(cm$r)[ij[, 1]],
        To = colnames(cm$r)[ij[, 2]],
        Rho = cm$r[ij],
        Abs_Rho = abs(Rho),
        P_FDR = cm$p_adj[ij],
        Sign = if_else(
          Rho >= 0,
          "Positive",
          "Negative"
        ),
        Group = group
      ) %>%
        mutate(
          Key = map2_chr(
            From,
            To,
            ~paste(
              sort(
                c(.x, .y)
              ),
              collapse = "||"
            )
          )
        )
  
    } else {
  
      edges <- tibble(
        From = character(),
        To = character(),
        Rho = double(),
        Abs_Rho = double(),
        P_FDR = double(),
        Sign = character(),
        Group = character(),
        Key = character()
      )
    }
  
    # Same retained family set is used as the vertex denominator in both networks.
    g <- graph_from_data_frame(
      edges,
      directed = FALSE,
      vertices = tibble(
        name = rownames(x)
      )
    )
  
    E(g)$weight <- edges$Abs_Rho
    E(g)$sign <- edges$Sign
  
    comps <- components(g)
  
    comm <- if (ecount(g) > 0) {
      cluster_louvain(
        g,
        weights = E(g)$weight
      )
    } else {
      NULL
    }
  
    nodes <- tibble(
      Family = V(g)$name,
      Degree = degree(g),
      Strength = strength(
        g,
        weights = E(g)$weight
      ),
      Community = if (is.null(comm)) {
        NA_integer_
      } else {
        as.integer(
          membership(comm)
        )
      },
      Component = comps$membership,
      Group = group
    ) %>%
      arrange(
        desc(Degree),
        desc(Strength),
        Family
      )
  
    stats <- tibble(
      Group = group,
      Nodes = vcount(g),
      Active_Nodes = sum(
        degree(g) > 0
      ),
      Edges = ecount(g),
      Avg_Degree = mean(
        degree(g)
      ),
      Density = edge_density(g),
      Modularity = if (is.null(comm)) {
        NA_real_
      } else {
        modularity(comm)
      },
      Components = comps$no,
      Largest_Component_N = max(
        comps$csize
      ),
      Largest_Component_Fraction =
        max(comps$csize) /
        vcount(g),
      Positive_Edges = sum(
        edges$Sign == "Positive"
      ),
      Negative_Edges = sum(
        edges$Sign == "Negative"
      ),
      Positive_Edge_Ratio = if (
        ecount(g) > 0
      ) {
        mean(
          edges$Sign == "Positive"
        )
      } else {
        NA_real_
      }
    )
  
    list(
      graph = g,
      edges = edges,
      nodes = nodes,
      stats = stats,
      cor = cm
    )
  }
  
  
  # Draw each group to exactly the Healthy sample-count vector at every location.
  sample_balanced_by_site <- function(
      meta,
      group,
      site_col = "Location"
  ) {
  
    target <- meta %>%
      filter(
        Network_Group == "Healthy"
      ) %>%
      count(
        .data[[site_col]],
        name = "Target_N"
      )
  
    pool <- meta %>%
      filter(
        Network_Group == group
      )
  
    out <- map_dfr(
      seq_len(
        nrow(target)
      ),
      function(i) {
  
        site <- target[[site_col]][i]
        n_take <- target$Target_N[i]
  
        available <- pool %>%
          filter(
            .data[[site_col]] == site
          )
  
        if (nrow(available) == 0) {
          stop(
            "No ",
            group,
            " sample is available at required site: ",
            site
          )
        }
  
        available[
          sample(
            seq_len(
              nrow(available)
            ),
            n_take,
            replace = TRUE
          ),
          ,
          drop = FALSE
        ]
      }
    )
  
    out
  }
  
  
  topology_from_graph <- function(
      g,
      group,
      resample_id,
      sampled_n
  ) {
  
    if (
      vcount(g) < 2 ||
        ecount(g) < 1
    ) {
  
      return(
        tibble(
          Group = group,
          Resample = resample_id,
          Sample_N = sampled_n,
          Nodes = vcount(g),
          Active_Nodes = sum(
            degree(g) > 0
          ),
          Edges = ecount(g),
          Avg_Degree = 0,
          Density = 0,
          Modularity = NA_real_,
          Components = ifelse(
            vcount(g) > 0,
            components(g)$no,
            0
          ),
          Largest_Component_Fraction = ifelse(
            vcount(g) > 0,
            max(
              components(g)$csize
            ) /
              vcount(g),
            0
          ),
          Positive_Edge_Ratio = NA_real_
        )
      )
    }
  
    comm <- cluster_louvain(
      g,
      weights = E(g)$weight
    )
  
    comps <- components(g)
  
    tibble(
      Group = group,
      Resample = resample_id,
      Sample_N = sampled_n,
      Nodes = vcount(g),
      Active_Nodes = sum(
        degree(g) > 0
      ),
      Edges = ecount(g),
      Avg_Degree = mean(
        degree(g)
      ),
      Density = edge_density(g),
      Modularity = modularity(comm),
      Components = comps$no,
      Largest_Component_Fraction =
        max(comps$csize) /
        vcount(g),
      Positive_Edge_Ratio = mean(
        E(g)$sign == "Positive"
      )
    )
  }
  
  
  # ==============================================================================
  # 2. LOAD INPUT DATA
  # ==============================================================================
  
  log_msg("Loading data")
  
  required_files <- c(
    "24083_metagenome_TPM.tsv",
    "dna_virus_taxonomy.tsv",
    "sample_info_metagenome.tsv"
  )
  
  missing_files <- required_files[
    !file.exists(
      file.path(
        INPUT_DIR,
        required_files
      )
    )
  ]
  
  if (length(missing_files) > 0) {
  
    stop(
      paste0(
        "Missing input files in:\n",
        INPUT_DIR,
        "\n\n",
        paste(
          missing_files,
          collapse = "\n"
        )
      )
    )
  }
  
  
  tpm <- read.table(
    file.path(
      INPUT_DIR,
      "24083_metagenome_TPM.tsv"
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    row.names = 1
  )
  
  tax <- read.table(
    file.path(
      INPUT_DIR,
      "dna_virus_taxonomy.tsv"
    ),
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    check.names = FALSE
  )
  
  meta <- read.table(
    file.path(
      INPUT_DIR,
      "sample_info_metagenome.tsv"
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = ""
  )
  
  colnames(tax)[1] <- "Genome"
  
  tax_fixed <- fix_taxonomy_levels(
    tax
  )
  
  
  # ==============================================================================
  # 3. AGGREGATE TPM TO FORMALLY NAMED VIRAL FAMILIES
  # ==============================================================================
  
  tpm_fam <- tpm %>%
  
    rownames_to_column(
      "Genome"
    ) %>%
  
    inner_join(
      tax_fixed %>%
  
        select(
          Genome,
          Family_Smart
        ),
      by = "Genome"
    ) %>%
  
    filter(
      !is.na(
        Family_Smart
      )
    ) %>%
  
    group_by(
      Family_Smart
    ) %>%
  
    summarise(
      across(
        where(is.numeric),
        sum
      ),
      .groups = "drop"
    ) %>%
  
    column_to_rownames(
      "Family_Smart"
    )
  
  
  # ==============================================================================
  # 4. PROJECT-AWARE SAMPLE FILTERING
  # ==============================================================================
  
  required_meta <- c(
    "Sample",
    "Project",
    "Ocean_Basin",
    "Coral_Host_Species",
    "Sample_Type",
    "Location",
    "Status"
  )
  
  missing_meta <- setdiff(
    required_meta,
    names(meta)
  )
  
  if (length(missing_meta) > 0) {
  
    stop(
      paste0(
        "Metadata is missing:\n",
        paste(
          missing_meta,
          collapse = "\n"
        )
      )
    )
  }
  
  
  # First export all projects that match the broad ecological context.
  context_project_qc <- meta %>%
  
    filter(
      Ocean_Basin == TARGET_OCEAN,
      Coral_Host_Species == TARGET_SPECIES,
      Sample_Type == TARGET_SAMPLE_TYPE
    ) %>%
  
    mutate(
      Network_Group = case_when(
  
        grepl(
          "Mix|Recovery|Treated|Unknown|Unkonwn",
          Status,
          ignore.case = TRUE
        ) ~ NA_character_,
  
        grepl(
          "SCTLD",
          Status,
          ignore.case = TRUE
        ) ~ "SCTLD",
  
        grepl(
          "Healthy|Normal",
          Status,
          ignore.case = TRUE
        ) ~ "Healthy",
  
        TRUE ~ NA_character_
      )
    ) %>%
  
    filter(
      !is.na(
        Network_Group
      )
    ) %>%
  
    count(
      Project,
      Location,
      Network_Group,
      name = "N"
    ) %>%
  
    arrange(
      Project,
      Location,
      Network_Group
    )
  
  
  write.csv(
    context_project_qc,
    file.path(
      OUTPUT_DIR,
      "00_project_context_QC_before_target_filter.csv"
    ),
    row.names = FALSE
  )
  
  
  meta_sub <- meta %>%
  
    filter(
      Project == TARGET_PROJECT,
      Ocean_Basin == TARGET_OCEAN,
      Coral_Host_Species == TARGET_SPECIES,
      Sample_Type == TARGET_SAMPLE_TYPE
    ) %>%
  
    mutate(
  
      Network_Group = case_when(
  
        grepl(
          "Mix|Recovery|Treated|Unknown|Unkonwn",
          Status,
          ignore.case = TRUE
        ) ~ NA_character_,
  
        grepl(
          "SCTLD",
          Status,
          ignore.case = TRUE
        ) ~ "SCTLD",
  
        grepl(
          "Healthy|Normal",
          Status,
          ignore.case = TRUE
        ) ~ "Healthy",
  
        TRUE ~ NA_character_
      ),
  
      In_TPM =
        Sample %in%
        colnames(
          tpm_fam
        )
    ) %>%
  
    filter(
      !is.na(
        Network_Group
      ),
      In_TPM
    )
  
  
  if (
    n_distinct(
      meta_sub$Project
    ) != 1L ||
      unique(
        meta_sub$Project
      ) != TARGET_PROJECT
  ) {
  
    stop(
      "Project-aware QC failed."
    )
  }
  
  
  samps_h <- meta_sub %>%
  
    filter(
      Network_Group == "Healthy"
    ) %>%
  
    pull(
      Sample
    )
  
  
  samps_d <- meta_sub %>%
  
    filter(
      Network_Group == "SCTLD"
    ) %>%
  
    pull(
      Sample
    )
  
  
  if (
    length(samps_h) < 3 ||
      length(samps_d) < 3
  ) {
  
    stop(
      "Too few Healthy or SCTLD samples after Project-aware filtering."
    )
  }
  
  
  site_counts <- meta_sub %>%
  
    count(
      Project,
      Network_Group,
      Location,
      name = "N"
    ) %>%
  
    arrange(
      Location,
      Network_Group
    )
  
  
  write.csv(
    meta_sub,
    file.path(
      OUTPUT_DIR,
      "01_filtered_metadata.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    site_counts,
    file.path(
      OUTPUT_DIR,
      "05_site_group_counts.csv"
    ),
    row.names = FALSE
  )
  
  
  cat("\nProject-aware sample counts:\n")
  print(
    site_counts
  )
  
  cat(
    "\nHealthy total = ",
    length(samps_h),
    "; SCTLD total = ",
    length(samps_d),
    "\n",
    sep = ""
  )
  
  
  # ==============================================================================
  # 5. FEATURE FILTERING
  # ==============================================================================
  
  fam_means <- rowMeans(
    tpm_fam[
      ,
      c(
        samps_h,
        samps_d
      ),
      drop = FALSE
    ]
  )
  
  
  top_fams <- names(
    sort(
      fam_means,
      decreasing = TRUE
    )
  )[
    seq_len(
      min(
        100,
        length(fam_means)
      )
    )
  ]
  
  
  mat <- tpm_fam[
    top_fams,
    c(
      samps_h,
      samps_d
    ),
    drop = FALSE
  ]
  
  
  prev_h <- rowMeans(
    mat[
      ,
      samps_h,
      drop = FALSE
    ] > 0
  )
  
  prev_d <- rowMeans(
    mat[
      ,
      samps_d,
      drop = FALSE
    ] > 0
  )
  
  prev_pooled <- rowMeans(
    mat > 0
  )
  
  
  keep_prev <- prev_pooled >= 0.20
  
  mat_f <- mat[
    keep_prev,
    ,
    drop = FALSE
  ]
  
  mat_clr <- clr_transform(
    mat_f
  )
  
  
  # Preserve the exact feature-by-sample matrices that enter the network model.
  write.csv(
    data.frame(
      Family = rownames(mat_f),
      mat_f,
      check.names = FALSE
    ),
    file.path(
      OUTPUT_DIR,
      "02b_network_input_family_TPM_matrix.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    data.frame(
      Family = rownames(mat_clr),
      mat_clr,
      check.names = FALSE
    ),
    file.path(
      OUTPUT_DIR,
      "02c_network_input_family_CLR_matrix.csv"
    ),
    row.names = FALSE
  )
  
  
  prevalence_table <- tibble(
    Family = rownames(mat),
    Healthy_Prevalence = prev_h,
    SCTLD_Prevalence = prev_d,
    Pooled_Prevalence = prev_pooled,
    Retained_Pooled_Prevalence_Filter = keep_prev,
    Candidate_Top100 = TRUE
  )
  
  
  write.csv(
    prevalence_table,
    file.path(
      OUTPUT_DIR,
      "02_prevalence_filter_table.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    tibble(
      Family = top_fams,
      Combined_Mean_TPM =
        fam_means[top_fams]
    ),
    file.path(
      OUTPUT_DIR,
      "03_top100_candidate_groups.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    tibble(
      Group = c(
        "Healthy",
        "SCTLD"
      ),
      N = c(
        length(samps_h),
        length(samps_d)
      )
    ),
    file.path(
      OUTPUT_DIR,
      "04_group_sample_sizes.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 6. REFERENCE NETWORKS
  # ==============================================================================
  
  log_msg(
    "Building reference networks"
  )
  
  ref_h <- make_network(
    mat_clr[
      ,
      samps_h,
      drop = FALSE
    ],
    "Healthy"
  )
  
  ref_d <- make_network(
    mat_clr[
      ,
      samps_d,
      drop = FALSE
    ],
    "SCTLD"
  )
  
  
  write.csv(
    bind_rows(
      ref_h$stats,
      ref_d$stats
    ),
    file.path(
      OUTPUT_DIR,
      "06_reference_topology.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    bind_rows(
      ref_h$edges,
      ref_d$edges
    ),
    file.path(
      OUTPUT_DIR,
      "07_reference_edges.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    bind_rows(
      ref_h$nodes,
      ref_d$nodes
    ),
    file.path(
      OUTPUT_DIR,
      "08_reference_nodes.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 7. SITE-BALANCED BOOTSTRAP
  # ==============================================================================
  
  log_msg(
    "Running site-balanced resampling: ",
    N_RESAMPLES,
    " replicates"
  )
  
  
  topology_list <- vector(
    "list",
    N_RESAMPLES * 2
  )
  
  edge_list <- vector(
    "list",
    N_RESAMPLES * 2
  )
  
  design_list <- vector(
    "list",
    N_RESAMPLES * 2
  )
  
  
  for (
    b in seq_len(
      N_RESAMPLES
    )
  ) {
  
    h_bal_meta <- sample_balanced_by_site(
      meta_sub,
      "Healthy"
    )
  
    d_bal_meta <- sample_balanced_by_site(
      meta_sub,
      "SCTLD"
    )
  
    h_bal <- h_bal_meta$Sample
    d_bal <- d_bal_meta$Sample
  
    nh <- make_network(
      mat_clr[
        ,
        h_bal,
        drop = FALSE
      ],
      "Healthy"
    )
  
    nd <- make_network(
      mat_clr[
        ,
        d_bal,
        drop = FALSE
      ],
      "SCTLD"
    )
  
    topology_list[[2 * b - 1]] <-
      topology_from_graph(
        nh$graph,
        "Healthy",
        b,
        length(h_bal)
      )
  
    topology_list[[2 * b]] <-
      topology_from_graph(
        nd$graph,
        "SCTLD",
        b,
        length(d_bal)
      )
  
    edge_list[[2 * b - 1]] <-
      nh$edges %>%
  
        select(
          Key,
          From,
          To,
          Rho,
          Abs_Rho,
          Sign
        ) %>%
  
        mutate(
          Group = "Healthy",
          Resample = b
        )
  
    edge_list[[2 * b]] <-
      nd$edges %>%
  
        select(
          Key,
          From,
          To,
          Rho,
          Abs_Rho,
          Sign
        ) %>%
  
        mutate(
          Group = "SCTLD",
          Resample = b
        )
  
    design_list[[2 * b - 1]] <-
      h_bal_meta %>%
  
        mutate(
          Group = "Healthy",
          Resample = b,
          Draw_Index = row_number()
        ) %>%
  
        select(
          Resample,
          Group,
          Location,
          Draw_Index,
          Sample
        )
  
    design_list[[2 * b]] <-
      d_bal_meta %>%
  
        mutate(
          Group = "SCTLD",
          Resample = b,
          Draw_Index = row_number()
        ) %>%
  
        select(
          Resample,
          Group,
          Location,
          Draw_Index,
          Sample
        )
  }
  
  
  topology_resamples <- bind_rows(
    topology_list
  )
  
  edges_resamples <- bind_rows(
    edge_list
  )
  
  resampling_design <- bind_rows(
    design_list
  )
  
  
  write.csv(
    edges_resamples,
    file.path(
      OUTPUT_DIR,
      "09c_site_balanced_edge_resamples.csv"
    ),
    row.names = FALSE
  )
  
  
  edge_stability <- edges_resamples %>%
  
    group_by(
      Group,
      Key,
      From,
      To,
      Sign
    ) %>%
  
    summarise(
      Retained_Replicates = n(),
      Mean_Abs_Rho = mean(
        Abs_Rho
      ),
      Median_Abs_Rho = median(
        Abs_Rho
      ),
      .groups = "drop"
    ) %>%
  
    mutate(
      Resample_N = N_RESAMPLES,
      Retention =
        Retained_Replicates /
        Resample_N
    ) %>%
  
    arrange(
      Group,
      desc(Retention),
      desc(Mean_Abs_Rho)
    )
  
  
  topology_summary <- topology_resamples %>%
  
    group_by(
      Group
    ) %>%
  
    summarise(
      across(
        c(
          Nodes,
          Active_Nodes,
          Edges,
          Avg_Degree,
          Density,
          Modularity,
          Components,
          Largest_Component_Fraction,
          Positive_Edge_Ratio
        ),
        list(
          Median = ~median(
            .x,
            na.rm = TRUE
          ),
          Q025 = ~quantile(
            .x,
            .025,
            na.rm = TRUE
          ),
          Q975 = ~quantile(
            .x,
            .975,
            na.rm = TRUE
          )
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  
  write.csv(
    topology_resamples,
    file.path(
      OUTPUT_DIR,
      "09_site_balanced_topology_resamples.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    resampling_design,
    file.path(
      OUTPUT_DIR,
      "09b_site_balanced_resampling_design.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    topology_summary,
    file.path(
      OUTPUT_DIR,
      "10_site_balanced_topology_summary.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    edge_stability,
    file.path(
      OUTPUT_DIR,
      "11_edge_stability_all.csv"
    ),
    row.names = FALSE
  )
  
  
  stable_edges <- edge_stability %>%
  
    filter(
      Retention >= 0.50
    ) %>%
  
    mutate(
      Threshold =
        "Retention >= 0.50"
    )
  
  
  write.csv(
    stable_edges,
    file.path(
      OUTPUT_DIR,
      "12_edges_retained_in_at_least_50pct.csv"
    ),
    row.names = FALSE
  )
  
  
  # The common stable core is a display layer only. Inferential topology remains
  # based on the prespecified full family set selected before network comparison.
  stable_nodes_by_group <- bind_rows(
    stable_edges %>% transmute(Group, Family = From),
    stable_edges %>% transmute(Group, Family = To)
  ) %>%
    distinct(Group, Family)
  
  common_stable_families <- intersect(
    stable_nodes_by_group %>% filter(Group == "Healthy") %>% pull(Family),
    stable_nodes_by_group %>% filter(Group == "SCTLD") %>% pull(Family)
  ) %>%
    sort()
  
  if (length(common_stable_families) != 8L ||
      !"Microviridae" %in% common_stable_families) {
    stop("Expected an eight-family common stable core containing Microviridae.")
  }
  
  write.csv(
    tibble(
      Family = common_stable_families,
      Display_Set = "Stable-edge endpoint in both Healthy and SCTLD",
      Used_For_Topology_Inference = FALSE
    ),
    file.path(OUTPUT_DIR, "12a_common_stable_core_families.csv"),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 8. PAIRED BOOTSTRAP DELTAS
  # ==============================================================================
  
  metric_tests <- topology_resamples %>%
  
    select(
      Resample,
      Group,
      Active_Nodes,
      Edges,
      Density,
      Avg_Degree,
      Modularity,
      Components,
      Largest_Component_Fraction,
      Positive_Edge_Ratio
    ) %>%
  
    pivot_wider(
      names_from = Group,
      values_from = -Resample
    ) %>%
  
    mutate(
      Delta_SCTLD_minus_Healthy_Active_Nodes =
        Active_Nodes_SCTLD -
        Active_Nodes_Healthy,
  
      Delta_SCTLD_minus_Healthy_Edges =
        Edges_SCTLD -
        Edges_Healthy,
  
      Delta_SCTLD_minus_Healthy_Density =
        Density_SCTLD -
        Density_Healthy,
  
      Delta_SCTLD_minus_Healthy_Avg_Degree =
        Avg_Degree_SCTLD -
        Avg_Degree_Healthy,
  
      Delta_SCTLD_minus_Healthy_Modularity =
        Modularity_SCTLD -
        Modularity_Healthy,
  
      Delta_SCTLD_minus_Healthy_Components =
        Components_SCTLD -
        Components_Healthy,
  
      Delta_SCTLD_minus_Healthy_Largest_Component =
        Largest_Component_Fraction_SCTLD -
        Largest_Component_Fraction_Healthy,
  
      Delta_SCTLD_minus_Healthy_Positive_Ratio =
        Positive_Edge_Ratio_SCTLD -
        Positive_Edge_Ratio_Healthy
    )
  
  
  write.csv(
    metric_tests,
    file.path(
      OUTPUT_DIR,
      "13_paired_resample_deltas.csv"
    ),
    row.names = FALSE
  )
  
  
  delta_summary <- metric_tests %>%
  
    summarise(
      across(
        starts_with(
          "Delta_"
        ),
        list(
          Median = ~median(
            .x,
            na.rm = TRUE
          ),
          Q025 = ~quantile(
            .x,
            .025,
            na.rm = TRUE
          ),
          Q975 = ~quantile(
            .x,
            .975,
            na.rm = TRUE
          ),
          Proportion_Delta_Negative =
            ~mean(
              .x < 0,
              na.rm = TRUE
            )
        ),
        .names = "{.col}_{.fn}"
      )
    )
  
  
  write.csv(
    delta_summary,
    file.path(
      OUTPUT_DIR,
      "14_paired_resample_delta_summary.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 9. WITHIN-LOCATION LABEL PERMUTATION
  # ==============================================================================
  
  log_msg(
    "Running ",
    N_PERMUTATIONS,
    " within-location label permutations"
  )
  
  
  perm_topology_list <- vector(
    "list",
    N_PERMUTATIONS * 2
  )
  
  
  perm_design_list <- vector(
    "list",
    N_PERMUTATIONS
  )
  
  
  for (
    b in seq_len(
      N_PERMUTATIONS
    )
  ) {
  
    perm_meta <- meta_sub %>%
  
      group_by(
        Location
      ) %>%
  
      mutate(
        Perm_Group =
          sample(
            Network_Group,
            size = n(),
            replace = FALSE
          )
      ) %>%
  
      ungroup()
  
  
    perm_design_list[[b]] <- perm_meta %>%
  
      transmute(
        Permutation = b,
        Sample,
        Location,
        Observed_Group = Network_Group,
        Permuted_Group = Perm_Group
      )
  
  
    h_perm <- perm_meta %>%
  
      filter(
        Perm_Group == "Healthy"
      ) %>%
  
      pull(
        Sample
      )
  
  
    d_perm <- perm_meta %>%
  
      filter(
        Perm_Group == "SCTLD"
      ) %>%
  
      pull(
        Sample
      )
  
  
    nh <- make_network(
      mat_clr[
        ,
        h_perm,
        drop = FALSE
      ],
      "Healthy"
    )
  
  
    nd <- make_network(
      mat_clr[
        ,
        d_perm,
        drop = FALSE
      ],
      "SCTLD"
    )
  
  
    perm_topology_list[[2 * b - 1]] <-
      topology_from_graph(
        nh$graph,
        "Healthy",
        b,
        length(h_perm)
      )
  
  
    perm_topology_list[[2 * b]] <-
      topology_from_graph(
        nd$graph,
        "SCTLD",
        b,
        length(d_perm)
      )
  }
  
  
  perm_topology <- bind_rows(
    perm_topology_list
  )
  
  
  perm_design <- bind_rows(
    perm_design_list
  )
  
  
  perm_delta <- perm_topology %>%
  
    select(
      Resample,
      Group,
      Active_Nodes,
      Edges,
      Density,
      Avg_Degree,
      Modularity,
      Components,
      Largest_Component_Fraction,
      Positive_Edge_Ratio
    ) %>%
  
    pivot_wider(
      names_from = Group,
      values_from = -Resample
    ) %>%
  
    transmute(
      Permutation = Resample,
  
      Delta_Active_Nodes =
        Active_Nodes_SCTLD -
        Active_Nodes_Healthy,
  
      Delta_Edges =
        Edges_SCTLD -
        Edges_Healthy,
  
      Delta_Density =
        Delta_Edges /
        choose(
          nrow(mat_clr),
          2
        ),
  
      Delta_Avg_Degree =
        2 *
        Delta_Edges /
        nrow(mat_clr),
  
      Delta_Modularity =
        Modularity_SCTLD -
        Modularity_Healthy,
  
      Delta_Components =
        Components_SCTLD -
        Components_Healthy,
  
      Delta_Largest_Component_Fraction =
        Largest_Component_Fraction_SCTLD -
        Largest_Component_Fraction_Healthy,
  
      Delta_Positive_Edge_Ratio =
        Positive_Edge_Ratio_SCTLD -
        Positive_Edge_Ratio_Healthy
    )
  
  
  observed_delta <- tibble(
  
    Metric = c(
      "Active_Nodes",
      "Edges",
      "Density",
      "Avg_Degree",
      "Modularity",
      "Components",
      "Largest_Component_Fraction",
      "Positive_Edge_Ratio"
    ),
  
    Observed_Delta_SCTLD_minus_Healthy = c(
      ref_d$stats$Active_Nodes -
        ref_h$stats$Active_Nodes,
  
      ref_d$stats$Edges -
        ref_h$stats$Edges,
  
      (
        ref_d$stats$Edges -
          ref_h$stats$Edges
      ) /
        choose(
          nrow(mat_clr),
          2
        ),
  
      2 *
        (
          ref_d$stats$Edges -
            ref_h$stats$Edges
        ) /
        nrow(mat_clr),
  
      ref_d$stats$Modularity -
        ref_h$stats$Modularity,
  
      ref_d$stats$Components -
        ref_h$stats$Components,
  
      ref_d$stats$Largest_Component_Fraction -
        ref_h$stats$Largest_Component_Fraction,
  
      ref_d$stats$Positive_Edge_Ratio -
        ref_h$stats$Positive_Edge_Ratio
    ),
  
    Permutation_Column = c(
      "Delta_Active_Nodes",
      "Delta_Edges",
      "Delta_Density",
      "Delta_Avg_Degree",
      "Delta_Modularity",
      "Delta_Components",
      "Delta_Largest_Component_Fraction",
      "Delta_Positive_Edge_Ratio"
    )
  )
  
  
  perm_summary <- observed_delta %>%
  
    rowwise() %>%
  
    mutate(
  
      Valid_Permutations =
        sum(
          is.finite(
            perm_delta[[as.character(Permutation_Column)]]
          )
        ),
  
      Null_Median =
        median(
          perm_delta[[as.character(Permutation_Column)]],
          na.rm = TRUE
        ),
  
      Null_Q025 =
        quantile(
          perm_delta[[as.character(Permutation_Column)]],
          0.025,
          na.rm = TRUE
        ),
  
      Null_Q975 =
        quantile(
          perm_delta[[as.character(Permutation_Column)]],
          0.975,
          na.rm = TRUE
        ),
  
      P_Lower =
        (
          1 +
            sum(
              perm_delta[[as.character(Permutation_Column)]] <=
                Observed_Delta_SCTLD_minus_Healthy,
              na.rm = TRUE
            )
        ) /
        (
          Valid_Permutations +
            1
        ),
  
      P_Upper =
        (
          1 +
            sum(
              perm_delta[[as.character(Permutation_Column)]] >=
                Observed_Delta_SCTLD_minus_Healthy,
              na.rm = TRUE
            )
        ) /
        (
          Valid_Permutations +
            1
        ),
  
      P_Two_Sided =
        min(
          1,
          2 *
            min(
              P_Lower,
              P_Upper
            )
        )
    ) %>%
  
    ungroup() %>%
  
    mutate(
      P_BH_Across_8_Topology_Metrics =
        p.adjust(
          P_Two_Sided,
          method = "BH"
        )
    )
  
  
  write.csv(
    perm_topology,
    file.path(
      OUTPUT_DIR,
      "18_within_location_permutation_topology.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    perm_design,
    file.path(
      OUTPUT_DIR,
      "18a_within_location_permutation_design.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    perm_delta,
    file.path(
      OUTPUT_DIR,
      "19_within_location_permutation_deltas.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    perm_summary,
    file.path(
      OUTPUT_DIR,
      "20_within_location_permutation_summary.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 10. PREVALENCE SENSITIVITY
  # ==============================================================================
  
  sensitivity <- map_dfr(
    c(
      0.10,
      0.20,
      0.30
    ),
    function(prev_cut) {
  
      keep <- prev_pooled >= prev_cut
  
      x <- clr_transform(
        mat[
          keep,
          ,
          drop = FALSE
        ]
      )
  
      h <- make_network(
        x[
          ,
          samps_h,
          drop = FALSE
        ],
        "Healthy"
      )
  
      d <- make_network(
        x[
          ,
          samps_d,
          drop = FALSE
        ],
        "SCTLD"
      )
  
      bind_rows(
        h$stats,
        d$stats
      ) %>%
  
        mutate(
          Prevalence_Cutoff =
            prev_cut,
  
          Features_Retained =
            sum(keep)
        )
    }
  )
  
  
  write.csv(
    sensitivity,
    file.path(
      OUTPUT_DIR,
      "21_prevalence_threshold_sensitivity.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==============================================================================
  # 11. AUDIT NOTE
  # ==============================================================================
  
  healthy_site_vector <- meta_sub %>%
  
    filter(
      Network_Group == "Healthy"
    ) %>%
  
    count(
      Location
    ) %>%
  
    arrange(
      Location
    )
  
  
  writeLines(
    c(
      "Project-aware robust viral family association network",
      "",
      paste0(
        "Project: ",
        TARGET_PROJECT
      ),
      paste0(
        "Sampling region: ",
        TARGET_OCEAN
      ),
      paste0(
        "Coral taxon: ",
        TARGET_SPECIES
      ),
      paste0(
        "Compartment: ",
        TARGET_SAMPLE_TYPE
      ),
      "",
      paste0(
        "Healthy n = ",
        length(samps_h)
      ),
      paste0(
        "SCTLD n = ",
        length(samps_d)
      ),
      paste0(
        "Locations = ",
        paste(
          healthy_site_vector$Location,
          collapse = "; "
        )
      ),
      paste0(
        "Healthy location-count vector = ",
        paste(
          healthy_site_vector$n,
          collapse = ","
        )
      ),
      "",
      "Preprocessing: formally named viral families only; top 100 by combined mean TPM; pooled prevalence >=20%; CLR transformation.",
      "Network: Spearman association; |rho| >= 0.7 and BH-adjusted P <= 0.05 across tested family pairs.",
      paste0(
        "Bootstrap: ",
        N_RESAMPLES,
        " site-balanced resamples; both groups sampled to the Healthy location-count vector."
      ),
      "Stable edge: bootstrap retention >= 0.50.",
      paste0(
        "Permutation: ",
        N_PERMUTATIONS,
        " health-label permutations within Location; observed location-specific group counts preserved."
      ),
      "Topology permutation P values are BH-adjusted across eight prespecified topology metrics, including Active_Nodes (connected families).",
      "Prevalence sensitivity: 10%, 20%, 30%.",
      "",
      "IMPORTANT: all samples entering the network are explicitly restricted to one Project."
    ),
    file.path(
      OUTPUT_DIR,
      "17_analysis_audit_notes.txt"
    )
  )
  
  
  writeLines(
    capture.output(
      sessionInfo()
    ),
    file.path(
      OUTPUT_DIR,
      "R_sessionInfo.txt"
    )
  )
  
  
  log_msg(
    "Complete. Output: ",
    OUTPUT_DIR
  )
}

run_network_figure <- function(stage_args) {
  
  # Fixed manuscript-dataset figure module with intentional result assertions.
  
  args <- stage_args
  usage <- paste(
    "Usage:",
    "  Rscript run_network.R plot-only <analysis_package_dir> <output_dir> <lefse_source_data.tsv>",
    "",
    "The analysis package must contain input_data/ and process_data/.",
    sep = "\n"
  )
  if (length(args) == 1L && args[[1]] %in% c("-h", "--help")) {
    cat(usage, "\n")
    quit(status = 0L)
  }
  if (length(args) != 3L) {
    stop(usage, call. = FALSE)
  }
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(ggrepel)
  })
  
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Package 'ragg' is required for high-resolution raster export.")
  }
  
  source_package_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  source_dir <- file.path(source_package_dir, "process_data")
  input_dir <- file.path(source_package_dir, "input_data")
  figure_dir <- file.path(output_dir, "figures")
  figure_source_dir <- file.path(output_dir, "figure_source_data")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)
  
  required_files <- c(
    "01_filtered_metadata.csv",
    "05_site_group_counts.csv",
    "06_reference_topology.csv",
    "08_reference_nodes.csv",
    "09_site_balanced_topology_resamples.csv",
    "10_site_balanced_topology_summary.csv",
    "11_edge_stability_all.csv",
    "12_edges_retained_in_at_least_50pct.csv",
    "12a_common_stable_core_families.csv",
    "20_within_location_permutation_summary.csv"
  )
  missing_files <- required_files[!file.exists(file.path(source_dir, required_files))]
  if (length(missing_files) > 0) {
    stop("Missing prepared source files: ", paste(missing_files, collapse = ", "))
  }
  
  healthy_col <- "#2878A8"
  sctld_col <- "#D0584B"
  candidate_fill <- "#F2C14E"
  other_fill <- "#FFFFFF"
  inactive_col <- "#BCC3C7"
  text_col <- "#202427"
  subtle_col <- "#657078"
  grid_col <- "#D6DBDE"
  
  stable_edges <- read.csv(
    file.path(source_dir, "12_edges_retained_in_at_least_50pct.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Retention = as.numeric(Retention),
      Mean_Abs_Rho = as.numeric(Mean_Abs_Rho),
      Median_Abs_Rho = as.numeric(Median_Abs_Rho)
    )
  
  edge_stability_all <- read.csv(
    file.path(source_dir, "11_edge_stability_all.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Retention = as.numeric(Retention),
      Mean_Abs_Rho = as.numeric(Mean_Abs_Rho),
      Median_Abs_Rho = as.numeric(Median_Abs_Rho)
    )
  
  topology <- read.csv(
    file.path(source_dir, "09_site_balanced_topology_resamples.csv"),
    check.names = FALSE
  )
  topology_summary_full <- read.csv(
    file.path(source_dir, "10_site_balanced_topology_summary.csv"),
    check.names = FALSE
  )
  permutation_summary <- read.csv(
    file.path(source_dir, "20_within_location_permutation_summary.csv"),
    check.names = FALSE
  )
  site_counts <- read.csv(
    file.path(source_dir, "05_site_group_counts.csv"),
    check.names = FALSE
  )
  reference_topology <- read.csv(
    file.path(source_dir, "06_reference_topology.csv"),
    check.names = FALSE
  )
  reference_nodes <- read.csv(
    file.path(source_dir, "08_reference_nodes.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata <- read.csv(
    file.path(source_dir, "01_filtered_metadata.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  lda_candidates_tbl <- read.csv(
    args[[3]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # The LEfSe row and network must use exactly the same samples in each health
  # group. Context labels and sample counts alone cannot establish this.
  lda_source_row <- read.delim(file.path(dirname(args[[3]]), "lefse_annotation_source_row.tsv"),
                               check.names = FALSE, stringsAsFactors = FALSE)
  for (health_group in c("Healthy", "SCTLD")) {
    membership_field <- if (health_group == "Healthy") "Healthy_Samples" else "Condition_Samples"
    lefse_samples <- strsplit(lda_source_row[[membership_field]][[1]], ";", fixed = TRUE)[[1]]
    network_samples <- as.character(metadata$Sample[metadata$Network_Group == health_group])
    if (anyDuplicated(lefse_samples) || !setequal(lefse_samples, network_samples))
      stop("LEfSe and network sample membership differs for ", health_group)
  }
  lda_candidates <- unique(lda_candidates_tbl$LDA_Family)
  union_layout_path <- file.path(input_dir, "network_union_layout.csv")
  if (!file.exists(union_layout_path)) {
    stop("Missing fixed 11-family union layout: ", union_layout_path)
  }
  
  EXPECTED_RESULTS <- list(
    analysis_families = 18L,
    union_edges = 17L,
    group_stable_edge_records = 25L,
    topology_tests = 8L,
    display_families = 11L,
    omitted_families = 7L,
    lda_candidate = "Microviridae"
  )
  
  # These are manuscript-dataset reproduction assertions, not tunable filters.
  if (!identical(lda_candidates, EXPECTED_RESULTS$lda_candidate)) {
    stop("Only Microviridae may be annotated as the named LDA candidate.")
  }
  if (n_distinct(reference_nodes$Family) != EXPECTED_RESULTS$analysis_families) {
    stop("The underlying analysis must retain all 18 family-level network nodes.")
  }
  if (n_distinct(stable_edges$Key) != EXPECTED_RESULTS$union_edges ||
      nrow(stable_edges) != EXPECTED_RESULTS$group_stable_edge_records) {
    stop("Expected 17 union edges and 25 group-specific stable edges.")
  }
  if (nrow(permutation_summary) != EXPECTED_RESULTS$topology_tests ||
      !"Active_Nodes" %in% permutation_summary$Metric) {
    stop("The corrected permutation table must contain eight topology metrics.")
  }
  
  stable_keys <- stable_edges %>% distinct(Key)
  edge_stability_union <- edge_stability_all %>%
    semi_join(stable_keys, by = "Key") %>%
    select(
      Group, Key, From, To, Sign, Retained_Replicates, Resample_N,
      Retention, Mean_Abs_Rho, Median_Abs_Rho
    )
  
  if (nrow(edge_stability_union) != 2L * EXPECTED_RESULTS$union_edges) {
    stop("Expected two health-state retention values for each of 17 union edges.")
  }
  
  network_families <- sort(unique(reference_nodes$Family))
  node_layout <- read.csv(
    union_layout_path, check.names = FALSE, stringsAsFactors = FALSE
  ) %>%
    select(Family, x, y) %>%
    mutate(x = as.numeric(x), y = as.numeric(y))
  display_families <- sort(unique(c(stable_edges$From, stable_edges$To)))
  hidden_families <- setdiff(network_families, display_families)
  if (length(display_families) != EXPECTED_RESULTS$display_families ||
      length(hidden_families) != EXPECTED_RESULTS$omitted_families ||
      !setequal(node_layout$Family, display_families) ||
      !EXPECTED_RESULTS$lda_candidate %in% display_families) {
    stop("Expected the fixed 11-family stable-edge union containing Microviridae.")
  }
  
  family_degree_audit <- expand_grid(
    Family = network_families, Group = c("Healthy", "SCTLD")
  ) %>%
    left_join(
      bind_rows(
        stable_edges %>% transmute(Group, Family = From),
        stable_edges %>% transmute(Group, Family = To)
      ) %>% count(Group, Family, name = "Stable_Degree"),
      by = c("Group", "Family")
    ) %>%
    mutate(
      Stable_Degree = replace_na(Stable_Degree, 0L),
      Shown_In_Network_Panels = Family %in% display_families,
      Display_Reason = if_else(
        Shown_In_Network_Panels,
        "Stable-edge endpoint in at least one health state; 11-family union display",
        "No stable edge in either group; excluded from a/b display only"
      )
    )
  write.csv(
    family_degree_audit,
    file.path(figure_source_dir, "panel_a_b_full18_display_audit.csv"),
    row.names = FALSE
  )
  
  hidden_family_summary <- family_degree_audit %>%
    filter(!Shown_In_Network_Panels) %>%
    select(Family, Group, Stable_Degree) %>%
    pivot_wider(names_from = Group, values_from = Stable_Degree,
                names_prefix = "Stable_Degree_") %>%
    left_join(
      edge_stability_all %>%
        pivot_longer(c(From, To), values_to = "Family") %>%
        group_by(Family, Group) %>%
        summarise(Max_Edge_Retention = max(Retention), .groups = "drop") %>%
        pivot_wider(names_from = Group, values_from = Max_Edge_Retention,
                    names_prefix = "Max_Retention_"),
      by = "Family"
    ) %>%
    mutate(Display_Exclusion_Only = TRUE)
  prevalence_path <- file.path(source_dir, "02_prevalence_filter_table.csv")
  if (file.exists(prevalence_path)) {
    hidden_family_summary <- hidden_family_summary %>%
      left_join(read.csv(prevalence_path, check.names = FALSE), by = "Family")
  }
  write.csv(
    hidden_family_summary,
    file.path(figure_source_dir, "panel_a_b_excluded7_evidence.csv"),
    row.names = FALSE
  )
  
  group_panel_data <- function(group_name) {
    edges <- stable_edges %>%
      filter(
        Group == group_name,
        From %in% display_families,
        To %in% display_families
      ) %>%
      left_join(node_layout %>% select(From = Family, x, y), by = "From") %>%
      rename(x_from = x, y_from = y) %>%
      left_join(node_layout %>% select(To = Family, x, y), by = "To") %>%
      rename(x_to = x, y_to = y)
  
    degree_tbl <- bind_rows(
      edges %>% transmute(Family = From),
      edges %>% transmute(Family = To)
    ) %>%
      count(Family, name = "Stable_Degree")
  
    nodes <- node_layout %>%
      left_join(degree_tbl, by = "Family") %>%
      mutate(
        Stable_Degree = replace_na(Stable_Degree, 0L),
        Activity = if_else(Stable_Degree > 0, "Connected", "No stable edge"),
        LDA_Status = if_else(
          Family %in% lda_candidates,
          "SCTLD-enriched LDA candidate",
          "Other family"
        ),
        Node_Class = case_when(
          Family %in% lda_candidates ~ "Microviridae LDA overlay",
          Stable_Degree == 0 ~ "No stable edge",
          TRUE ~ "Connected family"
        ),
        Border = if_else(
          Activity == "Connected",
          if_else(group_name == "Healthy", healthy_col, sctld_col),
          inactive_col
        ),
        Label_Colour = if_else(Activity == "Connected", text_col, "#899298"),
        Node_Alpha = if_else(Activity == "Connected", 1, 0.50)
      )
  
    list(edges = edges, nodes = nodes)
  }
  
  theme_network <- theme_void(base_size = 7, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9.2, face = "bold", colour = text_col, hjust = 0),
      plot.subtitle = element_text(
        size = 6.4, colour = subtle_col, hjust = 0, margin = margin(b = 2.5)
      ),
      plot.margin = margin(3, 4, 3, 4),
      panel.border = element_rect(colour = grid_col, fill = NA, linewidth = 0.35)
    )
  
  network_panel <- function(group_name) {
    dat <- group_panel_data(group_name)
    group_col <- if (group_name == "Healthy") healthy_col else sctld_col
    group_title <- if (group_name == "Healthy") {
      "a  Healthy / Normal (Exposed)"
    } else {
      "b  SCTLD"
    }
    group_subtitle <- paste0(
      "11-family stable-edge union\n", nrow(dat$edges), " edges | ",
      sum(dat$nodes$Stable_Degree > 0), " connected"
    )
  
    ggplot() +
      geom_segment(
        data = dat$edges,
        aes(
          x = x_from, y = y_from, xend = x_to, yend = y_to,
          linewidth = Retention
        ),
        colour = group_col, alpha = 0.80, lineend = "round"
      ) +
      geom_point(
        data = dat$nodes,
        aes(x = x, y = y, size = Stable_Degree, fill = Node_Class),
        shape = 21, colour = dat$nodes$Border, alpha = dat$nodes$Node_Alpha,
        stroke = 0.72
      ) +
      ggrepel::geom_text_repel(
        data = dat$nodes,
        aes(x = x, y = y, label = Family),
        colour = dat$nodes$Label_Colour, family = "Arial", size = 1.95,
        lineheight = 0.92, seed = 20260909, max.overlaps = Inf,
        box.padding = 0.24, point.padding = 0.28, force = 0.75,
        min.segment.length = 0, segment.colour = "#AAB1B5",
        segment.size = 0.28
      ) +
      scale_fill_manual(
        values = c(
          "Microviridae LDA overlay" = candidate_fill,
          "Connected family" = other_fill,
          "No stable edge" = "#E4E7E9"
        ),
        breaks = "Microviridae LDA overlay",
        labels = "Microviridae (LDA |score| >= 2.5)",
        name = NULL
      ) +
      scale_size_continuous(
        range = c(2.0, 5.4), breaks = c(0, 2, 4), limits = c(0, 4),
        name = "Stable degree"
      ) +
      scale_linewidth_continuous(
        range = c(0.55, 2.0), breaks = c(0.50, 0.75, 1.00),
        limits = c(0.50, 1.00), name = "Bootstrap retention"
      ) +
      guides(
        linewidth = guide_legend(order = 1, override.aes = list(colour = "#687078")),
        size = guide_legend(order = 2),
        fill = guide_legend(order = 3, override.aes = list(size = 3.6, alpha = 1))
      ) +
      coord_equal(xlim = c(-2.75, 2.75), ylim = c(-2.55, 2.55), clip = "off") +
      labs(title = group_title, subtitle = group_subtitle) +
      theme_network
  }
  
  p_healthy <- network_panel("Healthy")
  p_sctld <- network_panel("SCTLD")
  
  network_figure <- p_healthy + p_sctld +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.title = element_text(size = 6.0, face = "bold"),
      legend.text = element_text(size = 5.65),
      legend.key.height = unit(3.2, "mm")
    )
  
  edge_stability_plot_data <- edge_stability_union %>%
    mutate(
      Edge_Label = str_replace(Key, "\\|\\|", " - "),
      Group_Label = factor(
        recode(Group, Healthy = "Healthy / Normal", SCTLD = "SCTLD"),
        levels = c("Healthy / Normal", "SCTLD")
      ),
      Microviridae_Linked = str_detect(Key, "(^|\\|\\|)Microviridae(\\|\\||$)"),
      Retention_Percent = sprintf("%.1f%%", 100 * Retention),
      Stable_At_50pct = Retention >= 0.50
    )
  
  edge_stability_wide <- edge_stability_plot_data %>%
    select(Key, Edge_Label, Group_Label, Retention, Microviridae_Linked) %>%
    pivot_wider(names_from = Group_Label, values_from = Retention) %>%
    rename(
      Healthy_Retention = `Healthy / Normal`,
      SCTLD_Retention = SCTLD
    ) %>%
    mutate(
      Delta_Retention = SCTLD_Retention - Healthy_Retention,
      Delta_Percentage_Points = 100 * Delta_Retention,
      Change_Direction = case_when(
        Delta_Retention < 0 ~ "Lower in SCTLD",
        Delta_Retention > 0 ~ "Higher in SCTLD",
        TRUE ~ "No change"
      )
    ) %>%
    arrange(Delta_Retention, desc(Healthy_Retention), Edge_Label)
  
  # The first row in the ordered table should appear at the top of the plot.
  edge_levels <- rev(edge_stability_wide$Edge_Label)
  edge_stability_wide <- edge_stability_wide %>%
    mutate(Edge_Label = factor(Edge_Label, levels = edge_levels))
  edge_stability_plot_data <- edge_stability_plot_data %>%
    mutate(Edge_Label = factor(Edge_Label, levels = edge_levels))
  
  p_retention_heatmap <- ggplot(
    edge_stability_plot_data,
    aes(x = Group_Label, y = Edge_Label, fill = Retention)
  ) +
    geom_tile(width = 0.88, height = 0.79, colour = "white", linewidth = 0.65) +
    geom_text(
      aes(label = Retention_Percent), family = "Arial", fontface = "bold",
      size = 1.85, colour = text_col
    ) +
    geom_point(
      data = edge_stability_plot_data %>% filter(Stable_At_50pct),
      position = position_nudge(x = 0.31, y = 0.25),
      shape = 16, size = 0.65, colour = text_col
    ) +
    scale_fill_gradient2(
      low = "#F4F6F7", mid = "#B9CCD9", high = "#2D678F",
      midpoint = 0.50, limits = c(0, 1), breaks = c(0, 0.5, 1),
      labels = scales::label_percent(accuracy = 0.1), name = "Bootstrap retention",
      guide = guide_colourbar(display = "rectangles", nbin = 128)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 6.0, face = "bold", colour = text_col),
      axis.text.y = element_text(size = 5.55, colour = text_col, margin = margin(r = 3)),
      legend.position = "bottom", legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.5), legend.key.width = unit(17, "mm"),
      plot.margin = margin(2, 2, 3, 4)
    )
  
  p_retention_delta <- ggplot(
    edge_stability_wide,
    aes(x = Delta_Percentage_Points, y = Edge_Label, fill = Change_Direction)
  ) +
    geom_vline(xintercept = 0, colour = "#60666B", linewidth = 0.45) +
    geom_col(width = 0.58) +
    geom_text(
      aes(
        label = sprintf("%+.1f", Delta_Percentage_Points),
        hjust = if_else(Delta_Percentage_Points >= 0, -0.14, 1.14)
      ),
      family = "Arial", size = 1.85, colour = text_col
    ) +
    scale_fill_manual(values = c(
      "Lower in SCTLD" = "#3D78A0",
      "Higher in SCTLD" = "#D36A58",
      "No change" = "#BFC4C7"
    )) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),
      labels = function(x) paste0(x, " pp"),
      expand = expansion(mult = c(0.20, 0.20))
    ) +
    labs(x = "SCTLD - Healthy retention (percentage points)", y = NULL) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.35, colour = "#33383C"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 5.7, colour = text_col),
      axis.title.x = element_text(size = 6.0, colour = text_col, margin = margin(t = 3)),
      legend.position = "none", plot.margin = margin(2, 5, 3, 2)
    )
  
  edge_stability_figure <- p_retention_heatmap + p_retention_delta +
    plot_layout(widths = c(0.92, 1.08)) +
    plot_annotation(
      title = "c  Stability of the 17-edge union",
      subtitle = paste0(
        "Cells show retention across 500 site-balanced bootstrap resamples; dot, retention >= 50% (not a P value). ",
        "Bars show descriptive SCTLD - Healthy change."
      ),
      theme = theme(
        plot.title = element_text(family = "Arial", size = 8.6, face = "bold", colour = text_col),
        plot.subtitle = element_text(family = "Arial", size = 6.1, colour = subtle_col),
        plot.margin = margin(3, 4, 1, 4)
      )
    )
  
  micro_stability_data <- edge_stability_plot_data %>%
    filter(Microviridae_Linked) %>%
    mutate(
      Partner = if_else(From == "Microviridae", To, From),
      Partner = factor(
        Partner,
        levels = c("Nanoviridae", "Schizomimiviridae", "Straboviridae")
      )
    )
  
  micro_segments <- micro_stability_data %>%
    select(Partner, Group_Label, Retention) %>%
    pivot_wider(names_from = Group_Label, values_from = Retention) %>%
    rename(
      Healthy_Retention = `Healthy / Normal`,
      SCTLD_Retention = SCTLD
    )
  
  micro_stability_figure <- ggplot() +
    geom_vline(
      xintercept = 0.50, colour = "#8B9398", linewidth = 0.45, linetype = "22"
    ) +
    geom_segment(
      data = micro_segments,
      aes(
        y = Partner, yend = Partner,
        x = Healthy_Retention, xend = SCTLD_Retention
      ),
      colour = "#D2A21D", linewidth = 1.25, lineend = "round"
    ) +
    geom_point(
      data = micro_stability_data,
      aes(x = Retention, y = Partner, fill = Group_Label),
      shape = 21, colour = "white", stroke = 0.45, size = 3.2
    ) +
    scale_fill_manual(values = c(
      "Healthy / Normal" = healthy_col,
      "SCTLD" = sctld_col
    )) +
    scale_x_continuous(
      limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = "Microviridae-linked association stability",
      subtitle = paste0(
        "Descriptive bootstrap retention; no edge-wise differential test was performed.\n",
        "Dashed line, 50% stability threshold."
      ),
      x = "Bootstrap retention", y = "Partner family", fill = NULL
    ) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#33383C"),
      axis.ticks = element_line(linewidth = 0.35),
      axis.text = element_text(size = 6.3, colour = text_col),
      axis.title = element_text(size = 6.7, colour = text_col),
      plot.title = element_text(size = 9.0, face = "bold", colour = text_col),
      plot.subtitle = element_text(size = 6.2, colour = subtle_col, margin = margin(b = 4)),
      legend.position = "top",
      legend.justification = "right",
      legend.text = element_text(size = 6.0),
      plot.margin = margin(5, 6, 5, 5)
    )
  
  metric_map <- c(
    Edges = "d  Retained associations",
    Active_Nodes = "e  Connected families",
    Components = "f  Connected components"
  )
  
  metric_data <- topology %>%
    select(Group, Resample, all_of(names(metric_map))) %>%
    pivot_longer(all_of(names(metric_map)), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric_Label = factor(unname(metric_map[Metric]), levels = unname(metric_map)),
      Group_Label = factor(
        recode(Group, Healthy = "Healthy / Normal (Exposed)", SCTLD = "SCTLD"),
        levels = c("Healthy / Normal (Exposed)", "SCTLD")
      )
    )
  
  metric_annotation <- permutation_summary %>%
    filter(Metric %in% names(metric_map)) %>%
    transmute(
      Metric,
      Annotation = paste0(
        "Permutation P = ", formatC(P_Two_Sided, format = "f", digits = 3),
        intToUtf8(10), "BH-adjusted P = ",
        formatC(P_BH_Across_8_Topology_Metrics, format = "f", digits = 3)
      ),
      Metric_Label = factor(unname(metric_map[Metric]), levels = unname(metric_map))
    ) %>%
    left_join(
      metric_data %>%
        group_by(Metric_Label) %>%
        summarise(y = max(Value, na.rm = TRUE), .groups = "drop"),
      by = "Metric_Label"
    ) %>%
    mutate(x = 1.5, y = y * 1.08)
  
  metric_summary <- metric_data %>%
    group_by(Metric, Metric_Label, Group_Label) %>%
    summarise(
      N_Resamples = n(),
      Median = median(Value),
      Q025 = quantile(Value, 0.025),
      Q25 = quantile(Value, 0.25),
      Q75 = quantile(Value, 0.75),
      Q975 = quantile(Value, 0.975),
      .groups = "drop"
    )
  
  metric_figure <- ggplot(
    metric_data,
    aes(x = Group_Label, y = Value, fill = Group_Label, colour = Group_Label)
  ) +
    geom_violin(width = 0.82, trim = FALSE, alpha = 0.18, linewidth = 0.35) +
    geom_point(
      position = position_jitter(width = 0.075, height = 0, seed = 20260813),
      alpha = 0.085, size = 0.40, stroke = 0
    ) +
    geom_boxplot(
      width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.38,
      colour = "#303438"
    ) +
    geom_text(
      data = metric_annotation,
      aes(x = x, y = y, label = Annotation),
      inherit.aes = FALSE, family = "Arial", size = 2.0,
      colour = subtle_col, lineheight = 0.95
    ) +
    facet_wrap(~ Metric_Label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(
      "Healthy / Normal (Exposed)" = healthy_col,
      "SCTLD" = sctld_col
    )) +
    scale_colour_manual(values = c(
      "Healthy / Normal (Exposed)" = healthy_col,
      "SCTLD" = sctld_col
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
    labs(
      title = "Topology of the full 18-family candidate set",
      subtitle = paste0(
        "Distributions show 500 descriptive site-balanced bootstrap resamples.\n",
        "P values use 999 within-location label permutations and BH correction across eight metrics."
      ),
      x = NULL, y = NULL
    ) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#33383C"),
      axis.ticks = element_line(linewidth = 0.35),
      axis.text.x = element_text(size = 5.75, colour = text_col, margin = margin(t = 2)),
      axis.text.y = element_text(size = 5.8, colour = text_col),
      strip.background = element_rect(fill = "#F2F3F3", colour = "#D1D5D7", linewidth = 0.35),
      strip.text = element_text(
        size = 6.35, face = "bold", colour = text_col,
        margin = margin(2.5, 2.5, 2.5, 2.5)
      ),
      plot.title = element_text(size = 8.6, face = "bold", colour = text_col),
      plot.subtitle = element_text(size = 6.2, colour = subtle_col, margin = margin(b = 3)),
      legend.position = "none",
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(3, 4, 4, 5)
    )
  
  # Isolate the nested, guide-collecting network patchwork before stacking it.
  # `free(network_figure)` renders as a blank page with patchwork 1.3.2 on the
  # devices used here, whereas wrap_elements(full = ...) is stable across the
  # Cairo, ragg and base-SVG backends.
  final_figure <- wrap_elements(full = network_figure) /
    wrap_elements(full = edge_stability_figure) /
    wrap_elements(full = metric_figure) +
    plot_layout(heights = c(1.22, 1.25, 0.92)) +
    plot_annotation(
      title = "Viral family association networks across coral health states",
      subtitle = paste0(
        "Project-matched DNA viromes from Orbicella faveolata coral tissue at five Florida locations.\n",
        "Panels a,b show the 11-family stable-edge union; topology tests retain all 18 families."
      ),
      theme = theme(
        plot.title = element_text(
          family = "Arial", size = 10.8, face = "bold", colour = text_col
        ),
        plot.subtitle = element_text(
          family = "Arial", size = 6.8, colour = subtle_col, margin = margin(b = 4)
        ),
        plot.margin = margin(4, 5, 4, 5)
      )
    )
  
  save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4
  
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(
        file = file.path(figure_dir, paste0(stem, ".svg")),
        width = width_in, height = height_in,
        system_fonts = list(sans = "Arial")
      )
    } else {
      grDevices::svg(
        filename = file.path(figure_dir, paste0(stem, ".svg")),
        width = width_in, height = height_in, family = "Arial", onefile = TRUE
      )
    }
    print(plot)
    dev.off()
  
    grDevices::cairo_pdf(
      filename = file.path(figure_dir, paste0(stem, ".pdf")),
      width = width_in, height = height_in, family = "Arial"
    )
    print(plot)
    dev.off()
  
    ragg::agg_png(
      filename = file.path(figure_dir, paste0(stem, ".png")),
      width = width_in, height = height_in, units = "in", res = 300,
      background = "white"
    )
    print(plot)
    dev.off()
  
    ragg::agg_tiff(
      filename = file.path(figure_dir, paste0(stem, ".tiff")),
      width = width_in, height = height_in, units = "in", res = dpi,
      compression = "lzw", background = "white"
    )
    print(plot)
    dev.off()
  }
  
  save_figure(final_figure, "Figure_1_network_composite", 183, 245)
  save_figure(network_figure, "Figure_2_panels_ab_union11_networks", 183, 98)
  save_figure(edge_stability_figure, "Figure_3_panel_c_retention_heatmap_delta", 160, 120)
  save_figure(metric_figure, "Figure_4_panels_def_topology_tests", 183, 78)
  save_figure(micro_stability_figure, "Figure_5_Microviridae_annotation", 125, 78)
  
  panel_a <- group_panel_data("Healthy")
  panel_b <- group_panel_data("SCTLD")
  network_node_cols <- c(
    "Family", "x", "y", "Stable_Degree", "Activity", "Node_Class", "LDA_Status"
  )
  network_edge_cols <- c(
    "Group", "Key", "From", "To", "Sign", "Retained_Replicates",
    "Resample_N", "Retention", "Mean_Abs_Rho", "Median_Abs_Rho"
  )
  
  write.csv(
    panel_a$nodes %>% select(all_of(network_node_cols)),
    file.path(figure_source_dir, "panel_a_healthy_union11_nodes.csv"), row.names = FALSE
  )
  write.csv(
    panel_a$edges %>% select(all_of(network_edge_cols)),
    file.path(figure_source_dir, "panel_a_healthy_union11_edges.csv"), row.names = FALSE
  )
  write.csv(
    panel_b$nodes %>% select(all_of(network_node_cols)),
    file.path(figure_source_dir, "panel_b_sctld_union11_nodes.csv"), row.names = FALSE
  )
  write.csv(
    panel_b$edges %>% select(all_of(network_edge_cols)),
    file.path(figure_source_dir, "panel_b_sctld_union11_edges.csv"), row.names = FALSE
  )
  write.csv(
    edge_stability_plot_data %>% arrange(Edge_Label, Group_Label),
    file.path(figure_source_dir, "panel_c_retention_heatmap_long_full17_union.csv"), row.names = FALSE
  )
  write.csv(
    edge_stability_wide %>% arrange(Delta_Retention, Edge_Label),
    file.path(figure_source_dir, "panel_c_retention_delta_wide_full17_union.csv"), row.names = FALSE
  )
  write.csv(
    metric_data %>% filter(Metric == "Edges"),
    file.path(figure_source_dir, "panel_d_retained_associations_full18.csv"), row.names = FALSE
  )
  write.csv(
    metric_data %>% filter(Metric == "Active_Nodes"),
    file.path(figure_source_dir, "panel_e_connected_families_full18.csv"), row.names = FALSE
  )
  write.csv(
    metric_data %>% filter(Metric == "Components"),
    file.path(figure_source_dir, "panel_f_connected_components_full18.csv"), row.names = FALSE
  )
  write.csv(
    metric_summary,
    file.path(figure_source_dir, "panel_d_f_bootstrap_summary_full18.csv"), row.names = FALSE
  )
  write.csv(
    permutation_summary,
    file.path(figure_source_dir, "panel_d_f_permutation_tests_full18.csv"), row.names = FALSE
  )
  walk2(
    c("Edges", "Active_Nodes", "Components"),
    c("d", "e", "f"),
    function(metric_name, panel_letter) {
      write.csv(
        metric_summary %>% filter(Metric == metric_name),
        file.path(
          figure_source_dir,
          paste0("panel_", panel_letter, "_bootstrap_summary_full18.csv")
        ),
        row.names = FALSE
      )
      write.csv(
        permutation_summary %>% filter(Metric == metric_name),
        file.path(
          figure_source_dir,
          paste0("panel_", panel_letter, "_permutation_test_full18.csv")
        ),
        row.names = FALSE
      )
    }
  )
  write.csv(
    micro_stability_data %>% arrange(Partner, Group_Label),
    file.path(figure_source_dir, "panel_Microviridae_edge_stability.csv"), row.names = FALSE
  )
  write.csv(
    lda_candidates_tbl,
    file.path(figure_source_dir, "panel_LDA_Microviridae_annotation.csv"), row.names = FALSE
  )
  write.csv(
    site_counts,
    file.path(figure_source_dir, "panel_sample_counts_by_location.csv"), row.names = FALSE
  )
  
  lda_legend <- sprintf(
    "(%s-enriched; LDA = %.2f; raw Kruskal-Wallis P = %.4f; within-stratum BH-adjusted P = %.3f)",
    lda_candidates_tbl$Enriched_In[[1]], lda_candidates_tbl$LDA_Signed[[1]],
    lda_candidates_tbl$KW_P[[1]], lda_candidates_tbl$KW_P_BH_within_stratum[[1]])
  writeLines(c(
    "# Figure legend",
    "",
    sub("{LDA_ANNOTATION}", lda_legend, "Fig. X | Viral family association networks across coral health states. a,b, Healthy/normal (exposed) (blue; n = 20) and SCTLD (red; n = 69) DNA viromes from Orbicella faveolata coral tissue at five Florida locations within PRJNA925892. Both panels show the 11-family union of stable-edge endpoints at identical coordinates. Healthy contains 17 stable associations and SCTLD contains 8; grey nodes have no stable edge in that state. Edge width denotes sign-specific retention across 500 site-balanced bootstrap resamples (20 samples per group per resample); stable associations had retention >= 50% after satisfying |rho| >= 0.7 and BH-adjusted P <= 0.05 in each resample. Gold fill marks Microviridae, the only named network family meeting the separate |LDA score| >= 2.5 display criterion {LDA_ANNOTATION}. LDA was not used to construct or filter the network, and this annotation is exploratory. c, Absolute bootstrap retention for the 17-edge union and descriptive change in retention from Healthy to SCTLD, expressed in percentage points. Dots mark cells with retention >= 50% and are not P-value markers; no edge-wise differential test was performed. d-f, Bootstrap distributions of retained associations, connected families and connected components calculated on the prespecified full 18-family candidate set. Boxes show medians and interquartile ranges; whiskers extend to 1.5 times the interquartile range; violins and points show 500 resamples per group. Empirical two-sided P values were calculated from 999 health-label permutations within location on the full reference networks, with BH correction across eight topology metrics. The a,b union display is descriptive and does not redefine the inferential network universe. Figure source data are provided in the accompanying package.", fixed = TRUE)
  ), file.path(output_dir, "Figure_legend.md"))
  
  writeLines(c(
    "Nature-figure QA notes",
    "",
    "- Backend: R only (ggplot2/patchwork; cairo_pdf and ragg exports).",
    "- Final composite size: 183 x 245 mm.",
    "- Separate figures are supplied as PDF, SVG, PNG and 600-dpi TIFF.",
    paste0(
      "- SVG device: ",
      if (requireNamespace("svglite", quietly = TRUE)) {
        "svglite (editable text)"
      } else {
        "base R SVG fallback; use the Cairo PDF as the primary editable vector file"
      },
      "."
    ),
    "- Both network panels use the original fixed coordinate map and the same 11-family stable-edge union.",
    "- Panels a,b retain families that are stable-edge endpoints in at least one health state (18 -> 11 shown per group).",
    "- panel_a_b_full18_display_audit.csv retains the complete family set, degrees and exact display rule.",
    "- panel_a_b_excluded7_evidence.csv retains stability evidence for the seven families omitted from panels a,b.",
    "- Topology distributions and all eight tests remain calculated on the full 18-family set; no new statistical filtering is applied.",
    "- Network construction did not use LDA. Microviridae is overlaid only after network construction.",
    "- LDA annotation threshold: |LDA| >= 2.5; only Microviridae is mapped to a named family node.",
    "- The annotation is extracted from the explicitly supplied current LEfSe source data, with exact sample-membership checks; no unclassified feature is mapped to the network.",
    "- Panel e explicitly reports the Active_Nodes/connected-families permutation P value and BH-adjusted P value.",
    "- Formal topology annotations use 999 within-location permutations and BH adjustment across eight metrics.",
    "- The 500 bootstrap distributions and edge retentions are descriptive stability summaries.",
    "- Panel c uses one decimal place for both absolute retention labels and retention-difference labels; values are unchanged.",
    "- Panel c dots denote retention >= 50%, not statistical significance; delta bars are percentage-point differences.",
    "- No edge-wise differential P values were calculated; panels c and Figure 5 are explicitly descriptive.",
    "- Healthy/SCTLD colour is reinforced by panel labels; the candidate node uses a separate gold fill.",
    "- Source data for every quantitative panel are in figure_source_data/."
  ), file.path(output_dir, "Figure_QA_notes.txt"))
  
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo_plot.txt"))
  
  cat("Created updated figures in:", figure_dir, "\n")
}

write_master_results <- function(analysis_dir, output_dir, mode) {
  metric_panel <- c(
    Active_Nodes = "e", Edges = "d", Avg_Degree = "",
    Density = "", Modularity = "", Components = "f",
    Largest_Component_Fraction = "", Positive_Edge_Ratio = ""
  )
  metric_unit <- c(
    Active_Nodes = "families", Edges = "associations", Avg_Degree = "degree",
    Density = "proportion", Modularity = "index", Components = "components",
    Largest_Component_Fraction = "proportion", Positive_Edge_Ratio = "proportion"
  )
  
  blank_result <- function(n_rows) {
    tibble(
      section = rep(NA_character_, n_rows), panel = rep(NA_character_, n_rows),
      record_id = rep(NA_character_, n_rows), group = rep(NA_character_, n_rows),
      comparison = rep(NA_character_, n_rows), entity_1 = rep(NA_character_, n_rows),
      entity_2 = rep(NA_character_, n_rows), metric = rep(NA_character_, n_rows),
      estimate = rep(NA_real_, n_rows), lower_2.5 = rep(NA_real_, n_rows),
      upper_97.5 = rep(NA_real_, n_rows), unit = rep(NA_character_, n_rows),
      p_value = rep(NA_real_, n_rows), p_adjusted_bh = rep(NA_real_, n_rows),
      n = rep(NA_real_, n_rows), sign = rep(NA_character_, n_rows),
      threshold = rep(NA_character_, n_rows), note = rep(NA_character_, n_rows),
      source_file = rep(NA_character_, n_rows)
    )
  }
  
  coerce_result_schema <- function(x) {
    template <- blank_result(0)
    missing <- setdiff(names(template), names(x))
    for (col in missing) {
      x[[col]] <- vctrs::vec_init(template[[col]], nrow(x))
    }
    x %>% select(all_of(names(template)))
  }
  
  process_dir <- file.path(analysis_dir, "process_data")
  figure_source_dir <- file.path(output_dir, "figure_source_data")
  
  sample_sizes <- read_csv(
    file.path(process_dir, "04_group_sample_sizes.csv"), show_col_types = FALSE
  ) %>%
    transmute(
      section = "sample_size", panel = "a,b",
      record_id = paste0("sample_size:", Group), group = Group,
      metric = "sample_count", estimate = as.numeric(N), unit = "samples",
      n = as.numeric(N),
      note = "Samples in the project-, ocean-, host- and compartment-matched cohort.",
      source_file = "process_data/04_group_sample_sizes.csv"
    ) %>% coerce_result_schema()
  
  reference_topology <- read_csv(
    file.path(process_dir, "06_reference_topology.csv"), show_col_types = FALSE
  ) %>%
    pivot_longer(-Group, names_to = "metric", values_to = "estimate") %>%
    transmute(
      section = "reference_topology", panel = if_else(metric %in% names(metric_panel), metric_panel[metric], ""),
      record_id = paste0("reference_topology:", Group, ":", metric),
      group = Group, metric, estimate = as.numeric(estimate),
      unit = if_else(metric %in% names(metric_unit), metric_unit[metric], "value"),
      note = "Topology of the full 18-family reference network.",
      source_file = "process_data/06_reference_topology.csv"
    ) %>% coerce_result_schema()
  
  bootstrap_raw <- read_csv(
    file.path(process_dir, "10_site_balanced_topology_summary.csv"),
    show_col_types = FALSE
  )
  bootstrap_topology <- map_dfr(names(metric_panel), function(metric_name) {
    tibble(
      section = "bootstrap_topology_summary",
      panel = unname(metric_panel[[metric_name]]),
      record_id = paste0("bootstrap:", bootstrap_raw$Group, ":", metric_name),
      group = bootstrap_raw$Group,
      metric = metric_name,
      estimate = as.numeric(bootstrap_raw[[paste0(metric_name, "_Median")]]),
      lower_2.5 = as.numeric(bootstrap_raw[[paste0(metric_name, "_Q025")]]),
      upper_97.5 = as.numeric(bootstrap_raw[[paste0(metric_name, "_Q975")]]),
      unit = unname(metric_unit[[metric_name]]), n = 500,
      note = "Median and 2.5th-97.5th percentiles across site-balanced bootstrap resamples.",
      source_file = "process_data/10_site_balanced_topology_summary.csv"
    )
  }) %>% coerce_result_schema()
  
  edge_retention <- read_csv(
    file.path(figure_source_dir, "panel_c_retention_heatmap_long_full17_union.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      section = "edge_retention", panel = "c",
      record_id = paste0("edge_retention:", Group, ":", Key), group = Group,
      entity_1 = From, entity_2 = To, metric = "bootstrap_retention",
      estimate = 100 * as.numeric(Retention), unit = "percent",
      n = as.numeric(Resample_N), sign = Sign, threshold = ">=50%",
      note = if_else(
        as.logical(Stable_At_50pct),
        "Stable at the prespecified 50% retention threshold.",
        "Below the prespecified 50% retention threshold."
      ),
      source_file = "figure_source_data/panel_c_retention_heatmap_long_full17_union.csv"
    ) %>% coerce_result_schema()
  
  edge_delta <- read_csv(
    file.path(figure_source_dir, "panel_c_retention_delta_wide_full17_union.csv"),
    show_col_types = FALSE
  ) %>%
    separate(Key, into = c("entity_1", "entity_2"), sep = "\\|\\|", remove = FALSE) %>%
    transmute(
      section = "edge_retention_delta", panel = "c",
      record_id = paste0("edge_delta:", Key),
      comparison = "SCTLD minus Healthy", entity_1, entity_2,
      metric = "retention_difference", estimate = as.numeric(Delta_Percentage_Points),
      unit = "percentage points", note = paste0(Change_Direction, "; descriptive only."),
      source_file = "figure_source_data/panel_c_retention_delta_wide_full17_union.csv"
    ) %>% coerce_result_schema()
  
  permutation_tests <- read_csv(
    file.path(process_dir, "20_within_location_permutation_summary.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      section = "topology_permutation_test",
      panel = if_else(Metric %in% names(metric_panel), metric_panel[Metric], ""),
      record_id = paste0("permutation:", Metric),
      comparison = "SCTLD minus Healthy", metric = Metric,
      estimate = as.numeric(Observed_Delta_SCTLD_minus_Healthy),
      lower_2.5 = as.numeric(Null_Q025), upper_97.5 = as.numeric(Null_Q975),
      unit = if_else(Metric %in% names(metric_unit), metric_unit[Metric], "value"),
      p_value = as.numeric(P_Two_Sided),
      p_adjusted_bh = as.numeric(P_BH_Across_8_Topology_Metrics),
      n = as.numeric(Valid_Permutations),
      note = "Two-sided empirical P value from within-location label permutations; BH across eight topology metrics.",
      source_file = "process_data/20_within_location_permutation_summary.csv"
    ) %>% coerce_result_schema()
  
  lda_raw <- read_csv(
    file.path(figure_source_dir, "panel_LDA_Microviridae_annotation.csv"),
    show_col_types = FALSE
  )
  lda_result <- lda_raw %>%
    transmute(
      section = "LDA_annotation", panel = "a,b",
      record_id = paste0("LDA:", LDA_Family), group = Enriched_In,
      comparison = "SCTLD versus Healthy", entity_1 = LDA_Family,
      metric = "LDA_score", estimate = as.numeric(LDA_Signed), unit = "LDA score",
      p_value = as.numeric(KW_P), p_adjusted_bh = as.numeric(KW_P_BH_within_stratum),
      threshold = paste0("|LDA| >= ", LDA_Display_Threshold),
      note = "Exploratory post hoc annotation; LDA was not used to build or filter the network.",
      source_file = "figure_source_data/panel_LDA_Microviridae_annotation.csv"
    ) %>% coerce_result_schema()
  
  lda_abundance <- bind_rows(
    lda_raw %>% transmute(group = "Healthy", estimate = Healthy_mean_CPM, n = Healthy_n),
    lda_raw %>% transmute(group = Enriched_In, estimate = Condition_mean_CPM, n = Condition_n)
  ) %>%
    mutate(
      section = "LDA_annotation_abundance", panel = "a,b",
      record_id = paste0("LDA_mean_CPM:", group), comparison = "SCTLD versus Healthy",
      entity_1 = lda_raw$LDA_Family[[1]], metric = "mean_CPM", unit = "CPM",
      note = "Group mean abundance extracted from the current LEfSe source data.",
      source_file = "figure_source_data/panel_LDA_Microviridae_annotation.csv"
    ) %>% coerce_result_schema()
  
  master_results <- bind_rows(
    sample_sizes, reference_topology, bootstrap_topology, edge_retention,
    edge_delta, permutation_tests, lda_result, lda_abundance
  ) %>%
    arrange(section, panel, group, record_id)
  
  results_dir <- file.path(output_dir, "results")
  write_csv(master_results, file.path(results_dir, "network_analysis_master_results.csv"), na = "")
  write_tsv(master_results, file.path(results_dir, "network_analysis_master_results.tsv"), na = "")
  
  writeLines(
    c(
      "Complete viral network pipeline",
      paste0("Mode: ", mode),
      paste0("Analysis source: ", analysis_dir),
      paste0("Output: ", output_dir),
      "Panel c display: one decimal place on absolute retention and retention-difference labels.",
      "No numerical values are rounded in the source data or master result table."
    ),
    file.path(output_dir, "00_complete_pipeline_summary.txt")
  )
  
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo_complete.txt"))
  cat("Complete pipeline outputs created in:", output_dir, "\n")
}

usage <- paste(
  'Usage:',
  '  Rscript run_network.R plot-only <analysis_package_dir> <output_dir> <lefse_source_data.tsv>',
  '  Rscript run_network.R full <input_dir> <output_dir> <lefse_source_data.tsv>',
  '',
  'plot-only: redraw figures and rebuild the master table from saved process data.',
  'full: rerun the fixed manuscript analysis from four network inputs and current LEfSe output.',
  'First run: Rscript ../lefse_analysis/run_lefse.R ../lefse_analysis/inputs <lefse_output>',
  'Then supply <lefse_output>/results/source_data.tsv as the final argument.',
  'If omitted, only a same-run sibling ../lefse_analysis/results/source_data.tsv is accepted.',
  sep = '\n'
)
cli <- commandArgs(trailingOnly = TRUE)
if (length(cli) == 1L && cli[[1]] %in% c('-h', '--help')) {
  cat(usage, '\n')
  quit(status = 0L)
}
if (!length(cli) %in% c(3L, 4L)) stop(usage, call. = FALSE)
mode <- tolower(cli[[1]])
if (!mode %in% c('plot-only', 'full')) {
  stop('Mode must be plot-only or full.\n', usage, call. = FALSE)
}

explicit_source_dir <- normalizePath(cli[[2]], winslash = '/', mustWork = TRUE)
output_dir <- normalizePath(cli[[3]], winslash = '/', mustWork = FALSE)
# Run LEfSe first and supply its freshly computed source table explicitly.
# The default path also supports results stored beside this output directory.
lefse_source <- if (length(cli) == 4L) cli[[4]] else
  file.path(dirname(output_dir), 'lefse_analysis', 'results', 'source_data.tsv')
if (!file.exists(lefse_source))
  stop('Current LEfSe source data are required.\n', usage, call. = FALSE)
lefse_source <- normalizePath(lefse_source, winslash = '/', mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, 'results'), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({ library(tidyverse) })
annotation_path <- derive_lda_annotation(lefse_source, file.path(output_dir, 'provenance', 'lefse'))
analysis_dir <- explicit_source_dir
if (identical(mode, 'full')) {
  input_dir <- explicit_source_dir
  required_inputs <- c(
    '24083_metagenome_TPM.tsv', 'dna_virus_taxonomy.tsv', 'sample_info_metagenome.tsv',
    'network_union_layout.csv'
  )
  missing_inputs <- required_inputs[!file.exists(file.path(input_dir, required_inputs))]
  if (length(missing_inputs) > 0L) {
    stop('Missing input files: ', paste(missing_inputs, collapse = ', '))
  }
  snapshot_dir <- file.path(output_dir, 'input_data')
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
  source_paths <- normalizePath(file.path(input_dir, required_inputs), winslash = '/', mustWork = TRUE)
  destination_paths <- file.path(snapshot_dir, required_inputs)
  same_path <- tolower(source_paths) == tolower(
    normalizePath(destination_paths, winslash = '/', mustWork = FALSE)
  )
  if (any(!same_path)) {
    copied <- file.copy(source_paths[!same_path], destination_paths[!same_path], overwrite = TRUE)
    if (!all(copied)) stop('Failed to create the complete input snapshot.')
  }
  run_network_analysis(c(snapshot_dir, output_dir))
  analysis_dir <- output_dir
}

run_network_figure(c(analysis_dir, output_dir, annotation_path))
write_master_results(analysis_dir, output_dir, mode)
cat('Standalone workflow completed:', output_dir, '\n')
