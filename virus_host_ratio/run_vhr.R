#!/usr/bin/env Rscript
# Complete raw-input-to-statistics-to-figure VHR workflow.
# Usage: Rscript run_vhr.R [input_dir] [output_root]
# RStudio: source("/path/to/run_vhr.R") runs using the sibling inputs/ directory.
# Project-aware comparisons use only the input snapshot in inputs/.
# Status harmonisation is defined in this script; it is specific to this VHR analysis.
# No precomputed statistics or second R script are needed.

.vhr_script_dir <- local({
  frames <- sys.frames()
  ofiles <- lapply(frames, function(f) f$ofile)
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L, ofiles)
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  path <- if (length(ofiles)) tail(ofiles, 1L)[[1L]] else
          if (length(arg)) sub("^--file=", "", arg[1L]) else "run_vhr.R"
  dirname(normalizePath(path, winslash = "/", mustWork = FALSE))
})

run_vhr <- function(input_dir = file.path(.vhr_script_dir, "inputs"),
                    output_root = .vhr_script_dir) {
  # Preserve Unicode metadata when Windows starts in a non-UTF-8 C locale.
  # Restore the caller's locale when this function finishes.
  previous_ctype <- Sys.getlocale("LC_CTYPE")
  if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
    utf8_ctype <- suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
    if (!nzchar(utf8_ctype)) stop("Run this analysis in a UTF-8 R locale.")
    on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", previous_ctype)), add = TRUE)
  }
  # Write CSV text explicitly as UTF-8 on all platforms.
  write.csv <- function(..., fileEncoding = "UTF-8") {
    utils::write.csv(..., fileEncoding = fileEncoding)
  }
  run_started <- Sys.time()
#!/usr/bin/env Rscript

###############################################################################
# FINAL Project-aware total linkage-informed virus-to-host abundance ratio
#
# Analysis design and reproducibility rules
# -----------------------------------------------------------------------------
# 1. Metagenome only.
# 2. Statistical comparisons are strictly within:
#      Analysis_Project + Ocean_Basin + Species + Sample_Type
#    so independent Projects can never be pooled to reach n >= 3.
# 3. Mixed / recovery / unknown statuses are excluded explicitly.
# 4. Healthy/Normal is the reference group.
# 5. Only samples with linked-virus signal > 0 AND linked-host signal > 0
#    are retained for the ratio analysis.
# 6. Main statistical test:
#      two-sided Wilcoxon rank-sum test
#    on log10(total linkage-informed virus-to-host abundance ratio).
# 7. Multiple testing:
#      primary = BH correction WITHIN each
#      Analysis_Project x Ocean_Basin x Species x Sample_Type context;
#      supplementary = global BH across all pairwise status comparisons.
# 8. Boxplots and points use matched dodge widths. Point jitter is deliberately
#    tiny and reproducible, so points stay centered on their own box.
# 9. If the same Species + Ocean_Basin + Sample_Type is represented by multiple
#    Projects among eligible comparisons, the plotted x label automatically
#    becomes Project_Species to prevent visual merging.
#
# IMPORTANT INTERPRETATION
# -----------------------------------------------------------------------------
# This is an AGGREGATE linkage-informed viral-to-host abundance ratio:
#
#   sum(abundance of linked viruses)
#   --------------------------------
#   sum(abundance of linked hosts)
#
# It is NOT a pair-specific virus/host ratio.
#
# The biological meaning of the ratio also depends on the abundance scale of
# virus_abundance.tsv and host_abundance.tsv. If the two tables were independently
# compositional-normalized (e.g. each separately sums to 1 or 1e6), interpret
# this as a relative virus-to-host SIGNAL ratio rather than an absolute VHR.
###############################################################################


# ==============================================================================
# 0. Paths and parameters
# ==============================================================================

BASE_DIR <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
OUTDIR <- file.path(output_root, "results")
FIGDIR <- file.path(output_root, "figures")
PROVDIR <- file.path(output_root, "provenance")
for (d in c(OUTDIR, FIGDIR, PROVDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
previous_wd <- getwd()
on.exit(setwd(previous_wd), add = TRUE)
setwd(BASE_DIR)

PHAGE_OTU_FILE <- "virus_abundance.tsv"
HOST_OTU_FILE  <- "host_abundance.tsv"
META_FILE      <- "sample_metadata.tsv"
LINK_FILE      <- "virus_host_links.tsv"

REF_STATUS <- "Healthy/Normal"

MIN_SAMPLE_SIZE <- 3L

# Tight point jitter around the CENTER of each status-specific box.
# Set to 0 for points exactly on the box center line.
POINT_JITTER_WIDTH <- 0.025

# Box/point dodge must be identical.
DODGE_WIDTH <- 0.78

set.seed(123)

# OUTDIR is the portable results directory defined above.

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "ggh4x",
  "svglite",
  "ragg",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {

  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install them with install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), ")) and rerun.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ggh4x)
})


# ==============================================================================
# 2. Helper functions
# ==============================================================================

normalize_host_id <- function(x) {

  x <- stringr::str_trim(
    as.character(x)
  )

  x <- gsub(
    "^(GCA_[0-9]+)_([0-9]+)$",
    "\\1.\\2",
    x
  )

  x <- gsub(
    "^(GCF_[0-9]+)_([0-9]+)$",
    "\\1.\\2",
    x
  )

  x
}


read_otu_table <- function(file) {

  df <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    encoding = "UTF-8"
  )

  names(df)[1] <- "ID"

  df$ID <- as.character(
    df$ID
  )

  sample_cols <- setdiff(
    names(df),
    "ID"
  )

  df[sample_cols] <- lapply(
    df[sample_cols],
    function(x) {

      y <- suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )

      y[
        !is.finite(y)
      ] <- 0

      y
    }
  )

  df %>%

    group_by(
      ID
    ) %>%

    summarise(
      across(
        all_of(sample_cols),
        ~sum(
          .x,
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    )
}


# Conservative status normalization.
# Mixed/recovery/unknown categories are deliberately excluded.
normalize_status <- function(x) {

  s <- stringr::str_trim(
    as.character(x)
  )

  dplyr::case_when(

    is.na(s) |
      s == "" ~ NA_character_,

    stringr::str_detect(
      s,
      stringr::regex(
        "Mix|Recovery|Treated|Unknown|Unkonwn",
        ignore_case = TRUE
      )
    ) ~ NA_character_,

    stringr::str_detect(
      s,
      stringr::regex(
        "Healthy|Normal",
        ignore_case = TRUE
      )
    ) ~ REF_STATUS,

    stringr::str_detect(
      s,
      stringr::regex(
        "SCTLD",
        ignore_case = TRUE
      )
    ) ~ "Disease: SCTLD",

    stringr::str_detect(
      s,
      stringr::regex(
        "black band",
        ignore_case = TRUE
      )
    ) ~ "Disease: black band disease",

    stringr::str_detect(
      s,
      stringr::regex(
        "white plague",
        ignore_case = TRUE
      )
    ) ~ "Disease: white plague disease",

    stringr::str_detect(
      s,
      stringr::regex(
        "Dark Spot|DSS",
        ignore_case = TRUE
      )
    ) ~ "Disease: Dark Spot Syndrome",

    stringr::str_detect(
      s,
      stringr::regex(
        "Bleaching lesion",
        ignore_case = TRUE
      )
    ) ~ "Disease: Bleaching lesion",

    stringr::str_detect(
      s,
      stringr::regex(
        "Pathogen",
        ignore_case = TRUE
      )
    ) ~ "Disease: Pathogen-induced disease",

    stringr::str_detect(
      s,
      stringr::regex(
        "Bleached",
        ignore_case = TRUE
      )
    ) ~ "Bleached",

    stringr::str_detect(
      s,
      stringr::regex(
        "Heat",
        ignore_case = TRUE
      )
    ) &
      stringr::str_detect(
        s,
        stringr::regex(
          "Acidification",
          ignore_case = TRUE
        )
      ) ~ "Stressed: Combined",

    stringr::str_detect(
      s,
      stringr::regex(
        "Heat",
        ignore_case = TRUE
      )
    ) ~ "Stressed: Heat stress",

    stringr::str_detect(
      s,
      stringr::regex(
        "Acidification",
        ignore_case = TRUE
      )
    ) ~ "Stressed: Acidification",

    stringr::str_detect(
      s,
      stringr::regex(
        "Cold",
        ignore_case = TRUE
      )
    ) ~ "Stressed: Cold stress",

    TRUE ~ NA_character_
  )
}


short_status <- function(x) {

  dplyr::case_when(

    x == "Disease: Pathogen-induced disease" ~ "PID",
    x == "Disease: SCTLD" ~ "SCTLD",
    x == "Disease: Bleaching lesion" ~ "BL",
    x == "Disease: Dark Spot Syndrome" ~ "DSS",
    x == "Disease: black band disease" ~ "BBD",
    x == "Disease: white plague disease" ~ "WPD",
    x == "Bleached" ~ "Bleached",
    x == "Stressed: Heat stress" ~ "HS",
    x == "Stressed: Acidification" ~ "AS",
    x == "Stressed: Cold stress" ~ "CS",
    x == "Stressed: Combined" ~ "Combined",

    TRUE ~ x
  )
}


# ==============================================================================
# 3. Read data
# ==============================================================================

phage_otu <- read_otu_table(
  PHAGE_OTU_FILE
)

host_otu <- read_otu_table(
  HOST_OTU_FILE
) %>%

  mutate(
    ID = normalize_host_id(
      ID
    )
  )


meta <- read.delim(
  META_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  encoding = "UTF-8"
)


link_df <- read.delim(
  LINK_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  encoding = "UTF-8"
)


# ==============================================================================
# 4. Metadata: metagenome only + Project-aware context
# ==============================================================================

required_meta_columns <- c(
  "Sample",
  "Project",
  "Sample_Type",
  "Ocean_Basin",
  "Status"
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  names(meta)
)

if (length(missing_meta_columns) > 0) {

  stop(
    paste0(
      "Metadata is missing required columns:\n",
      paste(
        missing_meta_columns,
        collapse = "\n"
      )
    )
  )
}


if (!"Species" %in% names(meta)) {

  if (
    "Coral_Host_Species" %in%
      names(meta)
  ) {

    meta <- meta %>%

      mutate(
        Species =
          Coral_Host_Species
      )

  } else {

    stop(
      "Metadata contains neither Species nor Coral_Host_Species."
    )
  }
}


if ("Data_Type" %in% names(meta)) {

  meta <- meta %>%

    filter(
      stringr::str_to_lower(
        Data_Type
      ) == "metagenome"
    )
}


meta <- meta %>%

  mutate(

    Sample =
      stringr::str_trim(
        as.character(Sample)
      ),

    Project =
      stringr::str_trim(
        as.character(Project)
      ),

    Analysis_Project =
      Project,

    Species =
      stringr::str_trim(
        as.character(Species)
      ),

    Sample_Type =
      stringr::str_trim(
        as.character(Sample_Type)
      ),

    Ocean_Basin =
      stringr::str_trim(
        as.character(Ocean_Basin)
      ),

    Status_Original =
      as.character(Status),

    Status_Detail =
      normalize_status(
        Status
      )
  )


# Explicitly remove statuses that are not part of our inferential framework.
meta <- meta %>%

  filter(
    !is.na(
      Status_Detail
    ),
    !is.na(
      Analysis_Project
    ),
    Analysis_Project != "",
    !is.na(
      Ocean_Basin
    ),
    Ocean_Basin != "",
    !is.na(
      Species
    ),
    Species != "",
    !is.na(
      Sample_Type
    ),
    Sample_Type != ""
  )


# ==============================================================================
# 5. Match samples across virus table, host table and metadata
# ==============================================================================

sample_cols <- Reduce(
  intersect,
  list(
    setdiff(
      names(phage_otu),
      "ID"
    ),
    setdiff(
      names(host_otu),
      "ID"
    ),
    meta$Sample
  )
)


if (length(sample_cols) == 0) {

  stop(
    "No samples are shared by phage_OTU, Host_OTU and metadata."
  )
}


phage_otu <- phage_otu %>%

  select(
    ID,
    all_of(
      sample_cols
    )
  )


host_otu <- host_otu %>%

  select(
    ID,
    all_of(
      sample_cols
    )
  )


meta <- meta %>%

  filter(
    Sample %in%
      sample_cols
  )


cat(
  ">>> Matched metagenome samples: ",
  length(sample_cols),
  "\n",
  sep = ""
)


# ==============================================================================
# 6. QC: abundance-table sample sums
#
# This does NOT prove the abundance units, but helps diagnose whether each OTU
# table may have been independently normalized to a constant compositional sum.
# ==============================================================================

phage_total_all <- colSums(
  as.matrix(
    phage_otu[
      ,
      setdiff(
        names(phage_otu),
        "ID"
      ),
      drop = FALSE
    ]
  ),
  na.rm = TRUE
)


host_total_all <- colSums(
  as.matrix(
    host_otu[
      ,
      setdiff(
        names(host_otu),
        "ID"
      ),
      drop = FALSE
    ]
  ),
  na.rm = TRUE
)


abundance_scale_qc <- data.frame(
  Sample = sample_cols,
  Total_All_Phage_OTU = phage_total_all[
    sample_cols
  ],
  Total_All_Host_OTU = host_total_all[
    sample_cols
  ],
  stringsAsFactors = FALSE
)


write.csv(
  abundance_scale_qc,
  file.path(
    OUTDIR,
    "QC_abundance_table_sample_sums.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Parse virus-host linkages
# ==============================================================================

link_df_parsed <- link_df %>%

  transmute(

    Virus =
      as.character(
        Virus
      ),

    Host_ID =
      normalize_host_id(
        `Host genome`
      ),

    Virus_taxonomy =
      as.character(
        `Virus taxonomy`
      ),

    Host_taxonomy =
      as.character(
        `Host taxonomy`
      )
  ) %>%

  distinct(
    Virus,
    Host_ID,
    .keep_all = TRUE
  )


linked_viruses <- unique(
  link_df_parsed$Virus
)

linked_hosts <- unique(
  link_df_parsed$Host_ID
)


matched_linked_viruses <- linked_viruses[
  linked_viruses %in%
    phage_otu$ID
]

matched_linked_hosts <- linked_hosts[
  linked_hosts %in%
    host_otu$ID
]


diagnostic_summary <- data.frame(

  total_unique_linkages =
    nrow(
      link_df_parsed
    ),

  unique_viruses_in_linkage_table =
    length(
      linked_viruses
    ),

  unique_hosts_in_linkage_table =
    length(
      linked_hosts
    ),

  matched_linked_viruses =
    length(
      matched_linked_viruses
    ),

  matched_linked_hosts =
    length(
      matched_linked_hosts
    ),

  matched_virus_percent =
    length(
      matched_linked_viruses
    ) /
    length(
      linked_viruses
    ) *
    100,

  matched_host_percent =
    length(
      matched_linked_hosts
    ) /
    length(
      linked_hosts
    ) *
    100
)


write.csv(
  diagnostic_summary,
  file.path(
    OUTDIR,
    "QC_linkage_matching_diagnostic_summary.csv"
  ),
  row.names = FALSE
)


print(
  diagnostic_summary
)


# ==============================================================================
# 8. Calculate linked-virus / linked-host pool signal per sample
# ==============================================================================

virus_signal <- phage_otu %>%

  filter(
    ID %in%
      matched_linked_viruses
  ) %>%

  summarise(
    across(
      all_of(
        sample_cols
      ),
      ~sum(
        .x,
        na.rm = TRUE
      )
    )
  ) %>%

  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Linked_Virus_Signal"
  )


host_signal <- host_otu %>%

  filter(
    ID %in%
      matched_linked_hosts
  ) %>%

  summarise(
    across(
      all_of(
        sample_cols
      ),
      ~sum(
        .x,
        na.rm = TRUE
      )
    )
  ) %>%

  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Linked_Host_Signal"
  )


n_detected_viruses <- phage_otu %>%

  filter(
    ID %in%
      matched_linked_viruses
  ) %>%

  pivot_longer(
    cols = -ID,
    names_to = "Sample",
    values_to = "Virus_Abundance"
  ) %>%

  group_by(
    Sample
  ) %>%

  summarise(
    n_detected_linked_viruses =
      n_distinct(
        ID[
          Virus_Abundance > 0
        ]
      ),
    .groups = "drop"
  )


n_detected_hosts <- host_otu %>%

  filter(
    ID %in%
      matched_linked_hosts
  ) %>%

  pivot_longer(
    cols = -ID,
    names_to = "Sample",
    values_to = "Host_Abundance"
  ) %>%

  group_by(
    Sample
  ) %>%

  summarise(
    n_detected_linked_hosts =
      n_distinct(
        ID[
          Host_Abundance > 0
        ]
      ),
    .groups = "drop"
  )


# ==============================================================================
# 9. Optional but useful diagnostic:
#    number of predicted virus-host linkages co-detected in each sample
# ==============================================================================

matched_pairs <- link_df_parsed %>%

  filter(
    Virus %in%
      phage_otu$ID,
    Host_ID %in%
      host_otu$ID
  ) %>%

  distinct(
    Virus,
    Host_ID
  )


phage_matrix <- as.matrix(
  phage_otu[
    ,
    sample_cols,
    drop = FALSE
  ]
)

rownames(
  phage_matrix
) <- phage_otu$ID


host_matrix <- as.matrix(
  host_otu[
    ,
    sample_cols,
    drop = FALSE
  ]
)

rownames(
  host_matrix
) <- host_otu$ID


n_codetected_pairs <- vapply(

  sample_cols,

  function(s) {

    if (
      nrow(
        matched_pairs
      ) == 0
    ) {

      return(
        0L
      )
    }

    sum(
      phage_matrix[
        matched_pairs$Virus,
        s
      ] > 0 &
        host_matrix[
          matched_pairs$Host_ID,
          s
        ] > 0,
      na.rm = TRUE
    )
  },

  FUN.VALUE = integer(1)
)


codetection_df <- data.frame(
  Sample = sample_cols,
  n_codetected_predicted_pairs =
    as.integer(
      n_codetected_pairs[
        sample_cols
      ]
    ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 10. Merge signals and zero-status diagnosis
# ==============================================================================

total_signal_raw <- tibble(
  Sample = sample_cols
) %>%

  left_join(
    virus_signal,
    by = "Sample"
  ) %>%

  left_join(
    host_signal,
    by = "Sample"
  ) %>%

  left_join(
    n_detected_viruses,
    by = "Sample"
  ) %>%

  left_join(
    n_detected_hosts,
    by = "Sample"
  ) %>%

  left_join(
    codetection_df,
    by = "Sample"
  ) %>%

  mutate(

    Linked_Virus_Signal =
      replace_na(
        Linked_Virus_Signal,
        0
      ),

    Linked_Host_Signal =
      replace_na(
        Linked_Host_Signal,
        0
      ),

    n_detected_linked_viruses =
      replace_na(
        n_detected_linked_viruses,
        0
      ),

    n_detected_linked_hosts =
      replace_na(
        n_detected_linked_hosts,
        0
      ),

    n_codetected_predicted_pairs =
      replace_na(
        n_codetected_predicted_pairs,
        0
      ),

    Zero_Status =
      case_when(

        Linked_Virus_Signal > 0 &
          Linked_Host_Signal > 0 ~
          "Retained: virus>0 and host>0",

        Linked_Virus_Signal == 0 &
          Linked_Host_Signal > 0 ~
          "Filtered: virus=0",

        Linked_Virus_Signal > 0 &
          Linked_Host_Signal == 0 ~
          "Filtered: host=0",

        Linked_Virus_Signal == 0 &
          Linked_Host_Signal == 0 ~
          "Filtered: virus=0 and host=0",

        TRUE ~
          "Filtered: unknown"
      ),

    Pass_Nonzero_Filter =
      Linked_Virus_Signal > 0 &
      Linked_Host_Signal > 0
  ) %>%

  left_join(
    meta %>%

      select(
        Sample,
        Analysis_Project,
        Project,
        Ocean_Basin,
        Species,
        Sample_Type,
        Status_Original,
        Status_Detail,
        everything()
      ),
    by = "Sample"
  ) %>%

  arrange(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type,
    Status_Detail,
    Sample
  )


write.csv(
  total_signal_raw,
  file.path(
    OUTDIR,
    "Total_linkage_informed_VHR_each_sample_raw_with_zero_status.csv"
  ),
  row.names = FALSE
)


zero_status_summary <- total_signal_raw %>%

  count(
    Zero_Status,
    name = "n_samples"
  ) %>%

  mutate(
    percent =
      n_samples /
      sum(
        n_samples
      ) *
      100
  )


write.csv(
  zero_status_summary,
  file.path(
    OUTDIR,
    "QC_zero_status_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Nonzero ratio
# ==============================================================================

total_vhr_each_sample_nonzero <- total_signal_raw %>%

  filter(
    Pass_Nonzero_Filter
  ) %>%

  mutate(

    Linkage_informed_total_VHR =
      Linked_Virus_Signal /
      Linked_Host_Signal,

    log10_Linkage_informed_total_VHR =
      log10(
        Linkage_informed_total_VHR
      )
  ) %>%

  filter(
    is.finite(
      log10_Linkage_informed_total_VHR
    )
  )


write.csv(
  total_vhr_each_sample_nonzero,
  file.path(
    OUTDIR,
    "Total_linkage_informed_VHR_each_sample_nonzero_only.csv"
  ),
  row.names = FALSE
)


cat(
  ">>> Nonzero samples retained: ",
  nrow(
    total_vhr_each_sample_nonzero
  ),
  "\n",
  sep = ""
)


# ==============================================================================
# 12. Project-confounding QC
# ==============================================================================

cross_project_context_qc <- total_vhr_each_sample_nonzero %>%

  distinct(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type
  ) %>%

  group_by(
    Ocean_Basin,
    Species,
    Sample_Type
  ) %>%

  summarise(
    N_Analysis_Project =
      n_distinct(
        Analysis_Project
      ),
    Projects =
      paste(
        sort(
          unique(
            Analysis_Project
          )
        ),
        collapse = "; "
      ),
    .groups = "drop"
  ) %>%

  filter(
    N_Analysis_Project > 1
  )


write.csv(
  cross_project_context_qc,
  file.path(
    OUTDIR,
    "QC_same_ecological_context_across_multiple_projects.csv"
  ),
  row.names = FALSE
)


status_counts_projectaware <- total_vhr_each_sample_nonzero %>%

  count(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type,
    Status_Detail,
    name = "n_nonzero"
  )


write.csv(
  status_counts_projectaware,
  file.path(
    OUTDIR,
    "QC_status_counts_by_project_context.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Summary by Project + ecological context + status
# ==============================================================================

context_summary_nonzero <- total_vhr_each_sample_nonzero %>%

  group_by(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type,
    Status_Detail
  ) %>%

  summarise(

    n_samples_nonzero =
      n_distinct(
        Sample
      ),

    sum_linked_virus_signal =
      sum(
        Linked_Virus_Signal,
        na.rm = TRUE
      ),

    mean_linked_virus_signal =
      mean(
        Linked_Virus_Signal,
        na.rm = TRUE
      ),

    median_linked_virus_signal =
      median(
        Linked_Virus_Signal,
        na.rm = TRUE
      ),

    sum_linked_host_signal =
      sum(
        Linked_Host_Signal,
        na.rm = TRUE
      ),

    mean_linked_host_signal =
      mean(
        Linked_Host_Signal,
        na.rm = TRUE
      ),

    median_linked_host_signal =
      median(
        Linked_Host_Signal,
        na.rm = TRUE
      ),

    mean_total_VHR =
      mean(
        Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    median_total_VHR =
      median(
        Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    sd_total_VHR =
      sd(
        Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    mean_log10_total_VHR =
      mean(
        log10_Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    median_log10_total_VHR =
      median(
        log10_Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    sd_log10_total_VHR =
      sd(
        log10_Linkage_informed_total_VHR,
        na.rm = TRUE
      ),

    mean_n_detected_linked_viruses =
      mean(
        n_detected_linked_viruses,
        na.rm = TRUE
      ),

    mean_n_detected_linked_hosts =
      mean(
        n_detected_linked_hosts,
        na.rm = TRUE
      ),

    mean_n_codetected_predicted_pairs =
      mean(
        n_codetected_predicted_pairs,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%

  arrange(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type,
    Status_Detail
  )


write.csv(
  context_summary_nonzero,
  file.path(
    OUTDIR,
    "Total_linkage_informed_VHR_by_project_context_state_nonzero_only.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Eligibility:
#     each status n >= 3 WITHIN Project + region + species + compartment
# ==============================================================================

valid_subgroups <- context_summary_nonzero %>%

  filter(
    n_samples_nonzero >=
      MIN_SAMPLE_SIZE
  ) %>%

  select(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type,
    Status_Detail
  )


data_qc <- total_vhr_each_sample_nonzero %>%

  semi_join(
    valid_subgroups,
    by = c(
      "Analysis_Project",
      "Ocean_Basin",
      "Species",
      "Sample_Type",
      "Status_Detail"
    )
  )


keep_groups <- data_qc %>%

  group_by(
    Analysis_Project,
    Ocean_Basin,
    Species,
    Sample_Type
  ) %>%

  summarise(

    has_ref =
      REF_STATUS %in%
      Status_Detail,

    n_status =
      n_distinct(
        Status_Detail
      ),

    .groups = "drop"
  ) %>%

  filter(
    has_ref,
    n_status > 1
  )


plot_df <- data_qc %>%

  inner_join(
    keep_groups,
    by = c(
      "Analysis_Project",
      "Ocean_Basin",
      "Species",
      "Sample_Type"
    )
  )


write.csv(
  keep_groups,
  file.path(
    OUTDIR,
    "QC_eligible_project_contexts.csv"
  ),
  row.names = FALSE
)


write.csv(
  plot_df,
  file.path(
    OUTDIR,
    "Total_linkage_informed_VHR_QC_filtered_values_for_testing.csv"
  ),
  row.names = FALSE
)


cat(
  ">>> QC-filtered samples for testing: ",
  nrow(
    plot_df
  ),
  "\n",
  sep = ""
)


# ==============================================================================
# 15. Wilcoxon rank-sum tests
#
# Main correction:
#   BH within Analysis_Project + Ocean_Basin + Species + Sample_Type
#
# Supplementary correction:
#   global BH across all comparisons
# ==============================================================================

groups_to_test <- plot_df %>%

  distinct(
    Analysis_Project,
    Ocean_Basin,
    Sample_Type,
    Species
  )


results_list <- list()

counter <- 1L


if (
  nrow(
    groups_to_test
  ) > 0
) {

  for (
    i in seq_len(
      nrow(
        groups_to_test
      )
    )
  ) {

    cur_project <-
      groups_to_test$Analysis_Project[i]

    cur_ocean <-
      groups_to_test$Ocean_Basin[i]

    cur_type <-
      groups_to_test$Sample_Type[i]

    cur_species <-
      groups_to_test$Species[i]


    sub_df <- plot_df %>%

      filter(
        Analysis_Project ==
          cur_project,
        Ocean_Basin ==
          cur_ocean,
        Sample_Type ==
          cur_type,
        Species ==
          cur_species
      )


    ref_vals <- sub_df %>%

      filter(
        Status_Detail ==
          REF_STATUS
      ) %>%

      pull(
        log10_Linkage_informed_total_VHR
      )


    other_statuses <- setdiff(
      unique(
        as.character(
          sub_df$Status_Detail
        )
      ),
      REF_STATUS
    )


    for (
      stat in other_statuses
    ) {

      test_vals <- sub_df %>%

        filter(
          Status_Detail ==
            stat
        ) %>%

        pull(
          log10_Linkage_informed_total_VHR
        )


      if (
        length(
          ref_vals
        ) <
          MIN_SAMPLE_SIZE ||
          length(
            test_vals
          ) <
          MIN_SAMPLE_SIZE
      ) {

        next
      }


      p_val <- tryCatch(

        wilcox.test(
          ref_vals,
          test_vals,
          alternative = "two.sided",
          exact = FALSE
        )$p.value,

        error = function(e) {

          NA_real_
        }
      )


      if (
        is.na(
          p_val
        )
      ) {

        next
      }


      mean_ref <- mean(
        ref_vals,
        na.rm = TRUE
      )

      mean_test <- mean(
        test_vals,
        na.rm = TRUE
      )

      median_ref <- median(
        ref_vals,
        na.rm = TRUE
      )

      median_test <- median(
        test_vals,
        na.rm = TRUE
      )


      delta_median <-
        median_test -
        median_ref


      direction <- dplyr::case_when(

        delta_median > 0 ~
          "Higher in test group",

        delta_median < 0 ~
          "Lower in test group",

        TRUE ~
          "No change"
      )


      results_list[[counter]] <- data.frame(

        Analysis_Project =
          cur_project,

        Ocean_Basin =
          cur_ocean,

        Sample_Type =
          cur_type,

        Species =
          cur_species,

        Group_Ref =
          REF_STATUS,

        Group_Test =
          stat,

        Sample_Size_Ref =
          length(
            ref_vals
          ),

        Sample_Size_Test =
          length(
            test_vals
          ),

        Mean_Ref =
          mean_ref,

        Mean_Test =
          mean_test,

        Median_Ref =
          median_ref,

        Median_Test =
          median_test,

        Delta_Median =
          delta_median,

        Direction =
          direction,

        Test =
          "Two-sided Wilcoxon rank-sum",

        P_Value =
          p_val,

        stringsAsFactors = FALSE
      )


      counter <- counter + 1L
    }
  }
}


stats_table <- dplyr::bind_rows(
  results_list
)


if (
  nrow(
    stats_table
  ) > 0
) {

  stats_table <- stats_table %>%

    group_by(
      Analysis_Project,
      Ocean_Basin,
      Sample_Type,
      Species
    ) %>%

    mutate(

      # PRIMARY multiplicity correction
      P_Adjust_Within_Context =
        p.adjust(
          P_Value,
          method = "BH"
        ),

      Significance_Within_Context =
        case_when(

          P_Adjust_Within_Context <
            0.001 ~ "***",

          P_Adjust_Within_Context <
            0.01 ~ "**",

          P_Adjust_Within_Context <
            0.05 ~ "*",

          P_Adjust_Within_Context <
            0.1 ~ "trend",

          TRUE ~ "ns"
        )
    ) %>%

    ungroup() %>%

    mutate(

      # Supplementary sensitivity / reference only
      P_Adjust_Global =
        p.adjust(
          P_Value,
          method = "BH"
        ),

      Significance_Global =
        case_when(

          P_Adjust_Global <
            0.001 ~ "***",

          P_Adjust_Global <
            0.01 ~ "**",

          P_Adjust_Global <
            0.05 ~ "*",

          P_Adjust_Global <
            0.1 ~ "trend",

          TRUE ~ "ns"
        ),

      P_Value_Formatted =
        format.pval(
          P_Value,
          digits = 3
        ),

      P_Adjust_Within_Context_Formatted =
        format.pval(
          P_Adjust_Within_Context,
          digits = 3
        ),

      P_Adjust_Global_Formatted =
        format.pval(
          P_Adjust_Global,
          digits = 3
        )
    ) %>%

    arrange(
      P_Adjust_Within_Context,
      Analysis_Project,
      Ocean_Basin,
      Species,
      Sample_Type,
      Group_Test
    )
}


write.csv(
  stats_table,
  file.path(
    OUTDIR,
    "Total_linkage_informed_VHR_Wilcoxon_ProjectAware_within_context_BH.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Prevent visual merging of equivalent contexts from different Projects
# ==============================================================================

project_overlap <- plot_df %>%

  distinct(
    Analysis_Project,
    Ocean_Basin,
    Sample_Type,
    Species
  ) %>%

  group_by(
    Ocean_Basin,
    Sample_Type,
    Species
  ) %>%

  summarise(
    N_Project =
      n_distinct(
        Analysis_Project
      ),
    .groups = "drop"
  )


plot_df <- plot_df %>%

  left_join(
    project_overlap,
    by = c(
      "Ocean_Basin",
      "Sample_Type",
      "Species"
    )
  ) %>%

  mutate(

    Display_Species =
      if_else(
        N_Project > 1,
        paste0(
          Analysis_Project,
          "_",
          Species
        ),
        Species
      )
  )


if (
  nrow(
    stats_table
  ) > 0
) {

  stats_table <- stats_table %>%

    left_join(
      project_overlap,
      by = c(
        "Ocean_Basin",
        "Sample_Type",
        "Species"
      )
    ) %>%

    mutate(
      Display_Species =
        if_else(
          N_Project > 1,
          paste0(
            Analysis_Project,
            "_",
            Species
          ),
          Species
        )
    )
}


# ==============================================================================
# 17. Plot — explicit box/point x coordinates
#
# Explicit point and box alignment
# ------------------------------------------------------------------------------
# position_dodge2() and position_jitterdodge() calculate horizontal offsets
# independently. When different coral taxa contain different numbers of status
# groups, their apparent centres can differ slightly.
#
# Here we DO NOT use automatic dodge/jitter at all.
#
# For every:
#   Project + Ocean_Basin + Sample_Type + Species + Status
#
# the script calculates ONE explicit Box_X coordinate.
#
# Every replicate point is then placed around that exact Box_X using a small,
# deterministic and symmetric offset that is always INSIDE the corresponding
# box. Therefore a point can no longer "belong" visually to the wrong box.
# ==============================================================================

if (
  nrow(
    plot_df
  ) > 0
) {

  # ---------------------------------------------------------------------------
  # 17.1 Status order and colours
  # ---------------------------------------------------------------------------

  status_levels <- c(
    REF_STATUS,
    sort(
      setdiff(
        unique(
          as.character(
            plot_df$Status_Detail
          )
        ),
        REF_STATUS
      )
    )
  )


  plot_df$Status_Detail <- factor(
    plot_df$Status_Detail,
    levels = status_levels
  )


  status_palette <- c(
    "Healthy/Normal" = "#4E9F6D",
    "Bleached" = "#E69F00",
    "Disease: SCTLD" = "#C44E52",
    "Disease: Bleaching lesion" = "#D17A5D",
    "Disease: Pathogen-induced disease" = "#8C6BB1",
    "Disease: Dark Spot Syndrome" = "#B07AA1",
    "Disease: black band disease" = "#A65628",
    "Disease: white plague disease" = "#E15759",
    "Stressed: Heat stress" = "#4E79A7",
    "Stressed: Acidification" = "#76B7B2",
    "Stressed: Cold stress" = "#59A14F",
    "Stressed: Combined" = "#9C755F"
  )


  missing_status_colours <- setdiff(
    status_levels,
    names(
      status_palette
    )
  )


  if (
    length(
      missing_status_colours
    ) > 0
  ) {

    extra_cols <- grDevices::hcl.colors(
      length(
        missing_status_colours
      ),
      palette = "Dark 3"
    )

    names(
      extra_cols
    ) <- missing_status_colours

    status_palette <- c(
      status_palette,
      extra_cols
    )
  }


  status_palette <- status_palette[
    status_levels
  ]


  # ---------------------------------------------------------------------------
  # 17.2 Explicit x geometry
  # ---------------------------------------------------------------------------

  # Distance between coral-taxon centres.
  # Larger than 1 so even taxa with several status boxes cannot overlap.
  SPECIES_SPACING <- 1.35

  # Distance between adjacent status-specific boxes for the same coral taxon.
  STATUS_SPACING <- 0.24

  # Width of one box.
  BOX_WIDTH <- 0.18

  # Replicate points are kept within +/- this distance around Box_X.
  # BOX_WIDTH / 2 = 0.09, so +/-0.055 remains safely inside the box.
  POINT_SPREAD_MAX <- 0.055


  # ---------------------------------------------------------------------------
  # Unique species centres.
  #
  # X coordinates are made globally unique across facets. This lets a single
  # scale_x_continuous() map each numeric break to the correct coral-taxon label.
  # ---------------------------------------------------------------------------

  species_map <- plot_df %>%

    distinct(
      Analysis_Project,
      Ocean_Basin,
      Sample_Type,
      Species,
      Display_Species
    ) %>%

    arrange(
      Ocean_Basin,
      Sample_Type,
      Analysis_Project,
      Display_Species
    )


  panel_map <- species_map %>%

    count(
      Ocean_Basin,
      Sample_Type,
      name = "N_Species"
    ) %>%

    arrange(
      Ocean_Basin,
      Sample_Type
    ) %>%

    mutate(
      Panel_Order =
        row_number(),

      Previous_Species =
        cumsum(
          lag(
            N_Species + 2L,
            default = 0L
          )
        ),

      Panel_X_Start =
        Previous_Species *
        SPECIES_SPACING
    )


  species_map <- species_map %>%

    left_join(
      panel_map %>%

        select(
          Ocean_Basin,
          Sample_Type,
          Panel_Order,
          Panel_X_Start
        ),
      by = c(
        "Ocean_Basin",
        "Sample_Type"
      )
    ) %>%

    group_by(
      Ocean_Basin,
      Sample_Type
    ) %>%

    arrange(
      Analysis_Project,
      Display_Species,
      .by_group = TRUE
    ) %>%

    mutate(
      Species_Local_Order =
        row_number(),

      X_Base =
        Panel_X_Start +
        Species_Local_Order *
        SPECIES_SPACING
    ) %>%

    ungroup()


  plot_df <- plot_df %>%

    left_join(
      species_map %>%

        select(
          Analysis_Project,
          Ocean_Basin,
          Sample_Type,
          Species,
          Display_Species,
          X_Base
        ),
      by = c(
        "Analysis_Project",
        "Ocean_Basin",
        "Sample_Type",
        "Species",
        "Display_Species"
      )
    )


  # ---------------------------------------------------------------------------
  # Exact centre of every status-specific box
  # ---------------------------------------------------------------------------

  status_position_map <- plot_df %>%

    distinct(
      Analysis_Project,
      Ocean_Basin,
      Sample_Type,
      Species,
      Display_Species,
      X_Base,
      Status_Detail
    ) %>%

    mutate(
      Status_Order =
        match(
          as.character(
            Status_Detail
          ),
          status_levels
        )
    ) %>%

    group_by(
      Analysis_Project,
      Ocean_Basin,
      Sample_Type,
      Species,
      Display_Species,
      X_Base
    ) %>%

    arrange(
      Status_Order,
      .by_group = TRUE
    ) %>%

    mutate(
      N_Status =
        n(),

      Status_Rank =
        row_number(),

      Status_Offset =
        (
          Status_Rank -
          (
            N_Status + 1
          ) / 2
        ) *
        STATUS_SPACING,

      Box_X =
        X_Base +
        Status_Offset
    ) %>%

    ungroup()


  plot_df <- plot_df %>%

    left_join(
      status_position_map %>%

        select(
          Analysis_Project,
          Ocean_Basin,
          Sample_Type,
          Species,
          Display_Species,
          Status_Detail,
          Box_X
        ),
      by = c(
        "Analysis_Project",
        "Ocean_Basin",
        "Sample_Type",
        "Species",
        "Display_Species",
        "Status_Detail"
      )
    )


  # ---------------------------------------------------------------------------
  # Deterministic point locations INSIDE each box
  #
  # n=1:       0
  # n=2:      -a, +a
  # n=3:      -a, 0, +a
  # n=4:      evenly spaced from -a to +a
  #
  # There is no random jitter and no automatic dodge.
  # ---------------------------------------------------------------------------

  plot_df <- plot_df %>%

    group_by(
      Analysis_Project,
      Ocean_Basin,
      Sample_Type,
      Species,
      Display_Species,
      Status_Detail
    ) %>%

    arrange(
      log10_Linkage_informed_total_VHR,
      Sample,
      .by_group = TRUE
    ) %>%

    mutate(
      N_Points =
        n(),

      Point_Rank =
        row_number(),

      Point_Centered_Rank =
        Point_Rank -
        (
          N_Points + 1
        ) / 2,

      Point_Max_Abs_Rank =
        max(
          abs(
            Point_Centered_Rank
          )
        ),

      Point_Offset =
        if_else(
          N_Points <= 1 |
            Point_Max_Abs_Rank == 0,
          0,
          Point_Centered_Rank /
            Point_Max_Abs_Rank *
            POINT_SPREAD_MAX
        ),

      Point_X =
        Box_X +
        Point_Offset
    ) %>%

    ungroup()


  # ---------------------------------------------------------------------------
  # QC export: exact x positions used by boxes and points
  # ---------------------------------------------------------------------------

  write.csv(
    status_position_map,
    file.path(
      OUTDIR,
      "QC_plot_exact_box_x_positions.csv"
    ),
    row.names = FALSE
  )


  write.csv(
    plot_df %>%

      select(
        Sample,
        Analysis_Project,
        Ocean_Basin,
        Sample_Type,
        Species,
        Display_Species,
        Status_Detail,
        X_Base,
        Box_X,
        Point_X,
        Point_Offset,
        log10_Linkage_informed_total_VHR
      ),
    file.path(
      OUTDIR,
      "QC_plot_exact_point_x_positions.csv"
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 17.3 Significance labels
  #
  # Labels are placed above the EXACT centre of the corresponding test-status
  # box, not merely above the general species centre.
  # Primary significance = within-context BH.
  # ---------------------------------------------------------------------------

  if (
    nrow(
      stats_table
    ) > 0
  ) {

    stat_label_df <- stats_table %>%

      mutate(

        label =
          case_when(

            P_Adjust_Within_Context <
              0.001 ~ "***",

            P_Adjust_Within_Context <
              0.01 ~ "**",

            P_Adjust_Within_Context <
              0.05 ~ "*",

            P_Adjust_Within_Context <
              0.1 ~ "trend",

            TRUE ~ "ns"
          ),

        direction_symbol =
          case_when(

            Direction ==
              "Higher in test group" ~ "↑",

            Direction ==
              "Lower in test group" ~ "↓",

            TRUE ~ ""
          ),

        Status_short =
          short_status(
            Group_Test
          ),

        plot_label =
          paste0(
            Status_short,
            " ",
            direction_symbol,
            label
          )
      ) %>%

      filter(
        label != "ns"
      ) %>%

      left_join(
        status_position_map %>%

          transmute(
            Analysis_Project,
            Ocean_Basin,
            Sample_Type,
            Species,
            Display_Species,
            Group_Test =
              as.character(
                Status_Detail
              ),
            Label_X =
              Box_X
          ),
        by = c(
          "Analysis_Project",
          "Ocean_Basin",
          "Sample_Type",
          "Species",
          "Display_Species",
          "Group_Test"
        )
      ) %>%

      left_join(
        plot_df %>%

          group_by(
            Analysis_Project,
            Ocean_Basin,
            Sample_Type,
            Species,
            Display_Species
          ) %>%

          summarise(
            base_y =
              max(
                log10_Linkage_informed_total_VHR,
                na.rm = TRUE
              ),
            .groups = "drop"
          ),
        by = c(
          "Analysis_Project",
          "Ocean_Basin",
          "Sample_Type",
          "Species",
          "Display_Species"
        )
      ) %>%

      group_by(
        Analysis_Project,
        Ocean_Basin,
        Sample_Type,
        Species,
        Display_Species
      ) %>%

      arrange(
        Group_Test,
        .by_group = TRUE
      ) %>%

      mutate(
        y.position =
          base_y +
          0.23 +
          (
            row_number() - 1
          ) *
          0.22
      ) %>%

      ungroup()

  } else {

    stat_label_df <- tibble()
  }


  # ---------------------------------------------------------------------------
  # 17.4 Figure size
  # ---------------------------------------------------------------------------

  n_total_groups <- species_map %>%

    nrow()


  calc_width <- max(
    12,
    n_total_groups *
      1.0 +
      3.5
  )

  calc_height <- 7.5


  # ---------------------------------------------------------------------------
  # 17.5 Plot
  #
  # IMPORTANT:
  # no position_dodge()
  # no position_dodge2()
  # no position_jitterdodge()
  #
  # x is already final and explicit.
  # ---------------------------------------------------------------------------

  p_total <- ggplot() +

    geom_boxplot(
      data = plot_df,
      aes(
        x = Box_X,
        y = log10_Linkage_informed_total_VHR,
        fill = Status_Detail,
        group = interaction(
          Analysis_Project,
          Ocean_Basin,
          Sample_Type,
          Species,
          Status_Detail,
          drop = TRUE
        )
      ),
      outlier.shape = NA,
      alpha = 0.78,
      width = BOX_WIDTH,
      linewidth = 0.48,
      colour = "#333333"
    ) +

    geom_point(
      data = plot_df,
      aes(
        x = Point_X,
        y = log10_Linkage_informed_total_VHR,
        fill = Status_Detail
      ),
      size = 2.35,
      alpha = 0.78,
      colour = "#222222",
      shape = 21,
      stroke = 0.25
    ) +

    ggh4x::facet_nested(
      ~ Ocean_Basin +
        Sample_Type,
      scales = "free_x",
      space = "free_x",
      nest_line =
        element_line(
          linewidth = 0.45,
          colour = "grey55"
        )
    ) +

    scale_x_continuous(
      breaks =
        species_map$X_Base,
      labels =
        species_map$Display_Species,
      expand =
        expansion(
          add = 0.5,
          mult = c(
            0.04,
            0.04
          )
        )
    ) +

    scale_fill_manual(
      values =
        status_palette,
      drop =
        FALSE
    ) +

    labs(
      title =
        "Aggregate linkage-informed virus-to-host abundance ratio across coral status categories",
      subtitle =
        paste0(
          "Metagenome-only; comparisons restricted within Project × sampling region × coral taxon × compartment. ",
          "Two-sided Wilcoxon rank-sum tests; labels use within-context BH-adjusted P values."
        ),
      x =
        NULL,
      y =
        expression(
          Log[10]~
          "(aggregate linkage-informed virus-to-host abundance ratio)"
        ),
      fill =
        "Coral status"
    ) +

    theme_bw(
      base_size = 10,
      base_family = "Arial"
    ) +

    theme(

      strip.background =
        element_rect(
          fill = "white",
          colour = "#4D4D4D"
        ),

      strip.text =
        element_text(
          face = "bold",
          size = 9
        ),

      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          face = "italic",
          size = 9
        ),

      axis.text.y =
        element_text(
          size = 9
        ),

      axis.title.y =
        element_text(
          size = 10
        ),

      panel.grid.major.x =
        element_blank(),

      panel.grid.minor =
        element_blank(),

      panel.grid.major.y =
        element_line(
          linewidth = 0.30,
          colour = "#E6E6E6"
        ),

      legend.position =
        "bottom",

      plot.title =
        element_text(
          face = "bold",
          size = 12
        ),

      plot.subtitle =
        element_text(
          size = 8.5
        )
    )


  if (
    nrow(
      stat_label_df
    ) > 0
  ) {

    p_total <- p_total +

      geom_text(
        data = stat_label_df,
        aes(
          x = Label_X,
          y = y.position,
          label = plot_label
        ),
        inherit.aes = FALSE,
        size = 3.15,
        fontface = "bold"
      )
  }


  # Export only: Cairo embeds fonts and preserves Unicode direction arrows.
  grDevices::cairo_pdf(file.path(FIGDIR, "VHR.pdf"),
                      width = calc_width, height = calc_height,
                      family = "Arial", bg = "white")
  print(p_total)
  grDevices::dev.off()
  ggsave(file.path(FIGDIR, "VHR.png"), p_total,
         width = calc_width, height = calc_height, dpi = 300,
         device = ragg::agg_png, bg = "white")
  svglite::svglite(file.path(FIGDIR, "VHR.svg"),
                   width = calc_width, height = calc_height, bg = "white")
  print(p_total)
  grDevices::dev.off()
  ggsave(file.path(FIGDIR, "VHR.tiff"), p_total,
         width = calc_width, height = calc_height, dpi = 600,
         device = ragg::agg_tiff, compression = "lzw", bg = "white")
  write.csv(stat_label_df, file.path(OUTDIR, "QC_plot_significance_labels.csv"),
            row.names = FALSE)

}


# ==============================================================================
# 18. Interpretation note
# ==============================================================================

interpretation_note <- c(

  "Metric:",
  "Aggregate linkage-informed virus-to-host abundance ratio =",
  "sum abundance of all viruses represented in the virus-host linkage table /",
  "sum abundance of all hosts represented in the virus-host linkage table.",
  "",

  "This is NOT a pair-specific virus-to-host ratio.",
  "",

  "Statistical design:",
  "Comparisons are performed only within the same Analysis_Project, Ocean_Basin,",
  "Species and Sample_Type. Each retained status group contains at least 3",
  "nonzero samples, including a Healthy/Normal reference.",
  "",

  "Primary statistical test:",
  "two-sided Wilcoxon rank-sum test on log10-transformed aggregate ratios.",
  "",

  "Primary multiple-testing correction:",
  "Benjamini-Hochberg adjustment across Healthy/Normal-vs-other-status",
  "comparisons within each Analysis_Project x Ocean_Basin x Species x Sample_Type",
  "context.",
  "",

  "A global BH-adjusted P value is also exported as a supplementary reference.",
  "",

  "Abundance-scale caveat:",
  "If virus_abundance.tsv and host_abundance.tsv were independently normalized to separate",
  "compositional totals, interpret this metric as a relative viral-to-host signal",
  "ratio rather than an absolute virus-to-host abundance ratio."
)


writeLines(
  interpretation_note,
  file.path(
    OUTDIR,
    "README_metric_and_statistics_interpretation.txt"
  ),
  useBytes = TRUE
)


cat("\n")
cat("============================================================\n")
cat("PROJECT-AWARE VHR ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Test: two-sided Wilcoxon rank-sum\n"
)

cat(
  "Primary correction: BH within Project × region × species × compartment\n"
)

cat(
  "Point/box alignment: explicit x coordinates; no automatic dodge/jitter\n"
)

cat(
  "Output directory: ",
  OUTDIR,
  "\n",
  sep = ""
)

cat("============================================================\n")

# ============================================================================
# 19. Consolidated source data and machine-readable run record
# ============================================================================
# One sample row per matched, status-eligible accession (including zero signals).
# Comparison rows and group-summary rows have their own record_type; values are
# never copied across a sample x comparison Cartesian join.
sample_source <- total_signal_raw %>%
  left_join(total_vhr_each_sample_nonzero %>%
              select(Sample, Linkage_informed_total_VHR,
                     log10_Linkage_informed_total_VHR), by = "Sample") %>%
  mutate(In_Testing_And_Figure = Sample %in% plot_df$Sample,
         record_type = "sample", Record_ID = paste0("sample:", Sample))
geometry_cols <- intersect(c("Sample", "Display_Species", "X_Base", "Box_X",
                            "Point_X", "Point_Offset"), names(plot_df))
if (length(geometry_cols) > 1L) {
  sample_source <- left_join(sample_source, plot_df[, geometry_cols], by = "Sample")
}
comparison_source <- stats_table %>%
  mutate(record_type = "comparison", Record_ID = paste0("comparison:", row_number()))
group_source <- context_summary_nonzero %>%
  mutate(record_type = "group_summary", Record_ID = paste0("group:", row_number()))
source_data <- bind_rows(sample_source, comparison_source, group_source) %>%
  select(record_type, Record_ID, everything())
write.table(source_data, file.path(OUTDIR, "source_data.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
input_names <- c(PHAGE_OTU_FILE, HOST_OTU_FILE, META_FILE, LINK_FILE)
runtime_inputs <- data.frame(file = input_names,
                            bytes = file.info(input_names)$size,
                            md5 = unname(tools::md5sum(input_names)))
write.table(runtime_inputs, file.path(PROVDIR, "runtime_input_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- list(
  run_started_utc = format(run_started, tz = "UTC", usetz = TRUE),
  run_completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  viral_input = PHAGE_OTU_FILE,
  matched_samples = length(sample_cols),
  positive_signal_samples = nrow(total_vhr_each_sample_nonzero),
  plotted_samples = nrow(plot_df),
  eligible_contexts = nrow(keep_groups), comparisons = nrow(stats_table),
  within_context_q_below_0_05 = sum(stats_table$P_Adjust_Within_Context < 0.05),
  global_q_below_0_05 = sum(stats_table$P_Adjust_Global < 0.05),
  source_data_records = nrow(source_data),
  source_data_columns = ncol(source_data),
  source_data_record_counts = as.list(table(source_data$record_type)),
  matched_linked_virus_signal_sum = sum(total_signal_raw$Linked_Virus_Signal),
  matched_linked_host_signal_sum = sum(total_signal_raw$Linked_Host_Signal),
  figure_width_inches = if (exists("calc_width")) calc_width else NA_real_,
  figure_height_inches = if (exists("calc_height")) calc_height else NA_real_,
  status_rule = "VHR status harmonisation includes nutrient-enrichment Healthy/DSS labels",
  packages = as.list(vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
jsonlite::write_json(summary, file.path(PROVDIR, "run_summary.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 16)
writeLines(capture.output(sessionInfo()), file.path(PROVDIR, "sessionInfo.txt"))
cat("Consolidated source data: ", file.path(OUTDIR, "source_data.tsv"), "\n", sep = "")
invisible(summary)
}

.vhr_args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
if (length(.vhr_args) > 2L) stop("Usage: Rscript run_vhr.R [input_dir] [output_root]")
run_vhr(input_dir = if (length(.vhr_args) >= 1L) .vhr_args[1L] else file.path(.vhr_script_dir, "inputs"),
        output_root = if (length(.vhr_args) >= 2L) .vhr_args[2L] else .vhr_script_dir)
