#!/usr/bin/env Rscript
# Coral virome: one complete R pipeline, from SIX TPM/taxonomy/metadata inputs
# through normalization, project-aware LEfSe/KW calculations and the DNA/RNA plot.
# No precomputed LEfSe table, RDS object, second R script or Python is required.
# CLI: Rscript run_lefse.R [input_directory] [output_directory]
# Defaults: inputs/ beside this script; results/ and figures/ beside this script.
# RStudio: source("run_lefse.R"); run_pipeline()
# README.md documents the input schemas, methods, result fields and environment.
# Original scientific framework retained; DSS gates require exact original labels.

.script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[1]) else NULL
}
.package_dir <- if (is.null(.script_file)) getwd() else
  dirname(normalizePath(.script_file, winslash = "/", mustWork = TRUE))

run_pipeline <- function(input_dir = file.path(.package_dir, "inputs"),
                         output_dir = .package_dir) {
  # The UTF-8 locale is needed by labels and grid text on Windows.
  if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
    old_ctype <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", ".UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", old_ctype), add = TRUE)
  }
  required_packages <- c("Biobase", "BiocGenerics", "S4Vectors",
    "SummarizedExperiment", "lefser", "dplyr", "tidyr", "ggplot2", "cowplot",
    "svglite", "ragg", "systemfonts", "jsonlite")
  missing <- required_packages[!vapply(required_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages first (see README): ",
                            paste(missing, collapse = ", "))
  if (!capabilities("cairo")) stop("R with Cairo PDF support is required.")
  if (!any(tolower(systemfonts::system_fonts()$family) == "arial"))
    stop("Install Arial to preserve the verified figure layout, then rerun.")
  suppressPackageStartupMessages({
    library(SummarizedExperiment); library(S4Vectors); library(lefser)
    library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
    library(svglite); library(ragg); library(grid)
  })
  old_options <- options(stringsAsFactors = FALSE, scipen = 999)
  on.exit(options(old_options), add = TRUE)
  old_null_device <- cowplot::set_null_device("agg")
  on.exit(cowplot::set_null_device(old_null_device), add = TRUE)
  INPUT_DIR <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  RESULT_DIR <- file.path(output_root, "results")
  FIGURE_DIR <- file.path(output_root, "figures")
  # Keep one consolidated main table and the original audit exports separately.
  OUTPUT_DIR <- file.path(RESULT_DIR, "audit_tables")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  MIN_N <- 3L
  set.seed(123)
  cat("Input directory:", INPUT_DIR, "\nOutput directory:", output_root, "\n")
  # A font-aware scratch device measures grobs without an implicit Rplots.pdf.
  scratch_pdf <- tempfile(fileext = ".pdf")
  cairo_pdf(scratch_pdf, family = "Arial")
  scratch_device <- dev.cur()
  on.exit({
    if (scratch_device %in% dev.list()) dev.off(scratch_device)
    unlink(scratch_pdf)
  }, add = TRUE)

# 1. Required input files
# ==============================================================================

REQUIRED_FILES <- c(
  "24083_metagenome_TPM.tsv",
  "dna_virus_taxonomy.tsv",
  "sample_info_metagenome.tsv",
  "rna_virus_tpm.tsv",
  "rna_virus_taxonomy.tsv",
  "sample_info_metatranscriptome.tsv"
)

missing_files <- REQUIRED_FILES[
  !file.exists(file.path(INPUT_DIR, REQUIRED_FILES))
]

if (length(missing_files) > 0) {

  stop(
    paste0(
      "The following required input files are missing from:\n",
      INPUT_DIR,
      "\n\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat("All six required input files were found.\n\n")


# ==============================================================================
# 3. Analysis settings
# ==============================================================================

LDA_THRESHOLD <- 2
KW_THRESHOLD  <- 0.05

# Main-figure display threshold.
# LEfSe itself is run at the official/default |LDA| >= 2 threshold,
# while the final publication figure shows markers at |LDA| >= 2.5.
LDA_THRESHOLD_PLOT <- 2.5

FIG_WIDTH_MM  <- 205
FIG_HEIGHT_MM <- 235
BASE_FONT     <- "Arial"

# IMPORTANT:
# Location is intentionally NOT included.
STRATA_COLUMNS <- c(
  "Analysis_Project",
  "Ocean_Basin",
  "Coral_Host_Species",
  "Sample_Type"
)

# Preserve the comparison scope of the original LEfSe analysis.
DNA_CONDITIONS <- c(
  "SCTLD",
  "Bleached",
  "Heat_Stress",
  "DSS"
)

RNA_CONDITIONS <- c(
  "Heat_Stress"
)


# ==============================================================================
# 4. RNA special case: Avila-Magana et al. 2021
#
# These PRJNA accessions belong to the same study source and are therefore
# treated as one Analysis_Project.
# ==============================================================================

AVILA_PROJECTS <- c(

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


make_analysis_project <- function(project, virome_label) {

  project <- trimws(
    as.character(project)
  )

  project[
    is.na(project) |
      project == ""
  ] <- "Unknown_Project"


  # RNA special case only
  if (virome_label == "RNA") {

    project[
      project %in% AVILA_PROJECTS
    ] <- "Avila2021"
  }

  project
}


# ==============================================================================
# 5. Status normalization
# ==============================================================================

normalise_status <- function(x) {

  x <- trimws(
    as.character(x)
  )

  dplyr::case_when(

    is.na(x) |
      x == "" ~ NA_character_,

    # Explicitly exclude mixed, recovery and unknown categories
    grepl(
      "Mix",
      x,
      ignore.case = TRUE
    ) ~ NA_character_,

    grepl(
      "Recovery|Treated",
      x,
      ignore.case = TRUE
    ) ~ NA_character_,

    grepl(
      "Unknown|Unkonwn",
      x,
      ignore.case = TRUE
    ) ~ NA_character_,

    # Healthy
    grepl(
      "Healthy|Normal",
      x,
      ignore.case = TRUE
    ) ~ "Healthy",

    # Disease / stress groups used in this LEfSe analysis
    grepl(
      "SCTLD",
      x,
      ignore.case = TRUE
    ) ~ "SCTLD",

    grepl(
      "Bleach",
      x,
      ignore.case = TRUE
    ) ~ "Bleached",

    grepl(
      "Heat",
      x,
      ignore.case = TRUE
    ) ~ "Heat_Stress",

    grepl(
      "Dark Spot|DSS",
      x,
      ignore.case = TRUE
    ) ~ "DSS",

    TRUE ~ NA_character_
  )
}


# ==============================================================================
# 6. Taxonomy cleanup
# ==============================================================================

fix_taxonomy <- function(tax_df) {

  ranks <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family"
  )

  for (rank in ranks) {

    if (!rank %in% names(tax_df)) {

      tax_df[[rank]] <- "Unclassified"
    }

    tax_df[[rank]] <- as.character(
      tax_df[[rank]]
    )

    tax_df[[rank]][
      is.na(tax_df[[rank]]) |
        trimws(tax_df[[rank]]) == ""
    ] <- "Unclassified"
  }


  is_unknown <- function(x) {

    grepl(
      "Unclassified",
      x,
      ignore.case = TRUE
    ) |
      x == "NA"
  }


  tax_df$Family <- apply(

    tax_df[
      ,
      ranks,
      drop = FALSE
    ],

    1,

    function(row) {

      if (!is_unknown(row[["Family"]])) {

        return(
          row[["Family"]]
        )
      }


      if (!is_unknown(row[["Order"]])) {

        return(
          paste0(
            "Unclass. Order (",
            row[["Order"]],
            ")"
          )
        )
      }


      if (!is_unknown(row[["Class"]])) {

        return(
          paste0(
            "Unclass. Class (",
            row[["Class"]],
            ")"
          )
        )
      }


      if (!is_unknown(row[["Phylum"]])) {

        return(
          paste0(
            "Unclass. Phylum (",
            row[["Phylum"]],
            ")"
          )
        )
      }


      if (!is_unknown(row[["Kingdom"]])) {

        return(
          paste0(
            "Unclass. Kingdom (",
            row[["Kingdom"]],
            ")"
          )
        )
      }


      "Unclassified_Root"
    }
  )

  tax_df
}


# ==============================================================================
# 7. Read one virome dataset
# ==============================================================================

read_virome <- function(
    label,
    abundance_file,
    taxonomy_file,
    metadata_file
) {

  cat("\n------------------------------------------------------------\n")
  cat("Reading", label, "dataset\n")
  cat("------------------------------------------------------------\n")


  abundance <- read.delim(
    file.path(
      INPUT_DIR,
      abundance_file
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )


  taxonomy <- read.delim(
    file.path(
      INPUT_DIR,
      taxonomy_file
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )


  metadata <- read.delim(
    file.path(
      INPUT_DIR,
      metadata_file
    ),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )


  names(abundance)[1] <- "Genome"
  names(taxonomy)[1]  <- "Genome"


  # ---------------------------------------------------------------------------
  # Check metadata columns
  # ---------------------------------------------------------------------------

  required_metadata_columns <- c(
    "Sample",
    "Project",
    "Ocean_Basin",
    "Coral_Host_Species",
    "Sample_Type",
    "Status"
  )


  missing_columns <- setdiff(
    required_metadata_columns,
    names(metadata)
  )


  if (length(missing_columns) > 0) {

    stop(
      paste0(
        label,
        " metadata is missing:\n",
        paste(
          missing_columns,
          collapse = "\n"
        )
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Taxonomy
  # ---------------------------------------------------------------------------

  taxonomy <- fix_taxonomy(
    taxonomy
  )


  missing_tax <- setdiff(
    abundance$Genome,
    taxonomy$Genome
  )


  if (length(missing_tax) > 0) {

    stop(
      label,
      ": ",
      length(missing_tax),
      " genomes do not have taxonomy."
    )
  }


  # ---------------------------------------------------------------------------
  # Convert abundance columns to numeric
  # ---------------------------------------------------------------------------

  numeric_columns <- setdiff(
    names(abundance),
    "Genome"
  )


  abundance[numeric_columns] <- lapply(

    abundance[numeric_columns],

    function(x) {

      y <- suppressWarnings(
        as.numeric(x)
      )


      if (any(!is.finite(y))) {

        stop(
          label,
          ": non-finite abundance values found."
        )
      }


      if (any(y < 0)) {

        stop(
          label,
          ": negative abundance values found."
        )
      }


      y
    }
  )


  # ---------------------------------------------------------------------------
  # Aggregate genome TPM to viral-family abundance
  # ---------------------------------------------------------------------------

  tax_subset <- taxonomy[
    ,
    c(
      "Genome",
      "Family",
      "Phylum"
    )
  ]


  merged <- dplyr::left_join(
    abundance,
    tax_subset,
    by = "Genome"
  )


  family_df <- merged %>%

    dplyr::select(
      -Genome,
      -Phylum
    ) %>%

    dplyr::group_by(
      Family
    ) %>%

    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )


  family_mat <- as.matrix(
    family_df[
      ,
      setdiff(
        names(family_df),
        "Family"
      ),
      drop = FALSE
    ]
  )


  rownames(
    family_mat
  ) <- family_df$Family


  storage.mode(
    family_mat
  ) <- "double"


  family_mat[
    !is.finite(family_mat)
  ] <- 0


  # Remove fully unclassified root category
  family_mat <- family_mat[
    rownames(family_mat) != "Unclassified_Root",
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # Remove zero-sum samples and normalize to CPM
  # ---------------------------------------------------------------------------

  sample_totals <- colSums(
    family_mat
  )


  zero_sum_samples <- names(
    sample_totals[
      sample_totals <= 0
    ]
  )


  keep_samples <- names(
    sample_totals[
      sample_totals > 0
    ]
  )


  if (length(zero_sum_samples) > 0) {

    cat(
      label,
      ": excluding",
      length(zero_sum_samples),
      "zero-sum samples.\n"
    )
  }


  family_mat <- family_mat[
    ,
    keep_samples,
    drop = FALSE
  ]


  family_cpm <- sweep(
    family_mat,
    2,
    colSums(family_mat),
    "/"
  ) * 1e6


  # ---------------------------------------------------------------------------
  # Clean metadata
  # ---------------------------------------------------------------------------

  metadata <- metadata %>%

    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        as.character
      )
    )


  if (!"Location" %in% names(metadata)) {

    metadata$Location <- NA_character_
  }


  metadata <- metadata %>%

    dplyr::mutate(

      Sample = trimws(
        Sample
      ),

      Project = trimws(
        Project
      ),

      Ocean_Basin = trimws(
        Ocean_Basin
      ),

      Coral_Host_Species = trimws(
        Coral_Host_Species
      ),

      Sample_Type = trimws(
        Sample_Type
      ),

      Location = trimws(
        Location
      ),

      Status_Normalized = normalise_status(
        Status
      ),

      Analysis_Project = make_analysis_project(
        Project,
        label
      )
    )


  # ---------------------------------------------------------------------------
  # Match abundance and metadata samples
  # ---------------------------------------------------------------------------

  matched_samples <- intersect(
    colnames(family_cpm),
    metadata$Sample
  )


  if (length(matched_samples) == 0) {

    stop(
      label,
      ": no abundance samples match metadata."
    )
  }


  family_cpm <- family_cpm[
    ,
    matched_samples,
    drop = FALSE
  ]


  metadata <- metadata[
    match(
      matched_samples,
      metadata$Sample
    ),
    ,
    drop = FALSE
  ]


  rownames(
    metadata
  ) <- metadata$Sample


  # ---------------------------------------------------------------------------
  # Export raw Project -> Analysis_Project mapping
  # ---------------------------------------------------------------------------

  project_mapping <- metadata %>%

    dplyr::distinct(
      Project,
      Analysis_Project
    ) %>%

    dplyr::arrange(
      Analysis_Project,
      Project
    )


  write.csv(
    project_mapping,
    file.path(
      OUTPUT_DIR,
      paste0(
        label,
        "_Project_to_AnalysisProject.csv"
      )
    ),
    row.names = FALSE
  )


  cat(
    label,
    ": matched positive-abundance samples =",
    ncol(family_cpm),
    "\n"
  )


  list(
    label = label,
    family_cpm = family_cpm,
    metadata = metadata,
    taxonomy = taxonomy,
    zero_sum_samples = zero_sum_samples
  )
}


# ==============================================================================
# 8. Helper functions
# ==============================================================================

make_stratum_id <- function(
    analysis_project,
    ocean_basin,
    species,
    sample_type
) {

  paste(
    analysis_project,
    ocean_basin,
    species,
    sample_type,
    sep = " | "
  )
}


safe_kw_p <- function(
    values,
    class_factor
) {

  if (
    length(unique(class_factor)) != 2L |
      length(values) != length(class_factor)
  ) {

    return(
      NA_real_
    )
  }


  tryCatch(
    kruskal.test(
      values ~ class_factor
    )$p.value,
    error = function(e) NA_real_
  )
}


# ==============================================================================
# 9. Run one comparison
# ==============================================================================

run_one_comparison <- function(
    virome,
    condition,
    stratum_row
) {

  md <- virome$metadata

  # DSS only: exclude nutrient-enriched variants before selecting model samples.
  if (identical(condition, "DSS")) {
    md <- md[trimws(as.character(md$Status)) %in%
      c("Healthy/Normal", "Disease(Dark Spot Syndrome)"), , drop = FALSE]
  }


  # ---------------------------------------------------------------------------
  # Final filtering rule:
  #
  # same Analysis_Project
  # + same Ocean_Basin
  # + same Coral_Host_Species
  # + same Sample_Type
  #
  # Location is NOT used here.
  # ---------------------------------------------------------------------------

  sub_md <- md %>%

    dplyr::filter(

      Analysis_Project ==
        stratum_row$Analysis_Project,

      Ocean_Basin ==
        stratum_row$Ocean_Basin,

      Coral_Host_Species ==
        stratum_row$Coral_Host_Species,

      Sample_Type ==
        stratum_row$Sample_Type,

      Status_Normalized %in%
        c(
          "Healthy",
          condition
        )
    )


  ids <- intersect(
    sub_md$Sample,
    colnames(
      virome$family_cpm
    )
  )


  sub_md <- sub_md[
    match(
      ids,
      sub_md$Sample
    ),
    ,
    drop = FALSE
  ]


  sub_mat <- virome$family_cpm[
    ,
    ids,
    drop = FALSE
  ]


  class_factor <- factor(
    sub_md$Status_Normalized,
    levels = c(
      "Healthy",
      condition
    )
  )


  counts <- table(
    class_factor
  )


  # ---------------------------------------------------------------------------
  # HARD PROJECT QC
  # ---------------------------------------------------------------------------

  if (
    dplyr::n_distinct(
      sub_md$Analysis_Project
    ) != 1L
  ) {

    stop(
      "Internal error: multiple Analysis_Project values entered one comparison."
    )
  }


  current_project <- unique(
    sub_md$Analysis_Project
  )


  raw_projects <- sort(
    unique(
      sub_md$Project
    )
  )


  # Ordinary Analysis_Project:
  # exactly one raw Project accession must be present.
  if (
    current_project != "Avila2021" &&
      length(raw_projects) != 1L
  ) {

    stop(
      paste0(
        "Cross-project mixing detected in comparison:\n",
        current_project,
        "\nRaw Projects: ",
        paste(
          raw_projects,
          collapse = ", "
        )
      )
    )
  }


  # Avila2021:
  # multiple raw PRJNA accessions are allowed, but ALL must belong to the
  # explicitly approved PRJNA405xxx list.
  if (
    current_project == "Avila2021" &&
      any(
        !raw_projects %in%
          AVILA_PROJECTS
      )
  ) {

    stop(
      paste0(
        "Unexpected Project accession inside Avila2021:\n",
        paste(
          raw_projects[
            !raw_projects %in%
              AVILA_PROJECTS
          ],
          collapse = ", "
        )
      )
    )
  }


  stratum_id <- make_stratum_id(
    stratum_row$Analysis_Project,
    stratum_row$Ocean_Basin,
    stratum_row$Coral_Host_Species,
    stratum_row$Sample_Type
  )


  # Keep Location only for traceability
  locations <- sort(
    unique(
      sub_md$Location[
        !is.na(sub_md$Location) &
          sub_md$Location != ""
      ]
    )
  )


  location_text <- if (
    length(locations) > 0
  ) {

    paste(
      locations,
      collapse = "; "
    )

  } else {

    ""
  }


  raw_project_text <- paste(
    raw_projects,
    collapse = "; "
  )


  # ---------------------------------------------------------------------------
  # Family-level Kruskal-Wallis tests
  # ---------------------------------------------------------------------------

  feature_rows <- lapply(

    seq_len(
      nrow(sub_mat)
    ),

    function(i) {

      vals <- as.numeric(
        sub_mat[i, ]
      )


      data.frame(

        Virome = virome$label,

        Condition = condition,

        Stratum_ID = stratum_id,

        Analysis_Project =
          stratum_row$Analysis_Project,

        Raw_Projects =
          raw_project_text,

        Ocean_Basin =
          stratum_row$Ocean_Basin,

        Coral_Host_Species =
          stratum_row$Coral_Host_Species,

        Sample_Type =
          stratum_row$Sample_Type,

        Locations =
          location_text,

        Family =
          rownames(sub_mat)[i],

        Healthy_n =
          unname(
            counts[["Healthy"]]
          ),

        Condition_n =
          unname(
            counts[[condition]]
          ),

        Healthy_mean_CPM =
          mean(
            vals[
              class_factor ==
                "Healthy"
            ],
            na.rm = TRUE
          ),

        Condition_mean_CPM =
          mean(
            vals[
              class_factor ==
                condition
            ],
            na.rm = TRUE
          ),

        KW_P =
          safe_kw_p(
            vals,
            class_factor
          ),

        stringsAsFactors = FALSE
      )
    }
  )


  all_tests <- dplyr::bind_rows(
    feature_rows
  )


  all_tests$KW_P_BH_within_stratum <- p.adjust(
    all_tests$KW_P,
    method = "BH"
  )


  # ---------------------------------------------------------------------------
  # LEfSe
  # ---------------------------------------------------------------------------

  selected <- data.frame()

  warning_text <- character()


  if (
    all(counts >= MIN_N) &&
      ncol(sub_mat) >= 4L
  ) {

    se <- SummarizedExperiment::SummarizedExperiment(

      assays = list(
        CPM = sub_mat
      ),

      colData =
        S4Vectors::DataFrame(
          Class = class_factor
        )
    )


    set.seed(123)


    fit <- tryCatch(

      withCallingHandlers(

        lefser::lefser(

          se,

          classCol =
            "Class",

          kruskal.threshold =
            KW_THRESHOLD,

          lda.threshold =
            LDA_THRESHOLD,

          method =
            "none",

          checkAbundances =
            TRUE,

          trim.names =
            FALSE
        ),

        warning = function(w) {

          warning_text <<- c(
            warning_text,
            conditionMessage(w)
          )

          invokeRestart(
            "muffleWarning"
          )
        }
      ),

      error = function(e) {

        warning_text <<- c(
          warning_text,
          paste0(
            "ERROR: ",
            conditionMessage(e)
          )
        )

        NULL
      }
    )


    if (
      !is.null(fit) &&
        nrow(fit) > 0
    ) {

      selected <- as.data.frame(
        fit
      )


      names(selected)[
        names(selected) ==
          "features"
      ] <- "Family"


      names(selected)[
        names(selected) ==
          "scores"
      ] <- "LDA_Signed"


      selected <- selected %>%

        dplyr::mutate(

          LDA_Abs =
            abs(
              LDA_Signed
            ),

          Enriched_In =
            dplyr::if_else(
              LDA_Signed >= 0,
              condition,
              "Healthy"
            ),

          Virome =
            virome$label,

          Condition =
            condition,

          Stratum_ID =
            stratum_id,

          Analysis_Project =
            stratum_row$Analysis_Project,

          Raw_Projects =
            raw_project_text,

          Ocean_Basin =
            stratum_row$Ocean_Basin,

          Coral_Host_Species =
            stratum_row$Coral_Host_Species,

          Sample_Type =
            stratum_row$Sample_Type,

          Locations =
            location_text,

          Healthy_n =
            unname(
              counts[["Healthy"]]
            ),

          Condition_n =
            unname(
              counts[[condition]]
            ),

          LDA_Threshold =
            LDA_THRESHOLD,

          KW_Threshold =
            KW_THRESHOLD
        ) %>%

        dplyr::left_join(

          all_tests %>%

            dplyr::select(
              Stratum_ID,
              Family,
              Healthy_mean_CPM,
              Condition_mean_CPM,
              KW_P,
              KW_P_BH_within_stratum
            ),

          by = c(
            "Stratum_ID",
            "Family"
          )
        )


      selected$Direction_Matches_Group_Mean <- with(

        selected,

        (
          Enriched_In ==
            "Healthy" &
            Healthy_mean_CPM >=
              Condition_mean_CPM
        ) |
          (
            Enriched_In ==
              condition &
              Condition_mean_CPM >=
                Healthy_mean_CPM
          )
      )
    }
  }


  # ---------------------------------------------------------------------------
  # Add LEfSe fields to all tested families
  # ---------------------------------------------------------------------------

  all_tests <- all_tests %>%

    dplyr::mutate(

      LDA_Signed =
        NA_real_,

      LDA_Abs =
        NA_real_,

      Enriched_In =
        NA_character_,

      Retained =
        FALSE,

      LDA_Threshold =
        LDA_THRESHOLD,

      KW_Threshold =
        KW_THRESHOLD,

      Fit_Warnings =
        paste(
          unique(
            warning_text
          ),
          collapse = " | "
        )
    )


  if (
    nrow(selected) > 0
  ) {

    all_tests <- all_tests %>%

      dplyr::select(
        -LDA_Signed,
        -LDA_Abs,
        -Enriched_In
      ) %>%

      dplyr::left_join(

        selected %>%

          dplyr::select(
            Stratum_ID,
            Family,
            LDA_Signed,
            LDA_Abs,
            Enriched_In
          ),

        by = c(
          "Stratum_ID",
          "Family"
        )
      ) %>%

      dplyr::mutate(
        Retained =
          is.finite(
            LDA_Abs
          ) &
          LDA_Abs >=
            LDA_THRESHOLD
      )
  }


  # ---------------------------------------------------------------------------
  # Exact sample membership for QC
  # ---------------------------------------------------------------------------

  membership <- sub_md %>%

    dplyr::transmute(

      Virome =
        virome$label,

      Condition =
        condition,

      Stratum_ID =
        stratum_id,

      Analysis_Project,

      Raw_Project =
        Project,

      Ocean_Basin,

      Coral_Host_Species,

      Sample_Type,

      Location,

      Status_Original = Status,

      Status_Normalized,

      Sample
    )


  list(

    all_tests =
      all_tests,

    selected =
      selected,

    membership =
      membership,

    warning_text =
      warning_text
  )
}


# ==============================================================================
# 10. Identify eligible project-aware comparisons
# ==============================================================================

analyse_virome <- function(
    virome,
    conditions
) {

  md <- virome$metadata %>%

    dplyr::filter(

      !is.na(
        Status_Normalized
      ),

      dplyr::if_all(
        dplyr::all_of(
          STRATA_COLUMNS
        ),
        ~!is.na(.x) &
          .x != ""
      )
    )


  strata_results <-
    list()

  selected_results <-
    list()

  all_tests_results <-
    list()

  membership_results <-
    list()


  for (
    condition in conditions
  ) {

    # Use the same DSS-only rule for eligibility counts and actual model input.
    condition_md <- md
    if (identical(condition, "DSS")) {
      condition_md <- condition_md[trimws(as.character(condition_md$Status)) %in%
        c("Healthy/Normal", "Disease(Dark Spot Syndrome)"), , drop = FALSE]
    }

    candidate <- condition_md %>%

      dplyr::filter(
        Status_Normalized %in%
          c(
            "Healthy",
            condition
          )
      ) %>%

      dplyr::count(
        dplyr::across(
          dplyr::all_of(
            STRATA_COLUMNS
          )
        ),
        Status_Normalized,
        name = "n"
      ) %>%

      tidyr::pivot_wider(
        names_from =
          Status_Normalized,
        values_from =
          n,
        values_fill =
          0
      )


    if (
      !"Healthy" %in%
        names(candidate) ||
        !condition %in%
          names(candidate)
    ) {

      next
    }


    candidate <- candidate %>%

      dplyr::mutate(

        Virome =
          virome$label,

        Condition =
          condition,

        Healthy_n =
          .data[["Healthy"]],

        Condition_n =
          .data[[condition]],

        Eligible =
          Healthy_n >=
            MIN_N &
            Condition_n >=
              MIN_N
      )


    strata_results[[condition]] <- candidate


    eligible <- candidate %>%

      dplyr::filter(
        Eligible
      )


    if (
      nrow(eligible) == 0
    ) {

      next
    }


    for (
      i in seq_len(
        nrow(eligible)
      )
    ) {

      result <- run_one_comparison(
        virome,
        condition,
        eligible[i, ]
      )


      all_tests_results[[length(all_tests_results) + 1L]] <- result$all_tests


      membership_results[[length(membership_results) + 1L]] <- result$membership


      if (
        nrow(result$selected) > 0
      ) {

        selected_results[[length(selected_results) + 1L]] <- result$selected
      }
    }
  }


  list(

    strata =
      dplyr::bind_rows(
        strata_results
      ),

    all_tests =
      dplyr::bind_rows(
        all_tests_results
      ),

    selected =
      dplyr::bind_rows(
        selected_results
      ),

    membership =
      dplyr::bind_rows(
        membership_results
      )
  )
}


# ==============================================================================
# 11. Read DNA and RNA
# ==============================================================================

viromes <- list(

  DNA = read_virome(
    label =
      "DNA",
    abundance_file =
      "24083_metagenome_TPM.tsv",
    taxonomy_file =
      "dna_virus_taxonomy.tsv",
    metadata_file =
      "sample_info_metagenome.tsv"
  ),

  RNA = read_virome(
    label =
      "RNA",
    abundance_file =
      "rna_virus_tpm.tsv",
    taxonomy_file =
      "rna_virus_taxonomy.tsv",
    metadata_file =
      "sample_info_metatranscriptome.tsv"
  )
)


# ==============================================================================
# 12. Run analyses
# ==============================================================================

results <- list(

  DNA = analyse_virome(
    viromes$DNA,
    DNA_CONDITIONS
  ),

  RNA = analyse_virome(
    viromes$RNA,
    RNA_CONDITIONS
  )
)


all_strata <- dplyr::bind_rows(
  results$DNA$strata,
  results$RNA$strata
)


all_tests <- dplyr::bind_rows(
  results$DNA$all_tests,
  results$RNA$all_tests
)


selected <- dplyr::bind_rows(
  results$DNA$selected,
  results$RNA$selected
)


membership <- dplyr::bind_rows(
  results$DNA$membership,
  results$RNA$membership
)


if (
  nrow(all_tests) == 0
) {

  stop(
    "No eligible Healthy-vs-condition comparisons were found."
  )
}


# ==============================================================================
# 13. Taxonomy map
# ==============================================================================

tax_map <- dplyr::bind_rows(

  viromes$DNA$taxonomy %>%

    dplyr::transmute(
      Virome =
        "DNA",
      Family,
      Phylum
    ),

  viromes$RNA$taxonomy %>%

    dplyr::transmute(
      Virome =
        "RNA",
      Family,
      Phylum
    )
) %>%

  dplyr::distinct() %>%

  dplyr::group_by(
    Virome,
    Family
  ) %>%

  dplyr::summarise(
    Phylum =
      dplyr::first(
        Phylum
      ),
    .groups =
      "drop"
  )


# ==============================================================================
# 14. Comparison IDs
# ==============================================================================

comparison_key <- all_strata %>%

  dplyr::filter(
    Eligible
  ) %>%

  dplyr::distinct(

    Virome,
    Condition,
    Analysis_Project,
    Ocean_Basin,
    Coral_Host_Species,
    Sample_Type,
    Healthy_n,
    Condition_n
  ) %>%

  dplyr::mutate(

    Stratum_ID =
      make_stratum_id(
        Analysis_Project,
        Ocean_Basin,
        Coral_Host_Species,
        Sample_Type
      ),

    Condition_Order =
      match(
        Condition,
        c(
          "SCTLD",
          "Bleached",
          "Heat_Stress",
          "DSS"
        )
      )
  ) %>%

  dplyr::arrange(

    Virome,
    Condition_Order,
    Analysis_Project,
    Ocean_Basin,
    Coral_Host_Species,
    Sample_Type
  ) %>%

  dplyr::group_by(
    Virome
  ) %>%

  dplyr::mutate(

    Comparison_ID =
      paste0(

        dplyr::if_else(
          Virome ==
            "DNA",
          "D",
          "R"
        ),

        sprintf(
          "%02d",
          dplyr::row_number()
        )
      )
  ) %>%

  dplyr::ungroup() %>%

  dplyr::select(
    -Condition_Order
  )


# ==============================================================================
# 15. Comparison traceability
# ==============================================================================

comparison_trace <- membership %>%

  dplyr::group_by(

    Virome,
    Condition,
    Stratum_ID,
    Analysis_Project,
    Ocean_Basin,
    Coral_Host_Species,
    Sample_Type
  ) %>%

  dplyr::summarise(

    Raw_Projects =
      paste(
        sort(
          unique(
            Raw_Project
          )
        ),
        collapse = "; "
      ),

    N_Raw_Projects =
      dplyr::n_distinct(
        Raw_Project
      ),

    Locations =
      paste(
        sort(
          unique(
            Location[
              !is.na(Location) &
                Location != ""
            ]
          )
        ),
        collapse = "; "
      ),

    N_Locations =
      dplyr::n_distinct(
        Location[
          !is.na(Location) &
            Location != ""
        ]
      ),

    .groups =
      "drop"
  )


comparison_key <- comparison_key %>%

  dplyr::left_join(

    comparison_trace,

    by = c(
      "Virome",
      "Condition",
      "Stratum_ID",
      "Analysis_Project",
      "Ocean_Basin",
      "Coral_Host_Species",
      "Sample_Type"
    )
  )


# ==============================================================================
# 16. Project integrity QC
# ==============================================================================

project_integrity_qc <- comparison_key %>%

  dplyr::mutate(

    Project_integrity_OK =
      dplyr::case_when(

        Analysis_Project ==
          "Avila2021" &
          N_Raw_Projects >=
            1 ~ TRUE,

        Analysis_Project !=
          "Avila2021" &
          N_Raw_Projects ==
            1 ~ TRUE,

        TRUE ~ FALSE
      )
  )


if (
  any(
    !project_integrity_qc$Project_integrity_OK,
    na.rm = TRUE
  )
) {

  bad <- project_integrity_qc %>%

    dplyr::filter(
      !Project_integrity_OK
    )


  write.csv(

    bad,

    file.path(
      OUTPUT_DIR,
      "ERROR_cross_project_mixing_detected.csv"
    ),

    row.names = FALSE
  )


  stop(
    paste0(
      "Project integrity QC failed.\n",
      "See ERROR_cross_project_mixing_detected.csv"
    )
  )
}


# ==============================================================================
# 17. Add taxonomy and Comparison_ID
# ==============================================================================

all_tests <- all_tests %>%

  dplyr::left_join(

    tax_map,

    by = c(
      "Virome",
      "Family"
    )
  ) %>%

  dplyr::mutate(

    Phylum =
      dplyr::if_else(
        is.na(Phylum) |
          Phylum == "",
        "Unclassified",
        Phylum
      )
  ) %>%

  dplyr::left_join(

    comparison_key %>%

      dplyr::select(
        Virome,
        Condition,
        Stratum_ID,
        Comparison_ID
      ),

    by = c(
      "Virome",
      "Condition",
      "Stratum_ID"
    )
  )


if (
  nrow(selected) > 0
) {

  selected <- selected %>%

    dplyr::left_join(

      tax_map,

      by = c(
        "Virome",
        "Family"
      )
    ) %>%

    dplyr::mutate(

      Phylum =
        dplyr::if_else(
          is.na(Phylum) |
            Phylum == "",
          "Unclassified",
          Phylum
        )
    ) %>%

    dplyr::left_join(

      comparison_key %>%

        dplyr::select(
          Virome,
          Condition,
          Stratum_ID,
          Comparison_ID
        ),

      by = c(
        "Virome",
        "Condition",
        "Stratum_ID"
      )
    )
}


# ==============================================================================
# 18. Comparison summary
# ==============================================================================

comparison_summary <- all_tests %>%

  dplyr::group_by(

    Virome,
    Condition,
    Comparison_ID,
    Stratum_ID,
    Analysis_Project,
    Raw_Projects,
    Ocean_Basin,
    Coral_Host_Species,
    Sample_Type,
    Locations,
    Healthy_n,
    Condition_n
  ) %>%

  dplyr::summarise(

    Families_tested =
      dplyr::n(),

    Families_KW_P_lt_0_05 =
      sum(
        KW_P <
          KW_THRESHOLD,
        na.rm = TRUE
      ),

    Families_KW_BH_lt_0_05 =
      sum(
        KW_P_BH_within_stratum <
          0.05,
        na.rm = TRUE
      ),

    Families_LDA_abs_ge_2 =
      sum(
        Retained,
        na.rm = TRUE
      ),

    Fit_Warnings =
      paste(
        unique(
          Fit_Warnings[
            !is.na(Fit_Warnings) &
              nzchar(Fit_Warnings)
          ]
        ),
        collapse = " | "
      ),

    .groups =
      "drop"
  )


# ==============================================================================
# 19. Export all QC / source tables
# ==============================================================================

write.csv(
  project_integrity_qc,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_project_integrity_QC.csv"
  ),
  row.names = FALSE
)


write.csv(
  comparison_key,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_comparison_key.csv"
  ),
  row.names = FALSE
)


write.csv(
  membership,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_sample_membership.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_strata,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_all_strata.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_tests,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_all_tested.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  selected,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_selected.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  comparison_summary,
  file.path(
    OUTPUT_DIR,
    "ProjectAware_LEfSe_comparison_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ==============================================================================
# 20-25. FINAL SQUARE-CELL SIDE-BY-SIDE PUBLICATION FIGURE
# ==============================================================================
#
# IMPORTANT GEOMETRY
# ------------------------------------------------------------------------------
# The DNA and RNA panels contain very different numbers of viral families.
# Therefore three goals cannot all be achieved by simply stretching both
# matrices to the same height:
#   (i) identical column width,
#   (ii) square cells,
#   (iii) every RNA data row stretched to fill the full DNA matrix height.
#
# This version prioritizes:
#   * identical DNA/RNA column width;
#   * square bubble-matrix cells;
#   * identical physical height of Sampling region / Compartment / Status rows;
#   * DNA and RNA panels with the same total external height.
#
# The RNA data matrix therefore keeps square cells at the TOP of its panel and
# receives white padding below it to match the total DNA-panel height. This is
# visually cleaner than making the four RNA-family rows extremely tall.
#
# The final assembly below preserves these square cells and top alignment.
# ==============================================================================


# ==============================================================================
# 20. FIGURE SETTINGS
# ==============================================================================

# LEfSe analysis was performed at |LDA| >= 2.
# Final manuscript display retains markers with |LDA| >= 2.5.
LDA_THRESHOLD_PLOT <- 2.5

BASE_FONT <- "Arial"

# Physical geometry (millimetres)
CELL_MM          <- 5.6   # one comparison column = one cell width
LEFT_LABEL_MM    <- 45    # phylum strip + viral-family names
TITLE_MM         <- 8
X_LABEL_MM       <- 31    # room for rotated coral-taxon labels
PANEL_GAP_MM     <- 6
LEGEND_MM        <- 48
FIG_MARGIN_MM    <- 8

# Narrow phylum strip inside the left-label panel
PHYLUM_STRIP_FRAC <- 0.085


# ==============================================================================
# 21. FILTER MARKERS DISPLAYED IN THE FINAL FIGURE
# ==============================================================================

selected_plot <- selected %>%

  dplyr::filter(
    is.finite(LDA_Abs),
    LDA_Abs >= LDA_THRESHOLD_PLOT
  )

if (
  nrow(selected_plot) == 0
) {

  stop(
    paste0(
      "No viral families remain after |LDA| >= ",
      LDA_THRESHOLD_PLOT,
      ".\nLower LDA_THRESHOLD_PLOT if you want to display more LEfSe markers."
    )
  )
}


write.csv(
  selected_plot,
  file.path(
    OUTPUT_DIR,
    paste0(
      "ProjectAware_LEfSe_selected_LDA",
      LDA_THRESHOLD_PLOT,
      "_display.csv"
    )
  ),
  row.names = FALSE,
  na = ""
)


# ==============================================================================
# 22. GLOBAL ORDERS, COLOURS AND SHARED LDA SCALE
# ==============================================================================

CONDITION_LEVELS <- c(
  "SCTLD",
  "Bleached",
  "Heat_Stress",
  "DSS"
)

CONDITION_LABELS <- c(
  "SCTLD"       = "SCTLD",
  "Bleached"    = "Bleached",
  "Heat_Stress" = "Heat stress",
  "DSS"         = "DSS"
)


# ------------------------------------------------------------------------------
# Bubble colours
# ------------------------------------------------------------------------------

ENRICH_PALETTE <- c(
  "Healthy"     = "#3C9D63",
  "SCTLD"       = "#B2182B",
  "Bleached"    = "#E5A33D",
  "Heat_Stress" = "#6259A8",
  "DSS"         = "#9A4A87"
)

ENRICH_LABELS <- c(
  "Healthy"     = "Healthy",
  "SCTLD"       = "SCTLD",
  "Bleached"    = "Bleached",
  "Heat_Stress" = "Heat stress",
  "DSS"         = "DSS"
)


# Status-header fills: intentionally lighter than bubble colours
STATUS_HEADER_PALETTE <- c(
  "SCTLD"       = "#D8A3AC",
  "Bleached"    = "#F2C277",
  "Heat_Stress" = "#9895CE",
  "DSS"         = "#C896B7"
)


# ------------------------------------------------------------------------------
# Identical LDA-size scale for DNA and RNA
# ------------------------------------------------------------------------------

GLOBAL_LDA_MAX <- ceiling(
  max(
    selected_plot$LDA_Abs,
    na.rm = TRUE
  )
)

GLOBAL_LDA_MAX <- max(
  GLOBAL_LDA_MAX,
  LDA_THRESHOLD_PLOT
)

LDA_BREAKS <- pretty(
  c(
    LDA_THRESHOLD_PLOT,
    GLOBAL_LDA_MAX
  ),
  n = 4
)

LDA_BREAKS <- LDA_BREAKS[
  LDA_BREAKS >= LDA_THRESHOLD_PLOT &
    LDA_BREAKS <= GLOBAL_LDA_MAX
]

if (
  !LDA_THRESHOLD_PLOT %in% LDA_BREAKS
) {

  LDA_BREAKS <- sort(
    unique(
      c(
        LDA_THRESHOLD_PLOT,
        LDA_BREAKS
      )
    )
  )
}


# ------------------------------------------------------------------------------
# Palette helper
# ------------------------------------------------------------------------------

complete_named_palette <- function(
    values,
    preferred,
    fallback_palette = "Dynamic"
) {

  values <- unique(
    as.character(
      values
    )
  )

  values <- values[
    !is.na(values) &
      values != ""
  ]

  missing_values <- setdiff(
    values,
    names(preferred)
  )

  out <- preferred[
    intersect(
      names(preferred),
      values
    )
  ]

  if (
    length(missing_values) > 0
  ) {

    extra <- grDevices::hcl.colors(
      length(missing_values),
      palette = fallback_palette
    )

    names(extra) <- missing_values

    out <- c(
      out,
      extra
    )
  }

  out[
    values
  ]
}


# ------------------------------------------------------------------------------
# Virus phylum colours shared across DNA and RNA
# ------------------------------------------------------------------------------

PHYLUM_PREFERRED <- c(
  "Uroviricota"           = "#63B7A6",
  "Nucleocytoviricota"    = "#839CCB",
  "Phixviricota"          = "#9CCB5A",
  "Cressdnaviricota"      = "#EE8766",
  "Preplasmiviricota"     = "#D884B7",
  "Pisuviricota"          = "#E9C55D",
  "Kitrinoviricota"       = "#65B8B5",
  "Lenarviricota"         = "#A9C96B",
  "Negarnaviricota"       = "#72A8D0",
  "Duplornaviricota"      = "#B57C5C",
  "Hofneiviricota"        = "#D8B884",
  "Unclassified"          = "#B8B8B8"
)


ALL_PHYLA <- selected_plot %>%

  dplyr::mutate(
    Phylum = dplyr::if_else(
      is.na(Phylum) |
        Phylum == "",
      "Unclassified",
      Phylum
    )
  ) %>%

  dplyr::distinct(
    Phylum
  ) %>%

  dplyr::pull(
    Phylum
  )


PHYLUM_PALETTE <- complete_named_palette(
  ALL_PHYLA,
  PHYLUM_PREFERRED,
  "Set 3"
)


# ------------------------------------------------------------------------------
# Sampling-region colours
# ------------------------------------------------------------------------------

REGION_PREFERRED <- c(
  "Atlantic Ocean"    = "#4E79A7",
  "Caribbean Sea"     = "#59A14F",
  "Pacific Ocean"     = "#F28E2B",
  "Red Sea"           = "#E15759",
  "Arabian Gulf"      = "#B07AA1",
  "Mediterranean Sea" = "#76B7B2",
  "Indian Ocean"      = "#EDC948"
)


ALL_REGIONS <- comparison_key %>%

  dplyr::distinct(
    Ocean_Basin
  ) %>%

  dplyr::pull(
    Ocean_Basin
  )


REGION_PALETTE <- complete_named_palette(
  ALL_REGIONS,
  REGION_PREFERRED,
  "Dark 3"
)


# ------------------------------------------------------------------------------
# Compartment colours
# ------------------------------------------------------------------------------

COMPARTMENT_PREFERRED <- c(
  "Coral tissue"          = "#F1D1C1",
  "Coral fragment"        = "#C9B3D9",
  "Coral mucus"           = "#A6D9D5",
  "Coral skeleton"        = "#B8D9F0",
  "Coral endolithic band" = "#F1D29B",
  "Coral nubbins"         = "#C7E0CA",
  "Mucus-tissue slurry"   = "#E8C1D8"
)


ALL_COMPARTMENTS <- comparison_key %>%

  dplyr::distinct(
    Sample_Type
  ) %>%

  dplyr::pull(
    Sample_Type
  )


COMPARTMENT_PALETTE <- complete_named_palette(
  ALL_COMPARTMENTS,
  COMPARTMENT_PREFERRED,
  "Pastel 1"
)


# One fill scale for the two annotation rows inside the matrix plot
REGION_CONTEXT_KEYS <- paste0(
  "Region::",
  names(
    REGION_PALETTE
  )
)

COMPARTMENT_CONTEXT_KEYS <- paste0(
  "Compartment::",
  names(
    COMPARTMENT_PALETTE
  )
)

CONTEXT_COLOURS <- c(
  stats::setNames(
    unname(
      REGION_PALETTE
    ),
    REGION_CONTEXT_KEYS
  ),
  stats::setNames(
    unname(
      COMPARTMENT_PALETTE
    ),
    COMPARTMENT_CONTEXT_KEYS
  )
)


# ==============================================================================
# 23. ORDER PHYLA AND FAMILIES
# ==============================================================================

make_phylum_order <- function(
    df,
    virome_name
) {

  phylum_counts <- df %>%

    dplyr::distinct(
      Phylum,
      Family
    ) %>%

    dplyr::count(
      Phylum,
      name = "N_Families"
    )


  if (
    virome_name == "DNA"
  ) {

    phylum_counts <- phylum_counts %>%

      dplyr::mutate(
        Priority = dplyr::case_when(
          Phylum == "Uroviricota"        ~ 1L,
          Phylum == "Nucleocytoviricota" ~ 2L,
          TRUE                           ~ 3L
        )
      ) %>%

      dplyr::arrange(
        Priority,
        dplyr::desc(
          N_Families
        ),
        Phylum
      )

  } else {

    phylum_counts <- phylum_counts %>%

      dplyr::mutate(
        Priority = 1L
      ) %>%

      dplyr::arrange(
        dplyr::desc(
          N_Families
        ),
        Phylum
      )
  }


  phylum_counts %>%

    dplyr::mutate(
      Phylum_Order = dplyr::row_number()
    )
}


# Within each phylum:
#   1. formally named families ending in "-idae"
#   2. other named labels
#   3. "Unclass..." labels
#
# This implements the requested "known family -> unclassified" order.
make_family_map <- function(
    df,
    virome_name
) {

  phylum_order <- make_phylum_order(
    df,
    virome_name
  )


  family_map <- df %>%

    dplyr::distinct(
      Family,
      Phylum
    ) %>%

    dplyr::left_join(
      phylum_order %>%

        dplyr::select(
          Phylum,
          N_Families,
          Phylum_Order
        ),
      by = "Phylum"
    ) %>%

    dplyr::mutate(

      Family_Class = dplyr::case_when(

        grepl(
          "idae$",
          Family,
          ignore.case = TRUE
        ) ~ 1L,

        grepl(
          "^Unclass",
          Family,
          ignore.case = TRUE
        ) ~ 3L,

        TRUE ~ 2L
      )
    ) %>%

    dplyr::arrange(
      Phylum_Order,
      Family_Class,
      Family
    ) %>%

    dplyr::mutate(
      Display_Order = dplyr::row_number(),
      Y_Pos = dplyr::n() -
        Display_Order +
        1L
    )


  list(
    family_map = family_map,
    phylum_order = phylum_order
  )
}


# ==============================================================================
# 24. BUILD ONE PANEL AS:
#
#     LEFT: phylum strip + family labels
#     RIGHT: status / region / compartment / square bubble matrix
#
# No facets are used. All x columns are numeric and evenly spaced, so every
# comparison column has the same physical width.
# ==============================================================================

build_virome_panel <- function(
    virome_name
) {

  df <- selected_plot %>%

    dplyr::filter(
      Virome == virome_name
    ) %>%

    dplyr::mutate(
      Phylum = dplyr::if_else(
        is.na(Phylum) |
          Phylum == "",
        "Unclassified",
        Phylum
      )
    )


  if (
    nrow(df) == 0
  ) {

    stop(
      "No displayed ",
      virome_name,
      " markers."
    )
  }


  # ---------------------------------------------------------------------------
  # Comparison-column order
  # ---------------------------------------------------------------------------

  shown_comparisons <- unique(
    df$Comparison_ID
  )


  lookup <- comparison_key %>%

    dplyr::filter(
      Virome == virome_name,
      Comparison_ID %in% shown_comparisons
    ) %>%

    dplyr::distinct(
      Comparison_ID,
      Condition,
      Analysis_Project,
      Ocean_Basin,
      Coral_Host_Species,
      Sample_Type
    ) %>%

    dplyr::mutate(
      Condition = factor(
        Condition,
        levels = CONDITION_LEVELS
      )
    ) %>%

    dplyr::arrange(
      Condition,
      Comparison_ID
    ) %>%

    dplyr::mutate(
      X_Pos = dplyr::row_number()
    )


  if (
    nrow(lookup) == 0
  ) {

    stop(
      "No comparison-key rows matched ",
      virome_name,
      "."
    )
  }


  N_COMP <- nrow(
    lookup
  )


  # ---------------------------------------------------------------------------
  # Family / phylum order
  # ---------------------------------------------------------------------------

  family_objects <- make_family_map(
    df,
    virome_name
  )

  family_map <- family_objects$family_map
  phylum_order <- family_objects$phylum_order


  N_FAMILY <- nrow(
    family_map
  )


  # ---------------------------------------------------------------------------
  # Attach x/y coordinates
  # ---------------------------------------------------------------------------

  df <- df %>%

    dplyr::left_join(
      lookup %>%

        dplyr::select(
          Comparison_ID,
          X_Pos
        ),
      by = "Comparison_ID"
    ) %>%

    dplyr::left_join(
      family_map %>%

        dplyr::select(
          Family,
          Phylum,
          Y_Pos
        ),
      by = c(
        "Family",
        "Phylum"
      )
    )


  # ---------------------------------------------------------------------------
  # Three fixed-height top rows
  # ---------------------------------------------------------------------------

  Y_COMPARTMENT <- N_FAMILY + 1
  Y_REGION      <- N_FAMILY + 2
  Y_STATUS      <- N_FAMILY + 3


  region_df <- lookup %>%

    dplyr::transmute(
      X_Pos,
      Context_Key = paste0(
        "Region::",
        Ocean_Basin
      ),
      Y_Pos = Y_REGION
    )


  compartment_df <- lookup %>%

    dplyr::transmute(
      X_Pos,
      Context_Key = paste0(
        "Compartment::",
        Sample_Type
      ),
      Y_Pos = Y_COMPARTMENT
    )


  annotation_df <- dplyr::bind_rows(
    region_df,
    compartment_df
  )


  # ---------------------------------------------------------------------------
  # Status blocks
  # ---------------------------------------------------------------------------

  status_blocks <- lookup %>%

    dplyr::group_by(
      Condition
    ) %>%

    dplyr::summarise(
      X_Min = min(
        X_Pos
      ) - 0.5,
      X_Max = max(
        X_Pos
      ) + 0.5,
      X_Mid = mean(
        range(
          X_Pos
        )
      ),
      .groups = "drop"
    ) %>%

    dplyr::mutate(
      Condition_Character = as.character(
        Condition
      ),
      Status_Label = CONDITION_LABELS[
        Condition_Character
      ],
      Status_Fill = STATUS_HEADER_PALETTE[
        Condition_Character
      ]
    )


  # vertical separators between status blocks
  #
  # Do NOT use dplyr::n() inside slice_head(n = ...): `n` must be a constant.
  # Instead, first arrange the blocks and then explicitly remove the final block,
  # because a separator is only needed BETWEEN adjacent status blocks.
  status_boundaries <- status_blocks %>%

    dplyr::arrange(
      X_Min
    )


  if (
    nrow(status_boundaries) > 1
  ) {

    status_boundaries <- status_boundaries[
      seq_len(
        nrow(status_boundaries) - 1L
      ),
      ,
      drop = FALSE
    ]

  } else {

    status_boundaries <- status_boundaries[
      0,
      ,
      drop = FALSE
    ]
  }


  # ---------------------------------------------------------------------------
  # Bubble-grid cells
  # ---------------------------------------------------------------------------

  grid_df <- tidyr::expand_grid(
    X_Pos = seq_len(
      N_COMP
    ),
    Y_Pos = family_map$Y_Pos
  )


  # ---------------------------------------------------------------------------
  # Main square-cell matrix
  # ---------------------------------------------------------------------------

  p_matrix <- ggplot2::ggplot() +

    # family x comparison square cells
    ggplot2::geom_tile(
      data = grid_df,
      ggplot2::aes(
        x = X_Pos,
        y = Y_Pos
      ),
      width = 1,
      height = 1,
      fill = "#FFFFFF",
      colour = "#DDE2E7",
      linewidth = 0.32
    ) +

    # Sampling region / compartment: each exactly one cell high
    ggplot2::geom_tile(
      data = annotation_df,
      ggplot2::aes(
        x = X_Pos,
        y = Y_Pos,
        fill = Context_Key
      ),
      width = 1,
      height = 1,
      colour = "#FFFFFF",
      linewidth = 0.45
    ) +

    # Status header blocks
    ggplot2::geom_rect(
      data = status_blocks,
      ggplot2::aes(
        xmin = X_Min,
        xmax = X_Max,
        ymin = Y_STATUS - 0.5,
        ymax = Y_STATUS + 0.5
      ),
      fill = status_blocks$Status_Fill,
      colour = "#4C5258",
      linewidth = 0.55
    ) +

    ggplot2::geom_text(
      data = status_blocks,
      ggplot2::aes(
        x = X_Mid,
        y = Y_STATUS,
        label = Status_Label
      ),
      family = BASE_FONT,
      fontface = "bold",
      size = 2.45,
      colour = "#232629"
    ) +

    # separators between status blocks
    ggplot2::geom_vline(
      data = status_boundaries,
      ggplot2::aes(
        xintercept = X_Max
      ),
      linewidth = 0.65,
      colour = "#6F757A"
    ) +

    # separator between ecological annotations and viral-family matrix
    ggplot2::geom_hline(
      yintercept = N_FAMILY + 0.5,
      linewidth = 0.55,
      colour = "#8D949A"
    ) +

    # LDA bubbles
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(
        x = X_Pos,
        y = Y_Pos,
        size = LDA_Abs,
        colour = Enriched_In
      ),
      alpha = 0.94
    ) +

    ggplot2::scale_x_continuous(
      limits = c(
        0.5,
        N_COMP + 0.5
      ),
      breaks = lookup$X_Pos,
      labels = lookup$Coral_Host_Species,
      expand = c(
        0,
        0
      )
    ) +

    ggplot2::scale_y_continuous(
      limits = c(
        0.5,
        Y_STATUS + 0.5
      ),
      breaks = NULL,
      expand = c(
        0,
        0
      )
    ) +

    ggplot2::scale_fill_manual(
      values = CONTEXT_COLOURS,
      guide = "none"
    ) +

    ggplot2::scale_colour_manual(
      values = ENRICH_PALETTE,
      limits = names(
        ENRICH_PALETTE
      ),
      labels = ENRICH_LABELS,
      drop = FALSE,
      name = "Enriched in",
      guide = ggplot2::guide_legend(
        order = 1
      )
    ) +

    # exact same LDA-size limits for DNA and RNA
    ggplot2::scale_size_continuous(
      name = "|LDA score|",
      range = c(
        2.4,
        7.3
      ),
      limits = c(
        LDA_THRESHOLD_PLOT,
        GLOBAL_LDA_MAX
      ),
      breaks = LDA_BREAKS,
      guide = ggplot2::guide_legend(
        order = 2
      )
    ) +

    # one x-unit = one y-unit -> square matrix cells
    ggplot2::coord_fixed(
      ratio = 1,
      clip = "off"
    ) +

    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +

    ggplot2::theme_minimal(
      base_size = 7,
      base_family = BASE_FONT
    ) +

    ggplot2::theme(

      panel.grid = ggplot2::element_blank(),

      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),

      axis.text.x = ggplot2::element_text(
        angle = 50,
        hjust = 1,
        vjust = 1,
        face = "italic",
        size = 5.7,
        colour = "#263238"
      ),

      axis.ticks.x = ggplot2::element_blank(),

      legend.position = "right",

      legend.title = ggplot2::element_text(
        size = 7,
        face = "bold"
      ),

      legend.text = ggplot2::element_text(
        size = 6
      ),

      plot.margin = ggplot2::margin(
        0,
        1,
        1,
        0
      )
    )


  # ---------------------------------------------------------------------------
  # Left labels + phylum strip
  #
  # Family labels are sorted:
  #   known "-idae" families first,
  #   then other named labels,
  #   then "Unclass..." labels.
  # ---------------------------------------------------------------------------

  label_df <- family_map %>%

    dplyr::mutate(
      Label_Font = dplyr::if_else(
        grepl(
          "idae$",
          Family,
          ignore.case = TRUE
        ),
        "italic",
        "plain"
      )
    )


  left_track_df <- data.frame(
    Label = c(
      "Compartment",
      "Sampling region"
    ),
    Y_Pos = c(
      Y_COMPARTMENT,
      Y_REGION
    ),
    stringsAsFactors = FALSE
  )


  p_left <- ggplot2::ggplot() +

    # phylum strip
    ggplot2::geom_tile(
      data = label_df,
      ggplot2::aes(
        x = PHYLUM_STRIP_FRAC / 2,
        y = Y_Pos,
        fill = Phylum
      ),
      width = PHYLUM_STRIP_FRAC,
      height = 1
    ) +

    # family names
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(
        x = PHYLUM_STRIP_FRAC + 0.025,
        y = Y_Pos,
        label = Family,
        fontface = Label_Font
      ),
      hjust = 0,
      family = BASE_FONT,
      size = 2.15,
      colour = "#263238"
    ) +

    # top-row labels
    ggplot2::geom_text(
      data = left_track_df,
      ggplot2::aes(
        x = 0.99,
        y = Y_Pos,
        label = Label
      ),
      hjust = 1,
      family = BASE_FONT,
      fontface = "bold",
      size = 2.35,
      colour = "#30343B"
    ) +

    ggplot2::scale_x_continuous(
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +

    ggplot2::scale_y_continuous(
      limits = c(
        0.5,
        Y_STATUS + 0.5
      ),
      expand = c(
        0,
        0
      )
    ) +

    ggplot2::scale_fill_manual(
      values = PHYLUM_PALETTE,
      limits = names(
        PHYLUM_PALETTE
      ),
      drop = FALSE,
      name = "Virus phylum"
    ) +

    ggplot2::theme_void(
      base_family = BASE_FONT
    ) +

    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = 7,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 6
      ),
      plot.margin = ggplot2::margin(
        0,
        0,
        1,
        0
      )
    )


  # ---------------------------------------------------------------------------
  # Remove legends for the actual panel body, then align left labels and matrix.
  # ---------------------------------------------------------------------------

  p_left_body <- p_left +

    ggplot2::theme(
      legend.position = "none"
    )


  p_matrix_body <- p_matrix +

    ggplot2::theme(
      legend.position = "none"
    )


  aligned_body <- cowplot::align_plots(
    p_left_body,
    p_matrix_body,
    align = "h",
    axis = "tb"
  )


  MATRIX_WIDTH_MM <- N_COMP * CELL_MM


  content_row <- cowplot::plot_grid(
    aligned_body[[1]],
    aligned_body[[2]],
    nrow = 1,
    rel_widths = c(
      LEFT_LABEL_MM,
      MATRIX_WIDTH_MM
    ),
    align = "h",
    axis = "tb"
  )


  title_plot <- cowplot::ggdraw() +

    cowplot::draw_label(
      paste0(
        virome_name,
        " viruses"
      ),
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      vjust = 0.5,
      fontfamily = BASE_FONT,
      fontface = "bold",
      size = 10.5,
      colour = "#263238"
    )


  # Natural height:
  # all matrix rows are CELL_MM high, plus space for rotated x labels.
  NATURAL_CONTENT_MM <- (
    N_FAMILY +
      3
  ) * CELL_MM +
    X_LABEL_MM


  panel_natural <- cowplot::plot_grid(
    title_plot,
    content_row,
    ncol = 1,
    rel_heights = c(
      TITLE_MM,
      NATURAL_CONTENT_MM
    ),
    align = "v"
  )


  list(
    natural_panel = panel_natural,
    natural_height_mm = TITLE_MM +
      NATURAL_CONTENT_MM,
    panel_width_mm = LEFT_LABEL_MM +
      MATRIX_WIDTH_MM,
    n_family = N_FAMILY,
    n_comp = N_COMP,
    phylum_order = phylum_order,
    matrix_for_legend = p_matrix,
    left_for_legend = p_left
  )
}


# ==============================================================================
# 25. BUILD DNA + RNA
# ==============================================================================

dna_panel <- build_virome_panel(
  "DNA"
)

rna_panel <- build_virome_panel(
  "RNA"
)


# Export final ordering
phylum_order_output <- dplyr::bind_rows(

  dna_panel$phylum_order %>%

    dplyr::mutate(
      Virome = "DNA"
    ),

  rna_panel$phylum_order %>%

    dplyr::mutate(
      Virome = "RNA"
    )
) %>%

  dplyr::select(
    Virome,
    dplyr::everything()
  )


write.csv(
  phylum_order_output,
  file.path(
    OUTPUT_DIR,
    paste0(
      "ProjectAware_LEfSe_phylum_order_LDA",
      LDA_THRESHOLD_PLOT,
      ".csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================

# ==============================================================================
# 26. Consolidate freshly computed results (never read a precomputed result table)
# ==============================================================================
  key_fields <- c("Virome", "Comparison_ID", "Family")
  key_extra <- comparison_key %>%
    mutate(Comparison_Order = row_number()) %>%
    select(Virome, Comparison_ID, N_Raw_Projects, N_Locations, Comparison_Order)
  selected_extra <- selected %>%
    mutate(Selected_Order = row_number()) %>%
    select(all_of(key_fields), Selected_Order, Direction_Matches_Group_Mean)
  member_extra <- membership %>%
    group_by(Virome, Condition, Stratum_ID) %>%
    summarise(
      Healthy_Samples = paste(sort(Sample[Status_Normalized == "Healthy"]), collapse = ";"),
      Condition_Samples = paste(sort(Sample[Status_Normalized == first(Condition)]), collapse = ";"),
      Healthy_Original_Statuses = as.character(jsonlite::toJSON(
        sort(unique(Status_Original[Status_Normalized == "Healthy"])), auto_unbox = FALSE)),
      Condition_Original_Statuses = as.character(jsonlite::toJSON(
        sort(unique(Status_Original[Status_Normalized == first(Condition)])), auto_unbox = FALSE)),
      .groups = "drop")
  source_data <- all_tests %>%
    left_join(key_extra, by = c("Virome", "Comparison_ID")) %>%
    left_join(selected_extra, by = key_fields) %>%
    mutate(Display_Flag = Retained & is.finite(LDA_Abs) & LDA_Abs >= LDA_THRESHOLD_PLOT) %>%
    left_join(member_extra, by = c("Virome", "Condition", "Stratum_ID")) %>%
    select(all_of(names(all_tests)), N_Raw_Projects, N_Locations, Comparison_Order,
           Selected_Order, Display_Flag, Direction_Matches_Group_Mean, Healthy_Samples,
           Condition_Samples, Healthy_Original_Statuses, Condition_Original_Statuses)
  stopifnot(nrow(source_data) == nrow(all_tests),
            !anyDuplicated(source_data[key_fields]),
            sum(source_data$Retained) == nrow(selected),
            sum(source_data$Display_Flag) == nrow(selected_plot))
  # Quoted TSV uses standard doubled quotes, UTF-8, and empty fields for missing data.
  write.table(source_data, file.path(RESULT_DIR, "source_data.tsv"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = TRUE, qmethod = "double",
              na = "", fileEncoding = "UTF-8")

# ==============================================================================
# 27. Final aligned figure assembly, using panel objects calculated above
# ==============================================================================
obj <- list(dna_panel = dna_panel, rna_panel = rna_panel)
out <- OUTPUT_DIR
CELL <- 5.6
LABEL <- 45
MARGIN <- 8
GAP <- 6
LEGEND <- 48
TITLE <- 8
XLABEL <- 32
PHYLUM_STRIP_FRAC <- 0.085
W <- 2*MARGIN + 2*LABEL + (obj$dna_panel$n_comp+obj$rna_panel$n_comp)*CELL + 2*GAP + LEGEND
H <- 2*MARGIN + TITLE + (obj$dna_panel$n_family+3)*CELL + XLABEL
TOP <- H-MARGIN-TITLE
font <- "Arial"


extract_panel <- function(p) {
  g <- ggplotGrob(p + theme(legend.position="none"))
  g$grobs[[which(g$layout$name=="panel")]]
}

panel_data <- function(p) {
  pm <- p$matrix_for_legend
  # Split the crowded two-column header into two lines; retain its meaning.
  for (i in seq_along(pm$layers)) {
    d <- pm$layers[[i]]$data
    if (is.data.frame(d) && "Status_Label" %in% names(d)) {
      d$Status_Label[d$Status_Label == "Heat stress"] <- "Heat\nstress"
      pm$layers[[i]]$data <- d
    }
  }
  # Capture exact point coordinates and rendered aesthetics of the existing plot.
  built <- ggplot_build(pm)
  pts <- Filter(function(d) all(c("size","colour","shape","x","y") %in% names(d)), built$data)
  stopifnot(length(pts)==1L)
  original_pts <- Filter(function(d) all(c("size","colour","shape","x","y") %in% names(d)), ggplot_build(p$matrix_for_legend)$data)
  stopifnot(identical(pts[[1]], original_pts[[1]]))
  lab <- pm$scales$get_scales("x")$labels
  list(left=extract_panel(p$left_for_legend), matrix=extract_panel(pm), labels=lab, points=pts[[1]], n_comp=p$n_comp, n_family=p$n_family)
}
dna <- panel_data(obj$dna_panel)
rna <- panel_data(obj$rna_panel)
write.csv(rbind(transform(dna$points,Virome="DNA"),transform(rna$points,Virome="RNA")), file.path(out,"layout_repair_point_audit.csv"),row.names=FALSE)

place_grob <- function(g,x,y,w,h) {
  pushViewport(viewport(x=unit(x+w/2,"mm"),y=unit(y+h/2,"mm"),width=unit(w,"mm"),height=unit(h,"mm")))
  grid.draw(g)
  popViewport()
}
txt <- function(s,x,y,cex=7,face="plain",just="left",rot=0) {
  grid.text(s,x=unit(x,"mm"),y=unit(y,"mm"),just=just,rot=rot,gp=gpar(fontfamily=font,fontsize=cex,fontface=face,col="#263238"))
}
draw_panel <- function(p,x,title) {
  h <- (p$n_family+3)*CELL
  y <- TOP-h
  place_grob(p$left,x,y,LABEL,h)
  place_grob(p$matrix,x+LABEL,y,p$n_comp*CELL,h)
  txt(title,x+(LABEL+p$n_comp*CELL)/2,TOP+TITLE/2,10.5,"bold","centre")
  txt("Status",x+LABEL-0.6,TOP-CELL/2,6.7,"bold","right")
  for(i in seq_len(p$n_comp)) {
    txt(p$labels[i],x+LABEL+(i-.5)*CELL,y-1.5,6.5,"italic",just="right",rot=50)
  }
}

colour_scale <- obj$dna_panel$matrix_for_legend$scales$get_scales("colour")
enrich <- colour_scale$palette(1)
enrich_labels <- colour_scale$labels[names(enrich)]
fill_palette <- obj$dna_panel$matrix_for_legend$scales$get_scales("fill")$palette(1)
regions <- fill_palette[startsWith(names(fill_palette),"Region::")]
names(regions) <- sub("^Region::","",names(regions))
compartments <- fill_palette[startsWith(names(fill_palette),"Compartment::")]
names(compartments) <- sub("^Compartment::","",names(compartments))
phyla <- obj$dna_panel$left_for_legend$scales$get_scales("fill")$palette(1)
size_scale <- ggplot_build(obj$dna_panel$matrix_for_legend)$plot$scales$get_scales("size")
lda_breaks <- size_scale$get_breaks()
lda_sizes <- size_scale$map(lda_breaks)

draw_legend <- function(x) {
  y <- H-MARGIN-2
  txt("Enriched in",x,y,8,"bold")
  y <- y-5
  for(i in seq_along(enrich)) {
    grid.points(x=unit(x+2.2,"mm"),y=unit(y,"mm"),pch=19,gp=gpar(col=enrich[i],fontsize=8))
    txt(enrich_labels[i],x+7,y)
    y <- y-5
  }
  y <- y-3
  txt("|LDA score|",x,y,8,"bold")
  y <- y-7
  for(i in seq_along(lda_breaks)) {
    grid.points(x=unit(x+3,"mm"),y=unit(y,"mm"),pch=19,gp=gpar(col="#111111",fontsize=lda_sizes[i]*(72.27/25.4)+.5*(96/25.4)/2,lwd=.5*(96/25.4)/2))
    txt(format(lda_breaks[i],nsmall=1),x+9,y)
    y <- y-8
  }
  for(group in list(list(title="Sampling region",pal=regions),list(title="Compartment",pal=compartments),list(title="Virus phylum",pal=phyla))) {
    y <- y-2
    txt(group$title,x,y,8,"bold")
    y <- y-5
    for(i in seq_along(group$pal)) {
      grid.rect(x=unit(x+2.2,"mm"),y=unit(y,"mm"),width=unit(4.2,"mm"),height=unit(4.2,"mm"),gp=gpar(fill=group$pal[i],col=NA))
      txt(names(group$pal)[i],x+7,y)
      y <- y-5
    }
    y <- y-2
  }
  stopifnot(y > 0)
}
draw <- function() {
  grid.newpage()
  grid.rect(gp=gpar(fill="white",col=NA))
  draw_panel(dna,MARGIN,"DNA viruses")
  xr <- MARGIN+LABEL+dna$n_comp*CELL+GAP
  draw_panel(rna,xr,"RNA viruses")
  draw_legend(xr+LABEL+rna$n_comp*CELL+GAP)
}
stem <- file.path(FIGURE_DIR,"Coral_virome_LEfSe_DSS_exact")
cairo_pdf(paste0(stem,".pdf"),width=W/25.4,height=H/25.4,family=font)
draw(); dev.off()
svglite(paste0(stem,".svg"),width=W/25.4,height=H/25.4)
draw(); dev.off()
agg_png(paste0(stem,".png"),width=W,height=H,units="mm",res=300)
draw(); dev.off()
agg_tiff(paste0(stem,".tiff"),width=W,height=H,units="mm",res=600,compression="lzw")
draw(); dev.off()

# ==============================================================================
# 28. Reproducibility record and input manifest
# ==============================================================================
  input_paths <- file.path(INPUT_DIR, REQUIRED_FILES)
  input_manifest <- data.frame(File = REQUIRED_FILES,
    Bytes = file.info(input_paths)$size, MD5 = unname(tools::md5sum(input_paths)))
  write.table(input_manifest, file.path(RESULT_DIR, "input_manifest.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "sessionInfo.txt"))
  summary_lines <- c(
    "Complete coral virome pipeline: six input tables -> LEfSe/KW -> results -> DNA/RNA figure.",
    paste("UTC run time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("All tested feature/comparison rows:", nrow(all_tests)),
    paste("Eligible comparisons:", nrow(comparison_key)),
    paste("Selected at analysis threshold |LDA| >= 2:", nrow(selected)),
    paste("Displayed at |LDA| >= 2.5:", nrow(selected_plot)),
    paste("DNA / RNA bubbles:", nrow(dna$points), "/", nrow(rna$points)),
    paste("Figure width / height (mm):", W, "/", H),
    "Selection: unadjusted KW P < 0.05, |LDA| >= 2; BH is exported, not used for selection.",
    "DSS: exact trimws(Status) Healthy/Normal versus Disease(Dark Spot Syndrome).",
    "No enrichment-labelled samples enter DSS eligibility or DSS model fitting.",
    "Other comparisons retain their original filtering and statistical framework.",
    "PDF/SVG: vector; PNG: 300 dpi; TIFF: 600 dpi with LZW compression.",
    "Main results: results/source_data.tsv; original diagnostic exports: results/audit_tables/.")
  writeLines(summary_lines, file.path(RESULT_DIR, "run_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")
  invisible(list(all_tests = all_tests, selected = selected, displayed = selected_plot,
                 comparison_key = comparison_key, source_data = source_data,
                 membership = membership, output_directory = output_root))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L)
    stop("Usage: Rscript run_lefse.R [input_directory] [output_directory]")
  if (length(args) == 0L) run_pipeline()
  else if (length(args) == 1L) run_pipeline(input_dir = args[1])
  else run_pipeline(input_dir = args[1], output_dir = args[2])
}
