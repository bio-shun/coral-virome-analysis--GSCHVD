# Nested beta diversity and PCoA
# Standalone R workflow. See ../README.md for inputs, commands and methods.
# Original analyses and input files are never overwritten.
WORKFLOW_SCOPE <- "nested"
WORKFLOW_SCRIPT <- "run_pcoa.R"
WORKFLOW_TITLE <- "Nested beta diversity and PCoA"

# RStudio users: edit CONFIG below, source this file, then run main().
# Command-line users: Rscript <script.R> --help
CONFIG <- list(
  input_dir = "",                 # Folder containing the four input TSV/TXT files.
  output_dir = "",                # A new, empty folder for this workflow.
  mode = "all",                   # all, analysis, figures, or source-data.
  analysis_dir = "",              # Existing results for figures/source-data.
  figure_data_dir = "",           # Beta source-data: archived coordinates/manifest bundle, if separate.
  workers = 4L,
  formats = c("pdf", "svg", "png", "tiff"),
  individual = TRUE,              # Beta workflow: also export individual panels.
  legacy_bundle = NULL,           # Optional verified replay of archived native tests.
  cache_only = FALSE
)


# ---- Analysis functions ----
factor_order <- c("Coral_taxa", "Sampling_region", "Compartments", "Status")

factor_labels <- c(Coral_taxa = "Coral taxa", Sampling_region = "Sampling region", Compartments = "Compartments", 
    Status = "Status")

check_analysis_packages <- function() {
    packages <- c("vegan", "permute", "dplyr", "stringr", "digest")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Install required packages: ", paste(missing, collapse = ", "))
    suppressPackageStartupMessages(library(dplyr))
    invisible(TRUE)
}

clean_species <- function(x) {
    x <- stringr::str_squish(trimws(as.character(x)))
    sub("^p?ocillopora\\b", "Pocillopora", x, ignore.case = TRUE, perl = TRUE)
}

clean_compartment <- function(x) {
    x <- stringr::str_squish(trimws(as.character(x)))
    x <- sub("^Coral\\s+", "", x, ignore.case = TRUE)
    x[grepl("mucus.*tissue|tissue.*mucus", x, ignore.case = TRUE)] <- "Mucus-tissue slurry"
    key <- tolower(x)
    replacements <- c(tissue = "Tissue", mucus = "Mucus", skeleton = "Skeleton", fragment = "Fragment", polyp = "Polyp", 
        `endolithic band` = "Endolithic band", `surface and interior segments` = "Surface and interior segments", 
        `mucus-tissue slurry` = "Mucus-tissue slurry")
    for (k in names(replacements)) x[key == k] <- replacements[[k]]
    x
}

clean_status <- function(x) {
    x <- stringr::str_squish(trimws(as.character(x)))
    out <- rep("Others", length(x))
    out[grepl("^Healthy/Normal$|^Healthy/Normal\\(Exposed\\)$", x, ignore.case = TRUE)] <- "Healthy/Normal"
    out[grepl("^Bleached$", x, ignore.case = TRUE)] <- "Bleached"
    out[grepl("SCTLD", x, ignore.case = TRUE)] <- "Disease(SCTLD)"
    out[grepl("black band", x, ignore.case = TRUE)] <- "Diseased(black band disease)"
    out[grepl("white plague", x, ignore.case = TRUE)] <- "Diseased(white plague disease)"
    out[grepl("dark spot", x, ignore.case = TRUE)] <- "Disease(Dark Spot Syndrome)"
    out[grepl("bleaching lesion", x, ignore.case = TRUE)] <- "Disease (Bleaching lesion)"
    out[grepl("pathogen", x, ignore.case = TRUE)] <- "Disease (Pathogen-induced disease)"
    out[grepl("heat", x, ignore.case = TRUE) & grepl("acidification", x, ignore.case = TRUE)] <- "Stressed(Combined)"
    out[grepl("heat", x, ignore.case = TRUE) & !grepl("acidification", x, ignore.case = TRUE)] <- "Stressed(Heat stress)"
    out[grepl("acidification", x, ignore.case = TRUE) & !grepl("heat", x, ignore.case = TRUE)] <- "Stressed(Acidification stress)"
    out[grepl("cold", x, ignore.case = TRUE)] <- "Stressed(Cold stress)"
    out
}

read_dataset <- function(virus, input_dir) {
    stopifnot(virus %in% c("DNA", "RNA"))
    abundance_file <- file.path(input_dir, if (virus == "DNA") 
        "24083_metagenome_TPM.tsv"
    else "rna_virus_tpm.tsv")
    metadata_file <- file.path(input_dir, if (virus == "DNA") 
        "sample_info_metagenome.tsv"
    else "sample_info_metatranscriptome.tsv")
    if (!all(file.exists(c(abundance_file, metadata_file)))) 
        stop("Missing input files for ", virus)
    abundance <- read.delim(abundance_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "", 
        comment.char = "")
    abundance <- as.matrix(abundance)
    storage.mode(abundance) <- "numeric"
    if (any(!is.finite(abundance)) || any(abundance < 0)) 
        stop(virus, ": TPM must be finite and non-negative")
    if (anyDuplicated(colnames(abundance))) 
        stop(virus, ": duplicated abundance sample IDs")
    abundance <- abundance[rowSums(abundance, na.rm = TRUE) > 0, , drop = FALSE]
    metadata <- read.delim(metadata_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
    if (!"Sample" %in% names(metadata)) 
        names(metadata)[1] <- "Sample"
    required <- c("Sample", "Project", "Ocean_Basin", "Coral_Host_Species", "Sample_Type", "Status")
    if (!all(required %in% names(metadata))) 
        stop("Missing metadata columns: ", paste(setdiff(required, names(metadata)), collapse = ", "))
    metadata$Sample <- stringr::str_squish(trimws(as.character(metadata$Sample)))
    if (anyDuplicated(metadata$Sample)) 
        stop(virus, ": duplicated Sample identifiers in metadata")
    common <- intersect(colnames(abundance), metadata$Sample)
    if (length(common) < 2) 
        stop(virus, ": fewer than two samples after ID matching")
    abundance <- abundance[, common, drop = FALSE]
    metadata <- metadata[match(common, metadata$Sample), , drop = FALSE]
    nonzero_sample <- colSums(abundance, na.rm = TRUE) > 0
    abundance <- abundance[, nonzero_sample, drop = FALSE]
    metadata <- metadata[nonzero_sample, , drop = FALSE]
    metadata_clean <- data.frame(Sample = metadata$Sample, Project = stringr::str_squish(trimws(as.character(metadata$Project))), 
        Sampling_region = stringr::str_squish(trimws(as.character(metadata$Ocean_Basin))), Coral_taxa = clean_species(metadata$Coral_Host_Species), 
        Compartments = clean_compartment(metadata$Sample_Type), Status = clean_status(metadata$Status), stringsAsFactors = FALSE)
    valid <- complete.cases(metadata_clean) & metadata_clean$Project != "" & metadata_clean$Project != "Unknown" & 
        metadata_clean$Sampling_region != "" & metadata_clean$Sampling_region != "Unknown" & metadata_clean$Coral_taxa != 
        "" & metadata_clean$Coral_taxa != "Unknown" & metadata_clean$Compartments != "" & metadata_clean$Compartments != 
        "Unknown"
    abundance <- abundance[, valid, drop = FALSE]
    metadata_clean <- metadata_clean[valid, , drop = FALSE]
    rownames(metadata_clean) <- metadata_clean$Sample
    if (nrow(metadata_clean) < 3) 
        stop(virus, ": fewer than three valid samples")
    for (fac in c("Project", factor_order)) metadata_clean[[fac]] <- droplevels(factor(metadata_clean[[fac]]))
    community <- sqrt(t(abundance))
    distance <- vegan::vegdist(community, method = "bray")
    message(virus, ": ", nrow(metadata_clean), " samples, ", nrow(abundance), " viral features, ", nlevels(metadata_clean$Project), 
        " Projects")
    list(virus = virus, abundance = abundance, community = community, metadata = metadata_clean, distance = distance, 
        n_features = nrow(abundance))
}

make_project_perm <- function(project, permutations = 9999L) {
    ctrl <- permute::how(nperm = permutations)
    permute::setBlocks(ctrl) <- droplevels(factor(project))
    ctrl
}

within_project_variation <- function(metadata, factor_name) {
    split_values <- split(as.character(metadata[[factor_name]]), metadata$Project)
    sum(vapply(split_values, function(x) length(unique(x)) > 1, logical(1)))
}

canonical_checkpoint_value <- function(x) {
    if (is.null(x)) 
        return(NULL)
    if (is.environment(x)) 
        return(list(Environment_omitted = TRUE))
    if (is.function(x)) 
        stop("Unexpected function in checkpoint specification")
    attributes_clean <- function(z) {
        a <- attributes(z)
        if (is.null(a)) 
            return(NULL)
        a[[".Environment"]] <- NULL
        lapply(a, canonical_checkpoint_value)
    }
    if (is.language(x) || is.pairlist(x)) 
        return(list(Language = paste(deparse(x, width.cutoff = 500L), collapse = "\n"), Attributes = attributes_clean(x)))
    if (is.list(x)) {
        values <- lapply(unclass(x), canonical_checkpoint_value)
        attributes(values) <- attributes_clean(x)
        return(values)
    }
    if (is.atomic(x)) {
        values <- x
        attributes(values) <- attributes_clean(x)
        return(values)
    }
    stop("Unsupported checkpoint object: ", paste(class(x), collapse = ", "))
}

stable_checkpoint_signature <- function(specification) {
    stopifnot(identical(names(specification), c("Virus", "Type", "Design", "Requested", "Permutations")))
    list(Signature_version = "canonical_no_formula_environment_v1", Specification = canonical_checkpoint_value(specification))
}

stable_checkpoint_key <- function(specification) {
    digest::digest(stable_checkpoint_signature(specification), algo = "md5")
}

legacy_checkpoint_index <- local({
    indices <- new.env(parent = emptyenv())
    function(bundle) {
        bundle <- normalizePath(bundle, winslash = "/", mustWork = TRUE)
        if (exists(bundle, envir = indices, inherits = FALSE)) 
            return(get(bundle, envir = indices))
        index_path <- file.path(bundle, "Documentation/stable_checkpoint_cache_index.csv")
        if (file.exists(index_path)) {
            index <- read.csv(index_path, check.names = FALSE, stringsAsFactors = FALSE)
            result <- list(index = index, has_signature_files = TRUE)
        }
        else {
            execution_path <- file.path(bundle, "Data/Analysis_Results/10_permutation_execution_audit.csv")
            if (!file.exists(execution_path)) 
                stop("Missing legacy execution audit: ", execution_path)
            execution <- read.csv(execution_path, check.names = FALSE, stringsAsFactors = FALSE)
            if (anyDuplicated(execution$Test_ID)) 
                stop("Duplicate native legacy Test_ID values")
            rows <- lapply(seq_len(nrow(execution)), function(i) {
                recorded <- execution[i, , drop = FALSE]
                checkpoint <- file.path(bundle, recorded$Checkpoint[1])
                if (!file.exists(checkpoint)) 
                  stop("Missing legacy checkpoint: ", checkpoint)
                saved <- readRDS(checkpoint)
                native_id <- as.character(saved$audit$Test_ID[1])
                if (!identical(sub("\\.rds$", "", basename(checkpoint)), native_id) || !identical(native_id, as.character(recorded$Test_ID[1]))) 
                  stop("Legacy filename/Test_ID mismatch")
                if (!identical(as.integer(saved$audit$Seed[1]), as.integer(recorded$Seed[1])) || nrow(saved$permutation_matrix) != 
                  recorded$Actual_permutations[1] || nrow(saved$permutation_matrix) != saved$audit$Actual_permutations[1]) 
                  stop("Legacy audit/payload mismatch: ", native_id)
                data.frame(Stable_Key = stable_checkpoint_key(saved$specification), Native_Test_ID = native_id, 
                  Archived_seed = recorded$Seed[1], Actual_permutations = recorded$Actual_permutations[1], Checkpoint = recorded$Checkpoint[1], 
                  Checkpoint_file_MD5 = unname(tools::md5sum(checkpoint)), stringsAsFactors = FALSE)
            })
            index <- dplyr::bind_rows(rows)
            if (!nrow(index)) 
                stop("No archived native tests in legacy execution audit")
            if (anyDuplicated(index$Stable_Key)) 
                stop("Ambiguous canonical legacy checkpoint keys")
            result <- list(index = index, has_signature_files = FALSE)
            message("Indexed ", nrow(index), " native legacy tests in memory (source bundle unchanged)")
        }
        required <- c("Stable_Key", "Native_Test_ID", "Archived_seed", "Actual_permutations", "Checkpoint", "Checkpoint_file_MD5")
        if (!all(required %in% names(result$index))) 
            stop("Invalid legacy checkpoint index schema")
        assign(bundle, result, envir = indices)
        result
    }
})

load_legacy_checkpoint <- function(bundle, specification) {
    if (is.null(bundle)) 
        return(NULL)
    lookup <- legacy_checkpoint_index(bundle)
    index <- lookup$index
    key <- stable_checkpoint_key(specification)
    row <- index[index$Stable_Key == key, , drop = FALSE]
    if (!nrow(row)) 
        return(NULL)
    if (nrow(row) != 1L) 
        stop("Ambiguous legacy cache key: ", key)
    signature <- stable_checkpoint_signature(specification)
    if (lookup$has_signature_files) {
        stored_signature <- readRDS(file.path(bundle, "Documentation/Stable_checkpoint_signatures", paste0(key, 
            ".rds")))
        if (!identical(signature, stored_signature)) 
            stop("Legacy canonical signature mismatch: ", key)
    }
    checkpoint <- file.path(bundle, row$Checkpoint)
    if (!file.exists(checkpoint)) 
        stop("Missing legacy checkpoint: ", checkpoint)
    if (!identical(unname(tools::md5sum(checkpoint)), row$Checkpoint_file_MD5)) 
        stop("Legacy checkpoint checksum mismatch: ", key)
    saved <- readRDS(checkpoint)
    if (!identical(signature, stable_checkpoint_signature(saved$specification))) 
        stop("Legacy native payload mismatch: ", key)
    if (!identical(as.character(saved$audit$Test_ID[1]), row$Native_Test_ID)) 
        stop("Legacy Test_ID mismatch: ", key)
    if (!identical(as.integer(saved$audit$Seed[1]), as.integer(row$Archived_seed))) 
        stop("Legacy seed mismatch: ", key)
    if (nrow(saved$permutation_matrix) != row$Actual_permutations) 
        stop("Legacy permutation count mismatch: ", key)
    attr(saved, "verified_checkpoint_path") <- normalizePath(checkpoint, winslash = "/", mustWork = TRUE)
    saved
}

guard_checkpoint_read <- function(reader) {
    tryCatch(reader(), error = function(e) {
        class(e) <- unique(c("checkpoint_validation_error", class(e)))
        stop(e)
    })
}

must_rethrow_engine_error <- function(error, engine) {
    isTRUE(engine$cache_only) || inherits(error, "checkpoint_validation_error") ||
        grepl("No verified cached native test", conditionMessage(error), fixed = TRUE) ||
        grepl("publication permutation plan|Permutation plan design mismatch|Rebuilt permutation matrix differs", conditionMessage(error))
}

new_permutation_engine <- function(output_dir, workers = 1L, nperm = 9999L, seed = 20260902L, legacy_bundle = NULL,
    cache = TRUE, cache_only = FALSE, permutation_plan = NULL) {
    # The fixed publication snapshot uses the original recorded randomization
    # plan. Seeds rebuild the assignments; no archived statistic is read here.
    plan <- if (is.null(permutation_plan)) NULL else read.csv(permutation_plan, stringsAsFactors = FALSE)
    if (!is.null(plan)) {
        required <- c("Stable_key", "Original_test_ID", "Virus", "Type", "N", "Requested_permutations",
                      "Actual_permutations", "Seed", "RNG_kind", "Normal_kind", "Sample_kind", "Permutation_SHA256")
        if (!all(required %in% names(plan)) || anyDuplicated(plan$Stable_key)) stop("Invalid permutation plan")
    }
    state <- new.env(parent = emptyenv())
    state$virus <- ""
    state$context <- list(Analysis = "", Factor = "", Stratum = "", Project = "")
    state$audit <- list()
    state$mapping <- list()
    state$cluster <- NULL
    checkpoint_dir <- file.path(output_dir, "Permutation_Checkpoints")
    if (cache) 
        dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    set_context <- function(analysis, factor = "", stratum = "", project = "") {
        state$context <- list(Analysis = analysis, Factor = factor, Stratum = stratum, Project = project)
        invisible(NULL)
    }
    flush <- function() {
        if (length(state$audit)) 
            write.csv(dplyr::bind_rows(state$audit), file.path(output_dir, "10_permutation_execution_audit.csv"), 
                row.names = FALSE, na = "")
        if (length(state$mapping)) 
            write.csv(dplyr::bind_rows(state$mapping), file.path(output_dir, "11_test_context_map.csv"), row.names = FALSE, 
                na = "")
        invisible(NULL)
    }
    register <- function(saved, origin, checkpoint) {
        row <- saved$audit
        row$Checkpoint <- checkpoint
        row$Cache_origin <- origin
        state$audit[[as.character(row$Test_ID[1])]] <- row
        state$mapping[[length(state$mapping) + 1L]] <- cbind(data.frame(Virus = state$virus, Analysis = state$context$Analysis, 
            Factor = state$context$Factor, Stratum = state$context$Stratum, Project = state$context$Project, stringsAsFactors = FALSE), 
            row[, setdiff(names(row), c("Virus", "Model")), drop = FALSE])
        flush()
        saved$result
    }
    audit_test <- function(type, design, N, permutations, runner) {
        specification <- list(Virus = state$virus, Type = type, Design = design, Requested = nperm, Permutations = permutations)
        key <- stable_checkpoint_key(specification)
        signature <- stable_checkpoint_signature(specification)
        checkpoint <- file.path(checkpoint_dir, paste0(key, ".rds"))
        if (cache && file.exists(checkpoint)) {
            saved <- guard_checkpoint_read(function() {
                value <- readRDS(checkpoint)
                if (!identical(signature, stable_checkpoint_signature(value$specification))) 
                  stop("Checkpoint signature mismatch: ", key)
                if (nrow(value$permutation_matrix) != value$audit$Actual_permutations[1]) 
                  stop("Checkpoint permutation count mismatch: ", key)
                value
            })
            message("Cached ", state$virus, " ", type, " ", key)
            return(register(saved, "Current output", file.path("Permutation_Checkpoints", basename(checkpoint))))
        }
        saved <- guard_checkpoint_read(function() load_legacy_checkpoint(legacy_bundle, specification))
        if (!is.null(saved)) {
            legacy_path <- attr(saved, "verified_checkpoint_path")
            message("Legacy cache ", state$virus, " ", type, " ", saved$audit$Test_ID[1])
            return(register(saved, "Verified legacy bundle", legacy_path))
        }
        if (cache_only) 
            stop("No verified cached native test: ", state$virus, " ", type, " ", key)
        test_seed <- as.integer((as.double(strtoi(substr(key, 1, 7), base = 16L)) + as.double(seed))%%2147483646) + 
            1L
        planned <- NULL
        if (!is.null(plan)) {
            planned <- plan[plan$Stable_key == key, , drop = FALSE]
            if (nrow(planned) != 1L) stop("Test is absent or duplicated in the publication permutation plan: ", key)
            if (planned$Virus != state$virus || planned$Type != type || planned$N != N ||
                planned$Requested_permutations != nperm) stop("Permutation plan design mismatch: ", key)
            test_seed <- as.integer(planned$Seed)
            RNGkind(kind = planned$RNG_kind, normal.kind = planned$Normal_kind, sample.kind = planned$Sample_kind)
        }
        set.seed(test_seed)
        ctrl <- if (is.numeric(permutations) && length(permutations) == 1L) 
            permute::how(nperm = nperm)
        else permutations
        pm <- vegan:::getPermuteMatrix(ctrl, N)
        actual <- nrow(pm)
        permutation_sha256 <- digest::digest(list(dim = dim(pm), assignments = as.integer(pm)), algo = "sha256")
        if (!is.null(planned) && (actual != planned$Actual_permutations ||
            permutation_sha256 != planned$Permutation_SHA256)) stop("Rebuilt permutation matrix differs from the publication plan: ", key)
        complete <- isTRUE(permute::getComplete(attr(pm, "control")))
        if (actual < nperm && !complete) 
            stop("Fewer requested permutations without complete enumeration: ", key)
        if (actual < 1L) 
            stop("No non-identity permutations available: ", key)
        execution <- 1L
        if (N >= 40L && workers > 1L) {
            if (is.null(state$cluster)) {
                state$cluster <- parallel::makePSOCKcluster(workers, outfile = file.path(output_dir, "worker_log.txt"))
                parallel::clusterEvalQ(state$cluster, suppressPackageStartupMessages(library(vegan)))
            }
            execution <- state$cluster
        }
        message("Starting ", state$virus, " ", type, " N=", N, " permutations=", actual)
        started <- proc.time()[["elapsed"]]
        result <- runner(pm, execution)
        fperm <- if (type == "PERMDISP") 
            as.matrix(result$perm)
        else attr(result, "F.perm")
        stopifnot(is.matrix(fperm), nrow(fperm) == actual)
        row <- data.frame(Test_ID = key, Virus = state$virus, Type = type, N = N, Requested_permutations = nperm, 
            Actual_permutations = actual, Complete_enumeration = complete, Minimum_P_resolution = 1/(actual + 1), 
            Seed = test_seed, Workers = if (N >= 40L) 
                workers
            else 1L, Elapsed_seconds = proc.time()[["elapsed"]] - started, Model = design$Model, Checkpoint = if (cache) 
                file.path("Permutation_Checkpoints", basename(checkpoint))
            else "")
        row$Original_test_ID <- if (is.null(planned)) "" else planned$Original_test_ID
        row$Permutation_SHA256 <- permutation_sha256
        saved <- list(specification = specification, audit = row, permutation_matrix = pm, result = result)
        if (cache) 
            saveRDS(saved, checkpoint, compress = "gzip")
        register(saved, "Computed", row$Checkpoint[1])
    }
    adonis <- function(formula, data, permutations = nperm, by = NULL, ...) {
        dist <- eval(formula[[2]], envir = environment(formula))
        distance <- dist
        design <- list(Model = paste(deparse(formula), collapse = ""), Data = data, Distance = unclass(dist), By = by)
        execution_formula <- formula
        execution_formula[[2]] <- as.name("distance")
        environment(execution_formula) <- environment()
        audit_test("PERMANOVA", design, nrow(data), permutations, function(pm, execution) vegan::adonis2(execution_formula, 
            data = data, permutations = pm, by = by, parallel = execution, ...))
    }
    anova <- function(object, permutations = nperm, ...) {
        design <- list(Model = paste(deparse(formula(object)), collapse = ""), Fit = object)
        audit_test("partial dbRDA", design, nobs(object), permutations, function(pm, execution) vegan::anova.cca(object, 
            permutations = pm, parallel = execution, ...))
    }
    dispersion <- function(x, permutations = nperm, ...) {
        design <- list(Model = "Distance to group spatial median", Fit = x)
        audit_test("PERMDISP", design, length(x$distances), permutations, function(pm, execution) vegan::permutest(x, 
            permutations = pm, parallel = execution, ...))
    }
    close <- function() {
        if (!is.null(state$cluster)) {
            parallel::stopCluster(state$cluster)
            state$cluster <- NULL
        }
        invisible(NULL)
    }
    list(adonis = adonis, anova = anova, dispersion = dispersion, set_context = set_context, set_virus = function(virus) {
        state$virus <- virus
        invisible(NULL)
    }, close = close, flush = flush, nperm = nperm, cache_only = cache_only)
}

safe_dispersion <- function(distance, grouping, engine, project = NULL) {
    grouping <- droplevels(factor(grouping))
    if (nlevels(grouping) < 2 || any(table(grouping) < 2)) 
        return(c(F = NA_real_, P = NA_real_, Note = "Insufficient replication for PERMDISP"))
    bd <- tryCatch(vegan::betadisper(distance, grouping), error = function(e) e)
    if (inherits(bd, "error")) 
        return(c(F = NA_real_, P = NA_real_, Note = conditionMessage(bd)))
    perm <- if (is.null(project)) 
        engine$nperm
    else make_project_perm(project, engine$nperm)
    test <- tryCatch(engine$dispersion(bd, permutations = perm), error = function(e) e)
    if (inherits(test, "error")) {
        if (must_rethrow_engine_error(test, engine)) 
            stop(test)
        return(c(F = NA_real_, P = NA_real_, Note = conditionMessage(test)))
    }
    tab <- as.data.frame(test$tab)
    c(F = unname(tab$F[1]), P = unname(tab$`Pr(>F)`[1]), Note = "")
}

safe_adonis <- function(distance, metadata, factor_name, engine, blocked = FALSE) {
    if (nlevels(droplevels(factor(metadata[[factor_name]]))) < 2) 
        return(list(Df = NA_real_, R2 = NA_real_, F = NA_real_, P = NA_real_, Estimable = FALSE, Note = "Fewer than two factor levels"))
    perm <- if (blocked) 
        make_project_perm(metadata$Project, engine$nperm)
    else engine$nperm
    fit <- tryCatch(engine$adonis(as.formula(paste("distance ~", factor_name)), data = metadata, permutations = perm), 
        error = function(e) e)
    if (inherits(fit, "error")) {
        if (must_rethrow_engine_error(fit, engine)) 
            stop(fit)
        return(list(Df = NA_real_, R2 = NA_real_, F = NA_real_, P = NA_real_, Estimable = FALSE, Note = conditionMessage(fit)))
    }
    tab <- as.data.frame(fit)
    list(Df = tab$Df[1], R2 = tab$R2[1], F = tab$F[1], P = tab$`Pr(>F)`[1], Estimable = TRUE, Note = "")
}

run_nested_factor <- function(dataset, focal, engine, strata = character(), healthy_only = FALSE, require_healthy = FALSE) {
    metadata <- dataset$metadata
    if (healthy_only) 
        metadata <- metadata %>% filter(Status == "Healthy/Normal")
    if (nrow(metadata) < 3) 
        return(data.frame())
    if (length(strata) == 0) 
        groups <- tibble(.group_id = "Global")
    else groups <- metadata %>% distinct(across(all_of(strata)))
    rows <- lapply(seq_len(nrow(groups)), function(i) {
        if (!length(strata)) {
            sub_meta <- metadata
            stratum_label <- "Global"
        }
        else {
            sub_meta <- metadata
            for (s in strata) {
                value <- as.character(groups[[s]][i])
                sub_meta <- sub_meta %>% filter(as.character(.data[[s]]) == value)
            }
            stratum_label <- paste(vapply(strata, function(s) paste0(factor_labels[[s]], "=", as.character(groups[[s]][i])), 
                character(1)), collapse = " | ")
        }
        counts <- sub_meta %>% count(.data[[focal]], name = "n") %>% filter(n >= 3)
        valid_levels <- as.character(counts[[focal]])
        if (length(valid_levels) < 2 || (require_healthy && !"Healthy/Normal" %in% valid_levels)) 
            return(NULL)
        sub_meta <- sub_meta %>% filter(as.character(.data[[focal]]) %in% valid_levels) %>% droplevels()
        sub_meta[[focal]] <- droplevels(factor(sub_meta[[focal]]))
        dist_sub <- as.dist(as.matrix(dataset$distance)[rownames(sub_meta), rownames(sub_meta), drop = FALSE])
        engine$set_context("Nested", focal, stratum_label)
        pooled <- safe_adonis(dist_sub, sub_meta, focal, engine, blocked = FALSE)
        projects_varied <- within_project_variation(sub_meta, focal)
        blocked <- if (projects_varied > 0) 
            safe_adonis(dist_sub, sub_meta, focal, engine, blocked = TRUE)
        else list(Df = NA_real_, R2 = NA_real_, F = NA_real_, P = NA_real_, Estimable = FALSE, Note = "No within-Project variation in this stratum")
        disp_pool <- safe_dispersion(dist_sub, sub_meta[[focal]], engine)
        disp_block <- if (projects_varied > 0) 
            safe_dispersion(dist_sub, sub_meta[[focal]], engine, sub_meta$Project)
        else c(F = NA_real_, P = NA_real_, Note = "No within-Project variation in this stratum")
        data.frame(Virus = dataset$virus, Factor = focal, Stratum = stratum_label, N = nrow(sub_meta), Levels = nlevels(sub_meta[[focal]]), 
            Projects = nlevels(sub_meta$Project), Projects_within_factor_variation = projects_varied, R2_pooled = pooled$R2, 
            F_pooled = pooled$F, P_pooled = pooled$P, R2_blocked = blocked$R2, F_blocked = blocked$F, P_blocked = blocked$P, 
            Dispersion_F_pooled = as.numeric(disp_pool[["F"]]), Dispersion_P_pooled = as.numeric(disp_pool[["P"]]), 
            Dispersion_F_blocked = as.numeric(disp_block[["F"]]), Dispersion_P_blocked = as.numeric(disp_block[["P"]]), 
            Blocked_note = blocked$Note, stringsAsFactors = FALSE)
    })
    result <- bind_rows(rows)
    if (!nrow(result)) 
        return(result)
    result %>% group_by(Virus, Factor) %>% mutate(Q_pooled = p.adjust(P_pooled, method = "BH"), Q_blocked = p.adjust(P_blocked, 
        method = "BH"), Q_dispersion_pooled = p.adjust(Dispersion_P_pooled, method = "BH"), Q_dispersion_blocked = p.adjust(Dispersion_P_blocked, 
        method = "BH")) %>% ungroup()
}

analyse_nested <- function(dataset, engine) {
    bind_rows(run_nested_factor(dataset, "Sampling_region", engine, character(), healthy_only = TRUE), run_nested_factor(dataset, 
        "Coral_taxa", engine, "Sampling_region", healthy_only = TRUE), run_nested_factor(dataset, "Compartments", 
        engine, c("Sampling_region", "Coral_taxa"), healthy_only = TRUE), run_nested_factor(dataset, "Status", 
        engine, c("Sampling_region", "Coral_taxa", "Compartments"), healthy_only = FALSE, require_healthy = TRUE))
}

dataset_summary <- function(dataset) {
    metadata <- dataset$metadata
    data.frame(Virus = dataset$virus, Final_samples = nrow(metadata), Viral_features = dataset$n_features, Projects = nlevels(metadata$Project), 
        Sampling_region_levels = nlevels(metadata$Sampling_region), Coral_taxa_levels = nlevels(metadata$Coral_taxa), 
        Compartments_levels = nlevels(metadata$Compartments), Status_levels = nlevels(metadata$Status), stringsAsFactors = FALSE)
}

# This workflow writes a new analysis folder; input files are read-only.
run_analysis <- function(input_dir, output_dir, workers = 1L,
                         nperm = 9999L, seed = 20260902L,
                         legacy_bundle = NULL, cache_only = FALSE) {
  check_analysis_packages()
  workers <- as.integer(workers)
  nperm <- as.integer(nperm)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) stop("workers must be a positive integer")
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L) stop("nperm must be a positive integer")
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (identical(tolower(input_dir), tolower(output_dir))) stop("Input and output directories must differ")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Analysis folder is not empty; choose a new output folder: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(input_dir, c("24083_metagenome_TPM.tsv", "sample_info_metagenome.tsv",
                                 "rna_virus_tpm.tsv", "sample_info_metatranscriptome.tsv", "permutation_plan.csv"))
  if (!all(file.exists(files))) stop("Required abundance or metadata input is missing")
  checksums <- data.frame(File = basename(files), MD5 = unname(tools::md5sum(files)))
  write.csv(checksums, file.path(output_dir, "00_input_checksums.csv"), row.names = FALSE)
  engine <- new_permutation_engine(output_dir, workers, nperm, seed,
                                    legacy_bundle = legacy_bundle, cache_only = cache_only,
                                    permutation_plan = file.path(input_dir, "permutation_plan.csv"))
  on.exit(engine$close(), add = TRUE)
  results <- list()
  for (virus in c("DNA", "RNA")) {
    engine$set_virus(virus)
    dataset <- read_dataset(virus, input_dir)
    saveRDS(list(metadata = dataset$metadata, distance = dataset$distance, n_features = dataset$n_features),
      file.path(output_dir, paste0(virus, "_filtered_distance_metadata.rds")))
    write.csv(dataset$metadata, file.path(output_dir, paste0(virus, "_filtered_metadata.csv")), row.names = FALSE)
    result <- list(nested = analyse_nested(dataset, engine))
    result$sample_summary <- dataset_summary(dataset)
    results[[virus]] <- result
    saveRDS(result, file.path(output_dir, paste0(virus, "_completed_analysis.rds")))
    rm(dataset)
    invisible(gc(FALSE))
  }
  tables <- c(nested = "07_nested_factor_comparisons.csv", sample_summary = "08_sample_summary.csv")
  for (name in names(tables)) {
    values <- lapply(results, function(x) x[[name]])
    write.csv(dplyr::bind_rows(values), file.path(output_dir, tables[[name]]), row.names = FALSE, na = "")
  }
  nested <- dplyr::bind_rows(lapply(results, function(x) x$nested))
  if (nrow(nested)) {
    long <- dplyr::bind_rows(
      dplyr::transmute(nested, Virus, Factor, Stratum, R2 = R2_pooled, Inference = "Pooled"),
      dplyr::transmute(nested, Virus, Factor, Stratum, R2 = R2_blocked, Inference = "Project-blocked"))
    long <- long[!is.na(long$R2), , drop = FALSE]
    write.csv(long, file.path(output_dir, "07b_nested_effect_sizes_long.csv"), row.names = FALSE, na = "")
  }
  engine$flush()
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "09_sessionInfo.txt"))
  writeLines(c(
    "DNA and RNA are analysed separately using square-root TPM and Bray-Curtis distances.",
    paste0("Requested permutations: ", nperm, ". Native actual counts and finite-space enumeration are in 10_permutation_execution_audit.csv."),
    "Per-test identities and ecological strata are linked in 11_test_context_map.csv.",
    "Nested comparisons require at least two retained levels with at least three samples per level.",
    "Sampling region, Coral taxa and Compartments use Healthy/Normal samples; Status comparisons require a retained Healthy/Normal level.",
    "Blocked tests restrict permutations within Project but do not subtract Project sums of squares.",
    "No within-Project focal variation is recorded as not estimable for blocked tests.",
    "BH correction is within Virus x Factor, separately for pooled/blocked PERMANOVA/PERMDISP families.",
    "All eligible Status tests remain in tables; figure panels show nominal pooled P < 0.05 only.",
    "R-squared values from different nested subsets are not directly comparable independent effects.",
    "Sample libraries are not automatically independent biological replicates.",
    "Recorded per-test seeds reconstruct the publication permutation matrices; assignment SHA-256 values are checked before every newly computed test.",
    "Optional verified legacy replay preserves archived native tests without writing to the legacy folder."
  ), file.path(output_dir, "README_method_and_interpretation.txt"))
  message("Analysis completed: ", output_dir)
  invisible(results)
}

# ---- Figure functions ----
fig_factor_order <- c("Sampling_region", "Coral_taxa", "Compartments", "Status")

fig_factor_labels <- c(Sampling_region = "Sampling region", Coral_taxa = "Coral taxa", Compartments = "Compartments", 
    Status = "Status")

fig_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Missing R packages: ", paste(missing, collapse = ", "), ". Install these packages before plotting.")
}

fig_analysis_dir <- function(analysis_root) {
    root <- normalizePath(analysis_root, winslash = "/", mustWork = TRUE)
    nested <- file.path(root, "Data", "Analysis_Results")
    if (dir.exists(nested)) 
        nested
    else root
}

fig_read_table <- function(root, filename, columns = character()) {
    path <- file.path(root, filename)
    if (!file.exists(path)) 
        stop("Missing analysis table: ", path)
    out <- read.csv(path, check.names = FALSE, na.strings = c("", "NA"), stringsAsFactors = FALSE)
    absent <- setdiff(columns, names(out))
    if (length(absent)) 
        stop(filename, " is missing columns: ", paste(absent, collapse = ", "))
    out
}

fig_format_p <- function(p) {
    if (length(p) != 1L || is.na(p)) 
        return("NA")
    if (!is.numeric(p) || !is.finite(p) || p < 0 || p > 1) 
        stop("Invalid P/Q value")
    if (p == 0) 
        stop("Zero permutation P value; inspect the analysis before plotting")
    if (p < 1e-06) 
        return(format(p, digits = 15, scientific = TRUE, trim = TRUE))
    for (places in 0:15) {
        candidate <- round(p, places)
        crossed <- any((candidate < c(0.001, 0.01, 0.05)) != (p < c(0.001, 0.01, 0.05)))
        if (candidate > 0 && !crossed && abs(candidate - p) <= max(1e-16, abs(p) * 1e-14)) {
            return(formatC(candidate, format = "f", digits = places))
        }
    }
    format(p, digits = 16, scientific = FALSE, trim = TRUE)
}

fig_write_csv <- function(data, filename) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(data, filename, na = "")
}

fig_safe_name <- function(x) substr(gsub("[^A-Za-z0-9._-]+", "_", x), 1, 150)

fig_save <- function(plot, stem, width, height, formats, family = "sans") {
    formats <- unique(tolower(formats))
    if (length(setdiff(formats, c("pdf", "svg", "tiff", "png")))) {
        stop("formats must contain only pdf, svg, tiff and/or png")
    }
    dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
    files <- character()
    for (fmt in formats) {
        path <- paste0(stem, ".", fmt)
        switch(fmt, pdf = grDevices::cairo_pdf(path, width = width, height = height, family = family, onefile = TRUE), 
            svg = if (requireNamespace("svglite", quietly = TRUE)) {
                svglite::svglite(path, width = width, height = height)
            } else grDevices::svg(path, width = width, height = height, family = family), tiff = ragg::agg_tiff(path, 
                width = width, height = height, units = "in", res = 600, compression = "lzw"), png = ragg::agg_png(path, 
                width = width, height = height, units = "in", res = 300))
        tryCatch(print(plot), finally = grDevices::dev.off())
        files <- c(files, path)
    }
    invisible(files)
}

fig_palette <- function(focal, groups) {
    if (focal == "Coral_taxa") {
        cols <- grDevices::hcl(h = seq(15, 375, length.out = length(groups) + 1)[seq_along(groups)], c = 70, l = 55)
        return(stats::setNames(cols, groups))
    }
    base <- switch(focal, Sampling_region = c(`Arabian Gulf` = "#009E73", `Atlantic Ocean` = "#0072B2", `Caribbean Sea` = "#D55E00", 
        `Indian Ocean` = "#CC79A7", `Mediterranean Sea` = "#E69F00", `Pacific Ocean` = "#56B4E9", `Red Sea` = "#6A3D9A"), 
        Compartments = c(Tissue = "#D95F02", Mucus = "#1B9E77", Skeleton = "#8D40CA", Fragment = "#E72989", `Endolithic band` = "#377EB8", 
            Polyp = "#3EB855", `Surface and interior segments` = "#E6AB02", `Mucus-tissue slurry` = "#66C2A5"), 
        Status = c(`Healthy/Normal` = "#33A02C", Bleached = "#E31A1C", `Disease(SCTLD)` = "#542788", `Diseased(black band disease)` = "#40004B", 
            `Diseased(white plague disease)` = "#9970AB", `Disease(Dark Spot Syndrome)` = "#762A83", `Disease (Bleaching lesion)` = "#FD8D3C", 
            `Disease (Pathogen-induced disease)` = "#8073AC", `Stressed(Heat stress)` = "#FF7F00", `Stressed(Acidification stress)` = "#FDBF6F", 
            `Stressed(Combined)` = "#B15928", `Stressed(Cold stress)` = "#1F78B4", Others = "#999999"))
    extra <- setdiff(groups, names(base))
    if (length(extra)) {
        if (focal == "Sampling_region") 
            stop("No predefined sampling-region colour: ", paste(extra, collapse = ", "))
        add <- grDevices::hcl(h = seq(20, 340, length.out = length(extra) + 1)[seq_along(extra)], c = 60, l = 55)
        base <- c(base, stats::setNames(add, extra))
    }
    base[groups]
}

fig_read_dataset <- function(root, virus) {
    path <- file.path(root, paste0(virus, "_filtered_distance_metadata.rds"))
    if (!file.exists(path)) 
        stop("Missing saved distance/metadata: ", path)
    saved <- readRDS(path)
    if (!all(c("distance", "metadata", "n_features") %in% names(saved))) 
        stop("Incomplete saved dataset: ", virus)
    metadata <- as.data.frame(saved$metadata)
    if (!all(c("Sample", "Project", fig_factor_order) %in% names(metadata))) 
        stop("Missing metadata columns: ", virus)
    if (anyDuplicated(metadata$Sample)) 
        stop("Duplicated sample IDs: ", virus)
    rownames(metadata) <- as.character(metadata$Sample)
    ids <- attr(saved$distance, "Labels")
    if (is.null(ids) || !setequal(ids, rownames(metadata))) 
        stop("Distance and metadata sample IDs disagree: ", virus)
    metadata <- metadata[ids, , drop = FALSE]
    for (f in c("Project", fig_factor_order)) metadata[[f]] <- droplevels(factor(metadata[[f]]))
    list(virus = virus, distance = saved$distance, metadata = metadata, n_features = saved$n_features)
}

fig_nested_metadata <- function(dataset, focal, strata, healthy, require_healthy) {
    metadata <- dataset$metadata
    if (healthy) 
        metadata <- metadata[metadata$Status == "Healthy/Normal", , drop = FALSE]
    groups <- if (length(strata)) 
        unique(metadata[, strata, drop = FALSE])
    else data.frame(.group = "Global")
    out <- list()
    for (i in seq_len(nrow(groups))) {
        sub <- metadata
        for (s in strata) sub <- sub[as.character(sub[[s]]) == as.character(groups[[s]][i]), , drop = FALSE]
        counts <- table(sub[[focal]])
        levels <- names(counts)[counts >= 3]
        if (length(levels) < 2 || (require_healthy && !"Healthy/Normal" %in% levels)) 
            next
        sub <- droplevels(sub[as.character(sub[[focal]]) %in% levels, , drop = FALSE])
        sub[[focal]] <- droplevels(factor(sub[[focal]]))
        out[[length(out) + 1L]] <- sub
    }
    out
}

fig_stratum <- function(metadata, strata) {
    if (!length(strata)) 
        return("Global")
    paste(vapply(strata, function(s) paste0(fig_factor_labels[[s]], "=", unique(as.character(metadata[[s]]))[1]), 
        character(1)), collapse = " | ")
}

fig_panel_title <- function(focal, metadata, virus) {
    region <- unique(as.character(metadata$Sampling_region))[1]
    taxa <- unique(as.character(metadata$Coral_taxa))[1]
    compartment <- unique(as.character(metadata$Compartments))[1]
    if (focal == "Sampling_region") 
        return(bquote(.(virus) ~ " | Healthy/Normal samples by sampling region"))
    if (focal == "Coral_taxa") 
        return(bquote(.(virus) ~ " | " ~ .(region)))
    if (focal == "Compartments") 
        return(bquote(.(virus) ~ " | " ~ italic(.(taxa)) ~ " | " ~ .(region)))
    bquote(atop(.(virus) ~ " | " ~ italic(.(taxa)), .(region) ~ " | " ~ .(compartment)))
}

fig_significance <- function(p) {
    if (is.na(p)) 
        "NA"
    else if (p <= 0.001) 
        "***"
    else if (p <= 0.01) 
        "**"
    else if (p <= 0.05) 
        "*"
    else "ns"
}

fig_stats_subtitle <- function(r2, p) {
    r <- if (is.na(r2)) 
        "NA"
    else sprintf("%.1f%%", 100 * r2)
    label <- fig_format_p(p)
    bquote(italic(R)^2 == .(r) ~ ";" ~ italic(P) == .(label) ~ " " ~ .(fig_significance(p)))
}

fig_pcoa_theme <- function() {
    ggplot2::theme_classic(base_size = 7.5, base_family = "sans") + ggplot2::theme(panel.border = ggplot2::element_rect(linewidth = 0.35, 
        colour = "black", fill = NA), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_line(linewidth = 0.35, 
        colour = "black"), axis.text = ggplot2::element_text(colour = "black"), panel.grid = ggplot2::element_blank(), 
        plot.title = ggplot2::element_text(size = 8.6, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 2)), 
        plot.subtitle = ggplot2::element_text(size = 7.2, hjust = 0.5, colour = "#3F4650", margin = ggplot2::margin(b = 3)), 
        legend.title = ggplot2::element_text(size = 7, face = "bold"), legend.text = ggplot2::element_text(size = 6.5), 
        legend.key.height = grid::unit(0.28, "cm"), legend.key.width = grid::unit(0.35, "cm"), aspect.ratio = 1, 
        plot.margin = ggplot2::margin(4, 5, 4, 5))
}

fig_pcoa_panel <- function(dataset, metadata, focal, stats_row, colours, legend_levels) {
    ids <- rownames(metadata)
    full_distance <- as.matrix(dataset$distance)
    if (!all(ids %in% rownames(full_distance))) 
        stop("Missing distance sample IDs")
    pc <- ape::pcoa(stats::as.dist(full_distance[ids, ids, drop = FALSE]))
    if (ncol(pc$vectors) < 2L) 
        stop("PCoA has fewer than two axes: ", stats_row$Stratum[1])
    points <- as.data.frame(pc$vectors[, 1:2, drop = FALSE])
    names(points) <- c("Axis1", "Axis2")
    points$Sample <- rownames(points)
    metadata <- metadata[points$Sample, , drop = FALSE]
    points$Group <- as.character(metadata[[focal]])
    centroids <- stats::aggregate(points[, c("Axis1", "Axis2")], list(Group = points$Group), mean)
    names(centroids)[2:3] <- c("CentroidX", "CentroidY")
    points$CentroidX <- centroids$CentroidX[match(points$Group, centroids$Group)]
    points$CentroidY <- centroids$CentroidY[match(points$Group, centroids$Group)]
    counts <- table(points$Group)
    points$Group_n <- as.integer(counts[points$Group])
    hulls <- do.call(rbind, lapply(names(counts)[counts == 3], function(g) {
        d <- points[points$Group == g, , drop = FALSE]
        h <- grDevices::chull(d$Axis1, d$Axis2)
        d[c(h, h[1]), c("Axis1", "Axis2", "Group"), drop = FALSE]
    }))
    lines <- points[points$Group_n == 2, , drop = FALSE]
    p <- ggplot2::ggplot(points, ggplot2::aes(Axis1, Axis2, fill = Group)) + ggplot2::geom_segment(data = points[points$Group_n >= 
        2, , drop = FALSE], ggplot2::aes(x = Axis1, y = Axis2, xend = CentroidX, yend = CentroidY, group = Group), 
        inherit.aes = FALSE, linewidth = 0.25, colour = "#9AA0A6", alpha = 0.55)
    if (any(points$Group_n >= 4)) 
        p <- p + ggplot2::stat_ellipse(data = points[points$Group_n >= 4, , drop = FALSE], ggplot2::aes(x = Axis1, 
            y = Axis2, group = Group, fill = Group, colour = Group), geom = "polygon", type = "t", level = 0.8, 
            alpha = 0.13, linewidth = 0.35, inherit.aes = FALSE, show.legend = FALSE)
    if (!is.null(hulls) && nrow(hulls)) 
        p <- p + ggplot2::geom_polygon(data = hulls, ggplot2::aes(Axis1, Axis2, group = Group, fill = Group, colour = Group), 
            alpha = 0.13, linewidth = 0.35, inherit.aes = FALSE, show.legend = FALSE)
    if (nrow(lines)) 
        p <- p + ggplot2::geom_path(data = lines, ggplot2::aes(Axis1, Axis2, group = Group, colour = Group), linewidth = 0.45, 
            alpha = 0.65, inherit.aes = FALSE, show.legend = FALSE)
    p <- p + ggplot2::geom_point(ggplot2::aes(fill = Group), shape = 21, colour = "black", stroke = 0.35, size = 2.4, 
        alpha = 0.92)
    if (focal == "Coral_taxa") {
        centroids$label <- centroids$Group
        p <- p + ggrepel::geom_text_repel(data = centroids, ggplot2::aes(CentroidX, CentroidY, label = label), 
            parse = FALSE, fontface = "italic", inherit.aes = FALSE, size = 2.15, colour = "#20252B", box.padding = 0.3, 
            point.padding = 0.15, min.segment.length = 0, segment.colour = "#8C939B", show.legend = FALSE)
    }
    if (focal == "Status") {
        centroids$label <- stringr::str_wrap(centroids$Group, width = 20)
        p <- p + ggrepel::geom_text_repel(data = centroids, ggplot2::aes(CentroidX, CentroidY, label = label), 
            inherit.aes = FALSE, parse = FALSE, fontface = "plain", size = 1.9, colour = "#20252B", box.padding = 0.32, 
            point.padding = 0.15, min.segment.length = 0, force = 0.65, max.overlaps = Inf, direction = "both", 
            seed = 20260904, segment.colour = "#8C939B", segment.size = 0.22, show.legend = FALSE)
    }
    p1 <- pc$values$Relative_eig[1] * 100
    p2 <- pc$values$Relative_eig[2] * 100
    p <- p + ggplot2::scale_fill_manual(values = colours, breaks = legend_levels, drop = FALSE) + ggplot2::scale_colour_manual(values = colours, 
        breaks = legend_levels, guide = "none", drop = FALSE) + ggplot2::labs(title = fig_panel_title(focal, metadata, 
        dataset$virus), subtitle = fig_stats_subtitle(as.numeric(stats_row$R2_pooled[1]), as.numeric(stats_row$P_pooled[1])), 
        x = sprintf("PCoA1 (%.1f%%)", p1), y = sprintf("PCoA2 (%.1f%%)", p2), fill = fig_factor_labels[[focal]]) + 
        ggplot2::coord_equal(expand = TRUE) + fig_pcoa_theme()
    if (focal == "Status") 
        p <- p + ggplot2::theme(plot.title = ggplot2::element_text(size = 6.8, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 1)))
    if (focal == "Coral_taxa") 
        p <- p + ggplot2::theme(legend.position = "none")
    list(plot = p, points = points, meta = metadata, focal = focal, stats = stats_row, p1 = p1, p2 = p2)
}

fig_composite <- function(focal, records, colours, levels) {
    records <- records[order(vapply(records, function(x) x$virus, character(1)), vapply(records, function(x) x$stratum, 
        character(1)))]
    panels <- lapply(records, function(x) x$plot + ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 7.2, 
        face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(size = 7.2, hjust = 0.5), plot.margin = ggplot2::margin(3, 
        3, 3, 3)))
    ncol <- max(1L, min(if (focal == "Sampling_region") 2L else 4L, length(panels)))
    nrow <- ceiling(length(panels)/ncol)
    panels <- c(panels, rep(list(patchwork::plot_spacer()), nrow * ncol - length(panels)))
    body <- patchwork::wrap_plots(panels, ncol = ncol, byrow = TRUE)
    if (focal == "Status") 
        levels <- levels[levels %in% unique(unlist(lapply(records, function(x) x$points$Group)))]
    legend_data <- data.frame(Group = factor(levels, levels = levels), x = 1, y = 1)
    legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x, y, fill = Group)) + ggplot2::geom_point(shape = 21, 
        size = 2.7, colour = "black", stroke = 0.25) + ggplot2::scale_fill_manual(values = colours, breaks = levels, 
        limits = levels, drop = FALSE) + ggplot2::labs(fill = fig_factor_labels[[focal]]) + ggplot2::theme_void(base_family = "sans") + 
        ggplot2::theme(legend.position = "bottom", legend.box = "horizontal", legend.title = ggplot2::element_text(size = 7, 
            face = "bold"), legend.text = ggplot2::element_text(size = 6.2, face = if (focal == "Coral_taxa") 
            "italic"
        else "plain"), legend.key.height = grid::unit(0.25, "cm"), legend.key.width = grid::unit(0.32, "cm")) + 
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = if (focal == "Coral_taxa") 
            4
        else if (focal == "Status") 
            3
        else 2, byrow = TRUE, title.position = "top"))
    out <- body + patchwork::plot_annotation(title = paste0(fig_factor_labels[[focal]], "-associated coral viromes"), 
        subtitle = "DNA and RNA; sqrt(TPM)-Bray-Curtis PCoA. Labels show uncorrected pooled PERMANOVA P values.", 
        tag_levels = "a", theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5), 
            plot.subtitle = ggplot2::element_text(size = 6.8, colour = "#59616A", hjust = 0.5), plot.margin = ggplot2::margin(5, 
                5, 5, 5)))
    out <- cowplot::plot_grid(out, cowplot::get_legend(legend_plot), ncol = 1, align = "v", rel_heights = c(1, 
        switch(focal, Coral_taxa = 0.2, Status = 0.13, Sampling_region = 0.17, 0.1)))
    list(plot = out, width = 3.45 * ncol, height = 3.45 * nrow + switch(focal, Coral_taxa = 2, Status = 1.45, Sampling_region = 1.5, 
        1.25))
}

fig_export_stats <- function(stats, output_root, manifest, requested = NA_integer_) {
    tables <- list()
    for (focal in fig_factor_order) {
        data <- stats[stats$Factor == focal, , drop = FALSE]
        key <- paste(data$Virus, data$Factor, data$Stratum)
        matched <- manifest[match(key, paste(manifest$Virus, manifest$Factor, manifest$Stratum)), , drop = FALSE]
        if (anyNA(matched$Coordinates)) 
            stop("Missing coordinate manifest entry for ", focal)
        data$Figure_panel_included <- matched$Figure_panel_included
        data$Coordinates <- matched$Coordinates
        data$PCoA1_percent <- matched$PCoA1_percent
        data$PCoA2_percent <- matched$PCoA2_percent
        for (column in intersect(c("P_pooled", "P_blocked", "Dispersion_P_pooled", "Dispersion_P_blocked", "Q_pooled", 
            "Q_blocked", "Q_dispersion_pooled", "Q_dispersion_blocked"), names(data))) {
            data[[paste0(column, "_label")]] <- vapply(data[[column]], fig_format_p, character(1))
        }
        data$Permutations_requested <- requested
        data$Nominal_minimum_Monte_Carlo_P <- 1/(requested + 1)
        data$Permutation_count_note <- paste("Actual counts are recorded per test in", "Data/Analysis_Results/10_permutation_execution_audit.csv; small strata may use exhaustive restricted permutations")
        data$Figure_P_source <- "P_pooled (uncorrected pooled PERMANOVA)"
        fig_write_csv(data, file.path(output_root, "Data", "Statistics", paste0("Stats_", focal, ".csv")))
        tables[[focal]] <- data
    }
    all_stats <- do.call(rbind, tables)
    fig_write_csv(all_stats[all_stats$Figure_panel_included, , drop = FALSE], file.path(output_root, "Data", "Statistics", 
        "Figure_panel_statistics.csv"))
    invisible(all_stats)
}

render_nested_figures <- function(analysis_root, output_root, formats = c("pdf", "svg", "tiff", "png"), individual = TRUE) {
    fig_require(c("ape", "ggplot2", "patchwork", "cowplot", "ggrepel", "stringr", "readr", if (any(c("tiff", "png") %in% 
        formats)) "ragg"))
    set.seed(20260903)
    root <- fig_analysis_dir(analysis_root)
    output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
    stats <- fig_read_table(root, "07_nested_factor_comparisons.csv", c("Virus", "Factor", "Stratum", "N", "Levels", 
        "R2_pooled", "P_pooled", "R2_blocked", "P_blocked", "Q_pooled", "Q_blocked"))
    stats_key <- paste(stats$Virus, stats$Factor, stats$Stratum)
    if (anyDuplicated(stats_key)) 
        stop("Duplicated nested statistics rows")
    datasets <- stats::setNames(lapply(c("DNA", "RNA"), function(v) fig_read_dataset(root, v)), c("DNA", "RNA"))
    levels <- stats::setNames(lapply(fig_factor_order, function(f) sort(unique(unlist(lapply(datasets, function(d) as.character(d$metadata[[f]])))))), 
        fig_factor_order)
    colours <- stats::setNames(lapply(fig_factor_order, function(f) fig_palette(f, levels[[f]])), fig_factor_order)
    specs <- list(list(focal = "Sampling_region", strata = character(), healthy = TRUE, require_healthy = FALSE), 
        list(focal = "Coral_taxa", strata = "Sampling_region", healthy = TRUE, require_healthy = FALSE), list(focal = "Compartments", 
            strata = c("Sampling_region", "Coral_taxa"), healthy = TRUE, require_healthy = FALSE), list(focal = "Status", 
            strata = c("Sampling_region", "Coral_taxa", "Compartments"), healthy = FALSE, require_healthy = TRUE))
    records <- list()
    manifest <- list()
    files <- character()
    for (virus in c("DNA", "RNA")) for (spec in specs) {
        dataset <- datasets[[virus]]
        metas <- fig_nested_metadata(dataset, spec$focal, spec$strata, spec$healthy, spec$require_healthy)
        for (metadata in metas) {
            stratum <- fig_stratum(metadata, spec$strata)
            row <- stats[stats_key == paste(virus, spec$focal, stratum), , drop = FALSE]
            if (nrow(row) != 1L) 
                stop("No unique statistics row for ", virus, " / ", spec$focal, " / ", stratum)
            if (nrow(metadata) != row$N || nlevels(metadata[[spec$focal]]) != row$Levels) {
                stop("Plot sample set disagrees with analysed N/Levels: ", virus, " / ", stratum)
            }
            shown <- spec$focal != "Status" || (!is.na(row$P_pooled) && row$P_pooled < 0.05)
            rec <- fig_pcoa_panel(dataset, metadata, spec$focal, row, colours[[spec$focal]], levels[[spec$focal]])
            key <- fig_safe_name(paste(virus, spec$focal, gsub("Sampling region=|Coral taxa=|Compartments=", "", 
                stratum), sep = "__"))
            rec$key <- key
            rec$virus <- virus
            rec$stratum <- stratum
            rec$included <- shown
            coordinates <- paste0("Data/PCoA_coordinates/", key, "_coordinates.csv")
            individual_stem <- paste0("Figures/Individual/", virus, "/", spec$focal, "/", key)
            fig_write_csv(rec$points, file.path(output_root, coordinates))
            if (shown && individual) 
                files <- c(files, fig_save(rec$plot, file.path(output_root, individual_stem), 6.2, 4.8, formats))
            records[[length(records) + 1L]] <- rec
            manifest[[length(manifest) + 1L]] <- data.frame(Panel = length(records), Virus = virus, Factor = spec$focal, 
                Stratum = stratum, N = nrow(metadata), Levels = nlevels(metadata[[spec$focal]]), R2_pooled = row$R2_pooled, 
                P_pooled = row$P_pooled, R2_blocked = row$R2_blocked, P_blocked = row$P_blocked, P_label = fig_format_p(row$P_pooled), 
                Q_pooled = row$Q_pooled, Q_blocked = row$Q_blocked, Q_pooled_label = fig_format_p(row$Q_pooled), 
                Q_blocked_label = fig_format_p(row$Q_blocked), PCoA1_percent = rec$p1, PCoA2_percent = rec$p2, 
                Figure_panel_included = shown, Individual_plot = if (shown && individual) 
                  individual_stem
                else "", Coordinates = coordinates, stringsAsFactors = FALSE)
        }
    }
    if (!length(manifest)) 
        stop("No eligible nested comparisons")
    manifest <- do.call(rbind, manifest)
    manifest <- manifest[order(match(manifest$Factor, fig_factor_order), match(manifest$Virus, c("DNA", "RNA")), 
        manifest$Stratum), , drop = FALSE]
    if (nrow(manifest) != nrow(stats) || !setequal(stats_key, paste(manifest$Virus, manifest$Factor, manifest$Stratum))) 
        stop("Manifest does not cover every analysed stratum")
    fig_write_csv(manifest, file.path(output_root, "Data/Statistics/PCoA_exact_p_manifest.csv"))
    fig_write_csv(manifest[manifest$Figure_panel_included, , drop = FALSE], file.path(output_root, "Data/Statistics/Figure_panel_manifest.csv"))
    colour_manifest <- do.call(rbind, lapply(fig_factor_order, function(f) data.frame(Factor = f, Group = names(colours[[f]]), 
        Colour = unname(colours[[f]]))))
    fig_write_csv(colour_manifest, file.path(output_root, "Data/Statistics/PCoA_colour_manifest.csv"))
    audit_path <- file.path(root, "10_permutation_execution_audit.csv")
    requested <- if (file.exists(audit_path)) 
        unique(read.csv(audit_path)$Requested_permutations)
    else NA_integer_
    if (length(requested) != 1L) 
        requested <- NA_integer_
    fig_export_stats(stats, output_root, manifest, requested)
    composites <- list()
    for (focal in fig_factor_order) {
        chosen <- records[vapply(records, function(x) x$focal == focal && x$included, logical(1))]
        if (!length(chosen)) 
            next
        out <- fig_composite(focal, chosen, colours[[focal]], levels[[focal]])
        composites[[focal]] <- out$plot
        files <- c(files, fig_save(out$plot, file.path(output_root, "Figures", paste0("PCoA_FourFactor_", focal, 
            "_DNA_RNA")), out$width, out$height, formats))
    }
    message("PCoA figures completed; all ", nrow(manifest), " comparisons retained in statistics/coordinates.")
    invisible(list(manifest = manifest, colours = colour_manifest, composites = composites, records = records, 
        files = files))
}

# ---- Organised figure source data ----
# Source data are copied from the saved plotting coordinates and test results.
# This path never calls an ordination or permutation routine.
sdn_read <- function(path) {
  if (!file.exists(path)) stop("Missing saved figure data: ", path,
    ". Supply the completed figure bundle with --figure-data-dir; no ordination is silently recalculated.")
  as.data.frame(readr::read_csv(path, na = c("", "NA"), show_col_types = FALSE,
    name_repair = "minimal", progress = FALSE))
}

sdn_write <- function(data, stem) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, paste0(stem, ".csv"), na = "NA")
  readr::write_tsv(data, paste0(stem, ".tsv"), na = "NA")
}

sdn_portable_trace <- function(data) {
  for (column in names(data)) {
    if (!is.character(data[[column]]) || !grepl("path|checkpoint|directory|file", column, ignore.case = TRUE)) next
    absolute <- !is.na(data[[column]]) & grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\]{1,2})", data[[column]], perl = TRUE)
    data[[column]][absolute] <- basename(gsub("\\\\", "/", data[[column]][absolute]))
  }
  data
}

sdn_bundle <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (basename(path) == "Analysis_Results" && basename(dirname(path)) == "Data") {
    path <- dirname(dirname(path))
  }
  path
}

sdn_key <- function(data) paste(data$Virus, data$Factor, data$Stratum, sep = "\r")

sdn_equal <- function(x, y, what) {
  # Historical tables used different serializers (15 versus 17 digits).
  # Permit only serialization-scale numeric differences; retain the analysis value.
  equal <- if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
    identical(is.na(x), is.na(y)) && all(abs(x[!is.na(x)] - y[!is.na(y)]) <=
      pmax(1e-14, abs(x[!is.na(x)]) * 1e-13))
  } else isTRUE(all.equal(x, y, check.attributes = FALSE))
  if (!equal) {
    stop("Saved plotting data disagree with analysis results: ", what)
  }
}

sd_nested_export <- function(analysis_root, output_root, figure_data_dir = NULL) {
  fig_require(c("readr", "digest"))
  root <- fig_analysis_dir(analysis_root)
  bundle <- sdn_bundle(if (is.null(figure_data_dir) || !nzchar(figure_data_dir)) analysis_root else figure_data_dir)
  destination <- file.path(output_root, "Source_Data")
  if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE))) {
    stop("Source_Data is not empty: ", destination)
  }
  stats <- sdn_read(file.path(root, "07_nested_factor_comparisons.csv"))
  manifest <- sdn_read(file.path(bundle, "Data", "Statistics", "PCoA_exact_p_manifest.csv"))
  colours <- sdn_read(file.path(bundle, "Data", "Statistics", "PCoA_colour_manifest.csv"))
  keys <- sdn_key(stats)
  if (anyDuplicated(keys) || anyDuplicated(sdn_key(manifest)) ||
      !setequal(keys, sdn_key(manifest))) stop("Analysis/figure manifest keys are not unique and complete")
  datasets <- setNames(lapply(c("DNA", "RNA"), function(v) fig_read_dataset(root, v)), c("DNA", "RNA"))
  summary <- list()
  source_files <- c(file.path(root, "07_nested_factor_comparisons.csv"),
    file.path(bundle, "Data", "Statistics", c("PCoA_exact_p_manifest.csv", "PCoA_colour_manifest.csv")),
    file.path(root, paste0(c("DNA", "RNA"), "_filtered_distance_metadata.rds")))
  dictionary <- list()
  for (focal in fig_factor_order) {
    base <- stats[stats$Factor == focal, , drop = FALSE]
    base <- base[order(base$Virus, base$Stratum), , drop = FALSE]
    m <- manifest[match(sdn_key(base), sdn_key(manifest)), , drop = FALSE]
    shown <- focal != "Status" | (!is.na(base$P_pooled) & base$P_pooled < 0.05)
    sdn_equal(as.logical(m$Figure_panel_included), shown, paste(focal, "display rule"))
    for (column in intersect(names(base), names(m))) {
      sdn_equal(base[[column]], m[[column]], paste(focal, column))
    }
    # Preserve existing derived labels/per-test permutation counts, after validating
    # every original analysis column. Never recalculate a P or BH-adjusted Q.
    legacy_stats <- file.path(bundle, "Data", "Statistics", paste0("Stats_", focal, ".csv"))
    if (file.exists(legacy_stats)) {
      extra <- sdn_read(legacy_stats)
      if (anyDuplicated(sdn_key(extra)) || !setequal(sdn_key(base), sdn_key(extra))) stop("Incomplete Stats table: ", focal)
      extra <- extra[match(sdn_key(base), sdn_key(extra)), , drop = FALSE]
      for (column in names(base)) sdn_equal(base[[column]], extra[[column]], paste(focal, "Stats", column))
      base <- cbind(base, extra[, setdiff(names(extra), names(base)), drop = FALSE])
      source_files <- c(source_files, legacy_stats)
    }
    figure <- paste0("PCoA_FourFactor_", focal, "_DNA_RNA")
    order_in_figure <- rep(NA_integer_, nrow(base))
    order_in_figure[shown] <- seq_len(sum(shown))
    panel <- rep(NA_character_, nrow(base))
    panel[shown] <- letters[seq_len(sum(shown))]
    id <- data.frame(Figure = figure, Composite_panel_order = order_in_figure,
      Composite_panel_label = panel, Virus = base$Virus, Factor = base$Factor,
      Factor_label = unname(fig_factor_labels[focal]), Stratum = base$Stratum,
      Figure_panel_included = shown, stringsAsFactors = FALSE)
    base$Figure_panel_included <- shown
    base$PCoA1_percent <- m$PCoA1_percent
    base$PCoA2_percent <- m$PCoA2_percent
    base$Coordinates <- NULL
    base$Source_coordinate_file <- m$Coordinates
    base$Figure_P_source <- "P_pooled (uncorrected pooled PERMANOVA)"
    statistic <- cbind(id, base[, setdiff(names(base), names(id)), drop = FALSE])
    point_list <- vector("list", nrow(base))
    for (i in seq_len(nrow(base))) {
      relative <- m$Coordinates[i]
      if (is.na(relative) || grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", relative)) {
        stop("Coordinate path must be relative to the figure bundle")
      }
      path <- file.path(bundle, relative)
      point <- sdn_read(path)
      required <- c("Sample", "Group", "Axis1", "Axis2", "CentroidX", "CentroidY", "Group_n")
      if (!all(required %in% names(point))) stop("Incomplete coordinate columns: ", relative)
      if (anyDuplicated(point$Sample) || nrow(point) != base$N[i] ||
          length(unique(point$Group)) != base$Levels[i]) stop("Coordinate N/levels mismatch: ", relative)
      if (any(!is.finite(as.matrix(point[, c("Axis1", "Axis2", "CentroidX", "CentroidY")])))) {
        stop("Non-finite coordinates: ", relative)
      }
      metadata <- datasets[[base$Virus[i]]]$metadata
      if (!all(point$Sample %in% metadata$Sample)) stop("Coordinate sample absent from metadata")
      metadata <- metadata[match(point$Sample, metadata$Sample), , drop = FALSE]
      sdn_equal(as.character(metadata[[focal]]), point$Group, paste(relative, "groups"))
      specs <- switch(focal, Sampling_region = character(), Coral_taxa = "Sampling_region",
        Compartments = c("Sampling_region", "Coral_taxa"), Status = c("Sampling_region", "Coral_taxa", "Compartments"))
      if (any(vapply(specs, function(s) length(unique(as.character(metadata[[s]]))) != 1L, logical(1)))) {
        stop("Coordinate samples cross ecological strata: ", relative)
      }
      sdn_equal(fig_stratum(metadata, specs), base$Stratum[i], paste(relative, "stratum"))
      expected <- fig_nested_metadata(datasets[[base$Virus[i]]], focal, specs,
        healthy = focal != "Status", require_healthy = focal == "Status")
      expected <- expected[vapply(expected, function(d) fig_stratum(d, specs) == base$Stratum[i], logical(1))]
      if (length(expected) != 1L || !setequal(expected[[1]]$Sample, point$Sample)) {
        stop("Archived coordinates do not cover the exact analysed sample subset: ", relative)
      }
      if (focal != "Status" && any(as.character(metadata$Status) != "Healthy/Normal")) stop("Unhealthy sample in baseline plot")
      for (group in unique(point$Group)) {
        at <- point$Group == group
        sdn_equal(point$Group_n[at], rep(sum(at), sum(at)), paste(relative, "group size"))
        sdn_equal(point$CentroidX[at], rep(mean(point$Axis1[at]), sum(at)), paste(relative, "x centroid"))
        sdn_equal(point$CentroidY[at], rep(mean(point$Axis2[at]), sum(at)), paste(relative, "y centroid"))
      }
      colour <- colours$Colour[match(paste(focal, point$Group), paste(colours$Factor, colours$Group))]
      if (anyNA(colour)) stop("Missing palette entry: ", focal)
      annotation <- id[rep(i, nrow(point)), , drop = FALSE]
      rownames(annotation) <- NULL
      point_list[[i]] <- cbind(annotation, point[, required, drop = FALSE],
        data.frame(Project = as.character(metadata$Project),
          Sampling_region = as.character(metadata$Sampling_region), Coral_taxa = as.character(metadata$Coral_taxa),
          Compartments = as.character(metadata$Compartments), Status = as.character(metadata$Status),
          PCoA1_percent = m$PCoA1_percent[i], PCoA2_percent = m$PCoA2_percent[i], Colour = colour,
          Group_envelope = ifelse(point$Group_n >= 4, "80% t ellipse", ifelse(point$Group_n == 3, "convex hull", "line")),
          stringsAsFactors = FALSE))
      source_files <- c(source_files, path)
    }
    points <- do.call(rbind, point_list)
    folder <- file.path(destination, focal)
    sdn_write(points, file.path(folder, "points"))
    sdn_write(statistic, file.path(folder, "statistics"))
    summary[[focal]] <- data.frame(Figure = figure, Factor = focal, Comparisons = nrow(base),
      Displayed_panels = sum(shown), Sample_coordinate_rows = nrow(points))
    dictionary[[focal]] <- rbind(
      data.frame(Table = paste0(focal, "/points"), Column = names(points)),
      data.frame(Table = paste0(focal, "/statistics"), Column = names(statistic)))
  }
  sdn_write(do.call(rbind, summary), file.path(destination, "figure_index"))
  # This optional table retains finite-space counts and exact test identifiers.
  details <- file.path(bundle, "Data", "Statistics", "Nested_permutation_execution_details.csv")
  if (file.exists(details)) {
    sdn_write(sdn_portable_trace(sdn_read(details)), file.path(destination, "Supporting", "permutation_execution"))
    source_files <- c(source_files, details)
  } else {
    audit <- file.path(root, "10_permutation_execution_audit.csv")
    context <- file.path(root, "11_test_context_map.csv")
    if (file.exists(context) && file.exists(audit)) {
      map <- sdn_read(context)
      map <- map[map$Analysis == "Nested", , drop = FALSE]
      sdn_write(sdn_portable_trace(map), file.path(destination, "Supporting", "test_context"))
      execution <- sdn_read(audit)
      execution <- execution[execution$Test_ID %in% map$Test_ID, , drop = FALSE]
      sdn_write(sdn_portable_trace(execution), file.path(destination, "Supporting", "permutation_execution"))
      source_files <- c(source_files, context, audit)
    }
  }
  descriptions <- c(Figure = "Composite figure filename stem", Composite_panel_order = "Order among displayed panels, reset within each figure; NA for omitted comparisons",
    Composite_panel_label = "Lowercase panel tag in the composite; NA for omitted comparisons", Virus = "DNA or RNA viral dataset",
    Factor = "Internal focal-factor key", Factor_label = "Manuscript factor name", Stratum = "Ecological subset; key with Virus and Factor",
    Figure_panel_included = "TRUE if shown; Status requires nominal pooled P < 0.05, without BH filtering",
    Sample = "Sample accession/identifier; unique within each comparison", Group = "Focal-factor level used for colour and geometry",
    Axis1 = "Archived first PCoA coordinate, not recomputed", Axis2 = "Archived second PCoA coordinate, not recomputed",
    CentroidX = "Group mean of Axis1", CentroidY = "Group mean of Axis2", Group_n = "Number of retained samples in the group",
    Project = "Study identifier used for permutation blocks", Sampling_region = "Sampling region", Coral_taxa = "Coral taxon", Compartments = "Sample compartment", Status = "Harmonised sample status",
    PCoA1_percent = "First-axis percentage as used in the archived plot", PCoA2_percent = "Second-axis percentage as used in the archived plot",
    Colour = "Hexadecimal group colour in the archived palette", Group_envelope = "Geometric envelope rule; unrelated to statistical significance",
    N = "Retained samples in the comparison", Levels = "Retained focal-factor levels", Projects = "Represented Projects",
    Projects_within_factor_variation = "Projects containing eligible within-Project focal-factor variation",
    R2_pooled = "Pooled one-factor PERMANOVA R-squared, fraction", F_pooled = "Pooled PERMANOVA pseudo-F", P_pooled = "Unadjusted pooled PERMANOVA P",
    R2_blocked = "One-factor PERMANOVA R-squared; blocking changes permutation inference, not the fitted effect",
    F_blocked = "Project-blocked PERMANOVA pseudo-F", P_blocked = "Unadjusted PERMANOVA P from within-Project permutations",
    Q_pooled = "BH-adjusted pooled PERMANOVA P within Virus x Factor", Q_blocked = "BH-adjusted blocked PERMANOVA P within Virus x Factor",
    Dispersion_F_pooled = "Unblocked PERMDISP F", Dispersion_P_pooled = "Unadjusted unblocked PERMDISP P",
    Dispersion_F_blocked = "Project-blocked PERMDISP F", Dispersion_P_blocked = "Unadjusted Project-blocked PERMDISP P",
    Q_dispersion_pooled = "BH-adjusted pooled PERMDISP P within Virus x Factor", Q_dispersion_blocked = "BH-adjusted blocked PERMDISP P within Virus x Factor",
    Blocked_note = "Eligibility or failure note for blocked inference; NA is not a null result",
    Source_coordinate_file = "Relative archived coordinate filename", Figure_P_source = "Statistic used for displayed significance annotation")
  dictionary <- do.call(rbind, dictionary)
  dictionary$Description <- unname(descriptions[dictionary$Column])
  unlabelled <- is.na(dictionary$Description)
  dictionary$Description[unlabelled] <- ifelse(grepl("_label$", dictionary$Column[unlabelled]),
    "Display string only; use the corresponding numeric P/Q for analysis",
    "Retained archived permutation/provenance field; actual counts may differ from requested counts in finite spaces")
  sdn_write(dictionary, file.path(destination, "data_dictionary"))
  source_files <- unique(source_files)
  provenance <- data.frame(Source_file = vapply(source_files, function(p) {
    if (startsWith(p, paste0(bundle, "/"))) substring(p, nchar(bundle) + 2L)
    else paste0("Analysis_Results/", basename(p))
  }, character(1)), SHA256 = vapply(source_files, digest::digest, character(1), algo = "sha256", file = TRUE))
  sdn_write(provenance, file.path(destination, "source_file_checksums"))
  writeLines(c("# PCoA figure source data", "",
    "Each factor folder contains points and statistics in equivalent UTF-8 CSV and TSV formats. Choose either format; they contain the same values.",
    "Keys are Figure + Virus + Factor + Stratum (and Sample for points). Panel order/tags restart within each composite and include displayed panels only.",
    "Coordinates, centroids, axis percentages, colours and all numerical statistics are read from completed outputs; neither ordinations nor tests are rerun.",
    "All eligible comparisons are retained. Status panels use nominal pooled PERMANOVA P < 0.05; Figure_panel_included is a display flag, not evidence of BH significance.",
    "Numeric P/Q and R-squared values are not rounded for presentation. R-squared values are fractions; axis percentages are percentages. NA denotes unavailable/not estimable, never zero.",
    "BH families are Virus x Factor, separately for pooled/blocked PERMANOVA and PERMDISP. Blocking restricts permutations within Project and does not condition Project sums of squares.",
    "Group envelopes are 80% t ellipses for n >= 4 or convex hulls for n = 3, with spokes to centroids; they are descriptive geometry, not significance regions.",
    "Source_file_checksums identifies the read-only archived inputs. Supporting holds available per-test permutation audit information. No abundance matrices or native checkpoint binaries are duplicated.",
    "", "中文：每个因素文件夹含坐标与完整统计表，CSV/TSV内容相同。Status未入图的比较也予以保留；本次仅导出，不重算PCoA或置换检验。"),
    file.path(destination, "README.md"), useBytes = TRUE)
  files <- list.files(destination, recursive = TRUE, full.names = TRUE)
  checksum <- data.frame(File = substring(files, nchar(destination) + 2L), Bytes = file.info(files)$size,
    SHA256 = vapply(files, digest::digest, character(1), algo = "sha256", file = TRUE))
  readr::write_tsv(checksum, file.path(destination, "file_manifest.tsv"), na = "NA")
  message("Exported PCoA source data: ", destination)
  invisible(destination)
}

# ---- Entry point ----
workflow_help <- function() {
  cat(paste0(
    "Coral virome: ", WORKFLOW_TITLE, "\n\n",
    "Usage: Rscript ", WORKFLOW_SCRIPT, " [options]\n\n",
    "  --input-dir PATH       Four input abundance/metadata files\n",
    "  --output-dir PATH      New, empty results folder (required)\n",
    "  --mode all|analysis|figures|source-data  Default: all\n",
    "  --analysis-dir PATH    Existing bundle or Analysis_Results (figures/source-data)\n",
    if (WORKFLOW_SCOPE == "nested") "  --figure-data-dir PATH  Optional separate bundle with archived plotting coordinates\n" else "",
    "  --workers INTEGER      PSOCK workers for larger tests (default: 4)\n",
    "  --formats pdf,svg,png,tiff   Figure export formats\n",
    if (WORKFLOW_SCOPE == "nested") "  --individual true|false  Export individual PCoA panels (default: true)\n" else "",
    "  --legacy-bundle PATH   Optional read-only native checkpoint replay\n",
    "  --cache-only           Require verified cached tests; never rerun permutations\n",
    "  --help                 Show this help\n\n",
    "Analysis requests 9,999 permutations; finite spaces may be enumerated.\n",
    "DNA and RNA use the same rules and are analysed separately.\n",
    "source-data only organises completed results; no tests or figures are rerun.\n",
    "With RStudio: edit CONFIG, source the script, then run main(character()).\n"
  ))
  invisible(NULL)
}

parse_options <- function(args, config) {
  value_options <- c("input-dir", "output-dir", "mode", "analysis-dir", "workers",
                     "formats", "individual", "legacy-bundle", "figure-data-dir")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) return(NULL)
    if (arg == "--cache-only") { config$cache_only <- TRUE; i <- i + 1L; next }
    if (!startsWith(arg, "--")) stop("Unexpected argument: ", arg)
    name <- substring(arg, 3L)
    if (!name %in% value_options) stop("Unknown option: ", arg)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) stop("Missing value for ", arg)
    config[[gsub("-", "_", name, fixed = TRUE)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  config$mode <- match.arg(config$mode, c("all", "analysis", "figures", "source-data"))
  workers <- suppressWarnings(as.numeric(config$workers))
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers)) {
    stop("workers must be a positive integer")
  }
  config$workers <- as.integer(workers)
  if (is.character(config$formats) && length(config$formats) == 1L) {
    config$formats <- strsplit(config$formats, ",", fixed = TRUE)[[1]]
  }
  config$formats <- unique(tolower(trimws(config$formats)))
  if (!length(config$formats) || length(setdiff(config$formats, c("pdf", "svg", "png", "tiff")))) {
    stop("formats must contain pdf, svg, png and/or tiff")
  }
  if (is.character(config$individual)) {
    if (!tolower(config$individual) %in% c("true", "false")) stop("individual must be true or false")
    config$individual <- tolower(config$individual) == "true"
  }
  config
}

main <- function(args = commandArgs(trailingOnly = TRUE), config = CONFIG) {
  options <- parse_options(args, config)
  if (is.null(options) || (!length(args) && !nzchar(options$output_dir))) return(workflow_help())
  if (!nzchar(options$output_dir)) stop("Set --output-dir to a new, empty folder")
  reuse_results <- options$mode %in% c("figures", "source-data")
  if (reuse_results) {
    if (!nzchar(options$analysis_dir)) stop(options$mode, " mode requires --analysis-dir")
    source_dir <- fig_analysis_dir(options$analysis_dir)
  } else {
    if (!nzchar(options$input_dir)) stop("Set --input-dir to your inputs folder")
    options$input_dir <- normalizePath(options$input_dir, winslash = "/", mustWork = TRUE)
    files <- c("24083_metagenome_TPM.tsv", "rna_virus_tpm.tsv", "sample_info_metagenome.tsv",
               "sample_info_metatranscriptome.tsv", "permutation_plan.csv")
    missing <- files[!file.exists(file.path(options$input_dir, files))]
    if (length(missing)) stop("Missing input files: ", paste(missing, collapse = ", "))
    check_analysis_packages()
  }
  # Validate plotting dependencies before starting an expensive analysis.
  if (options$mode %in% c("all", "figures")) {
    packages <- if (WORKFLOW_SCOPE == "nested") {
      c("ape", "ggplot2", "patchwork", "cowplot", "ggrepel", "stringr", "readr")
    } else c("ggplot2", "patchwork", "scales")
    if (any(c("png", "tiff") %in% options$formats)) packages <- c(packages, "ragg")
    fig_require(packages)
  }
  output <- normalizePath(options$output_dir, winslash = "/", mustWork = FALSE)
  if (file.exists(output) && !dir.exists(output)) stop("Output path is an existing file: ", output)
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Output folder is not empty. Choose a new folder: ", output)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output)) stop("Cannot create output folder: ", output)
  if (!reuse_results) {
    source_dir <- file.path(output, "Data", "Analysis_Results")
    run_analysis(options$input_dir, source_dir, workers = options$workers,
      nperm = 9999L, legacy_bundle = options$legacy_bundle, cache_only = options$cache_only)
  }
  if (options$mode %in% c("all", "figures")) {
    if (WORKFLOW_SCOPE == "nested") {
      render_nested_figures(source_dir, output, formats = options$formats, individual = options$individual)
    } else render_global_figure(source_dir, output, formats = options$formats)
  }
  if (WORKFLOW_SCOPE == "global") {
    sd_global_export(source_dir, output)
  } else if (options$mode != "analysis") {
    # After rendering, export the exact point records just written by that render.
    # Export-only requires an existing figure bundle and never recomputes PCoA.
    figures <- if (options$mode == "source-data") {
      if (nzchar(options$figure_data_dir)) options$figure_data_dir else options$analysis_dir
    } else output
    sd_nested_export(source_dir, output, figure_data_dir = figures)
  }
  provenance <- c(
    paste("Workflow:", WORKFLOW_TITLE), paste("Mode:", options$mode),
    paste("Generated UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Input directory:", options$input_dir), paste("Analysis directory:", source_dir),
    paste("Requested permutations:", if (reuse_results) "not rerun" else "9999"),
    paste("Legacy replay:", if (is.null(options$legacy_bundle)) "none" else options$legacy_bundle),
    "", capture.output(sessionInfo()))
  writeLines(provenance, file.path(output, "run_information.txt"))
  message("Completed ", WORKFLOW_TITLE, ": ", output)
  invisible(output)
}

# Sourcing defines functions and CONFIG without launching an expensive run.
if (sys.nframe() == 0L) main()
