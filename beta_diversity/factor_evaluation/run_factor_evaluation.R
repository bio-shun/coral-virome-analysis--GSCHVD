# Global four-factor evaluation and summary
# Standalone R workflow. See ../README.md for inputs, commands and methods.
# Original analyses and input files are never overwritten.
WORKFLOW_SCOPE <- "global"
WORKFLOW_SCRIPT <- "run_factor_evaluation.R"
WORKFLOW_TITLE <- "Global four-factor evaluation and summary"

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

safe_partial_dbrda <- function(distance, metadata, focal, engine, project_aware = FALSE) {
    other <- setdiff(factor_order, focal)
    conditions <- if (project_aware) 
        c("Project", other)
    else other
    formula_text <- paste("distance ~", focal, "+ Condition(", paste(conditions, collapse = " + "), ")")
    fit <- tryCatch(vegan::dbrda(as.formula(formula_text), data = metadata), error = function(e) e)
    analysis <- if (project_aware) 
        "Project-conditioned partial dbRDA"
    else "Global partial dbRDA"
    if (inherits(fit, "error") || is.null(fit$CCA) || is.null(fit$CCA$rank) || is.na(fit$CCA$rank) || fit$CCA$rank < 
        1) {
        note <- if (inherits(fit, "error")) 
            conditionMessage(fit)
        else "No independent constrained rank"
        return(data.frame(Analysis = analysis, Factor = focal, Df = 0, R2 = NA_real_, Adj_R2 = NA_real_, F = NA_real_, 
            P = NA_real_, Estimable = FALSE, Note = note))
    }
    permutation <- if (project_aware) 
        make_project_perm(metadata$Project, engine$nperm)
    else engine$nperm
    test <- tryCatch(engine$anova(fit, permutations = permutation), error = function(e) e)
    adj <- tryCatch(vegan::RsquareAdj(fit), error = function(e) NULL)
    r2 <- if (is.null(adj)) 
        NA_real_
    else unname(adj$r.squared)
    adj_r2 <- if (is.null(adj)) 
        NA_real_
    else unname(adj$adj.r.squared)
    if (inherits(test, "error")) {
        if (must_rethrow_engine_error(test, engine)) 
            stop(test)
        return(data.frame(Analysis = analysis, Factor = focal, Df = fit$CCA$rank, R2 = r2, Adj_R2 = adj_r2, F = NA_real_, 
            P = NA_real_, Estimable = TRUE, Note = conditionMessage(test)))
    }
    tab <- as.data.frame(test)
    data.frame(Analysis = analysis, Factor = focal, Df = tab$Df[1], R2 = r2, Adj_R2 = adj_r2, F = tab$F[1], P = tab$`Pr(>F)`[1], 
        Estimable = TRUE, Note = "")
}

project_coverage <- function(metadata, factor_name, virus) {
    counts <- metadata %>% count(Project, .data[[factor_name]], name = "n_level") %>% filter(n_level >= 2)
    eligible <- counts %>% count(Project, name = "Retained_levels") %>% filter(Retained_levels >= 2)
    data.frame(Virus = virus, Factor = factor_name, Eligible_projects = nrow(eligible), Total_projects = nlevels(metadata$Project), 
        Eligible_samples = if (!nrow(eligible)) 
            0
        else sum(counts$n_level[counts$Project %in% eligible$Project]), Coverage = if (!nrow(eligible)) 
            "Not estimable within Project"
        else if (nrow(eligible) < 3) 
            "Case-study coverage only"
        else "Replicated across Projects")
}

project_specific_tests <- function(dataset, factor_name, engine) {
    metadata <- dataset$metadata
    level_counts <- metadata %>% count(Project, .data[[factor_name]], name = "Group_n")
    eligible <- level_counts %>% filter(Group_n >= 2) %>% count(Project, name = "Levels") %>% filter(Levels >= 
        2)
    if (!nrow(eligible)) 
        return(data.frame())
    rows <- lapply(as.character(eligible$Project), function(project_id) {
        valid_levels <- level_counts %>% filter(Project == project_id, Group_n >= 2) %>% pull(!!rlang::sym(factor_name)) %>% 
            as.character()
        keep <- metadata$Project == project_id & as.character(metadata[[factor_name]]) %in% valid_levels
        meta_sub <- droplevels(metadata[keep, , drop = FALSE])
        dist_sub <- as.dist(as.matrix(dataset$distance)[rownames(meta_sub), rownames(meta_sub), drop = FALSE])
        engine$set_context("Project-specific", factor_name, project = project_id)
        ad <- safe_adonis(dist_sub, meta_sub, factor_name, engine, blocked = FALSE)
        disp <- safe_dispersion(dist_sub, meta_sub[[factor_name]], engine)
        data.frame(Virus = dataset$virus, Factor = factor_name, Project = project_id, N = nrow(meta_sub), Levels = nlevels(meta_sub[[factor_name]]), 
            R2 = ad$R2, F = ad$F, P = ad$P, Dispersion_F = as.numeric(disp[["F"]]), Dispersion_P = as.numeric(disp[["P"]]), 
            Note = ad$Note, stringsAsFactors = FALSE)
    })
    bind_rows(rows)
}

analyse_global <- function(dataset, engine) {
    virus <- dataset$virus
    metadata <- dataset$metadata
    engine$set_context("Global marginal PERMANOVA")
    global_fit <- engine$adonis(dataset$distance ~ Coral_taxa + Sampling_region + Compartments + Status, data = metadata, 
        by = "margin", permutations = engine$nperm)
    global_stats <- as.data.frame(global_fit)
    global_stats$Factor <- rownames(global_stats)
    rownames(global_stats) <- NULL
    global_stats <- global_stats %>% filter(Factor %in% factor_order) %>% transmute(Virus = virus, Analysis = "Global marginal PERMANOVA", 
        Factor, Df, SumOfSqs, R2, F, P = `Pr(>F)`) %>% mutate(Q = p.adjust(P, method = "BH"))
    global_partial <- bind_rows(lapply(factor_order, function(f) {
        engine$set_context("Global partial dbRDA", f)
        safe_partial_dbrda(dataset$distance, metadata, f, engine, FALSE)
    })) %>% mutate(Virus = virus, .before = 1)
    project_partial <- bind_rows(lapply(factor_order, function(f) {
        engine$set_context("Project-conditioned partial dbRDA", f)
        safe_partial_dbrda(dataset$distance, metadata, f, engine, TRUE)
    })) %>% mutate(Virus = virus, .before = 1)
    dispersion <- bind_rows(lapply(factor_order, function(f) {
        engine$set_context("Global PERMDISP", f)
        x <- safe_dispersion(dataset$distance, metadata[[f]], engine)
        data.frame(Virus = virus, Factor = f, Dispersion_F = as.numeric(x[["F"]]), Dispersion_P = as.numeric(x[["P"]]), 
            Note = unname(x[["Note"]]))
    })) %>% mutate(Dispersion_Q = p.adjust(Dispersion_P, method = "BH"))
    coverage <- bind_rows(lapply(factor_order, function(f) project_coverage(metadata, f, virus)))
    project_partial <- project_partial %>% left_join(coverage %>% select(Virus, Factor, Eligible_projects, Total_projects, 
        Coverage), by = c("Virus", "Factor")) %>% mutate(Model_estimable = Estimable, Estimable = Estimable & Eligible_projects > 
        0, Note = ifelse(Eligible_projects == 0, "No Project with >=2 samples per retained level across >=2 levels", 
        Note))
    per_project <- bind_rows(lapply(factor_order, function(f) project_specific_tests(dataset, f, engine)))
    if (nrow(per_project)) 
        per_project <- per_project %>% group_by(Virus, Factor) %>% mutate(Q = p.adjust(P, method = "BH"), Dispersion_Q = p.adjust(Dispersion_P, 
            method = "BH")) %>% ungroup()
    list(global_stats = global_stats, global_partial = global_partial, project_partial = project_partial, dispersion = dispersion, 
        coverage = coverage, per_project = per_project)
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
    result <- analyse_global(dataset, engine)
    result$sample_summary <- dataset_summary(dataset)
    results[[virus]] <- result
    saveRDS(result, file.path(output_dir, paste0(virus, "_completed_analysis.rds")))
    rm(dataset)
    invisible(gc(FALSE))
  }
  tables <- c(global_stats = "01_global_marginal_permanova.csv", global_partial = "02_global_partial_dbrda.csv", project_partial = "03_project_conditioned_partial_dbrda.csv", dispersion = "04_global_permdisp.csv", coverage = "05_within_project_coverage.csv", per_project = "06_project_specific_permanova_permdisp.csv", sample_summary = "08_sample_summary.csv")
  for (name in names(tables)) {
    values <- lapply(results, function(x) x[[name]])
    write.csv(dplyr::bind_rows(values), file.path(output_dir, tables[[name]]), row.names = FALSE, na = "")
  }
  
  engine$flush()
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "09_sessionInfo.txt"))
  writeLines(c(
    "DNA and RNA are analysed separately using square-root TPM and Bray-Curtis distances.",
    paste0("Requested permutations: ", nperm, ". Native actual counts and finite-space enumeration are in 10_permutation_execution_audit.csv."),
    "Per-test identities and ecological strata are linked in 11_test_context_map.csv.",
    "Global marginal PERMANOVA tests each factor while accounting for the other three.",
    "Partial dbRDA conditions on the other factors; Project-conditioned models additionally condition on Project and use within-Project permutations.",
    "Project coverage and Project-specific tests require at least two retained levels with at least two samples per level.",
    "Project-specific BH correction is within Virus x Factor; global marginal and PERMDISP corrections are within Virus.",
    "PERMDISP uses distances to spatial medians.",
    "Sample libraries are not automatically independent biological replicates.",
    "Recorded per-test seeds reconstruct the publication permutation matrices; assignment SHA-256 values are checked before every newly computed test.",
    "Optional verified legacy replay preserves archived native tests without writing to the legacy folder."
  ), file.path(output_dir, "README_method_and_interpretation.txt"))
  message("Analysis completed: ", output_dir)
  invisible(results)
}

# ---- Figure functions ----
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

render_global_figure <- function(analysis_root, output_root, formats = c("pdf", "svg", "tiff", "png")) {
    fig_require(c("ggplot2", "patchwork", "scales", if (any(c("tiff", "png") %in% formats)) "ragg"))
    root <- fig_analysis_dir(analysis_root)
    order <- c("Coral_taxa", "Sampling_region", "Compartments", "Status")
    labels <- fig_factor_labels
    colours <- c(Coral_taxa = "#3B6FB6", Sampling_region = "#2A9D8F", Compartments = "#E9A23B", Status = "#C65A64")
    global <- fig_read_table(root, "01_global_marginal_permanova.csv", c("Virus", "Factor", "R2"))
    partial <- fig_read_table(root, "03_project_conditioned_partial_dbrda.csv", c("Virus", "Factor", "Adj_R2", 
        "Estimable", "Eligible_projects"))
    coverage <- fig_read_table(root, "05_within_project_coverage.csv", c("Virus", "Factor", "Eligible_projects", 
        "Total_projects"))
    projects <- fig_read_table(root, "06_project_specific_permanova_permdisp.csv", c("Virus", "Factor", "R2"))
    theme <- ggplot2::theme_classic(base_size = 8.5, base_family = "Arial") + ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.35, 
        colour = "black"), axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"), axis.text = ggplot2::element_text(colour = "black"), 
        panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"), 
        plot.title = ggplot2::element_text(face = "bold", size = 9), plot.subtitle = ggplot2::element_text(size = 7, 
            colour = "#444444"), legend.position = "none")
    global$Factor <- factor(global$Factor, levels = rev(order))
    global$Label <- sprintf("%.1f%%", 100 * global$R2)
    a <- ggplot2::ggplot(global, ggplot2::aes(Factor, 100 * R2, fill = Factor)) + ggplot2::geom_col(width = 0.68, 
        colour = "black", linewidth = 0.2) + ggplot2::geom_text(ggplot2::aes(label = Label), hjust = -0.1, size = 2.3) + 
        ggplot2::coord_flip() + ggplot2::facet_wrap(~Virus, ncol = 1, scales = "free_y") + ggplot2::scale_fill_manual(values = colours, 
        labels = labels) + ggplot2::scale_x_discrete(labels = labels) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
        0.2))) + ggplot2::labs(title = "Global pooled association", subtitle = "Marginal PERMANOVA; sqrt(TPM) Bray-Curtis", 
        x = NULL, y = "Variation explained (%)") + theme
    partial$Factor <- factor(partial$Factor, levels = rev(order))
    partial$Value <- ifelse(!partial$Estimable | is.na(partial$Adj_R2), 0, pmax(0, 100 * partial$Adj_R2))
    partial$Label <- ifelse(!partial$Estimable | is.na(partial$Adj_R2), "not estimable", sprintf("%.2f%% (%s)", 
        100 * partial$Adj_R2, paste0(partial$Eligible_projects, ifelse(partial$Eligible_projects == 1, " Project", 
            " Projects"))))
    b <- ggplot2::ggplot(partial, ggplot2::aes(Factor, Value, fill = Factor)) + ggplot2::geom_col(width = 0.68, 
        colour = "black", linewidth = 0.2) + ggplot2::geom_text(ggplot2::aes(label = Label), hjust = -0.05, size = 2.1) + 
        ggplot2::coord_flip() + ggplot2::facet_wrap(~Virus, ncol = 1, scales = "free_y") + ggplot2::scale_fill_manual(values = colours, 
        labels = labels) + ggplot2::scale_x_discrete(labels = labels) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
        0.55))) + ggplot2::labs(title = "Project-conditioned association", subtitle = "Adjusted partial R-squared; within-Project permutations", 
        x = NULL, y = "Adjusted unique variation (%)") + theme
    coverage$Factor <- factor(coverage$Factor, levels = order)
    coverage$Fraction <- coverage$Eligible_projects/coverage$Total_projects
    coverage$Label <- paste0(coverage$Eligible_projects, "/", coverage$Total_projects)
    c <- ggplot2::ggplot(coverage, ggplot2::aes(Factor, Fraction, fill = Factor)) + ggplot2::geom_col(width = 0.68, 
        colour = "black", linewidth = 0.2) + ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -0.35, size = 2.3) + 
        ggplot2::facet_wrap(~Virus, ncol = 1) + ggplot2::scale_fill_manual(values = colours, labels = labels) + 
        ggplot2::scale_x_discrete(labels = labels) + ggplot2::scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) + 
        ggplot2::labs(title = "Within-Project evidence coverage", subtitle = "Projects with >=2 replicated levels", 
            x = NULL, y = "Eligible Projects") + theme
    if (nrow(projects)) {
        projects$Factor <- factor(projects$Factor, levels = order)
        d <- ggplot2::ggplot(projects, ggplot2::aes(Factor, 100 * R2, colour = Factor)) + ggplot2::geom_boxplot(outlier.shape = NA, 
            width = 0.55, colour = "#555555", linewidth = 0.3) + ggplot2::geom_jitter(size = 1.1, alpha = 0.75, 
            position = ggplot2::position_jitter(width = 0.12, seed = 20260903)) + ggplot2::facet_wrap(~Virus, ncol = 1, 
            scales = "free_y") + ggplot2::scale_colour_manual(values = colours, labels = labels) + ggplot2::scale_x_discrete(labels = labels) + 
            ggplot2::labs(title = "Project-specific association", subtitle = "Each point is one eligible Project", 
                x = NULL, y = "PERMANOVA R-squared (%)") + theme
    }
    else d <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No eligible Project-specific comparisons")
    plot <- patchwork::wrap_plots(a, b, c, d, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
    files <- fig_save(plot, file.path(output_root, "Figures", "Global_summary", "Figure_Global_Four_Factor_Associations"), 
        11.5, 8.2, formats, family = "Arial")
    message("Global four-factor summary completed.")
    invisible(list(plot = plot, files = files))
}

# ---- Organised figure source data ----
# Portable Source Data export for the saved global four-factor figure.
# This file defines functions only. No model fitting, permutations or plotting
# are invoked. All original statistical columns are retained without rounding.

sdg_root <- function(analysis_root) {
  root <- normalizePath(analysis_root, winslash = "/", mustWork = TRUE)
  nested <- file.path(root, "Data", "Analysis_Results")
  if (dir.exists(nested)) nested else root
}

sdg_read <- function(root, filename, required = character()) {
  path <- file.path(root, filename)
  if (!file.exists(path)) stop("Missing saved analysis table: ", filename)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"))
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(filename, " is missing: ", paste(absent, collapse = ", "))
  x
}

sdg_write_pair <- function(x, stem) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  # Keep missing values explicit in source-data files; this distinguishes NA
  # from an intentionally empty display label (for example Panel D).
  readr::write_csv(x, paste0(stem, ".csv"), na = "NA")
  readr::write_tsv(x, paste0(stem, ".tsv"), na = "NA")
  c(paste0(stem, ".csv"), paste0(stem, ".tsv"))
}

sdg_safe_trace <- function(x) {
  # Checkpoint provenance can refer to a private archive outside this bundle.
  # Keep its portable basename and hash identity, never an absolute local path.
  for (column in names(x)) {
    if (!is.character(x[[column]])) next
    if (grepl("path|checkpoint|directory|file", column, ignore.case = TRUE)) {
      absolute <- !is.na(x[[column]]) & grepl("^(?:[A-Za-z]:[/\\\\]|[/\\\\]{1,2})", x[[column]], perl = TRUE)
      x[[column]][absolute] <- basename(gsub("\\\\", "/", x[[column]][absolute]))
    }
  }
  x
}

sdg_annotate <- function(x, filename, panel, labels) {
  if (any(c("Figure", "Panel", "Factor_label", "Source_table", "Source_row") %in% names(x))) {
    stop("Saved statistical table already contains reserved Source Data columns: ", filename)
  }
  x$Figure <- rep("Figure_Global_Four_Factor_Associations", nrow(x))
  x$Panel <- rep(panel, nrow(x))
  x$Factor_label <- if ("Factor" %in% names(x)) unname(labels[as.character(x$Factor)]) else rep(NA_character_, nrow(x))
  if ("Factor" %in% names(x) && anyNA(x$Factor_label)) stop("Unrecognised factor in ", filename)
  x$Source_table <- rep(filename, nrow(x))
  x$Source_row <- seq_len(nrow(x))
  first <- c("Figure", "Panel", intersect("Virus", names(x)), intersect("Factor", names(x)), "Factor_label")
  x[, c(first, setdiff(names(x), first)), drop = FALSE]
}

sdg_factor_display <- function(x, factor_order, colours) {
  if (!"Factor" %in% names(x)) return(x)
  f <- as.character(x$Factor)
  x$Factor_order <- unname(match(f, factor_order))
  x$Colour <- unname(colours[f])
  if (anyNA(x$Factor_order[!is.na(f)])) {
    stop("Unrecognised factor in display metadata: ",
         paste(unique(f[is.na(x$Factor_order) & !is.na(f)]), collapse = ", "))
  }
  first <- c("Figure", "Panel", intersect("Virus", names(x)),
             intersect("Factor", names(x)), intersect("Factor_label", names(x)),
             "Factor_order", "Colour")
  x[, c(first, setdiff(names(x), first)), drop = FALSE]
}

sdg_dictionary <- function(tables) {
  descriptions <- c(
    Figure = "Composite figure filename stem; no extension.",
    Panel = "Composite panel A-D, or Supporting for saved results not directly plotted.",
    Virus = "DNA and RNA communities are analysed separately.",
    Factor = "Original machine-readable ecological factor identifier.",
    Factor_label = "Human-readable factor name, without replacing the original Factor.",
    Factor_order = "Display order (1 Coral_taxa, 2 Sampling_region, 3 Compartments, 4 Status): top to bottom for horizontal Panels A/B and left to right for Panels C/D. Panels A/B reverse the internal factor levels before coord_flip to achieve this order.",
    Colour = "Factor colour used by the renderer: Coral_taxa #3B6FB6, Sampling_region #2A9D8F, Compartments #E9A23B, Status #C65A64.",
    Analysis = "Saved analysis/model label.",
    Project = "Study identifier; each Panel D row is one eligible Project-factor comparison.",
    N = "Number of analysed samples in the corresponding saved comparison.",
    Levels = "Number of retained levels of the focal factor.",
    Df = "Saved degrees of freedom for the focal term.",
    SumOfSqs = "Saved distance-based sum of squares for the focal term.",
    R2 = "Saved unadjusted R-squared, on its original 0-1 scale; not rounded.",
    Adj_R2 = "Saved adjusted partial R-squared on its original scale; negative/nonestimable values are preserved.",
    F = "Saved PERMANOVA pseudo-F or dbRDA F statistic, as identified by Analysis.",
    P = "Saved unrounded permutation P value; no new test or threshold substitution.",
    Q = "Saved Benjamini-Hochberg-adjusted P value; no new correction is performed.",
    Dispersion_F = "Saved PERMDISP F statistic using distances to group spatial medians.",
    Dispersion_P = "Saved unrounded PERMDISP permutation P value.",
    Dispersion_Q = "Saved Benjamini-Hochberg-adjusted PERMDISP P value.",
    Estimable = "Saved analysis reporting flag; preserved exactly and used by the existing Panel B renderer.",
    Model_estimable = "Saved model-rank estimability flag, distinct from study-coverage/reporting estimability; NA if absent from the source table.",
    Note = "Saved analysis or test-status note, including reasons a value is unavailable.",
    Eligible_projects = "Projects containing at least two retained focal-factor levels, each represented by at least two samples.",
    Total_projects = "Total study Projects in the corresponding filtered viral dataset.",
    Eligible_samples = "Samples in retained levels of eligible within-Project comparisons.",
    Coverage = "Saved within-Project evidence-coverage category.",
    Final_samples = "Total samples retained after the saved filtering rules.",
    Viral_features = "Number of nonzero viral features retained by the saved workflow.",
    Projects = "Number of study Projects in the filtered dataset.",
    Sampling_region_levels = "Number of retained Sampling region levels.",
    Coral_taxa_levels = "Number of retained Coral taxa levels.",
    Compartments_levels = "Number of retained Compartments levels.",
    Status_levels = "Number of retained Status levels.",
    Source_table = "Portable filename of the original saved analysis table.",
    Source_row = "One-based data-row index in Source_table, excluding its header.",
    Plotted_value_percent = "Panel A/D: 100*R2. Panel B: 0 when !Estimable or Adj_R2 is missing, otherwise max(0,100*Adj_R2), matching the renderer. This is the data value before decorative point jitter.",
    Plotted_value_fraction = "Panel C bar height: Eligible_projects/Total_projects, before percent axis formatting.",
    Raw_Adj_R2_percent = "100*Adj_R2 without clamping, retaining negative values and missing values.",
    Plotted_label = "Exact displayed data label used by the renderer; empty for Panel D points, which have no point-specific text label.",
    Plotted_value_rule = "Explicit transformation used for this panel's display data.",
    Test_ID = "Saved native permutation-test identifier; links execution audit and context map.",
    Type = "Saved native test type: PERMANOVA, partial dbRDA or PERMDISP.",
    Requested_permutations = "Requested number of permutations for the saved native test.",
    Actual_permutations = "Actual number of nonidentity permutations executed; finite spaces can be enumerated.",
    Complete_enumeration = "Saved flag indicating exhaustive finite-space permutation enumeration.",
    Minimum_P_resolution = "Saved numerical permutation P resolution, based on the actual execution count.",
    Seed = "Saved per-test random seed.",
    Workers = "Saved number of execution workers.",
    Elapsed_seconds = "Saved native-test elapsed runtime in seconds.",
    Model = "Saved model text, rather than a newly fitted model.",
    Checkpoint = "Saved relative checkpoint reference, or basename when an absolute local path was recorded; checkpoint files are not copied.",
    Cache_origin = "Saved origin of the native result, for example computed or verified archive replay.",
    Stratum = "Saved ecological-stratum context, when applicable.",
    Audit_scope = "Whether trace records are explicitly selected by global context or are the complete legacy audit without inferred test mappings.",
    File = "Relative file path inside this Source_Data folder.",
    Bytes = "File size in bytes.",
    SHA256 = "SHA-256 checksum of the file bytes.",
    Rows = "Number of data rows; NA for text documentation.",
    Columns = "Number of data columns; NA for text documentation."
  )
  rows <- lapply(names(tables), function(name) {
    x <- tables[[name]]
    fields <- names(x)
    description <- unname(descriptions[fields])
    description[is.na(description)] <- "Additional original saved field; retained without modification."
    data.frame(Table = name, Field = fields,
      Type = vapply(x, function(z) paste(class(z), collapse = "/"), character(1)),
      Description = description, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

sd_global_export <- function(analysis_root, output_root) {
  packages <- c("readr", "digest")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Source Data export requires: ", paste(missing, collapse = ", "))
  root <- sdg_root(analysis_root)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  dest <- file.path(output_root, "Source_Data")
  if (dir.exists(dest) && length(list.files(dest, all.files = TRUE, no.. = TRUE))) {
    stop("Source_Data is not empty; choose a new output folder: ", dest)
  }
  labels <- c(Coral_taxa = "Coral taxa", Sampling_region = "Sampling region",
              Compartments = "Compartments", Status = "Status")
  factor_order <- names(labels)
  colours <- c(Coral_taxa = "#3B6FB6", Sampling_region = "#2A9D8F",
               Compartments = "#E9A23B", Status = "#C65A64")
  specs <- list(
    Panel_A_global_marginal = list(file = "01_global_marginal_permanova.csv", panel = "A",
      required = c("Virus", "Factor", "R2")),
    Panel_B_project_conditioned = list(file = "03_project_conditioned_partial_dbrda.csv", panel = "B",
      required = c("Virus", "Factor", "Adj_R2", "Estimable", "Eligible_projects")),
    Panel_C_within_project_coverage = list(file = "05_within_project_coverage.csv", panel = "C",
      required = c("Virus", "Factor", "Eligible_projects", "Total_projects")),
    Panel_D_project_specific = list(file = "06_project_specific_permanova_permdisp.csv", panel = "D",
      required = c("Virus", "Factor", "Project", "R2")))
  # Read and validate all figure tables before creating any output files.
  tables <- lapply(specs, function(s) sdg_factor_display(
    sdg_annotate(sdg_read(root, s$file, s$required), s$file, s$panel, labels),
    factor_order, colours))
  a <- tables$Panel_A_global_marginal
  a$Plotted_value_percent <- 100 * a$R2
  a$Plotted_label <- sprintf("%.1f%%", a$Plotted_value_percent)
  a$Plotted_value_rule <- "100 * R2"
  tables$Panel_A_global_marginal <- a
  b <- tables$Panel_B_project_conditioned
  if (!"Model_estimable" %in% names(b)) b$Model_estimable <- NA
  b$Raw_Adj_R2_percent <- 100 * b$Adj_R2
  b$Plotted_value_percent <- ifelse(!b$Estimable | is.na(b$Adj_R2), 0, pmax(0, 100 * b$Adj_R2))
  b$Plotted_label <- ifelse(!b$Estimable | is.na(b$Adj_R2), "not estimable",
    sprintf("%.2f%% (%s)", 100 * b$Adj_R2,
      paste0(b$Eligible_projects, ifelse(b$Eligible_projects == 1, " Project", " Projects"))))
  b$Plotted_value_rule <- "ifelse(!Estimable | is.na(Adj_R2), 0, pmax(0, 100 * Adj_R2))"
  tables$Panel_B_project_conditioned <- b
  c <- tables$Panel_C_within_project_coverage
  c$Plotted_value_fraction <- c$Eligible_projects / c$Total_projects
  c$Plotted_value_percent <- 100 * c$Plotted_value_fraction
  c$Plotted_label <- paste0(c$Eligible_projects, "/", c$Total_projects)
  c$Plotted_value_rule <- "Eligible_projects / Total_projects (fraction; axis formatted as percent)"
  tables$Panel_C_within_project_coverage <- c
  d <- tables$Panel_D_project_specific
  d$Plotted_value_percent <- 100 * d$R2
  d$Plotted_label <- rep("", nrow(d))
  d$Plotted_value_rule <- rep("100 * R2 (point statistic before display jitter; boxplots summarize these values)", nrow(d))
  tables$Panel_D_project_specific <- d
  supporting <- c("02_global_partial_dbrda.csv", "04_global_permdisp.csv", "08_sample_summary.csv")
  missing_support <- supporting[!file.exists(file.path(root, supporting))]
  for (filename in setdiff(supporting, missing_support)) {
    name <- file.path("Supporting", sub("\\.csv$", "", filename))
    tables[[name]] <- sdg_annotate(sdg_read(root, filename), filename, "Supporting", labels)
  }
  context_file <- "11_test_context_map.csv"
  audit_file <- "10_permutation_execution_audit.csv"
  context <- NULL
  if (file.exists(file.path(root, context_file))) {
    context <- sdg_read(root, context_file, c("Analysis", "Test_ID"))
    keep <- context$Analysis %in% c("Global marginal PERMANOVA", "Global partial dbRDA",
      "Project-conditioned partial dbRDA", "Global PERMDISP", "Project-specific")
    context <- context[keep, , drop = FALSE]
    context$Audit_scope <- rep("Explicit global analysis contexts from the saved context map", nrow(context))
    tables[["Supporting/test_context_map"]] <- sdg_safe_trace(context)
  }
  if (file.exists(file.path(root, audit_file))) {
    audit <- sdg_read(root, audit_file, "Test_ID")
    if (!is.null(context)) audit <- audit[audit$Test_ID %in% context$Test_ID, , drop = FALSE]
    audit$Audit_scope <- rep(if (is.null(context))
      "Complete saved legacy audit; context map absent, so ambiguous test-to-panel mappings are not inferred"
      else "Tests selected using the explicit global context map", nrow(audit))
    tables[["Supporting/permutation_execution_audit"]] <- sdg_safe_trace(audit)
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  for (name in names(tables)) files <- c(files, sdg_write_pair(tables[[name]], file.path(dest, name)))
  dictionary <- sdg_dictionary(tables)
  files <- c(files, sdg_write_pair(dictionary, file.path(dest, "data_dictionary")))
  readme <- c(
    "# Source Data: global four-factor evaluation", "",
    "These files contain saved results and the numeric values used by Figure_Global_Four_Factor_Associations.",
    "The exporter reads completed analysis tables only. It does not fit models, run permutations, correct P values or draw figures.", "",
    "- Panel_A_global_marginal: all saved marginal PERMANOVA terms and 100*R2 bar heights.",
    "- Panel_B_project_conditioned: all saved Project-conditioned partial dbRDA terms, raw adjusted R2 and display values/labels.",
    "- Panel_C_within_project_coverage: Project eligibility counts, denominators and displayed fractions.",
    "- Panel_D_project_specific: every eligible saved Project-specific result, including nonsignificant results; 100*R2 are the values summarized by the boxes and points.", "",
    "CSV and TSV are parallel exports of the same table. Numeric P/Q and effect-size columns are not rounded; scientific notation denotes a number, not a threshold. Missing values are written as NA. An intentionally empty Plotted_label (Panel D) means that the renderer has no point-specific text label; it is not a missing statistic.",
    "Panel B preserves negative Adj_R2 and separate Estimable/Model_estimable flags. As in the renderer, the bar height is zero for nonestimable/missing values and is clamped at zero for a negative adjusted R2. The saved raw result is never replaced by that display value.",
    "Panel D display jitter is decorative and is not included in the exported R2-derived values; the source statistic remains exact. Boxplots are generated from these point statistics using the unchanged renderer.",
    "Source_table and Source_row link each exported result to its original saved table. Original statistical columns are retained; the new Figure, Panel, Factor_label, Factor_order, Colour and Plotted_* columns document figure mapping.",
    "Factor_order is top-to-bottom in horizontal Panels A/B and left-to-right in Panels C/D. Internal factor levels are reversed in Panels A/B before coord_flip, matching the renderer. Colour records the exact factor palette.", "",
    "Saved Benjamini-Hochberg families: Panel A Q was adjusted across the four ecological factors separately within DNA and RNA. Supporting global PERMDISP Dispersion_Q was adjusted across the four factors separately within each viral dataset. Panel D Q and Dispersion_Q were adjusted across eligible Project-specific tests separately within each Virus x Factor combination, with PERMANOVA and PERMDISP corrected as separate families. The partial dbRDA source tables do not contain a Q column; no adjusted value is inferred or calculated by this export.", "",
    "Supporting/ contains global partial dbRDA, global PERMDISP and sample totals when available. Execution audit and context records are included when available; absolute checkpoint paths are reduced to basenames. Checkpoint payloads and original input matrices are not copied.",
    "An explicit saved test-context map is used to select global audit records. If this map is absent in a legacy bundle, the complete saved audit is retained with an explicit scope note; no ambiguous mapping is inferred.",
    "Missing supporting tables: ",
    if (length(missing_support)) paste(missing_support, collapse = ", ") else "None.", "",
    "Data dictionary: data_dictionary.csv / data_dictionary.tsv.",
    "SHA256_manifest.csv / .tsv list every generated file except the manifests themselves, with relative paths and SHA-256 hashes.",
    "Raw table provenance is in source_table_checksums.csv / .tsv; it records original filenames and hashes without private directory paths.", "",
    "Project is the study/blocking variable. The four ecological factors remain Sampling region, Coral taxa, Compartments and Status.",
    "The exporter preserves the existing replication, exchangeability and modelling assumptions; it does not validate biological independence or remove study confounding."
  )
  readme_path <- file.path(dest, "README.md")
  writeLines(readme, readme_path, useBytes = TRUE)
  files <- c(files, readme_path)
  input_files <- unique(c(vapply(specs, `[[`, character(1), "file"), setdiff(supporting, missing_support),
    intersect(c(context_file, audit_file), list.files(root))))
  source_checksums <- data.frame(Source_table = input_files,
    Bytes = as.numeric(file.info(file.path(root, input_files))$size),
    SHA256 = vapply(file.path(root, input_files), function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    stringsAsFactors = FALSE)
  files <- c(files, sdg_write_pair(source_checksums, file.path(dest, "source_table_checksums")))
  manifest <- data.frame(
    File = substring(gsub("\\\\", "/", files), nchar(gsub("\\\\", "/", dest)) + 2L),
    Bytes = as.numeric(file.info(files)$size),
    SHA256 = vapply(files, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    stringsAsFactors = FALSE)
  files <- c(files, sdg_write_pair(manifest, file.path(dest, "SHA256_manifest")))
  message("Global Source Data exported: ", length(tables), " tables; no statistical tests rerun.")
  invisible(list(directory = dest, tables = tables, files = files, manifest = manifest,
                  dictionary = dictionary, source_checksums = source_checksums))
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
