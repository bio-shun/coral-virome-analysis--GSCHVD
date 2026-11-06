#!/usr/bin/env Rscript
# Install the exact dependency versions into a project-local renv library.
# Run from any working directory: Rscript environment/restore_environment.R
# Internet access is needed; source-only packages may need Rtools44 on Windows.
# Optional: --library PATH, --bootstrap-archive renv_1.1.5.tar.gz,
#           --packages NAME,NAME (repair only these locked packages).
# This explicit setup never writes to a personal or system R package library.
flag <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(flag)) dirname(dirname(normalizePath(sub("^--file=", "", flag[1]), winslash = "/"))) else getwd()
if (as.character(getRversion()) != "4.4.1") stop("Install R 4.4.1 before restoring this release.")
source(file.path(root, "environment", "project_library.R"))
options <- environment_arguments(commandArgs(trailingOnly = TRUE), restore = TRUE)
library <- if (is.null(options$library)) project_library(root) else options$library
dir.create(library, recursive = TRUE, showWarnings = FALSE)
library <- normalizePath(library, winslash = "/", mustWork = TRUE)
bootstrap <- file.path(root, "environment", "bootstrap", paste0("R-", getRversion()), R.version$platform)
dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
# renv's user cache and sandbox would otherwise write outside the project.
Sys.setenv(RENV_CONFIG_CACHE_ENABLED = "FALSE", RENV_CONFIG_SANDBOX_ENABLED = "FALSE",
           RENV_PATHS_ROOT = file.path(root, "environment", ".renv"))
record <- read.delim(file.path(root, "environment", "bootstrap_package.tsv"), stringsAsFactors = FALSE)
present <- tryCatch(as.character(packageVersion("renv", lib.loc = bootstrap)), error = function(e) "")
if (!identical(present, record$Version)) {
  archive <- options$bootstrap_archive
  if (is.null(archive)) {
    archive <- tempfile(fileext = ".tar.gz")
    download.file(record$Source_URL, archive, mode = "wb", quiet = FALSE)
  }
  if (!file.exists(archive) || tolower(unname(tools::md5sum(archive))) != record$MD5)
    stop("The renv bootstrap archive is missing or failed its recorded checksum.")
  install.packages(archive, lib = bootstrap, repos = NULL, type = "source")
}
if (as.character(packageVersion("renv", lib.loc = bootstrap)) != record$Version)
  stop("Could not install the pinned project-local renv bootstrap.")
.libPaths(c(bootstrap, .Library), include.site = FALSE)
packages <- if (is.null(options$packages)) NULL else strsplit(options$packages, ",", fixed = TRUE)[[1L]]
renv::restore(project = root, library = library, lockfile = file.path(root, "renv.lock"),
              packages = packages, prompt = FALSE)
if (is.null(packages)) {
  cat("Restore completed. Set R_LIBS to this library before running the module scripts.\n")
} else {
  cat("Partial repair completed for: ", paste(packages, collapse = ", "),
      ". Check the complete environment before analysis.\n", sep = "")
}
