# Shared path/argument rules for the explicit restore and launch commands.
# Installed packages stay inside the unpacked archive unless --library is used.
project_library <- function(root) {
  file.path(root, "environment", "library", paste0("R-", getRversion()), R.version$platform)
}
environment_arguments <- function(args, restore = FALSE) {
  options <- list(library = NULL, bootstrap_archive = NULL, packages = NULL,
                  check = FALSE, forwarded = character())
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    named <- c("--library", if (restore) c("--bootstrap-archive", "--packages"))
    if (key %in% named) {
      if (i == length(args) || !nzchar(args[i + 1L])) stop("Missing value after ", key)
      slot <- switch(key, "--library" = "library", "--bootstrap-archive" = "bootstrap_archive", "--packages" = "packages")
      options[[slot]] <- args[i + 1L]
      i <- i + 2L
    } else if (!restore && identical(key, "--check")) {
      options$check <- TRUE
      i <- i + 1L
    } else if (restore) {
      stop("Unknown setup option: ", key)
    } else {
      options$forwarded <- c(options$forwarded, key)
      i <- i + 1L
    }
  }
  options
}
