#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

dir.create("results/spatial/phase28b_loading_rescue", recursive = TRUE, showWarnings = FALSE)

attempted <- FALSE
success <- requireNamespace("hdf5r", quietly = TRUE)
err <- ""

if (!success) {
  attempted <- TRUE
  res <- tryCatch({
    install.packages("hdf5r", repos = "https://cloud.r-project.org")
    TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    FALSE
  }, warning = function(w) {
    err <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  success <- requireNamespace("hdf5r", quietly = TRUE)
  if (!isTRUE(res) && !success && !nzchar(err)) err <- "install.packages returned non-TRUE result"
}

if (success) {
  suppressPackageStartupMessages(library(hdf5r))
}

out <- data.table(
  attempted = attempted,
  success = success,
  error_message = err,
  fallback_required = !success,
  notes = ifelse(success, "hdf5r is available; direct Seurat HDF5 route can be attempted.",
                 "hdf5r unavailable; proceed to Python h5-to-mtx fallback.")
)
fwrite(out, "results/spatial/phase28b_loading_rescue/hdf5r_install_attempt_v0.1.tsv", sep = "\t", na = "")
message("hdf5r_status=", ifelse(success, "available", "unavailable"))
