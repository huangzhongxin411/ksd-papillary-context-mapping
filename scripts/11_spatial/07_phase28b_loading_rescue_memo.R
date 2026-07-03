#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

out_dir <- "results/spatial/phase28b_loading_rescue"
memo_path <- "docs/spatial/phase28b_loading_rescue_completion_memo_v0.1.md"
dir.create(dirname(memo_path), recursive = TRUE, showWarnings = FALSE)

read_or_empty <- function(path) {
  if (file.exists(path)) fread(path) else data.table()
}

dep <- read_or_empty(file.path(out_dir, "spatial_loading_dependency_audit_v0.1.tsv"))
install <- read_or_empty(file.path(out_dir, "hdf5r_install_attempt_v0.1.tsv"))
conv <- read_or_empty(file.path(out_dir, "h5_to_mtx_conversion_qc_v0.1.tsv"))
load <- read_or_empty(file.path(out_dir, "spatial_load_qc_v0.2.tsv"))
gate <- read_or_empty(file.path(out_dir, "spatial_projection_gate_v0.1.tsv"))

hdf5r_installed <- if (nrow(dep) && "hdf5r" %in% dep$package) dep[package == "hdf5r", installed][1] else FALSE
hdf5r_install_success <- if (nrow(install)) install$success[1] else FALSE
fallback_used <- nrow(conv) > 0 && any(conv$conversion_attempted == TRUE | conv$conversion_success == TRUE)
loaded <- if (nrow(load)) load[load_success == TRUE] else data.table()
pass <- if (nrow(load)) load[qc_status == "pass"] else data.table()
projection_allowed <- nrow(gate) > 0 && any(gate$projection_allowed == TRUE)

fmt_loaded <- function(dt) {
  if (!nrow(dt)) return("None.")
  paste(dt[, paste0("- ", sample_id, ": route=", load_route,
                    ", spots=", n_spots_matrix,
                    ", median_features=", round(n_features_median, 1),
                    ", median_counts=", round(n_counts_median, 1),
                    ", image_loaded=", image_loaded,
                    ", coordinates_loaded=", coordinates_loaded,
                    ", qc_status=", qc_status)], collapse = "\n")
}

memo <- c(
  "# Phase 28B-L GSE206306 spatial loading rescue completion memo",
  "",
  "## 1. Was hdf5r installed?",
  "",
  paste0("Initial hdf5r availability: ", hdf5r_installed, ". Install success: ", hdf5r_install_success, "."),
  "",
  "## 2. Was h5-to-mtx fallback used?",
  "",
  paste0(fallback_used, "."),
  "",
  "## 3. Which samples loaded successfully?",
  "",
  fmt_loaded(loaded),
  "",
  "## 4. Were images and coordinates loaded?",
  "",
  if (nrow(loaded)) paste0("Loaded samples with images: ", sum(loaded$image_loaded), "/", nrow(loaded),
                           "; coordinates: ", sum(loaded$coordinates_loaded), "/", nrow(loaded), ".") else "No samples loaded successfully.",
  "",
  "## 5. How many spots and genes were retained?",
  "",
  if (nrow(loaded)) paste(loaded[, paste0("- ", sample_id, ": spots=", n_spots_matrix,
                                          ", median detected genes=", round(n_features_median, 1))], collapse = "\n") else "Not available because no sample loaded successfully.",
  "",
  "## 6. Is spatial projection now allowed?",
  "",
  if (projection_allowed) "Yes for samples passing the projection gate; proceed only to Phase 28C signature projection, without biological claims until projection QC passes." else "No. Spatial projection is not allowed.",
  "",
  "## 7. Should manuscript or Figure 1-5 change?",
  "",
  "No. This was a loading/QC rescue phase only.",
  "",
  "## 8. Manuscript-safe wording",
  "",
  if (projection_allowed) {
    "GSE206306 spatial transcriptomic resources were recovered and loaded for QC; downstream projection was then evaluated separately."
  } else {
    "GSE206306 Visium-like resources were recovered, but spatial object loading remained limited by local software/input compatibility; spatial validation was not used as evidence."
  },
  "",
  "## Claim boundary",
  "",
  "Do not write spatial validation, plaque-specific localization, Loop/TAL spatial enrichment confirmed, or risk-gene spatial localization unless Phase 28C projection is completed and QC-passed."
)
writeLines(memo, memo_path)
message("wrote ", memo_path)
