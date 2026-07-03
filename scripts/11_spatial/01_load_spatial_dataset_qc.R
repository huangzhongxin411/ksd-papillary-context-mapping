#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

manifest_path <- "results/spatial/phase28a/spatial_dataset_manifest_v0.1.tsv"
out_path <- "results/spatial/phase28a/spatial_load_qc_v0.1.tsv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manifest_path)) {
  stop("Missing manifest: run scripts/11_spatial/00_phase28a_spatial_resource_audit.R first.")
}

manifest <- fread(manifest_path)
ready <- manifest[analysis_ready == TRUE]

if (nrow(ready) == 0) {
  qc <- data.table(
    dataset_id = manifest$dataset_id,
    sample_id = NA_character_,
    n_spots = NA_integer_,
    n_genes_detected_median = NA_real_,
    umi_median = NA_real_,
    percent_mito_median = NA_real_,
    image_loaded = FALSE,
    coordinates_loaded = FALSE,
    qc_status = "not_attempted_not_analysis_ready",
    notes = "No dataset satisfied matrix + spatial coordinates + scalefactors/image + metadata criteria."
  )
  fwrite(qc, out_path, sep = "\t", na = "")
  message("No analysis-ready spatial dataset; wrote not-attempted QC.")
  quit(save = "no", status = 0)
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  qc <- data.table(
    dataset_id = ready$dataset_id,
    sample_id = NA_character_,
    n_spots = NA_integer_,
    n_genes_detected_median = NA_real_,
    umi_median = NA_real_,
    percent_mito_median = NA_real_,
    image_loaded = FALSE,
    coordinates_loaded = FALSE,
    qc_status = "not_attempted_seurat_missing",
    notes = "Seurat package is not installed; cannot call Load10X_Spatial."
  )
  fwrite(qc, out_path, sep = "\t", na = "")
  message("Seurat unavailable; wrote QC blocker.")
  quit(save = "no", status = 0)
}

load_one <- function(ds, local_path) {
  roots <- trimws(unlist(strsplit(local_path, ";", fixed = TRUE)))
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) {
    return(data.table(
      dataset_id = ds, sample_id = NA_character_, n_spots = NA_integer_,
      n_genes_detected_median = NA_real_, umi_median = NA_real_,
      percent_mito_median = NA_real_, image_loaded = FALSE, coordinates_loaded = FALSE,
      qc_status = "load_failed_missing_local_root",
      notes = "No existing local root for analysis-ready manifest row."
    ))
  }
  root <- roots[[1]]
  obj <- tryCatch(Seurat::Load10X_Spatial(data.dir = root), error = function(e) e)
  if (inherits(obj, "error")) {
    return(data.table(
      dataset_id = ds, sample_id = basename(root), n_spots = NA_integer_,
      n_genes_detected_median = NA_real_, umi_median = NA_real_,
      percent_mito_median = NA_real_, image_loaded = FALSE, coordinates_loaded = FALSE,
      qc_status = "load_failed",
      notes = conditionMessage(obj)
    ))
  }
  counts <- obj$nCount_Spatial
  feats <- obj$nFeature_Spatial
  mt <- tryCatch(Seurat::PercentageFeatureSet(obj, pattern = "^MT-"), error = function(e) rep(NA_real_, ncol(obj)))
  data.table(
    dataset_id = ds,
    sample_id = basename(root),
    n_spots = ncol(obj),
    n_genes_detected_median = median(feats, na.rm = TRUE),
    umi_median = median(counts, na.rm = TRUE),
    percent_mito_median = median(mt, na.rm = TRUE),
    image_loaded = length(obj@images) > 0,
    coordinates_loaded = length(obj@images) > 0,
    qc_status = "loaded_pending_manual_review",
    notes = "Loaded with Seurat::Load10X_Spatial; spot-level QC requires manual review before claims."
  )
}

qc <- rbindlist(lapply(seq_len(nrow(ready)), function(i) {
  load_one(ready$dataset_id[[i]], ready$local_path[[i]])
}), fill = TRUE)
fwrite(qc, out_path, sep = "\t", na = "")
message("Spatial load QC complete: ", out_path)
