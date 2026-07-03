#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("data/raw/spatial", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/spatial", recursive = TRUE, showWarnings = FALSE)
dir.create("results/spatial/phase28a", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/spatial", recursive = TRUE, showWarnings = FALSE)
dir.create("scripts/11_spatial", recursive = TRUE, showWarnings = FALSE)

out_dir <- "results/spatial/phase28a"

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", na = "")
}

file_size_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  round(file.info(path)$size / 1024^2, 3)
}

first_line <- function(path) {
  if (!file.exists(path)) return("")
  if (file.info(path)$size > 5e6) return("[large file; first line not read]")
  paste(readLines(path, n = 1, warn = FALSE), collapse = "")
}

roots <- c("data/raw", "data/processed", "results", "config")
roots <- roots[dir.exists(roots)]
all_files <- list.files(roots, recursive = TRUE, full.names = TRUE, all.files = FALSE)

patterns <- c(
  filtered_h5 = "filtered_feature_bc_matrix\\.h5$",
  matrix_mtx = "matrix\\.mtx\\.gz$",
  barcodes = "barcodes\\.tsv\\.gz$",
  features = "features\\.tsv\\.gz$",
  positions_csv = "tissue_positions(_list)?\\.csv$",
  scalefactors = "scalefactors_json\\.json$",
  lowres = "tissue_lowres_image\\.png$",
  hires = "tissue_hires_image\\.png$",
  h5ad = "\\.h5ad$",
  rds = "\\.rds$|\\.RDS$",
  loom = "\\.loom$",
  h5 = "\\.h5$",
  tsv_gz = "\\.tsv\\.gz$",
  tar = "\\.tar$",
  zip = "\\.zip$",
  png = "\\.png$",
  tif = "\\.tif$|\\.tiff$"
)

candidate_files <- rbindlist(lapply(names(patterns), function(tp) {
  hit <- all_files[grepl(patterns[[tp]], basename(all_files), ignore.case = TRUE) |
                     (tp %in% c("positions_csv", "scalefactors", "lowres", "hires") &
                        grepl(patterns[[tp]], all_files, ignore.case = TRUE))]
  if (!length(hit)) {
    return(data.table())
  }
  data.table(file_path = hit, file_type = tp)
}), fill = TRUE)

if (nrow(candidate_files) == 0) {
  candidate_files <- data.table(file_path = character(), file_type = character())
}

assign_dataset <- function(path) {
  if (grepl("GSE206306", path, ignore.case = TRUE)) return("GSE206306")
  if (grepl("GSE231630", path, ignore.case = TRUE)) return("GSE231630")
  if (grepl("HuBMAP|hubmap", path, ignore.case = TRUE)) return("HuBMAP_or_kidney_spatial_referenced")
  if (grepl("GSE231569", path, ignore.case = TRUE)) return("GSE231569_snRNA_not_spatial")
  if (grepl("spatial", path, ignore.case = TRUE)) return("other_spatial_reference")
  "other_nonspatial_project_file"
}

required_types <- c(
  "filtered_h5", "matrix_mtx", "barcodes", "features",
  "positions_csv", "scalefactors", "lowres", "hires", "h5ad", "rds"
)

inventory <- copy(candidate_files)
if (nrow(inventory) > 0) {
  inventory[, dataset_id := vapply(file_path, assign_dataset, character(1))]
  inventory[, file_name := basename(file_path)]
  inventory[, size_mb := vapply(file_path, file_size_mb, numeric(1))]
  inventory[, is_required := file_type %in% required_types]
  inventory[, detected := TRUE]
  inventory[, notes := fifelse(
    dataset_id %in% c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced"),
    "Candidate spatial-resource file category detected.",
    "Detected by extension pattern but not necessarily a spatial-resource file for Phase 28A."
  )]
  setcolorder(inventory, c("dataset_id", "file_path", "file_name", "file_type", "size_mb", "is_required", "detected", "notes"))
} else {
  inventory <- data.table(
    dataset_id = character(), file_path = character(), file_name = character(),
    file_type = character(), size_mb = numeric(), is_required = logical(),
    detected = logical(), notes = character()
  )
}

candidate_spatial_ids <- c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced")
required_visium_types <- c("filtered_h5_or_matrix_mtx", "barcodes", "features", "positions_csv", "scalefactors", "lowres_or_hires_image")
missing_rows <- rbindlist(lapply(candidate_spatial_ids, function(ds) {
  rbindlist(lapply(required_visium_types, function(tp) {
    has_type <- if (nrow(inventory) == 0) FALSE else any(inventory$dataset_id == ds & inventory$file_type == tp)
    broad_has <- if (nrow(inventory) == 0) FALSE else {
      if (tp == "filtered_h5_or_matrix_mtx") any(inventory$dataset_id == ds & inventory$file_type %in% c("filtered_h5", "matrix_mtx", "h5ad"))
      else if (tp == "lowres_or_hires_image") any(inventory$dataset_id == ds & inventory$file_type %in% c("lowres", "hires", "png", "tif"))
      else any(inventory$dataset_id == ds & inventory$file_type == tp)
    }
    if (has_type || broad_has) return(data.table())
    data.table(
      dataset_id = ds,
      file_path = "",
      file_name = "",
      file_type = tp,
      size_mb = NA_real_,
      is_required = TRUE,
      detected = FALSE,
      notes = "Required Visium-compatible file category not detected locally."
    )
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(missing_rows) > 0) {
  inventory <- rbind(inventory, missing_rows, fill = TRUE)
}
write_tsv(inventory, file.path(out_dir, "spatial_local_file_inventory_v0.1.tsv"))

dataset_info <- data.table(
  dataset_id = c("GSE206306", "GSE231630", "GSE231569_snRNA_not_spatial", "HuBMAP_or_kidney_spatial_referenced"),
  source = c("GEO", "GEO", "GEO/local project reference", "project text search"),
  species = c("Homo sapiens", "Homo sapiens", "Homo sapiens", "unknown"),
  tissue = c("human kidney papilla", "human kidney papilla / kidney medulla", "human kidney papilla nuclei", "kidney papilla/medulla candidate if present"),
  disease_context = c("calcium oxalate stone disease and healthy papilla context", "stone disease superseries / spatially anchored papilla atlas", "single-nucleus reference; not spatial", "unknown"),
  technology = c("spatial transcriptomics plus snRNA/CODEX according to GEO summary", "spatially anchored transcriptomic atlas; GEO supplement includes MTX/TAR/TIFF/TSV", "snRNA-seq", "Visium-compatible candidate if files exist"),
  local_path = c("data/raw/GSE206306; data/raw/spatial/GSE206306", "data/raw/GSE231630; data/raw/spatial/GSE231630", "data/raw/GSE231569; data/processed/gse231569_audited_seurat.rds", "not detected locally")
)

has_any <- function(ds, patt) {
  if (nrow(inventory) == 0) return(FALSE)
  any(inventory$dataset_id == ds & inventory$detected == TRUE & grepl(patt, inventory$file_type, ignore.case = TRUE))
}

has_path_file <- function(ds) {
  if (nrow(inventory) == 0) return(FALSE)
  any(inventory$dataset_id == ds & inventory$detected == TRUE)
}

manifest <- copy(dataset_info)
manifest[, `:=`(
  has_matrix = vapply(dataset_id, function(ds) has_any(ds, "filtered_h5|matrix_mtx|h5ad|rds"), logical(1)),
  has_barcodes = vapply(dataset_id, function(ds) has_any(ds, "barcodes"), logical(1)),
  has_features = vapply(dataset_id, function(ds) has_any(ds, "features"), logical(1)),
  has_positions = vapply(dataset_id, function(ds) has_any(ds, "positions_csv"), logical(1)),
  has_scalefactors = vapply(dataset_id, function(ds) has_any(ds, "scalefactors"), logical(1)),
  has_image = vapply(dataset_id, function(ds) has_any(ds, "lowres|hires|png|tif"), logical(1)),
  has_metadata = vapply(dataset_id, function(ds) has_path_file(ds) && ds %in% c("GSE206306", "GSE231630", "GSE231569_snRNA_not_spatial"), logical(1))
)]

manifest[dataset_id %in% c("GSE206306", "GSE231630") & !has_metadata, has_metadata := dir.exists(file.path("data/raw", dataset_id))]
manifest[dataset_id == "GSE231569_snRNA_not_spatial", `:=`(
  has_matrix = TRUE,
  has_barcodes = TRUE,
  has_features = TRUE,
  has_metadata = TRUE,
  has_positions = FALSE,
  has_scalefactors = FALSE,
  has_image = FALSE
)]

manifest[, download_status := fifelse(
  dataset_id == "GSE231569_snRNA_not_spatial", "local_snRNA_available_not_spatial",
  fifelse(has_matrix | has_positions | has_image, "partial_local_files_detected", "resource_limited_no_required_spatial_files_detected")
)]

manifest[, readiness_class := fifelse(
  dataset_id == "GSE231569_snRNA_not_spatial", "not_spatial",
  fifelse((has_matrix | (has_barcodes & has_features)) & has_positions & has_scalefactors & has_image & has_metadata, "analysis_ready_visium",
          fifelse((has_matrix | (has_barcodes & has_features)) & !(has_positions & has_image), "expression_only",
                  fifelse(has_image & !(has_matrix | (has_barcodes & has_features)), "image_only",
                          fifelse(has_metadata & !(has_matrix | has_positions | has_image), "metadata_only", "missing_required_files"))))
)]
manifest[, analysis_ready := readiness_class == "analysis_ready_visium"]
manifest[, claim_allowed := FALSE]
manifest[, notes := fifelse(
  analysis_ready,
  "Potentially analysis-ready; Seurat loading and spot-level QC are still required before any claim.",
  fifelse(readiness_class == "not_spatial",
          "Reference snRNA dataset; not a spatial transcriptomics resource.",
          "Spatial projection not allowed: matrix, coordinates, scalefactors and/or image are missing locally.")
)]

setcolorder(manifest, c(
  "dataset_id", "source", "species", "tissue", "disease_context", "technology",
  "has_matrix", "has_barcodes", "has_features", "has_positions", "has_scalefactors",
  "has_image", "has_metadata", "local_path", "download_status", "analysis_ready",
  "claim_allowed", "notes"
))
write_tsv(manifest, file.path(out_dir, "spatial_dataset_manifest_v0.1.tsv"))

ready <- manifest[analysis_ready == TRUE]
load_qc <- data.table(
  dataset_id = manifest$dataset_id,
  sample_id = NA_character_,
  n_spots = NA_integer_,
  n_genes_detected_median = NA_real_,
  umi_median = NA_real_,
  percent_mito_median = NA_real_,
  image_loaded = FALSE,
  coordinates_loaded = FALSE,
  qc_status = fifelse(manifest$analysis_ready, "not_attempted_by_manifest_script", "not_attempted_not_analysis_ready"),
  notes = fifelse(manifest$analysis_ready,
                  "Use scripts/11_spatial/01_load_spatial_dataset_qc.R for Seurat Load10X_Spatial QC.",
                  paste0("Readiness class: ", manifest$readiness_class, "; no Seurat loading attempted."))
)
write_tsv(load_qc, file.path(out_dir, "spatial_load_qc_v0.1.tsv"))

signature_files <- c(
  "results/supplementary_tables_v0.2/S5_audited_cell_type_markers.tsv",
  "results/tables/gse231569_marker_audit.tsv",
  "results/gene_sets/spatial_TAL_marker_set.txt",
  "results/gene_sets/magma_top50.txt",
  "results/gene_sets/magma_top100.txt",
  "results/gene_sets/magma_fdr05.txt"
)
signature_available <- any(file.exists(signature_files))
tal_genes <- if (file.exists("results/gene_sets/spatial_TAL_marker_set.txt")) unique(readLines("results/gene_sets/spatial_TAL_marker_set.txt", warn = FALSE)) else character()
risk_genes <- unique(unlist(lapply(c("results/gene_sets/magma_top50.txt", "results/gene_sets/magma_top100.txt", "results/gene_sets/magma_fdr05.txt"), function(p) {
  if (file.exists(p)) readLines(p, warn = FALSE) else character()
})))

projection <- data.table(
  dataset_id = manifest$dataset_id,
  sample_id = NA_character_,
  celltype_signature_available = signature_available,
  gene_overlap_n = NA_integer_,
  gene_overlap_fraction = NA_real_,
  loop_tal_signature_overlap = ifelse(length(tal_genes) > 0, NA_integer_, NA_integer_),
  risk_module_overlap = ifelse(length(risk_genes) > 0, NA_integer_, NA_integer_),
  projection_ready = FALSE,
  notes = fifelse(
    manifest$analysis_ready,
    "Spatial object must load successfully before gene overlap can be computed.",
    "Projection not reliable/not allowed because no analysis-ready spatial object is available."
  )
)
write_tsv(projection, file.path(out_dir, "spatial_projection_readiness_v0.1.tsv"))

memo <- c(
  "# Phase 28A spatial transcriptomics resource audit memo",
  "",
  "## Datasets checked",
  "",
  "- GSE206306: GEO metadata indicates human kidney papilla / stone disease spatial transcriptomics context; local required spatial files were not detected.",
  "- GSE231630: GEO SuperSeries for a spatially anchored human kidney papilla atlas; local required spatial files were not detected.",
  "- GSE231569_snRNA_not_spatial: local snRNA reference exists but is not a spatial transcriptomics dataset.",
  "- HuBMAP_or_kidney_spatial_referenced: no local Visium-compatible HuBMAP kidney papilla/medulla files were detected.",
  "",
  "## Required files present",
  "",
  "No candidate spatial dataset has the complete required combination of expression matrix, spatial coordinates, scalefactors/image and sample metadata. Existing local GSE231569 files are snRNA matrices and do not provide Visium coordinates or histology images.",
  "",
  "## Analysis readiness",
  "",
  "No dataset is analysis-ready for spatial projection. `analysis_ready` is FALSE for all candidate spatial datasets and `claim_allowed` is FALSE.",
  "",
  "## Seurat loading",
  "",
  "Seurat loading was not attempted because no dataset satisfied the analysis-ready Visium criteria. The loader script is provided for future use if complete Visium inputs are landed.",
  "",
  "## Spatial projection",
  "",
  "Spatial projection is not allowed in the current workspace. Existing GSE231569 cell-type signatures and MAGMA gene sets are available, but no spatial object is available for gene-overlap or spot-level QC.",
  "",
  "## Manuscript and figure claims",
  "",
  "No manuscript, main Figure 1-5, candidate gene tier or integrated evidence claim should change.",
  "",
  "## Manuscript-safe wording",
  "",
  "Spatial transcriptomic resources were audited for future projection analyses, but no complete analysis-ready spatial dataset was available locally; therefore spatial validation was not used as an evidence layer."
)
writeLines(memo, "docs/phase28a_spatial_resource_audit_memo_v0.1.md")

message("Phase 28A spatial resource audit complete: ", out_dir)
