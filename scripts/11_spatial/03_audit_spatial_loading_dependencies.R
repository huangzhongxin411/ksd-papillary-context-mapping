#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

dir.create("results/spatial/phase28b_loading_rescue", recursive = TRUE, showWarnings = FALSE)

pkgs <- data.table(
  package = c("Seurat", "SeuratObject", "Matrix", "hdf5r", "jsonlite", "png", "ggplot2", "patchwork", "dplyr", "readr"),
  required_for = c(
    "direct Load10X_Spatial and fallback object construction",
    "Seurat object infrastructure",
    "sparse matrix operations",
    "direct 10x HDF5 matrix loading",
    "scalefactors_json parsing",
    "spatial image reading",
    "QC plotting",
    "QC contact sheet composition",
    "metadata/QC table manipulation",
    "robust delimited file reading"
  )
)

pkgs[, installed := vapply(package, requireNamespace, logical(1), quietly = TRUE)]
pkgs[, version := vapply(package, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(p))
}, character(1))]
pkgs[, status := fifelse(installed, "PASS", fifelse(package == "hdf5r", "WARNING", "FAIL"))]
pkgs[, notes := fifelse(
  installed,
  "Package available.",
  fifelse(package == "hdf5r",
          "Missing; direct HDF5 route unavailable, use Python h5-to-mtx fallback if install fails.",
          "Missing package may limit loading/QC plotting.")
)]

fwrite(pkgs, "results/spatial/phase28b_loading_rescue/spatial_loading_dependency_audit_v0.1.tsv", sep = "\t", na = "")
message("wrote spatial loading dependency audit")
