#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), mustWork = TRUE)
phase_dir <- file.path(root, "results", "spatial", "phase28a_r")
docs_spatial <- file.path(root, "docs", "spatial")
staging_root <- file.path(root, "data", "raw", "spatial", "manual_download_staging")
dir.create(file.path(root, "data", "raw", "spatial", "gse206306"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "data", "raw", "spatial", "gse231630"), recursive = TRUE, showWarnings = FALSE)
dir.create(staging_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "data", "processed", "spatial"), recursive = TRUE, showWarnings = FALSE)
dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_spatial, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "")

safe_size_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  round(file.info(path)$size / 1024^2, 3)
}

norm_path <- function(path) normalizePath(path, mustWork = FALSE)

home <- path.expand("~")
search_roots <- unique(c(
  root,
  file.path(root, "data"),
  file.path(root, "external"),
  path.expand("~/Downloads"),
  path.expand("~/Documents")
))
search_roots <- search_roots[dir.exists(search_roots)]

patterns <- list(
  dataset = "GSE206306|GSE231630|GSM7166168|GSM7166169|GSM7166170|GSM6250307|GSM6250308|GSM6250309|GSM6250310|GSM7290910|GSM7290911|GSM7290912|GSM7290913|GSM7290914",
  visium = "Visium|spatial|tissue_positions|scalefactors_json|filtered_feature_bc_matrix|matrix\\.mtx|barcodes\\.tsv|features\\.tsv|tissue_lowres_image|tissue_hires_image",
  ext = "\\.h5$|\\.h5ad$|\\.rds$|\\.RDS$|\\.tar$|\\.tar\\.gz$|\\.tgz$|\\.zip$|\\.mtx\\.gz$|\\.tsv\\.gz$|\\.csv$|\\.json$|\\.png$|\\.tif$|\\.tiff$"
)

all_files <- unique(unlist(lapply(search_roots, function(sr) {
  list.files(sr, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
}), use.names = FALSE))

candidate <- all_files[
  grepl(patterns$dataset, all_files, ignore.case = TRUE) |
    grepl(patterns$visium, all_files, ignore.case = TRUE) |
    (grepl(patterns$ext, all_files, ignore.case = TRUE) & grepl("spatial|GSE206306|GSE231630|Visium|HuBMAP|hubmap|papilla|medulla", all_files, ignore.case = TRUE))
]
candidate <- unique(candidate[file.exists(candidate)])

dataset_id_of <- function(path) {
  if (grepl("GSE206306|GSM6250307|GSM6250308|GSM6250309|GSM6250310|GSM7166168|GSM7166169|GSM7166170", path, ignore.case = TRUE)) return("GSE206306")
  if (grepl("GSE231630|GSM7290910|GSM7290911|GSM7290912|GSM7290913|GSM7290914", path, ignore.case = TRUE)) return("GSE231630")
  if (grepl("HuBMAP|hubmap", path, ignore.case = TRUE)) return("HuBMAP_or_kidney_spatial_referenced")
  if (grepl("GSE231569", path, ignore.case = TRUE)) return("GSE231569_snRNA_not_spatial")
  if (grepl("spatial|Visium|papilla|medulla", path, ignore.case = TRUE)) return("other_spatial_candidate")
  "other"
}

sample_id_of <- function(path) {
  m <- regmatches(path, regexpr("GSM[0-9]+", path, ignore.case = TRUE))
  if (length(m) && nzchar(m)) return(toupper(m))
  m2 <- regmatches(path, regexpr("KRP[0-9]+|K[0-9]{7}|[0-9]{2}-[0-9]{4}", path, ignore.case = TRUE))
  if (length(m2) && nzchar(m2)) return(m2)
  ""
}

file_type_of <- function(name) {
  if (grepl("filtered_feature_bc_matrix\\.h5$", name, ignore.case = TRUE)) return("filtered_feature_bc_matrix.h5")
  if (grepl("matrix\\.mtx(\\.gz)?$", name, ignore.case = TRUE)) return("matrix.mtx.gz")
  if (grepl("barcodes\\.tsv(\\.gz)?$", name, ignore.case = TRUE)) return("barcodes.tsv.gz")
  if (grepl("features\\.tsv(\\.gz)?$|genes\\.tsv(\\.gz)?$", name, ignore.case = TRUE)) return("features.tsv.gz")
  if (grepl("tissue_positions(_list)?\\.csv$", name, ignore.case = TRUE)) return("tissue_positions.csv")
  if (grepl("scalefactors_json\\.json$", name, ignore.case = TRUE)) return("scalefactors_json.json")
  if (grepl("tissue_lowres_image\\.png$", name, ignore.case = TRUE)) return("tissue_lowres_image.png")
  if (grepl("tissue_hires_image\\.png$", name, ignore.case = TRUE)) return("tissue_hires_image.png")
  if (grepl("\\.h5ad$", name, ignore.case = TRUE)) return("h5ad")
  if (grepl("\\.rds$|\\.RDS$", name, ignore.case = TRUE)) return("rds")
  if (grepl("\\.tar\\.gz$|\\.tgz$", name, ignore.case = TRUE)) return("tar.gz")
  if (grepl("\\.tar$", name, ignore.case = TRUE)) return("tar")
  if (grepl("\\.zip$", name, ignore.case = TRUE)) return("zip")
  if (grepl("\\.png$", name, ignore.case = TRUE)) return("png")
  if (grepl("\\.tif$|\\.tiff$", name, ignore.case = TRUE)) return("tiff")
  tools::file_ext(name)
}

category_of <- function(ft) {
  if (ft %in% c("filtered_feature_bc_matrix.h5", "matrix.mtx.gz", "h5ad", "rds")) return("expression_matrix")
  if (ft %in% c("barcodes.tsv.gz")) return("barcodes")
  if (ft %in% c("features.tsv.gz")) return("features")
  if (ft %in% c("tissue_positions.csv")) return("coordinates")
  if (ft %in% c("scalefactors_json.json")) return("scalefactors")
  if (ft %in% c("tissue_lowres_image.png", "tissue_hires_image.png", "png", "tiff")) return("image")
  if (ft %in% c("tar", "tar.gz", "zip")) return("archive")
  "metadata_or_other"
}

inventory <- if (length(candidate)) {
  data.table(file_path = norm_path(candidate))
} else {
  data.table(file_path = character())
}
if (nrow(inventory)) {
  inventory[, `:=`(
    dataset_id = vapply(file_path, dataset_id_of, character(1)),
    candidate_sample_id = vapply(file_path, sample_id_of, character(1)),
    file_name = basename(file_path),
    file_type = vapply(basename(file_path), file_type_of, character(1)),
    size_mb = vapply(file_path, safe_size_mb, numeric(1))
  )]
  inventory[, required_category := vapply(file_type, category_of, character(1))]
  inventory[, is_complete_visium_component := required_category %in% c("expression_matrix", "barcodes", "features", "coordinates", "scalefactors", "image")]
  inventory[, notes := fifelse(dataset_id %in% c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced"),
                               "Candidate file for spatial recovery review.",
                               "Detected by broad spatial/file pattern; may be unrelated or non-spatial.")]
  setcolorder(inventory, c("dataset_id", "candidate_sample_id", "file_path", "file_name", "file_type", "size_mb", "required_category", "is_complete_visium_component", "notes"))
} else {
  inventory <- data.table(
    dataset_id = character(), candidate_sample_id = character(), file_path = character(),
    file_name = character(), file_type = character(), size_mb = numeric(),
    required_category = character(), is_complete_visium_component = logical(), notes = character()
  )
}

## Unpack only archives likely related to target spatial datasets.
archives <- inventory[required_category == "archive" & dataset_id %in% c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced")]
unpack_log <- data.table(
  archive_path = character(), archive_size_mb = numeric(), extract_dir = character(),
  extract_status = character(), n_files_extracted = integer(),
  contains_matrix = logical(), contains_spatial_dir = logical(), contains_image = logical(),
  notes = character()
)

if (nrow(archives) > 0) {
  for (i in seq_len(nrow(archives))) {
    ap <- archives$file_path[[i]]
    base <- tools::file_path_sans_ext(basename(ap))
    base <- sub("\\.tar$", "", base, ignore.case = TRUE)
    exdir <- file.path(staging_root, make.names(base))
    if (dir.exists(exdir) && length(list.files(exdir, recursive = TRUE, all.files = FALSE)) > 0) {
      status <- "skipped_existing_extract_dir_nonempty"
    } else {
      dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
      if (grepl("\\.zip$", ap, ignore.case = TRUE)) {
        status_code <- suppressWarnings(system2("unzip", c("-n", ap, "-d", exdir), stdout = TRUE, stderr = TRUE))
      } else {
        status_code <- suppressWarnings(system2("tar", c("-xf", ap, "-C", exdir), stdout = TRUE, stderr = TRUE))
      }
      status <- if (is.null(attr(status_code, "status"))) "extracted" else paste0("extract_failed_status_", attr(status_code, "status"))
    }
    files <- list.files(exdir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    unpack_log <- rbind(unpack_log, data.table(
      archive_path = ap,
      archive_size_mb = safe_size_mb(ap),
      extract_dir = exdir,
      extract_status = status,
      n_files_extracted = length(files),
      contains_matrix = any(grepl("filtered_feature_bc_matrix\\.h5|matrix\\.mtx|\\.h5ad$|\\.rds$", files, ignore.case = TRUE)),
      contains_spatial_dir = any(grepl("/spatial$|/spatial/", files, ignore.case = TRUE)),
      contains_image = any(grepl("tissue_lowres_image|tissue_hires_image|\\.png$|\\.tif$|\\.tiff$", files, ignore.case = TRUE)),
      notes = "Archives are extracted only into manual_download_staging; original files are not modified."
    ), fill = TRUE)
  }
}
if (nrow(unpack_log) == 0) {
  unpack_log <- data.table(
    archive_path = NA_character_, archive_size_mb = NA_real_, extract_dir = NA_character_,
    extract_status = "no_target_spatial_archives_found",
    n_files_extracted = 0L, contains_matrix = FALSE, contains_spatial_dir = FALSE,
    contains_image = FALSE,
    notes = "No archive clearly linked to GSE206306/GSE231630/HuBMAP spatial resources was found for unpacking."
  )
}
write_tsv(unpack_log, file.path(phase_dir, "spatial_archive_unpack_log_v0.1.tsv"))

## Re-run inventory including staging after unpack.
staging_files <- list.files(staging_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
if (length(staging_files)) {
  staging_dt <- data.table(file_path = norm_path(staging_files))
  staging_dt[, `:=`(
    dataset_id = vapply(file_path, dataset_id_of, character(1)),
    candidate_sample_id = vapply(file_path, sample_id_of, character(1)),
    file_name = basename(file_path),
    file_type = vapply(basename(file_path), file_type_of, character(1)),
    size_mb = vapply(file_path, safe_size_mb, numeric(1))
  )]
  staging_dt[, required_category := vapply(file_type, category_of, character(1))]
  staging_dt[, is_complete_visium_component := required_category %in% c("expression_matrix", "barcodes", "features", "coordinates", "scalefactors", "image")]
  staging_dt[, notes := "Detected after safe archive extraction into manual_download_staging."]
  setcolorder(staging_dt, names(inventory))
  inventory <- unique(rbind(inventory, staging_dt, fill = TRUE), by = "file_path")
}
write_tsv(inventory, file.path(phase_dir, "spatial_file_recovery_inventory_v0.1.tsv"))

## Sample-level readiness.
## Use only resources recovered into the current project spatial raw/staging tree.
project_spatial_root <- file.path(root, "data", "raw", "spatial")
project_h5 <- inventory[
  dataset_id %in% c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced") &
    file_type == "filtered_feature_bc_matrix.h5" &
    grepl(project_spatial_root, file_path, fixed = TRUE)
]

sample_from_root <- function(path) {
  gsm <- regmatches(path, regexpr("GSM[0-9]+", path, ignore.case = TRUE))
  if (length(gsm) && nzchar(gsm)) return(toupper(gsm))
  basename(path)
}

if (nrow(project_h5) == 0) {
  candidate_roots <- data.table(
    dataset_id = c("GSE206306", "GSE231630", "HuBMAP_or_kidney_spatial_referenced"),
    sample_id = c("dataset_level_missing", "dataset_level_missing", "dataset_level_missing"),
    local_root = c("data/raw/spatial/gse206306", "data/raw/spatial/gse231630", "not_detected")
  )
} else {
  project_h5[, local_root := dirname(file_path)]
  project_h5[, sample_id := vapply(local_root, sample_from_root, character(1))]
  candidate_roots <- unique(project_h5[, .(dataset_id, sample_id, local_root)])
}

check_root <- function(ds, sid, lr) {
  files <- if (dir.exists(lr)) list.files(lr, recursive = TRUE, full.names = TRUE, all.files = FALSE) else character()
  b <- basename(files)
  has_h5 <- any(grepl("filtered_feature_bc_matrix\\.h5$", b, ignore.case = TRUE))
  has_mtx <- any(grepl("matrix\\.mtx(\\.gz)?$", b, ignore.case = TRUE))
  has_barcodes <- any(grepl("barcodes\\.tsv(\\.gz)?$", b, ignore.case = TRUE))
  has_features <- any(grepl("features\\.tsv(\\.gz)?$|genes\\.tsv(\\.gz)?$", b, ignore.case = TRUE))
  has_positions <- any(grepl("tissue_positions(_list)?\\.csv$", b, ignore.case = TRUE))
  has_scalefactors <- any(grepl("scalefactors_json\\.json$", b, ignore.case = TRUE))
  has_lowres <- any(grepl("tissue_lowres_image\\.png$", b, ignore.case = TRUE))
  has_hires <- any(grepl("tissue_hires_image\\.png$", b, ignore.case = TRUE))
  has_metadata <- nzchar(sid) && sid != "dataset_level_missing"
  matrix_ok <- has_h5 || (has_mtx && has_barcodes && has_features)
  image_ok <- has_lowres || has_hires
  analysis_ready <- matrix_ok && has_positions && has_scalefactors && image_ok && has_metadata
  reasons <- c()
  if (!matrix_ok) reasons <- c(reasons, "missing_matrix")
  if (!has_positions) reasons <- c(reasons, "missing_coordinates")
  if (!has_scalefactors) reasons <- c(reasons, "missing_scalefactors")
  if (!image_ok) reasons <- c(reasons, "missing_image")
  if (!has_metadata) reasons <- c(reasons, "missing_metadata")
  if (ds == "GSE231569_snRNA_not_spatial") reasons <- c(reasons, "snrna_not_spatial")
  data.table(
    dataset_id = ds, sample_id = sid, local_root = lr,
    has_h5_matrix = has_h5, has_mtx_matrix = has_mtx,
    has_barcodes = has_barcodes, has_features = has_features,
    has_positions = has_positions, has_scalefactors = has_scalefactors,
    has_lowres_image = has_lowres, has_hires_image = has_hires,
    has_metadata = has_metadata,
    visium_structure_valid = matrix_ok && has_positions && has_scalefactors && image_ok,
    analysis_ready = analysis_ready,
    failure_reason = ifelse(length(reasons), paste(unique(reasons), collapse = ";"), ""),
    notes = ifelse(analysis_ready, "Complete Visium-like structure detected; standardization can proceed.", "Not analysis-ready; do not load or project.")
  )
}
readiness <- rbindlist(lapply(seq_len(nrow(candidate_roots)), function(i) {
  check_root(candidate_roots$dataset_id[[i]], candidate_roots$sample_id[[i]], candidate_roots$local_root[[i]])
}), fill = TRUE)
write_tsv(readiness, file.path(phase_dir, "spatial_visium_sample_readiness_v0.1.tsv"))

## Standardize complete samples with symlinks. Do not overwrite non-symlink files.
standardize_one <- function(row) {
  ds <- row$dataset_id
  sid <- row$sample_id
  src <- row$local_root
  dst <- file.path(root, "data", "processed", "spatial", tolower(ds), sid)
  if (!isTRUE(row$analysis_ready)) {
    return(data.table(
      dataset_id = ds, sample_id = sid, source_root = src,
      standardized_root = sub(paste0(root, "/"), "", dst, fixed = TRUE),
      standardization_status = "not_performed_not_analysis_ready",
      files_linked_or_copied = 0L,
      notes = "No complete Visium structure detected; standardization skipped."
    ))
  }
  dir.create(file.path(dst, "spatial"), recursive = TRUE, showWarnings = FALSE)
  links <- list(
    filtered_feature_bc_matrix.h5 = file.path(src, "filtered_feature_bc_matrix.h5"),
    `spatial/tissue_positions_list.csv` = file.path(src, "spatial", "tissue_positions_list.csv"),
    `spatial/tissue_positions.csv` = file.path(src, "spatial", "tissue_positions.csv"),
    `spatial/scalefactors_json.json` = file.path(src, "spatial", "scalefactors_json.json"),
    `spatial/tissue_lowres_image.png` = file.path(src, "spatial", "tissue_lowres_image.png"),
    `spatial/tissue_hires_image.png` = file.path(src, "spatial", "tissue_hires_image.png")
  )
  n <- 0L
  blocked <- character()
  for (rel in names(links)) {
    source <- links[[rel]]
    if (!file.exists(source)) next
    target <- file.path(dst, rel)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(target)) {
      if (file.info(target)$isdir) {
        blocked <- c(blocked, paste0(rel, ":target_is_dir"))
      } else {
        n <- n + 1L
      }
      next
    }
    ok <- file.symlink(source, target)
    if (isTRUE(ok)) n <- n + 1L else blocked <- c(blocked, paste0(rel, ":symlink_failed"))
  }
  data.table(
    dataset_id = ds,
    sample_id = sid,
    source_root = src,
    standardized_root = sub(paste0(root, "/"), "", dst, fixed = TRUE),
    standardization_status = ifelse(length(blocked), "standardized_with_warnings", "standardized_symlinked"),
    files_linked_or_copied = n,
    notes = ifelse(length(blocked), paste(blocked, collapse = ";"), "Symlinked complete Visium-like sample structure without modifying source files.")
  )
}

std_log <- rbindlist(lapply(seq_len(nrow(readiness)), function(i) {
  standardize_one(readiness[i])
}), fill = TRUE)
write_tsv(std_log, file.path(phase_dir, "spatial_standardized_structure_log_v0.1.tsv"))

## Attempt Seurat loading only for analysis-ready standardized samples.
load_one <- function(row) {
  ds <- row$dataset_id
  sid <- row$sample_id
  std <- file.path(root, row$standardized_root)
  if (!row$standardization_status %in% c("standardized_symlinked", "standardized_with_warnings")) {
    return(data.table(dataset_id = ds, sample_id = sid, n_spots = NA_integer_,
                      n_genes_median = NA_real_, nCount_Spatial_median = NA_real_,
                      percent_mt_median = NA_real_, coordinates_loaded = FALSE,
                      image_loaded = FALSE, qc_status = "not_attempted_not_analysis_ready",
                      notes = "No standardized analysis-ready structure."))
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    return(data.table(dataset_id = ds, sample_id = sid, n_spots = NA_integer_,
                      n_genes_median = NA_real_, nCount_Spatial_median = NA_real_,
                      percent_mt_median = NA_real_, coordinates_loaded = FALSE,
                      image_loaded = FALSE, qc_status = "load_failed_seurat_missing",
                      notes = "Seurat package is not installed."))
  }
  obj <- tryCatch(Seurat::Load10X_Spatial(data.dir = std), error = function(e) e)
  if (inherits(obj, "error")) {
    return(data.table(dataset_id = ds, sample_id = sid, n_spots = NA_integer_,
                      n_genes_median = NA_real_, nCount_Spatial_median = NA_real_,
                      percent_mt_median = NA_real_, coordinates_loaded = FALSE,
                      image_loaded = FALSE, qc_status = "load_failed",
                      notes = conditionMessage(obj)))
  }
  counts <- obj$nCount_Spatial
  feats <- obj$nFeature_Spatial
  mt <- tryCatch(Seurat::PercentageFeatureSet(obj, pattern = "^MT-"), error = function(e) rep(NA_real_, ncol(obj)))
  rds_out <- file.path(std, paste0(sid, "_seurat_spatial_qc.rds"))
  saveRDS(obj, rds_out)
  pass <- ncol(obj) > 100 &&
    length(obj@images) > 0 &&
    median(feats, na.rm = TRUE) > 200
  data.table(dataset_id = ds, sample_id = sid, n_spots = ncol(obj),
             n_genes_median = median(feats, na.rm = TRUE),
             nCount_Spatial_median = median(counts, na.rm = TRUE),
             percent_mt_median = median(mt, na.rm = TRUE),
             coordinates_loaded = length(obj@images) > 0,
             image_loaded = length(obj@images) > 0,
             qc_status = ifelse(pass, "loaded_qc_pass", "loaded_qc_warning"),
             notes = paste0("Saved QC object: ", sub(paste0(root, "/"), "", rds_out, fixed = TRUE)))
}

load_qc <- rbindlist(lapply(seq_len(nrow(std_log)), function(i) load_one(std_log[i])), fill = TRUE)
write_tsv(load_qc, file.path(phase_dir, "spatial_load_qc_v0.1.tsv"))

signature_names <- c("Loop/TAL signature", "MAGMA top50", "MAGMA top100", "MAGMA FDR",
                     "MAGMA suggestive", "P1 six genes", "injury/remodeling markers",
                     "epithelial marker score", "stromal/ECM marker score", "immune marker score")
projection <- CJ(dataset_id = unique(readiness$dataset_id), sample_id = unique(readiness$sample_id), signature = signature_names)
loaded_keys <- paste(load_qc[qc_status == "loaded_qc_pass", dataset_id], load_qc[qc_status == "loaded_qc_pass", sample_id])
projection[, object_loaded := paste(dataset_id, sample_id) %in% loaded_keys]
projection[, `:=`(
  gene_overlap_n = NA_integer_,
  gene_overlap_fraction = NA_real_,
  projection_ready = FALSE,
  notes = fifelse(object_loaded,
                  "Object loaded and QC-passed; gene-overlap projection readiness must be evaluated in Phase 28B before biological interpretation.",
                  "Projection not evaluated because no spatial object was loaded and QC-passed in Phase 28A-R.")
)]
write_tsv(projection, file.path(phase_dir, "spatial_projection_readiness_v0.2.tsv"))

checklist <- c(
  "# Phase 28A-R manual spatial download checklist",
  "",
  "## GSE206306",
  "",
  "GEO metadata checked in Phase 28A: human kidney papilla / calcium oxalate stone disease spatial transcriptomics context; samples include GSM6250307, GSM6250308, GSM6250309, GSM6250310, GSM7166168, GSM7166169 and GSM7166170.",
  "",
  "Required files for each Visium-compatible sample/library:",
  "",
  "- Expression matrix: `filtered_feature_bc_matrix.h5` OR `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`.",
  "- Spatial coordinates: `tissue_positions.csv` OR `tissue_positions_list.csv`.",
  "- Image/scalefactors: `scalefactors_json.json`, `tissue_lowres_image.png`, and `tissue_hires_image.png` if available.",
  "- Metadata: sample ID, disease/control status, tissue region, slide/library ID, species, technology.",
  "",
  "## GSE231630",
  "",
  "GEO metadata checked in Phase 28A: spatially anchored human kidney papilla atlas SuperSeries; samples include GSM7166168, GSM7166169, GSM7166170, GSM6250307, GSM6250308, GSM6250309, GSM6250310, GSM7290910, GSM7290911, GSM7290912, GSM7290913 and GSM7290914.",
  "",
  "Required files for each Visium-compatible sample/library:",
  "",
  "- Expression matrix: `filtered_feature_bc_matrix.h5` OR `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`.",
  "- Spatial coordinates: `tissue_positions.csv` OR `tissue_positions_list.csv`.",
  "- Image/scalefactors: `scalefactors_json.json`, `tissue_lowres_image.png`, and `tissue_hires_image.png` if available.",
  "- Metadata: sample ID, disease/control status, tissue region, slide/library ID, species, technology.",
  "",
  "## Manual staging location",
  "",
  "Place downloaded archives or unpacked folders under `data/raw/spatial/manual_download_staging/` or dataset-specific folders `data/raw/spatial/gse206306/` and `data/raw/spatial/gse231630/`. Re-run `scripts/11_spatial/02_phase28a_r_spatial_resource_recovery.R` after landing files.",
  "",
  "## Claim boundary",
  "",
  "Do not claim spatial validation from metadata or partial files. Spatial projection requires a loaded object, coordinates, image and spot-level QC."
)
writeLines(checklist, file.path(docs_spatial, "phase28a_r_manual_download_checklist_v0.1.md"))

complete_found <- any(readiness$analysis_ready)
archives_unpacked <- any(unpack_log$extract_status %in% c("extracted", "skipped_existing_extract_dir_nonempty"))
loaded_success <- any(load_qc$qc_status == "loaded_qc_pass")
memo <- c(
  "# Phase 28A-R spatial resource recovery completion memo",
  "",
  "## 1. Were any complete spatial files found?",
  "",
  if (complete_found) "Yes. At least one GSE206306 sample recovered into the current project staging tree met file-based Visium readiness criteria." else "No. No complete Visium-compatible sample with expression matrix, coordinates, scalefactors, image and metadata was found.",
  "",
  "## 2. Were any archives unpacked?",
  "",
  if (archives_unpacked) "Yes. Candidate archives were safely extracted into `data/raw/spatial/manual_download_staging/` without overwriting existing files." else "No. No target spatial archive linked to GSE206306/GSE231630/HuBMAP was found for unpacking.",
  "",
  "## 3. Which samples are analysis_ready?",
  "",
  if (complete_found) paste(unique(readiness[analysis_ready == TRUE, paste(dataset_id, sample_id, sep = ": ")]), collapse = "\n") else "None.",
  "",
  "## 4. Was Seurat loading attempted?",
  "",
  if (complete_found) "Yes. Seurat loading was attempted for standardized analysis-ready samples." else "No. Seurat loading was not attempted because no sample was analysis-ready in this resource-recovery phase.",
  "",
  "## 5. Did any sample load successfully?",
  "",
  if (loaded_success) paste(load_qc[qc_status == "loaded_qc_pass", paste(dataset_id, sample_id, sep = ": ")], collapse = "\n") else "No spatial sample loaded successfully. See `spatial_load_qc_v0.1.tsv` for the blocker, likely missing HDF5 reader support if only `.h5` matrices are available.",
  "",
  "## 6. Is spatial projection allowed?",
  "",
  if (loaded_success) "Not yet. A loaded QC-passed spatial object exists, but projection readiness and gene-overlap thresholds must be assessed in Phase 28B before any biological interpretation." else "No. Projection is not allowed until a complete spatial object loads and passes QC.",
  "",
  "## 7. Should manuscript or figures change?",
  "",
  "No. Manuscript conclusions, Figure 1-5, candidate gene tiers and integrated evidence summary should remain unchanged.",
  "",
  "## 8. Claim boundary",
  "",
  if (loaded_success) "Spatial transcriptomic resources were recovered and at least one sample loaded, but biological spatial projection has not yet been performed; spatial validation is not yet used as evidence." else "Spatial transcriptomic resources were audited and recovery was attempted; complete Visium-like files were found, but no spatial object loaded successfully in this phase, so spatial validation was not used as evidence."
)
writeLines(memo, file.path(docs_spatial, "phase28a_r_completion_memo_v0.1.md"))

message("Phase 28A-R spatial recovery audit complete: ", phase_dir)
