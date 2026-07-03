#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

samples <- c("GSM6250307", "GSM6250308", "GSM6250309", "GSM6250310", "GSM7166170")
base_dir <- "data/processed/spatial/gse206306"
obj_dir <- file.path(base_dir, "seurat_objects")
out_dir <- "results/spatial/phase28b_loading_rescue"
fig_dir <- "results/figures/supplementary/spatial_qc"
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "")

has_file <- function(...) file.exists(file.path(...))

ensure_positions_csv <- function(sample_dir) {
  sp <- file.path(sample_dir, "spatial")
  pos_csv <- file.path(sp, "tissue_positions.csv")
  pos_list <- file.path(sp, "tissue_positions_list.csv")
  if (file.exists(pos_csv)) return("tissue_positions.csv")
  if (!file.exists(pos_list)) return("missing")
  lines <- readLines(pos_list, warn = FALSE)
  header <- "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"
  has_header <- length(lines) > 0 && grepl("barcode|in_tissue|array_row", lines[[1]], ignore.case = TRUE)
  if (has_header) {
    file.copy(pos_list, pos_csv, overwrite = FALSE)
  } else {
    writeLines(c(header, lines), pos_csv)
  }
  "created_from_tissue_positions_list.csv"
}

structure_rows <- rbindlist(lapply(samples, function(sid) {
  d <- file.path(base_dir, sid)
  pos_format <- ensure_positions_csv(d)
  data.table(
    sample_id = sid,
    has_h5 = has_file(d, "filtered_feature_bc_matrix.h5"),
    has_mtx = has_file(d, "filtered_feature_bc_matrix", "matrix.mtx.gz"),
    has_barcodes = has_file(d, "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
    has_features = has_file(d, "filtered_feature_bc_matrix", "features.tsv.gz"),
    has_positions = has_file(d, "spatial", "tissue_positions.csv") || has_file(d, "spatial", "tissue_positions_list.csv"),
    positions_format = pos_format,
    has_scalefactors = has_file(d, "spatial", "scalefactors_json.json"),
    has_lowres_image = has_file(d, "spatial", "tissue_lowres_image.png"),
    has_hires_image = has_file(d, "spatial", "tissue_hires_image.png"),
    structure_ready = (has_file(d, "filtered_feature_bc_matrix.h5") ||
                         (has_file(d, "filtered_feature_bc_matrix", "matrix.mtx.gz") &&
                            has_file(d, "filtered_feature_bc_matrix", "barcodes.tsv.gz") &&
                            has_file(d, "filtered_feature_bc_matrix", "features.tsv.gz"))) &&
      (has_file(d, "spatial", "tissue_positions.csv") || has_file(d, "spatial", "tissue_positions_list.csv")) &&
      has_file(d, "spatial", "scalefactors_json.json") &&
      (has_file(d, "spatial", "tissue_lowres_image.png") || has_file(d, "spatial", "tissue_hires_image.png")),
    notes = "Final structure check before Seurat loading."
  )
}), fill = TRUE)
write_tsv(structure_rows, file.path(out_dir, "spatial_structure_final_check_v0.1.tsv"))

read_positions_barcodes <- function(sample_dir) {
  pos <- file.path(sample_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(pos)) return(character())
  dt <- tryCatch(fread(pos), error = function(e) data.table())
  if (!nrow(dt)) return(character())
  if ("barcode" %in% names(dt)) return(dt$barcode)
  dt[[1]]
}

load_one <- function(sid) {
  sample_dir <- file.path(base_dir, sid)
  route <- "not_attempted"
  attempted <- FALSE
  success <- FALSE
  failure <- ""
  obj <- NULL

  hdf5r_available <- requireNamespace("hdf5r", quietly = TRUE)
  seurat_available <- requireNamespace("Seurat", quietly = TRUE)
  if (!seurat_available) {
    return(data.table(
      sample_id = sid, load_route = "none", load_attempted = FALSE, load_success = FALSE,
      n_spots_matrix = NA_integer_, n_spots_tissue = length(read_positions_barcodes(sample_dir)),
      barcode_overlap_n = NA_integer_, barcode_overlap_fraction = NA_real_,
      image_loaded = FALSE, coordinates_loaded = FALSE,
      n_features_median = NA_real_, n_counts_median = NA_real_, percent_mt_median = NA_real_,
      object_saved = FALSE, qc_status = "fail", failure_reason = "Seurat_missing",
      notes = "Seurat package is not installed."
    ))
  }

  if (hdf5r_available && file.exists(file.path(sample_dir, "filtered_feature_bc_matrix.h5"))) {
    attempted <- TRUE
    route <- "direct_h5_Load10X_Spatial"
    obj <- tryCatch(Seurat::Load10X_Spatial(data.dir = sample_dir), error = function(e) e)
    if (!inherits(obj, "error")) {
      success <- TRUE
    } else {
      failure <- conditionMessage(obj)
      obj <- NULL
    }
  }

  if (!success && file.exists(file.path(sample_dir, "filtered_feature_bc_matrix", "matrix.mtx.gz"))) {
    attempted <- TRUE
    route <- ifelse(nzchar(failure), "mtx_fallback_after_h5_failure", "mtx_fallback")
    obj <- tryCatch({
      counts <- Seurat::Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))
      if (is.list(counts)) {
        counts <- counts[[1]]
      }
      so <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial", project = sid)
      img <- tryCatch(Seurat::Read10X_Image(image.dir = file.path(sample_dir, "spatial"), filter.matrix = TRUE),
                      error = function(e) e)
      if (!inherits(img, "error")) {
        Seurat::DefaultAssay(img) <- "Spatial"
        so[["slice1"]] <- img
      } else {
        attr(so, "image_error") <- conditionMessage(img)
      }
      so
    }, error = function(e) e)
    if (!inherits(obj, "error")) {
      success <- TRUE
    } else {
      failure <- paste(c(failure, conditionMessage(obj)), collapse = " | ")
      obj <- NULL
    }
  }

  if (!success && file.exists(file.path(sample_dir, "filtered_feature_bc_matrix.h5"))) {
    attempted <- TRUE
    route <- ifelse(nzchar(failure), "manual_h5_matrix_image_after_Load10X_Spatial_failure", "manual_h5_matrix_image")
    obj <- tryCatch({
      counts <- Seurat::Read10X_h5(file.path(sample_dir, "filtered_feature_bc_matrix.h5"))
      if (is.list(counts)) counts <- counts[[1]]
      so <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial", project = sid)

      ## Seurat v5 Read10X_Image glob fails if both tissue_positions.csv and
      ## tissue_positions_list.csv are present. Use a temporary image dir with
      ## exactly one tissue_positions file without modifying the standardized tree.
      tmp_img <- file.path(tempdir(), paste0("spatial_image_", sid))
      if (dir.exists(tmp_img)) unlink(tmp_img, recursive = TRUE, force = TRUE)
      dir.create(tmp_img, recursive = TRUE, showWarnings = FALSE)
      for (fn in c("scalefactors_json.json", "tissue_lowres_image.png", "tissue_hires_image.png")) {
        src <- file.path(sample_dir, "spatial", fn)
        if (file.exists(src)) file.symlink(normalizePath(src), file.path(tmp_img, fn))
      }
      pos_src <- file.path(sample_dir, "spatial", "tissue_positions.csv")
      if (!file.exists(pos_src)) pos_src <- file.path(sample_dir, "spatial", "tissue_positions_list.csv")
      if (file.exists(pos_src)) file.symlink(normalizePath(pos_src), file.path(tmp_img, "tissue_positions.csv"))

      img <- tryCatch(Seurat::Read10X_Image(image.dir = tmp_img, filter.matrix = TRUE),
                      error = function(e) e)
      if (!inherits(img, "error")) {
        Seurat::DefaultAssay(img) <- "Spatial"
        so[["slice1"]] <- img
      } else {
        attr(so, "image_error") <- conditionMessage(img)
      }
      so
    }, error = function(e) e)
    if (!inherits(obj, "error")) {
      success <- TRUE
    } else {
      failure <- paste(c(failure, conditionMessage(obj)), collapse = " | ")
      obj <- NULL
    }
  }

  tissue_barcodes <- read_positions_barcodes(sample_dir)
  n_tissue <- length(unique(tissue_barcodes))
  if (success) {
    matrix_barcodes <- colnames(obj)
    overlap_n <- length(intersect(matrix_barcodes, tissue_barcodes))
    overlap_frac <- ifelse(length(matrix_barcodes) > 0, overlap_n / length(matrix_barcodes), NA_real_)
    image_loaded <- length(obj@images) > 0
    coordinates_loaded <- image_loaded
    nfeat <- obj$nFeature_Spatial
    ncount <- obj$nCount_Spatial
    pct_mt <- tryCatch(Seurat::PercentageFeatureSet(obj, pattern = "^MT-"), error = function(e) rep(NA_real_, ncol(obj)))
    object_path <- file.path(obj_dir, paste0(sid, "_spatial_qc.rds"))
    saveRDS(obj, object_path)
    qc_status <- if (overlap_frac >= 0.80 && coordinates_loaded && image_loaded &&
                       ncol(obj) > 100 && median(nfeat, na.rm = TRUE) > 200) {
      "pass"
    } else if (length(matrix_barcodes) > 0 && overlap_n > 0) {
      "partial_load"
    } else {
      "fail"
    }
    return(data.table(
      sample_id = sid, load_route = route, load_attempted = attempted, load_success = TRUE,
      n_spots_matrix = ncol(obj), n_spots_tissue = n_tissue,
      barcode_overlap_n = overlap_n, barcode_overlap_fraction = overlap_frac,
      image_loaded = image_loaded, coordinates_loaded = coordinates_loaded,
      n_features_median = median(nfeat, na.rm = TRUE),
      n_counts_median = median(ncount, na.rm = TRUE),
      percent_mt_median = median(pct_mt, na.rm = TRUE),
      object_saved = TRUE, qc_status = qc_status, failure_reason = "",
      notes = paste0("Saved: ", object_path)
    ))
  }

  data.table(
    sample_id = sid, load_route = route, load_attempted = attempted, load_success = FALSE,
    n_spots_matrix = NA_integer_, n_spots_tissue = n_tissue,
    barcode_overlap_n = NA_integer_, barcode_overlap_fraction = NA_real_,
    image_loaded = FALSE, coordinates_loaded = FALSE,
    n_features_median = NA_real_, n_counts_median = NA_real_, percent_mt_median = NA_real_,
    object_saved = FALSE, qc_status = "fail", failure_reason = failure,
    notes = ifelse(requireNamespace("hdf5r", quietly = TRUE), "Both direct and fallback routes failed.", "hdf5r unavailable; fallback route attempted if mtx exists.")
  )
}

qc <- rbindlist(lapply(samples, load_one), fill = TRUE)
write_tsv(qc, file.path(out_dir, "spatial_load_qc_v0.2.tsv"))

gate <- qc[, .(
  sample_id,
  object_loaded = load_success,
  qc_status,
  projection_allowed = qc_status == "pass",
  reason = fifelse(qc_status == "pass", "object_loaded_and_qc_pass",
                   fifelse(qc_status == "partial_load", "partial_load_not_projection_ready", failure_reason)),
  next_phase = fifelse(qc_status == "pass", "Phase 28C spatial signature projection",
                       fifelse(qc_status == "partial_load", "image/coordinate rescue, not projection", "freeze spatial as loading-limited or fix loader dependency"))
)]
write_tsv(gate, file.path(out_dir, "spatial_projection_gate_v0.1.tsv"))

## QC contact sheet only for loaded samples. This is loading QC only.
loaded <- qc[load_success == TRUE]
pdf_path <- file.path(fig_dir, "gse206306_spatial_loading_qc_contact_sheet_v0.1.pdf")
png_path <- file.path(fig_dir, "gse206306_spatial_loading_qc_contact_sheet_v0.1.png")
if (nrow(loaded) > 0 && requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
  plots <- lapply(loaded$sample_id, function(sid) {
    obj <- readRDS(file.path(obj_dir, paste0(sid, "_spatial_qc.rds")))
    row <- loaded[sample_id == sid][1]
    p <- tryCatch(
      Seurat::SpatialDimPlot(obj, images = "slice1", pt.size.factor = 1.2) +
        ggplot2::ggtitle(paste0(sid, " | spots=", row$n_spots_matrix,
                                " | med genes=", round(row$n_features_median, 0))) +
        ggplot2::labs(caption = "Loading QC only; no biological projection"),
      error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0, y = 0, label = paste(sid, conditionMessage(e), sep = "\n")) +
          ggplot2::theme_void()
      }
    )
    p + ggplot2::theme(plot.title = ggplot2::element_text(size = 9), plot.caption = ggplot2::element_text(size = 7))
  })
  sheet <- patchwork::wrap_plots(plots, ncol = 3) +
    patchwork::plot_annotation(title = "GSE206306 spatial loading QC contact sheet",
                               subtitle = "Spot/image loading QC only; no biological projection or spatial validation claim")
  ggplot2::ggsave(pdf_path, sheet, width = 11, height = 8)
  ggplot2::ggsave(png_path, sheet, width = 11, height = 8, dpi = 300)
} else {
  pdf(pdf_path, width = 11, height = 8)
  plot.new()
  text(0.5, 0.6, "GSE206306 spatial loading QC", cex = 1.2)
  text(0.5, 0.45, "No samples loaded; no biological projection performed.", cex = 0.9)
  dev.off()
  png(png_path, width = 3300, height = 2400, res = 300)
  plot.new()
  text(0.5, 0.6, "GSE206306 spatial loading QC", cex = 1.2)
  text(0.5, 0.45, "No samples loaded; no biological projection performed.", cex = 0.9)
  dev.off()
}

message("wrote spatial_load_qc_v0.2 and projection gate")
