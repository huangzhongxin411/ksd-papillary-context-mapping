suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
})

fig_dir <- "results/figures/supplementary/spatial_projection"
out_dir <- "results/spatial/phase28c_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal <- list(
  deep_teal = "#005A64",
  slate = "#6F929B",
  sand = "#C99A2E",
  terracotta = "#A85B4B",
  pale_bg = "#EAF1F2",
  light_grey = "#D3DADC",
  ink = "#22313B"
)

load_qc <- read_tsv("results/spatial/phase28b_loading_rescue/spatial_load_qc_v0.2.tsv", show_col_types = FALSE)
objects <- list.files(obj_dir, pattern = "_spatial_scores_v0\\.1\\.rds$", full.names = TRUE)
if (!length(objects)) stop("No scored objects found.")

plot_missing <- function(sample_id, signature) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.55, label = sample_id, fontface = "bold", size = 4, color = pal$ink) +
    annotate("text", x = 0.5, y = 0.42, label = paste(signature, "not scored"), size = 3, color = pal$slate) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.background = element_rect(fill = "white", color = pal$light_grey, linewidth = 0.4))
}

plot_one <- function(path, signature, cols) {
  sample_id <- sub("_spatial_scores_v0\\.1\\.rds$", "", basename(path))
  obj <- readRDS(path)
  col <- paste0("score_", signature)
  qc <- load_qc %>% filter(sample_id == !!sample_id)
  n_spots <- if (nrow(qc) && "n_spots_matrix" %in% colnames(qc)) qc$n_spots_matrix[1] else NA_integer_
  median_genes <- if (nrow(qc) && "n_features_median" %in% colnames(qc)) qc$n_features_median[1] else NA_real_
  subtitle <- if (nrow(qc)) paste0("spots=", n_spots, "; median genes=", round(as.numeric(median_genes), 0)) else ""
  if (!col %in% colnames(obj@meta.data)) return(plot_missing(sample_id, signature))
  p <- tryCatch({
    SpatialFeaturePlot(
      obj, features = col, images = names(obj@images)[1],
      pt.size.factor = 1.6, alpha = c(0.35, 1),
      combine = FALSE
    )[[1]] +
      scale_fill_gradient(low = cols[1], high = cols[2])
  }, error = function(e) plot_missing(sample_id, paste(signature, "plot failed")))
  p +
    labs(title = sample_id, subtitle = subtitle, fill = "score") +
    theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 11, color = pal$ink),
      plot.subtitle = element_text(family = "Arial", size = 8.5, color = pal$slate),
      legend.position = "right",
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8)
    )
}

save_contact_sheet <- function(signatures, filename_stub, title, cols) {
  plots <- list()
  for (sig in signatures) {
    plots <- c(plots, lapply(objects, plot_one, signature = sig, cols = cols))
  }
  ncol <- length(objects)
  nrow <- length(signatures)
  sheet <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = title,
      subtitle = "Spatial projection; not causal or plaque-specific validation",
      theme = theme(
        plot.title = element_text(family = "Arial", face = "bold", size = 15, color = pal$ink),
        plot.subtitle = element_text(family = "Arial", size = 10, color = pal$terracotta)
      )
    )
  width <- 3.1 * ncol
  height <- 3.2 * nrow + 0.7
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), sheet, width = width, height = height, device = "pdf")
  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), sheet, width = width, height = height, dpi = 600)
}

risk_cols <- c(pal$pale_bg, pal$deep_teal)
context_cols <- c(pal$pale_bg, pal$slate)
p1_cols <- c(pal$pale_bg, pal$sand)

save_contact_sheet("MAGMA_top50", "gse206306_spatial_overlay_MAGMA_top50_v0.1",
                   "GSE206306 spatial overlay: MAGMA_top50", risk_cols)
save_contact_sheet("MAGMA_top100", "gse206306_spatial_overlay_MAGMA_top100_v0.1",
                   "GSE206306 spatial overlay: MAGMA_top100", risk_cols)
save_contact_sheet("MAGMA_FDR", "gse206306_spatial_overlay_MAGMA_FDR_v0.1",
                   "GSE206306 spatial overlay: MAGMA_FDR", risk_cols)
save_contact_sheet("P1_core", "gse206306_spatial_overlay_P1_core_v0.1",
                   "GSE206306 spatial overlay: P1_core", p1_cols)
save_contact_sheet("Loop_TAL", "gse206306_spatial_overlay_Loop_TAL_v0.1",
                   "GSE206306 spatial overlay: Loop_TAL", context_cols)
save_contact_sheet("Injury_remodeling", "gse206306_spatial_overlay_injury_remodeling_v0.1",
                   "GSE206306 spatial overlay: Injury_remodeling", context_cols)
save_contact_sheet(c("Fibrosis_ECM", "Immune_myeloid"), "gse206306_spatial_overlay_ecm_immune_v0.1",
                   "GSE206306 spatial overlay: ECM and immune/myeloid context", context_cols)
