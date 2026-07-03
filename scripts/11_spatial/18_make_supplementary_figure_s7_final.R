suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
})

fig_dir <- "results/figures/supplementary/spatial_projection/final"
out_dir <- "results/spatial/phase28c_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal <- list(
  deep_teal = "#245A64",
  loop_teal = "#0F4C5C",
  bluegrey = "#7F9DA6",
  sand_gold = "#B99B5A",
  terracotta = "#9B5C4D",
  pale_grey = "#E6E9EA",
  light_grey = "#D3DADC",
  ink = "#22313B"
)

load_qc <- read_tsv("results/spatial/phase28b_loading_rescue/spatial_load_qc_v0.2.tsv", show_col_types = FALSE)
summary <- read_tsv(file.path(out_dir, "spatial_projection_summary_v0.1.tsv"), show_col_types = FALSE)
p1 <- read_tsv(file.path(out_dir, "p1_gene_spatial_expression_audit_v0.1.tsv"), show_col_types = FALSE)

rep_sample <- "GSM7166170"
rep_obj <- readRDS(file.path(obj_dir, paste0(rep_sample, "_spatial_scores_v0.1.rds")))

spatial_panel <- function(feature, title, tag, cols) {
  p <- SpatialFeaturePlot(rep_obj, features = feature, images = names(rep_obj@images)[1],
                          pt.size.factor = 1.75, alpha = c(0.35, 1), combine = FALSE)[[1]] +
    scale_fill_gradient(low = cols[1], high = cols[2]) +
    labs(tag = tag, title = title, subtitle = rep_sample, fill = "score") +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 15, color = pal$ink),
      plot.title = element_text(family = "Arial", face = "bold", size = 12, color = pal$ink),
      plot.subtitle = element_text(family = "Arial", size = 9, color = pal$bluegrey),
      legend.title = element_text(family = "Arial", size = 8.5),
      legend.text = element_text(family = "Arial", size = 8)
    )
  p
}

panel_a <- ggplot() +
  annotate("rect", xmin = 0.03, xmax = 0.97, ymin = 0.08, ymax = 0.92,
           fill = "white", color = pal$bluegrey, linewidth = 0.45) +
  annotate("text", x = 0.10, y = 0.80, hjust = 0, label = "A. Loading / projection gate",
           family = "Arial", fontface = "bold", size = 4.7, color = pal$ink) +
  annotate("text", x = 0.10, y = 0.61, hjust = 0,
           label = "5/5 samples loaded\nImages and coordinates loaded\nScore projection completed\nDeterministic mean-z scoring",
           family = "Arial", size = 3.65, color = pal$ink, lineheight = 1.08) +
  annotate("text", x = 0.10, y = 0.20, hjust = 0,
           label = "Boundary: spatial context mapping only",
           family = "Arial", fontface = "italic", size = 3.2, color = pal$terracotta) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

panel_b <- spatial_panel("score_MAGMA_top50", "MAGMA_top50", "B", c(pal$pale_grey, pal$deep_teal))
panel_c <- spatial_panel("score_Loop_TAL", "Loop/TAL context", "C", c(pal$pale_grey, pal$bluegrey))
panel_d <- spatial_panel("score_Ion_mineral_handling", "Ion/mineral handling", "D", c(pal$pale_grey, pal$bluegrey))

pairs <- tibble(
  risk_signature = c("MAGMA_top50", "MAGMA_top100", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_FDR", "MAGMA_top50", "MAGMA_top100", "MAGMA_top50", "MAGMA_top100"),
  context_signature = c("Loop_TAL", "Loop_TAL", "Ion_mineral_handling", "Ion_mineral_handling", "Loop_TAL", "Ion_mineral_handling", "Injury_remodeling", "Injury_remodeling", "Fibrosis_ECM", "Fibrosis_ECM")
) %>%
  mutate(pair = paste(risk_signature, context_signature, sep = " | "))

panel_e_dat <- pairs %>%
  left_join(summary, by = c("risk_signature", "context_signature")) %>%
  mutate(
    pair = factor(pair, levels = rev(pair)),
    support_display = ifelse(support_class %in% c("moderate_spatial_context_support", "strong_spatial_context_support"),
                             "moderate support", "no consistent support")
  )

panel_e <- ggplot(panel_e_dat, aes(x = median_rho_adjusted, y = pair)) +
  geom_vline(xintercept = 0, color = pal$light_grey, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = median_rho_adjusted, yend = pair),
               color = pal$loop_teal, linewidth = 0.75) +
  geom_point(aes(fill = support_display), shape = 21, size = 3.4, color = "white", stroke = 0.35) +
  scale_fill_manual(values = c("moderate support" = pal$sand_gold, "no consistent support" = pal$light_grey)) +
  labs(tag = "E", title = "Across-sample adjusted co-distribution",
       subtitle = "Median residual Spearman rho across five samples; spot-level associations are descriptive",
       x = "Median residual Spearman rho", y = NULL, fill = NULL) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 15, color = pal$ink),
    plot.title = element_text(face = "bold", size = 12, color = pal$ink),
    plot.subtitle = element_text(size = 8.5, color = pal$bluegrey),
    axis.text.y = element_text(size = 8.2, color = pal$ink),
    legend.position = "bottom"
  )

gene_order <- c("CLDN10", "HIBADH", "PKD2", "UMOD", "CASR", "CLDN14")
p1_tile <- p1 %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    sample_id = factor(sample_id, levels = sort(unique(sample_id))),
    status = ifelse(spatial_plot_allowed, "resolved", "low detection")
  )

panel_f <- ggplot(p1_tile, aes(x = sample_id, y = gene, fill = status)) +
  geom_tile(color = "white", linewidth = 0.55) +
  scale_fill_manual(values = c("resolved" = pal$sand_gold, "low detection" = pal$light_grey)) +
  labs(tag = "F", title = "P1 expression resolvability boundary",
       subtitle = "Detection is interpreted as assay resolvability, not disease validation",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold", size = 15, color = pal$ink),
    plot.title = element_text(face = "bold", size = 12, color = pal$ink),
    plot.subtitle = element_text(size = 8.5, color = pal$terracotta),
    axis.text.x = element_text(angle = 35, hjust = 1, color = pal$ink, size = 8),
    axis.text.y = element_text(color = pal$ink, size = 8.5),
    legend.position = "bottom"
  )

footer <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "Spatial projection supports papillary context mapping only; causal, plaque-specific and gene-validation claims remain unresolved.",
           family = "Arial", size = 3.6, color = pal$terracotta) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

top <- panel_a + panel_b + panel_c + panel_d + plot_layout(widths = c(1.05, 1, 1, 1))
bottom <- panel_e + panel_f + plot_layout(widths = c(1.35, 1))
fig <- top / bottom / footer + plot_layout(heights = c(1.1, 1.05, 0.10)) +
  plot_annotation(
    title = "GSE206306 spatial projection provides papillary context support for MAGMA-prioritized modules",
    subtitle = "Spatial context mapping in QC-passed human renal papillary sections; not plaque-specific or causal validation",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 16, color = pal$ink),
      plot.subtitle = element_text(family = "Arial", size = 10, color = pal$bluegrey)
    )
  )

pdf_path <- file.path(fig_dir, "supplementary_figure_s7_gse206306_spatial_projection_final_v0.1.pdf")
png_path <- file.path(fig_dir, "supplementary_figure_s7_gse206306_spatial_projection_final_v0.1.png")
ggsave(png_path, fig, width = 15.5, height = 9.7, dpi = 600)
ggsave(pdf_path, fig, width = 15.5, height = 9.7, device = "pdf")

# The local R PDF device may omit raster spatial layers on some macOS setups.
# Save a robust image-backed PDF from the generated 600 dpi PNG.
if (requireNamespace("png", quietly = TRUE)) {
  img <- png::readPNG(png_path)
  grDevices::pdf(pdf_path, width = 15.5, height = 9.7)
  op <- par(mar = c(0, 0, 0, 0))
  plot.new()
  rasterImage(img, 0, 0, 1, 1)
  par(op)
  dev.off()
}
