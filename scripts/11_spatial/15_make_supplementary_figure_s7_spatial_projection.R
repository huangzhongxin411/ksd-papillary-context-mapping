suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
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
score_qc <- read_tsv(file.path(out_dir, "spatial_module_score_qc_v0.1.tsv"), show_col_types = FALSE)
summary <- read_tsv(file.path(out_dir, "spatial_projection_summary_v0.1.tsv"), show_col_types = FALSE)
p1 <- read_tsv(file.path(out_dir, "p1_gene_spatial_expression_audit_v0.1.tsv"), show_col_types = FALSE)

rep_sample <- load_qc %>% arrange(desc(n_spots_matrix)) %>% slice(1) %>% pull(sample_id)
rep_obj <- readRDS(file.path(obj_dir, paste0(rep_sample, "_spatial_scores_v0.1.rds")))

panel_label <- function(p, label) {
  p + labs(tag = label) +
    theme(plot.tag = element_text(family = "Arial", face = "bold", size = 15, color = pal$ink))
}

plot_spatial <- function(obj, feature, title, cols) {
  if (!feature %in% colnames(obj@meta.data)) {
    return(ggplot() + annotate("text", 0.5, 0.5, label = paste(title, "not scored"), size = 4) +
             xlim(0, 1) + ylim(0, 1) + theme_void())
  }
  SpatialFeaturePlot(obj, features = feature, images = names(obj@images)[1],
                     pt.size.factor = 1.8, alpha = c(0.35, 1),
                     combine = FALSE)[[1]] +
    scale_fill_gradient(low = cols[1], high = cols[2]) +
    labs(title = title, subtitle = rep_sample, fill = "score") +
    theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 12, color = pal$ink),
      plot.subtitle = element_text(family = "Arial", size = 9, color = pal$slate),
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8)
    )
}

panel_a <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.08, ymax = 0.92,
           fill = "white", color = pal$slate, linewidth = 0.4) +
  annotate("text", x = 0.08, y = 0.80, hjust = 0, label = "A. Loading and projection gate",
           family = "Arial", fontface = "bold", size = 4.7, color = pal$ink) +
  annotate("text", x = 0.08, y = 0.62, hjust = 0,
           label = paste0("5/5 GSE206306 samples loaded\nImages and coordinates loaded\nQC pass; projection allowed\nScored by deterministic mean-z"),
           family = "Arial", size = 3.6, color = pal$ink, lineheight = 1.12) +
  annotate("text", x = 0.08, y = 0.20, hjust = 0,
           label = "Boundary: spatial context mapping only",
           family = "Arial", fontface = "italic", size = 3.2, color = pal$terracotta) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

panel_b <- panel_label(plot_spatial(rep_obj, "score_MAGMA_top50", "MAGMA_top50", c(pal$pale_bg, pal$deep_teal)), "B")
panel_c <- panel_label(plot_spatial(rep_obj, "score_Loop_TAL", "Loop/TAL context", c(pal$pale_bg, pal$slate)), "C")
panel_d <- panel_label(plot_spatial(rep_obj, "score_Injury_remodeling", "Injury/remodeling", c(pal$pale_bg, pal$slate)), "D")

summary_e <- summary %>%
  filter(risk_signature %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR"),
         context_signature %in% c("Loop_TAL", "Injury_remodeling", "Fibrosis_ECM", "Immune_myeloid")) %>%
  mutate(pair = paste(risk_signature, context_signature, sep = " | "),
         context_signature = factor(context_signature, levels = c("Loop_TAL", "Injury_remodeling", "Fibrosis_ECM", "Immune_myeloid")))

panel_e <- ggplot(summary_e, aes(x = median_rho_adjusted, y = reorder(pair, median_rho_adjusted))) +
  geom_vline(xintercept = 0, color = pal$light_grey, linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = median_rho_adjusted, yend = reorder(pair, median_rho_adjusted)),
               color = pal$slate, linewidth = 0.7) +
  geom_point(aes(fill = support_class), shape = 21, size = 3.2, color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    strong_spatial_context_support = pal$deep_teal,
    moderate_spatial_context_support = pal$sand,
    no_consistent_support = pal$light_grey
  ), drop = FALSE) +
  labs(tag = "E", title = "Across-sample adjusted co-distribution",
       x = "Median residual Spearman rho", y = NULL, fill = "Support") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 15, color = pal$ink),
        plot.title = element_text(face = "bold", size = 12, color = pal$ink),
        axis.text.y = element_text(size = 7.8, color = pal$ink),
        legend.position = "bottom")

p1_tile <- p1 %>%
  mutate(status = case_when(
    spatial_plot_allowed ~ "resolved",
    gene_present ~ "low detection",
    TRUE ~ "not present"
  ))

panel_f <- ggplot(p1_tile, aes(x = sample_id, y = gene, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(resolved = pal$sand, `low detection` = pal$light_grey, `not present` = "white")) +
  labs(tag = "F", title = "P1 expression resolvability boundary",
       subtitle = "Low detection is not interpreted as biological absence",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.tag = element_text(face = "bold", size = 15, color = pal$ink),
        plot.title = element_text(face = "bold", size = 12, color = pal$ink),
        plot.subtitle = element_text(size = 8.5, color = pal$terracotta),
        axis.text.x = element_text(angle = 35, hjust = 1, color = pal$ink),
        axis.text.y = element_text(color = pal$ink),
        legend.position = "bottom")

top <- panel_a + panel_b + panel_c + panel_d + plot_layout(widths = c(1.05, 1, 1, 1))
bottom <- panel_e + panel_f + plot_layout(widths = c(1.35, 1))
fig <- top / bottom +
  plot_annotation(
    title = "GSE206306 spatial projection supports papillary context mapping but not causal or plaque-specific validation",
    subtitle = "Supplementary Figure S7 candidate; scores are spatial projections in QC-passed human renal papillary sections",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 16, color = pal$ink),
      plot.subtitle = element_text(family = "Arial", size = 10, color = pal$slate)
    )
  )

ggsave(file.path(fig_dir, "supplementary_figure_s7_gse206306_spatial_projection_v0.1.pdf"),
       fig, width = 15.5, height = 9.5, device = "pdf")
ggsave(file.path(fig_dir, "supplementary_figure_s7_gse206306_spatial_projection_v0.1.png"),
       fig, width = 15.5, height = 9.5, dpi = 600)
