suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- c(
  P1_core = "#B59A5B",
  candidate_TAL_driver = "#3E6672",
  supporting_TAL_expressed_gene = "#6F8F98",
  non_TAL_or_low_detection = "#C9C9C9",
  not_detected_in_seurat = "#E1E1E1"
)

scale01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(1, length(x)))
  (x - rng[1]) / diff(rng)
}

drivers <- fread("results/tables/magma_top50_tal_driver_genes.tsv")
p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")

drivers <- drivers[present_in_seurat == TRUE]
drivers[, magma_strength := -log10(pmax(magma_p, .Machine$double.xmin))]
drivers[, tal_specificity_score := fifelse(
  is.na(all_celltype_max) | all_celltype_max <= 0,
  NA_real_,
  TAL_avg_expression / all_celltype_max
)]
drivers[, donor_detection := TAL_pct_expressed]
drivers[, magma_strength_scaled := scale01(magma_strength)]
drivers[, tal_specificity_scaled := scale01(tal_specificity_score)]
drivers[, donor_detection_scaled := scale01(donor_detection)]
drivers[, contribution_score := magma_strength_scaled * tal_specificity_scaled * donor_detection_scaled]
drivers[, candidate_role := fcase(
  gene %in% p1$gene, "P1_core",
  driver_class == "candidate_TAL_driver", "candidate_TAL_driver",
  driver_class == "supporting_TAL_expressed_gene", "supporting_TAL_expressed_gene",
  default = "non_TAL_or_low_detection"
)]
drivers[, gene_set := "MAGMA_top50"]
drivers[, contribution_rank := frank(-contribution_score, ties.method = "first")]
drivers[, gene_label := factor(gene, levels = rev(gene[order(contribution_rank)]))]

out <- drivers[, .(
  gene,
  MAGMA_p = magma_p,
  MAGMA_rank = magma_rank,
  MAGMA_fdr = magma_fdr,
  gene_set,
  locus_id,
  Loop_TAL_avg_expr = TAL_avg_expression,
  Loop_TAL_pct_expressed = TAL_pct_expressed,
  Loop_TAL_specificity = tal_specificity_score,
  donor_detection,
  leave_one_gene_delta,
  contribution_score,
  contribution_rank,
  candidate_role,
  top_celltype,
  claim_boundary = "Influential expression-context contributor; not a causal driver gene"
)]
setorder(out, contribution_rank)
fwrite(out, "results/tables/loop_tal_influential_magma_genes.tsv", sep = "\t")

top_plot <- out[is.finite(contribution_score)][order(contribution_rank)][1:min(.N, 18)]
top_plot[, gene_label := factor(gene, levels = rev(gene))]
top_plot[, label_p1 := ifelse(gene %in% p1$gene, "P1", "")]

fig <- ggplot(top_plot, aes(contribution_score, gene_label)) +
  geom_segment(aes(x = 0, xend = contribution_score, yend = gene_label),
               color = "#A8A8A8", linewidth = 0.45) +
  geom_point(aes(fill = candidate_role, size = donor_detection),
             shape = 21, color = "#444444", stroke = 0.25) +
  geom_text(aes(label = label_p1), nudge_x = 0.025, size = 2.6, fontface = "bold", color = "#4F4F4F") +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_size_continuous(range = c(2.4, 5.4)) +
  labs(
    title = "Loop/TAL influential MAGMA genes",
    subtitle = "Contribution score = scaled MAGMA strength x scaled Loop/TAL specificity x donor detection",
    x = "Contribution score",
    y = NULL,
    fill = "Role",
    size = "TAL detection"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "#555555"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/figure2_loop_tal_influential_genes_lollipop.pdf", fig,
       width = 8.2, height = 6.4, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure2_loop_tal_influential_genes_lollipop.png", fig,
       width = 8.2, height = 6.4, units = "in", dpi = 260, bg = "white")

heat <- melt(
  drivers[gene %in% top_plot$gene,
          .(gene, magma_strength_scaled, tal_specificity_scaled, donor_detection_scaled, contribution_score)],
  id.vars = "gene",
  variable.name = "metric",
  value.name = "scaled_value"
)
heat[, metric := factor(metric,
                        levels = c("magma_strength_scaled", "tal_specificity_scaled", "donor_detection_scaled", "contribution_score"),
                        labels = c("MAGMA strength", "Loop/TAL specificity", "Donor detection", "Combined score"))]
heat[, gene := factor(gene, levels = rev(top_plot$gene))]

fig_heat <- ggplot(heat, aes(metric, gene, fill = scaled_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", scaled_value)), size = 2.3, color = "#303030") +
  scale_fill_gradient(low = "#EEF3F4", high = "#3E6672", na.value = "#F1F1F1") +
  labs(title = "Contribution-score components", x = NULL, y = NULL, fill = "Scaled\nvalue") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        panel.grid = element_blank())

ggsave("results/figures/figure2_loop_tal_gene_contribution_heatmap.pdf", fig_heat,
       width = 7.8, height = 6.2, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure2_loop_tal_gene_contribution_heatmap.png", fig_heat,
       width = 7.8, height = 6.2, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# Loop/TAL Influential MAGMA Gene Analysis v0.1",
  "",
  "This analysis ranks detected MAGMA top-50 genes by a conservative expression-context contribution score:",
  "",
  "contribution score = scaled MAGMA strength x scaled Loop/TAL specificity x donor detection.",
  "",
  "The score is intended to identify genes contributing to the observed Loop/TAL-associated signal. It should not be interpreted as evidence of causal driver genes."
), "docs/loop_tal_influential_gene_analysis_v0.1.md")

message("wrote Loop/TAL influential gene analysis")
