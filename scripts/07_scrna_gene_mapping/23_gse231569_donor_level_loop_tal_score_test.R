suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

scores <- fread("results/tables/gse231569_donor_celltype_magma_scores_v0.2.tsv")
module_order <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")
scores <- scores[module_name %in% module_order]

label_cell <- function(x) {
  map <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like",
    Pericyte_smooth_muscle = "Perivascular/mural-like"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}
scores[, cell_label := label_cell(audited_broad_cell_type)]
scores[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
scores[, module_label := factor(module_name, levels = module_order,
                                labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]

rank_dt <- scores[, {
  d <- .SD[order(-mean_score)]
  d[, rank_desc := seq_len(.N)]
  d
}, by = .(module_name, donor_id)]
tal_rank <- rank_dt[is_tal == TRUE, .(
  n_donors = .N,
  median_tal_rank = median(rank_desc, na.rm = TRUE),
  max_tal_rank = max(rank_desc, na.rm = TRUE),
  n_rank1 = sum(rank_desc == 1, na.rm = TRUE),
  n_rank_top2 = sum(rank_desc <= 2, na.rm = TRUE),
  tal_mean_score = mean(mean_score, na.rm = TRUE)
), by = module_name]

contrast <- scores[, {
  tal <- mean_score[is_tal == TRUE]
  other <- mean_score[is_tal == FALSE]
  data.table(
    n_donors_with_tal = length(tal),
    tal_median = median(tal, na.rm = TRUE),
    other_median = median(other, na.rm = TRUE),
    tal_minus_other_median = median(tal, na.rm = TRUE) - median(other, na.rm = TRUE),
    wilcox_p = if (length(tal) >= 2 && length(other) >= 2) suppressWarnings(wilcox.test(tal, other)$p.value) else NA_real_
  )
}, by = module_name]
stats <- merge(tal_rank, contrast, by = "module_name", all = TRUE)
stats[, fdr := p.adjust(wilcox_p, method = "BH")]
stats[, interpretation := fcase(
  median_tal_rank <= 2 & tal_minus_other_median > 0, "donor-level Loop/TAL score support",
  default = "limited donor-level Loop/TAL support"
)]
fwrite(stats, "results/tables/gse231569_donor_celltype_magma_score_stats.tsv", sep = "\t")

p <- ggplot(scores, aes(cell_label, mean_score, fill = is_tal)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "#555555", linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 1.7,
             shape = 21, color = "#555555", stroke = 0.2) +
  facet_wrap(~ module_label, ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#3E6672", "FALSE" = "#D8D8D8"), labels = c("Other", "Loop/TAL")) +
  labs(title = "Donor-level MAGMA module score support for Loop/TAL context",
       subtitle = "Each point represents one donor-cell-type summary; cell-level UMAP is visualization only",
       x = NULL, y = "Mean per-cell module score", fill = NULL) +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555555"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave("results/figures/figure2_donor_level_magma_score_boxplot.pdf", p,
       width = 9.2, height = 6.5, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure2_donor_level_magma_score_boxplot.png", p,
       width = 9.2, height = 6.5, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# GSE231569 Donor-Level Loop/TAL MAGMA Score Test v0.1",
  "",
  "Cell-level UMAPs are used for visualization. Donor-cell-type summaries provide the analysis unit for Loop/TAL module-score support."
), "docs/gse231569_donor_level_magma_score_test_v0.1.md", useBytes = TRUE)

message("wrote GSE231569 donor-level MAGMA score outputs")
