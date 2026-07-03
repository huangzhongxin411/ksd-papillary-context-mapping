suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
p1 <- fread(file.path(table_dir, "gse73680_patient_level_p1_gene_response.tsv"))$gene
module_scores <- fread(file.path(table_dir, "gse73680_patient_level_module_score_matrix.tsv"))
module_w <- dcast(module_scores, module_name + patient_id ~ group_curated, value.var = "patient_level_module_score")
module_delta <- module_w[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla),
                         .(module_or_gene = module_name, patient_id,
                           delta = plaque_or_stone_papilla - control_or_adjacent, feature_type = "module")]

mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[gene %in% p1]
sample_ids <- intersect(meta$sample_id, names(mat))
expr <- as.matrix(mat[, ..sample_ids])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
logexpr <- log2(expr + 1)
gene_long <- melt(as.data.table(logexpr, keep.rownames = "module_or_gene"), id.vars = "module_or_gene",
                  variable.name = "sample_id", value.name = "score")
gene_long <- merge(gene_long, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
gene_w <- dcast(gene_long, module_or_gene + patient_id ~ group_curated, value.var = "score", fun.aggregate = mean)
gene_delta <- gene_w[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla),
                     .(module_or_gene, patient_id, delta = plaque_or_stone_papilla - control_or_adjacent, feature_type = "gene")]
all_delta <- rbind(module_delta, gene_delta, fill = TRUE)
res <- all_delta[, {
  n_pos <- sum(delta > 0, na.rm = TRUE)
  n_neg <- sum(delta < 0, na.rm = TRUE)
  n <- sum(delta != 0 & is.finite(delta))
  p <- if (n > 0) binom.test(n_pos, n, p = 0.5, alternative = "greater")$p.value else NA_real_
  .(n_paired_patients = n, n_positive_delta = n_pos, n_negative_delta = n_neg,
    positive_fraction = n_pos / max(1, n), median_delta = median(delta, na.rm = TRUE),
    sign_test_p = p)
}, by = .(feature_type, module_or_gene)]
res[, interpretation := fcase(
  positive_fraction >= 0.70, "directionally_consistent",
  positive_fraction >= 0.60, "moderate_direction_consistency",
  default = "weak_or_heterogeneous"
)]
fwrite(res, file.path(table_dir, "gse73680_paired_direction_consistency.tsv"), sep = "\t")

plot_modules <- c("MAGMA_top50", "MAGMA_FDR", "MAGMA_suggestive")
sp <- module_scores[module_name %in% plot_modules]
sp[, group_curated := factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))]
p <- ggplot(sp, aes(group_curated, patient_level_module_score, group = patient_id)) +
  geom_line(color = "#6A7378", alpha = 0.45, linewidth = 0.35) +
  geom_point(aes(color = group_curated), size = 1.8, alpha = 0.85) +
  facet_wrap(~ module_name, nrow = 1) +
  scale_color_manual(values = c(control_or_adjacent = "#8AA0A8", plaque_or_stone_papilla = "#B08A45")) +
  labs(x = NULL, y = "Patient-level module score", color = NULL,
       title = "Paired patient module-score shifts") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "gse73680_magma_module_paired_delta_spaghetti.pdf"), p, width = 8.5, height = 3.6)
message("wrote paired direction consistency")
