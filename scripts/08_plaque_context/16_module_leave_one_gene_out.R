suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_set <- function(path) unique(fread(path, header = FALSE)[[1]])
sets <- list(
  MAGMA_top50 = read_set("results/gene_sets/magma_top50.txt"),
  MAGMA_top100 = read_set("results/gene_sets/magma_top100.txt"),
  MAGMA_FDR = read_set("results/gene_sets/magma_fdr05.txt"),
  MAGMA_suggestive = read_set("results/gene_sets/magma_suggestive_p1e4.txt")
)

meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
sample_ids <- intersect(meta$sample_id, names(mat))
meta <- meta[match(sample_ids, sample_id)]
expr <- as.matrix(mat[, ..sample_ids])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
zexpr <- t(scale(t(log2(expr + 1))))
zexpr[!is.finite(zexpr)] <- NA_real_
paired <- meta[, .N, by = .(patient_id, group_curated)][, .N, by = patient_id][N == 2, patient_id]

paired_delta_test <- function(genes) {
  genes <- intersect(genes, rownames(zexpr))
  sc <- colMeans(zexpr[genes, , drop = FALSE], na.rm = TRUE)
  dt <- data.table(sample_id = sample_ids, score = sc)
  dt <- merge(dt, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
  pw <- dcast(dt[patient_id %in% paired], patient_id ~ group_curated, value.var = "score", fun.aggregate = mean)
  pw <- pw[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
  delta <- mean(pw$plaque_or_stone_papilla - pw$control_or_adjacent, na.rm = TRUE)
  p <- if (nrow(pw) >= 3) t.test(pw$plaque_or_stone_papilla, pw$control_or_adjacent, paired = TRUE)$p.value else NA_real_
  list(delta = delta, p = p, n_pairs = nrow(pw))
}

orig <- fread(file.path(table_dir, "gse73680_patient_level_module_response.tsv"))
res <- rbindlist(lapply(names(sets), function(nm) {
  genes <- intersect(sets[[nm]], rownames(zexpr))
  original_effect <- orig[module_name == nm, paired_delta][1]
  rbindlist(lapply(genes, function(g) {
    remaining <- setdiff(genes, g)
    tt <- paired_delta_test(remaining)
    data.table(module_name = nm, removed_gene = g, n_genes_remaining = length(remaining),
               effect_size_after_removal = tt$delta, p_value_after_removal = tt$p,
               original_effect = original_effect,
               effect_retention_fraction = ifelse(abs(original_effect) > 0, tt$delta / original_effect, NA_real_),
               direction_preserved = sign(tt$delta) == sign(original_effect))
  }))
}), fill = TRUE)
res[, fdr_after_removal := p.adjust(p_value_after_removal, method = "BH"), by = module_name]
res[, interpretation := fcase(
  direction_preserved == TRUE & effect_retention_fraction >= 0.80, "robust",
  direction_preserved == TRUE & effect_retention_fraction >= 0.50, "moderate",
  default = "gene-driven_or_fragile"
)]
fwrite(res[, .(module_name, removed_gene, n_genes_remaining, effect_size_after_removal,
               p_value_after_removal, fdr_after_removal, effect_retention_fraction,
               direction_preserved, interpretation)],
       file.path(table_dir, "gse73680_module_leave_one_gene_out.tsv"), sep = "\t")

plot_dt <- res[, .(min_retention = min(effect_retention_fraction, na.rm = TRUE),
                   median_retention = median(effect_retention_fraction, na.rm = TRUE),
                   frac_direction_preserved = mean(direction_preserved, na.rm = TRUE)), by = module_name]
p <- ggplot(plot_dt, aes(module_name, min_retention, fill = frac_direction_preserved)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#777777") +
  geom_col(width = 0.65, color = "#555555", linewidth = 0.2) +
  geom_text(aes(label = sprintf("dir %.0f%%", 100 * frac_direction_preserved)), vjust = -0.35, size = 3) +
  scale_fill_gradient(low = "#C9C9C9", high = "#3E6672", limits = c(0, 1)) +
  labs(x = NULL, y = "Minimum effect retention after removing one gene", fill = "Direction preserved",
       title = "MAGMA module leave-one-gene-out sensitivity") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "gse73680_module_leave_one_gene_out.pdf"), p, width = 7.2, height = 4.8)
message("wrote module leave-one-gene-out sensitivity")
