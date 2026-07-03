suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(73682)
table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_set <- function(path, fallback = character()) if (file.exists(path)) unique(fread(path, header = FALSE)[[1]]) else fallback
sets <- list(
  P1_core_TAL_candidates = c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2"),
  MAGMA_top50 = read_set("results/gene_sets/magma_top50.txt"),
  MAGMA_top100 = read_set("results/gene_sets/magma_top100.txt"),
  MAGMA_FDR = read_set("results/gene_sets/magma_fdr05.txt"),
  MAGMA_suggestive = read_set("results/gene_sets/magma_suggestive_p1e4.txt")
)

meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
paired_ids <- meta[, .N, by = .(patient_id, group_curated)][, .N, by = patient_id][N == 2, patient_id]
paired_meta <- meta[patient_id %in% paired_ids]
control_cols <- paired_meta[group_curated == "control_or_adjacent"][order(patient_id), sample_id]
plaque_cols <- paired_meta[group_curated == "plaque_or_stone_papilla"][order(patient_id), sample_id]
paired_ids <- paired_meta[group_curated == "control_or_adjacent"][order(patient_id), patient_id]
all_cols <- meta$sample_id

mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
expr <- as.matrix(mat[, ..all_cols])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
logexpr <- log2(expr + 1)
zexpr <- t(scale(t(logexpr)))
zexpr[!is.finite(zexpr)] <- NA_real_

delta_mat <- zexpr[, plaque_cols, drop = FALSE] - zexpr[, control_cols, drop = FALSE]
colnames(delta_mat) <- paired_ids
gene_delta <- rowMeans(delta_mat, na.rm = TRUE)
gene_mean <- rowMeans(logexpr, na.rm = TRUE)
background <- names(gene_delta)[is.finite(gene_delta) & is.finite(gene_mean)]
qs <- unique(quantile(gene_mean[background], probs = seq(0, 1, length.out = 21), na.rm = TRUE))
bin_dt <- data.table(gene = background,
                     bin = cut(gene_mean[background], breaks = qs, include.lowest = TRUE, labels = FALSE))

module_delta <- function(genes) {
  genes <- intersect(genes, names(gene_delta))
  mean(gene_delta[genes], na.rm = TRUE)
}
sample_matched <- function(genes) {
  genes <- intersect(genes, bin_dt$gene)
  sampled <- vapply(genes, function(g) {
    b <- bin_dt[gene == g, bin][1]
    pool <- setdiff(bin_dt[bin == b, gene], genes)
    if (!length(pool)) pool <- bin_dt[bin == b, gene]
    sample(pool, 1)
  }, character(1))
  unique(sampled)
}

res <- rbindlist(lapply(names(sets), function(nm) {
  genes <- intersect(sets[[nm]], names(gene_delta))
  obs <- module_delta(genes)
  rand <- replicate(1000, module_delta(sample_matched(genes)))
  pct <- mean(rand <= obs)
  emp <- mean(abs(rand) >= abs(obs))
  data.table(module_name = nm, n_genes_detected = length(genes), observed_delta = obs,
             n_random = length(rand), empirical_p = emp,
             expression_matched_percentile = pct,
             random_delta_mean = mean(rand), random_delta_sd = sd(rand),
             interpretation = fcase(
               pct >= 0.95 | emp < 0.05, "robust_beyond_expression_level_background",
               pct >= 0.90, "moderate_expression_matched_signal",
               default = "may_reflect_expression_level_background"
             ))
}), fill = TRUE)
fwrite(res, file.path(table_dir, "gse73680_expression_matched_random_benchmark.tsv"), sep = "\t")

p <- ggplot(res, aes(expression_matched_percentile, module_name, fill = interpretation)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "#888888") +
  geom_col(width = 0.65, color = "#555555", linewidth = 0.2) +
  geom_text(aes(label = fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))),
            hjust = 1.05, size = 2.7) +
  coord_cartesian(xlim = c(0, 1.05)) +
  scale_fill_manual(values = c(robust_beyond_expression_level_background = "#3E6672",
                               moderate_expression_matched_signal = "#B08A45",
                               may_reflect_expression_level_background = "#C9C9C9")) +
  labs(x = "Observed paired delta percentile vs expression-matched random modules",
       y = NULL, fill = "Interpretation", title = "Expression-matched random benchmark") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "gse73680_expression_matched_random_benchmark.pdf"), p, width = 8, height = 5)
message("wrote expression-matched random benchmark")
