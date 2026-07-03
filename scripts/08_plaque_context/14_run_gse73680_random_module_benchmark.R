suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(73680)
table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

target_modules <- c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")
meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
sample_ids <- intersect(meta$sample_id, names(mat))
meta <- meta[match(sample_ids, sample_id)]
expr <- as.matrix(mat[, ..sample_ids])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
expr <- log2(expr + 1)
zexpr <- t(scale(t(expr)))
zexpr[!is.finite(zexpr)] <- NA_real_
background <- rownames(zexpr)[rowSums(is.finite(zexpr)) >= ceiling(0.8 * ncol(zexpr))]

scores <- fread(file.path(table_dir, "gse73680_module_score_matrix.tsv"))
resp <- fread(file.path(table_dir, "gse73680_module_score_response.tsv"))
res <- rbindlist(lapply(target_modules, function(mn) {
  row <- resp[module_name == mn]
  n <- row$n_genes_detected[1]
  obs <- row$delta[1]
  if (!length(n) || is.na(n) || n < 2 || n > length(background) || is.na(obs)) {
    return(data.table(module_name = mn, n_genes_detected = n, observed_delta = obs,
                      n_random = 0L, empirical_p = NA_real_, percentile = NA_real_,
                      benchmark_interpretation = "not_interpretable"))
  }
  rand <- replicate(1000, {
    genes <- sample(background, n)
    sc <- colMeans(zexpr[genes, , drop = FALSE], na.rm = TRUE)
    mean(sc[meta$group_curated == "plaque_or_stone_papilla"], na.rm = TRUE) -
      mean(sc[meta$group_curated == "control_or_adjacent"], na.rm = TRUE)
  })
  empirical_p <- mean(abs(rand) >= abs(obs))
  percentile <- mean(rand <= obs)
  data.table(module_name = mn, n_genes_detected = n, observed_delta = obs, n_random = length(rand),
             empirical_p = empirical_p, percentile = percentile,
             benchmark_interpretation = fcase(
               empirical_p < 0.05 | percentile >= 0.95 | percentile <= 0.05, "module response exceeds random expectation",
               percentile >= 0.90 | percentile <= 0.10, "moderate exploratory random-benchmark signal",
               default = "not stronger than random background"
             ),
             random_delta_mean = mean(rand), random_delta_sd = sd(rand))
}), fill = TRUE)
fwrite(res, file.path(table_dir, "gse73680_random_module_benchmark.tsv"), sep = "\t")

plot_dt <- res[!is.na(observed_delta)]
p <- ggplot(plot_dt, aes(module_name, observed_delta, fill = benchmark_interpretation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#8A8A8A") +
  geom_col(width = 0.62, color = "#555555", linewidth = 0.2) +
  geom_errorbar(aes(ymin = random_delta_mean - 1.96 * random_delta_sd,
                    ymax = random_delta_mean + 1.96 * random_delta_sd), width = 0.2, color = "#333333") +
  scale_fill_manual(values = c("module response exceeds random expectation" = "#3E6672",
                               "moderate exploratory random-benchmark signal" = "#B08A45",
                               "not stronger than random background" = "#C9C9C9",
                               "not_interpretable" = "#E0E0E0")) +
  labs(x = NULL, y = "Observed disease-context delta", fill = "Benchmark",
       title = "Random gene-set benchmark for GSE73680 module responses",
       subtitle = "Error bars show mean +/- 1.96 SD of random module deltas") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
ggsave(file.path(fig_dir, "gse73680_random_module_benchmark.pdf"), p, width = 8.2, height = 5.2)
message("wrote GSE73680 random module benchmark")
