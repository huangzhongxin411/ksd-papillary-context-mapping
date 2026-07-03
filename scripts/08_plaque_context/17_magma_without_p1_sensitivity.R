suppressPackageStartupMessages(library(data.table))

set.seed(73681)
table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
p1 <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
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
background <- rownames(zexpr)[rowSums(is.finite(zexpr)) >= ceiling(0.8 * ncol(zexpr))]
paired_delta <- function(genes) {
  genes <- intersect(genes, rownames(zexpr))
  sc <- colMeans(zexpr[genes, , drop = FALSE], na.rm = TRUE)
  dt <- merge(data.table(sample_id = sample_ids, score = sc), meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
  pw <- dcast(dt[patient_id %in% paired], patient_id ~ group_curated, value.var = "score", fun.aggregate = mean)
  pw <- pw[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
  delta <- mean(pw$plaque_or_stone_papilla - pw$control_or_adjacent, na.rm = TRUE)
  p <- if (nrow(pw) >= 3) t.test(pw$plaque_or_stone_papilla, pw$control_or_adjacent, paired = TRUE)$p.value else NA_real_
  list(delta = delta, p = p, n = length(genes))
}
orig <- fread(file.path(table_dir, "gse73680_patient_level_module_response.tsv"))
res <- rbindlist(lapply(names(sets), function(nm) {
  original_genes <- intersect(sets[[nm]], rownames(zexpr))
  genes <- setdiff(original_genes, p1)
  tt <- paired_delta(genes)
  rand <- replicate(1000, paired_delta(sample(background, length(genes)))$delta)
  data.table(module_name = paste0(nm, "_without_P1"),
             original_n_genes = length(original_genes),
             n_p1_removed = length(intersect(original_genes, p1)),
             n_genes_remaining = length(genes),
             original_effect = orig[module_name == nm, paired_delta][1],
             effect_without_p1 = tt$delta,
             p_value_without_p1 = tt$p,
             random_benchmark_percentile_without_p1 = mean(rand <= tt$delta),
             empirical_p_without_p1 = mean(abs(rand) >= abs(tt$delta)))
}), fill = TRUE)
res[, fdr_without_p1 := p.adjust(p_value_without_p1, method = "BH")]
res[, interpretation := fcase(
  fdr_without_p1 <= 0.05 & random_benchmark_percentile_without_p1 >= 0.95, "robust_without_P1",
  p_value_without_p1 < 0.05, "nominal_without_P1",
  default = "not_supported_without_P1"
)]
fwrite(res, file.path(table_dir, "gse73680_magma_without_p1_sensitivity.tsv"), sep = "\t")
message("wrote MAGMA without-P1 sensitivity")
