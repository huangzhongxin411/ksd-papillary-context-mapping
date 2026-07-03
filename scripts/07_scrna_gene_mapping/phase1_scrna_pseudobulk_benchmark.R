suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
table_dir <- "results/tables"
figure_dir <- "results/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260616)
obj <- readRDS(obj_path)
mat <- GetAssayData(obj, assay = "RNA", layer = "data")
meta <- obj@meta.data

top50 <- scan("results/tables/phase1_candidate_genes_top50.txt", what = character(), quiet = TRUE)
top100 <- scan("results/tables/phase1_candidate_genes_top100.txt", what = character(), quiet = TRUE)
ranked <- fread("results/tables/phase1_candidate_genes_ranked_protein_coding.tsv")

score_genes <- function(genes) {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[present, , drop = FALSE])
}

meta$risk_top50_score <- score_genes(top50)
meta$risk_top100_score <- score_genes(top100)

donor_score <- as.data.table(meta)[
  ,
  .(
    n_cells = .N,
    mean_top50_score = mean(risk_top50_score, na.rm = TRUE),
    median_top50_score = median(risk_top50_score, na.rm = TRUE),
    mean_top100_score = mean(risk_top100_score, na.rm = TRUE),
    median_top100_score = median(risk_top100_score, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE)
  ),
  by = .(sample_id, disease_status, phase1_cell_type)
]
setnames(donor_score, "phase1_cell_type", "cell_type")
fwrite(donor_score, file.path(table_dir, "scrna_module_score_by_donor_celltype.tsv"), sep = "\t")

candidate <- ranked[rank_in_candidate_list <= 100, .(gene, rank_in_candidate_list)]
candidate <- candidate[gene %in% rownames(mat)]
det_rows <- lapply(candidate$gene, function(g) {
  vals <- mat[g, ] > 0
  dt <- as.data.table(meta[, c("sample_id", "disease_status", "phase1_cell_type")])
  dt[, detected := as.numeric(vals)]
  out <- dt[, .(pct_expressed = mean(detected), n_cells = .N), by = .(sample_id, disease_status, phase1_cell_type)]
  out[, gene := g]
  out
})
det <- rbindlist(det_rows)
setnames(det, "phase1_cell_type", "cell_type")
det <- merge(det, candidate, by = "gene", all.x = TRUE)
setcolorder(det, c("gene", "rank_in_candidate_list", "sample_id", "disease_status", "cell_type", "pct_expressed", "n_cells"))
fwrite(det, file.path(table_dir, "scrna_candidate_gene_detection_rate.tsv"), sep = "\t")

cell_type <- factor(meta$phase1_cell_type)
groups <- levels(cell_type)
group_mat <- sparseMatrix(
  i = seq_along(cell_type),
  j = as.integer(cell_type),
  x = 1,
  dims = c(length(cell_type), length(groups)),
  dimnames = list(colnames(mat), groups)
)
group_counts <- Matrix::colSums(group_mat)
gene_by_group <- mat %*% group_mat
gene_by_group <- t(t(gene_by_group) / as.numeric(group_counts))

global_mean <- Matrix::rowMeans(mat)
genes_universe <- names(global_mean)[is.finite(global_mean) & global_mean > 0]
bins <- cut(
  global_mean[genes_universe],
  breaks = unique(quantile(global_mean[genes_universe], probs = seq(0, 1, length.out = 21), na.rm = TRUE)),
  include.lowest = TRUE
)
gene_bins <- split(genes_universe, bins)
bin_lookup <- setNames(as.character(bins), genes_universe)

sample_matched <- function(real_genes) {
  real_genes <- intersect(real_genes, names(bin_lookup))
  sampled <- character()
  for (g in real_genes) {
    pool <- setdiff(gene_bins[[bin_lookup[[g]]]], real_genes)
    if (length(pool) == 0) pool <- setdiff(genes_universe, real_genes)
    sampled <- c(sampled, sample(pool, 1))
  }
  unique(sampled)
}

benchmark_one <- function(real_genes, label, n_iter = 1000) {
  real_genes <- intersect(real_genes, rownames(gene_by_group))
  real_scores <- Matrix::colMeans(gene_by_group[real_genes, , drop = FALSE])
  rand_scores <- matrix(NA_real_, nrow = n_iter, ncol = length(groups), dimnames = list(NULL, groups))
  for (i in seq_len(n_iter)) {
    rg <- sample_matched(real_genes)
    rg <- intersect(rg, rownames(gene_by_group))
    rand_scores[i, ] <- Matrix::colMeans(gene_by_group[rg, , drop = FALSE])
  }
  rbindlist(lapply(groups, function(ct) {
    rs <- rand_scores[, ct]
    data.table(
      module_group = label,
      cell_type = ct,
      n_real_genes = length(real_genes),
      real_score = as.numeric(real_scores[ct]),
      random_mean = mean(rs, na.rm = TRUE),
      random_sd = sd(rs, na.rm = TRUE),
      random_p95 = quantile(rs, 0.95, na.rm = TRUE),
      percentile = mean(rs <= as.numeric(real_scores[ct]), na.rm = TRUE)
    )
  }))
}

bench <- rbind(
  benchmark_one(top50, "top50", n_iter = 1000),
  benchmark_one(top100, "top100", n_iter = 1000)
)
fwrite(bench, file.path(table_dir, "scrna_random_gene_set_benchmark.tsv"), sep = "\t")

pdf(file.path(figure_dir, "scrna_random_benchmark_celltype.pdf"), width = 9, height = 4.8)
print(
  ggplot(bench, aes(x = reorder(cell_type, percentile), y = percentile, fill = module_group)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(x = NULL, y = "Expression-matched random set percentile", fill = "Module") +
    theme_bw(base_size = 10)
)
dev.off()

cat("wrote\t", file.path(table_dir, "scrna_module_score_by_donor_celltype.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "scrna_candidate_gene_detection_rate.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "scrna_random_gene_set_benchmark.tsv"), "\n", sep = "")

