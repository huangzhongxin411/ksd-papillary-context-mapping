suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

set.seed(20260617)

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
cfg_path <- "config/gse231569_cluster_to_celltype_audited.tsv"
benchmark_path <- "results/tables/audited_locus_scrna_random_benchmark.tsv"
cell_score_path <- "results/tables/audited_locus_scrna_module_score_by_celltype.tsv"
table_dir <- "results/tables"
figure_dir <- "results/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- fread(cfg_path)
cfg[, cluster := as.character(cluster)]
bench <- fread(benchmark_path)
cell_scores <- fread(cell_score_path)

cell_meta <- cfg[
  ,
  .(
    n_cells = sum(n_cluster_cells),
    n_donors = max(n_celltype_donors, na.rm = TRUE),
    annotation_confidence = paste(sort(unique(annotation_confidence)), collapse = ";"),
    immune_review_flag = any(immune_review_flag),
    mural_review_flag = any(mural_review_flag)
  ),
  by = .(cell_type = broad_cell_type)
]

score_long <- melt(
  cell_scores,
  id.vars = c(
    "audited_broad_cell_type",
    "n_cells",
    "n_donors",
    "annotation_confidence",
    "immune_review_cells",
    "mural_review_cells"
  ),
  measure.vars = patterns("^mean_audited_"),
  variable.name = "score_name",
  value.name = "mean_module_score"
)
score_long[, gene_set := sub("^mean_audited_(.*)_score$", "\\1", score_name)]
score_long <- score_long[, .(gene_set, cell_type = audited_broad_cell_type, mean_module_score)]

evidence <- merge(
  bench[, .(
    gene_set,
    cell_type = audited_broad_cell_type,
    benchmark_percentile,
    n_real_genes
  )],
  score_long,
  by = c("gene_set", "cell_type"),
  all.x = TRUE
)
evidence <- merge(evidence, cell_meta, by = "cell_type", all.x = TRUE)

evidence[, evidence_class := fifelse(
  mural_review_flag | annotation_confidence == "low_or_exploratory" | n_cells < 200 | n_donors < 4,
  "exploratory_low_abundance",
  fifelse(
    immune_review_flag,
    "review_flagged",
    fifelse(
      benchmark_percentile >= 0.95,
      "strong_primary",
      fifelse(benchmark_percentile >= 0.90, "moderate_primary", "not_supported")
    )
  )
)]

evidence[cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95,
  evidence_class := "strong_primary"
]
evidence[cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.90 & benchmark_percentile < 0.95,
  evidence_class := "moderate_primary"
]
evidence[immune_review_flag == TRUE & benchmark_percentile >= 0.90 & !mural_review_flag,
  evidence_class := "review_flagged"
]

evidence[, interpretation := fifelse(
  evidence_class == "strong_primary",
  "stable high-percentile signal after audited annotation",
  fifelse(
    evidence_class == "moderate_primary",
    "moderate audited signal requiring confirmation by MAGMA-prioritized genes",
    fifelse(
      evidence_class == "exploratory_low_abundance",
      "high or notable signal but limited by low cell count, donor coverage, or exploratory annotation",
      fifelse(
        evidence_class == "review_flagged",
        "primary label retained but immune-review clusters require caution",
        "not supported by expression-matched benchmark"
      )
    )
  )
)]

setorder(evidence, gene_set, -benchmark_percentile)
fwrite(
  evidence,
  file.path(table_dir, "audited_locus_celltype_evidence_summary.tsv"),
  sep = "\t"
)

obj <- readRDS(obj_path)
cluster_to_broad <- setNames(cfg$broad_cell_type, cfg$cluster)
cluster_to_conf <- setNames(cfg$annotation_confidence, cfg$cluster)
cluster_to_immune <- setNames(cfg$immune_review_flag, cfg$cluster)
cluster_to_mural <- setNames(cfg$mural_review_flag, cfg$cluster)
obj$audited_broad_cell_type <- unname(cluster_to_broad[as.character(obj$seurat_clusters)])
obj$audited_annotation_confidence <- unname(cluster_to_conf[as.character(obj$seurat_clusters)])
obj$audited_immune_review_flag <- unname(cluster_to_immune[as.character(obj$seurat_clusters)])
obj$audited_mural_review_flag <- unname(cluster_to_mural[as.character(obj$seurat_clusters)])
if (!"donor_id" %in% colnames(obj@meta.data)) obj$donor_id <- obj$sample_id

keep_cells <- colnames(obj)[
  obj$audited_annotation_confidence != "low_or_exploratory" &
    !(obj$audited_immune_review_flag %in% TRUE)
]
obj_cons <- subset(obj, cells = keep_cells)
mat <- GetAssayData(obj_cons, assay = "RNA", layer = "data")

gene_sets <- list(
  locus_top50 = scan("results/gene_sets/locus_top50.txt", what = character(), quiet = TRUE),
  locus_top100 = scan("results/gene_sets/locus_top100.txt", what = character(), quiet = TRUE)
)

cell_type <- factor(obj_cons$audited_broad_cell_type)
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
bin_breaks <- unique(quantile(global_mean[genes_universe], probs = seq(0, 1, length.out = 21), na.rm = TRUE))
bins <- cut(global_mean[genes_universe], breaks = bin_breaks, include.lowest = TRUE)
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

score_genes <- function(genes) {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[present, , drop = FALSE])
}

meta_cons <- as.data.table(obj_cons@meta.data)
for (nm in names(gene_sets)) {
  meta_cons[[paste0(nm, "_score")]] <- score_genes(gene_sets[[nm]])
}

sensitivity_module <- meta_cons[
  ,
  .(
    n_cells = .N,
    n_donors = uniqueN(donor_id),
    mean_locus_top50_score = mean(locus_top50_score, na.rm = TRUE),
    mean_locus_top100_score = mean(locus_top100_score, na.rm = TRUE)
  ),
  by = .(audited_broad_cell_type)
]

benchmark_one <- function(real_genes, label, n_iter = 1000) {
  real_genes <- intersect(real_genes, rownames(gene_by_group))
  real_scores <- Matrix::colMeans(gene_by_group[real_genes, , drop = FALSE])
  rand_scores <- matrix(NA_real_, nrow = n_iter, ncol = length(groups), dimnames = list(NULL, groups))
  for (i in seq_len(n_iter)) {
    rg <- intersect(sample_matched(real_genes), rownames(gene_by_group))
    rand_scores[i, ] <- Matrix::colMeans(gene_by_group[rg, , drop = FALSE])
  }
  rbindlist(lapply(groups, function(ct) {
    rs <- rand_scores[, ct]
    data.table(
      analysis_version = "conservative_exclude_low_or_exploratory_and_immune_review",
      gene_set = label,
      audited_broad_cell_type = ct,
      n_real_genes = length(real_genes),
      real_score = as.numeric(real_scores[ct]),
      random_mean = mean(rs, na.rm = TRUE),
      random_sd = sd(rs, na.rm = TRUE),
      random_p95 = as.numeric(quantile(rs, 0.95, na.rm = TRUE)),
      benchmark_percentile = mean(rs <= as.numeric(real_scores[ct]), na.rm = TRUE)
    )
  }))
}

sensitivity <- rbindlist(lapply(names(gene_sets), function(nm) {
  benchmark_one(gene_sets[[nm]], nm, n_iter = 1000)
}))
sensitivity <- merge(
  sensitivity,
  sensitivity_module,
  by = "audited_broad_cell_type",
  all.x = TRUE
)
sensitivity[, interpretation := fifelse(
  audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95,
  "TAL signal persists after excluding low-confidence and immune-review clusters",
  fifelse(
    benchmark_percentile >= 0.95,
    "high benchmark percentile in conservative analysis",
    fifelse(benchmark_percentile >= 0.90, "moderate conservative signal", "not supported in conservative analysis")
  )
)]
setorder(sensitivity, gene_set, -benchmark_percentile)
fwrite(
  sensitivity,
  file.path(table_dir, "audited_locus_scrna_sensitivity_summary.tsv"),
  sep = "\t"
)

pdf(file.path(figure_dir, "audited_locus_scrna_sensitivity_benchmark.pdf"), width = 9, height = 4.8)
print(
  ggplot(sensitivity, aes(
    x = reorder(audited_broad_cell_type, benchmark_percentile),
    y = benchmark_percentile,
    fill = gene_set
  )) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(x = NULL, y = "Conservative benchmark percentile", fill = "Gene set") +
    theme_bw(base_size = 10)
)
dev.off()

cat("wrote\t", file.path(table_dir, "audited_locus_celltype_evidence_summary.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "audited_locus_scrna_sensitivity_summary.tsv"), "\n", sep = "")
cat("wrote\t", file.path(figure_dir, "audited_locus_scrna_sensitivity_benchmark.pdf"), "\n", sep = "")
