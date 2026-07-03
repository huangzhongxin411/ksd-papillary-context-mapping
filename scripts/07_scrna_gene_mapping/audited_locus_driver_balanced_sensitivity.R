suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

set.seed(20260617)

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
cfg_path <- "config/gse231569_cluster_to_celltype_audited.tsv"
ranked_path <- "results/tables/phase1_candidate_genes_ranked_protein_coding.tsv"
gene_set_dir <- "results/gene_sets"
table_dir <- "results/tables"
figure_dir <- "results/figures"
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ranked <- fread(ranked_path)
top50 <- scan(file.path(gene_set_dir, "locus_top50.txt"), what = character(), quiet = TRUE)
top50_meta <- ranked[gene %in% top50]
setorder(top50_meta, rank_in_candidate_list)

obj <- readRDS(obj_path)
cfg <- fread(cfg_path)
cfg[, cluster := as.character(cluster)]
cluster_to_broad <- setNames(cfg$broad_cell_type, cfg$cluster)
cluster_to_conf <- setNames(cfg$annotation_confidence, cfg$cluster)
cluster_to_immune <- setNames(cfg$immune_review_flag, cfg$cluster)
cluster_to_mural <- setNames(cfg$mural_review_flag, cfg$cluster)
obj$audited_broad_cell_type <- unname(cluster_to_broad[as.character(obj$seurat_clusters)])
obj$audited_annotation_confidence <- unname(cluster_to_conf[as.character(obj$seurat_clusters)])
obj$audited_immune_review_flag <- unname(cluster_to_immune[as.character(obj$seurat_clusters)])
obj$audited_mural_review_flag <- unname(cluster_to_mural[as.character(obj$seurat_clusters)])
if (!"donor_id" %in% colnames(obj@meta.data)) obj$donor_id <- obj$sample_id

make_context <- function(obj_in, analysis_version) {
  mat <- GetAssayData(obj_in, assay = "RNA", layer = "data")
  cell_type <- factor(obj_in$audited_broad_cell_type)
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

  list(
    obj = obj_in,
    mat = mat,
    gene_by_group = gene_by_group,
    groups = groups,
    genes_universe = genes_universe,
    gene_bins = gene_bins,
    bin_lookup = bin_lookup,
    analysis_version = analysis_version
  )
}

sample_matched <- function(real_genes, ctx) {
  real_genes <- intersect(real_genes, names(ctx$bin_lookup))
  sampled <- character()
  for (g in real_genes) {
    pool <- setdiff(ctx$gene_bins[[ctx$bin_lookup[[g]]]], real_genes)
    if (length(pool) == 0) pool <- setdiff(ctx$genes_universe, real_genes)
    sampled <- c(sampled, sample(pool, 1))
  }
  unique(sampled)
}

benchmark_one <- function(real_genes, label, ctx, n_iter = 1000) {
  real_genes <- intersect(unique(real_genes), rownames(ctx$gene_by_group))
  real_scores <- Matrix::colMeans(ctx$gene_by_group[real_genes, , drop = FALSE])
  rand_scores <- matrix(NA_real_, nrow = n_iter, ncol = length(ctx$groups), dimnames = list(NULL, ctx$groups))
  for (i in seq_len(n_iter)) {
    rg <- intersect(sample_matched(real_genes, ctx), rownames(ctx$gene_by_group))
    rand_scores[i, ] <- Matrix::colMeans(ctx$gene_by_group[rg, , drop = FALSE])
  }
  rbindlist(lapply(ctx$groups, function(ct) {
    rs <- rand_scores[, ct]
    data.table(
      analysis_version = ctx$analysis_version,
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

full_ctx <- make_context(obj, "full_audited")
keep_cells <- colnames(obj)[
  obj$audited_annotation_confidence != "low_or_exploratory" &
    !(obj$audited_immune_review_flag %in% TRUE)
]
cons_ctx <- make_context(subset(obj, cells = keep_cells), "conservative_exclude_low_or_exploratory_and_immune_review")

top50_present <- intersect(top50_meta$gene, rownames(full_ctx$gene_by_group))
tal_cell_type <- "Loop_of_Henle_TAL"
top50_tal_score <- mean(full_ctx$gene_by_group[top50_present, tal_cell_type])

driver_rows <- lapply(top50_meta$gene, function(g) {
  if (!g %in% rownames(full_ctx$gene_by_group)) {
    return(data.table(gene = g, present_in_seurat = FALSE))
  }
  vals <- as.numeric(full_ctx$gene_by_group[g, ])
  names(vals) <- colnames(full_ctx$gene_by_group)
  tal_avg <- vals[tal_cell_type]
  top_ct <- names(vals)[which.max(vals)]
  specificity_rank <- rank(-vals, ties.method = "min")[tal_cell_type]
  cells_tal <- colnames(full_ctx$obj)[full_ctx$obj$audited_broad_cell_type == tal_cell_type]
  gene_vec <- full_ctx$mat[g, cells_tal]
  without <- setdiff(top50_present, g)
  without_score <- mean(full_ctx$gene_by_group[without, tal_cell_type])
  meta <- top50_meta[gene == g][1]
  data.table(
    gene = g,
    present_in_seurat = TRUE,
    locus_id = meta$locus_id,
    rank = meta$rank_in_candidate_list,
    mapping_method = meta$mapping_method,
    evidence_level = meta$evidence_level,
    nearest_lead_snp = meta$nearest_lead_snp,
    lead_snp_p = meta$lead_snp_p,
    distance_to_lead_snp = meta$distance_to_lead_snp,
    TAL_avg_expression = as.numeric(tal_avg),
    TAL_pct_expressed = mean(gene_vec > 0),
    TAL_specificity_rank = as.integer(specificity_rank),
    all_celltype_max = max(vals, na.rm = TRUE),
    top_celltype = top_ct,
    leave_one_gene_TAL_score = without_score,
    leave_one_gene_delta = top50_tal_score - without_score
  )
})
drivers <- rbindlist(driver_rows, fill = TRUE)
drivers[, driver_class := fifelse(
  present_in_seurat == FALSE,
  "not_detected_in_seurat",
  fifelse(
    TAL_specificity_rank == 1 & TAL_pct_expressed >= 0.10 & leave_one_gene_delta >= quantile(leave_one_gene_delta, 0.75, na.rm = TRUE),
    "candidate_TAL_driver",
    fifelse(TAL_specificity_rank <= 2 & TAL_pct_expressed >= 0.05, "supporting_TAL_expressed_gene", "non_TAL_or_low_detection")
  )
)]
setorder(drivers, -leave_one_gene_delta, -TAL_avg_expression)
fwrite(drivers, file.path(table_dir, "audited_locus_top50_tal_driver_genes.tsv"), sep = "\t")

balanced <- ranked[
  gene_type == "protein_coding"
][
  order(rank_in_candidate_list)
][
  ,
  .SD[1],
  by = locus_id
][
  order(rank_in_candidate_list)
][
  1:min(.N, 50)
]
writeLines(balanced$gene, file.path(gene_set_dir, "locus_balanced_top50.txt"))
fwrite(
  balanced[, .(
    gene,
    locus_id,
    rank_in_candidate_list,
    mapping_method,
    evidence_level,
    nearest_lead_snp,
    lead_snp_p,
    distance_to_lead_snp
  )],
  file.path(table_dir, "locus_balanced_top50_genes.tsv"),
  sep = "\t"
)

balanced_bench <- rbind(
  benchmark_one(balanced$gene, "locus_balanced_top50", full_ctx),
  benchmark_one(balanced$gene, "locus_balanced_top50", cons_ctx)
)
fwrite(balanced_bench, file.path(table_dir, "audited_locus_balanced_scrna_benchmark.tsv"), sep = "\t")

top50_loci <- unique(top50_meta$locus_id)
base_tal <- benchmark_one(top50_meta$gene, "locus_top50", cons_ctx)[audited_broad_cell_type == tal_cell_type]
loo <- rbindlist(lapply(top50_loci, function(locus) {
  removed <- top50_meta[locus_id == locus]
  remaining <- top50_meta[locus_id != locus, gene]
  bench <- benchmark_one(remaining, paste0("remove_", locus), cons_ctx)[audited_broad_cell_type == tal_cell_type]
  data.table(
    analysis_version = cons_ctx$analysis_version,
    removed_locus_id = locus,
    removed_genes = paste(removed$gene, collapse = ","),
    n_removed_genes = nrow(removed),
    n_remaining_genes = length(intersect(remaining, rownames(cons_ctx$gene_by_group))),
    TAL_percentile_before_removal = base_tal$benchmark_percentile,
    TAL_percentile_after_removal = bench$benchmark_percentile,
    delta_percentile = bench$benchmark_percentile - base_tal$benchmark_percentile,
    TAL_score_after_removal = bench$real_score
  )
}))
loo[, interpretation := fifelse(
  TAL_percentile_after_removal >= 0.90,
  "TAL_signal_retained_above_0.90",
  fifelse(TAL_percentile_after_removal >= 0.80, "TAL_signal_weakened_but_retained_above_0.80", "TAL_signal_sensitive_to_removed_locus")
)]
setorder(loo, TAL_percentile_after_removal)
fwrite(loo, file.path(table_dir, "audited_locus_leave_one_locus_out.tsv"), sep = "\t")

pdf(file.path(figure_dir, "audited_locus_leave_one_locus_out_tal.pdf"), width = 9, height = 5)
print(
  ggplot(loo, aes(x = reorder(removed_locus_id, TAL_percentile_after_removal), y = TAL_percentile_after_removal)) +
    geom_col(fill = "#4E79A7") +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "#E15759") +
    coord_flip() +
    labs(x = "Removed locus", y = "TAL percentile after removal") +
    theme_bw(base_size = 10)
)
dev.off()

cat("wrote\t", file.path(table_dir, "audited_locus_top50_tal_driver_genes.tsv"), "\n", sep = "")
cat("wrote\t", file.path(gene_set_dir, "locus_balanced_top50.txt"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "audited_locus_balanced_scrna_benchmark.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "audited_locus_leave_one_locus_out.tsv"), "\n", sep = "")
