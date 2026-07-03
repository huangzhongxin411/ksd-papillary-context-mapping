suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

set.seed(20260617)

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
cfg_path <- "config/gse231569_cluster_to_celltype_audited.tsv"
magma_path <- "results/tables/magma_genes.tsv"
candidate_path <- "results/tables/phase1_candidate_genes.tsv"
loci_path <- "results/tables/phase1_2025_loci.tsv"
gene_set_dir <- "results/gene_sets"
table_dir <- "results/tables"
figure_dir <- "results/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

magma <- fread(magma_path)
candidate <- fread(candidate_path)
loci <- fread(loci_path)
top50 <- scan(file.path(gene_set_dir, "magma_top50.txt"), what = character(), quiet = TRUE)

top50_meta <- magma[gene_symbol %in% top50]
setorder(top50_meta, rank)
cand_gene_map <- candidate[
  ,
  .(
    locus_id = paste(sort(unique(locus_id)), collapse = ";"),
    nearest_lead_snp = paste(sort(unique(nearest_lead_snp)), collapse = ";"),
    lead_snp_p = min(lead_snp_p, na.rm = TRUE),
    mapping_method = paste(sort(unique(mapping_method)), collapse = ";"),
    distance_to_lead_snp = min(distance_to_lead_snp, na.rm = TRUE)
  ),
  by = .(gene_symbol = gene)
]
top50_meta <- merge(
  top50_meta,
  cand_gene_map,
  by = "gene_symbol",
  all.x = TRUE
)
top50_meta <- top50_meta[!duplicated(gene_symbol)]

assign_locus_by_position <- function(chr, start, stop) {
  hit <- loci[CHR == chr & start <= end & end >= start]
  if (nrow(hit) == 0) return(NA_character_)
  paste(hit$locus_id, collapse = ";")
}
missing_locus <- is.na(top50_meta$locus_id) | top50_meta$locus_id == ""
if (any(missing_locus)) {
  top50_meta[missing_locus, locus_id := mapply(assign_locus_by_position, chr, start, stop)]
}
top50_meta[is.na(locus_id) | locus_id == "", locus_group := paste0("MAGMA_ONLY_", gene_symbol)]
top50_meta[!is.na(locus_id) & locus_id != "", locus_group := locus_id]

obj <- readRDS(obj_path)
cfg <- fread(cfg_path)
cfg[, cluster := as.character(cluster)]
obj$audited_broad_cell_type <- unname(setNames(cfg$broad_cell_type, cfg$cluster)[as.character(obj$seurat_clusters)])
obj$audited_annotation_confidence <- unname(setNames(cfg$annotation_confidence, cfg$cluster)[as.character(obj$seurat_clusters)])
obj$audited_immune_review_flag <- unname(setNames(cfg$immune_review_flag, cfg$cluster)[as.character(obj$seurat_clusters)])
obj$audited_mural_review_flag <- unname(setNames(cfg$mural_review_flag, cfg$cluster)[as.character(obj$seurat_clusters)])
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
  bins <- cut(
    global_mean[genes_universe],
    breaks = unique(quantile(global_mean[genes_universe], probs = seq(0, 1, length.out = 21), na.rm = TRUE)),
    include.lowest = TRUE
  )
  list(
    obj = obj_in,
    mat = mat,
    gene_by_group = gene_by_group,
    groups = groups,
    genes_universe = genes_universe,
    gene_bins = split(genes_universe, bins),
    bin_lookup = setNames(as.character(bins), genes_universe),
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

tal_cell_type <- "Loop_of_Henle_TAL"
top50_present <- intersect(top50_meta$gene_symbol, rownames(full_ctx$gene_by_group))
top50_tal_score <- mean(full_ctx$gene_by_group[top50_present, tal_cell_type])

driver_rows <- lapply(top50_meta$gene_symbol, function(g) {
  meta <- top50_meta[gene_symbol == g][1]
  if (!g %in% rownames(full_ctx$gene_by_group)) {
    return(data.table(gene = g, present_in_seurat = FALSE, magma_rank = meta$rank, magma_p = meta$p, locus_group = meta$locus_group, locus_id = meta$locus_id))
  }
  vals <- as.numeric(full_ctx$gene_by_group[g, ])
  names(vals) <- colnames(full_ctx$gene_by_group)
  cells_tal <- colnames(full_ctx$obj)[full_ctx$obj$audited_broad_cell_type == tal_cell_type]
  gene_vec <- full_ctx$mat[g, cells_tal]
  without <- setdiff(top50_present, g)
  without_score <- mean(full_ctx$gene_by_group[without, tal_cell_type])
  data.table(
    gene = g,
    present_in_seurat = TRUE,
    locus_group = meta$locus_group,
    locus_id = meta$locus_id,
    magma_rank = meta$rank,
    magma_p = meta$p,
    magma_fdr = meta$fdr,
    nearest_lead_snp = meta$nearest_lead_snp,
    lead_snp_p = meta$lead_snp_p,
    TAL_avg_expression = as.numeric(vals[tal_cell_type]),
    TAL_pct_expressed = mean(gene_vec > 0),
    TAL_specificity_rank = as.integer(rank(-vals, ties.method = "min")[tal_cell_type]),
    all_celltype_max = max(vals, na.rm = TRUE),
    top_celltype = names(vals)[which.max(vals)],
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
fwrite(drivers, file.path(table_dir, "magma_top50_tal_driver_genes.tsv"), sep = "\t")

balanced <- top50_meta[order(rank)][, .SD[1], by = locus_group][order(rank)]
fwrite(
  balanced[, .(gene = gene_symbol, locus_group, locus_id, magma_rank = rank, magma_p = p, magma_fdr = fdr, nearest_lead_snp, lead_snp_p)],
  file.path(table_dir, "magma_locus_balanced_top50_genes.tsv"),
  sep = "\t"
)
balanced_bench <- rbind(
  benchmark_one(balanced$gene_symbol, "magma_locus_balanced_top50", full_ctx),
  benchmark_one(balanced$gene_symbol, "magma_locus_balanced_top50", cons_ctx)
)
fwrite(balanced_bench, file.path(table_dir, "magma_locus_balanced_scrna_benchmark.tsv"), sep = "\t")

base_tal <- benchmark_one(top50_meta$gene_symbol, "magma_top50", cons_ctx)[audited_broad_cell_type == tal_cell_type]
loo <- rbindlist(lapply(unique(top50_meta$locus_group), function(lg) {
  removed <- top50_meta[locus_group == lg]
  remaining <- top50_meta[locus_group != lg, gene_symbol]
  bench <- benchmark_one(remaining, paste0("remove_", lg), cons_ctx)[audited_broad_cell_type == tal_cell_type]
  data.table(
    analysis_version = cons_ctx$analysis_version,
    removed_locus_group = lg,
    removed_locus_id = paste(unique(na.omit(removed$locus_id)), collapse = ";"),
    removed_genes = paste(removed$gene_symbol, collapse = ","),
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
fwrite(loo, file.path(table_dir, "magma_leave_one_locus_out.tsv"), sep = "\t")

pdf(file.path(figure_dir, "magma_leave_one_locus_out_tal.pdf"), width = 9, height = 5)
print(
  ggplot(loo, aes(x = reorder(removed_locus_group, TAL_percentile_after_removal), y = TAL_percentile_after_removal)) +
    geom_col(fill = "#4E79A7") +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "#E15759") +
    coord_flip() +
    labs(x = "Removed locus/group", y = "TAL percentile after removal") +
    theme_bw(base_size = 10)
)
dev.off()

cat("wrote\t", file.path(table_dir, "magma_top50_tal_driver_genes.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "magma_locus_balanced_top50_genes.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "magma_locus_balanced_scrna_benchmark.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "magma_leave_one_locus_out.tsv"), "\n", sep = "")
