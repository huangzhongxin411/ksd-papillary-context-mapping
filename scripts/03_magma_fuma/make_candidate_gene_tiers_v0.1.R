suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

magma <- fread("results/tables/magma_genes.tsv")
candidate <- fread("results/tables/phase1_candidate_genes.tsv")
overlap <- fread("results/tables/magma_vs_locus_overlap.tsv")
drivers <- fread("results/tables/magma_top50_tal_driver_genes.tsv")
evidence <- fread("results/tables/magma_scrna_evidence_summary.tsv")
cfg <- fread("config/gse231569_cluster_to_celltype_audited.tsv")

gene_sets <- list(
  magma_top50 = scan("results/gene_sets/magma_top50.txt", what = character(), quiet = TRUE),
  magma_top100 = scan("results/gene_sets/magma_top100.txt", what = character(), quiet = TRUE),
  magma_top200 = scan("results/gene_sets/magma_top200.txt", what = character(), quiet = TRUE),
  magma_suggestive_p1e4 = scan("results/gene_sets/magma_suggestive_p1e4.txt", what = character(), quiet = TRUE),
  magma_fdr05 = scan("results/gene_sets/magma_fdr05.txt", what = character(), quiet = TRUE)
)

genes <- sort(unique(c(gene_sets$magma_top200, overlap[is_tal_driver_gene == TRUE, gene], "UMOD", "HIBADH", "FAM13A")))
out <- data.table(gene = genes)

cand_map <- candidate[
  ,
  .(
    locus_id = paste(sort(unique(locus_id)), collapse = ";"),
    lead_snp = paste(sort(unique(nearest_lead_snp)), collapse = ";"),
    lead_snp_p = min(lead_snp_p, na.rm = TRUE),
    mapping_method = paste(sort(unique(mapping_method)), collapse = ";")
  ),
  by = .(gene)
]
out <- merge(out, cand_map, by = "gene", all.x = TRUE)
out <- merge(
  out,
  magma[, .(gene = gene_symbol, magma_rank = rank, magma_p = p, magma_fdr = fdr, magma_bonferroni = bonferroni_significant, magma_suggestive = suggestive)],
  by = "gene",
  all.x = TRUE
)
out <- merge(
  out,
  unique(drivers[, .(gene, locus_driver_gene = driver_class %in% c("candidate_TAL_driver", "supporting_TAL_expressed_gene"), magma_tal_driver_class = driver_class)]),
  by = "gene",
  all.x = TRUE
)
out[is.na(locus_driver_gene), locus_driver_gene := FALSE]

for (nm in names(gene_sets)) {
  out[, (paste0("in_", nm)) := gene %in% gene_sets[[nm]]]
}
out[, gene_set_basis := fifelse(
  in_magma_top50, "magma_top50",
  fifelse(in_magma_top100, "magma_top100",
    fifelse(in_magma_top200, "magma_top200",
      fifelse(in_magma_suggestive_p1e4, "magma_suggestive_p1e4",
        fifelse(in_magma_fdr05, "magma_fdr05", "locus_only")
      )
    )
  )
)]

obj <- readRDS("data/processed/GSE231569/phase1_gse231569_quick_seurat.rds")
cfg[, cluster := as.character(cluster)]
obj$audited_broad_cell_type <- unname(setNames(cfg$broad_cell_type, cfg$cluster)[as.character(obj$seurat_clusters)])
obj$audited_annotation_confidence <- unname(setNames(cfg$annotation_confidence, cfg$cluster)[as.character(obj$seurat_clusters)])
mat <- GetAssayData(obj, assay = "RNA", layer = "data")

gene_expr_summary <- rbindlist(lapply(out$gene, function(g) {
  if (!g %in% rownames(mat)) {
    return(data.table(gene = g, scrna_top_celltype = NA_character_, scrna_avg_expression = NA_real_, scrna_detection_rate = NA_real_))
  }
  rows <- rbindlist(lapply(sort(unique(obj$audited_broad_cell_type)), function(ct) {
    cells <- colnames(obj)[obj$audited_broad_cell_type == ct]
    vals <- mat[g, cells]
    data.table(
      gene = g,
      scrna_top_celltype = ct,
      scrna_avg_expression = mean(vals),
      scrna_detection_rate = mean(vals > 0)
    )
  }))
  rows[order(-scrna_avg_expression, -scrna_detection_rate)][1]
}))
out <- merge(out, gene_expr_summary, by = "gene", all.x = TRUE)

ev_map <- evidence[
  ,
  .(
    gene_set_basis = gene_set,
    scrna_top_celltype = audited_broad_cell_type,
    scrna_benchmark_percentile = benchmark_percentile,
    scrna_evidence_class = evidence_class,
    annotation_confidence
  )
]
out <- merge(out, ev_map, by = c("gene_set_basis", "scrna_top_celltype"), all.x = TRUE)

out[, evidence_level_pre_twas := fifelse(
  !is.na(magma_rank) & magma_rank <= 200 & scrna_top_celltype == "Loop_of_Henle_TAL" & scrna_benchmark_percentile >= 0.90,
  "MAGMA_supported_TAL_candidate",
  fifelse(
    !is.na(magma_rank) & magma_rank <= 200,
    "MAGMA_supported_non_TAL_candidate",
    fifelse(locus_driver_gene == TRUE & scrna_top_celltype == "Loop_of_Henle_TAL", "locus_only_TAL_driver", "exploratory_context_gene")
  )
)]
out[, tier := fifelse(
  evidence_level_pre_twas == "MAGMA_supported_TAL_candidate",
  "pre_TWAS_Tier2_candidate",
  fifelse(
    evidence_level_pre_twas == "MAGMA_supported_non_TAL_candidate",
    "pre_TWAS_Tier3_context",
    "pre_TWAS_exploratory"
  )
)]
out[, interpretation := fifelse(
  evidence_level_pre_twas == "MAGMA_supported_TAL_candidate",
  "MAGMA-prioritized gene with audited TAL localization support; requires TWAS/SMR-coloc before Tier 1.",
  fifelse(
    evidence_level_pre_twas == "MAGMA_supported_non_TAL_candidate",
    "MAGMA-prioritized gene with non-TAL or heterogeneous single-cell localization; retain for TWAS/SMR-coloc.",
    fifelse(
      evidence_level_pre_twas == "locus_only_TAL_driver",
      "Locus-based TAL driver without MAGMA top200 support; supplementary exploratory signal.",
      "Exploratory context gene pending additional genetic prioritization."
    )
  )
)]

setcolorder(out, c(
  "gene",
  "locus_id",
  "lead_snp",
  "lead_snp_p",
  "mapping_method",
  "locus_driver_gene",
  "magma_rank",
  "magma_p",
  "magma_fdr",
  "magma_bonferroni",
  "gene_set_basis",
  "scrna_top_celltype",
  "scrna_benchmark_percentile",
  "scrna_detection_rate",
  "annotation_confidence",
  "evidence_level_pre_twas",
  "tier",
  "interpretation"
))
setorder(out, magma_rank, gene)

fwrite(out, "results/tables/candidate_gene_tiers_v0.1.tsv", sep = "\t")
cat("wrote\tresults/tables/candidate_gene_tiers_v0.1.tsv\n")
