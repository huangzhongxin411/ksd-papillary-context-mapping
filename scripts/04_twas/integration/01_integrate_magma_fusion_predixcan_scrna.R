suppressPackageStartupMessages(library(data.table))

read_or_empty <- function(path, cols) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(as.data.table(setNames(replicate(length(cols), logical(), simplify = FALSE), cols)))
  }
  fread(path, fill = TRUE)
}

candidate <- fread("results/tables/candidate_gene_tiers_v0.1.tsv")

fusion <- read_or_empty(
  "results/tables/twas_fusion_best_per_gene.tsv",
  c("gene", "tissue", "fusion_p", "fusion_fdr")
)
if (nrow(fusion)) {
  fusion_best <- fusion[, .SD[1], by = gene]
  setnames(fusion_best, "tissue", "fusion_best_tissue", skip_absent = TRUE)
} else {
  fusion_best <- data.table(gene = character(), fusion_best_tissue = character(), fusion_p = numeric(), fusion_fdr = numeric())
}
setnames(fusion_best, c("fusion_p", "fusion_fdr"), c("fusion_best_p", "fusion_best_fdr"), skip_absent = TRUE)

predixcan <- read_or_empty(
  "results/tables/spredixcan_best_per_gene.tsv",
  c("gene", "tissue", "predixcan_p", "predixcan_fdr")
)
if (nrow(predixcan)) {
  predixcan_best <- predixcan[, .SD[1], by = gene]
  setnames(predixcan_best, "tissue", "predixcan_best_tissue", skip_absent = TRUE)
} else {
  predixcan_best <- data.table(gene = character(), predixcan_best_tissue = character(), predixcan_p = numeric(), predixcan_fdr = numeric())
}
setnames(predixcan_best, c("predixcan_p", "predixcan_fdr"), c("predixcan_best_p", "predixcan_best_fdr"), skip_absent = TRUE)

smulti <- read_or_empty(
  "results/tables/smultixcan_results.tsv",
  c("gene", "smultixcan_p", "smultixcan_fdr")
)

dt <- merge(candidate, fusion_best[, .(gene, fusion_best_tissue, fusion_best_p, fusion_best_fdr)], by = "gene", all.x = TRUE)
dt <- merge(dt, predixcan_best[, .(gene, predixcan_best_tissue, predixcan_best_p, predixcan_best_fdr)], by = "gene", all.x = TRUE)
if (nrow(smulti)) {
  dt <- merge(dt, smulti[, .(gene, smultixcan_p, smultixcan_fdr)], by = "gene", all.x = TRUE)
} else {
  dt[, `:=`(smultixcan_p = NA_real_, smultixcan_fdr = NA_real_)]
}

dt[, fusion_support := fifelse(!is.na(fusion_best_fdr) & fusion_best_fdr < 0.05, "strong",
  fifelse(!is.na(fusion_best_p) & fusion_best_p < 1e-4, "suggestive",
    fifelse(!is.na(fusion_best_p), "not_significant", "not_available")
  )
)]
dt[, predixcan_support := fifelse(!is.na(predixcan_best_fdr) & predixcan_best_fdr < 0.05, "strong",
  fifelse(!is.na(predixcan_best_p) & predixcan_best_p < 1e-4, "suggestive",
    fifelse(!is.na(predixcan_best_p), "not_significant", "not_available")
  )
)]
dt[, magma_support := fifelse(!is.na(magma_fdr) & magma_fdr < 0.05, "fdr_significant",
  fifelse(!is.na(magma_rank) & magma_rank <= 200, "top200", "not_significant_or_context")
)]

dt[, has_twas_support := fusion_support %in% c("strong", "suggestive") | predixcan_support %in% c("strong", "suggestive")]
dt[, has_tal_support := scrna_top_celltype == "Loop_of_Henle_TAL" & scrna_evidence_class == "strong_primary"]
dt[, twas_status := fifelse(
  fusion_support != "not_available" & predixcan_support != "not_available",
  "fusion_predixcan_completed",
  fifelse(
    fusion_support != "not_available",
    "fusion_only",
    fifelse(predixcan_support != "not_available", "predixcan_only", "not_run_resource_missing")
  )
)]

dt[, current_tier := fifelse(
  magma_support %in% c("fdr_significant", "top200") & has_twas_support & has_tal_support,
  "Tier1_candidate_pending_coloc",
  fifelse(
    magma_support %in% c("fdr_significant", "top200") & has_tal_support,
    "Tier2_moderate_MAGMA_scRNA",
    fifelse(
      has_twas_support | magma_support %in% c("fdr_significant", "top200"),
      "Tier3_context",
      "Exploratory"
    )
  )
)]

dt[, priority_reason := fifelse(
  current_tier == "Tier1_candidate_pending_coloc",
  "MAGMA plus TWAS plus TAL scRNA support; requires SMR/coloc before high-confidence causal interpretation.",
  fifelse(
    current_tier == "Tier2_moderate_MAGMA_scRNA",
    "MAGMA plus TAL scRNA support; TWAS unavailable or not yet supportive.",
    fifelse(current_tier == "Tier3_context", "Partial MAGMA/TWAS/scRNA support.", "Exploratory or context-only evidence.")
  )
)]

priority_genes <- unique(c(
  "UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2", "FAM13A",
  dt[!is.na(magma_rank) & magma_rank <= 50, gene],
  dt[has_twas_support == TRUE, gene],
  dt[current_tier == "Tier1_candidate_pending_coloc", gene]
))

priority <- dt[gene %in% priority_genes, .(
  gene,
  locus_id,
  lead_snp,
  lead_snp_p,
  magma_rank,
  magma_p,
  fusion_best_p,
  predixcan_best_p,
  scrna_top_celltype,
  priority_reason,
  recommended_for_smr = TRUE,
  recommended_for_coloc = TRUE
)]
priority <- unique(priority, by = "gene")
setorder(priority, magma_rank, gene)

cross <- dt[, .(
  gene,
  locus_id,
  lead_snp,
  lead_snp_p,
  magma_rank,
  magma_p,
  magma_fdr,
  magma_support,
  fusion_best_tissue,
  fusion_best_p,
  fusion_best_fdr,
  fusion_support,
  predixcan_best_tissue,
  predixcan_best_p,
  predixcan_best_fdr,
  predixcan_support,
  smultixcan_p,
  smultixcan_fdr,
  twas_status,
  scrna_top_celltype,
  scrna_benchmark_percentile,
  scrna_evidence_class,
  current_tier,
  priority_reason
)]

fwrite(cross, "results/tables/twas_cross_method_overlap.tsv", sep = "\t")
fwrite(cross[order(current_tier, magma_rank)], "results/tables/twas_cross_method_priority_summary.tsv", sep = "\t")
fwrite(cross, "results/tables/candidate_gene_tiers_v0.2.tsv", sep = "\t")
fwrite(priority, "results/tables/priority_loci_for_smr_coloc.tsv", sep = "\t")

message("wrote\tresults/tables/twas_cross_method_overlap.tsv")
message("wrote\tresults/tables/candidate_gene_tiers_v0.2.tsv")
message("wrote\tresults/tables/priority_loci_for_smr_coloc.tsv")
