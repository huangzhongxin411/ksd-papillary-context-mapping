suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final_main_figures_v0.2", recursive = TRUE, showWarnings = FALSE)

magma <- fread("results/tables/magma_genes.tsv")
p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")
plaque <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
coupling <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")

genes <- unique(c(magma[rank <= 100, gene_symbol], p1$gene, infl[contribution_rank <= 30, gene]))
tiers <- data.table(gene = genes)
tiers <- merge(tiers, magma[, .(gene = gene_symbol, MAGMA_P = p, MAGMA_FDR = fdr, MAGMA_rank = rank)], by = "gene", all.x = TRUE)
tiers <- merge(tiers, p1[, .(gene, snRNA_top_celltype = scrna_top_celltype,
                             Loop_TAL_score = TAL_avg_expression,
                             TAL_specificity = specificity_ratio_avg,
                             interpretation_role = manuscript_role)], by = "gene", all.x = TRUE)
tiers <- merge(tiers, infl[, .(gene, contribution_score, contribution_rank, candidate_role)], by = "gene", all.x = TRUE)
tiers <- merge(tiers, plaque[, .(gene, plaque_response_FDR = fdr, plaque_response_P = p_value, paired_delta)], by = "gene", all.x = TRUE)

inj <- coupling[analysis == "Paired delta" & injury_module == "injury_remodeling",
                .(module_name, injury_coupling = rho, injury_coupling_FDR = fdr, robustness_summary)]
module_coupling <- inj[module_name == "MAGMA_top50"]
tiers[, injury_coupling := module_coupling$injury_coupling[1]]
tiers[, injury_coupling_FDR := module_coupling$injury_coupling_FDR[1]]

tiers[, `:=`(
  locus_id = NA_character_,
  GWAS_nearest_or_mapped = "MAGMA_prioritized",
  TWAS_min_P = NA_real_,
  TWAS_FDR = NA_real_,
  TWAS_tissues = NA_character_,
  SMR_P = NA_real_,
  HEIDI_P = NA_real_,
  coloc_PP4 = NA_real_,
  spatial_context_support = "pending_resource_landing"
)]
idx <- is.na(tiers$snRNA_top_celltype) & tiers$gene %in% infl$gene
tiers[idx, snRNA_top_celltype := infl[match(tiers[idx, gene], gene), top_celltype]]
tiers[, snRNA_support := !is.na(snRNA_top_celltype) & (snRNA_top_celltype == "Loop_of_Henle_TAL" | !is.na(contribution_score))]
tiers[, plaque_support := (!is.na(plaque_response_FDR) & plaque_response_FDR < 0.10) | (!is.na(injury_coupling_FDR) & injury_coupling_FDR < 0.10)]
tiers[, gwas_magma_support := !is.na(MAGMA_FDR) & MAGMA_FDR < 0.05]
tiers[, final_tier := fcase(
  gwas_magma_support & snRNA_support & plaque_support & gene %in% p1$gene, "Tier 1 current integrated candidate",
  gwas_magma_support & (snRNA_support | plaque_support), "Tier 2 MAGMA-context candidate",
  !is.na(contribution_score) | plaque_support, "Tier 3 functional/context candidate",
  default = "Context-only or pending"
)]
tiers[, claim_boundary := "Prioritization/context tier only; not causal validation, TWAS convergence, colocalization or spatial validation"]
setorder(tiers, final_tier, MAGMA_rank)

out_cols <- c("gene", "locus_id", "GWAS_nearest_or_mapped", "MAGMA_P", "MAGMA_FDR",
              "TWAS_min_P", "TWAS_FDR", "TWAS_tissues", "SMR_P", "HEIDI_P", "coloc_PP4",
              "snRNA_top_celltype", "Loop_TAL_score", "TAL_specificity", "spatial_context_support",
              "plaque_response_FDR", "paired_delta", "injury_coupling", "injury_coupling_FDR",
              "final_tier", "interpretation_role", "claim_boundary")
fwrite(tiers[, ..out_cols], "results/tables/candidate_gene_tiers_v1.0.tsv", sep = "\t")

supp_tables <- data.table(
  table_id = paste0("Table S", 1:12),
  file = c(
    "results/phase1_v0.1/tables/phase1_gwas_qc_report.tsv",
    "results/phase1_v0.1/tables/phase1_2025_lead_snps.tsv",
    "results/tables/magma_genes.tsv",
    "results/tables/gse231569_celllevel_magma_scores.tsv",
    "results/tables/gse231569_donor_celltype_magma_score_stats.tsv",
    "results/tables/loop_tal_influential_magma_genes.tsv",
    "results/gse73680/tables/gse73680_patient_level_module_response.tsv",
    "results/tables/gse73680_risk_injury_correlation_robustness.tsv",
    "results/tables/magma_pathway_enrichment.tsv",
    "results/tables/go_bp_redundancy_reduced_terms.tsv",
    "results/twas/twas_resource_status.tsv",
    "results/tables/candidate_gene_tiers_v1.0.tsv"
  ),
  purpose = c(
    "GWAS QC", "Lead loci", "MAGMA gene-level results", "Cell-level snRNA module scores",
    "Donor-level Loop/TAL score support", "Loop/TAL contributing genes",
    "GSE73680 module response", "Risk-injury coupling robustness",
    "GO enrichment", "GO redundancy audit", "TWAS resource status", "Integrated candidate tiers"
  )
)
fwrite(supp_tables, "docs/supplementary_table_plan_v0.2.tsv", sep = "\t")

supp_figs <- data.table(
  figure_id = paste0("Figure S", 1:8),
  file = c(
    "results/figures/figure2_umap_audited_celltypes_v0.1.pdf",
    "results/figures/figure2_magma_score_umap_v0.1.pdf",
    "results/figures/figure2_donor_level_magma_score_boxplot.pdf",
    "results/figures/figure2_loop_tal_gene_contribution_heatmap.pdf",
    "results/figures/figure5_risk_injury_correlation_robustness.pdf",
    "results/figures/figure5_go_redundancy_reduced_dotplot.pdf",
    "results/figures/gse73680_magma_vs_injury_scatter.pdf",
    "results/figures/figure4_paired_module_spaghetti_v0.1.pdf"
  ),
  purpose = c(
    "Audited UMAP source", "MAGMA score UMAP source", "Donor-level snRNA score support",
    "Contribution-score components", "Risk-injury robustness", "GO audit", "MAGMA-injury scatter",
    "Standalone paired spaghetti"
  )
)
fwrite(supp_figs, "docs/supplementary_figure_plan_v0.2.tsv", sep = "\t")

boundary <- data.table(
  claim = c("Loop/TAL-associated cellular context", "MAGMA module-level disease-context association",
            "Functional context and injury coupling", "TWAS/SMR/spatial extensions"),
  supported_language = c("supports/maps to/converges on", "supports module-level association",
                         "functional interpretation and disease-context coupling", "resource landing / pending evidence"),
  forbidden_language = c("causal cell type; causal mediation", "disease-gene validation; causal injury mechanism",
                         "mechanistic validation; therapeutic target", "causal validation; spatial validation; TWAS convergence"),
  status = c("allowed", "allowed", "allowed_with_boundary", "not_claimed")
)
fwrite(boundary, "docs/claim_boundary_audit_v0.2.tsv", sep = "\t")

figs <- c(
  "figure1_integrative_framework_v0.4",
  "figure2_magma_scrna_localization_v0.5",
  "figure3_p1_gene_evidence_v0.7",
  "figure4_gse73680_disease_context_v1.0",
  "figure5_functional_context_v0.2"
)
for (stem in figs) {
  for (ext in c("pdf", "png")) {
    src <- file.path("results/figures", paste0(stem, ".", ext))
    if (file.exists(src)) file.copy(src, file.path("results/figures/final_main_figures_v0.2", basename(src)), overwrite = TRUE)
  }
}

writeLines(c(
  "# Supervisor Review Cover Memo v0.2",
  "",
  "This package contains the Phase 10 converged manuscript and five main candidate figures.",
  "",
  "Core claim: MAGMA-prioritized KSD genes support a Loop/TAL-associated renal papillary single-nucleus context and MAGMA module-level plaque/stone papilla disease-context association.",
  "",
  "New anti-reviewer hardening includes donor-level snRNA score summaries, paired-delta and residual risk-injury coupling robustness, GO background/redundancy audit and candidate_gene_tiers_v1.0.",
  "",
  "TWAS, SMR/coloc and spatial transcriptomics remain resource-landing extensions and are not used as manuscript claims in v0.9."
), "docs/supervisor_review_cover_memo_v0.2.md", useBytes = TRUE)

message("wrote candidate tiers and supervisor package plans")
