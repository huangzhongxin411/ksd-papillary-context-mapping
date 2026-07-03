suppressPackageStartupMessages(library(data.table))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

integrated <- data.table(
  evidence_layer = c(
    "GWAS/MAGMA",
    "MAGMA to GSE231569 scRNA localization",
    "P1 six-gene evidence spectrum",
    "GSE73680 disease-context analysis",
    "TWAS/coloc/spatial extensions"
  ),
  analysis = c(
    "Locked KSD GWAS summary statistics followed by MAGMA gene-level prioritization.",
    "MAGMA-prioritized gene sets projected onto audited GSE231569 renal papillary single-nucleus cell types.",
    "Six P1 candidate genes evaluated across TAL expression, donor detection, specificity and gene-role assignment.",
    "Patient-aware GSE73680 plaque/stone papilla versus control/adjacent disease-context expression analysis.",
    "TWAS, SMR/coloc and spatial transcriptomics preparation/resource checks."
  ),
  main_result = c(
    "MAGMA prioritized KSD-associated genes from the locked GWAS reconstruction.",
    "MAGMA-prioritized gene sets converged on a Loop of Henle/TAL-associated cellular context.",
    "P1 genes formed an interpretable TAL/transport/calcium-handling expression spectrum rather than a uniform TAL marker set.",
    "GSE73680 supported MAGMA module-level disease-context expression association, not uniform P1 single-gene differential expression.",
    "Required external TWAS weights, eQTL/coloc resources and spatial matrix/image resources were not fully available."
  ),
  support_strength = c("strong", "strong", "moderate", "moderate", "resource_limited"),
  claim_supported = c(
    "KSD GWAS signals can be summarized as MAGMA-prioritized candidate genes and modules.",
    "MAGMA-prioritized KSD genes localize to a TAL-associated single-nucleus context.",
    "P1 genes support a MAGMA + scRNA-based TAL/transport/calcium expression spectrum.",
    "GSE73680 supports disease-context expression association for MAGMA-prioritized modules.",
    "These analyses define pending extensions and resource gaps."
  ),
  claim_not_supported = c(
    "Causal mediation; gene-level mechanism proof.",
    "Causal cell type mediation; spatial validation; TWAS or coloc convergence.",
    "Uniform TAL specificity for all P1 genes; causal P1 validation.",
    "Uniform P1 single-gene differential expression; causal plaque formation; TAL mechanism proof.",
    "Completed TWAS support; completed colocalization; completed spatial validation."
  ),
  figure = c("Figure 1/2", "Figure 2", "Figure 3", "Figure 4", "Limitations/Pending extensions"),
  table = c(
    "results/tables/magma_qc_summary.tsv; results/tables/magma_genes.tsv",
    "results/tables/magma_gene_set_summary.tsv; results/tables/magma_scrna_evidence_summary.tsv",
    "results/tables/p1_tal_gene_interpretation_summary.tsv",
    "results/gse73680/tables/gse73680_disease_context_summary_v0.2.tsv",
    "docs/twas_resource_acquisition_log.md; docs/coloc_run_plan.md; docs/phase4_spatial_prep_status.md"
  ),
  notes = c(
    "Main genetic prioritization layer.",
    "Single-nucleus localization layer; not causal mediation.",
    "Gene-centric interpretation layer for manuscript Result 3.",
    "Disease-context expression layer; module-level support only.",
    "Do not treat unavailable-resource extensions as negative results."
  )
)
fwrite(integrated, "results/tables/integrated_evidence_summary_v0.1.tsv", sep = "\t")

p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
single <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
p1 <- merge(p1, single[, .(gene, GSE73680_single_gene_response = interpretation,
                           GSE73680_single_gene_p = p_value,
                           GSE73680_single_gene_fdr = fdr)], by = "gene", all.x = TRUE)
p1[, GSE73680_module_context := "MAGMA module-level disease-context support; not gene-specific validation"]
p1[, TWAS_status := "resource_limited_not_completed"]
p1[, coloc_status := "resource_limited_not_completed"]
p1[, MAGMA_support := paste0("MAGMA rank ", magma_rank, "; P=", magma_p)]
p1[, scRNA_TAL_context := fcase(
  overall_evidence_class == "P1_strong_TAL_context", "strong TAL-associated scRNA context",
  overall_evidence_class == "P1_moderate_TAL_context", "moderate/broad TAL-associated scRNA context",
  default = "contextual scRNA support"
)]
p1[, P1_role := manuscript_role]
p1[, final_tier_v0.3 := fcase(
  gene %in% c("UMOD", "CLDN10", "CLDN14", "CASR"), "Tier1B_MAGMA_scRNA_context_candidate",
  gene == "HIBADH", "Tier2_MAGMA_scRNA_supporting_context_gene",
  gene == "PKD2", "Tier2_broad_epithelial_context_nominal_GSE73680",
  default = "Tier2_MAGMA_scRNA_candidate"
)]
p1[, interpretation := fcase(
  final_tier_v0.3 == "Tier1B_MAGMA_scRNA_context_candidate",
  "Prioritized as a MAGMA + scRNA-supported contextual candidate; GSE73680 supports the broader MAGMA module but not this gene as individually disease-validated.",
  gene == "HIBADH",
  "MAGMA + scRNA-supported contextual candidate with supporting rather than leading biological interpretation.",
  gene == "PKD2",
  "Broad renal epithelial context gene with nominal GSE73680 paired response; not FDR-supported as a single gene.",
  default = "MAGMA + scRNA-supported candidate."
)]
p1[, claim_boundary := "Do not describe as disease-validated, causal, TWAS-supported or colocalized; GSE73680 support is module-level unless single-gene FDR support is shown."]
tiers <- p1[, .(gene, current_tier, MAGMA_support, scRNA_TAL_context, P1_role,
                GSE73680_single_gene_response, GSE73680_module_context,
                TWAS_status, coloc_status, final_tier_v0.3, interpretation, claim_boundary)]
tiers[, gene_order := match(gene, c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"))]
setorder(tiers, gene_order)
tiers[, gene_order := NULL]
fwrite(tiers, "results/tables/candidate_gene_tiers_v0.3.tsv", sep = "\t")

message("wrote Phase 6 integration tables")
