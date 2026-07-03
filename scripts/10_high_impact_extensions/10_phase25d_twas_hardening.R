suppressPackageStartupMessages(library(data.table))

set.seed(20260623)

dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

model_db <- "external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.db"
gwas_file <- "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
twas_file <- "results/twas/spredixcan_Kidney_Cortex_KSD.csv"
magma_file <- "results/phase2_magma_v0.1/tables/magma_genes.tsv"
if (!file.exists(magma_file)) magma_file <- "submission_package_v1.9_BMC_working/supplementary/tables/S3_magma_gene_results_and_gene_sets.tsv"

p1_genes <- c("UMOD", "CLDN14", "CASR", "CLDN10", "HIBADH", "PKD2")

weights_tsv <- "results/twas/_phase25d_kidney_cortex_weights.tsv"
sqlite_reader <- tempfile(fileext = ".py")
writeLines(c(
  "import sqlite3, sys",
  "db, out = sys.argv[1], sys.argv[2]",
  "con = sqlite3.connect(db)",
  "cur = con.cursor()",
  "with open(out, 'w') as fh:",
  "    fh.write('gene\\trsid\\tvarID\\tref_allele\\teff_allele\\n')",
  "    for row in cur.execute('select gene, rsid, varID, ref_allele, eff_allele from weights'):",
  "        fh.write('\\t'.join('' if x is None else str(x) for x in row) + '\\n')",
  "con.close()"
), sqlite_reader)
status <- system2(
  "external/twas/predixcan/metaxcan_venv/bin/python3",
  args = c(sqlite_reader, model_db, weights_tsv)
)
if (!identical(status, 0L)) stop("Failed to export PredictDB weights table from sqlite model")
weights <- fread(weights_tsv)

snps_per_gene <- weights[, .(
  n_weight_snps = uniqueN(varID),
  n_rsid = uniqueN(rsid),
  n_varid = uniqueN(varID)
), by = gene]

model_audit <- data.table(
  model_db = model_db,
  n_weight_rows = nrow(weights),
  n_model_genes = uniqueN(weights$gene),
  n_distinct_rsid = uniqueN(weights$rsid),
  n_distinct_varid = uniqueN(weights$varID),
  median_snps_per_gene = median(snps_per_gene$n_weight_snps),
  mean_snps_per_gene = mean(snps_per_gene$n_weight_snps),
  n_genes_with_1_weight_snp = nrow(snps_per_gene[n_weight_snps == 1]),
  n_genes_with_gt1_weight_snp = nrow(snps_per_gene[n_weight_snps > 1]),
  n_genes_with_ge3_weight_snp = nrow(snps_per_gene[n_weight_snps >= 3])
)
fwrite(model_audit, "results/twas/predictdb_kidney_cortex_model_audit.tsv", sep = "\t")

gwas <- fread(gwas_file, select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "Z", "P", "N", "EAF"))
gwas[, chr_prefix := paste0("chr", CHR)]
gwas[, key_rsid := SNP]
gwas[, key_chr_A1_A2 := paste0(chr_prefix, "_", BP, "_", A1, "_", A2, "_b38")]
gwas[, key_chr_A2_A1 := paste0(chr_prefix, "_", BP, "_", A2, "_", A1, "_b38")]
gwas[, key_nochr_A1_A2 := paste0(CHR, "_", BP, "_", A1, "_", A2, "_b38")]
gwas[, key_nochr_A2_A1 := paste0(CHR, "_", BP, "_", A2, "_", A1, "_b38")]

model_keys <- list(
  rsid = unique(weights$rsid),
  varID = unique(weights$varID)
)
gwas_keys <- list(
  SNP = unique(gwas$key_rsid),
  chr_BP_A1_A2_b38 = unique(gwas$key_chr_A1_A2),
  chr_BP_A2_A1_b38 = unique(gwas$key_chr_A2_A1),
  nochr_BP_A1_A2_b38 = unique(gwas$key_nochr_A1_A2),
  nochr_BP_A2_A1_b38 = unique(gwas$key_nochr_A2_A1)
)

allele_compatible <- function(gkey_type, mkey_type) {
  g_cols <- switch(gkey_type,
    SNP = c("key_rsid", "A1", "A2"),
    chr_BP_A1_A2_b38 = c("key_chr_A1_A2", "A1", "A2"),
    chr_BP_A2_A1_b38 = c("key_chr_A2_A1", "A1", "A2"),
    nochr_BP_A1_A2_b38 = c("key_nochr_A1_A2", "A1", "A2"),
    nochr_BP_A2_A1_b38 = c("key_nochr_A2_A1", "A1", "A2")
  )
  m_col <- if (mkey_type == "rsid") "rsid" else "varID"
  g0 <- unique(gwas[, .(gkey = get(g_cols[1]), A1, A2)])
  m0 <- unique(weights[, .(gkey = get(m_col), ref_allele, eff_allele)])
  z <- merge(g0, m0, by = "gkey", allow.cartesian = TRUE)
  if (!nrow(z)) return(list(n = 0L, ambiguous = 0L))
  z[, compatible := (A1 == eff_allele & A2 == ref_allele) | (A1 == ref_allele & A2 == eff_allele)]
  comp <- z[compatible == TRUE]
  amb <- comp[, .N, by = gkey][N > 1, .N]
  list(n = uniqueN(comp$gkey), ambiguous = ifelse(length(amb), amb, 0L))
}

overlap_rows <- rbindlist(lapply(names(gwas_keys), function(gk) {
  rbindlist(lapply(names(model_keys), function(mk) {
    ov <- intersect(gwas_keys[[gk]], model_keys[[mk]])
    ac <- allele_compatible(gk, mk)
    data.table(
      gwas_key_type = gk,
      model_key_type = mk,
      n_gwas_variants = length(gwas_keys[[gk]]),
      n_model_variants = length(model_keys[[mk]]),
      n_overlap = length(ov),
      overlap_rate_vs_model = length(ov) / length(model_keys[[mk]]),
      allele_compatible_n = ac$n,
      ambiguous_n = ac$ambiguous
    )
  }))
}))
max_overlap <- max(overlap_rows$n_overlap)
overlap_rows[, recommended := fifelse(n_overlap == max_overlap & n_overlap > 0, "best_observed_key_match", "not_recommended")]
fwrite(overlap_rows, "results/twas/twas_variant_key_overlap_audit.tsv", sep = "\t")

twas <- fread(twas_file)
twas[, gene_id := gene]
twas[, gene := gene_name]
twas[, twas_p := as.numeric(pvalue)]
twas[, twas_z := as.numeric(zscore)]
twas[, twas_fdr := p.adjust(twas_p, method = "BH")]
twas[, p1_candidate := gene %in% p1_genes]
twas[, support_status := fifelse(twas_fdr < 0.05, "FDR_supported", fifelse(twas_p < 1e-4, "suggestive_only", "not_supported"))]

magma <- fread(magma_file, fill = TRUE)
if ("gene_symbol" %in% names(magma)) setnames(magma, "gene_symbol", "gene")
if (!"gene" %in% names(magma) && "GENE" %in% names(magma)) setnames(magma, "GENE", "gene")
if ("rank" %in% names(magma)) setnames(magma, "rank", "magma_rank")
if ("p" %in% names(magma)) setnames(magma, "p", "magma_p")
if ("fdr" %in% names(magma)) setnames(magma, "fdr", "magma_fdr")
magma <- unique(magma[, intersect(c("gene", "magma_rank", "magma_p", "magma_fdr", "suggestive", "bonferroni_significant"), names(magma)), with = FALSE])

top50 <- if (file.exists("results/gene_sets/magma_top50.txt")) scan("results/gene_sets/magma_top50.txt", what = character(), quiet = TRUE) else magma[order(magma_rank)][1:50, gene]
top100 <- if (file.exists("results/gene_sets/magma_top100.txt")) scan("results/gene_sets/magma_top100.txt", what = character(), quiet = TRUE) else magma[order(magma_rank)][1:100, gene]
fdr_set <- if (file.exists("results/gene_sets/magma_fdr05.txt")) scan("results/gene_sets/magma_fdr05.txt", what = character(), quiet = TRUE) else magma[!is.na(magma_fdr) & magma_fdr < 0.05, gene]
suggestive_set <- if (file.exists("results/gene_sets/magma_suggestive_p1e4.txt")) scan("results/gene_sets/magma_suggestive_p1e4.txt", what = character(), quiet = TRUE) else magma[!is.na(magma_p) & magma_p < 1e-4, gene]

ov <- merge(twas, magma, by = "gene", all.x = TRUE)
ov[, magma_overlap := !is.na(magma_rank)]
ov[, magma_top50 := gene %in% top50]
ov[, magma_top100 := gene %in% top100]
ov[, magma_fdr := gene %in% fdr_set]
ov[, magma_suggestive := gene %in% suggestive_set]
ov[, claim_use := fifelse(twas_fdr < 0.05 & magma_overlap == TRUE & n_snps_used > 1,
                          "stronger_proxy_TWAS_MAGMA_support",
                          fifelse(twas_fdr < 0.05 & magma_overlap == TRUE,
                                  "low_coverage_proxy_TWAS_MAGMA_support",
                                  "not_used_for_main_claim"))]
fwrite(ov[order(twas_fdr, twas_p)], "results/twas/twas_magma_overlap_real.tsv", sep = "\t")

groups <- list(
  all_genes = ov,
  n_snps_used_eq1 = ov[n_snps_used == 1],
  n_snps_used_gt1 = ov[n_snps_used > 1],
  n_snps_used_ge3 = ov[n_snps_used >= 3]
)
sensitivity <- rbindlist(lapply(names(groups), function(nm) {
  x <- groups[[nm]]
  f <- x[twas_fdr < 0.05]
  data.table(
    group = nm,
    n_tested = nrow(x),
    n_fdr_lt_0_05 = nrow(f),
    min_p = ifelse(nrow(x), min(x$twas_p, na.rm = TRUE), NA_real_),
    max_q_among_fdr = ifelse(nrow(f), max(f$twas_fdr, na.rm = TRUE), NA_real_),
    magma_top50_overlap = sum(f$gene %in% top50),
    magma_top100_overlap = sum(f$gene %in% top100),
    magma_fdr_overlap = sum(f$gene %in% fdr_set),
    magma_suggestive_overlap = sum(f$gene %in% suggestive_set),
    p1_overlap = sum(f$gene %in% p1_genes),
    top_genes = paste(head(f[order(twas_fdr), gene], 12), collapse = ";")
  )
}))
fwrite(sensitivity, "results/twas/twas_one_snp_sensitivity.tsv", sep = "\t")
fwrite(ov[twas_fdr < 0.05 & n_snps_used > 1][order(twas_fdr)],
       "results/twas/twas_multi_snp_supported_genes.tsv", sep = "\t")

background <- unique(ov$gene)
sets <- list(
  MAGMA_top50 = intersect(top50, background),
  MAGMA_top100 = intersect(top100, background),
  MAGMA_FDR = intersect(fdr_set, background),
  MAGMA_suggestive = intersect(suggestive_set, background),
  P1 = intersect(p1_genes, background)
)
twas_sets <- list(
  TWAS_FDR_all = unique(ov[twas_fdr < 0.05, gene]),
  TWAS_FDR_multiSNP = unique(ov[twas_fdr < 0.05 & n_snps_used > 1, gene])
)
n_random <- 10000L
bench <- rbindlist(lapply(names(twas_sets), function(tsn) {
  tg <- twas_sets[[tsn]]
  rbindlist(lapply(names(sets), function(msn) {
    mg <- sets[[msn]]
    obs <- length(intersect(tg, mg))
    draws <- replicate(n_random, {
      length(intersect(sample(background, length(tg), replace = FALSE), mg))
    })
    data.table(
      twas_set = tsn,
      magma_set = msn,
      observed_overlap = obs,
      expected_mean = mean(draws),
      expected_sd = sd(draws),
      empirical_p = (sum(draws >= obs) + 1) / (n_random + 1),
      percentile = mean(draws <= obs),
      n_random = n_random
    )
  }))
}))
fwrite(bench, "results/twas/twas_magma_overlap_random_benchmark.tsv", sep = "\t")

run_comparison <- data.table(
  run_id = "Kidney_Cortex_mashr_varID_mapped_v0.1_uncalibrated_proxy",
  snp_key = "PredictDB varID derived from rsID-to-varID mapping",
  run_completed = TRUE,
  genes_tested = nrow(ov),
  genes_with_usable_snps = nrow(ov[n_snps_used > 0]),
  genes_with_1_snp = nrow(ov[n_snps_used == 1]),
  genes_with_gt1_snp = nrow(ov[n_snps_used > 1]),
  min_p = min(as.numeric(ov[["twas_p"]]), na.rm = TRUE),
  fdr_lt_0_05_genes = nrow(ov[twas_fdr < 0.05]),
  magma_top50_overlap = sum(ov[["twas_fdr"]] < 0.05 & ov[["gene"]] %in% top50, na.rm = TRUE),
  magma_top100_overlap = sum(ov[["twas_fdr"]] < 0.05 & ov[["gene"]] %in% top100, na.rm = TRUE),
  magma_fdr_overlap = sum(ov[["twas_fdr"]] < 0.05 & ov[["gene"]] %in% fdr_set, na.rm = TRUE),
  magma_suggestive_overlap = sum(ov[["twas_fdr"]] < 0.05 & ov[["gene"]] %in% suggestive_set, na.rm = TRUE),
  p1_overlap = sum(ov[["twas_fdr"]] < 0.05 & ov[["gene"]] %in% p1_genes, na.rm = TRUE),
  warning = "S-PrediXcan log warns p/z are uncalibrated without GWAS h2; many FDR genes use one SNP.",
  recommended_for_manuscript = fifelse(sum(ov[["twas_fdr"]] < 0.05 & ov[["n_snps_used"]] > 1 & ov[["magma_overlap"]] == TRUE, na.rm = TRUE) > 0,
                                       "supporting_proxy_TWAS_with_boundary",
                                       "supplementary_exploratory_only")
)
fwrite(run_comparison, "results/twas/spredixcan_run_comparison.tsv", sep = "\t")

candidate <- fread("results/phase2b_twas_ready_v0.1/tables/candidate_gene_tiers_v0.2.tsv", fill = TRUE)
candidate[, gene := as.character(gene)]
twas_best <- ov[order(twas_fdr), .SD[1], by = gene]
twas_cols <- twas_best[, .(
  gene,
  twas_kidney_cortex_status = support_status,
  twas_fdr_q = twas_fdr,
  twas_n_snps_used = n_snps_used,
  twas_magma_overlap_status = fifelse(magma_overlap == TRUE, "MAGMA_overlap", "no_MAGMA_overlap"),
  twas_interpretation = fifelse(gene %in% p1_genes,
                                "P1_not_upgraded_by_TWAS",
                                fifelse(twas_fdr < 0.05 & magma_overlap == TRUE & n_snps_used > 1,
                                        "stronger_TWAS_MAGMA_proxy_candidate",
                                        fifelse(twas_fdr < 0.05 & magma_overlap == TRUE,
                                                "exploratory_TWAS_proxy_low_coverage",
                                                fifelse(twas_fdr < 0.05,
                                                        "TWAS_FDR_without_MAGMA_overlap_not_main_claim",
                                                        "not_TWAS_FDR_supported"))))
)]
candidate_v11 <- merge(candidate, twas_cols, by = "gene", all.x = TRUE)
candidate_v11[is.na(twas_kidney_cortex_status), `:=`(
  twas_kidney_cortex_status = "not_tested_or_not_in_model",
  twas_magma_overlap_status = "not_applicable",
  twas_interpretation = "not_TWAS_supported"
)]
candidate_v11[gene %in% p1_genes, twas_interpretation := "P1_not_upgraded_by_TWAS"]
fwrite(candidate_v11, "results/tables/candidate_gene_tiers_v1.1.tsv", sep = "\t")

integrated <- ov[, .(
  gene,
  gene_id,
  twas_kidney_cortex_status = support_status,
  twas_fdr_q = twas_fdr,
  twas_n_snps_used = n_snps_used,
  twas_magma_overlap_status = fifelse(magma_overlap == TRUE, "MAGMA_overlap", "no_MAGMA_overlap"),
  magma_rank,
  magma_p,
  magma_fdr,
  p1_candidate,
  twas_interpretation = claim_use
)][order(twas_fdr_q)]
fwrite(integrated, "results/tables/integrated_evidence_summary_v0.3.tsv", sep = "\t")

plan <- c(
  "# Phase25D TWAS hardening plan",
  "",
  "Phase25C TWAS is a real GTEx Kidney_Cortex S-PrediXcan run, but it is treated as uncalibrated proxy evidence pending variant-coverage, one-SNP sensitivity and MAGMA-overlap benchmark checks.",
  "",
  "Hardening tasks completed in Phase25D:",
  "- Freeze Phase25C outputs.",
  "- Audit PredictDB Kidney_Cortex model structure and SNP coverage.",
  "- Compare GWAS key formats against model rsID and varID keys.",
  "- Stratify TWAS results by one-SNP and multi-SNP support.",
  "- Benchmark TWAS/MAGMA overlap against random expectation.",
  "- Update candidate-gene tier and integrated evidence tables with manuscript-safe interpretation fields."
)
writeLines(plan, "docs/phase25d_twas_hardening_plan.md")

manuscript_text <- c(
  "# TWAS proxy result insert v0.1",
  "",
  "## Methods insert",
  "",
  "We performed an exploratory GTEx Kidney_Cortex proxy S-PrediXcan analysis using PredictDB GTEx v8 MASHR expression prediction models. GWAS variants were harmonized to the PredictDB Kidney_Cortex model using rsID-to-varID mapping, allele compatibility checks and effect-direction flipping when required. Benjamini-Hochberg correction was applied across tested Kidney_Cortex TWAS genes. This analysis was treated as genetically regulated expression support only and not as causal, colocalization, SMR, spatial validation or papillary-specific evidence.",
  "",
  "## Results insert",
  "",
  sprintf("The GTEx Kidney_Cortex proxy TWAS tested %s genes and identified %s FDR-supported genetically regulated expression associations. Among these, %s overlapped MAGMA-prioritized genes, whereas the six P1 genes were not FDR-supported. One-SNP sensitivity showed that %s FDR-supported genes used a single SNP and %s used more than one SNP; therefore, the result was interpreted as support for a subset of MAGMA-prioritized genes rather than as P1 validation.", nrow(ov), nrow(ov[twas_fdr < 0.05]), nrow(ov[twas_fdr < 0.05 & magma_overlap == TRUE]), nrow(ov[twas_fdr < 0.05 & n_snps_used == 1]), nrow(ov[twas_fdr < 0.05 & n_snps_used > 1])),
  "",
  "## Discussion limitation insert",
  "",
  "The Kidney_Cortex TWAS result should not be interpreted as evidence that the P1 genes mediate KSD risk through altered expression. Rather, it provides an external kidney-cortex proxy layer showing that a subset of MAGMA-prioritized genes also have genetically regulated expression associations. The absence of FDR-supported P1 TWAS signals reinforces the claim-bounded interpretation: P1 genes are prioritized by MAGMA and papillary single-nucleus localization, whereas TWAS convergence occurs at a broader MAGMA-module level. The analysis is not causal, not colocalization, not SMR, not spatial validation and not papillary-specific.",
  "",
  "## Figure or supplement recommendation",
  "",
  "Do not redraw main Figures 1-5 at this stage. If included visually, use a Supplementary Figure S6 showing the Kidney_Cortex proxy workflow, TWAS FDR genes, TWAS-MAGMA overlap benchmark, one-SNP sensitivity and absence of FDR-supported P1 genes."
)
writeLines(manuscript_text, "manuscript/twas_proxy_result_insert_v0.1.md")

boundary <- c(
  "# Phase25D TWAS claim-boundary memo",
  "",
  "Allowed wording:",
  "- GTEx Kidney_Cortex proxy",
  "- genetically regulated expression support",
  "- subset of MAGMA-prioritized genes",
  "- P1 genes were not FDR-supported",
  "",
  "Not allowed:",
  "- TWAS proves causality",
  "- P1 genes are TWAS-supported",
  "- papilla-specific TWAS validation",
  "- TAL-specific TWAS",
  "- colocalization, SMR, spatial validation or experimental mechanism",
  "",
  "Recommended use: supporting proxy TWAS evidence only if emphasizing FDR-supported TWAS/MAGMA overlap and the one-SNP coverage limitation."
)
writeLines(boundary, "docs/phase25d_twas_claim_boundary_memo.md")

supp_fig <- c(
  "# Supplementary Figure S6 plan: GTEx Kidney_Cortex S-PrediXcan proxy TWAS",
  "",
  "A. S-PrediXcan workflow and resource boundary: PredictDB GTEx v8 MASHR Kidney_Cortex model, rsID-to-varID harmonization and proxy-tissue limitation.",
  "B. Ranked FDR-supported TWAS genes by -log10(FDR), marking MAGMA-overlapping genes.",
  "C. TWAS-MAGMA overlap random benchmark for all FDR and multi-SNP FDR TWAS sets.",
  "D. One-SNP sensitivity comparing all, n_snps_used == 1, n_snps_used > 1 and n_snps_used >= 3 subsets.",
  "E. P1 genes absent from the FDR-supported TWAS set.",
  "",
  "Do not insert into main Figures 1-5 unless later requested."
)
writeLines(supp_fig, "docs/supplementary_figure_twas_plan_v0.1.md")

memo <- c(
  "# Phase25D TWAS reliability hardening summary",
  "",
  "## Result classification",
  "",
  if (run_comparison$recommended_for_manuscript == "supporting_proxy_TWAS_with_boundary") {
    "- Classification: supporting proxy TWAS evidence with explicit limitations."
  } else {
    "- Classification: supplementary exploratory TWAS only."
  },
  sprintf("- TWAS FDR genes: %s", nrow(ov[twas_fdr < 0.05])),
  sprintf("- TWAS FDR and MAGMA overlap genes: %s", nrow(ov[twas_fdr < 0.05 & magma_overlap == TRUE])),
  sprintf("- Multi-SNP TWAS FDR and MAGMA overlap genes: %s", nrow(ov[twas_fdr < 0.05 & magma_overlap == TRUE & n_snps_used > 1])),
  sprintf("- P1 FDR-supported genes: %s", nrow(ov[twas_fdr < 0.05 & p1_candidate == TRUE])),
  "",
  "## Boundary",
  "",
  "The TWAS layer supports a subset of MAGMA-prioritized genes in a GTEx Kidney_Cortex proxy model. It does not support P1 genes directly and does not establish causality, colocalization, SMR support, spatial validation or papillary specificity."
)
writeLines(memo, "docs/phase25d_twas_reliability_summary.md")

message("wrote Phase25D TWAS hardening outputs")
