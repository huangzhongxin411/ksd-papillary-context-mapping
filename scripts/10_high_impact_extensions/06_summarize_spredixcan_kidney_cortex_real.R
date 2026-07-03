suppressPackageStartupMessages(library(data.table))

out_file <- "results/twas/spredixcan_Kidney_Cortex_KSD.csv"
dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)

p1 <- c("UMOD", "CLDN14", "CASR", "CLDN10", "HIBADH", "PKD2")

if (!file.exists(out_file) || file.info(out_file)$size == 0) {
  fwrite(data.table(status = "not_run_or_empty", tissue = "Kidney_Cortex"),
         "results/twas/twas_results_real_summary.tsv", sep = "\t")
  quit(save = "no", status = 0)
}

x <- fread(out_file, fill = TRUE)
gene_col <- intersect(c("gene", "gene_name", "gene_name_in_model", "GENE"), names(x))[1]
gene_symbol_col <- intersect(c("gene_name", "gene_name_in_model", "GENE_NAME", "symbol"), names(x))[1]
p_col <- intersect(c("pvalue", "p_value", "P", "p"), names(x))[1]
z_col <- intersect(c("zscore", "z_score", "z"), names(x))[1]
nsnp_col <- intersect(c("n_snps_used", "n_snps_in_model", "n_snps"), names(x))[1]

if (is.na(gene_col) || is.na(p_col)) {
  fwrite(data.table(status = "failed_parse", tissue = "Kidney_Cortex",
                    columns = paste(names(x), collapse = ",")),
         "results/twas/twas_results_real_summary.tsv", sep = "\t")
  quit(save = "no", status = 0)
}

res <- data.table(
  tissue = "Kidney_Cortex",
  gene_id = as.character(x[[gene_col]]),
  gene = if (!is.na(gene_symbol_col)) as.character(x[[gene_symbol_col]]) else as.character(x[[gene_col]]),
  twas_p = as.numeric(x[[p_col]]),
  twas_z = if (!is.na(z_col)) as.numeric(x[[z_col]]) else NA_real_,
  n_snps_used = if (!is.na(nsnp_col)) as.numeric(x[[nsnp_col]]) else NA_real_
)
res <- res[!is.na(twas_p)]
res[, twas_fdr := p.adjust(twas_p, method = "BH")]
res[, p1_candidate := gene %in% p1]
res[, support_status := fifelse(twas_fdr < 0.05, "FDR_supported", fifelse(twas_p < 1e-4, "suggestive_only", "not_supported"))]
res[, claim_use := fifelse(twas_fdr < 0.05, "eligible_as_TWAS_enhancement", "not_used_for_main_claim")]
fwrite(res[order(twas_fdr, twas_p)], "results/twas/twas_results.tsv", sep = "\t")

magma_path <- "results/phase2_magma_v0.1/tables/magma_genes.tsv"
if (!file.exists(magma_path)) magma_path <- "submission_package_v1.9_BMC_working/supplementary/tables/S3_magma_gene_results_and_gene_sets.tsv"
magma <- fread(magma_path, fill = TRUE)
if ("gene_symbol" %in% names(magma)) setnames(magma, "gene_symbol", "gene")
keep <- intersect(c("gene", "rank", "p", "fdr"), names(magma))
magma <- unique(magma[, ..keep])
if ("rank" %in% names(magma)) setnames(magma, "rank", "magma_rank")
if ("p" %in% names(magma)) setnames(magma, "p", "magma_p")
if ("fdr" %in% names(magma)) setnames(magma, "fdr", "magma_fdr")

ov <- merge(res, magma, by = "gene", all.x = TRUE)
ov[, magma_overlap := !is.na(magma_rank)]
ov[, claim_use := fifelse(twas_fdr < 0.05 & magma_overlap, "eligible_as_TWAS_MAGMA_enhancement", "not_used_for_main_claim")]
fwrite(ov[order(twas_fdr, twas_p)], "results/twas/twas_magma_overlap_real.tsv", sep = "\t")

summary <- data.table(
  status = "completed",
  tissue = "Kidney_Cortex",
  n_genes_tested = nrow(res),
  n_genes_with_n_snps_used_gt0 = if ("n_snps_used" %in% names(res) && any(!is.na(res$n_snps_used))) nrow(res[n_snps_used > 0]) else NA_integer_,
  min_p = min(res$twas_p, na.rm = TRUE),
  n_fdr_lt_0_05 = nrow(res[twas_fdr < 0.05]),
  n_magma_overlap_fdr_lt_0_05 = nrow(ov[twas_fdr < 0.05 & magma_overlap == TRUE]),
  n_p1_fdr_lt_0_05 = nrow(res[twas_fdr < 0.05 & p1_candidate == TRUE]),
  claim_use = ifelse(nrow(res[twas_fdr < 0.05]) > 0, "FDR_supported_results_available_for_review", "not_FDR_supported_do_not_promote"),
  blocking_reason = ""
)
fwrite(summary, "results/twas/twas_results_real_summary.tsv", sep = "\t")

memo <- c(
  "# Phase25B Kidney_Cortex S-PrediXcan real run memo",
  "",
  paste0("- Status: ", summary$status),
  paste0("- Genes tested: ", summary$n_genes_tested),
  paste0("- Minimum P: ", signif(summary$min_p, 4)),
  paste0("- FDR-supported genes: ", summary$n_fdr_lt_0_05),
  paste0("- FDR-supported MAGMA overlaps: ", summary$n_magma_overlap_fdr_lt_0_05),
  "",
  "Claim rule: only FDR-supported TWAS results, especially MAGMA-overlapping genes, may be described as enhanced genetic evidence."
)
writeLines(memo, "docs/phase25b_twas_real_run_memo.md")
