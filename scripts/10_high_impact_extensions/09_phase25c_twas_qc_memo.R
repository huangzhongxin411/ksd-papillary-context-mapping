suppressPackageStartupMessages(library(data.table))

dir.create("docs", showWarnings = FALSE)
dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)

summary <- fread("results/twas/twas_results_real_summary.tsv", fill = TRUE)
res <- fread("results/twas/twas_results.tsv", fill = TRUE)
ov <- fread("results/twas/twas_magma_overlap_real.tsv", fill = TRUE)
map_qc <- fread("results/twas/ksd_2025_for_twas.Kidney_Cortex_varID_qc.tsv", fill = TRUE)
resource <- fread("results/twas/twas_resource_status.tsv", fill = TRUE)

get_metric <- function(name) {
  x <- map_qc[metric == name, value]
  if (length(x)) x[1] else NA_character_
}

warnings <- c()
log_path <- "logs/twas/spredixcan_Kidney_Cortex_KSD.log"
if (file.exists(log_path)) {
  log_txt <- readLines(log_path, warn = FALSE)
  warnings <- grep("WARNING|IMPORTANT|Missing", log_txt, value = TRUE)
}

qc <- data.table(
  item = c(
    "run_completed",
    "predictdb_tar_present",
    "kidney_cortex_model_found",
    "kidney_cortex_covariance_found",
    "gwas_rows",
    "model_rsid_overlap",
    "allele_compatible_overlap",
    "mapped_varid_rows",
    "flipped_rows",
    "ambiguous_multi_mapping_skipped",
    "n_genes_tested",
    "n_genes_with_usable_snps",
    "n_genes_with_single_snp_model_used",
    "n_genes_with_gt1_snp_used",
    "min_p",
    "n_fdr_lt_0_05",
    "n_fdr_magma_overlap",
    "n_fdr_p1_overlap",
    "harmonization_warning",
    "promotion_status"
  ),
  value = c(
    summary$status[1] == "completed",
    resource[resource == "PredictDB MASHR tar", status][1],
    resource[resource == "PredictDB Kidney_Cortex model db", status][1],
    resource[resource == "PredictDB Kidney_Cortex covariance", status][1],
    get_metric("gwas_rows"),
    get_metric("rsid_overlap"),
    get_metric("allele_compatible"),
    get_metric("written_rows"),
    get_metric("flipped_rows"),
    get_metric("ambiguous_multi_mapping_skipped"),
    summary$n_genes_tested[1],
    summary$n_genes_with_n_snps_used_gt0[1],
    nrow(res[n_snps_used == 1]),
    nrow(res[n_snps_used > 1]),
    summary$min_p[1],
    summary$n_fdr_lt_0_05[1],
    summary$n_magma_overlap_fdr_lt_0_05[1],
    summary$n_p1_fdr_lt_0_05[1],
    if (length(warnings)) paste(warnings, collapse = " | ") else "none",
    if (summary$status[1] == "completed" && summary$n_fdr_lt_0_05[1] > 0) {
      "TWAS_FDR_supported_available_for_cautious_enhancement"
    } else {
      "not_promoted"
    }
  ),
  note = c(
    "S-PrediXcan process completed and produced non-empty output.",
    "Manual tar was validated as POSIX tar and unpacked.",
    "Kidney_Cortex model db found after unpacking.",
    "Kidney_Cortex covariance found after unpacking.",
    "Original TWAS GWAS input rows.",
    "GWAS rsIDs overlapping Kidney_Cortex model rsIDs before allele filtering.",
    "Overlapping rows with compatible allele sets.",
    "Mapped varID GWAS rows used for S-PrediXcan.",
    "Rows whose beta/Z/EAF were flipped to model effect allele.",
    "Ambiguous rsID-to-varID mappings excluded.",
    "Rows in parsed S-PrediXcan result table.",
    "Genes with n_snps_used > 0.",
    "Single-SNP models require cautious interpretation.",
    "Multi-SNP models with >1 SNP used.",
    "Minimum TWAS P value.",
    "Benjamini-Hochberg FDR across Kidney_Cortex TWAS genes.",
    "FDR-supported TWAS genes overlapping MAGMA gene results.",
    "FDR-supported TWAS genes among six P1 candidates.",
    "Log warnings retained for interpretation boundary.",
    "Promotion means enhanced genetic evidence only, not causality, colocalization or validation."
  )
)
fwrite(qc, "results/twas/phase25c_spredixcan_kidney_cortex_qc.tsv", sep = "\t")

top <- ov[twas_fdr < 0.05][order(twas_fdr)]
top_magma <- top[magma_overlap == TRUE]
top_p1 <- top[p1_candidate == TRUE]

top_lines <- if (nrow(top_magma)) {
  paste0(
    "- ", head(top_magma$gene, 15),
    " (TWAS FDR=", signif(head(top_magma$twas_fdr, 15), 3),
    "; MAGMA rank=", head(top_magma$magma_rank, 15), ")"
  )
} else {
  "- No FDR-supported TWAS/MAGMA overlap."
}

memo <- c(
  "# Phase25C Kidney_Cortex S-PrediXcan real run memo",
  "",
  "## Resource verification",
  "",
  "- `mashr_eqtl.tar` exists, has plausible size (~250 MB), is a POSIX tar archive, and contains `mashr_Kidney_Cortex.db` plus `mashr_Kidney_Cortex.txt.gz`.",
  "- The archive was unpacked under `external/twas/predixcan/predictdb_gtex_v8_mashr/`.",
  "- S-PrediXcan ran with the project-local Python environment and PredictDB MASHR Kidney_Cortex files.",
  "",
  "## Harmonization",
  "",
  paste0("- Original GWAS rows: ", get_metric("gwas_rows")),
  paste0("- Model rsID overlap: ", get_metric("rsid_overlap")),
  paste0("- Allele-compatible overlap: ", get_metric("allele_compatible")),
  paste0("- Mapped varID rows written: ", get_metric("written_rows")),
  paste0("- Flipped rows: ", get_metric("flipped_rows")),
  paste0("- Ambiguous mappings skipped: ", get_metric("ambiguous_multi_mapping_skipped")),
  "",
  "## TWAS result summary",
  "",
  paste0("- Run completed: ", summary$status[1]),
  paste0("- Genes tested: ", summary$n_genes_tested[1]),
  paste0("- Genes with usable SNPs: ", summary$n_genes_with_n_snps_used_gt0[1]),
  paste0("- Minimum P: ", signif(summary$min_p[1], 4)),
  paste0("- FDR q < 0.05 genes: ", summary$n_fdr_lt_0_05[1]),
  paste0("- FDR-supported MAGMA overlaps: ", summary$n_magma_overlap_fdr_lt_0_05[1]),
  paste0("- FDR-supported P1 overlaps: ", summary$n_p1_fdr_lt_0_05[1]),
  paste0("- Genes with one SNP used: ", nrow(res[n_snps_used == 1])),
  paste0("- Genes with >1 SNP used: ", nrow(res[n_snps_used > 1])),
  "",
  "## Leading FDR-supported TWAS/MAGMA overlaps",
  "",
  top_lines,
  "",
  "## Warnings and interpretation boundary",
  "",
  if (length(warnings)) paste0("- ", warnings) else "- No warning lines were detected in the S-PrediXcan log.",
  "- TWAS support may be described only as GTEx Kidney_Cortex genetically regulated expression support for FDR-supported genes, preferably emphasizing the TWAS/MAGMA-overlapping subset.",
  "- This does not establish causality, colocalization, SMR support, spatial validation, cell-type-specific disease expression or experimental mechanism.",
  "- P1 genes were not FDR-supported in this Kidney_Cortex TWAS run."
)
writeLines(memo, "docs/phase25c_spredixcan_kidney_cortex_real_run_memo.md")

message("wrote Phase25C QC and memo")
