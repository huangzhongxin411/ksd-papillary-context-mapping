suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

exact <- fread(file.path(table_dir, "gse73680_exact_statistics_summary.tsv"))
loo <- fread(file.path(table_dir, "gse73680_module_leave_one_gene_out.tsv"))
without_p1 <- fread(file.path(table_dir, "gse73680_magma_without_p1_sensitivity.tsv"))
expr_match <- fread(file.path(table_dir, "gse73680_expression_matched_random_benchmark.tsv"))
direction <- fread(file.path(table_dir, "gse73680_paired_direction_consistency.tsv"))

loo_sum <- loo[, .(
  min_retention = min(effect_retention_fraction, na.rm = TRUE),
  frac_direction_preserved = mean(direction_preserved, na.rm = TRUE),
  frac_robust = mean(interpretation == "robust", na.rm = TRUE)
), by = module_name]

summary <- rbindlist(list(
  data.table(
    evidence_item = "single_gene_P1_response",
    result = "No P1 gene reached FDR support; PKD2 was nominal only in paired sensitivity.",
    strength = "negative_or_exploratory",
    interpretation = "P1 single-gene responses are heterogeneous and not FDR-supported.",
    claim_boundary = "Do not claim P1 genes are differentially expressed in plaque."
  ),
  data.table(
    evidence_item = "MAGMA_module_response",
    result = paste(exact[analysis_type == "patient_level_module_paired_sensitivity" &
                           feature_or_module %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive"),
                         paste0(feature_or_module, ":q=", fdr_rounded)], collapse = "; "),
    strength = "strong_module_context_support",
    interpretation = "MAGMA modules show paired patient-level disease-context shifts with q <= 0.05.",
    claim_boundary = "Module-level expression association only; no causal genetic mechanism."
  ),
  data.table(
    evidence_item = "leave_one_gene_out",
    result = paste(loo_sum[, paste0(module_name, ":min_retention=", sprintf("%.2f", min_retention),
                                    ",dir=", sprintf("%.0f%%", 100 * frac_direction_preserved))], collapse = "; "),
    strength = ifelse(all(loo_sum$min_retention >= 0.80 & loo_sum$frac_direction_preserved == 1), "strong_module_context_support", "moderate_module_context_support"),
    interpretation = "MAGMA module shifts are preserved after removing individual genes.",
    claim_boundary = "Robustness check does not prove causality."
  ),
  data.table(
    evidence_item = "MAGMA_without_P1",
    result = paste(without_p1[, paste0(module_name, ":q=", sprintf("%.3f", fdr_without_p1),
                                       ",pct=", sprintf("%.3f", random_benchmark_percentile_without_p1))], collapse = "; "),
    strength = ifelse(all(without_p1$interpretation == "robust_without_P1"), "strong_module_context_support", "moderate_module_context_support"),
    interpretation = "MAGMA module signal remains after removing the six P1 genes.",
    claim_boundary = "Disease-context module signal is not solely driven by P1 genes."
  ),
  data.table(
    evidence_item = "expression_matched_benchmark",
    result = paste(expr_match[, paste0(module_name, ":pct=", sprintf("%.3f", expression_matched_percentile),
                                       ",emp.P=", sprintf("%.3f", empirical_p))], collapse = "; "),
    strength = ifelse(any(expr_match[module_name %in% c("MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]$interpretation ==
                            "robust_beyond_expression_level_background"), "strong_module_context_support", "boundary_check_not_supported"),
    interpretation = ifelse(any(expr_match[module_name %in% c("MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]$interpretation ==
                                "robust_beyond_expression_level_background"),
                            "Selected MAGMA modules exceed expression-matched random expectations.",
                            "MAGMA modules do not exceed the stricter expression-matched random benchmark under the current implementation."),
    claim_boundary = "Expression-matched benchmark should be reported as a boundary check, not primary support, if not exceeded."
  ),
  data.table(
    evidence_item = "direction_consistency",
    result = paste(direction[feature_type == "module" & module_or_gene %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive"),
                             paste0(module_or_gene, ":positive=", sprintf("%.0f%%", 100 * positive_fraction))], collapse = "; "),
    strength = "moderate_module_context_support",
    interpretation = "MAGMA_top50, top100 and FDR show directionally consistent positive paired deltas; MAGMA_suggestive is moderate.",
    claim_boundary = "Direction consistency is supportive but still observational."
  ),
  data.table(
    evidence_item = "claim_boundary",
    result = "GSE73680 supports disease-context expression association for MAGMA-prioritized modules rather than uniform P1 single-gene differential expression.",
    strength = "moderate_module_context_support",
    interpretation = "Suitable as Result 4 candidate with careful wording.",
    claim_boundary = "Do not claim causal validation, P1 gene validation, TAL mechanism confirmation, spatial validation, TWAS convergence or colocalization."
  )
), fill = TRUE)

fwrite(summary, file.path(table_dir, "gse73680_disease_context_summary_v0.2.tsv"), sep = "\t")

writeLines(c(
  "# GSE73680 Disease-Context Results Memo v0.2",
  "",
  "## Main Conclusion",
  "",
  "GSE73680 provides disease-context expression support at the MAGMA module level, rather than uniform P1 single-gene differential expression.",
  "",
  "## Exact Statistics",
  "",
  "- Patient-level MAGMA module sensitivity reached q <= 0.05 for MAGMA_top50, MAGMA_top100, MAGMA_FDR and MAGMA_suggestive.",
  "- Sample-level module analysis was borderline rather than conventionally FDR-supported.",
  "- P1 single-gene analysis was not FDR-supported; PKD2 was nominal only in paired sensitivity.",
  "",
  "## Robustness",
  "",
  "- Leave-one-gene-out sensitivity preserved MAGMA module effect direction.",
  "- MAGMA-without-P1 sensitivity remained supported, indicating the module signal is not solely driven by the six P1 genes.",
  "- Size-matched random benchmarking supported MAGMA modules, whereas the stricter expression-matched benchmark did not provide primary support under the current implementation.",
  "- Direction consistency was strongest for MAGMA_top50, MAGMA_top100 and MAGMA_FDR.",
  "",
  "## Recommended Result 4 Wording",
  "",
  "GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression.",
  "",
  "## Claim Boundary",
  "",
  "Allowed: disease-context expression association for MAGMA-prioritized modules.",
  "",
  "Not allowed: causal validation, P1 gene validation, TAL mechanism confirmation, spatial validation, TWAS convergence or colocalization."
), "docs/gse73680_disease_context_results_memo_v0.2.md")

message("wrote GSE73680 robustness summary v0.2")
