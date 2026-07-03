suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

label_sig <- function(p, fdr) {
  fcase(
    !is.na(fdr) & fdr <= 0.05, "FDR_supported",
    !is.na(fdr) & fdr > 0.05 & fdr <= 0.10, "borderline_FDR",
    !is.na(p) & p < 0.05, "nominal_only",
    default = "not_detectable"
  )
}

rows <- list()
add_stats <- function(path, analysis_type, feature_col, effect_col, p_col, fdr_col, method_col = "method", notes = "") {
  if (!file.exists(path)) return(NULL)
  x <- fread(path)
  if (!nrow(x)) return(NULL)
  data.table(
    analysis_type = analysis_type,
    feature_or_module = as.character(x[[feature_col]]),
    method = if (method_col %in% names(x)) as.character(x[[method_col]]) else analysis_type,
    effect_size = as.numeric(x[[effect_col]]),
    p_value = if (p_col %in% names(x)) as.numeric(x[[p_col]]) else NA_real_,
    fdr_exact = if (fdr_col %in% names(x)) as.numeric(x[[fdr_col]]) else NA_real_,
    fdr_rounded = if (fdr_col %in% names(x)) sprintf("%.3f", as.numeric(x[[fdr_col]])) else NA_character_,
    notes = notes
  )
}

rows[[length(rows) + 1]] <- add_stats(file.path(table_dir, "gse73680_p1_gene_response.tsv"),
                                      "sample_level_P1_limma_duplicateCorrelation", "gene", "logFC", "p_value", "fdr",
                                      notes = "Sample-level limma duplicateCorrelation P1 gene response.")
rows[[length(rows) + 1]] <- add_stats(file.path(table_dir, "gse73680_patient_level_p1_gene_response.tsv"),
                                      "patient_level_P1_paired_sensitivity", "gene", "paired_delta", "p_value", "fdr",
                                      notes = "Patient-level paired P1 sensitivity.")
rows[[length(rows) + 1]] <- add_stats(file.path(table_dir, "gse73680_module_score_response.tsv"),
                                      "sample_level_module_limma_duplicateCorrelation", "module_name", "delta", "p_value", "fdr",
                                      notes = "Sample-level module score limma duplicateCorrelation.")
rows[[length(rows) + 1]] <- add_stats(file.path(table_dir, "gse73680_patient_level_module_response.tsv"),
                                      "patient_level_module_paired_sensitivity", "module_name", "paired_delta", "p_value", "fdr",
                                      notes = "Patient-level paired module sensitivity.")
if (file.exists(file.path(table_dir, "gse73680_random_module_benchmark.tsv"))) {
  x <- fread(file.path(table_dir, "gse73680_random_module_benchmark.tsv"))
  rows[[length(rows) + 1]] <- data.table(
    analysis_type = "size_matched_random_benchmark",
    feature_or_module = x$module_name,
    method = "1000_size_matched_random_gene_sets",
    effect_size = x$observed_delta,
    p_value = x$empirical_p,
    fdr_exact = NA_real_,
    fdr_rounded = NA_character_,
    notes = paste0("percentile=", sprintf("%.3f", x$percentile), "; ", x$benchmark_interpretation)
  )
}
if (file.exists(file.path(table_dir, "gse73680_expression_matched_random_benchmark.tsv"))) {
  x <- fread(file.path(table_dir, "gse73680_expression_matched_random_benchmark.tsv"))
  rows[[length(rows) + 1]] <- data.table(
    analysis_type = "expression_matched_random_benchmark",
    feature_or_module = x$module_name,
    method = "1000_expression_matched_random_gene_sets",
    effect_size = x$observed_delta,
    p_value = x$empirical_p,
    fdr_exact = NA_real_,
    fdr_rounded = NA_character_,
    notes = paste0("percentile=", sprintf("%.3f", x$expression_matched_percentile), "; ", x$interpretation)
  )
}

out <- rbindlist(rows, fill = TRUE)
out[, significance_label := label_sig(p_value, fdr_exact)]
out[analysis_type %like% "random_benchmark|matched_random", significance_label := fifelse(
  !is.na(p_value) & p_value < 0.05, "FDR_supported",
  fifelse(!is.na(p_value) & p_value <= 0.10, "borderline_FDR", "not_detectable")
)]
setcolorder(out, c("analysis_type", "feature_or_module", "method", "effect_size", "p_value",
                   "fdr_exact", "fdr_rounded", "significance_label", "notes"))
fwrite(out, file.path(table_dir, "gse73680_exact_statistics_summary.tsv"), sep = "\t")
message("wrote exact statistics summary")
