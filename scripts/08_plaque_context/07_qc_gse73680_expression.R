suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
processed_dir <- "data/processed/gse73680"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_matrix_path <- file.path(processed_dir, "gse73680_gene_expression_matrix.tsv.gz")
meta_path <- "config/gse73680_sample_metadata_curated.tsv"
if (!file.exists(gene_matrix_path) || !file.exists(meta_path)) {
  fwrite(data.table(n_samples = 0L, n_genes = 0L, n_samples_included = 0L,
                    n_genes_detected_in_80pct_samples = 0L, overall_missing_rate = NA_real_,
                    scale_status = "not_available", normalization_status = "not_available",
                    n_outlier_samples = NA_integer_, pca_group_separation_visible = NA,
                    batch_effect_suspected = NA, analysis_ready_qc = FALSE,
                    status = "fail", notes = "Gene expression matrix or curated metadata missing."),
         file.path(table_dir, "gse73680_expression_qc_v2.tsv"), sep = "\t")
  fwrite(data.table(), file.path(table_dir, "gse73680_sample_qc_metrics.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
mat <- fread(gene_matrix_path)
meta <- fread(meta_path)
sample_cols <- intersect(setdiff(names(mat), "gene"), meta[include_in_analysis == TRUE]$sample_id)
if (!length(sample_cols)) {
  fwrite(data.table(n_samples = length(setdiff(names(mat), "gene")), n_genes = nrow(mat), n_samples_included = 0L,
                    n_genes_detected_in_80pct_samples = 0L, overall_missing_rate = NA_real_,
                    scale_status = "not_available", normalization_status = "not_available",
                    n_outlier_samples = NA_integer_, pca_group_separation_visible = NA,
                    batch_effect_suspected = NA, analysis_ready_qc = FALSE, status = "fail",
                    notes = "No samples included after metadata curation."),
         file.path(table_dir, "gse73680_expression_qc_v2.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
vals <- as.matrix(mat[, ..sample_cols])
sample_metrics <- rbindlist(lapply(sample_cols, function(s) {
  v <- mat[[s]]
  data.table(sample_id = s, n_genes_detected = sum(!is.na(v)), missing_rate = mean(is.na(v)),
             median_expression = median(v, na.rm = TRUE), iqr_expression = IQR(v, na.rm = TRUE),
             mean_expression = mean(v, na.rm = TRUE), sd_expression = sd(v, na.rm = TRUE))
}))
sample_metrics <- merge(sample_metrics, meta[, .(sample_id, group_curated)], by = "sample_id", all.x = TRUE)
sample_metrics[, outlier_flag := abs(scale(median_expression)) > 3]
sample_metrics[, include_in_qc := TRUE]
sample_metrics[, exclude_reason := NA_character_]
fwrite(sample_metrics, file.path(table_dir, "gse73680_sample_qc_metrics.tsv"), sep = "\t")
overall_missing <- mean(is.na(vals))
scale_status <- if (max(vals, na.rm = TRUE) <= 25 && min(vals, na.rm = TRUE) >= -5) "likely_log2_or_normalized" else if (all(abs(vals - round(vals)) < 1e-8, na.rm = TRUE) && max(vals, na.rm = TRUE) > 100) "raw_count_or_unlogged" else "normalized_continuous_or_unknown"
analysis_ready <- ncol(vals) >= 6 && nrow(vals) >= 8000 && overall_missing <= 0.15 && scale_status != "raw_count_or_unlogged"
status <- if (analysis_ready && overall_missing <= 0.05) "pass" else if (ncol(vals) >= 4 && nrow(vals) >= 5000 && overall_missing <= 0.15) "warning" else "fail"
fwrite(data.table(n_samples = length(setdiff(names(mat), "gene")), n_genes = nrow(mat), n_samples_included = ncol(vals),
                  n_genes_detected_in_80pct_samples = sum(rowMeans(!is.na(vals)) >= 0.8),
                  overall_missing_rate = overall_missing, scale_status = scale_status,
                  normalization_status = "not_independently_verified", n_outlier_samples = sum(sample_metrics$outlier_flag),
                  pca_group_separation_visible = NA, batch_effect_suspected = NA,
                  analysis_ready_qc = status %in% c("pass", "warning"), status = status,
                  notes = "QC summary generated; visual plots only if matrix is available."),
       file.path(table_dir, "gse73680_expression_qc_v2.tsv"), sep = "\t")
pdf(file.path(fig_dir, "gse73680_boxplot.pdf")); boxplot(vals, las = 2, main = "GSE73680 expression distribution"); dev.off()
pdf(file.path(fig_dir, "gse73680_density_plot.pdf")); plot(density(vals[,1], na.rm = TRUE), main = "GSE73680 density preview"); dev.off()
pdf(file.path(fig_dir, "gse73680_pca.pdf")); pc <- prcomp(t(na.omit(vals)), scale. = TRUE); plot(pc$x[,1:2], main = "GSE73680 PCA"); dev.off()
pdf(file.path(fig_dir, "gse73680_sample_correlation_heatmap.pdf")); heatmap(cor(vals, use = "pairwise.complete.obs"), main = "Sample correlation"); dev.off()
message("wrote GSE73680 expression QC outputs")
