suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
processed_dir <- "data/processed/gse73680"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(processed_dir, "gse73680_expression_matrix.tsv.gz")
if (!file.exists(matrix_path)) {
  fwrite(data.table(n_features_input = 0L, n_features_mapped = 0L, n_features_unmapped = 0L,
                    mapping_rate = 0, n_unique_genes = 0L, n_multi_probe_genes = 0L,
                    collapse_rule = "max_mean_expression_probe", status = "fail",
                    notes = "Expression matrix missing; gene mapping not run."),
         file.path(table_dir, "gse73680_gene_mapping_qc.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
mat <- fread(matrix_path)
feature_ids <- mat$feature_id
is_symbol <- grepl("^[A-Z0-9][A-Z0-9.-]{1,15}$", feature_ids)
mapping <- data.table(feature_id = feature_ids,
                      gene_symbol = ifelse(is_symbol, feature_ids, NA_character_),
                      mapping_source = ifelse(is_symbol, "original_feature_id", "unmapped"),
                      mapping_status = ifelse(is_symbol, "direct_gene_symbol", "unmapped"),
                      is_ambiguous = FALSE,
                      notes = ifelse(is_symbol, "Feature ID treated as gene symbol by pattern.", "No platform annotation available."))
fwrite(mapping, file.path(table_dir, "gse73680_feature_gene_mapping.tsv"), sep = "\t")
mapped <- mapping[!is.na(gene_symbol)]
if (!nrow(mapped)) {
  fwrite(data.table(n_features_input = nrow(mat), n_features_mapped = 0L, n_features_unmapped = nrow(mat),
                    mapping_rate = 0, n_unique_genes = 0L, n_multi_probe_genes = 0L,
                    collapse_rule = "max_mean_expression_probe", status = "fail",
                    notes = "No features mapped to gene symbols; platform annotation required."),
         file.path(table_dir, "gse73680_gene_mapping_qc.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
expr <- merge(mapping[, .(feature_id, gene_symbol)], mat, by = "feature_id")
sample_cols <- setdiff(names(expr), c("feature_id", "gene_symbol"))
expr[, mean_expression := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_cols]
setorder(expr, gene_symbol, -mean_expression)
gene_expr <- expr[!duplicated(gene_symbol), c("gene_symbol", sample_cols), with = FALSE]
setnames(gene_expr, "gene_symbol", "gene")
fwrite(gene_expr, file.path(processed_dir, "gse73680_gene_expression_matrix.tsv.gz"), sep = "\t")
mapping_rate <- nrow(mapped) / nrow(mat)
n_unique <- uniqueN(mapped$gene_symbol)
status <- if (mapping_rate >= 0.60 && n_unique >= 8000) "pass" else if (mapping_rate >= 0.40 && n_unique >= 5000) "warning" else "fail"
fwrite(data.table(n_features_input = nrow(mat), n_features_mapped = nrow(mapped),
                  n_features_unmapped = nrow(mat) - nrow(mapped), mapping_rate = mapping_rate,
                  n_unique_genes = n_unique, n_multi_probe_genes = sum(duplicated(mapped$gene_symbol)),
                  collapse_rule = "max_mean_expression_probe", status = status,
                  notes = "Direct gene-symbol mapping attempted; use platform annotation if status is warning/fail."),
       file.path(table_dir, "gse73680_gene_mapping_qc.tsv"), sep = "\t")
message("wrote GSE73680 feature-to-gene mapping outputs")
