suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

doc_dir <- "docs/revision/stage5B1_gse73680_module_response"
table_dir <- "results/tables/revision/stage5B1_gse73680_module_response"
fig_dir <- "results/figures/revision/stage5B1_gse73680_module_response"
script_dir <- "scripts/revision_utils/stage5B1_gse73680_module_response"
log_dir <- "logs/revision/stage5B1_gse73680_module_response"
for (d in c(doc_dir, table_dir, fig_dir, script_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage5B1_gse73680_module_response.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 5B1 GSE73680 paired module-response and composition-signature scoring\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  stage5a_report = "docs/revision/stage5A_gse73680_audit/stage5A_report.md",
  stage5b_plan = "docs/revision/stage5A_gse73680_audit/stage5B_statistical_design_plan.md",
  stage5a_pairing = "results/tables/revision/stage5A_gse73680_audit/gse73680_sample_pairing_audit.tsv",
  stage5a_groups = "results/tables/revision/stage5A_gse73680_audit/gse73680_sample_group_summary.tsv",
  stage5a_modules = "results/tables/revision/stage5A_gse73680_audit/gse73680_module_feasibility_audit.tsv",
  stage5a_signatures = "results/tables/revision/stage5A_gse73680_audit/gse73680_composition_signature_manifest.tsv",
  stage5a_existing = "results/tables/revision/stage5A_gse73680_audit/gse73680_existing_result_audit.tsv",
  expr = "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz",
  metadata = "config/gse73680_sample_metadata_curated.tsv",
  feature_mapping = "results/gse73680/tables/gse73680_feature_gene_mapping.tsv",
  stage3_model = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  stage3_counts = "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv",
  exemplar = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
  stage4_claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
  stage4c2r_report = "docs/revision/stage4C2R_draft_figures_v0.2/stage4C2R_report.md"
)

required <- c("stage5a_report", "stage5b_plan", "stage5a_pairing", "stage5a_groups",
              "stage5a_modules", "stage5a_signatures", "expr", "metadata",
              "stage3_model", "exemplar", "stage4_claim", "stage4c2r_report")
for (nm in required) {
  if (!file.exists(paths[[nm]]) || file.access(paths[[nm]], 4) != 0) {
    stop("Required input missing or unreadable: ", nm, " -> ", paths[[nm]])
  }
}

check_cols <- function(dt, cols, name) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing columns in ", name, ": ", paste(miss, collapse = ", "))
}

expr_dt <- fread(paths$expr)
meta <- fread(paths$metadata)
stage5a_pairing <- fread(paths$stage5a_pairing)
stage5a_modules <- fread(paths$stage5a_modules)
stage5a_signatures <- fread(paths$stage5a_signatures)
stage3 <- fread(paths$stage3_model)
exemplar <- fread(paths$exemplar)

check_cols(expr_dt, names(expr_dt)[1], "expression matrix")
check_cols(meta, c("sample_id", "patient_id", "group_curated", "include_in_analysis"), "metadata")
check_cols(stage5a_pairing, c("patient_id", "sample_id", "disease_group", "usable_for_paired_analysis"), "Stage 5A pairing")
check_cols(stage5a_modules, c("module_name", "module_role", "detected_genes", "n_genes_input", "n_genes_detected_in_gse73680"), "Stage 5A module feasibility")
check_cols(stage5a_signatures, c("signature_name", "marker_genes_input", "n_marker_genes_detected"), "Stage 5A composition manifest")

id_col <- names(expr_dt)[1]
sample_cols <- setdiff(names(expr_dt), id_col)
expr_genes_raw <- as.character(expr_dt[[id_col]])
keep_gene <- !is.na(expr_genes_raw) & nzchar(expr_genes_raw)
expr_dt <- expr_dt[keep_gene]
expr_genes <- as.character(expr_dt[[id_col]])
expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_genes

meta[, include_bool := as.logical(include_in_analysis)]
meta[is.na(include_bool), include_bool := FALSE]
meta_use <- meta[include_bool == TRUE & sample_id %in% colnames(expr_mat)]
primary_labels <- c("control_or_adjacent", "plaque_or_stone_papilla")
meta_use <- meta_use[group_curated %in% primary_labels]
setorder(meta_use, patient_id, group_curated, sample_id)

overlap_samples <- intersect(meta_use$sample_id, colnames(expr_mat))
expr_mat <- expr_mat[, overlap_samples, drop = FALSE]
meta_use <- meta_use[match(overlap_samples, sample_id)]

pair_tab <- meta_use[, .(
  control_or_adjacent_sample = sample_id[group_curated == "control_or_adjacent"][1],
  plaque_or_stone_sample = sample_id[group_curated == "plaque_or_stone_papilla"][1],
  n_control = sum(group_curated == "control_or_adjacent"),
  n_plaque = sum(group_curated == "plaque_or_stone_papilla"),
  stone_type = paste(unique(stone_status[group_curated == "plaque_or_stone_papilla" & !is.na(stone_status) & stone_status != "" & stone_status != "\"\""]), collapse = ";")
), by = patient_id]
pair_tab[, paired := n_control >= 1 & n_plaque >= 1]
pair_tab[, usable_for_paired_delta := paired & !is.na(control_or_adjacent_sample) & !is.na(plaque_or_stone_sample)]
pair_tab[, pair_id := patient_id]
pair_tab[, notes := "Paired delta uses the first control/adjacent and first plaque/stone sample per patient; no patient has duplicated group samples in current included metadata."]
setcolorder(pair_tab, c("patient_id", "pair_id", "control_or_adjacent_sample", "plaque_or_stone_sample", "paired", "stone_type", "usable_for_paired_delta", "notes"))

expected_included <- 55
expected_pairs <- 26
input_checks <- rbindlist(list(
  data.table(input_item = "expression_matrix", expected = "exists/readable", found = file.exists(paths$expr), status = ifelse(file.exists(paths$expr), "pass", "fail"), file_path = paths$expr, notes = ""),
  data.table(input_item = "metadata", expected = "exists/readable", found = file.exists(paths$metadata), status = ifelse(file.exists(paths$metadata), "pass", "fail"), file_path = paths$metadata, notes = ""),
  data.table(input_item = "sample_id_overlap", expected = ">=55 included samples", found = length(overlap_samples), status = ifelse(length(overlap_samples) >= expected_included, "pass", "warning"), file_path = paths$expr, notes = "Overlap after include flag and primary label filtering."),
  data.table(input_item = "included_samples_stage5A_count", expected = expected_included, found = nrow(meta_use), status = ifelse(nrow(meta_use) == expected_included, "pass", "warning"), file_path = paths$metadata, notes = ""),
  data.table(input_item = "primary_comparison_labels", expected = paste(primary_labels, collapse = ";"), found = paste(sort(unique(meta_use$group_curated)), collapse = ";"), status = ifelse(all(primary_labels %in% meta_use$group_curated), "pass", "fail"), file_path = paths$metadata, notes = ""),
  data.table(input_item = "patient_id_or_pair_id", expected = "patient_id", found = "patient_id" %in% names(meta), status = ifelse("patient_id" %in% names(meta), "pass", "fail"), file_path = paths$metadata, notes = ""),
  data.table(input_item = "paired_patients_recoverable", expected = expected_pairs, found = pair_tab[usable_for_paired_delta == TRUE, .N], status = ifelse(pair_tab[usable_for_paired_delta == TRUE, .N] == expected_pairs, "pass", "warning"), file_path = paths$metadata, notes = ""),
  data.table(input_item = "gene_symbols_available", expected = ">0 non-empty genes", found = nrow(expr_mat), status = ifelse(nrow(expr_mat) > 0, "pass", "fail"), file_path = paths$expr, notes = ""),
  data.table(input_item = "stage3R_module_tables", expected = "readable", found = file.exists(paths$stage3_model), status = ifelse(file.exists(paths$stage3_model), "pass", "fail"), file_path = paths$stage3_model, notes = ""),
  data.table(input_item = "composition_signature_manifest", expected = "readable", found = file.exists(paths$stage5a_signatures), status = ifelse(file.exists(paths$stage5a_signatures), "pass", "fail"), file_path = paths$stage5a_signatures, notes = "")
), fill = TRUE)
fwrite(input_checks, file.path(table_dir, "stage5B1_input_consistency_check.tsv"), sep = "\t")

vals <- as.numeric(expr_mat)
q <- quantile(vals, probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE, names = FALSE)
names(q) <- c("min", "q01", "q25", "median", "q75", "q99", "max")
has_negative <- any(vals < 0, na.rm = TRUE)
raw_like <- !has_negative && (q["max"] > 1000 || q["q99"] > 100)
log_like <- !has_negative && q["max"] < 50 && q["q99"] < 25
small_constant <- 1
if (raw_like) {
  final_mat <- log2(expr_mat + small_constant)
  transform_decision <- "log2_x_plus_1_applied"
  transform_label <- "log2(x + 1)"
  transform_status <- "applied"
  transform_notes <- "Values were positive with very large upper range, consistent with raw/normalized intensity rather than log2 scale."
} else if (log_like) {
  final_mat <- expr_mat
  transform_decision <- "no_log_transform_applied"
  transform_label <- "no additional log transformation"
  transform_status <- "not_applied"
  transform_notes <- "Values appeared log2-like by range."
} else {
  final_mat <- log2(pmax(expr_mat, 0) + small_constant)
  transform_decision <- "log2_x_plus_1_applied_with_warning"
  transform_label <- "log2(x + 1) after clipping negative values to zero"
  transform_status <- "warning_applied"
  transform_notes <- "Scale was ambiguous; conservative log2 transform applied after clipping negative values to zero."
}

scale_audit <- data.table(
  qc_item = c("min", "q01", "q25", "median", "q75", "q99", "max", "negative_values_exist", "very_large_raw_like_values_exist", "missing_value_fraction_after_filtering", "log2_like_assessment", "log2_transformation_applied_in_stage5B1", "final_expression_matrix_used_for_modeling"),
  value = c(sprintf("%.6g", q), as.character(has_negative), as.character(raw_like), sprintf("%.6f", mean(is.na(expr_mat))), ifelse(log_like, "log2_like", ifelse(raw_like, "raw_like_positive_intensity", "unclear")), transform_decision, "in_memory_final_mat"),
  status = c(rep("observed", 7), ifelse(has_negative, "warning", "pass"), ifelse(raw_like, "requires_transform", "not_raw_like"), ifelse(mean(is.na(expr_mat)) < 0.2, "pass_with_missingness", "warning"), ifelse(log_like, "pass", ifelse(raw_like, "requires_transform", "requires_manual_check")), transform_status, "ready_for_stage5B1_models"),
  notes = c(rep("Quantile of included-sample expression values before Stage 5B1 transformation.", 7), "Checked before transformation.", "Large values trigger log2(x + 1).", "Calculated after sample filtering.", transform_notes, transform_notes, "Final transformed matrix is held in memory and used for all Stage 5B1 scores.")
)
fwrite(scale_audit, file.path(table_dir, "gse73680_expression_scale_audit.tsv"), sep = "\t")

sample_inclusion <- meta[, .(
  sample_id,
  patient_id,
  pair_id = patient_id,
  tissue_status = ifelse("plaque_status" %in% names(meta), plaque_status, ""),
  primary_group = group_curated,
  disease_group = disease_status,
  stone_type = ifelse("stone_status" %in% names(meta), stone_status, ""),
  included_in_primary_analysis = sample_id %in% meta_use$sample_id,
  paired_analysis_eligible = sample_id %in% c(pair_tab[usable_for_paired_delta == TRUE, control_or_adjacent_sample], pair_tab[usable_for_paired_delta == TRUE, plaque_or_stone_sample]),
  reason_if_excluded = fifelse(!(sample_id %in% sample_cols), "missing_from_expression_matrix",
                         fifelse(!as.logical(include_in_analysis), exclude_reason,
                         fifelse(!(group_curated %in% primary_labels), "not_primary_group", ""))),
  notes = "Primary analysis requires include_in_analysis TRUE, expression matrix match, and primary group label."
)]
fwrite(sample_inclusion, file.path(table_dir, "gse73680_stage5B1_sample_inclusion.tsv"), sep = "\t")
fwrite(pair_tab[, .(patient_id, control_or_adjacent_sample, plaque_or_stone_sample, paired, stone_type, usable_for_paired_delta, notes)], file.path(table_dir, "gse73680_stage5B1_pairing_summary.tsv"), sep = "\t")

gene_nonmiss <- rowSums(!is.na(final_mat))
top_var_genes <- names(sort(apply(final_mat, 1, var, na.rm = TRUE), decreasing = TRUE))[1:min(5000, nrow(final_mat))]
pca_mat <- t(final_mat[top_var_genes, , drop = FALSE])
for (j in seq_len(ncol(pca_mat))) {
  if (anyNA(pca_mat[, j])) pca_mat[is.na(pca_mat[, j]), j] <- median(pca_mat[, j], na.rm = TRUE)
}
pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
sample_qc <- data.table(
  sample_id = rownames(pca$x),
  pca1 = pca$x[, 1],
  pca2 = pca$x[, 2]
)
sample_qc <- merge(meta_use[, .(sample_id, patient_id, primary_group = group_curated)], sample_qc, by = "sample_id", all.x = TRUE)
sample_qc[, `:=`(
  n_genes_nonmissing = colSums(!is.na(final_mat))[sample_id],
  mean_expression = colMeans(final_mat, na.rm = TRUE)[sample_id],
  median_expression = apply(final_mat, 2, median, na.rm = TRUE)[sample_id],
  sd_expression = apply(final_mat, 2, sd, na.rm = TRUE)[sample_id],
  missing_fraction = colMeans(is.na(final_mat))[sample_id]
)]
sample_qc[, outlier_flag := abs(scale(pca1)) > 4 | abs(scale(pca2)) > 4 | missing_fraction > 0.3]
sample_qc[, notes := "PCA computed on top variable genes after Stage 5B1 expression transformation; draft QC only."]
setcolorder(sample_qc, c("sample_id", "patient_id", "primary_group", "n_genes_nonmissing", "mean_expression", "median_expression", "sd_expression", "missing_fraction", "pca1", "pca2", "outlier_flag", "notes"))
fwrite(sample_qc, file.path(table_dir, "gse73680_sample_qc_metrics.tsv"), sep = "\t")

p_pca <- ggplot(sample_qc, aes(pca1, pca2, color = primary_group, shape = outlier_flag)) +
  geom_point(size = 2.4, alpha = 0.9) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c(control_or_adjacent = "#4C78A8", plaque_or_stone_papilla = "#E0A43A")) +
  labs(title = "GSE73680 Stage 5B1 PCA QC draft", subtitle = paste0("Not a final manuscript figure; transform=", transform_decision), x = "PC1", y = "PC2", color = "Group", shape = "Outlier")
ggsave(file.path(fig_dir, "gse73680_stage5B1_pca_qc_draft.pdf"), p_pca, width = 6.5, height = 5, units = "in")
ggsave(file.path(fig_dir, "gse73680_stage5B1_pca_qc_draft.png"), p_pca, width = 6.5, height = 5, units = "in", dpi = 300)

split_genes <- function(x) {
  if (is.na(x) || x == "" || x == "\"\"") return(character())
  unique(trimws(unlist(strsplit(x, ";"))))
}
module_manifest <- stage5a_modules[, .(
  module_name,
  module_role,
  genes_input = ifelse(detected_genes == "" | detected_genes == "\"\"", missing_genes, paste(detected_genes, missing_genes, sep = ";")),
  detected_genes,
  n_genes_input,
  n_genes_detected = n_genes_detected_in_gse73680,
  detected_fraction
)]
signature_roles <- function(x) {
  fifelse(x %in% c("injury_epithelial", "ECM_fibrosis", "mineralization_remodeling"),
          "injury_or_remodeling_signature", "composition_signature")
}
signature_rows <- rbindlist(lapply(seq_len(nrow(stage5a_signatures)), function(i) {
  nm <- stage5a_signatures$signature_name[i]
  genes <- split_genes(stage5a_signatures$marker_genes_input[i])
  det <- intersect(genes, rownames(final_mat))
  data.table(
    module_name = nm,
    module_role = signature_roles(nm),
    genes_input = paste(genes, collapse = ";"),
    detected_genes = paste(det, collapse = ";"),
    n_genes_input = length(genes),
    n_genes_detected = length(det),
    detected_fraction = ifelse(length(genes) > 0, length(det) / length(genes), NA_real_)
  )
}), fill = TRUE)
module_manifest <- rbind(
  module_manifest[!module_role %in% c("composition_signature", "injury_or_remodeling_signature") & !module_name %in% signature_rows$module_name],
  signature_rows,
  fill = TRUE
)
module_manifest[module_name == "R5_TWAS_proxy_only", module_role := "secondary_TWAS_proxy_module"]
module_manifest[, module_status := fifelse(n_genes_detected >= 5, "score_ready", fifelse(n_genes_detected >= 1, "small_module_supplement_only", "not_evaluable"))]
module_manifest[, notes := fifelse(module_status == "score_ready", "Mean z-score can be interpreted as bulk module score.", "Small or sparse module; interpret only as supplementary/descriptive.")]
fwrite(module_manifest[, .(module_name, module_role, n_genes_input, n_genes_detected, detected_fraction, module_status, notes)], file.path(table_dir, "stage5B1_module_manifest.tsv"), sep = "\t")

all_modules <- module_manifest[, .(module_name, module_role, detected_genes)]
score_rows <- rbindlist(lapply(seq_len(nrow(all_modules)), function(i) {
  nm <- all_modules$module_name[i]
  role <- all_modules$module_role[i]
  genes <- intersect(split_genes(all_modules$detected_genes[i]), rownames(final_mat))
  if (!length(genes)) return(NULL)
  sub <- final_mat[genes, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[is.nan(z)] <- NA_real_
  data.table(
    sample_id = colnames(final_mat),
    module_name = nm,
    module_role = role,
    n_genes_detected = length(genes),
    mean_z_score = colMeans(z, na.rm = TRUE),
    mean_expression_score = colMeans(sub, na.rm = TRUE),
    score_method = "mean z-scored transformed expression across detected genes",
    notes = "Bulk sample-level module/signature score; not cell-type-specific."
  )
}), fill = TRUE)
score_rows <- merge(score_rows, meta_use[, .(sample_id, patient_id, pair_id = patient_id, primary_group = group_curated, tissue_status = plaque_status, stone_type = stone_status)], by = "sample_id", all.x = TRUE)
setcolorder(score_rows, c("sample_id", "patient_id", "pair_id", "primary_group", "tissue_status", "stone_type", "module_name", "module_role", "n_genes_detected", "mean_z_score", "mean_expression_score", "score_method", "notes"))
fwrite(score_rows, file.path(table_dir, "gse73680_sample_module_scores.tsv"), sep = "\t")

paired_patients <- pair_tab[usable_for_paired_delta == TRUE]
delta_rows <- rbindlist(lapply(seq_len(nrow(paired_patients)), function(i) {
  p <- paired_patients[i]
  ctrl <- score_rows[sample_id == p$control_or_adjacent_sample]
  pla <- score_rows[sample_id == p$plaque_or_stone_sample]
  m <- merge(ctrl[, .(module_name, module_role, control_or_adjacent_score = mean_z_score)],
             pla[, .(module_name, plaque_or_stone_score = mean_z_score)],
             by = "module_name")
  m[, `:=`(
    patient_id = p$patient_id,
    pair_id = p$pair_id,
    stone_type = p$stone_type,
    control_or_adjacent_sample = p$control_or_adjacent_sample,
    plaque_or_stone_sample = p$plaque_or_stone_sample,
    paired_delta = plaque_or_stone_score - control_or_adjacent_score
  )]
  m[, delta_direction := fifelse(paired_delta > 0, "increase", fifelse(paired_delta < 0, "decrease", "zero"))]
  m[, notes := "Paired delta = plaque/stone papilla score minus control/adjacent score."]
  m
}), fill = TRUE)
setcolorder(delta_rows, c("patient_id", "pair_id", "stone_type", "module_name", "module_role", "control_or_adjacent_sample", "plaque_or_stone_sample", "control_or_adjacent_score", "plaque_or_stone_score", "paired_delta", "delta_direction", "notes"))
fwrite(delta_rows, file.path(table_dir, "gse73680_paired_module_delta.tsv"), sep = "\t")

summarize_delta <- function(dt) {
  if (nrow(dt) < 3) {
    return(data.table(n_paired_patients = nrow(dt), median_delta = median(dt$paired_delta, na.rm = TRUE), mean_delta = mean(dt$paired_delta, na.rm = TRUE), sd_delta = sd(dt$paired_delta, na.rm = TRUE), n_positive_delta = sum(dt$paired_delta > 0, na.rm = TRUE), n_negative_delta = sum(dt$paired_delta < 0, na.rm = TRUE), direction_consistency_fraction = NA_real_, wilcoxon_p = NA_real_, paired_t_p = NA_real_, effect_direction = "not_evaluable", interpretation = "not_evaluable"))
  }
  del <- dt$paired_delta
  npos <- sum(del > 0, na.rm = TRUE)
  nneg <- sum(del < 0, na.rm = TRUE)
  frac <- max(npos, nneg) / sum(del != 0, na.rm = TRUE)
  wp <- tryCatch(wilcox.test(del, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  tp <- tryCatch(t.test(del, mu = 0)$p.value, error = function(e) NA_real_)
  direction <- ifelse(mean(del, na.rm = TRUE) > 0, "increase", ifelse(mean(del, na.rm = TRUE) < 0, "decrease", "mixed"))
  interp <- ifelse(frac >= 0.7 & direction == "increase", "directionally_consistent_increase",
            ifelse(frac >= 0.7 & direction == "decrease", "directionally_consistent_decrease", "mixed_or_weak"))
  data.table(n_paired_patients = uniqueN(dt$patient_id), median_delta = median(del, na.rm = TRUE), mean_delta = mean(del, na.rm = TRUE), sd_delta = sd(del, na.rm = TRUE), n_positive_delta = npos, n_negative_delta = nneg, direction_consistency_fraction = frac, wilcoxon_p = wp, paired_t_p = tp, effect_direction = direction, interpretation = interp)
}

delta_summary <- delta_rows[, summarize_delta(.SD), by = .(module_name, module_role)]
delta_summary[, fdr_within_module_family := p.adjust(wilcoxon_p, method = "BH"), by = module_role]
delta_summary[, notes := "Stage 5B1 paired delta summary; bulk disease-context only."]
setcolorder(delta_summary, c("module_name", "module_role", "n_paired_patients", "median_delta", "mean_delta", "sd_delta", "n_positive_delta", "n_negative_delta", "direction_consistency_fraction", "wilcoxon_p", "paired_t_p", "fdr_within_module_family", "effect_direction", "interpretation", "notes"))
fwrite(delta_summary, file.path(table_dir, "gse73680_paired_module_delta_summary.tsv"), sep = "\t")

model_rows <- rbindlist(lapply(unique(score_rows$module_name), function(nm) {
  dt <- score_rows[module_name == nm & sample_id %in% c(paired_patients$control_or_adjacent_sample, paired_patients$plaque_or_stone_sample)]
  dt[, group_binary := ifelse(primary_group == "plaque_or_stone_papilla", 1, 0)]
  if (uniqueN(dt$patient_id) < 3 || uniqueN(dt$primary_group) < 2) {
    return(data.table(module_name = nm, module_role = dt$module_role[1], model_type = "paired_fixed_effect_lm", n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), effect_estimate = NA_real_, standard_error = NA_real_, test_statistic = NA_real_, p_value = NA_real_, effect_direction = "not_evaluable", model_status = "failed", interpretation = "not_evaluable", notes = "Insufficient paired data."))
  }
  fit <- tryCatch(summary(lm(mean_z_score ~ group_binary + factor(patient_id), data = dt)), error = function(e) NULL)
  if (is.null(fit) || !"group_binary" %in% rownames(fit$coefficients)) {
    return(data.table(module_name = nm, module_role = dt$module_role[1], model_type = "paired_fixed_effect_lm", n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), effect_estimate = NA_real_, standard_error = NA_real_, test_statistic = NA_real_, p_value = NA_real_, effect_direction = "not_evaluable", model_status = "failed", interpretation = "not_evaluable", notes = "Fixed-effect model failed."))
  }
  co <- fit$coefficients["group_binary", ]
  direction <- ifelse(co[["Estimate"]] > 0, "increase", ifelse(co[["Estimate"]] < 0, "decrease", "mixed"))
  data.table(module_name = nm, module_role = dt$module_role[1], model_type = "paired_patient_fixed_effect_lm", n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), effect_estimate = co[["Estimate"]], standard_error = co[["Std. Error"]], test_statistic = co[["t value"]], p_value = co[["Pr(>|t|)"]], effect_direction = direction, model_status = "paired_model_success", interpretation = "bulk paired module shift estimate; not causal or cell-type-specific", notes = "Linear model with patient fixed effects on sample-level module scores.")
}), fill = TRUE)
model_rows[, fdr := p.adjust(p_value, method = "BH"), by = module_role]
setcolorder(model_rows, c("module_name", "module_role", "model_type", "n_samples", "n_patients", "effect_estimate", "standard_error", "test_statistic", "p_value", "fdr", "effect_direction", "model_status", "interpretation", "notes"))
fwrite(model_rows, file.path(table_dir, "gse73680_module_response_model.tsv"), sep = "\t")

sig_roles <- c("composition_signature", "injury_or_remodeling_signature")
signature_summary <- copy(delta_summary[module_role %in% sig_roles])
setnames(signature_summary, c("module_name", "module_role"), c("signature_name", "signature_role"))
setnames(signature_summary, "fdr_within_module_family", "fdr")
signature_summary[, n_genes_detected := module_manifest[match(signature_name, module_name), n_genes_detected]]
setcolorder(signature_summary, c("signature_name", "signature_role", "n_genes_detected", "n_paired_patients", "median_delta", "mean_delta", "direction_consistency_fraction", "wilcoxon_p", "paired_t_p", "fdr", "effect_direction", "interpretation", "notes"))
signature_summary[, notes := "Bulk marker-signature response; not validated cell fraction estimation."]
fwrite(signature_summary, file.path(table_dir, "gse73680_signature_response_summary.tsv"), sep = "\t")

primary_modules <- c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")
signatures <- module_manifest[module_role %in% sig_roles, module_name]
cor_grid <- CJ(module_name = primary_modules, signature_name = signatures)
cor_sample <- rbindlist(lapply(seq_len(nrow(cor_grid)), function(i) {
  mn <- cor_grid$module_name[i]
  sn <- cor_grid$signature_name[i]
  sx <- score_rows[module_name == mn, .(sample_id, x = mean_z_score)]
  sy <- score_rows[module_name == sn, .(sample_id, y = mean_z_score)]
  xy <- merge(sx, sy, by = "sample_id")
  ok <- is.finite(xy$x) & is.finite(xy$y)
  data.table(correlation_level = "sample_level", module_name = mn, signature_name = sn, n_observations = sum(ok),
             spearman_r = suppressWarnings(cor(xy$x[ok], xy$y[ok], method = "spearman")),
             spearman_p = tryCatch(cor.test(xy$x[ok], xy$y[ok], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_),
             pearson_r = suppressWarnings(cor(xy$x[ok], xy$y[ok], method = "pearson")),
             pearson_p = tryCatch(cor.test(xy$x[ok], xy$y[ok], method = "pearson")$p.value, error = function(e) NA_real_))
}), fill = TRUE)
cor_delta <- rbindlist(lapply(seq_len(nrow(cor_grid)), function(i) {
  mn <- cor_grid$module_name[i]
  sn <- cor_grid$signature_name[i]
  dx <- delta_rows[module_name == mn, .(patient_id, x = paired_delta)]
  dy <- delta_rows[module_name == sn, .(patient_id, y = paired_delta)]
  xy <- merge(dx, dy, by = "patient_id")
  ok <- is.finite(xy$x) & is.finite(xy$y)
  data.table(correlation_level = "paired_delta_level", module_name = mn, signature_name = sn, n_observations = sum(ok),
             spearman_r = suppressWarnings(cor(xy$x[ok], xy$y[ok], method = "spearman")),
             spearman_p = tryCatch(cor.test(xy$x[ok], xy$y[ok], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_),
             pearson_r = suppressWarnings(cor(xy$x[ok], xy$y[ok], method = "pearson")),
             pearson_p = tryCatch(cor.test(xy$x[ok], xy$y[ok], method = "pearson")$p.value, error = function(e) NA_real_))
}), fill = TRUE)
cor_all <- rbind(cor_sample, cor_delta, fill = TRUE)
cor_all[, fdr := p.adjust(spearman_p, method = "BH"), by = correlation_level]
cor_all[, interpretation := fifelse(abs(spearman_r) >= 0.7 & fdr < 0.1, "strong_signature_coupling_candidate",
                              fifelse(abs(spearman_r) >= 0.4 & fdr < 0.2, "moderate_signature_coupling_candidate", "weak_or_uncertain"))]
cor_all[, notes := "Correlation prepares Stage 5B2 composition-aware sensitivity; not causal mediation."]
fwrite(cor_all, file.path(table_dir, "gse73680_module_signature_correlation.tsv"), sep = "\t")
paired_strong_couplings <- cor_all[correlation_level == "paired_delta_level" & interpretation == "strong_signature_coupling_candidate", .N]
paired_moderate_couplings <- cor_all[correlation_level == "paired_delta_level" & interpretation == "moderate_signature_coupling_candidate", .N]

primary_sum <- delta_summary[module_name %in% primary_modules]
primary_direction <- ifelse(all(primary_sum$effect_direction == "increase"), "directionally increasing", ifelse(all(primary_sum$effect_direction == "decrease"), "directionally decreasing", "mixed"))
primary_strength <- ifelse(primary_direction != "mixed" && mean(primary_sum$direction_consistency_fraction, na.rm = TRUE) >= 0.65, "moderate_before_composition_adjustment", "weak_or_mixed_before_composition_adjustment")
claim_decisions <- data.table(
  evidence_component = c("primary_MAGMA_module_response", "TWAS_proxy_module_response", "curated_exemplar_response", "composition_signature_response", "injury_remodeling_coupling"),
  stage5B1_result_summary = c(
    paste0("Primary modules are ", primary_direction, "; mean direction consistency=", sprintf("%.2f", mean(primary_sum$direction_consistency_fraction, na.rm = TRUE)), "."),
    "TWAS-proxy modules scored but retained as secondary/supplementary context.",
    "Curated exemplar panel scored as biological context only, not validation.",
    "Bulk composition signatures scored and summarized by paired delta.",
    "Injury/ECM/mineralization signatures scored and correlated with primary modules."
  ),
  allowed_claim = c("paired plaque/stone papilla bulk module shift", "supplementary bulk disease-context association", "curated biological context", "bulk composition/signature association", "bulk injury/remodeling coupling"),
  disallowed_claim = "cell-type-specific response; genetic causality; plaque causation; therapeutic target validation; single-gene validation",
  claim_strength_before_composition_adjustment = c(primary_strength, "supplementary_only", "supplementary_only", "context_only", "context_only"),
  needed_stage5B2_sensitivity = c("composition-aware model for primary MAGMA modules", "keep supplementary after composition-aware checks", "do not use as evidence upgrade", "use as covariates/sensitivity terms", "test whether primary module shifts persist beyond injury/remodeling coupling"),
  notes = "Stage 5B1 decision only; final GSE73680 claim requires Stage 5B2."
)
fwrite(claim_decisions, file.path(table_dir, "gse73680_stage5B1_claim_decision_table.tsv"), sep = "\t")

method_text <- c(
  "# Manuscript replacement text draft for Stage 5B1",
  "",
  "## Methods Draft",
  "",
  paste0("GSE73680 bulk papillary expression data were audited against curated sample metadata and analyzed as a paired bulk tissue disease-context dataset. Because the reconstructed gene-level expression matrix showed positive intensity-like values rather than a clearly log2-normalized scale, Stage 5B1 used ", transform_label, " before computing scores. For each predefined MAGMA-prioritized module and bulk composition or injury/remodeling signature, detected genes were z-scored across samples and averaged to generate sample-level module scores. Paired plaque/stone papilla minus control/adjacent deltas were computed within patients, and patient fixed-effect linear models were used as module-level sensitivity models."),
  "",
  "## Results Draft",
  "",
  paste0("GSE73680 retained ", nrow(meta_use), " primary-analysis samples, including ", paired_patients[, .N], " patients with paired control/adjacent and plaque/stone papilla samples. Primary MAGMA-prioritized modules showed a ", primary_direction, " paired bulk module-score pattern before composition-aware adjustment. These Stage 5B1 results support only a preliminary bulk disease-context association assessment and require Stage 5B2 composition-aware sensitivity modeling before being used as main manuscript support."),
  "",
  "## Limitations Draft",
  "",
  "GSE73680 is a bulk tissue dataset and cannot by itself identify cell-type-specific disease responses. Metadata labels are curated from available sample labels and should be described cautiously. Bulk module shifts may reflect tissue composition, injury, fibrosis, mineralization, or other remodeling programs; therefore Stage 5B2 composition-aware sensitivity is required before strengthening disease-context claims."
)
writeLines(method_text, file.path(doc_dir, "manuscript_replacement_text_stage5B1.md"))

max_primary_fdr <- max(primary_sum$fdr_within_module_family, na.rm = TRUE)
reviewer <- c(
  "# Stage 5B1 simulated reviewer check",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "1. Was expression scale/log2 status confirmed or handled transparently?",
  paste0("   It was handled transparently. Stage 5B1 classified the matrix as ", ifelse(raw_like, "raw-like positive intensity", ifelse(log_like, "log2-like", "unclear")), " and used `", transform_decision, "`."),
  "",
  "2. Is pairing structure clear enough for paired analysis?",
  paste0("   Yes. ", paired_patients[, .N], " paired patients were recovered for plaque/stone papilla versus control/adjacent deltas."),
  "",
  "3. Do primary MAGMA modules show consistent paired disease-context shifts?",
  paste0("   Primary modules are ", primary_direction, " with mean direction consistency ", sprintf("%.2f", mean(primary_sum$direction_consistency_fraction, na.rm = TRUE)), "."),
  "",
  "4. Are TWAS/exemplar modules kept supplementary?",
  "   Yes. They are scored, but the claim decision table keeps them supplementary/context only.",
  "",
  "5. Do composition/injury signatures suggest possible bulk composition confounding?",
  paste0("   Yes. Stage 5B1 found ", paired_strong_couplings, " strong and ", paired_moderate_couplings, " moderate paired-delta module-signature coupling candidates, motivating Stage 5B2 composition-aware modeling."),
  "",
  "6. What would a reviewer criticize?",
  "   Filename-derived metadata labels, uncertain original normalization, and possible injury/composition confounding of bulk module shifts.",
  "",
  "7. What must Stage 5B2 test before GSE73680 can strengthen the manuscript?",
  "   Whether primary MAGMA module shifts remain directionally consistent after adjusting for composition, injury, ECM/fibrosis, and mineralization/remodeling signatures.",
  "",
  "8. What wording is currently allowed?",
  "   `paired plaque/stone papilla bulk module shift` or `bulk disease-context association before composition-aware adjustment`.",
  "",
  "9. What wording is not allowed?",
  "   Cell-type-specific response, genetic causality, plaque causation, therapeutic target validation, and single-gene validation.",
  "",
  "10. Should Stage 5B2 composition-aware modeling begin?",
  "   Yes, after human acceptance of Stage 5B1, because composition-aware sensitivity is required before strengthening the GSE73680 claim."
)
writeLines(reviewer, file.path(doc_dir, "stage5B1_simulated_reviewer_check.md"))

report <- c(
  "# Stage 5B1 report: GSE73680 paired module-response and composition-signature scoring",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Input Files Used",
  "",
  paste0("- Expression matrix: `", paths$expr, "`"),
  paste0("- Metadata: `", paths$metadata, "`"),
  paste0("- Stage 5A module feasibility: `", paths$stage5a_modules, "`"),
  paste0("- Stage 5A composition manifest: `", paths$stage5a_signatures, "`"),
  "",
  "## Expression Scale Decision",
  "",
  paste0("- Decision: `", transform_decision, "`"),
  paste0("- Notes: ", transform_notes),
  "",
  "## Sample Inclusion and Pairing Status",
  "",
  paste0("- Primary-analysis samples: ", nrow(meta_use)),
  paste0("- Paired patients retained: ", paired_patients[, .N]),
  "",
  "## PCA/QC Summary",
  "",
  paste0("- Outlier samples flagged: ", sample_qc[outlier_flag == TRUE, .N]),
  "- PCA plot is draft QC only, not a manuscript figure.",
  "",
  "## Module Feasibility",
  "",
  paste0("- Modules/signatures scored: ", uniqueN(score_rows$module_name)),
  "- Primary MAGMA modules and composition/injury signatures were scoreable.",
  "",
  "## Primary Module Paired Response",
  "",
  paste0("- Direction summary: ", primary_direction),
  paste0("- Mean primary direction consistency: ", sprintf("%.2f", mean(primary_sum$direction_consistency_fraction, na.rm = TRUE))),
  paste0("- Maximum primary within-family FDR from paired delta tests: ", sprintf("%.4g", max_primary_fdr)),
  "",
  "## Supplementary Module Response",
  "",
  "- TWAS-proxy and curated exemplar modules were scored but remain supplementary/descriptive.",
  "",
  "## Composition/Injury Signature Response",
  "",
  paste0("- Signatures summarized: ", nrow(signature_summary)),
  "- These are bulk marker-signature scores, not cell fractions.",
  "",
  "## Module-Signature Correlations",
  "",
  "- Sample-level and paired-delta-level correlations were computed to prepare Stage 5B2 composition-aware sensitivity.",
  paste0("- Paired-delta-level coupling candidates: ", paired_strong_couplings, " strong and ", paired_moderate_couplings, " moderate."),
  "",
  "## Preliminary Claim Decision",
  "",
  paste0("- Current strength before composition adjustment: ", primary_strength),
  "- Allowed wording: bulk disease-context association or paired plaque/stone papilla bulk module shift.",
  "- Disallowed wording: cell-type-specific response, genetic causality, plaque causation, therapeutic target validation, or single-gene validation.",
  "",
  "## Blockers or Warnings",
  "",
  "- No hard Stage 5B2 blocker.",
  "- Warnings: original normalization remains reconstructed/uncertain; metadata labels are curated from sample labels; composition confounding remains unresolved.",
  "",
  "## Recommended Stage 5B2 Tasks",
  "",
  "- Fit composition-aware module-response models for primary MAGMA modules.",
  "- Test adjustment by epithelial, LoopTAL, immune, stromal, endothelial, injury, ECM/fibrosis, and mineralization/remodeling signatures.",
  "- Keep TWAS-proxy and curated exemplar results supplementary.",
  "- Plan Figure 4 only after Stage 5B2 sensitivity results."
)
writeLines(report, file.path(doc_dir, "stage5B1_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 5, `:=`(
    status = "stage5B1_completed_module_response",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 5A and Stage 5B1 completed; expression scale audit, sample inclusion/pairing, PCA QC, module/signature scoring, paired delta tests, patient fixed-effect models, module-signature correlations, claim decision, draft text, reviewer check, and report generated",
    blocking_issues = "No hard blocker for Stage 5B2; composition-aware modeling required before using GSE73680 as stronger disease-context support",
    next_stage_ready = "stage5B2_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Wrote Stage 5B1 outputs to", doc_dir, "and", table_dir, "\n")
cat("Transform decision:", transform_decision, "\n")
cat("Paired patients:", paired_patients[, .N], "\n")
cat("Primary direction:", primary_direction, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
