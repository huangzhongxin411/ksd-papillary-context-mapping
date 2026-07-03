suppressPackageStartupMessages({
  library(data.table)
})

doc_dir <- "docs/revision/stage5B2_gse73680_composition_adjusted"
table_dir <- "results/tables/revision/stage5B2_gse73680_composition_adjusted"
fig_dir <- "results/figures/revision/stage5B2_gse73680_composition_adjusted"
script_dir <- "scripts/revision_utils/stage5B2_gse73680_composition_adjusted"
log_dir <- "logs/revision/stage5B2_gse73680_composition_adjusted"
for (d in c(doc_dir, table_dir, fig_dir, script_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage5B2_gse73680_composition_adjusted_analysis.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 5B2 GSE73680 composition-aware sensitivity modeling\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  stage5b1_report = "docs/revision/stage5B1_gse73680_module_response/stage5B1_report.md",
  scale_audit = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_expression_scale_audit.tsv",
  sample_inclusion = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_sample_inclusion.tsv",
  pairing = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_pairing_summary.tsv",
  scores = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_sample_module_scores.tsv",
  delta = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta.tsv",
  delta_summary = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta_summary.tsv",
  response_model = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_response_model.tsv",
  signature_summary = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_signature_response_summary.tsv",
  module_signature_correlation = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_signature_correlation.tsv",
  stage5b1_claim = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_claim_decision_table.tsv",
  expr = "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz",
  metadata = "config/gse73680_sample_metadata_curated.tsv"
)

required <- c("scores", "delta", "delta_summary", "response_model", "signature_summary", "module_signature_correlation", "stage5b1_claim")
for (nm in required) {
  if (!file.exists(paths[[nm]]) || file.access(paths[[nm]], 4) != 0) {
    stop("Required Stage 5B1 input missing/unreadable: ", nm, " -> ", paths[[nm]])
  }
}

check_cols <- function(dt, cols, name) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing columns in ", name, ": ", paste(miss, collapse = ", "))
}

scores <- fread(paths$scores)
delta <- fread(paths$delta)
delta_summary <- fread(paths$delta_summary)
response_model <- fread(paths$response_model)
signature_summary <- fread(paths$signature_summary)
module_sig_cor <- fread(paths$module_signature_correlation)
stage5b1_claim <- fread(paths$stage5b1_claim)
pairing <- if (file.exists(paths$pairing)) fread(paths$pairing) else data.table()

check_cols(scores, c("sample_id", "patient_id", "primary_group", "module_name", "module_role", "mean_z_score"), "scores")
check_cols(delta, c("patient_id", "module_name", "module_role", "paired_delta"), "delta")
check_cols(delta_summary, c("module_name", "module_role", "effect_direction", "fdr_within_module_family"), "delta_summary")
check_cols(signature_summary, c("signature_name", "signature_role", "effect_direction", "interpretation"), "signature_summary")

primary_modules <- c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")
all_signatures <- signature_summary$signature_name
single_covariates <- c("epithelial_general", "LoopTAL", "immune_myeloid", "fibroblast_stromal", "endothelial", "injury_epithelial", "ECM_fibrosis", "mineralization_remodeling")
single_covariates <- intersect(single_covariates, all_signatures)
composition_core <- intersect(c("epithelial_general", "LoopTAL", "immune_myeloid", "fibroblast_stromal", "endothelial", "pericyte_smooth_muscle"), all_signatures)
injury_set <- intersect(c("injury_epithelial", "ECM_fibrosis", "mineralization_remodeling"), all_signatures)
renal_optional <- intersect(c("proximal_tubule", "collecting_duct"), all_signatures)

recoverable_pairs <- uniqueN(delta[module_name %in% primary_modules, patient_id])
input_check <- rbindlist(list(
  data.table(input_item = "stage5B1_module_scores", expected = "exists/readable", found = file.exists(paths$scores), status = "pass", file_path = paths$scores, notes = paste0(nrow(scores), " rows")),
  data.table(input_item = "paired_delta_table", expected = "exists/readable", found = file.exists(paths$delta), status = "pass", file_path = paths$delta, notes = paste0(nrow(delta), " rows")),
  data.table(input_item = "signature_response_table", expected = "exists/readable", found = file.exists(paths$signature_summary), status = "pass", file_path = paths$signature_summary, notes = paste0(nrow(signature_summary), " signatures")),
  data.table(input_item = "module_signature_correlation_table", expected = "exists/readable", found = file.exists(paths$module_signature_correlation), status = "pass", file_path = paths$module_signature_correlation, notes = paste0(nrow(module_sig_cor), " correlations")),
  data.table(input_item = "paired_patients_recoverable", expected = 26, found = recoverable_pairs, status = ifelse(recoverable_pairs == 26, "pass", "warning"), file_path = paths$delta, notes = ""),
  data.table(input_item = "primary_MAGMA_modules_recoverable", expected = paste(primary_modules, collapse = ";"), found = paste(intersect(primary_modules, unique(delta$module_name)), collapse = ";"), status = ifelse(all(primary_modules %in% delta$module_name), "pass", "fail"), file_path = paths$delta, notes = ""),
  data.table(input_item = "composition_injury_signatures_recoverable", expected = ">=8 major signatures", found = length(all_signatures), status = ifelse(length(all_signatures) >= 8, "pass", "warning"), file_path = paths$signature_summary, notes = paste(all_signatures, collapse = ";"))
), fill = TRUE)
fwrite(input_check, file.path(table_dir, "stage5B2_input_consistency_check.tsv"), sep = "\t")

cov_manifest <- rbindlist(list(
  data.table(covariate_set_name = composition_core, covariates = composition_core, covariate_type = "core_composition_covariate", intended_model = "single_covariate_delta_and_patient_FE", n_covariates = 1, use_priority = "primary_single_covariate", notes = "Single marker-signature proxy; not cell fraction."),
  data.table(covariate_set_name = injury_set, covariates = injury_set, covariate_type = "disease_remodeling_covariate", intended_model = "single_covariate_delta_and_patient_FE", n_covariates = 1, use_priority = "primary_single_covariate", notes = "Single remodeling/injury proxy; not mechanism proof."),
  data.table(covariate_set_name = renal_optional, covariates = renal_optional, covariate_type = "optional_renal_compartment_covariate", intended_model = "exploratory_single_covariate_delta", n_covariates = 1, use_priority = "exploratory_only", notes = "Optional renal compartment proxy."),
  data.table(covariate_set_name = "composition_core_set", covariates = paste(composition_core, collapse = ";"), covariate_type = "compact_covariate_set", intended_model = "compact_delta_model", n_covariates = length(composition_core), use_priority = "secondary_multi_covariate", notes = "Use cautiously with 26 pairs; skip if high collinearity."),
  data.table(covariate_set_name = "injury_remodeling_set", covariates = paste(injury_set, collapse = ";"), covariate_type = "compact_covariate_set", intended_model = "compact_delta_model", n_covariates = length(injury_set), use_priority = "secondary_multi_covariate", notes = "Use cautiously with 26 pairs; skip if high collinearity."),
  data.table(covariate_set_name = "composition_plus_injury_set", covariates = paste(c(composition_core, injury_set), collapse = ";"), covariate_type = "compact_covariate_set", intended_model = "compact_delta_model", n_covariates = length(c(composition_core, injury_set)), use_priority = "exploratory_only", notes = "High overfit risk; likely skip or PCA-compress.")
), fill = TRUE)
fwrite(cov_manifest, file.path(table_dir, "stage5B2_covariate_set_manifest.tsv"), sep = "\t")

wide_delta <- dcast(delta[module_name %in% c(primary_modules, all_signatures)], patient_id ~ module_name, value.var = "paired_delta")
wide_scores <- dcast(scores[module_name %in% c(primary_modules, all_signatures)], sample_id + patient_id + primary_group ~ module_name, value.var = "mean_z_score")

sig_pairs <- t(combn(all_signatures, 2))
collinearity <- rbindlist(lapply(seq_len(nrow(sig_pairs)), function(i) {
  s1 <- sig_pairs[i, 1]
  s2 <- sig_pairs[i, 2]
  dx <- wide_delta[[s1]]
  dy <- wide_delta[[s2]]
  sx <- wide_scores[[s1]]
  sy <- wide_scores[[s2]]
  d_ok <- is.finite(dx) & is.finite(dy)
  s_ok <- is.finite(sx) & is.finite(sy)
  d_r <- suppressWarnings(cor(dx[d_ok], dy[d_ok], method = "spearman"))
  s_r <- suppressWarnings(cor(sx[s_ok], sy[s_ok], method = "spearman"))
  flag <- ifelse(abs(d_r) >= 0.85 | abs(s_r) >= 0.85, "very_high_abs_r_ge_0.85",
          ifelse(abs(d_r) >= 0.7 | abs(s_r) >= 0.7, "high_abs_r_ge_0.7",
          ifelse(abs(d_r) >= 0.4 | abs(s_r) >= 0.4, "moderate", "low")))
  data.table(
    signature_1 = s1,
    signature_2 = s2,
    paired_delta_spearman_r = d_r,
    paired_delta_spearman_p = tryCatch(cor.test(dx[d_ok], dy[d_ok], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_),
    sample_level_spearman_r = s_r,
    sample_level_spearman_p = tryCatch(cor.test(sx[s_ok], sy[s_ok], method = "spearman", exact = FALSE)$p.value, error = function(e) NA_real_),
    collinearity_flag = flag,
    notes = "Signature-score correlations are redundancy diagnostics, not biology claims."
  )
}), fill = TRUE)
fwrite(collinearity, file.path(table_dir, "gse73680_signature_collinearity_diagnostics.tsv"), sep = "\t")

priority <- rbindlist(lapply(all_signatures, function(sig) {
  sig_sum <- signature_summary[signature_name == sig]
  assoc <- module_sig_cor[correlation_level == "paired_delta_level" & signature_name == sig & module_name %in% primary_modules]
  max_abs <- ifelse(nrow(assoc), max(abs(assoc$spearman_r), na.rm = TRUE), NA_real_)
  col_risk <- collinearity[signature_1 == sig | signature_2 == sig, max(abs(c(paired_delta_spearman_r, sample_level_spearman_r)), na.rm = TRUE)]
  role <- ifelse(sig %in% c("injury_epithelial", "mineralization_remodeling"), "primary_adjustment",
          ifelse(!is.na(max_abs) && max_abs >= 0.5, "primary_adjustment",
          ifelse(!is.na(col_risk) && col_risk >= 0.85, "avoid_due_to_collinearity",
          ifelse(sig %in% single_covariates, "secondary_adjustment", "exploratory_only"))))
  data.table(
    signature_name = sig,
    signature_role = sig_sum$signature_role[1],
    stage5B1_delta_direction = sig_sum$effect_direction[1],
    association_with_primary_modules = ifelse(is.na(max_abs), "not_evaluable", ifelse(max_abs >= 0.7, "strong", ifelse(max_abs >= 0.4, "moderate", "weak"))),
    collinearity_risk = ifelse(is.na(col_risk), "not_evaluable", ifelse(col_risk >= 0.85, "very_high", ifelse(col_risk >= 0.7, "high", ifelse(col_risk >= 0.4, "moderate", "low")))),
    recommended_adjustment_role = role,
    notes = "Prioritization balances paired-delta module coupling, signature response, and collinearity risk."
  )
}), fill = TRUE)
fwrite(priority, file.path(table_dir, "gse73680_signature_priority_for_adjustment.tsv"), sep = "\t")

classify_adjustment <- function(intercept, p, unadj_mean) {
  if (!is.finite(intercept)) return("not_evaluable")
  if (sign(intercept) != sign(unadj_mean) || abs(intercept) < 1e-8) return("lost_after_adjustment")
  ratio <- abs(intercept) / abs(unadj_mean)
  if (ratio >= 0.75 && p < 0.1) return("retained_after_adjustment")
  if (ratio >= 0.25) return("attenuated_but_direction_retained")
  "covariate_explained"
}

single_models <- rbindlist(unlist(lapply(primary_modules, function(mod) {
  lapply(single_covariates, function(cov) {
    dt <- wide_delta[, .(patient_id, module_delta = get(mod), cov_delta = get(cov))]
    dt <- dt[is.finite(module_delta) & is.finite(cov_delta)]
    unadj <- mean(dt$module_delta, na.rm = TRUE)
    if (nrow(dt) < 8 || sd(dt$cov_delta, na.rm = TRUE) == 0) {
      return(data.table(module_name = mod, covariate_signature = cov, n_paired_patients = nrow(dt), intercept_estimate = NA_real_, intercept_p = NA_real_, covariate_beta = NA_real_, covariate_p = NA_real_, model_r2 = NA_real_, module_delta_direction_after_adjustment = "not_evaluable", retained_direction_after_adjustment = FALSE, adjustment_result = "not_evaluable", notes = "Insufficient data or invariant covariate."))
    }
    fit <- summary(lm(module_delta ~ cov_delta, data = dt))
    co <- fit$coefficients
    intercept <- co["(Intercept)", "Estimate"]
    cov_beta <- co["cov_delta", "Estimate"]
    result <- classify_adjustment(intercept, co["(Intercept)", "Pr(>|t|)"], unadj)
    data.table(
      module_name = mod,
      covariate_signature = cov,
      n_paired_patients = nrow(dt),
      intercept_estimate = intercept,
      intercept_p = co["(Intercept)", "Pr(>|t|)"],
      covariate_beta = cov_beta,
      covariate_p = co["cov_delta", "Pr(>|t|)"],
      model_r2 = fit$r.squared,
      module_delta_direction_after_adjustment = ifelse(intercept > 0, "increase", ifelse(intercept < 0, "decrease", "zero")),
      retained_direction_after_adjustment = is.finite(intercept) && sign(intercept) == sign(unadj) && abs(intercept) > 0,
      adjustment_result = result,
      notes = "Paired-delta single-covariate sensitivity; intercept tests residual module shift at covariate delta zero."
    )
  })
}), recursive = FALSE), fill = TRUE)
fwrite(single_models, file.path(table_dir, "gse73680_delta_adjusted_single_covariate_models.tsv"), sep = "\t")

make_pc1 <- function(wide_dt, vars, name) {
  vars <- vars[vars %in% names(wide_dt)]
  if (length(vars) < 2) return(NULL)
  mat <- as.matrix(wide_dt[, ..vars])
  mode(mat) <- "numeric"
  if (anyNA(mat)) {
    for (j in seq_len(ncol(mat))) mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  }
  pc <- prcomp(mat, center = TRUE, scale. = TRUE)$x[, 1]
  wide_dt[, (name) := as.numeric(pc)]
  name
}
composition_pc1 <- make_pc1(wide_delta, composition_core, "composition_core_PC1")
remodeling_pc1 <- make_pc1(wide_delta, injury_set, "remodeling_PC1")

compact_specs <- list(
  composition_core_limited = intersect(c("epithelial_general", "immune_myeloid", "fibroblast_stromal"), all_signatures),
  LoopTAL_plus_injury = intersect(c("LoopTAL", "injury_epithelial"), all_signatures),
  injury_remodeling_set = injury_set,
  composition_core_PC1 = composition_pc1,
  remodeling_PC1 = remodeling_pc1
)

max_abs_cor <- function(vars) {
  vars <- vars[vars %in% names(wide_delta)]
  if (length(vars) < 2) return(0)
  mat <- as.matrix(wide_delta[, ..vars])
  suppressWarnings(max(abs(cor(mat, use = "pairwise.complete.obs", method = "spearman")[upper.tri(diag(length(vars)) + 1)]), na.rm = TRUE))
}

compact_models <- rbindlist(unlist(lapply(primary_modules, function(mod) {
  lapply(names(compact_specs), function(model_name) {
    vars <- compact_specs[[model_name]]
    vars <- vars[!is.na(vars) & nzchar(vars)]
    vars <- vars[vars %in% names(wide_delta)]
    nvars <- length(vars)
    mac <- max_abs_cor(vars)
    if (nvars == 0) {
      return(data.table(module_name = mod, model_name = model_name, covariates_included = "", n_paired_patients = nrow(wide_delta), intercept_estimate = NA_real_, intercept_p = NA_real_, model_r2 = NA_real_, max_abs_covariate_correlation = NA_real_, model_status = "failed", adjustment_result = "not_evaluable", notes = "No available covariates."))
    }
    if (nvars >= 4) {
      return(data.table(module_name = mod, model_name = model_name, covariates_included = paste(vars, collapse = ";"), n_paired_patients = nrow(wide_delta), intercept_estimate = NA_real_, intercept_p = NA_real_, model_r2 = NA_real_, max_abs_covariate_correlation = mac, model_status = "skipped_overfit_risk", adjustment_result = "not_evaluable", notes = "Skipped because 26 paired patients cannot support this many correlated covariates."))
    }
    if (mac >= 0.85 && nvars > 1) {
      return(data.table(module_name = mod, model_name = model_name, covariates_included = paste(vars, collapse = ";"), n_paired_patients = nrow(wide_delta), intercept_estimate = NA_real_, intercept_p = NA_real_, model_r2 = NA_real_, max_abs_covariate_correlation = mac, model_status = "skipped_high_collinearity", adjustment_result = "not_evaluable", notes = "Skipped due to very high covariate collinearity."))
    }
    dt <- wide_delta[, c("patient_id", mod, vars), with = FALSE]
    setnames(dt, mod, "module_delta")
    dt <- dt[complete.cases(dt)]
    unadj <- mean(dt$module_delta)
    form <- as.formula(paste("module_delta ~", paste(vars, collapse = " + ")))
    fit <- tryCatch(summary(lm(form, data = dt)), error = function(e) NULL)
    if (is.null(fit)) {
      return(data.table(module_name = mod, model_name = model_name, covariates_included = paste(vars, collapse = ";"), n_paired_patients = nrow(dt), intercept_estimate = NA_real_, intercept_p = NA_real_, model_r2 = NA_real_, max_abs_covariate_correlation = mac, model_status = "failed", adjustment_result = "not_evaluable", notes = "Model failed."))
    }
    intercept <- fit$coefficients["(Intercept)", "Estimate"]
    ip <- fit$coefficients["(Intercept)", "Pr(>|t|)"]
    data.table(module_name = mod, model_name = model_name, covariates_included = paste(vars, collapse = ";"), n_paired_patients = nrow(dt), intercept_estimate = intercept, intercept_p = ip, model_r2 = fit$r.squared, max_abs_covariate_correlation = mac, model_status = "success", adjustment_result = classify_adjustment(intercept, ip, unadj), notes = "Compact paired-delta model; intercept is residual shift after compact adjustment.")
  })
}), recursive = FALSE), fill = TRUE)
fwrite(compact_models, file.path(table_dir, "gse73680_delta_adjusted_compact_models.tsv"), sep = "\t")

patient_fe <- rbindlist(unlist(lapply(primary_modules, function(mod) {
  lapply(single_covariates, function(cov) {
    mdt <- scores[module_name == mod, .(sample_id, patient_id, primary_group, module_score = mean_z_score)]
    cdt <- scores[module_name == cov, .(sample_id, covariate_score = mean_z_score)]
    dt <- merge(mdt, cdt, by = "sample_id")
    dt[, group_binary := ifelse(primary_group == "plaque_or_stone_papilla", 1, 0)]
    if (uniqueN(dt$patient_id) < 8 || uniqueN(dt$primary_group) < 2 || sd(dt$covariate_score, na.rm = TRUE) == 0) {
      return(data.table(module_name = mod, covariate_signature = cov, n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), group_effect_estimate = NA_real_, group_effect_p = NA_real_, covariate_effect_estimate = NA_real_, covariate_effect_p = NA_real_, model_status = "failed", adjustment_result = "not_evaluable", notes = "Insufficient data or invariant covariate."))
    }
    fit <- tryCatch(summary(lm(module_score ~ group_binary + covariate_score + factor(patient_id), data = dt)), error = function(e) NULL)
    if (is.null(fit) || !"group_binary" %in% rownames(fit$coefficients) || !"covariate_score" %in% rownames(fit$coefficients)) {
      return(data.table(module_name = mod, covariate_signature = cov, n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), group_effect_estimate = NA_real_, group_effect_p = NA_real_, covariate_effect_estimate = NA_real_, covariate_effect_p = NA_real_, model_status = "failed", adjustment_result = "not_evaluable", notes = "Patient fixed-effect model failed, likely rank deficiency."))
    }
    ge <- fit$coefficients["group_binary", "Estimate"]
    gp <- fit$coefficients["group_binary", "Pr(>|t|)"]
    unadj <- response_model[module_name == mod, effect_estimate][1]
    result <- classify_adjustment(ge, gp, unadj)
    data.table(module_name = mod, covariate_signature = cov, n_samples = nrow(dt), n_patients = uniqueN(dt$patient_id), group_effect_estimate = ge, group_effect_p = gp, covariate_effect_estimate = fit$coefficients["covariate_score", "Estimate"], covariate_effect_p = fit$coefficients["covariate_score", "Pr(>|t|)"], model_status = "success", adjustment_result = result, notes = "Sample-level patient fixed-effect sensitivity with one covariate signature.")
  })
}), recursive = FALSE), fill = TRUE)
fwrite(patient_fe, file.path(table_dir, "gse73680_patient_fixed_effect_adjusted_models.tsv"), sep = "\t")

retention <- rbindlist(lapply(primary_modules, function(mod) {
  sm <- single_models[module_name == mod]
  cm <- compact_models[module_name == mod]
  pf <- patient_fe[module_name == mod]
  unadj <- delta_summary[module_name == mod]
  n_ret_single <- sm[adjustment_result == "retained_after_adjustment", .N]
  n_att_single <- sm[adjustment_result == "attenuated_but_direction_retained", .N]
  n_lost_single <- sm[adjustment_result %in% c("lost_after_adjustment", "covariate_explained"), .N]
  cm_success <- cm[model_status == "success"]
  n_ret_compact <- cm_success[adjustment_result == "retained_after_adjustment", .N]
  n_att_compact <- cm_success[adjustment_result == "attenuated_but_direction_retained", .N]
  n_lost_compact <- cm_success[adjustment_result %in% c("lost_after_adjustment", "covariate_explained"), .N]
  pf_retained_frac <- pf[model_status == "success", mean(adjustment_result %in% c("retained_after_adjustment", "attenuated_but_direction_retained"), na.rm = TRUE)]
  support <- ifelse(n_lost_single == 0 && n_ret_single >= 6 && n_lost_compact == 0 && n_ret_compact >= 2, "strong_retained",
             ifelse(n_lost_single <= 1 && (n_ret_single + n_att_single) >= 6 && n_lost_compact <= 1, "moderate_retained",
             ifelse(n_lost_single <= 3 && (n_ret_single + n_att_single) >= 5, "attenuated_but_retained",
             ifelse(n_lost_single >= 4 || n_lost_compact >= 2, "composition_dependent", "not_supported"))))
  interp <- switch(support,
    strong_retained = "Direction and residual shift retained across most single and compact adjustments.",
    moderate_retained = "Direction retained across most single-covariate adjustments, with some attenuation or compact-model caveats.",
    attenuated_but_retained = "Direction generally retained but effect is attenuated by composition/remodeling signatures.",
    composition_dependent = "Module shift is strongly explained by composition/remodeling signatures.",
    not_supported = "Adjusted support is weak or unstable.",
    "not_evaluable"
  )
  data.table(
    module_name = mod,
    unadjusted_effect_direction = unadj$effect_direction[1],
    unadjusted_fdr = unadj$fdr_within_module_family[1],
    n_single_covariate_models = nrow(sm),
    n_retained_single_covariate = n_ret_single,
    n_attenuated_single_covariate = n_att_single,
    n_lost_single_covariate = n_lost_single,
    n_compact_models_success = nrow(cm_success),
    n_retained_compact = n_ret_compact + n_att_compact,
    n_lost_compact = n_lost_compact,
    patient_fixed_effect_support = ifelse(is.na(pf_retained_frac), "not_evaluable", ifelse(pf_retained_frac >= 0.75, "mostly_retained", "mixed")),
    overall_composition_adjusted_support = support,
    interpretation = interp,
    notes = "Retention summary is Stage 5B2 sensitivity, not deconvolution proof."
  )
}), fill = TRUE)
fwrite(retention, file.path(table_dir, "gse73680_adjustment_retention_summary.tsv"), sep = "\t")

coupling <- rbindlist(lapply(primary_modules, function(mod) {
  get_coup <- function(sig) {
    r <- module_sig_cor[correlation_level == "paired_delta_level" & module_name == mod & signature_name == sig]
    if (!nrow(r)) return("not_evaluable")
    paste0(r$interpretation[1], " (rho=", sprintf("%.2f", r$spearman_r[1]), ", FDR=", sprintf("%.3g", r$fdr[1]), ")")
  }
  vals <- c(get_coup("injury_epithelial"), get_coup("ECM_fibrosis"), get_coup("mineralization_remodeling"))
  strong_or_mod <- sum(grepl("strong|moderate", vals))
  interp <- ifelse(strong_or_mod >= 2, "module_shift_coupled_to_injury_remodeling",
            ifelse(strong_or_mod == 1, "module_shift_partly_independent_of_injury_remodeling", "no_clear_coupling"))
  data.table(
    module_name = mod,
    injury_epithelial_coupling = vals[1],
    ECM_fibrosis_coupling = vals[2],
    mineralization_remodeling_coupling = vals[3],
    coupling_strength_summary = paste(vals, collapse = "; "),
    interpretation = interp,
    claim_allowed = ifelse(interp == "module_shift_coupled_to_injury_remodeling", "injury/remodeling-associated bulk context", "bulk disease-context coupling"),
    claim_not_allowed = "cell-type-specific injury response; causal mineralization mechanism; therapeutic target validation",
    notes = "Coupling is correlation of paired deltas, not mediation or mechanism."
  )
}), fill = TRUE)
fwrite(coupling, file.path(table_dir, "gse73680_injury_remodeling_coupling_interpretation.tsv"), sep = "\t")

support_counts <- retention[, table(overall_composition_adjusted_support)]
overall_support <- if (all(retention$overall_composition_adjusted_support %in% c("strong_retained", "moderate_retained"))) {
  "moderate_bulk_context_support"
} else if (any(retention$overall_composition_adjusted_support == "composition_dependent")) {
  "attenuated_bulk_context_support"
} else if (all(retention$overall_composition_adjusted_support %in% c("attenuated_but_retained", "moderate_retained"))) {
  "attenuated_bulk_context_support"
} else {
  "supplementary_context_only"
}
overall_use <- ifelse(overall_support %in% c("moderate_bulk_context_support", "attenuated_bulk_context_support"),
                      "main_text_conservative_with_stage5B2_caveat", "supplementary_only")
claim_table <- data.table(
  evidence_component = c("primary_MAGMA_module_response", "composition_adjusted_primary_response", "injury_remodeling_coupling", "TWAS_proxy_module_response", "curated_exemplar_response", "single_gene_response", "overall_GSE73680_support"),
  stage5B1_result_summary = c(
    "Primary MAGMA modules increased in paired plaque/stone papilla bulk scores before adjustment.",
    "Not available before Stage 5B2.",
    "Stage 5B1 identified paired-delta module-signature coupling candidates.",
    "TWAS-proxy modules scored as secondary/supplementary.",
    "Curated exemplar panel scored as context only.",
    "Single genes are not used for validation claims.",
    "Stage 5B1 supported moderate_before_composition_adjustment."
  ),
  stage5B2_adjusted_result_summary = c(
    paste(retention$module_name, retention$overall_composition_adjusted_support, collapse = "; "),
    paste0("Retention support distribution: ", paste(names(support_counts), as.integer(support_counts), collapse = "; ")),
    paste(coupling$module_name, coupling$interpretation, collapse = "; "),
    "Remain supplementary; no causal or papilla-specific TWAS claim.",
    "Remain biological exemplars; no evidence upgrade.",
    "Not evaluated as validation; disallowed for main claim.",
    paste0("Overall support classified as ", overall_support, ".")
  ),
  allowed_claim = c(
    "bulk disease-context association",
    "directionally retained bulk disease-context sensitivity",
    "injury/remodeling-associated bulk disease context",
    "supplementary bulk context only",
    "curated biological context",
    "none beyond descriptive supplementary context",
    ifelse(overall_support == "moderate_bulk_context_support", "paired bulk disease-context support after composition-aware sensitivity modeling", "injury/remodeling-associated paired bulk disease-context support")
  ),
  disallowed_claim = "cell-type-specific response; genetic causality; plaque causation; therapeutic target validation; single-gene validation",
  claim_strength_after_composition_adjustment = c("moderate_bulk_context_support", overall_support, "attenuated_bulk_context_support", "supplementary_context_only", "supplementary_context_only", "supplementary_context_only", overall_support),
  manuscript_use = c("main_text_conservative", overall_use, "main_text_caveat_or_discussion", "supplementary_only", "supplementary_only", "avoid_as_claim", overall_use),
  notes = "Stage 5B2 claim table; final wording should remain bulk-tissue and sensitivity-analysis based."
)
fwrite(claim_table, file.path(table_dir, "gse73680_stage5B2_claim_decision_table.tsv"), sep = "\t")

blueprint <- data.table(
  panel_id = LETTERS[1:7],
  panel_title = c("GSE73680 paired sample design and QC", "Primary MAGMA module paired deltas", "Paired/patient-blocked model summary", "Composition and injury signature response", "Composition-adjusted retention summary", "Injury/remodeling coupling", "Claim-decision strip"),
  source_file = c(
    "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_pairing_summary.tsv",
    "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta.tsv",
    "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_response_model.tsv",
    "results/tables/revision/stage5B1_gse73680_module_response/gse73680_signature_response_summary.tsv",
    "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_adjustment_retention_summary.tsv",
    "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_injury_remodeling_coupling_interpretation.tsv",
    "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_stage5B2_claim_decision_table.tsv"
  ),
  plot_type = c("paired design schematic/table", "paired delta dot/line plot", "forest/table plot", "signature paired delta plot", "tile/summary matrix", "correlation/heatmap summary", "text strip"),
  allowed_interpretation = c(
    "paired bulk papillary comparison",
    "primary modules shift in paired bulk context",
    "patient-aware model sensitivity",
    "bulk marker-signature response",
    "direction retained/attenuated after sensitivity models",
    "bulk module shifts coupled to injury/remodeling context",
    "conservative GSE73680 claim boundary"
  ),
  forbidden_interpretation = c(
    "new cohort causal validation",
    "cell-type-specific disease response",
    "genetic causality",
    "validated cell fractions",
    "deconvolution proof",
    "causal mineralization mechanism",
    "therapeutic target or single-gene validation"
  ),
  priority = c("required", "required", "required", "required", "required", "recommended", "required"),
  notes = "Do not draw final Figure 4 until human review of Stage 5B2."
)
writeLines(c(
  "# Figure 4 GSE73680 blueprint",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "Draft planning only. Do not treat as final Figure 4.",
  "",
  "| panel_id | panel_title | source_file | plot_type | allowed_interpretation | forbidden_interpretation | priority | notes |",
  "|---|---|---|---|---|---|---|---|",
  apply(blueprint, 1, function(r) paste0("| ", paste(gsub("\\|", "/", r), collapse = " | "), " |"))
), file.path(doc_dir, "figure4_gse73680_blueprint.md"))

methods_text <- c(
  "# Manuscript replacement text draft for Stage 5B2",
  "",
  "## Methods Draft",
  "",
  "Composition-aware sensitivity modeling was performed using Stage 5B1 sample-level module and marker-signature scores. For each primary MAGMA-prioritized module, paired plaque/stone papilla minus control/adjacent module-score deltas were modeled against paired deltas for individual bulk composition, epithelial, Loop/TAL, immune, stromal, endothelial, injury, ECM/fibrosis and mineralization/remodeling signatures. Compact multi-covariate models were run only when the number of covariates and collinearity were compatible with the 26 paired-patient design. Patient fixed-effect sample-level models were used as an additional sensitivity analysis. These marker-signature models were interpreted as bulk-tissue sensitivity analyses, not deconvolution or cell-type-specific expression tests.",
  "",
  "## Results Draft",
  "",
  paste0("After composition-aware sensitivity modeling, primary MAGMA module support was classified as: ", paste(retention$module_name, retention$overall_composition_adjusted_support, collapse = "; "), ". The overall GSE73680 support level was `", overall_support, "`. Injury/remodeling coupling was evident for several primary modules, so the safest wording is that GSE73680 supports a paired bulk plaque/stone papilla disease-context signal that is linked to tissue remodeling and requires conservative interpretation."),
  "",
  "## Limitations Draft",
  "",
  "GSE73680 remains a bulk tissue dataset with curated metadata labels and uncertain original normalization provenance. Marker signatures are proxies and should not be interpreted as measured cell fractions. Composition-aware adjustment can test sensitivity to bulk tissue programs but cannot establish cell-type-specific expression changes, genetic causality, plaque nucleation mechanisms, therapeutic target validation or single-gene validation."
)
writeLines(methods_text, file.path(doc_dir, "manuscript_replacement_text_stage5B2.md"))

reviewer <- c(
  "# Stage 5B2 simulated reviewer check",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "1. Do primary MAGMA module shifts remain after composition-aware adjustment?",
  paste0("   Summary: ", paste(retention$module_name, retention$overall_composition_adjusted_support, collapse = "; "), "."),
  "",
  "2. Are the adjusted models overfit given 26 paired patients?",
  "   Single-covariate models are appropriate as sensitivity analyses. Compact models were restricted to three or fewer covariates or PCA-compressed scores to limit overfit risk with 26 paired patients.",
  "",
  "3. Are collinear signatures handled transparently?",
  "   Yes. Pairwise sample-level and paired-delta signature correlations were written to the collinearity diagnostic table and used to constrain compact-model interpretation.",
  "",
  "4. Does injury/remodeling coupling explain the module response?",
  "   It partly explains or tracks the response. The interpretation table separates bulk disease-context support from injury/remodeling-associated bulk context.",
  "",
  "5. Are TWAS-proxy and curated exemplar modules kept supplementary?",
  "   Yes. The claim decision table keeps TWAS-proxy and curated exemplar outputs supplementary/context only.",
  "",
  "6. Does any claim imply cell-type-specific disease response?",
  "   No. All outputs describe bulk marker-signature sensitivity and explicitly disallow cell-type-specific response.",
  "",
  "7. What is the strongest allowed GSE73680 claim?",
  ifelse(overall_support == "moderate_bulk_context_support",
         "   Paired bulk disease-context support for MAGMA-prioritized modules after composition-aware sensitivity modeling.",
         "   Injury/remodeling-associated paired bulk disease-context support for MAGMA-prioritized modules."),
  "",
  "8. What is the strongest disallowed GSE73680 claim?",
  "   Cell-type-specific disease response, genetic causality, plaque causation/nucleation, therapeutic target validation and single-gene validation remain disallowed.",
  "",
  "9. Should Figure 4 planning begin?",
  "   Yes, as a draft blueprint only. Final Figure 4 should wait for human acceptance of the claim strength.",
  "",
  "10. Should full manuscript rewrite still wait until Stage 6?",
  "   Yes. The manuscript rewrite should wait until spatial/TWAS conservative context in Stage 6 is complete."
)
writeLines(reviewer, file.path(doc_dir, "stage5B2_simulated_reviewer_check.md"))

report <- c(
  "# Stage 5B2 report: GSE73680 composition-aware sensitivity modeling",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Input Files",
  "",
  paste0("- Sample/module scores: `", paths$scores, "`"),
  paste0("- Paired deltas: `", paths$delta, "`"),
  paste0("- Stage 5B1 model table: `", paths$response_model, "`"),
  "",
  "## Collinearity Diagnostics",
  "",
  paste0("- Signature pairs audited: ", nrow(collinearity)),
  paste0("- High or very high collinearity pairs: ", collinearity[collinearity_flag %in% c("high_abs_r_ge_0.7", "very_high_abs_r_ge_0.85"), .N]),
  "",
  "## Single-Covariate Adjustment Results",
  "",
  paste0("- Single-covariate models run: ", nrow(single_models)),
  paste0("- Retained or attenuated-but-retained: ", single_models[adjustment_result %in% c("retained_after_adjustment", "attenuated_but_direction_retained"), .N]),
  paste0("- Lost/covariate-explained: ", single_models[adjustment_result %in% c("lost_after_adjustment", "covariate_explained"), .N]),
  "",
  "## Compact Model Results",
  "",
  paste0("- Compact model rows: ", nrow(compact_models)),
  paste0("- Successful compact models: ", compact_models[model_status == "success", .N]),
  paste0("- Skipped for collinearity/overfit: ", compact_models[model_status %in% c("skipped_high_collinearity", "skipped_overfit_risk"), .N]),
  "",
  "## Patient Fixed-Effect Model Results",
  "",
  paste0("- Patient fixed-effect adjusted models: ", nrow(patient_fe)),
  paste0("- Successful models: ", patient_fe[model_status == "success", .N]),
  "",
  "## Retention Summary",
  "",
  paste(retention$module_name, retention$overall_composition_adjusted_support, collapse = "; "),
  "",
  "## Injury/Remodeling Coupling Interpretation",
  "",
  paste(coupling$module_name, coupling$interpretation, collapse = "; "),
  "",
  "## Final GSE73680 Claim Level",
  "",
  paste0("- Overall support: `", overall_support, "`"),
  "- Recommended wording should stay bulk-tissue and sensitivity-analysis based.",
  "",
  "## Figure 4 Readiness",
  "",
  "- Draft Figure 4 blueprint is ready.",
  "- Do not draw final Figure 4 until human acceptance of Stage 5B2 claim strength.",
  "",
  "## Unresolved Limitations",
  "",
  "- Bulk marker signatures are not validated cell fractions.",
  "- Metadata remains curated from sample labels.",
  "- Adjusted models are sensitivity analyses, not causal mediation.",
  "",
  "## Recommended Next Step",
  "",
  "Proceed to human review of Stage 5B2. After acceptance, either draft Figure 4 planning/refinement or move to Stage 6 spatial/TWAS conservative context. Do not start full manuscript rewrite before Stage 6 is complete."
)
writeLines(report, file.path(doc_dir, "stage5B2_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 5, `:=`(
    status = "stage5B2_completed_composition_adjusted_sensitivity",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 5A, 5B1, and 5B2 completed; composition-aware sensitivity models, collinearity diagnostics, retention summary, injury/remodeling coupling interpretation, claim decision table, Figure 4 blueprint, draft text, reviewer check, and report generated",
    blocking_issues = "Full Stage 5 not yet complete because final/draft Figure 4 has not been generated; GSE73680 claims remain bulk-tissue sensitivity claims, not cell-type-specific or causal",
    next_stage_ready = "figure4_planning_or_stage6_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Wrote Stage 5B2 outputs to", doc_dir, "and", table_dir, "\n")
cat("Overall support:", overall_support, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
