suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

sample_cor <- fread("results/tables/gse73680_risk_injury_module_correlations.tsv")
delta_cor <- fread("results/tables/gse73680_risk_injury_paired_delta_correlations.tsv")
resid_cor <- fread("results/tables/gse73680_risk_injury_partial_correlations.tsv")

patient_scores <- fread("results/gse73680/tables/gse73680_patient_level_module_score_matrix.tsv")
inj_scores <- fread("results/gse73680/tables/gse73680_injury_program_score_matrix.tsv")
meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
inj_patient <- merge(inj_scores, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
inj_patient <- inj_patient[, .(patient_level_injury_score = mean(injury_score, na.rm = TRUE)),
                           by = .(injury_module, patient_id, group_curated)]

risk_modules <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
risk_w <- dcast(patient_scores[module_name %in% risk_modules],
                module_name + patient_id ~ group_curated,
                value.var = "patient_level_module_score")
risk_w <- risk_w[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
risk_w[, risk_delta := plaque_or_stone_papilla - control_or_adjacent]
inj_w <- dcast(inj_patient, injury_module + patient_id ~ group_curated,
               value.var = "patient_level_injury_score")
inj_w <- inj_w[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
inj_w[, injury_delta := plaque_or_stone_papilla - control_or_adjacent]
delta_dt <- merge(risk_w[, .(module_name, patient_id, risk_delta)],
                  inj_w[, .(injury_module, patient_id, injury_delta)],
                  by = "patient_id", allow.cartesian = TRUE)

cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

loo <- delta_dt[, {
  pats <- unique(patient_id)
  rbindlist(lapply(pats, function(pid) {
    d <- .SD[patient_id != pid]
    z <- cor_one(d$risk_delta, d$injury_delta)
    data.table(left_out_patient = pid, n_paired_patients = z[["n"]],
               rho = z[["rho"]], p_value = z[["p"]])
  }))
}, by = .(module_name, injury_module)]
loo[, direction_retained := rho > 0]
loo_summary <- loo[, .(
  n_leave_one_runs = .N,
  min_rho = min(rho, na.rm = TRUE),
  median_rho = median(rho, na.rm = TRUE),
  max_rho = max(rho, na.rm = TRUE),
  fraction_positive = mean(direction_retained, na.rm = TRUE)
), by = .(module_name, injury_module)]
loo_summary[, robustness_call := fcase(
  fraction_positive >= 0.95 & min_rho > 0, "direction_stable",
  fraction_positive >= 0.80, "mostly_direction_stable",
  default = "sensitive_to_patient_removal"
)]
fwrite(loo, "results/tables/gse73680_risk_injury_leave_one_patient_out.tsv", sep = "\t")
fwrite(loo_summary, "results/tables/gse73680_risk_injury_leave_one_patient_out_summary.tsv", sep = "\t")

sample_cor[, analysis := "Sample-level"]
delta_cor[, analysis := "Paired delta"]
resid_cor[, analysis := "Patient/group residual"]
common_cols <- c("module_name", "injury_module", "analysis", "rho", "p_value", "fdr", "interpretation")
robust <- rbindlist(list(sample_cor[, ..common_cols], delta_cor[, ..common_cols], resid_cor[, ..common_cols]), fill = TRUE)
robust <- merge(robust, loo_summary, by = c("module_name", "injury_module"), all.x = TRUE)
robust[, robustness_summary := fcase(
  analysis == "Paired delta" & rho > 0 & fdr < 0.10 & robustness_call %in% c("direction_stable", "mostly_direction_stable"), "robust paired module-level coupling",
  analysis == "Patient/group residual" & rho > 0 & fdr < 0.10, "residual coupling support",
  analysis == "Sample-level" & rho > 0 & fdr < 0.10, "sample-level coupling support",
  default = "limited support"
)]
fwrite(robust, "results/tables/gse73680_risk_injury_correlation_robustness.tsv", sep = "\t")

plot_dt <- robust[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates") &
                    injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
plot_dt[, module_label := factor(module_name,
                                 levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                                 labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
plot_dt[, injury_label := factor(injury_module,
                                 levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                                 labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
plot_dt[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual", "Sample-level"))]

p <- ggplot(plot_dt[analysis != "Sample-level"], aes(injury_label, module_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
  facet_wrap(~ analysis, ncol = 1) +
  scale_fill_gradient2(low = "#8AA0A8", mid = "white", high = "#9A5F52", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Risk-injury coupling robustness in GSE73680",
       subtitle = "Primary robustness views: paired patient delta and patient/group residual correlations",
       x = NULL, y = NULL, fill = "rho") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555555"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid = element_blank())
ggsave("results/figures/figure5_risk_injury_correlation_robustness.pdf", p,
       width = 8.8, height = 8.0, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure5_risk_injury_correlation_robustness.png", p,
       width = 8.8, height = 8.0, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# GSE73680 Risk-Injury Correlation Robustness v0.1",
  "",
  "This analysis separates sample-level correlations from paired patient delta correlations and patient/group residual correlations.",
  "Main text should emphasize paired-delta and residual support when present. Sample-level correlations are useful descriptively but are more vulnerable to disease-group composition effects."
), "docs/gse73680_risk_injury_correlation_robustness_v0.1.md", useBytes = TRUE)

message("wrote GSE73680 risk-injury robustness outputs")
