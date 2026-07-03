suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

module_scores <- fread(file.path(table_dir, "gse73680_module_score_matrix.tsv"))
patient_scores <- fread(file.path(table_dir, "gse73680_patient_level_module_score_matrix.tsv"))
meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
expr <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
expr <- expr[!is.na(gene) & gene != ""]

marker_sets <- list(
  injury_remodeling = c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "FN1", "VIM", "HAVCR1", "LCN2"),
  inflammation_immune = c("CCL2", "CCL7", "CXCL8", "IL1B", "TNF", "CD68", "LST1"),
  fibrosis_ecm = c("COL1A1", "COL3A1", "DCN", "LUM", "ACTA2", "TAGLN"),
  epithelial_injury = c("HAVCR1", "LCN2", "VCAM1", "SOX9", "PROM1")
)

sample_ids <- intersect(meta$sample_id, names(expr))
meta <- meta[match(sample_ids, sample_id)]
emat <- as.matrix(expr[, ..sample_ids])
rownames(emat) <- expr$gene
mode(emat) <- "numeric"
emat <- log2(emat + 1)
zmat <- t(scale(t(emat)))
zmat[!is.finite(zmat)] <- NA_real_

inj_scores <- rbindlist(lapply(names(marker_sets), function(nm) {
  genes <- unique(marker_sets[[nm]])
  detected <- intersect(genes, rownames(zmat))
  scores <- if (length(detected) > 0) colMeans(zmat[detected, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(sample_ids))
  data.table(
    injury_module = nm,
    sample_id = sample_ids,
    injury_score = as.numeric(scores),
    n_genes_input = length(genes),
    n_genes_detected = length(detected),
    detected_fraction = length(detected) / length(genes)
  )
}), fill = TRUE)
fwrite(inj_scores, file.path(table_dir, "gse73680_injury_program_score_matrix.tsv"), sep = "\t")

risk_modules <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
sample_risk <- module_scores[module_name %in% risk_modules, .(module_name, sample_id, module_score)]
sample_dt <- merge(sample_risk, inj_scores[, .(injury_module, sample_id, injury_score)], by = "sample_id", allow.cartesian = TRUE)
sample_dt <- merge(sample_dt, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")

cor_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

sample_cor <- sample_dt[, {
  z <- cor_one(module_score, injury_score)
  .(n_samples = z$n, cor_method = "Spearman_sample_level", rho = z$rho, p_value = z$p)
}, by = .(module_name, injury_module)]
sample_cor[, fdr := p.adjust(p_value, method = "BH")]
sample_cor[, interpretation := fcase(
  rho >= 0.50 & fdr < 0.05, "strong disease-context coupling",
  rho >= 0.30 & fdr < 0.10, "moderate disease-context coupling",
  default = "weak_or_no_coupling"
)]
fwrite(sample_cor, file.path(table_dir, "gse73680_risk_injury_module_correlations.tsv"), sep = "\t")
fwrite(sample_cor, "results/tables/gse73680_risk_injury_module_correlations.tsv", sep = "\t")

inj_patient <- merge(inj_scores, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
inj_patient <- inj_patient[, .(patient_level_injury_score = mean(injury_score, na.rm = TRUE)),
                           by = .(injury_module, patient_id, group_curated)]
risk_patient <- patient_scores[module_name %in% risk_modules, .(module_name, patient_id, group_curated, patient_level_module_score)]

risk_wide <- dcast(risk_patient, module_name + patient_id ~ group_curated, value.var = "patient_level_module_score")
risk_wide <- risk_wide[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
risk_wide[, risk_delta := plaque_or_stone_papilla - control_or_adjacent]
inj_wide <- dcast(inj_patient, injury_module + patient_id ~ group_curated, value.var = "patient_level_injury_score")
inj_wide <- inj_wide[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
inj_wide[, injury_delta := plaque_or_stone_papilla - control_or_adjacent]

delta_dt <- merge(risk_wide[, .(module_name, patient_id, risk_delta)],
                  inj_wide[, .(injury_module, patient_id, injury_delta)],
                  by = "patient_id", allow.cartesian = TRUE)
delta_cor <- delta_dt[, {
  z <- cor_one(risk_delta, injury_delta)
  .(n_paired_patients = z$n, cor_method = "Spearman_paired_delta", rho = z$rho, p_value = z$p)
}, by = .(module_name, injury_module)]
delta_cor[, fdr := p.adjust(p_value, method = "BH")]
delta_cor[, interpretation := fcase(
  rho >= 0.50 & fdr < 0.05, "strong paired-delta coupling",
  rho >= 0.30 & fdr < 0.10, "moderate paired-delta coupling",
  default = "weak_or_no_paired_delta_coupling"
)]
fwrite(delta_cor, file.path(table_dir, "gse73680_risk_injury_paired_delta_correlations.tsv"), sep = "\t")
fwrite(delta_cor, "results/tables/gse73680_risk_injury_paired_delta_correlations.tsv", sep = "\t")

resid_cor <- sample_dt[, {
  d <- .SD[is.finite(module_score) & is.finite(injury_score)]
  if (nrow(d) < 10 || uniqueN(d$patient_id) < 3) {
    .(n_samples = nrow(d), cor_method = "patient_group_residual_spearman", rho = NA_real_, p_value = NA_real_)
  } else {
    d[, patient_factor := factor(patient_id)]
    fit_x <- lm(module_score ~ group_curated + patient_factor, data = d)
    fit_y <- lm(injury_score ~ group_curated + patient_factor, data = d)
    z <- cor_one(residuals(fit_x), residuals(fit_y))
    .(n_samples = z$n, cor_method = "patient_group_residual_spearman", rho = z$rho, p_value = z$p)
  }
}, by = .(module_name, injury_module)]
resid_cor[, fdr := p.adjust(p_value, method = "BH")]
resid_cor[, interpretation := fcase(
  rho >= 0.50 & fdr < 0.05, "strong residual coupling",
  rho >= 0.30 & fdr < 0.10, "moderate residual coupling",
  default = "weak_or_no_residual_coupling"
)]
fwrite(resid_cor, file.path(table_dir, "gse73680_risk_injury_partial_correlations.tsv"), sep = "\t")
fwrite(resid_cor, "results/tables/gse73680_risk_injury_partial_correlations.tsv", sep = "\t")

heat <- copy(sample_cor)
heat[, module_label := factor(module_name, levels = rev(risk_modules),
                              labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
heat[, injury_label := factor(injury_module,
                              levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                              labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
heat[, sig_label := fifelse(fdr < 0.05, "*", "")]
p_heat <- ggplot(heat, aes(injury_label, module_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sig_label)), size = 2.8, color = "#303030") +
  scale_fill_gradient2(low = "#8AA0A8", mid = "white", high = "#9A5F52", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "GSE73680 risk-module coupling with papillary injury programs",
       x = NULL, y = NULL, fill = "Spearman\nrho") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "gse73680_risk_injury_correlation_heatmap.pdf"), p_heat,
       width = 7.4, height = 4.8, units = "in", device = "pdf", bg = "white")
ggsave(file.path(fig_dir, "gse73680_risk_injury_correlation_heatmap.png"), p_heat,
       width = 7.4, height = 4.8, units = "in", dpi = 260, bg = "white")
ggsave("results/figures/gse73680_risk_injury_correlation_heatmap.pdf", p_heat,
       width = 7.4, height = 4.8, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/gse73680_risk_injury_correlation_heatmap.png", p_heat,
       width = 7.4, height = 4.8, units = "in", dpi = 260, bg = "white")

scatter_dt <- sample_dt[module_name == "MAGMA_top50" & injury_module == "injury_remodeling"]
p_scatter <- ggplot(scatter_dt, aes(injury_score, module_score, fill = group_curated)) +
  geom_point(shape = 21, size = 2.5, color = "#555555", stroke = 0.25, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "#3E6672", linewidth = 0.55) +
  scale_fill_manual(values = c(control_or_adjacent = "#8AA0A8", plaque_or_stone_papilla = "#B08A45"),
                    labels = c("Control/adjacent", "Plaque/stone papilla")) +
  labs(title = "MAGMA top 50 module versus injury/remodeling program",
       x = "Injury/remodeling program score", y = "MAGMA top 50 module score", fill = NULL) +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "gse73680_magma_vs_injury_scatter.pdf"), p_scatter,
       width = 5.8, height = 4.6, units = "in", device = "pdf", bg = "white")
ggsave(file.path(fig_dir, "gse73680_magma_vs_injury_scatter.png"), p_scatter,
       width = 5.8, height = 4.6, units = "in", dpi = 260, bg = "white")
ggsave("results/figures/gse73680_magma_vs_injury_scatter.pdf", p_scatter,
       width = 5.8, height = 4.6, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/gse73680_magma_vs_injury_scatter.png", p_scatter,
       width = 5.8, height = 4.6, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# GSE73680 Risk-Injury Coupling v0.1",
  "",
  "This analysis evaluates whether MAGMA-prioritized risk modules covary with papillary injury/remodeling programs in GSE73680.",
  "",
  "Reported correlations are module-level disease-context coupling analyses. They should not be interpreted as evidence that genetic risk modules causally drive injury or remodeling."
), "docs/gse73680_risk_injury_coupling_v0.1.md", useBytes = TRUE)

message("wrote GSE73680 risk-injury coupling outputs")
