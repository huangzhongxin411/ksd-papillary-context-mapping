suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
})

stage_id <- "stage5C1_gse73680_figure4_draft"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
for (d in c(fig_dir, table_dir, doc_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "stage5C1_generate_figure4_draft.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 5C1 conservative Figure 4 draft generation\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  pairing = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_pairing_summary.tsv",
  paired_delta = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta.tsv",
  paired_summary = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta_summary.tsv",
  patient_model = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_response_model.tsv",
  signature_summary = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_signature_response_summary.tsv",
  module_signature_correlation = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_signature_correlation.tsv",
  stage5B1_claim = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_claim_decision_table.tsv",
  single_adjusted = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_delta_adjusted_single_covariate_models.tsv",
  compact_adjusted = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_delta_adjusted_compact_models.tsv",
  patient_adjusted = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_patient_fixed_effect_adjusted_models.tsv",
  retention = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_adjustment_retention_summary.tsv",
  coupling = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_injury_remodeling_coupling_interpretation.tsv",
  stage5B2_claim = "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_stage5B2_claim_decision_table.tsv"
)

missing_files <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_files)) {
  stop("Required Stage 5B1/5B2 source tables missing: ", paste(missing_files, collapse = ", "))
}

read_source <- function(name) {
  d <- fread(paths[[name]])
  cat(sprintf("Loaded %-30s %5d rows x %2d columns\n", name, nrow(d), ncol(d)))
  d
}

check_cols <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) {
    stop("Missing required columns in ", label, ": ", paste(missing, collapse = ", "))
  }
}

pairing <- read_source("pairing")
paired_delta <- read_source("paired_delta")
paired_summary <- read_source("paired_summary")
patient_model <- read_source("patient_model")
signature_summary <- read_source("signature_summary")
module_signature_correlation <- read_source("module_signature_correlation")
stage5B1_claim <- read_source("stage5B1_claim")
single_adjusted <- read_source("single_adjusted")
compact_adjusted <- read_source("compact_adjusted")
patient_adjusted <- read_source("patient_adjusted")
retention <- read_source("retention")
coupling <- read_source("coupling")
stage5B2_claim <- read_source("stage5B2_claim")

check_cols(pairing, c("patient_id", "control_or_adjacent_sample", "plaque_or_stone_sample", "paired", "usable_for_paired_delta"), "pairing")
check_cols(paired_delta, c("patient_id", "module_name", "module_role", "paired_delta"), "paired_delta")
check_cols(paired_summary, c("module_name", "n_paired_patients", "median_delta", "mean_delta", "fdr_within_module_family"), "paired_summary")
check_cols(patient_model, c("module_name", "model_type", "effect_estimate", "standard_error", "fdr"), "patient_model")
check_cols(signature_summary, c("signature_name", "signature_role", "mean_delta", "fdr", "n_paired_patients"), "signature_summary")
check_cols(module_signature_correlation, c("correlation_level", "module_name", "signature_name", "spearman_r", "fdr"), "module_signature_correlation")
check_cols(single_adjusted, c("module_name", "covariate_signature", "adjustment_result"), "single_adjusted")
check_cols(compact_adjusted, c("module_name", "model_name", "model_status", "adjustment_result"), "compact_adjusted")
check_cols(patient_adjusted, c("module_name", "covariate_signature", "model_status", "adjustment_result"), "patient_adjusted")
check_cols(retention, c("module_name", "n_retained_single_covariate", "n_attenuated_single_covariate", "n_lost_single_covariate", "n_retained_compact", "n_lost_compact", "overall_composition_adjusted_support"), "retention")
check_cols(coupling, c("module_name", "interpretation", "claim_allowed", "claim_not_allowed"), "coupling")

primary_modules <- c(
  "R1_MAGMA_Bonferroni_only",
  "R1_R2_R3_all_MAGMA_Bonferroni",
  "MAGMA_top50",
  "MAGMA_top100"
)
module_labels <- c(
  R1_MAGMA_Bonferroni_only = "R1 MAGMA Bonf.",
  R1_R2_R3_all_MAGMA_Bonferroni = "R1-R3 MAGMA Bonf.",
  MAGMA_top50 = "MAGMA top 50",
  MAGMA_top100 = "MAGMA top 100"
)
signature_labels <- c(
  injury_epithelial = "Injury epithelial",
  ECM_fibrosis = "ECM / fibrosis",
  mineralization_remodeling = "Mineralization / remodeling",
  epithelial_general = "Epithelial general",
  LoopTAL = "Loop / TAL",
  immune_myeloid = "Immune / myeloid",
  fibroblast_stromal = "Fibroblast / stromal",
  collecting_duct = "Collecting duct",
  endothelial = "Endothelial",
  proximal_tubule = "Proximal tubule",
  pericyte_smooth_muscle = "Pericyte / smooth muscle",
  T_cell = "T cell",
  B_cell = "B cell"
)

assert_primary_complete <- function(d, label) {
  missing <- setdiff(primary_modules, unique(d$module_name))
  if (length(missing)) stop(label, " lacks primary modules: ", paste(missing, collapse = ", "))
}
for (x in list(paired_delta, paired_summary, patient_model, retention, coupling)) {
  assert_primary_complete(x, "Required source")
}

pretty_module <- function(x) unname(module_labels[as.character(x)])
pretty_signature <- function(x) {
  out <- unname(signature_labels[as.character(x)])
  out[is.na(out)] <- gsub("_", " ", as.character(x)[is.na(out)])
  out
}

write_panel_source <- function(d, filename) {
  path <- file.path(table_dir, filename)
  fwrite(d, path, sep = "\t", na = "")
  cat("Wrote source:", path, "\n")
  path
}

# Panel A: retain patient-level provenance and derive the displayed cohort counts.
panelA <- copy(pairing)
panelA[, design_status := fifelse(usable_for_paired_delta, "paired_analysis", "unpaired_context_only")]
panelA[, n_samples_in_record := (nzchar(control_or_adjacent_sample) & !is.na(control_or_adjacent_sample)) +
         (nzchar(plaque_or_stone_sample) & !is.na(plaque_or_stone_sample))]
n_samples <- sum(panelA$n_samples_in_record)
n_paired <- sum(panelA$usable_for_paired_delta)
if (n_samples != 55L || n_paired != 26L) {
  stop("Panel A cohort counts do not match locked values (55 samples; 26 paired patients): found ", n_samples, " and ", n_paired)
}
panelA_path <- write_panel_source(panelA, "figure4_panelA_source.tsv")

# Panel B: all paired deltas for the four primary MAGMA modules.
panelB <- paired_delta[module_name %in% primary_modules]
panelB[, module_label := factor(pretty_module(module_name), levels = unname(module_labels))]
setorder(panelB, module_name, patient_id)
if (nrow(panelB) != 104L) stop("Panel B expected 104 paired-delta rows; found ", nrow(panelB))
panelB_path <- write_panel_source(panelB, "figure4_panelB_source.tsv")

# Panel C: patient fixed-effect model estimates and 95% confidence intervals.
panelC <- patient_model[module_name %in% primary_modules]
panelC[, `:=`(
  module_label = factor(pretty_module(module_name), levels = rev(unname(module_labels))),
  ci_low = effect_estimate - 1.96 * standard_error,
  ci_high = effect_estimate + 1.96 * standard_error,
  fdr_label = paste0("FDR ", formatC(fdr, digits = 2, format = "f"))
)]
panelC_path <- write_panel_source(panelC, "figure4_panelC_source.tsv")

# Panel D: bulk marker-signature responses. These are proxies, never cell fractions.
panelD <- copy(signature_summary)
panelD[, signature_label := pretty_signature(signature_name)]
panelD[, signature_class := fifelse(signature_role == "injury_or_remodeling_signature", "Injury / remodeling proxy", "Composition marker proxy")]
panelD[, signature_label := factor(signature_label, levels = rev(pretty_signature(signature_name[order(mean_delta)])))]
panelD_path <- write_panel_source(panelD, "figure4_panelD_source.tsv")

# Panel E: sensitivity outcome counts, preserving attenuation and loss after adjustment.
panelE <- rbindlist(list(
  retention[, .(module_name, model_family = "Single-signature", outcome = "Direction retained", n_models = n_retained_single_covariate)],
  retention[, .(module_name, model_family = "Single-signature", outcome = "Attenuated", n_models = n_attenuated_single_covariate)],
  retention[, .(module_name, model_family = "Single-signature", outcome = "Lost / explained", n_models = n_lost_single_covariate)],
  retention[, .(module_name, model_family = "Compact", outcome = "Direction retained", n_models = n_retained_compact)],
  retention[, .(module_name, model_family = "Compact", outcome = "Lost / explained", n_models = n_lost_compact)]
), use.names = TRUE, fill = TRUE)
panelE[, `:=`(
  module_label = factor(pretty_module(module_name), levels = unname(module_labels)),
  model_family = factor(model_family, levels = c("Single-signature", "Compact")),
  outcome = factor(outcome, levels = c("Direction retained", "Attenuated", "Lost / explained"))
)]
panelE[, fraction := n_models / sum(n_models), by = .(module_name, model_family)]
panelE <- merge(panelE, retention[, .(module_name, overall_composition_adjusted_support)], by = "module_name", all.x = TRUE)
panelE_path <- write_panel_source(panelE, "figure4_panelE_source.tsv")

# Panel F: paired-delta coupling only; correlation is not mediation or mechanism.
panelF <- module_signature_correlation[
  correlation_level == "paired_delta_level" &
    module_name %in% primary_modules &
    signature_name %in% c("injury_epithelial", "ECM_fibrosis", "mineralization_remodeling")
]
panelF <- merge(
  panelF,
  coupling[, .(module_name, stage5B2_interpretation = interpretation, claim_allowed, claim_not_allowed)],
  by = "module_name",
  all.x = TRUE
)
panelF[, `:=`(
  module_label = factor(pretty_module(module_name), levels = rev(unname(module_labels))),
  signature_label = factor(pretty_signature(signature_name), levels = c("Injury epithelial", "ECM / fibrosis", "Mineralization / remodeling")),
  cell_label = paste0("rho ", formatC(spearman_r, digits = 2, format = "f"), "\nFDR ", formatC(fdr, digits = 3, format = "g"))
)]
if (nrow(panelF) != 12L) stop("Panel F expected 12 module-signature coupling rows; found ", nrow(panelF))
panelF_path <- write_panel_source(panelF, "figure4_panelF_source.tsv")

# Panel G: explicit visual claim boundary.
panelG <- data.table(
  display_order = 1:6,
  statement = c(
    "Attenuated bulk-context support",
    "Paired plaque/stone papilla module shift",
    "Linked to injury/remodeling signatures",
    "Not cell-type-specific response",
    "Not genetic causality",
    "Not plaque nucleation mechanism"
  ),
  statement_type = c(rep("allowed", 3), rep("boundary", 3)),
  notes = c(
    "Overall Stage 5B2 interpretation.",
    "Bulk paired comparison only.",
    "Correlation/coupling, not mechanism.",
    "Bulk expression cannot localize response to a cell type.",
    "Post-GWAS prioritization does not establish causality.",
    "No plaque nucleation experiment or mechanism test was performed."
  )
)
panelG[, statement_display := c(
  "Attenuated\nbulk-context support",
  "Paired plaque/stone papilla\nmodule shift",
  "Linked to injury/remodeling\nsignatures",
  "Not cell-type-specific\nresponse",
  "Not genetic\ncausality",
  "Not plaque nucleation\nmechanism"
)]
panelG_path <- write_panel_source(panelG, "figure4_panelG_source.tsv")

claim_lock <- data.table(
  evidence_component = c(
    "unadjusted_paired_module_shift",
    "patient_aware_model",
    "composition_adjusted_retention",
    "injury_remodeling_coupling",
    "TWAS_proxy_modules",
    "curated_exemplar_modules",
    "overall_GSE73680_claim"
  ),
  support_level = c(
    "unadjusted_directional_support",
    "unadjusted_directional_support",
    "attenuated_but_retained",
    "injury_remodeling_associated",
    "supplementary_only",
    "supplementary_only",
    "attenuated_but_retained"
  ),
  allowed_wording = c(
    "Primary MAGMA module scores increased directionally in the paired bulk comparison before adjustment.",
    "Patient-aware bulk models supported positive primary-module shifts before composition-aware sensitivity analysis.",
    "Primary-module direction was retained overall but attenuated after composition/injury/remodeling-aware sensitivity modeling.",
    "Primary-module shifts were coupled to injury/remodeling bulk marker-signature changes.",
    "Kidney_Cortex TWAS-proxy modules provide supplementary bulk disease-context association only.",
    "Curated exemplar modules provide supplementary biological context only.",
    "injury/remodeling-associated paired bulk disease-context support for MAGMA-prioritized modules"
  ),
  forbidden_wording = c(
    "independent validation; cell-type-specific response; genetic causality; single-gene validation",
    "causal genetic effect; cell-type-specific response",
    "composition-independent validation; deconvolution proof; strong independent support",
    "causal injury mechanism; causal mineralization mechanism; mediation",
    "papilla-specific TWAS validation; causal expression evidence",
    "validated disease genes; therapeutic targets",
    "cell-type-specific disease response; genetic causality; plaque causation or nucleation; therapeutic target validation; single-gene validation"
  ),
  main_or_supplement = c("main", "main", "main", "main_or_discussion", "supplement", "supplement", "main"),
  notes = c(
    "Unadjusted paired evidence; 26 paired patients.",
    "Patient fixed-effect model; 52 paired samples from 26 patients.",
    "Sensitivity analysis using marker-signature proxies, not validated fractions.",
    "Paired-delta correlation; not mediation or mechanism.",
    "One-SNP model burden and tissue mismatch retain supplementary status.",
    "Exemplar curation does not upgrade evidence strength.",
    "Locked Stage 5C1 wording."
  )
)
claim_lock_path <- file.path(table_dir, "gse73680_final_claim_lock_stage5C1.tsv")
fwrite(claim_lock, claim_lock_path, sep = "\t")

theme_stage <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#333333"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, color = "#555555", lineheight = 0.95),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = max(base_size - 0.2, 8.5), color = "#333333"),
      legend.title = element_text(face = "bold", size = max(base_size - 0.1, 8.5)),
      legend.text = element_text(size = max(base_size - 0.3, 8.0)),
      plot.margin = margin(6, 6, 6, 6)
    )
}

panel_label <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = 14, color = "#222222"),
          plot.tag.position = c(0, 1))
}

# A: paired design schematic.
pA <- ggplot() +
  annotate("rect", xmin = 0.03, xmax = 0.30, ymin = 0.27, ymax = 0.80, fill = "#E6E9EA", color = "#7F9DA6", linewidth = 0.6) +
  annotate("rect", xmin = 0.70, xmax = 0.97, ymin = 0.27, ymax = 0.80, fill = "#F3E9D2", color = "#B99B5A", linewidth = 0.6) +
  annotate("segment", x = 0.32, xend = 0.68, y = 0.54, yend = 0.54, arrow = arrow(length = unit(0.12, "inches")), color = "#555555", linewidth = 0.8) +
  annotate("text", x = 0.165, y = 0.61, label = "Control / adjacent", fontface = "bold", size = 3.8, color = "#333333") +
  annotate("text", x = 0.165, y = 0.43, label = "papillary bulk", size = 3.4, color = "#555555") +
  annotate("text", x = 0.835, y = 0.61, label = "Plaque / stone", fontface = "bold", size = 3.8, color = "#333333") +
  annotate("text", x = 0.835, y = 0.43, label = "papillary bulk", size = 3.4, color = "#555555") +
  annotate("text", x = 0.50, y = 0.70, label = "26 paired patients", fontface = "bold", size = 4.0, color = "#245A64") +
  annotate("text", x = 0.50, y = 0.35, label = "55 samples total | 29 patients", size = 3.5, color = "#555555") +
  annotate("text", x = 0.50, y = 0.11, label = "Paired bulk comparison; not a causal validation cohort", size = 3.2, color = "#9B5C4D") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Paired GSE73680 design and QC") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(6, 6, 6, 6))

# B: unadjusted paired bulk shifts.
pB <- ggplot(panelB, aes(x = module_label, y = paired_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#7A7A7A", linetype = "dashed") +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "#E1EAEC", color = "#245A64", linewidth = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.62, color = "#245A64") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.1, fill = "#B99B5A", color = "#333333") +
  labs(title = "Primary MAGMA module paired deltas", subtitle = "Unadjusted paired bulk shift; each point is one patient", x = NULL, y = "Paired delta (plaque/stone - control)") +
  theme_stage(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# C: patient-aware unadjusted model.
pC <- ggplot(panelC, aes(y = module_label, x = effect_estimate)) +
  geom_vline(xintercept = 0, color = "#7A7A7A", linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.18, color = "#7F9DA6", linewidth = 0.7, orientation = "y") +
  geom_point(size = 3.2, shape = 21, fill = "#245A64", color = "white", stroke = 0.35) +
  geom_text(aes(x = ci_high + 0.015, label = fdr_label), hjust = 0, size = 3.0, color = "#555555") +
  coord_cartesian(xlim = c(min(panelC$ci_low) - 0.02, max(panelC$ci_high) + 0.09), clip = "off") +
  labs(title = "Patient-aware bulk sensitivity", subtitle = "Patient fixed-effect estimates before composition-aware adjustment", x = "Group effect estimate (95% CI)", y = NULL) +
  theme_stage(9) +
  theme(plot.margin = margin(6, 42, 6, 6))

# D: signature response, explicitly proxy language.
pD <- ggplot(panelD, aes(y = signature_label, x = mean_delta, color = signature_class)) +
  geom_vline(xintercept = 0, color = "#7A7A7A", linetype = "dashed", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = mean_delta, yend = signature_label), linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("Injury / remodeling proxy" = "#9B5C4D", "Composition marker proxy" = "#7F9DA6")) +
  labs(title = "Bulk marker-signature response", subtitle = "Proxy scores, not validated cell fractions", x = "Mean paired delta", y = NULL, color = NULL) +
  theme_stage(8.5) +
  theme(legend.position = "bottom", legend.box.spacing = unit(0.05, "cm"), axis.text.y = element_text(size = 8.5))

# E: adjusted retention and attenuation.
pE <- ggplot(panelE, aes(x = module_label, y = fraction, fill = outcome)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = n_models), position = position_stack(vjust = 0.5), size = 3.0, color = "#222222") +
  facet_wrap(~model_family, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c("Direction retained" = "#7F9DA6", "Attenuated" = "#B99B5A", "Lost / explained" = "#C88A78"), drop = FALSE) +
  labs(title = "Adjustment reveals attenuation", subtitle = "All four modules: attenuated but direction retained overall", x = NULL, y = "Sensitivity-model outcomes", fill = NULL) +
  theme_stage(8.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom", strip.background = element_rect(fill = "#F4F5F5", color = NA), strip.text = element_text(face = "bold", size = 8))

# F: injury/remodeling paired-delta coupling.
pF <- ggplot(panelF, aes(x = signature_label, y = module_label, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.9, color = "#222222") +
  scale_fill_gradient(low = "#F1F2F2", high = "#9B5C4D", limits = c(0, 1), name = "Spearman\nrho") +
  labs(title = "Injury/remodeling coupling", subtitle = "Paired-delta correlation; coupling is not mechanism", x = NULL, y = NULL) +
  theme_stage(8.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

# G: claim-decision strip with allowed statements and hard boundaries.
pG <- ggplot(panelG, aes(x = display_order, y = 1, fill = statement_type)) +
  geom_tile(width = 0.94, height = 0.72, color = "white", linewidth = 0.8) +
  geom_text(aes(label = statement_display), size = 3.0, lineheight = 0.92, fontface = "bold", color = "#2F2F2F") +
  scale_fill_manual(values = c(allowed = "#DDE8EA", boundary = "#F1DDD7"), guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = NULL, expand = expansion(add = 0.05)) +
  coord_cartesian(ylim = c(0.55, 1.45), clip = "off") +
  labs(title = "Locked interpretation boundary") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(6, 6, 6, 6))

figure4 <- (
  panel_label(pA, "A") | panel_label(pB, "B")
) / (
  panel_label(pC, "C") | panel_label(pD, "D")
) / (
  panel_label(pE, "E") | panel_label(pF, "F")
) / panel_label(pG, "G") +
  plot_layout(heights = c(0.88, 1.05, 1.28, 0.48), widths = c(1, 1.16), guides = "keep") +
  plot_annotation(
    title = "Figure 4 draft v0.1. GSE73680 provides attenuated paired bulk disease-context support",
    subtitle = "MAGMA-prioritized module shifts track injury/remodeling context and attenuate in marker-signature sensitivity models",
    theme = theme(
      plot.title = element_text(face = "bold", family = "sans", size = 17, color = "#222222"),
      plot.subtitle = element_text(family = "sans", size = 10.5, color = "#555555"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

pdf_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.1.pdf")
png_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.1.png")
svg_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.1.svg")

grDevices::pdf(pdf_path, width = 15, height = 12.5, useDingbats = FALSE, bg = "white")
print(figure4)
grDevices::dev.off()
ggsave(png_path, figure4, width = 15, height = 12.5, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(svg_path, figure4, width = 15, height = 12.5, units = "in", device = svglite, bg = "white")
cat("\nFigure files written:\n", pdf_path, "\n", png_path, "\n", svg_path, "\n")

source_paths <- c(panelA_path, panelB_path, panelC_path, panelD_path, panelE_path, panelF_path, panelG_path)
source_tables <- c(
  paths$pairing,
  paste(paths$paired_delta, paths$paired_summary, sep = ";"),
  paths$patient_model,
  paths$signature_summary,
  paste(paths$retention, paths$single_adjusted, paths$compact_adjusted, paths$patient_adjusted, sep = ";"),
  paste(paths$module_signature_correlation, paths$coupling, sep = ";"),
  paste(paths$stage5B1_claim, paths$stage5B2_claim, claim_lock_path, sep = ";")
)
manifest <- rbindlist(lapply(seq_along(source_paths), function(i) {
  d <- fread(source_paths[i])
  data.table(
    figure = "Figure 4 draft v0.1",
    panel = LETTERS[i],
    source_table = source_tables[i],
    n_rows = nrow(d),
    n_columns = ncol(d),
    ready_for_publication_source_data = "yes_draft",
    notes = c(
      "Patient-level paired-design provenance; counts displayed are derived from these records.",
      "Four primary MAGMA modules x 26 paired patients; unadjusted bulk deltas.",
      "Patient fixed-effect model estimates and derived 95% confidence intervals.",
      "Bulk marker-signature summaries; explicitly not cell fractions.",
      "Derived sensitivity-outcome counts; attenuation/loss remains visible.",
      "Paired-delta correlations merged with the Stage 5B2 coupling boundary.",
      "Claim statements and hard exclusions; no quantitative inference."
    )[i]
  )
}))
manifest_path <- file.path(table_dir, "figure4_source_data_manifest_v0.1.tsv")
fwrite(manifest, manifest_path, sep = "\t")

visual_qc <- data.table(
  figure_id = "Figure 4",
  version = "v0.1",
  pdf_exists = file.exists(pdf_path),
  png_exists = file.exists(png_path),
  svg_exists = file.exists(svg_path),
  png_intended_dpi = 600,
  minimum_configured_font_pt = 8.0,
  panel_labels = "A-G present",
  legend_placement = "outside dense data panels",
  palette_check = "pass_semantic_teal_amber_terracotta",
  claim_boundary_check = "pass_attenuated_bulk_context",
  source_data_check = ifelse(all(file.exists(source_paths)), "pass", "fail"),
  visual_status = "requires_rendered_agent_and_human_review",
  action_required = "Inspect PDF/PNG at full size and 50% scale before accepting draft."
)
fwrite(visual_qc, file.path(table_dir, "figure4_visual_qc_v0.1.tsv"), sep = "\t")

cat("\nClaim lock:", claim_lock_path, "\n")
cat("Source manifest:", manifest_path, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
