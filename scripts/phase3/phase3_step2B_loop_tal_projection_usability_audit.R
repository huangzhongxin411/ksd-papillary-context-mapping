#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

root <- getwd()
dir.create(file.path(root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "scripts/08_spatial_analysis"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
bool_chr <- function(x) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))

required <- c(
  "notes/phase3_step2_report.md",
  "results/tables/phase3_step2_label_transfer_summary.tsv",
  "results/tables/phase3_step2_reference_compartment_summary.tsv",
  "results/tables/phase3_step2_condition_context_summary.tsv",
  "notes/phase3_step2_spatial_label_transfer_claim_boundary.md",
  "notes/phase3_step2_limitations_and_next_steps.md",
  "config/spatial_sample_metadata_curated_phase3.tsv"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "))

score_files <- list.files(
  file.path(root, "results/spatial/phase3_step2_label_transfer"),
  pattern = "^label_transfer_scores.tsv.gz$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(score_files)) stop("No label_transfer_scores.tsv.gz files found.")

scores <- rbindlist(lapply(score_files, function(path) {
  dt <- fread(path)
  dt[, source_file := sub(paste0("^", root, "/?"), "", path)]
  dt
}), fill = TRUE)
scores[, in_tissue := as.integer(in_tissue)]
scores[, Loop_TAL_prediction_score := as.numeric(Loop_TAL_prediction_score)]
scores[, prediction_score_max := as.numeric(prediction_score_max)]

meta <- fread(file.path(root, "config/spatial_sample_metadata_curated_phase3.tsv"))
meta_small <- meta[, .(sample_id, section_id, control_or_disease, metadata_confidence)]
scores <- merge(scores, meta_small, by = c("sample_id", "section_id"), all.x = TRUE)
in_tissue <- scores[in_tissue == 1]

top_q_stats <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(list(threshold = NA_real_, n_above = NA_integer_, n_equal = NA_integer_, tie_fraction = NA_real_))
  }
  threshold <- as.numeric(quantile(x, 0.75, na.rm = TRUE, names = FALSE, type = 7))
  list(
    threshold = threshold,
    n_above = sum(x > threshold),
    n_equal = sum(x == threshold),
    tie_fraction = sum(x == threshold) / length(x)
  )
}

usability <- in_tissue[, {
  x <- Loop_TAL_prediction_score
  tq <- top_q_stats(x)
  pct_nonzero <- 100 * mean(x > 0, na.rm = TRUE)
  sparse <- pct_nonzero < 5 || uniqueN(x[!is.na(x)]) < 10
  continuous <- pct_nonzero >= 5 && uniqueN(x[!is.na(x)]) >= 20 && max(x, na.rm = TRUE) > 0
  rank_ok <- !is.na(tq$tie_fraction) && tq$tie_fraction <= 0.35 && pct_nonzero >= 5
  rec <- if (continuous) {
    "continuous_Loop_TAL_prediction_score"
  } else if (pct_nonzero >= 5) {
    "within_section_Loop_TAL_rank_percentile"
  } else if (pct_nonzero > 0) {
    "nonzero_Loop_TAL_indicator"
  } else {
    "not_usable_for_claim_grade_correlation"
  }
  note <- if (tq$threshold == 0 && tq$tie_fraction > 0.35) {
    "Top-quartile threshold is zero and highly tie-sensitive; use continuous scores, nonzero indicator, or rank only as descriptive support."
  } else if (continuous) {
    "Continuous Loop/TAL score has nonzero dynamic range but remains a label-transfer projection, not a cell fraction."
  } else {
    "Loop/TAL score is sparse; avoid claim-grade correlation."
  }
  .(
    n_spots_in_tissue = .N,
    n_nonzero_Loop_TAL_score_spots = sum(x > 0, na.rm = TRUE),
    percent_nonzero_Loop_TAL_score_spots = pct_nonzero,
    median_Loop_TAL_score = median(x, na.rm = TRUE),
    mean_Loop_TAL_score = mean(x, na.rm = TRUE),
    max_Loop_TAL_score = max(x, na.rm = TRUE),
    n_unique_Loop_TAL_score_values = uniqueN(x[!is.na(x)]),
    top_quartile_threshold = tq$threshold,
    n_spots_above_top_quartile_threshold = tq$n_above,
    n_spots_equal_to_top_quartile_threshold = tq$n_equal,
    top_quartile_tie_fraction = tq$tie_fraction,
    continuous_score_usable = continuous,
    rank_top_quartile_usable = rank_ok,
    recommended_Loop_TAL_metric_for_step3 = rec,
    notes = note
  )
}, by = .(sample_id, section_id, control_or_disease, metadata_confidence)]
setorder(usability, sample_id, section_id)
write_tsv(usability, file.path(root, "results/tables/phase3_step2B_loop_tal_projection_usability.tsv"))

composition <- in_tissue[, .(
  n_spots = .N,
  median_prediction_score_for_predicted_compartment = median(prediction_score_max, na.rm = TRUE)
), by = .(sample_id, section_id, control_or_disease, metadata_confidence, predicted_compartment)]
composition[, percent_spots := 100 * n_spots / sum(n_spots), by = .(sample_id, section_id)]
composition[, notes := ifelse(
  predicted_compartment == "Loop/TAL",
  "Loop/TAL is max-predicted for this fraction of in-tissue spots.",
  "Broad-compartment max-predicted label from Seurat TransferData."
)]
setcolorder(composition, c(
  "sample_id", "section_id", "control_or_disease", "metadata_confidence",
  "predicted_compartment", "n_spots", "percent_spots",
  "median_prediction_score_for_predicted_compartment", "notes"
))
setorder(composition, sample_id, section_id, -percent_spots)
write_tsv(composition, file.path(root, "results/tables/phase3_step2B_predicted_compartment_composition.tsv"))

step2_summary <- fread(file.path(root, "results/tables/phase3_step2_label_transfer_summary.tsv"))
quality <- merge(
  step2_summary[, .(sample_id, section_id, n_spots_in_tissue, median_prediction_score_max, percent_low_confidence_spots)],
  usability[, .(
    sample_id, section_id, metadata_confidence,
    percent_nonzero_Loop_TAL_score_spots, n_unique_Loop_TAL_score_values,
    recommended_Loop_TAL_metric_for_step3
  )],
  by = c("sample_id", "section_id"),
  all.x = TRUE
)
quality[, Loop_TAL_score_sparsity_flag := percent_nonzero_Loop_TAL_score_spots < 5 | n_unique_Loop_TAL_score_values < 10]
quality[, low_confidence_flag := percent_low_confidence_spots > 10]
quality[, include_for_step3 := fifelse(
  low_confidence_flag & Loop_TAL_score_sparsity_flag, "human_review",
  fifelse(low_confidence_flag | Loop_TAL_score_sparsity_flag, "yes_with_caution", "yes")
)]
quality[, recommended_use := fifelse(
  low_confidence_flag, "descriptive_overlay_only",
  fifelse(Loop_TAL_score_sparsity_flag, "descriptive_overlay_only", "within_section_context_projection")
)]
quality[, notes := fifelse(
  include_for_step3 == "yes",
  "Retain for Step 3 as contextual projection; avoid disease/control and ROI claims.",
  "Retain only with caution for descriptive context; avoid claim-grade Loop/TAL correlation."
)]
quality[, `:=`(
  Loop_TAL_score_sparsity_flag = bool_chr(Loop_TAL_score_sparsity_flag),
  low_confidence_flag = bool_chr(low_confidence_flag)
)]
setcolorder(quality, c(
  "sample_id", "section_id", "n_spots_in_tissue", "median_prediction_score_max",
  "percent_low_confidence_spots", "Loop_TAL_score_sparsity_flag",
  "low_confidence_flag", "metadata_confidence", "include_for_step3",
  "recommended_use", "notes"
))
setorder(quality, sample_id, section_id)
quality_out <- quality[, .(
  sample_id,
  section_id,
  n_spots_in_tissue,
  median_prediction_score_max,
  percent_low_confidence_spots,
  Loop_TAL_score_sparsity_flag,
  low_confidence_flag,
  metadata_confidence,
  include_for_step3,
  recommended_use,
  notes
)]
write_tsv(quality_out, file.path(root, "results/tables/phase3_step2B_section_quality_flags.tsv"))

in_tissue[, section_label := paste(sample_id, section_id, sep = "\n")]
in_tissue[, Loop_TAL_nonzero := Loop_TAL_prediction_score > 0]
in_tissue[, Loop_TAL_score_for_plot := pmax(Loop_TAL_prediction_score, 1e-8)]

p_dist <- ggplot(in_tissue, aes(x = Loop_TAL_score_for_plot)) +
  geom_histogram(bins = 50, fill = "#4C78A8", color = "white", linewidth = 0.12) +
  geom_vline(xintercept = 1e-8, linetype = "dashed", linewidth = 0.25, color = "grey35") +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~ section_id, scales = "free_y", ncol = 2) +
  labs(
    title = "Loop/TAL prediction-score distribution by section",
    subtitle = "Zero scores are shown at 1e-8 on the log-scale axis",
    x = "Loop/TAL prediction score (log10 scale; zeros offset)",
    y = "In-tissue spots"
  ) +
  theme_bw(base_size = 8)
ggsave(file.path(root, "results/figures/phase3_step2B_Loop_TAL_score_distribution_by_section.pdf"),
       p_dist, width = 8.5, height = 11, units = "in", device = "pdf")

p_nonzero <- ggplot(usability, aes(x = reorder(section_id, percent_nonzero_Loop_TAL_score_spots), y = percent_nonzero_Loop_TAL_score_spots, fill = control_or_disease)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = paste0(metadata_confidence, "\n", fmt(percent_nonzero_Loop_TAL_score_spots, 1), "%")), hjust = -0.05, size = 2.2) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Nonzero Loop/TAL prediction-score fraction by section",
    x = NULL,
    y = "Spots with Loop/TAL score > 0 (%)",
    fill = "Context"
  ) +
  theme_bw(base_size = 8) +
  theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
ggsave(file.path(root, "results/figures/phase3_step2B_Loop_TAL_nonzero_fraction_by_section.pdf"),
       p_nonzero, width = 8.5, height = 5.8, units = "in", device = "pdf")

p_comp <- ggplot(composition, aes(x = section_id, y = percent_spots, fill = predicted_compartment)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Max-predicted broad-compartment composition",
    x = NULL,
    y = "In-tissue spots (%)",
    fill = "Predicted compartment"
  ) +
  theme_bw(base_size = 8)
ggsave(file.path(root, "results/figures/phase3_step2B_predicted_compartment_composition.pdf"),
       p_comp, width = 8.5, height = 6.2, units = "in", device = "pdf")

p_overview <- ggplot(in_tissue, aes(x = x_coord, y = -y_coord)) +
  geom_point(aes(color = Loop_TAL_prediction_score), size = 0.45, alpha = 0.9) +
  scale_color_viridis_c(
    option = "magma",
    trans = "sqrt",
    limits = c(0, max(in_tissue$Loop_TAL_prediction_score, na.rm = TRUE)),
    oob = squish,
    name = "Loop/TAL score"
  ) +
  facet_wrap(~ section_id, scales = "free", ncol = 2) +
  labs(
    title = "Refined Loop/TAL label-transfer projection overview",
    subtitle = "label-transfer projection; no ROI annotation; sparse Loop/TAL score",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 8) +
  theme(
    strip.text = element_text(size = 6),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
ggsave(file.path(root, "results/figures/phase3_step2B_Loop_TAL_projection_overview_refined.pdf"),
       p_overview, width = 8.5, height = 11, units = "in", device = "pdf")

mean_nonzero <- mean(usability$percent_nonzero_Loop_TAL_score_spots, na.rm = TRUE)
range_nonzero <- range(usability$percent_nonzero_Loop_TAL_score_spots, na.rm = TRUE)
max_tie <- max(usability$top_quartile_tie_fraction, na.rm = TRUE)
continuous_sections <- usability[continuous_score_usable == TRUE, .N]
rank_sections <- usability[rank_top_quartile_usable == TRUE, .N]
low_conf_sections <- quality[low_confidence_flag == "TRUE", .N]
sparse_sections <- quality[Loop_TAL_score_sparsity_flag == "TRUE", .N]
include_counts <- quality[, .N, by = include_for_step3]
loop_max_spots <- composition[predicted_compartment == "Loop/TAL", sum(n_spots)]
if (!is.finite(loop_max_spots)) loop_max_spots <- 0L

primary_decision <- if (continuous_sections >= 8 && max_tie <= 0.35) {
  "A. Use continuous Loop/TAL prediction score for within-section co-distribution."
} else if (continuous_sections >= 8) {
  "A. Use continuous Loop/TAL prediction score for within-section co-distribution; keep top-quartile rank only as a descriptive overlay."
} else if (rank_sections >= 8) {
  "B. Use within-section Loop/TAL rank percentile only."
} else {
  "D. Do not use Loop/TAL score for claim-grade co-distribution; proceed with broad compartment projection only."
}
primary_metric <- if (grepl("^A", primary_decision)) "continuous_Loop_TAL_prediction_score" else if (grepl("^B", primary_decision)) "within_section_Loop_TAL_rank_percentile" else "not_usable_for_claim_grade_correlation"

metric_note <- c(
  "# Phase 3-Step 2B Step 3 metric decision",
  "",
  paste0("- Continuous Loop/TAL prediction score suitability: ", continuous_sections, " of ", nrow(usability), " sections met the pragmatic nonzero/dynamic-range criteria."),
  paste0("- Top-quartile rank usability: ", rank_sections, " of ", nrow(usability), " sections passed the tie-sensitivity screen; maximum threshold tie fraction was ", fmt(max_tie, 3), "."),
  paste0("- Loop/TAL was the max-predicted compartment for ", loop_max_spots, " in-tissue spots across all audited sections."),
  "- Top-quartile rank should be used only as a descriptive overlay because Loop/TAL medians are zero and threshold ties can distort high/low calls.",
  "- Retain sections with caution as tissue-context projections; do not use them for disease/control inference or ROI-level localization.",
  paste0("- Primary Step 3 Loop/TAL projection variable: `", primary_metric, "`."),
  "- Supplementary variables: nonzero Loop/TAL indicator and within-section rank percentile, reported as descriptive sensitivity/context only.",
  "- Disease/control comparison should be avoided; the single control/reference section is insufficient for claim-grade testing.",
  "",
  "## Decision",
  primary_decision
)
writeLines(metric_note, file.path(root, "notes/phase3_step2B_step3_metric_decision.md"))

wording <- c(
  "# Phase 3-Step 2B manuscript-safe spatial interpretation wording",
  "",
  "## Safe Results wording",
  "Seurat label transfer projected broad renal papillary compartments onto the locked Visium sections. Loop/TAL prediction scores were sparse across sections; therefore, spatial analyses used the label-transfer output as a supplementary tissue-context projection rather than a claim-grade localization test.",
  "",
  "## Safe Limitations wording",
  "The spatial projection is limited by the absence of plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation, limited control representation, and the sparsity of Loop/TAL prediction scores. Label-transfer scores should not be interpreted as true cell fractions or independent spatial validation.",
  "",
  "## Safe Figure legend wording",
  "Spatial overlays show Seurat label-transfer prediction scores for broad renal papillary compartments. Loop/TAL scores are displayed as sparse contextual projections; no plaque, mineral, lesion, or fibrosis ROI annotation was available.",
  "",
  "## Avoided wording",
  "- Loop/TAL-enriched spots",
  "- spatial validation",
  "- plaque localization",
  "- lesion localization",
  "- causal niche"
)
writeLines(wording, file.path(root, "notes/phase3_step2B_spatial_interpretation_wording.md"))

include_summary <- paste(include_counts$include_for_step3, include_counts$N, sep = "=", collapse = "; ")
report <- c(
  "# Phase 3-Step 2B report",
  "",
  paste0("- Label-transfer files read: ", length(score_files), "."),
  paste0("- Loop/TAL nonzero fraction ranged from ", fmt(range_nonzero[1], 2), "% to ", fmt(range_nonzero[2], 2), "% across sections; mean was ", fmt(mean_nonzero, 2), "%."),
  paste0("- Loop/TAL was the max-predicted compartment for ", loop_max_spots, " in-tissue spots across all audited sections."),
  paste0("- Sparse Loop/TAL sections flagged: ", sparse_sections, " of ", nrow(quality), "."),
  paste0("- Low-confidence sections flagged using >10% low-confidence spots: ", low_conf_sections, " of ", nrow(quality), "."),
  paste0("- Section inclusion flags for Step 3: ", include_summary, "."),
  paste0("- Recommended Step 3 Loop/TAL metric: `", primary_metric, "`."),
  "- Top-quartile Loop/TAL rank is retained only as a descriptive overlay, not as a claim-grade high/low classifier.",
  "- Disease/control comparison should be avoided because control representation is limited.",
  "- No MAGMA spatial module scoring, spatial correlations, RCTD/Cell2location, TWAS/bulk analyses, or manuscript edits were performed.",
  "",
  "## Step 3 readiness",
  if (primary_metric == "not_usable_for_claim_grade_correlation") {
    "Step 3 can proceed only if restricted to broad-compartment/contextual overlays, or after human review accepts this limitation."
  } else {
    "Step 3 can proceed using the selected Loop/TAL metric as a cautious within-section contextual projection."
  },
  "",
  "## Recommended next action",
  if (primary_metric == "not_usable_for_claim_grade_correlation") {
    "C. restrict spatial layer to descriptive overlays only."
  } else {
    "A. proceed to Phase 3-Step 3 using selected Loop/TAL metric."
  }
)
writeLines(report, file.path(root, "notes/phase3_step2B_report.md"))

checklist <- data.table(
  task_id = sprintf("P3S2B-%02d", 1:9),
  task_name = c(
    "main_qc_script_created",
    "loop_tal_usability_table_created",
    "predicted_compartment_composition_created",
    "section_quality_flags_created",
    "loop_tal_qc_figures_created",
    "step3_metric_decision_note_created",
    "manuscript_safe_wording_created",
    "phase3_step2B_report_created",
    "guardrails_observed"
  ),
  completed = "yes",
  output_file = c(
    "scripts/08_spatial_analysis/phase3_step2B_loop_tal_projection_usability_audit.R",
    "results/tables/phase3_step2B_loop_tal_projection_usability.tsv",
    "results/tables/phase3_step2B_predicted_compartment_composition.tsv",
    "results/tables/phase3_step2B_section_quality_flags.tsv",
    "results/figures/phase3_step2B_*.pdf",
    "notes/phase3_step2B_step3_metric_decision.md",
    "notes/phase3_step2B_spatial_interpretation_wording.md",
    "notes/phase3_step2B_report.md",
    "codex_tasks/phase3_step2B_completion_checklist.tsv"
  ),
  blocking_issue = c(rep("none", 8), "none"),
  manual_review_needed = c("no", "no", "no", "no", "no", "yes", "yes", "yes", "no"),
  notes = c(
    "Reads existing label_transfer_scores.tsv.gz files only.",
    "Audits sparsity, ties, nonzero fraction, and metric usability.",
    "Shows whether Loop/TAL is ever the max-predicted compartment.",
    "Flags low-confidence and Loop/TAL sparsity concerns.",
    "Figures emphasize zero inflation and contextual projection.",
    "Human review should accept selected metric before Step 3.",
    "Safe wording is provided for later manuscript revision only.",
    "Report recommends next action without proceeding.",
    "No MAGMA scoring, spatial correlation, RCTD/Cell2location, TWAS/bulk, or manuscript modification."
  )
)
write_tsv(checklist, file.path(root, "codex_tasks/phase3_step2B_completion_checklist.tsv"))

message("Phase 3-Step 2B complete: ", length(score_files), " label-transfer files audited.")
