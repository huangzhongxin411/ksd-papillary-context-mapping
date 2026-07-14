#!/usr/bin/env Rscript

# Phase 3-Step 4: assemble existing spatial evidence only. No new transfer,
# scoring, correlations, disease/control testing, or ROI inference is performed.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

root <- normalizePath(getwd())
fig_dir <- file.path(root, "results", "figures")
tab_dir <- file.path(root, "results", "tables")
note_dir <- file.path(root, "notes")
task_dir <- file.path(root, "codex_tasks")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(note_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_quiet <- function(path) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
write_tsv <- function(x, path) readr::write_tsv(x, path, na = "")

required <- c(
  "notes/phase3_step1_report.md",
  "notes/phase3_step1B_report.md",
  "notes/phase3_step2_report.md",
  "notes/phase3_step2B_report.md",
  "notes/phase3_step3_report.md",
  "config/spatial_sample_metadata_curated_phase3.tsv",
  "results/tables/phase3_step1B_section_inclusion_decision.tsv",
  "results/tables/phase3_step2_label_transfer_summary.tsv",
  "results/tables/phase3_step2B_loop_tal_projection_usability.tsv",
  "results/tables/phase3_step2B_predicted_compartment_composition.tsv",
  "results/tables/phase3_step2B_section_quality_flags.tsv",
  "results/tables/phase3_step3_spatial_magma_module_gene_mapping.tsv",
  "results/tables/phase3_step3_predicted_compartment_module_summary.tsv"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing) > 0L) stop("Missing required Step 4 input(s): ", paste(missing, collapse = "; "))

label_qc <- read_tsv_quiet(file.path(tab_dir, "phase3_step2_label_transfer_summary.tsv"))
loop <- read_tsv_quiet(file.path(tab_dir, "phase3_step2B_loop_tal_projection_usability.tsv"))
quality <- read_tsv_quiet(file.path(tab_dir, "phase3_step2B_section_quality_flags.tsv"))
mapping <- read_tsv_quiet(file.path(tab_dir, "phase3_step3_spatial_magma_module_gene_mapping.tsv"))
inclusion <- read_tsv_quiet(file.path(tab_dir, "phase3_step1B_section_inclusion_decision.tsv"))

n_sections <- nrow(label_qc)
n_complete <- sum(inclusion$complete_visium_input == "complete")
n_success <- sum(label_qc$label_transfer_status == "success")
mean_low_conf <- mean(label_qc$percent_low_confidence_spots)
loop_nonzero_mean <- mean(loop$percent_nonzero_Loop_TAL_score_spots)
loop_max_spots <- sum(loop$n_nonzero_Loop_TAL_score_spots > 0L) # all are zero here; the audited maximum-predicted count is separately documented.
n_human_review <- sum(quality$include_for_step3 == "human_review")
module_min_present <- min(mapping$n_present_in_spatial)
module_max_present <- max(mapping$n_present_in_spatial)

# 1. Integrated evidence summary
evidence <- tibble::tibble(
  evidence_component = c(
    "spatial input completeness", "curated metadata", "no ROI annotation", "Seurat label transfer",
    "prediction-score QC", "Loop/TAL sparsity", "predicted broad-compartment composition",
    "spatial MAGMA module score overlays", "predicted-compartment module summaries",
    "disease/control comparison limitation", "final spatial claim boundary"
  ),
  input_or_analysis = c(
    "Phase 3-Step 1/1B section audit and inclusion table", "Phase 3-Step 1B curated metadata table",
    "Phase 3-Step 1 ROI annotation search", "Phase 3-Step 2 Seurat anchor-based label transfer",
    "Phase 3-Step 2/2B maximum prediction-score and low-confidence QC", "Phase 3-Step 2B Loop/TAL usability audit",
    "Phase 3-Step 2B max-predicted broad-compartment composition", "Phase 3-Step 3 canonical MAGMA module overlays",
    "Phase 3-Step 3 broad predicted-compartment module summaries", "Locked spatial metadata and Phase 3-Step 2/2B reports",
    "Integrated Phase 3 evidence and claim-boundary notes"
  ),
  key_result = c(
    sprintf("%d sections were identified; %d had complete Visium inputs and all were retained.", n_sections, n_complete),
    "Curated section metadata were locked; confidence was high for 4 sections and moderate for 6 sections.",
    "No plaque-, mineral-, lesion-, or fibrosis-resolved spot-level ROI annotation was available.",
    sprintf("Seurat label transfer completed successfully for %d/%d sections using six broad reference compartments.", n_success, n_sections),
    sprintf("Median section-level maximum prediction score was 0.8664; mean low-confidence spot fraction was %.4f%%. %d sections exceeded the >10%% low-confidence flag.", mean_low_conf, n_human_review),
    sprintf("Loop/TAL nonzero-score fraction was %.2f%% across every section; Loop/TAL was max-predicted for 0 in-tissue spots.", loop_nonzero_mean),
    "Collecting duct and injured/undifferentiated epithelial labels dominated most max-predicted spot compositions; sparse stromal contexts occurred in only some sections.",
    sprintf("Five canonical MAGMA modules were scored; present genes per module-section ranged from %d to %d.", module_min_present, module_max_present),
    "Module scores were summarized only across broad label-transfer-predicted epithelial/stromal contexts.",
    "Control representation is limited; no claim-grade disease/control spatial hypothesis test was performed.",
    "Spatial results are supplementary anatomical/tissue-context projections, not independent validation or localized causal evidence."
  ),
  support_strength = c("high for technical completeness", "moderate", "high for absence within available files", "moderate", "moderate", "high for non-usability", "descriptive", "descriptive", "descriptive", "insufficient for inferential comparison", "boundary statement"),
  limitation = c(
    "Completeness does not establish biological comparability.", "Six sections have moderate metadata confidence.",
    "No anatomical lesion/plaque target can be assessed.", "Label transfer is not true deconvolution.",
    "Two sections had >10% low-confidence spots.", "No usable Loop/TAL spatial score variation remains.",
    "Predicted labels are contextual assignments, not cell fractions.", "Module scores are expression summaries, not causal or spatial-validation evidence.",
    "No ROI or claim-grade cell-type enrichment inference.", "Limited controls and heterogeneous section context.",
    "No plaque-, lesion-, causal-niche, or disease/control spatial claim."
  ),
  allowed_claim = c(
    "Ten complete spatial sections provided an available projection resource.", "Metadata-supported descriptive grouping and section traceability.",
    "ROI-specific localization could not be evaluated with the available resource.", "Broad tissue-context projection was feasible.",
    "Prediction confidence can be reported as technical QC.", "Loop/TAL co-distribution was not used for claim-grade analysis.",
    "Broad predicted tissue contexts can be displayed descriptively.", "MAGMA modules can be displayed as descriptive anatomical/tissue-context overlays.",
    "Module scores can be summarized across predicted broad contexts.", "Disease/control status is reported as a limitation only.", "Spatial evidence is supplementary context only."
  ),
  not_allowed_claim = c(
    "Independent spatial validation.", "Disease/control spatial inference.", "Plaque-specific, mineral-specific, lesion-stage, or fibrosis-resolved localization.",
    "True cell fractions or deconvolution.", "Biological validation from prediction confidence.", "Loop/TAL spatial enrichment or correlation.",
    "Cell-type abundance or lesion localization.", "Plaque-specific localization, causal niche, or causal gene evidence.",
    "Claim-grade compartment enrichment.", "Disease/control spatial difference.", "Spatial validation, plaque localization, lesion-stage localization, causal niche."
  ),
  recommended_location = c("Supplementary Figure S-spatial-1/Table", "Supplementary Table", "Supplementary Figure S-spatial-1", "Supplementary Figure S-spatial-1", "Supplementary Figure S-spatial-1", "Supplementary Figure S-spatial-1 and Methods", "Supplementary Figure S-spatial-1", "Supplementary Figure S-spatial-2", "Supplementary Figure S-spatial-2/Table", "Methods and Limitations", "Results, Methods, Limitations, and figure legends")
)
write_tsv(evidence, file.path(tab_dir, "phase3_step4_spatial_integrated_evidence_summary.tsv"))

# Render existing PDFs to temporary PNGs. This only converts prior figure outputs;
# it does not rerun any spatial statistical analysis.
tmp_dir <- file.path(tempdir(), "phase3_step4_pdf_rasters")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
render_pdf <- function(src, key) {
  out <- file.path(tmp_dir, paste0(key, ".png"))
  status <- system2("sips", c("-s", "format", "png", src, "--out", out), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out)) stop("Unable to rasterize source PDF: ", src)
  out
}
as_grob <- function(path) {
  img <- png::readPNG(path)
  grid::rasterGrob(img, interpolate = TRUE)
}
panel <- function(label, title, grob = NULL, text = NULL, subtitle = NULL) {
  multiline_subtitle <- !is.null(subtitle) && grepl("\\n", subtitle, fixed = TRUE)
  base <- ggdraw() +
    draw_label(label, x = 0.015, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 15, colour = "#333333") +
    draw_label(title, x = 0.09, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 11.5, colour = "#245A64")
  if (!is.null(subtitle)) base <- base + draw_label(subtitle, x = 0.09, y = 0.925, hjust = 0, vjust = 1, size = if (multiline_subtitle) 7 else 8.5, lineheight = 1.05, colour = "#555555")
  if (!is.null(grob)) base <- base + draw_grob(grob, x = 0.04, y = 0.03, width = 0.92, height = if (is.null(subtitle)) 0.88 else if (multiline_subtitle) 0.67 else 0.82)
  if (!is.null(text)) base <- base + draw_label(text, x = 0.10, y = 0.52, hjust = 0, vjust = 0.5, size = 13, lineheight = 1.35, colour = "#333333")
  base
}
fig_title <- function(title, subtitle) {
  ggdraw() +
    draw_label(title, x = 0, y = 1, hjust = 0, vjust = 1, fontface = "bold", size = 18, colour = "#245A64") +
    draw_label(subtitle, x = 0, y = 0.60, hjust = 0, vjust = 1, size = 9.5, colour = "#555555")
}
save_figure <- function(plot, stem, width = 13, height = 10.5) {
  # The base PDF device keeps the workflow functional on this macOS setup,
  # where the optional Cairo module depends on an unavailable X11 library.
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

qc_pdf <- render_pdf(file.path(fig_dir, "phase3_step2_prediction_score_distribution.pdf"), "prediction_score_distribution")
composition_pdf <- render_pdf(file.path(fig_dir, "phase3_step2B_predicted_compartment_composition.pdf"), "predicted_compartment_composition")
bonf_pdf <- render_pdf(file.path(fig_dir, "phase3_step3_all_sections_Bonferroni_module_overlay.pdf"), "bonferroni_overlay")
top100_pdf <- render_pdf(file.path(fig_dir, "phase3_step3_all_sections_Top100_module_overlay.pdf"), "top100_overlay")
sugg_pdf <- render_pdf(file.path(fig_dir, "phase3_step3_all_sections_Suggestive_module_overlay.pdf"), "suggestive_overlay")
module_summary_pdf <- render_pdf(file.path(fig_dir, "phase3_step3_module_score_by_predicted_compartment.pdf"), "module_by_compartment")

resource_text <- paste0(
  "10 Visium sections\n", n_complete, " complete inputs\n", n_success, "/", n_sections, " label-transfer successes\n\n",
  "Curated metadata: 4 high-confidence; 6 moderate-confidence sections\n\n",
  "No plaque, mineral, lesion, or fibrosis ROI annotation available"
)
loop_text <- paste0(
  "Loop/TAL usability boundary\n\n",
  "Nonzero prediction-score fraction: 0%\n",
  "Loop/TAL max-predicted in-tissue spots: 0\n\n",
  "Not used for claim-grade MAGMA x Loop/TAL co-distribution.\n",
  "Not evidence of Loop/TAL spatial enrichment."
)

supp1 <- plot_grid(
  fig_title("Supplementary Spatial Figure S-spatial-1 | Spatial resource and label-transfer QC",
            "Supplementary anatomical/tissue-context projection. Seurat label transfer is not true deconvolution."),
  plot_grid(
    panel("A", "Spatial resource summary and ROI boundary", text = resource_text),
    panel("B", "Prediction-score distribution", grob = as_grob(qc_pdf), subtitle = "Existing Phase 3-Step 2 QC output"),
    panel("C", "Predicted broad-compartment composition", grob = as_grob(composition_pdf), subtitle = "Existing Phase 3-Step 2B descriptive summary"),
    panel("D", "Loop/TAL projection usability boundary", text = loop_text),
    ncol = 2, labels = NULL, align = "hv"
  ), ncol = 1, rel_heights = c(0.12, 0.88)
)
save_figure(supp1, "phase3_step4_SuppFig_spatial_label_transfer_QC")

module_boundary <- "Descriptive spatial module projection; no ROI annotation\nNot plaque-specific localization; no disease/control test; not independent validation"
supp2 <- plot_grid(
  fig_title("Supplementary Spatial Figure S-spatial-2 | Descriptive MAGMA spatial module projections",
            "Canonical MAGMA modules displayed as broad tissue-context overlays; no claim-grade Loop/TAL co-distribution."),
  plot_grid(
    panel("A", "Bonferroni module overlay overview", grob = as_grob(bonf_pdf), subtitle = "Descriptive spatial module projection"),
    panel("B", "Top100 module overlay overview", grob = as_grob(top100_pdf), subtitle = "Descriptive spatial module projection"),
    panel("C", "Suggestive module overlay overview", grob = as_grob(sugg_pdf), subtitle = "Descriptive spatial module projection"),
    panel("D", "Module scores by predicted broad compartment", grob = as_grob(module_summary_pdf), subtitle = module_boundary),
    ncol = 2, labels = NULL, align = "hv"
  ), ncol = 1, rel_heights = c(0.12, 0.88)
)
save_figure(supp2, "phase3_step4_SuppFig_spatial_MAGMA_module_overlays")

# 5. Figure/source-data manifest
manifest <- tibble::tribble(
  ~figure, ~panel, ~panel_title, ~source_figure_or_table, ~source_data_file, ~main_message, ~interpretation_boundary, ~manual_polishing_needed, ~notes,
  "S-spatial-1", "A", "Spatial resource summary and ROI boundary", "phase3_step1_report.md; phase3_step1B_section_inclusion_decision.tsv", "results/tables/phase3_step1B_section_inclusion_decision.tsv", "Ten complete sections were available; ROI annotation was absent.", "No plaque/mineral/lesion/fibrosis localization can be claimed.", "no", "Text panel assembled from accepted audits.",
  "S-spatial-1", "B", "Prediction-score distribution", "phase3_step2_prediction_score_distribution.pdf", "source_data/figures/phase3_step2_label_transfer_qc_source_data.tsv", "Prediction-score QC supports technical interpretation of broad context projection.", "Not true deconvolution or biological validation.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production.",
  "S-spatial-1", "C", "Predicted broad-compartment composition", "phase3_step2B_predicted_compartment_composition.pdf", "results/tables/phase3_step2B_predicted_compartment_composition.tsv", "Broad max-predicted compartment composition is descriptive.", "Predicted labels are not cell fractions or disease comparisons.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production.",
  "S-spatial-1", "D", "Loop/TAL projection usability boundary", "phase3_step2B_loop_tal_projection_usability.tsv", "results/tables/phase3_step2B_loop_tal_projection_usability.tsv", "Loop/TAL score sparsity precluded claim-grade co-distribution.", "Do not claim Loop/TAL spatial enrichment or correlation.", "no", "Text panel reports audited 0% nonzero score fraction and zero max-predicted spots.",
  "S-spatial-2", "A", "Bonferroni module overlay overview", "phase3_step3_all_sections_Bonferroni_module_overlay.pdf", "results/tables/phase3_step3_spatial_module_scores.tsv.gz", "Bonferroni module is shown as a descriptive tissue-context overlay.", "Not plaque-specific localization or causal evidence.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production.",
  "S-spatial-2", "B", "Top100 module overlay overview", "phase3_step3_all_sections_Top100_module_overlay.pdf", "results/tables/phase3_step3_spatial_module_scores.tsv.gz", "Top100 module is shown as a descriptive tissue-context overlay.", "Not independent spatial validation or Loop/TAL enrichment.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production.",
  "S-spatial-2", "C", "Suggestive module overlay overview", "phase3_step3_all_sections_Suggestive_module_overlay.pdf", "results/tables/phase3_step3_spatial_module_scores.tsv.gz", "Suggestive module is shown as a descriptive tissue-context overlay.", "Not lesion-stage localization or disease/control evidence.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production.",
  "S-spatial-2", "D", "Module scores by predicted broad compartment", "phase3_step3_module_score_by_predicted_compartment.pdf", "results/tables/phase3_step3_predicted_compartment_module_summary.tsv", "Module scores can be summarized across broad predicted contexts.", "Predicted compartments are contextual labels, not true fractions or ROIs.", "yes", "Existing PDF rasterized for assembly; retain source PDF for final production."
)
write_tsv(manifest, file.path(tab_dir, "phase3_step4_spatial_figure_manifest.tsv"))

# Scientific figure design QC required for the assembled figure artifacts.
qc <- tibble::tribble(
  ~figure, ~check, ~status, ~evidence, ~manual_action,
  "S-spatial-1", "PDF and 600 dpi PNG exist", ifelse(all(file.exists(file.path(fig_dir, paste0("phase3_step4_SuppFig_spatial_label_transfer_QC.", c("pdf", "png"))))), "pass", "fail"), "Generated by Step 4 assembly script.", "none",
  "S-spatial-1", "panel labels and boundary statement", "pass", "Panels A-D include clear labels; panel A/D state no ROI and Loop/TAL non-usability.", "none",
  "S-spatial-1", "source-panel rendering", "manual review", "Existing source PDFs were rasterized for composite assembly.", "Replace embedded rasters with vector panels during final journal production if required.",
  "S-spatial-2", "PDF and 600 dpi PNG exist", ifelse(all(file.exists(file.path(fig_dir, paste0("phase3_step4_SuppFig_spatial_MAGMA_module_overlays.", c("pdf", "png"))))), "pass", "fail"), "Generated by Step 4 assembly script.", "none",
  "S-spatial-2", "claim-boundary wording", "pass", "Title/subtitle and panel D mark descriptive projection, no ROI, no plaque-specific localization.", "none",
  "S-spatial-2", "source-panel rendering", "manual review", "Existing source PDFs were rasterized for composite assembly.", "Replace embedded rasters with vector panels during final journal production if required.",
  "Both", "PNG readability at 50% scale", "manual review", "600 dpi PNGs generated; visual QA recorded after render.", "Confirm exact journal column width after manuscript layout is locked.")
write_tsv(qc, file.path(tab_dir, "phase3_step4_spatial_figure_qc.tsv"))

writeLines(c(
  "# Phase 3-Step 4 Spatial Figure Routing Plan", "",
  "## Main manuscript", "",
  "Do not use the spatial layer as a strong main-evidence figure. If spatial context is retained in a later integrated figure, use only a small boundary/context panel stating that it is supplementary broad tissue-context projection, not deconvolution or spatial validation.", "",
  "## Supplementary Figure S-spatial-1", "",
  "Spatial input and label-transfer QC: (A) 10-section resource and no-ROI boundary; (B) prediction-score distribution; (C) predicted broad-compartment composition; (D) Loop/TAL usability boundary. This figure documents why Loop/TAL co-distribution was not advanced to an inferential claim.", "",
  "## Supplementary Figure S-spatial-2", "",
  "Descriptive MAGMA spatial module projections: (A) Bonferroni; (B) Top100; (C) Suggestive; (D) module scores by predicted broad compartment. Use the explicit labels 'descriptive spatial module projection', 'no ROI annotation', and 'not plaque-specific localization'.", "",
  "## Supplementary tables", "",
  "Retain the curated section metadata, section inclusion decision, label-transfer QC, Loop/TAL usability audit, predicted-compartment composition, module gene mapping, and broad-context summaries as traceable supplementary tables.", "",
  "## Production note", "",
  "The composite PDFs use rasterized copies of existing source PDFs because vector-PDF embedding was unavailable in this environment. The source PDFs remain the authoritative plot files; perform a final vector-panel substitution if journal production requires it."
), file.path(note_dir, "phase3_step4_spatial_figure_routing_plan.md"))

writeLines(c(
  "# Phase 3-Step 4 Results Wording", "",
  "To provide supplementary spatial tissue-context information, we analyzed ten Visium sections with complete expression, position, scale-factor, and tissue-image inputs. Curated section metadata were available, although confidence was moderate for six sections. No plaque-, mineral-, lesion-, or fibrosis-resolved spot-level ROI annotation was available; consequently, the spatial resource was used for broad anatomical/tissue-context projection only.", "",
  "Using Seurat anchor-based label transfer from the single-nucleus reference, all ten sections received broad-compartment prediction scores. Section-level maximum prediction scores supported technical use of these projections as contextual annotations, while two sections exceeded the pre-specified low-confidence-spot flag and were retained only with caution. These transferred labels are not estimates of true cell fractions.", "",
  "Loop/TAL prediction scores were zero across the audited in-tissue spots in every section, and Loop/TAL was not the maximum-predicted compartment for any in-tissue spot. We therefore did not perform or interpret MAGMA x Loop/TAL spatial co-distribution as a claim-grade analysis.", "",
  "Canonical MAGMA gene modules were displayed as descriptive spatial expression overlays and summarized across broad predicted epithelial/stromal contexts. These supplementary displays provide tissue-context visualization only and do not establish spatial validation, plaque localization, lesion-stage enrichment, or causal spatial niches."
), file.path(note_dir, "phase3_step4_results_wording.md"))

writeLines(c(
  "# Phase 3-Step 4 Methods Wording", "",
  "Ten complete Visium sections were included after review of matrix, spatial-coordinate, scale-factor, and tissue-image inputs. Section-level metadata were locked before analysis; no plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation was available.", "",
  "Broad spatial tissue-context projection used Seurat anchor-based label transfer from the audited GSE231569 single-nucleus reference. Reference labels were the `phase1_cell_type` broad compartments: Collecting_duct_principal, Endothelial, Fibroblast_stromal, Injured_undifferentiated_epithelial, Loop_of_Henle_TAL, and Pericyte_smooth_muscle. We summarized section-level maximum prediction scores and the proportion of spots with maximum score below 0.5 as technical QC. Transferred labels and scores were treated as contextual predictions, not cell fractions or formal deconvolution estimates.", "",
  "Loop/TAL usability was audited from the transferred Loop/TAL prediction score across each section. Because the score was zero/nonzero-sparse across the full audited resource and Loop/TAL was not max-predicted for any in-tissue spot, no claim-grade Loop/TAL co-distribution or correlation analysis was performed.", "",
  "For the previously defined canonical MAGMA modules, per-spot module scores were calculated as the arithmetic mean of log-normalized expression across module genes present in each section. Module scores were displayed as descriptive overlays and summarized by broad predicted compartment. No ROI-resolved analysis and no disease/control spatial hypothesis testing were performed."
), file.path(note_dir, "phase3_step4_methods_wording.md"))

writeLines(c(
  "# Phase 3-Step 4 Limitations Wording", "",
  "The spatial analysis is supplementary and should be interpreted as broad anatomical/tissue-context projection. The available sections lacked plaque-, mineral-, lesion-, and fibrosis-resolved ROI annotation, precluding plaque-specific or lesion-stage localization. Seurat label transfer provides predicted context rather than true deconvolution or cell-fraction estimation. Loop/TAL prediction scores were zero/sparse across the audited sections, so Loop/TAL spatial co-distribution was not interpretable as a claim-grade analysis. Control representation was limited and no disease/control spatial hypothesis test was conducted. Finally, spatial MAGMA module overlays are descriptive expression summaries and do not constitute independent spatial validation, causal localization, or causal gene evidence."
), file.path(note_dir, "phase3_step4_limitations_wording.md"))

writeLines(c(
  "# Supplementary Spatial Figure Legends", "",
  "## Supplementary Spatial Figure S-spatial-1 | Spatial resource and broad-compartment label-transfer quality control.", "",
  "(A) Spatial resource summary for ten Visium sections with complete inputs, with the explicit boundary that plaque-, mineral-, lesion-, and fibrosis-resolved ROI annotations were unavailable. (B) Distribution of maximum Seurat label-transfer prediction scores across sections. (C) Composition of max-predicted broad compartments by section. (D) Loop/TAL projection usability boundary: the Loop/TAL prediction score had a 0% nonzero fraction across audited sections and Loop/TAL was not max-predicted for any in-tissue spot. Seurat label transfer is a broad tissue-context projection, not true deconvolution. Accordingly, this figure does not support Loop/TAL spatial enrichment, plaque-specific localization, lesion-stage localization, independent spatial validation, or disease/control spatial inference.", "",
  "## Supplementary Spatial Figure S-spatial-2 | Descriptive spatial projections of canonical MAGMA modules.", "",
  "(A-C) All-section overviews of Bonferroni, Top100, and suggestive MAGMA module scores, respectively. (D) Module scores summarized across broad Seurat label-transfer-predicted compartments. Module scores are arithmetic means of log-normalized expression across module genes present in each section. These panels are descriptive spatial module projections with no ROI annotation and are not plaque-specific localization. They do not provide claim-grade Loop/TAL co-distribution, lesion-stage localization, independent spatial validation, causal spatial-niche evidence, or disease/control spatial differences."
), file.path(note_dir, "phase3_step4_supplementary_spatial_figure_legends.md"))

reviewer <- tibble::tribble(
  ~reviewer_concern, ~action_taken, ~result, ~remaining_limitation, ~manuscript_or_supplement_location,
  "need for spatial analysis beyond correlation", "Completed resource audit, Seurat broad-compartment projection, QC, and descriptive MAGMA overlays.", "Spatial layer provides supplementary tissue-context visualization.", "No claim-grade spatial association test or spatial validation.", "Supplementary Figures S-spatial-1 and S-spatial-2; Methods/Limitations wording",
  "label transfer feasibility", "Locked 10 complete sections and used Seurat anchor-based transfer because RCTD/spacexr was unavailable.", "All 10 sections completed label transfer.", "Predictions are not true deconvolution or fractions.", "Supplementary Figure S-spatial-1; Methods wording",
  "no ROI annotation", "Audited spatial inputs and documented absence of plaque/mineral/lesion/fibrosis ROI annotation.", "ROI-specific analyses were not performed.", "Plaque- and lesion-specific localization cannot be assessed.", "Supplementary Figure S-spatial-1A; Limitations wording",
  "no plaque-specific claim", "Added synchronized claim boundaries to figures, legends, Results, Methods, and Limitations.", "Spatial context retained as supplementary only.", "No plaque-specific evidence exists.", "All Step 4 narrative outputs and both supplementary legends",
  "Loop/TAL prediction-score sparsity", "Audited continuous and rank usability across all sections.", "0% nonzero fraction; zero max-predicted Loop/TAL spots; co-distribution not performed.", "Loop/TAL spatial enrichment cannot be inferred.", "Supplementary Figure S-spatial-1D; Results/Methods wording",
  "descriptive module overlays", "Assembled canonical MAGMA overlay and broad-context summary panels.", "Overlays are available as descriptive tissue-context displays.", "No causal, ROI-specific, or independent-validation interpretation.", "Supplementary Figure S-spatial-2; Results/Limitations wording",
  "disease/control limitation", "Preserved condition labels only as metadata; did not run a hypothesis test.", "No disease/control spatial claim is made.", "Limited control representation prevents robust comparison.", "Integrated evidence table; Methods/Limitations wording",
  "source-data traceability", "Created panel-level figure/source-data manifest and QC table.", "Each panel maps to an accepted source plot/table and source-data file.", "Composite panels use rasterized existing PDFs pending optional vector substitution.", "results/tables/phase3_step4_spatial_figure_manifest.tsv"
)
write_tsv(reviewer, file.path(task_dir, "phase3_step4_reviewer_response_spatial_issue_table.tsv"))

writeLines(c(
  "# Phase 3-Step 4 Closure Report", "",
  "## What Phase 3 accomplished", "",
  "Phase 3 locked and audited ten complete Visium sections; curated spatial metadata; documented the absence of plaque/mineral/lesion/fibrosis ROI annotation; completed Seurat broad-compartment label transfer; audited prediction confidence and Loop/TAL usability; calculated canonical MAGMA spatial module scores; and assembled the accepted spatial evidence as two supplementary figures with source traceability.", "",
  "## Accepted outputs", "",
  "Accepted outputs include the Step 1/1B input and metadata audits, Step 2 label-transfer summaries and QC, the Step 2B Loop/TAL non-usability audit, Step 3 canonical module mapping/overlays/broad-context summaries, and the present Step 4 integrated evidence table, supplementary figures, legends, wording, and reviewer-response table.", "",
  "## Remaining limitations", "",
  "No plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation was available. Seurat label transfer is not true deconvolution. Loop/TAL prediction scores were zero/sparse across the resource, so no claim-grade MAGMA x Loop/TAL spatial co-distribution was performed. Control representation is limited, and no disease/control spatial hypothesis test was conducted. The assembled composite figures rasterize existing source PDFs and should receive optional vector-panel substitution at final journal production.", "",
  "## Figure placement", "",
  "Place the spatial results in the supplement as S-spatial-1 (resource and label-transfer QC) and S-spatial-2 (descriptive MAGMA module projections). Do not use spatial analyses as strong main-text evidence; at most retain a small boundary/context panel in a later integrated figure.", "",
  "## Final safe spatial claim", "",
  "The available Visium resource supports supplementary descriptive anatomical/tissue-context projection of canonical MAGMA modules across broad label-transfer-predicted compartments.", "",
  "## Final unsafe spatial claims", "",
  "Do not claim Loop/TAL spatial enrichment, plaque-specific localization, lesion-stage localization, causal spatial niche, independent spatial validation, true deconvolution/cell fractions, or disease/control spatial differences.", "",
  "## Closure decision", "",
  "Phase 3 can close with the spatial layer retained as supplementary descriptive context.", "",
  "## Recommended next step", "",
  "A. close Phase 3 and proceed to Phase 4: TWAS downgrading and candidate evidence-model repair. This choice best matches the accepted spatial limitations and avoids treating the spatial layer as validation."
), file.path(note_dir, "phase3_step4_phase3_closure_report.md"))

checklist <- tibble::tribble(
  ~task_id, ~task_name, ~completed, ~output_file, ~blocking_issue, ~manual_review_needed, ~notes,
  "P3S4-01", "Integrated spatial evidence summary", "yes", "results/tables/phase3_step4_spatial_integrated_evidence_summary.tsv", "none", "no", "Includes all required evidence components.",
  "P3S4-02", "Spatial figure routing plan", "yes", "notes/phase3_step4_spatial_figure_routing_plan.md", "none", "no", "Routes spatial evidence to supplement.",
  "P3S4-03", "Supplementary Spatial Figure S-spatial-1 PDF/PNG", "yes", "results/figures/phase3_step4_SuppFig_spatial_label_transfer_QC.pdf; results/figures/phase3_step4_SuppFig_spatial_label_transfer_QC.png", "none", "yes", "Optional vector-panel substitution at final production.",
  "P3S4-04", "Supplementary Spatial Figure S-spatial-2 PDF/PNG", "yes", "results/figures/phase3_step4_SuppFig_spatial_MAGMA_module_overlays.pdf; results/figures/phase3_step4_SuppFig_spatial_MAGMA_module_overlays.png", "none", "yes", "Optional vector-panel substitution at final production.",
  "P3S4-05", "Figure/source-data manifest", "yes", "results/tables/phase3_step4_spatial_figure_manifest.tsv", "none", "no", "All panels mapped to source figures/tables.",
  "P3S4-06", "Figure visual QC", "yes", "results/tables/phase3_step4_spatial_figure_qc.tsv", "none", "yes", "Check final journal-scale raster/text rendering.",
  "P3S4-07", "Results wording", "yes", "notes/phase3_step4_results_wording.md", "none", "no", "Conservative claim boundary.",
  "P3S4-08", "Methods wording", "yes", "notes/phase3_step4_methods_wording.md", "none", "no", "No new analysis described.",
  "P3S4-09", "Limitations wording", "yes", "notes/phase3_step4_limitations_wording.md", "none", "no", "Documents ROI, deconvolution, Loop/TAL, and control limits.",
  "P3S4-10", "Supplementary figure legends", "yes", "notes/phase3_step4_supplementary_spatial_figure_legends.md", "none", "no", "Both legends contain claim boundaries.",
  "P3S4-11", "Reviewer-response spatial issue table", "yes", "codex_tasks/phase3_step4_reviewer_response_spatial_issue_table.tsv", "none", "no", "Eight required reviewer issues included.",
  "P3S4-12", "Phase 3 closure report", "yes", "notes/phase3_step4_phase3_closure_report.md", "none", "no", "Recommends Phase 3 closure and Phase 4 option A.",
  "P3S4-13", "Stop rule", "yes", "N/A", "none", "no", "No new label transfer, scoring, correlation, RCTD/Cell2location, TWAS, bulk analysis, or manuscript DOCX edit performed.")
write_tsv(checklist, file.path(task_dir, "phase3_step4_completion_checklist.tsv"))

message("Phase 3-Step 4 outputs written successfully.")
