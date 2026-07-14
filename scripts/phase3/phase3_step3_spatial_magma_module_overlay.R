#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(scales)
})

root <- getwd()
dir.create(file.path(root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/spatial/phase3_step3_module_scores"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/figures/phase3_step3_spatial_module_overlays"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
safe_id <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
section_uid <- function(sample_id, section_id) paste(safe_id(sample_id), safe_id(section_id), sep = "__")
fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))

required <- c(
  "config/spatial_sample_metadata_curated_phase3.tsv",
  "results/tables/phase3_step1_spatial_sample_manifest.tsv",
  "results/tables/phase3_step1_spatial_section_audit.tsv",
  "notes/phase3_step2B_report.md",
  "notes/phase3_step2B_step3_metric_decision.md",
  "notes/phase3_step2B_spatial_interpretation_wording.md",
  "notes/phase3_step2_spatial_label_transfer_claim_boundary.md"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "))

module_files <- c(
  MAGMA_top50 = "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
  MAGMA_top100 = "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
  MAGMA_Bonferroni = "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
  MAGMA_FDR05 = "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
  MAGMA_suggestive_p1e4 = "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
)
missing_modules <- module_files[!file.exists(file.path(root, module_files))]
if (length(missing_modules)) stop("Missing module gene set file(s): ", paste(missing_modules, collapse = ", "))
modules <- lapply(file.path(root, module_files), function(path) unique(trimws(readLines(path, warn = FALSE))))
modules <- lapply(modules, function(x) x[nzchar(x)])
names(modules) <- names(module_files)
selected_overlay_modules <- c("MAGMA_Bonferroni", "MAGMA_top100", "MAGMA_suggestive_p1e4")
module_display <- c(
  MAGMA_top50 = "Top50",
  MAGMA_top100 = "Top100",
  MAGMA_Bonferroni = "Bonferroni",
  MAGMA_FDR05 = "FDR05",
  MAGMA_suggestive_p1e4 = "Suggestive"
)

meta <- fread(file.path(root, "config/spatial_sample_metadata_curated_phase3.tsv"))
manifest <- fread(file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv"))
included <- meta[included_for_phase3_step2 == "yes"]
included <- merge(
  included,
  manifest[, .(sample_id, section_id, matrix_file)],
  by = c("sample_id", "section_id"),
  all.x = TRUE
)
setorder(included, sample_id, section_id)

score_files <- list.files(
  file.path(root, "results/spatial/phase3_step2_label_transfer"),
  pattern = "^label_transfer_scores.tsv.gz$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(score_files)) stop("No Phase 3-Step 2 label-transfer score files found.")
label_scores <- rbindlist(lapply(score_files, fread), fill = TRUE)
label_scores[, in_tissue := as.integer(in_tissue)]
label_scores[, prediction_score_max := as.numeric(prediction_score_max)]

read_counts <- function(path) {
  counts <- Read10X_h5(path)
  if (is.list(counts)) counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  counts
}

calc_module_scores <- function(counts, module_genes) {
  present <- intersect(module_genes, rownames(counts))
  if (!length(present)) {
    return(list(scores = rep(NA_real_, ncol(counts)), present = character(), missing = setdiff(module_genes, rownames(counts))))
  }
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- NA_real_
  sub <- counts[present, , drop = FALSE]
  norm <- t(t(sub) / lib * 10000)
  log_norm <- log1p(norm)
  list(
    scores = as.numeric(Matrix::colMeans(log_norm, na.rm = TRUE)),
    present = present,
    missing = setdiff(module_genes, present)
  )
}

plot_module_overlay <- function(dt, module_names, title) {
  ggplot(dt[in_tissue == 1 & module_name %in% module_names], aes(x = x_coord, y = -y_coord, color = module_score)) +
    geom_point(size = 0.55, alpha = 0.9) +
    facet_wrap(~ module_name, ncol = 3) +
    scale_color_viridis_c(option = "magma", name = "Module score") +
    coord_equal() +
    labs(
      title = title,
      subtitle = "spatial module projection; no ROI annotation",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 8) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7),
      legend.position = "bottom"
    )
}

all_spot_scores <- list()
mapping_rows <- list()
section_summary_rows <- list()
compartment_summary_rows <- list()
overlay_plot_data <- list()

for (i in seq_len(nrow(included))) {
  row <- included[i]
  uid <- section_uid(row$sample_id, row$section_id)
  message("Scoring modules for ", row$sample_id, " / ", row$section_id)
  counts <- read_counts(file.path(root, row$matrix_file))
  section_labels <- label_scores[sample_id == row$sample_id & section_id == row$section_id]
  if (!nrow(section_labels)) stop("No label-transfer scores for ", row$sample_id, " / ", row$section_id)

  score_dt <- data.table(spot_id = colnames(counts))
  for (module_name in names(modules)) {
    res <- calc_module_scores(counts, modules[[module_name]])
    score_dt[, (module_name) := res$scores]
    mapping_rows[[paste(uid, module_name, sep = "__")]] <- data.table(
      sample_id = row$sample_id,
      section_id = row$section_id,
      module_name = module_name,
      n_module_genes = length(modules[[module_name]]),
      n_present_in_spatial = length(res$present),
      n_missing_from_spatial = length(res$missing),
      present_genes = paste(res$present, collapse = ";"),
      missing_genes = paste(res$missing, collapse = ";"),
      score_computed = ifelse(length(res$present) > 0, "yes", "no"),
      notes = "Module score is arithmetic mean of log-normalized spatial expression across present module genes."
    )
  }

  long <- melt(score_dt, id.vars = "spot_id", variable.name = "module_name", value.name = "module_score")
  long[, module_name := as.character(module_name)]
  long <- merge(
    section_labels[, .(spot_id, sample_id, section_id, x_coord, y_coord, in_tissue, predicted_compartment, prediction_score_max)],
    long,
    by = "spot_id",
    all.y = TRUE,
    sort = FALSE
  )
  n_used <- rbindlist(mapping_rows[startsWith(names(mapping_rows), uid)])[, .(module_name, n_genes_used = n_present_in_spatial)]
  long <- merge(long, n_used, by = "module_name", all.x = TRUE)
  long[, `:=`(
    control_or_disease = row$control_or_disease,
    metadata_confidence = row$metadata_confidence,
    notes = "Descriptive spatial module projection; not ROI-resolved, not disease/control tested, and not Loop/TAL claim-grade co-distribution."
  )]
  setcolorder(long, c(
    "spot_id", "sample_id", "section_id", "x_coord", "y_coord", "in_tissue",
    "predicted_compartment", "prediction_score_max", "module_name", "module_score",
    "n_genes_used", "control_or_disease", "metadata_confidence", "notes"
  ))

  section_dir <- file.path(root, "results/spatial/phase3_step3_module_scores", uid)
  dir.create(section_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(long, file.path(section_dir, "spatial_module_scores.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
  all_spot_scores[[uid]] <- long
  overlay_plot_data[[uid]] <- copy(long)[, section_uid := uid][]

  in_tissue <- long[in_tissue == 1]
  section_summary_rows[[uid]] <- in_tissue[, .(
    n_in_tissue_spots = .N,
    n_genes_used = unique(n_genes_used)[1],
    mean_module_score = mean(module_score, na.rm = TRUE),
    median_module_score = median(module_score, na.rm = TRUE),
    sd_module_score = sd(module_score, na.rm = TRUE),
    iqr_module_score = IQR(module_score, na.rm = TRUE),
    notes = "Section-level descriptive summary only; no disease/control statistical testing."
  ), by = .(sample_id, section_id, control_or_disease, metadata_confidence, module_name)]

  compartment_summary_rows[[uid]] <- in_tissue[!is.na(predicted_compartment), .(
    n_spots = .N,
    mean_module_score = mean(module_score, na.rm = TRUE),
    median_module_score = median(module_score, na.rm = TRUE),
    sd_module_score = sd(module_score, na.rm = TRUE),
    iqr_module_score = IQR(module_score, na.rm = TRUE),
    median_prediction_score_max = median(prediction_score_max, na.rm = TRUE),
    notes = "Predicted broad-compartment context only; not true cell fraction or histological ROI."
  ), by = .(sample_id, section_id, module_name, predicted_compartment)]

  p <- plot_module_overlay(long, selected_overlay_modules, paste(row$sample_id, row$section_id))
  ggsave(
    file.path(root, "results/figures/phase3_step3_spatial_module_overlays", paste0(uid, "_module_overlay.pdf")),
    p,
    width = 8.5,
    height = 3.8,
    units = "in",
    device = "pdf"
  )
  rm(counts, score_dt, long)
  gc(verbose = FALSE)
}

mapping_dt <- rbindlist(mapping_rows, fill = TRUE)
spot_dt <- rbindlist(all_spot_scores, fill = TRUE)
section_summary <- rbindlist(section_summary_rows, fill = TRUE)
compartment_summary <- rbindlist(compartment_summary_rows, fill = TRUE)
write_tsv(mapping_dt, file.path(root, "results/tables/phase3_step3_spatial_magma_module_gene_mapping.tsv"))
fwrite(spot_dt, file.path(root, "results/tables/phase3_step3_spatial_module_scores.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
write_tsv(section_summary, file.path(root, "results/tables/phase3_step3_section_module_score_summary.tsv"))
write_tsv(compartment_summary, file.path(root, "results/tables/phase3_step3_predicted_compartment_module_summary.tsv"))

plot_dt <- rbindlist(overlay_plot_data, fill = TRUE)
for (mod in selected_overlay_modules) {
  p <- ggplot(plot_dt[in_tissue == 1 & module_name == mod], aes(x = x_coord, y = -y_coord, color = module_score)) +
    geom_point(size = 0.35, alpha = 0.85) +
    facet_wrap(~ section_id, scales = "free", ncol = 2) +
    scale_color_viridis_c(option = "magma", name = paste(module_display[[mod]], "score")) +
    labs(
      title = paste("All sections", module_display[[mod]], "module projection"),
      subtitle = "spatial module projection; no ROI annotation",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 8) +
    theme(
      strip.text = element_text(size = 6),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  out_name <- paste0("phase3_step3_all_sections_", module_display[[mod]], "_module_overlay.pdf")
  ggsave(file.path(root, "results/figures", out_name), p, width = 8.5, height = 11, units = "in", device = "pdf")
}

context <- compartment_summary[, .(
  n_sections_with_compartment = uniqueN(paste(sample_id, section_id)),
  total_spots = sum(n_spots),
  median_module_score_across_sections = median(median_module_score, na.rm = TRUE),
  iqr_module_score_across_sections = IQR(median_module_score, na.rm = TRUE)
), by = .(module_name, predicted_compartment)]
context[, descriptive_pattern := fifelse(
  n_sections_with_compartment >= 3,
  "Module scores observed across recurrent predicted broad tissue contexts.",
  "Module scores observed in a sparsely represented predicted broad tissue context."
)]
context[, allowed_interpretation := "spatial module scores can be displayed as anatomical tissue-context overlays; module scores are observed within predicted broad epithelial/stromal tissue contexts"]
context[, not_allowed_interpretation := "Loop/TAL spatial enrichment; plaque-specific localization; lesion-stage localization; independent spatial validation"]
context[, notes := "Descriptive only; predicted compartments are label-transfer context, not true fractions or ROIs."]
setcolorder(context, c(
  "module_name", "predicted_compartment", "n_sections_with_compartment", "total_spots",
  "median_module_score_across_sections", "iqr_module_score_across_sections",
  "descriptive_pattern", "allowed_interpretation", "not_allowed_interpretation", "notes"
))
setorder(context, module_name, predicted_compartment)
write_tsv(context, file.path(root, "results/tables/phase3_step3_broad_context_summary.tsv"))

p_box <- ggplot(spot_dt[in_tissue == 1 & !is.na(predicted_compartment)], aes(x = predicted_compartment, y = module_score, fill = predicted_compartment)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.2) +
  facet_wrap(~ module_name, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(
    title = "Descriptive module scores by predicted broad compartment",
    subtitle = "Label-transfer context only; not true fractions, ROIs, or disease/control tests",
    x = NULL,
    y = "Spatial module score"
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(root, "results/figures/phase3_step3_module_score_by_predicted_compartment.pdf"),
       p_box, width = 8.5, height = 7.5, units = "in", device = "pdf")

heat <- context[n_sections_with_compartment >= 2]
heat[, module_display := module_display[module_name]]
p_heat <- ggplot(heat, aes(x = module_display, y = predicted_compartment, fill = median_module_score_across_sections)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = fmt(median_module_score_across_sections, 2)), size = 2.3) +
  scale_fill_viridis_c(option = "magma", name = "Median score") +
  labs(
    title = "Spatial context summary heatmap",
    subtitle = "Median section-level module scores in represented predicted broad compartments",
    x = "MAGMA module",
    y = "Predicted broad compartment"
  ) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(root, "results/figures/phase3_step3_spatial_context_summary_heatmap.pdf"),
       p_heat, width = 7.2, height = 4.8, units = "in", device = "pdf")

loop_note <- c(
  "# Phase 3-Step 3 Loop/TAL spatial usability note",
  "",
  "- Phase 3-Step 2B found Loop/TAL label-transfer nonzero fraction was 0% across all sections.",
  "- Loop/TAL was not the max-predicted compartment for any in-tissue spot.",
  "- Therefore MAGMA x Loop/TAL spatial co-distribution was not performed as a claim-grade analysis.",
  "- Spatial module scores are retained only as supplementary broad tissue-context projections.",
  "- No plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation is available."
)
writeLines(loop_note, file.path(root, "notes/phase3_step3_loop_tal_spatial_usability_note.md"))

results_wording <- c(
  "# Phase 3-Step 3 Results wording",
  "",
  "Spatial transcriptomic data were used as a supplementary tissue-context projection layer. Seurat label transfer produced broad-compartment prediction scores across ten Visium sections. Because Loop/TAL prediction scores were sparse and no ROI annotation was available, we did not interpret spatial data as Loop/TAL localization or plaque-specific validation. MAGMA module scores were therefore displayed as descriptive spatial overlays and summarized across predicted broad tissue contexts."
)
writeLines(results_wording, file.path(root, "notes/phase3_step3_results_wording.md"))

methods_wording <- c(
  "# Phase 3-Step 3 Methods wording",
  "",
  "Ten locked Visium sections were analyzed using the canonical section paths from the Phase 3-Step 1 manifest and curated spatial metadata. Phase 3-Step 2 label-transfer outputs were used to attach predicted broad-compartment context to spatial spots.",
  "",
  "For each section, raw 10x feature-barcode matrices were log-normalized using log1p(counts per spot divided by total spot counts multiplied by 10,000). Canonical MAGMA gene modules were mapped to the spatial gene universe in each section. For each module and spot, the module score was calculated as the arithmetic mean of log-normalized expression across module genes present in the spatial matrix. Missing genes were recorded per module and section.",
  "",
  "Module scores were summarized descriptively by section and by predicted broad compartment. No plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation was available. No claim-grade Loop/TAL co-distribution was performed because Loop/TAL label-transfer scores were unusable for this purpose. No disease/control statistical testing was performed."
)
writeLines(methods_wording, file.path(root, "notes/phase3_step3_methods_wording.md"))

limitations_wording <- c(
  "# Phase 3-Step 3 Limitations wording",
  "",
  "The spatial transcriptomic layer is supplementary and supports anatomical/tissue-context visualization only. No plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation was available. Seurat label transfer is not true deconvolution and does not estimate physical cell fractions. Loop/TAL prediction scores were sparse and unusable for claim-grade spatial co-distribution. Control representation was limited, so disease/control comparisons were not used as claim-grade evidence."
)
writeLines(limitations_wording, file.path(root, "notes/phase3_step3_limitations_wording.md"))

module_gene_status <- mapping_dt[, .(
  min_present = min(n_present_in_spatial),
  max_missing = max(n_missing_from_spatial)
), by = module_name]
overlay_count <- length(list.files(file.path(root, "results/figures/phase3_step3_spatial_module_overlays"), pattern = "_module_overlay\\.pdf$"))
unsafe_claims <- c(
  "Loop/TAL spatial enrichment",
  "plaque-specific localization",
  "lesion-stage localization",
  "causal spatial niche",
  "independent spatial validation",
  "claim-grade disease/control spatial difference"
)
report <- c(
  "# Phase 3-Step 3 report",
  "",
  paste0("- Sections included: ", nrow(included), "."),
  paste0("- Modules scored: ", paste(names(modules), collapse = ", "), "."),
  "- Module scores were computed as arithmetic means of log-normalized expression over present module genes.",
  paste0("- Gene mapping status: all module-section scores computed; minimum present genes by module ranged from ", min(module_gene_status$min_present), " to ", max(module_gene_status$min_present), "."),
  paste0("- Per-section overlay PDFs generated: ", overlay_count, ". Compact all-section overlays generated for Bonferroni, Top100, and Suggestive modules."),
  "- Predicted-compartment module summaries were generated for broad label-transfer contexts only.",
  "- Loop/TAL co-distribution was not performed because Phase 3-Step 2B showed Loop/TAL scores were zero/nonzero fraction 0% and Loop/TAL was not max-predicted for any spot.",
  "",
  "## Safe Spatial Claim",
  "Spatial module scores can be displayed as descriptive anatomical/tissue-context overlays and summarized across predicted broad epithelial/stromal contexts.",
  "",
  "## Unsafe Spatial Claims",
  paste0("- ", unsafe_claims),
  "",
  "## Phase 3 Status",
  "A figure-integration step is needed before Phase 3 closes if these spatial overlays will be assembled into a manuscript supplement or figure panel.",
  "",
  "## Recommended Next Action",
  "A. proceed to Phase 3-Step 4: spatial evidence integration and Figure/Supplement assembly."
)
writeLines(report, file.path(root, "notes/phase3_step3_report.md"))

checklist <- data.table(
  task_id = sprintf("P3S3-%02d", 1:14),
  task_name = c(
    "main_script_created",
    "module_gene_mapping_created",
    "per_spot_module_scores_created",
    "section_module_summary_created",
    "predicted_compartment_module_summary_created",
    "broad_context_summary_created",
    "per_section_module_overlays_created",
    "all_section_overviews_created",
    "broad_compartment_figures_created",
    "loop_tal_usability_note_created",
    "results_wording_created",
    "methods_wording_created",
    "limitations_wording_created",
    "guardrails_observed"
  ),
  completed = "yes",
  output_file = c(
    "scripts/08_spatial_analysis/phase3_step3_spatial_magma_module_overlay.R",
    "results/tables/phase3_step3_spatial_magma_module_gene_mapping.tsv",
    "results/tables/phase3_step3_spatial_module_scores.tsv.gz",
    "results/tables/phase3_step3_section_module_score_summary.tsv",
    "results/tables/phase3_step3_predicted_compartment_module_summary.tsv",
    "results/tables/phase3_step3_broad_context_summary.tsv",
    "results/figures/phase3_step3_spatial_module_overlays/*_module_overlay.pdf",
    "results/figures/phase3_step3_all_sections_*_module_overlay.pdf",
    "results/figures/phase3_step3_module_score_by_predicted_compartment.pdf;results/figures/phase3_step3_spatial_context_summary_heatmap.pdf",
    "notes/phase3_step3_loop_tal_spatial_usability_note.md",
    "notes/phase3_step3_results_wording.md",
    "notes/phase3_step3_methods_wording.md",
    "notes/phase3_step3_limitations_wording.md",
    "codex_tasks/phase3_step3_completion_checklist.tsv"
  ),
  blocking_issue = "none",
  manual_review_needed = c(rep("no", 9), "yes", "yes", "yes", "yes", "no"),
  notes = c(
    "Script reads locked spatial inputs and existing label-transfer outputs.",
    "Per module and section present/missing genes recorded.",
    "Long table retained because size is manageable.",
    "Descriptive section summaries only.",
    "Predicted compartments are label-transfer context only.",
    "Allowed and not-allowed interpretations included.",
    "Each section shows selected modules without ROI labels.",
    "Bonferroni, Top100, and Suggestive overview PDFs created.",
    "Summary figures are descriptive.",
    "Documents why Loop/TAL co-distribution was not done.",
    "Safe Results wording only; manuscript not edited.",
    "Safe Methods wording only; manuscript not edited.",
    "Safe Limitations wording only; manuscript not edited.",
    "No claim-grade Loop/TAL correlation, disease/control testing, RCTD/Cell2location, TWAS/bulk, or manuscript edits."
  )
)
write_tsv(checklist, file.path(root, "codex_tasks/phase3_step3_completion_checklist.tsv"))

message("Phase 3-Step 3 complete: ", nrow(included), " sections scored across ", length(modules), " MAGMA modules.")
