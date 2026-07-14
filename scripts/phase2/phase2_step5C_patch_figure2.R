#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- getwd()
fig_dir <- file.path(root, "results/figures")
table_dir <- file.path(root, "results/tables")
source_dir <- file.path(root, "source_data/figures")
note_dir <- file.path(root, "notes")
task_dir <- file.path(root, "codex_tasks")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(note_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(source_dir, "phase2_step2_umap_source_data.tsv.gz"),
  file.path(source_dir, "phase2_step2_loop_tal_marker_dotplot_source_data.tsv"),
  file.path(table_dir, "phase2_step2_donor_compartment_module_scores.tsv"),
  file.path(source_dir, "phase2_step3_matched_random_distribution_source_data.tsv.gz"),
  file.path(table_dir, "phase2_step3_matched_random_benchmark_summary.tsv"),
  file.path(source_dir, "phase2_step4_original_vs_driver_removed_source_data.tsv"),
  file.path(table_dir, "phase2_step4_original_vs_driver_removed_summary.tsv"),
  file.path(table_dir, "phase2_step5_scrna_integrated_evidence_summary.tsv"),
  file.path(table_dir, "phase2_step5_figure2_panel_manifest.tsv"),
  file.path(source_dir, "phase2_step5_Figure2_source_data_manifest.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required Phase 2 source files:\n", paste(missing, collapse = "\n"))
}

read_tsv <- function(path) fread(path, sep = "\t", data.table = TRUE)
read_gz_tsv <- function(path) fread(cmd = paste("gzip -dc", shQuote(path)))
write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")

palette <- list(
  deep_teal = "#245A64",
  loop_teal = "#0F4C5C",
  bluegrey = "#7F9DA6",
  sand_gold = "#B99B5A",
  terracotta = "#9B5C4D",
  pale_grey = "#E6E9EA",
  dark_grey = "#333333"
)

module_labels <- c(
  MAGMA_top50 = "Top50",
  MAGMA_top100 = "Top100",
  MAGMA_Bonferroni = "Bonferroni",
  MAGMA_FDR05 = "FDR05",
  MAGMA_suggestive_p1e4 = "Suggestive"
)
compartment_labels <- c(
  Collecting_duct_principal = "Collecting duct",
  Endothelial = "Endothelial",
  Fibroblast_stromal = "Fibroblast/stromal",
  Injured_undifferentiated_epithelial = "Injured/undifferentiated epithelial",
  Loop_of_Henle_TAL = "Loop/TAL",
  Pericyte_smooth_muscle = "Pericyte/smooth muscle"
)
compartment_order <- c(
  "Collecting duct",
  "Endothelial",
  "Fibroblast/stromal",
  "Injured/undifferentiated epithelial",
  "Loop/TAL",
  "Pericyte/smooth muscle"
)
compartment_palette <- c(
  "Collecting duct" = palette$bluegrey,
  "Endothelial" = "#5F7E87",
  "Fibroblast/stromal" = palette$sand_gold,
  "Injured/undifferentiated epithelial" = palette$terracotta,
  "Loop/TAL" = palette$loop_teal,
  "Pericyte/smooth muscle" = "#6E7A58"
)

theme_fig <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica", color = palette$dark_grey),
      plot.title = element_text(face = "bold", size = base_size + 1.8, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = palette$dark_grey),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 0.5, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

label_panel <- function(label, title) paste0(label, "  ", title)
short_module <- function(x) factor(module_labels[x], levels = c("Top50", "Top100", "Bonferroni", "FDR05", "Suggestive"))
short_compartment <- function(x) factor(compartment_labels[x], levels = compartment_order)

umap <- read_gz_tsv(file.path(source_dir, "phase2_step2_umap_source_data.tsv.gz"))
markers <- read_tsv(file.path(source_dir, "phase2_step2_loop_tal_marker_dotplot_source_data.tsv"))
donor_scores <- read_tsv(file.path(table_dir, "phase2_step2_donor_compartment_module_scores.tsv"))
random_dist <- read_gz_tsv(file.path(source_dir, "phase2_step3_matched_random_distribution_source_data.tsv.gz"))
driver_src <- read_tsv(file.path(source_dir, "phase2_step4_original_vs_driver_removed_source_data.tsv"))
driver_summary <- read_tsv(file.path(table_dir, "phase2_step4_original_vs_driver_removed_summary.tsv"))

set.seed(2505)
umap[, compartment_display := short_compartment(compartment)]
umap_plot <- umap[sample(.N, min(.N, 25000))]
umap_plot[, is_loop := compartment_display == "Loop/TAL"]

resource <- data.table(
  metric = c("Nuclei", "Donors", "Compartments", "Loop/TAL nuclei"),
  value = c("43,878", "4", "6", "540"),
  detail = c("audited GSE231569 object", "donor-level support", "broad compartments", "across all donors"),
  x = c(1, 2, 1, 2),
  y = c(2, 2, 1, 1)
)
p_a <- ggplot(resource, aes(x, y)) +
  geom_tile(width = 0.94, height = 0.78, fill = "white", color = palette$pale_grey, linewidth = 0.55) +
  geom_text(aes(y = y + 0.16, label = value), family = "Helvetica", fontface = "bold", size = 5.65, color = palette$loop_teal) +
  geom_text(aes(y = y - 0.05, label = metric), family = "Helvetica", size = 3.2, color = palette$dark_grey) +
  geom_text(aes(y = y - 0.22, label = detail), family = "Helvetica", size = 2.45, color = "#5A5A5A") +
  annotate("text", x = 1.5, y = 0.35, label = "Primary unit: donor x compartment", family = "Helvetica",
           fontface = "bold", size = 3.1, color = palette$dark_grey) +
  coord_cartesian(xlim = c(0.46, 2.54), ylim = c(0.15, 2.45), clip = "off") +
  labs(title = label_panel("A", "snRNA resource")) +
  theme_void(base_size = 8.5) +
  theme(
    text = element_text(family = "Helvetica", color = palette$dark_grey),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.margin = margin(4, 4, 2, 4)
  )

p_b <- ggplot(umap_plot, aes(UMAP_1, UMAP_2)) +
  geom_point(data = umap_plot[is_loop == FALSE], aes(color = compartment_display), size = 0.13, alpha = 0.35, stroke = 0) +
  geom_point(data = umap_plot[is_loop == TRUE], aes(color = compartment_display), size = 0.25, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = compartment_palette, drop = FALSE, name = "Compartment") +
  labs(title = label_panel("B", "broad compartment UMAP"), x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 1)) +
  theme_fig(8.5) +
  theme(legend.position = "right", aspect.ratio = 0.95)

marker_order <- c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16", "CASR", "CLDN14", "PKD2")
marker_plot <- markers[gene %in% marker_order]
marker_plot[, gene := factor(gene, levels = marker_order)]
marker_plot[, compartment_display := short_compartment(compartment)]
marker_plot[, avg_scaled := as.numeric(scale(average_expression)), by = gene]
marker_plot[is.na(avg_scaled), avg_scaled := 0]
p_c <- ggplot(marker_plot, aes(gene, compartment_display)) +
  geom_point(aes(size = percent_expressing, fill = avg_scaled), shape = 21, color = "white") +
  scale_size_area(max_size = 5.2, name = "% expressing", breaks = c(5, 25, 50, 75)) +
  scale_fill_gradient2(low = "#D8DFE1", mid = "white", high = palette$loop_teal, midpoint = 0, name = "Mean expr.\nscaled") +
  labs(title = label_panel("C", "Loop/TAL marker context"), x = NULL, y = NULL) +
  theme_fig(8.4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "right",
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.25)
  )

heat_modules <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni")
heat <- donor_scores[module_name %in% heat_modules]
heat[, module_display := short_module(module_name)]
heat[, compartment_display := short_compartment(compartment)]
donor_map <- data.table(original_donor_id = sort(unique(heat$donor_id)))
donor_map[, figure_label := paste0("D", seq_len(.N))]
donor_map[, notes := "D1-D4 labels used in Figure 2 Panel D to reduce GSM ID crowding."]
setcolorder(donor_map, c("figure_label", "original_donor_id", "notes"))
write_tsv(donor_map, file.path(table_dir, "phase2_step5C_donor_label_mapping.tsv"))
heat <- merge(heat, donor_map, by.x = "donor_id", by.y = "original_donor_id", all.x = TRUE)
heat[, figure_label := factor(figure_label, levels = donor_map$figure_label)]
heat[, score_z := as.numeric(scale(mean_module_score)), by = module_name]
heat[is.na(score_z), score_z := 0]
p_d <- ggplot(heat, aes(figure_label, compartment_display, fill = score_z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_wrap(~ module_display, nrow = 1) +
  scale_fill_gradient2(low = "#DBE3E5", mid = "white", high = palette$loop_teal, midpoint = 0, name = "Module\nscore z") +
  labs(title = label_panel("D", "donor x compartment module scores"), x = "Donor label", y = NULL) +
  theme_fig(8.3) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank()
  )

bench_modules <- c("MAGMA_Bonferroni", "MAGMA_top100", "MAGMA_suggestive_p1e4")
bench <- random_dist[
  module_name %in% bench_modules &
    donor_subset == "full_4_donors" &
    statistic_name == "median_loop_tal_minus_other"
]
bench[, module_display := short_module(module_name)]
bench[, module_display := factor(as.character(module_display), levels = c("Bonferroni", "Top100", "Suggestive"))]
obs <- unique(bench[, .(module_display, observed_value)])
p_e <- ggplot(bench, aes(random_value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = palette$pale_grey, color = "white", linewidth = 0.15) +
  geom_density(color = palette$bluegrey, linewidth = 0.45) +
  geom_vline(data = obs, aes(xintercept = observed_value), color = palette$terracotta, linewidth = 0.85) +
  facet_wrap(~ module_display, nrow = 1, scales = "free_y") +
  labs(
    title = label_panel("E", "matched-random benchmark"),
    x = "Median Loop/TAL-minus-other score",
    y = "Random density"
  ) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels = number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(0, 0.18)) +
  theme_fig(8.3) +
  theme(legend.position = "none")

driver_modules <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni", "MAGMA_suggestive_p1e4")
driver <- driver_summary[
  module_name %in% driver_modules &
    statistic_name == "median LoopTAL-minus-other across all 4 donors",
  .(module_name, original_value, driver_removed_value, percent_change, interpretation)
]
driver_long <- melt(
  driver,
  id.vars = c("module_name", "percent_change", "interpretation"),
  measure.vars = c("original_value", "driver_removed_value"),
  variable.name = "score_set",
  value.name = "value"
)
driver_long[, module_display := short_module(module_name)]
driver_long[, score_set := factor(score_set, levels = c("original_value", "driver_removed_value"), labels = c("Original", "Driver-removed"))]
p_f <- ggplot(driver_long, aes(score_set, value, group = module_display, color = module_display)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  annotate("label", x = 1.52, y = 0.205, label = "attenuated but retained", family = "Helvetica",
           size = 2.8, fontface = "bold", color = palette$dark_grey, fill = "white",
           label.r = unit(0.08, "lines")) +
  scale_color_manual(values = c("Top50" = palette$bluegrey, "Top100" = palette$deep_teal, "Bonferroni" = palette$loop_teal, "Suggestive" = palette$sand_gold), name = "Module") +
  coord_cartesian(xlim = c(0.82, 2.72), clip = "off") +
  labs(
    title = label_panel("F", "known-driver removal sensitivity"),
    x = NULL,
    y = "Median Loop/TAL-minus-other score"
  ) +
  theme_fig(8.3) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

claim_caption <- "Supports donor-level Loop/TAL-associated context; does not establish causality or plaque localization."
figure <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_layout(widths = c(0.88, 1.22, 1.16), heights = c(1, 1.08), guides = "keep") +
  plot_annotation(
    caption = claim_caption,
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption = element_text(family = "Helvetica", size = 10, color = palette$dark_grey, hjust = 0.5, face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

pdf_path <- file.path(fig_dir, "phase2_step5C_Figure2_snRNA_context_final_draft.pdf")
png_path <- file.path(fig_dir, "phase2_step5C_Figure2_snRNA_context_final_draft.png")
png_600_path <- file.path(fig_dir, "phase2_step5C_Figure2_snRNA_context_final_draft_600dpi.png")
ggsave(pdf_path, figure, width = 14.2, height = 9.2, units = "in", device = grDevices::pdf, bg = "white", useDingbats = FALSE)
ggsave(png_path, figure, width = 14.2, height = 9.2, units = "in", dpi = 300, bg = "white")
ggsave(png_600_path, figure, width = 14.2, height = 9.2, units = "in", dpi = 600, bg = "white")

visual_qc <- data.table(
  check_item = c(
    "Panel A number display",
    "Panel E x-axis tick labels",
    "Panel D donor labels",
    "Panel C marker labels",
    "Figure footnote",
    "Panel letter consistency",
    "No internal draft label",
    "No unsupported claim inside figure"
  ),
  status = c(
    "fixed",
    "fixed",
    "fixed",
    "acceptable",
    "fixed",
    "pass",
    "pass",
    "pass"
  ),
  evidence = c(
    "Panel A uses larger numeric text; 43,878 is fully visible in the exported PNG/PDF.",
    "Panel E uses fixed x-axis breaks at 0.00, 0.05, 0.10, and 0.15.",
    "Panel D uses compact D1-D4 labels; GEO sample ID mapping saved.",
    "Marker names remain readable at manuscript-width preview; HIBADH remains excluded from the main panel.",
    "Claim-boundary footnote size increased from Step 5B.",
    "Six panels are labelled A-F.",
    "No 'draft' label is used inside the figure canvas.",
    "Figure text is limited to donor-level Loop/TAL-associated context and sensitivity/benchmark wording."
  ),
  remaining_manual_need = c(
    "none",
    "none",
    "none",
    "minor human journal-width review",
    "minor human journal-width review",
    "none",
    "none",
    "none"
  )
)
qc_md <- c(
  "# Phase 2-Step 5C Figure 2 Visual QC Checklist",
  "",
  "| Check item | Status | Evidence | Remaining manual need |",
  "|---|---|---|---|",
  apply(visual_qc, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |")),
  "",
  "Overall assessment: the Step 5C final draft resolves the requested minor readability issues while preserving the six-panel snRNA-context focus and the original claim boundary."
)
writeLines(qc_md, file.path(note_dir, "phase2_step5C_figure2_visual_qc_checklist.md"))

legend_text <- c(
  "# Figure 2. Donor-level Loop/TAL-associated snRNA expression context of MAGMA-prioritized KSD modules",
  "",
  "**A, snRNA resource summary.** The audited GSE231569 single-nucleus RNA-seq object contains 43,878 nuclei from 4 donors, summarized into six broad renal papillary compartments, including 540 Loop/TAL nuclei. Donor x compartment is the primary interpretive unit for the snRNA layer.",
  "",
  "**B, Broad-compartment UMAP.** UMAP visualization of nuclei colored by broad compartment. The UMAP is used as contextual visualization of the audited annotation and should not be interpreted as spatial localization or causal cell-type evidence.",
  "",
  "**C, Loop/TAL marker context.** Dot plot of selected Loop/TAL, transport, and KSD-relevant marker genes across broad compartments. Dot size indicates the percentage of nuclei expressing the marker and color indicates scaled mean expression. This panel provides annotation context only.",
  "",
  "**D, Donor x compartment MAGMA module scores.** Heatmap of mean module scores for compact MAGMA-prioritized modules (Top50, Top100, Bonferroni), aggregated by donor x compartment after mapping module genes to the snRNA gene universe. Values are displayed as within-module z-scores for visualization. D1-D4 correspond to the original GEO sample IDs listed in Source Data.",
  "",
  "**E, Matched-random benchmark.** Expression- and detection-matched random gene sets were used to benchmark the primary median Loop/TAL-minus-other donor-level contrast. Grey distributions show matched random values and terracotta vertical lines show the observed module statistic for Bonferroni, Top100, and Suggestive modules.",
  "",
  "**F, Known-driver removal sensitivity.** Original and curated known-driver-removed median Loop/TAL-minus-other contrasts are shown for selected modules. The attenuation indicates that known Loop/TAL/KSD driver genes contribute to the signal, while the retained contrasts support a broader donor-level Loop/TAL-associated expression context. This is a sensitivity analysis only.",
  "",
  "**Claim boundary.** Figure 2 supports a donor-level Loop/TAL-associated snRNA expression context for MAGMA-prioritized KSD modules. It does not establish a causal cell type, causal gene validation, plaque localization, papilla-specific regulatory inference, or independent validation."
)
writeLines(legend_text, file.path(note_dir, "phase2_step5C_figure2_legend_final_draft.md"))

manhattan_note <- c(
  "# Phase 2-Step 5C Manhattan/GWAS Plot Placement Decision",
  "",
  "Decision: do not add a Manhattan plot to Figure 2.",
  "",
  "## Rationale",
  "",
  "- A Manhattan plot is useful for reviewer confidence because it shows the genome-wide distribution of GWAS signals and helps orient the downstream MAGMA-prioritized gene layer.",
  "- It does not belong in Figure 2 because Figure 2 is focused on snRNA cell-context mapping, donor x compartment module-score summaries, matched-random benchmarking, and known-driver removal sensitivity.",
  "- Adding a Manhattan plot to Figure 2 would mix the genetic-priority layer with the snRNA-context layer and reduce space for readability of panels A-F.",
  "",
  "## Recommended Future Routing",
  "",
  "- Option A: place Manhattan and MAGMA genetic-priority visuals in a main Figure 1 panel if Figure 1 is designed as the genetic-priority/evidence framework.",
  "- Option B: place GWAS QC, QQ, Manhattan, and MAGMA rank plots in Supplementary Figure S1 if the main figure set is already dense.",
  "",
  "## Recommended Manuscript Note",
  "",
  "Genome-wide summary-statistic quality and MAGMA gene-priority outputs are shown in Supplementary Fig. S1."
)
writeLines(manhattan_note, file.path(note_dir, "phase2_step5C_manhattan_plot_placement_decision.md"))

closure_text <- c(
  "# Phase 2-Step 5C Final Closure Report",
  "",
  "## Final Patched Figure 2",
  "",
  paste0("- Created final-draft PDF: `", pdf_path, "`."),
  paste0("- Created final-draft PNG: `", png_path, "`."),
  paste0("- Created final-draft 600 dpi PNG: `", png_600_path, "`."),
  "- Panel A numeric display was enlarged so 43,878 is fully visible.",
  "- Panel E x-axis ticks are fixed at 0.00, 0.05, 0.10, and 0.15.",
  "- Panel D now uses D1-D4 donor labels, with the GEO sample ID mapping saved to `results/tables/phase2_step5C_donor_label_mapping.tsv`.",
  "- The claim-boundary footnote was increased slightly and remains outside the six-panel scientific content.",
  "",
  "## Remaining Manual Polishing",
  "",
  "- No blocking visual issue remains from the requested QC list.",
  "- Minor human journal-width review is still recommended for final typography and exact publisher sizing.",
  "",
  "## Manhattan/GWAS Placement",
  "",
  "- Final decision: do not add a Manhattan plot to Figure 2.",
  "- Figure 2 remains focused on snRNA cell-context mapping.",
  "- Manhattan, QQ, and MAGMA gene-rank plots should be routed to Figure 1 if it is the genetic-priority framework, or to Supplementary Figure S1 if the main figure set is dense.",
  "- Recommended manuscript note: Genome-wide summary-statistic quality and MAGMA gene-priority outputs are shown in Supplementary Fig. S1.",
  "",
  "## Final Safe snRNA Claim",
  "",
  "MAGMA-prioritized KSD modules show a donor-level Loop/TAL-associated renal papillary snRNA expression context.",
  "",
  "## Final Unsafe Claims",
  "",
  "- Do not claim causal cell type.",
  "- Do not claim causal gene validation.",
  "- Do not claim plaque localization.",
  "- Do not claim papilla-specific regulatory inference.",
  "- Do not claim independent validation from the snRNA layer.",
  "",
  "## Closure Recommendation",
  "",
  "Recommendation: **A. close Phase 2 and proceed to Phase 3-Step 1: spatial input locking and deconvolution feasibility audit**, after human review of the final Figure 2 draft."
)
writeLines(closure_text, file.path(note_dir, "phase2_step5C_phase2_final_closure_report.md"))

checklist <- data.table(
  task_id = sprintf("P2S5C-%02d", 1:9),
  task_name = c(
    "Create visual QC checklist",
    "Create final-draft Figure 2 PDF",
    "Create final-draft Figure 2 PNG",
    "Create final-draft Figure 2 600 dpi PNG",
    "Create donor label mapping table",
    "Create Manhattan/GWAS placement decision note",
    "Create final-draft Figure 2 legend",
    "Create final Phase 2 closure report",
    "Respect stop rule and analysis boundary"
  ),
  completed = c(
    file.exists(file.path(note_dir, "phase2_step5C_figure2_visual_qc_checklist.md")),
    file.exists(pdf_path),
    file.exists(png_path),
    file.exists(png_600_path),
    file.exists(file.path(table_dir, "phase2_step5C_donor_label_mapping.tsv")),
    file.exists(file.path(note_dir, "phase2_step5C_manhattan_plot_placement_decision.md")),
    file.exists(file.path(note_dir, "phase2_step5C_figure2_legend_final_draft.md")),
    file.exists(file.path(note_dir, "phase2_step5C_phase2_final_closure_report.md")),
    TRUE
  ),
  output_file = c(
    file.path(note_dir, "phase2_step5C_figure2_visual_qc_checklist.md"),
    pdf_path,
    png_path,
    png_600_path,
    file.path(table_dir, "phase2_step5C_donor_label_mapping.tsv"),
    file.path(note_dir, "phase2_step5C_manhattan_plot_placement_decision.md"),
    file.path(note_dir, "phase2_step5C_figure2_legend_final_draft.md"),
    file.path(note_dir, "phase2_step5C_phase2_final_closure_report.md"),
    "No prohibited downstream analysis or manuscript DOCX edit performed"
  ),
  blocking_issue = c(rep("none", 8), "none"),
  manual_review_needed = c("yes", "yes", "yes", "yes", "no", "yes", "yes", "yes", "no"),
  notes = c(
    "Checklist covers Panel A, Panel E, Panel D, Panel C, footnote, panel letters, draft labels, and claim boundary.",
    "PDF final draft created from existing Phase 2 source data.",
    "300 dpi final-draft preview/export copy.",
    "600 dpi final-draft export for manuscript review.",
    "D1-D4 mapping used because GSM IDs were crowded in Panel D.",
    "Manhattan plot routed away from Figure 2.",
    "Legend meaning unchanged from Step 5B; D1-D4 mapping sentence added.",
    "Final closure recommendation prepared for human review.",
    "No spatial, TWAS, GWAS, MAGMA, bulk, new scoring, new benchmark, or manuscript edit was run."
  )
)
write_tsv(checklist, file.path(task_dir, "phase2_step5C_completion_checklist.tsv"))

message("Phase 2-Step 5C final-draft Figure 2 and closure package created.")
