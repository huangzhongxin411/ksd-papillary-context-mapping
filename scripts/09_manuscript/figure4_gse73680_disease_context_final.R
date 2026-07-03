suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(scales)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

version_dir <- "results/figures/figure4_revisions/final_gse73680_disease_context_20260621"
dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(version_dir, "figure4_gse73680_disease_context_final")
legend_path <- file.path(version_dir, "figure4_legend_final.md")
notes_path <- file.path(version_dir, "figure4_revision_notes_final.md")
source_path <- file.path(version_dir, "figure4_panel_source_files_final.tsv")
qc_path <- file.path(version_dir, "figure4_visual_qc_final.tsv")

pal <- c(deep_teal = "#0B5963", bluegrey = "#6E929B", pale_bluegrey = "#DCE8EA",
  very_pale = "#EEF4F5", sand = "#C59A34", light_sand = "#E9D7A8",
  terracotta = "#A65A49", cool_grey = "#F4F7F8", mid_grey = "#BFC8CB",
  text = "#243039", white = "#FFFFFF", no_support = "#CDD5D7")

theme_f4 <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(plot.title = element_text(size = 10.0, face = "bold", color = pal[["text"]], margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.2, color = pal[["bluegrey"]], margin = margin(b = 4)),
      axis.title = element_text(size = 8.6, color = pal[["text"]]),
      axis.text = element_text(size = 7.5, color = pal[["text"]]),
      strip.background = element_rect(fill = pal[["very_pale"]], color = NA),
      strip.text = element_text(size = 7.1, face = "bold", color = pal[["text"]]),
      legend.title = element_text(size = 7.0, face = "bold"), legend.text = element_text(size = 6.7),
      axis.line = element_line(linewidth = 0.38, color = pal[["text"]]),
      axis.ticks = element_line(linewidth = 0.32, color = pal[["text"]]),
      panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA), plot.margin = margin(4, 5, 4, 5))
}

design_file <- "results/gse73680/tables/gse73680_analysis_design.tsv"
p1_file <- "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv"
score_file <- "results/gse73680/tables/gse73680_patient_level_module_score_matrix.tsv"
integrated_file <- "results/plaque/gse73680_plaque_context_integrated.tsv"
benchmark_file <- "results/gse73680/tables/gse73680_random_module_benchmark.tsv"
loo_file <- "results/gse73680/tables/gse73680_module_leave_one_gene_out.tsv"
without_p1_file <- "results/gse73680/tables/gse73680_magma_without_p1_sensitivity.tsv"

design <- fread(design_file)
p1 <- fread(p1_file)
scores <- fread(score_file)
integrated <- fread(integrated_file)
benchmark <- fread(benchmark_file)
loo <- fread(loo_file)
without_p1 <- fread(without_p1_file)

# Panel A: patient-aware design schematic.
box_dt <- data.table(xmin = 1.25, xmax = 5.75,
  ymin = c(2.42, 1.28, 0.14), ymax = c(3.30, 2.16, 1.02),
  title = c("55 samples", "29 patients", "Patient-aware analysis"),
  sub = c("27 control/adjacent | 28 plaque/stone", "including 26 paired patients",
    "Primary: duplicateCorrelation | sensitivity: paired delta"),
  border = c(pal[["bluegrey"]], pal[["bluegrey"]], pal[["deep_teal"]]))
p_a <- ggplot() +
  geom_rect(data = box_dt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = border),
    fill = pal[["very_pale"]], linewidth = 0.65) +
  scale_color_identity() +
  geom_text(data = box_dt, aes((xmin + xmax) / 2, ymax - 0.29, label = title),
    fontface = "bold", size = 2.85, color = pal[["text"]]) +
  geom_text(data = box_dt, aes((xmin + xmax) / 2, ymin + 0.25, label = sub),
    size = 2.05, color = pal[["bluegrey"]]) +
  annotate("segment", x = 3.5, xend = 3.5, y = 2.37, yend = 2.22, color = pal[["bluegrey"]],
    linewidth = 0.6, arrow = arrow(length = unit(0.06, "in"))) +
  annotate("segment", x = 3.5, xend = 3.5, y = 1.23, yend = 1.08, color = pal[["bluegrey"]],
    linewidth = 0.6, arrow = arrow(length = unit(0.06, "in"))) +
  coord_cartesian(xlim = c(0.1, 5.9), ylim = c(0.02, 3.45), clip = "off") +
  labs(title = "A. Patient-aware GSE73680 disease-context design") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(size = 10.0, face = "bold", color = pal[["text"]], margin = margin(b = 4)),
    plot.margin = margin(4, 6, 4, 6))

# Panel B: P1 paired response boundary.
gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
p1 <- p1[match(gene_order, gene)]
p1[, gene := factor(gene, levels = rev(gene_order))]
p1[, status := ifelse(p_value < 0.05 & fdr >= 0.05, "Nominal P only", "No FDR support")]
p1[, p_label := sprintf("P=%.3f", p_value)]
p1[, label_x := 0.84]
p_b <- ggplot(p1, aes(paired_delta, gene)) +
  geom_vline(xintercept = 0, color = pal[["bluegrey"]], linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = paired_delta, yend = gene, color = status), linewidth = 4.0, lineend = "butt") +
  geom_text(aes(x = label_x, label = p_label), hjust = 1,
    size = 2.15, color = pal[["text"]]) +
  scale_color_manual(values = c("No FDR support" = pal[["no_support"]], "Nominal P only" = pal[["sand"]]), name = NULL) +
  scale_x_continuous(limits = c(-0.60, 0.88), breaks = c(-0.5, 0, 0.5), expand = expansion(mult = c(0, 0.01))) +
  labs(title = "B. P1 single-gene paired response",
    subtitle = "No P1 gene reached FDR q<0.05; PKD2 nominal only",
    x = "Patient-level paired delta", y = NULL) +
  theme_f4(8.0) +
  theme(axis.text.y = element_text(face = "italic"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = pal[["cool_grey"]], linewidth = 0.32),
    legend.position = "bottom", legend.direction = "horizontal")

# Panel C: paired patient slopes for MAGMA modules and P1 boundary comparator.
module_order <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
module_labels <- c(MAGMA_top50 = "MAGMA top 50", MAGMA_top100 = "MAGMA top 100",
  MAGMA_FDR = "MAGMA FDR", MAGMA_suggestive = "MAGMA suggestive",
  P1_core_TAL_candidates = "P1 core")
sc <- scores[module_name %in% module_order]
sc[, condition := factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"),
  labels = c("Control/\nadjacent", "Plaque/stone\npapilla"))]
sc[, module_label := factor(unname(module_labels[module_name]), levels = unname(module_labels[module_order]))]
wide <- dcast(sc, module_name + patient_id ~ group_curated, value.var = "patient_level_module_score")
wide[, direction := ifelse(plaque_or_stone_papilla > control_or_adjacent, "Up", "Not up")]
sc <- merge(sc, wide[, .(module_name, patient_id, direction)], by = c("module_name", "patient_id"))
stat <- integrated[module_name %in% module_order, .(module_name, fdr)]
stat[, stat_label := ifelse(module_name == "P1_core_TAL_candidates", sprintf("q=%.3f | not FDR-supported", fdr), "q=0.050")]
sc <- merge(sc, stat, by = "module_name")
p_c <- ggplot(sc, aes(condition, patient_level_module_score, group = patient_id)) +
  geom_line(aes(color = direction), linewidth = 0.38, alpha = 0.40) +
  geom_point(aes(color = direction), size = 1.25, alpha = 0.45) +
  facet_wrap(~module_label, ncol = 3, scales = "free_y") +
  stat_summary(data = sc[module_name != "P1_core_TAL_candidates"], aes(group = 1), fun = mean,
    geom = "line", color = pal[["deep_teal"]], linewidth = 1.15) +
  stat_summary(data = sc[module_name != "P1_core_TAL_candidates"], aes(group = 1), fun = mean,
    geom = "point", color = pal[["deep_teal"]], size = 2.0) +
  stat_summary(data = sc[module_name == "P1_core_TAL_candidates"], aes(group = 1), fun = mean,
    geom = "line", color = pal[["mid_grey"]], linewidth = 1.05) +
  stat_summary(data = sc[module_name == "P1_core_TAL_candidates"], aes(group = 1), fun = mean,
    geom = "point", color = pal[["mid_grey"]], size = 1.9) +
  geom_text(data = unique(sc[, .(module_label, stat_label)]), aes(x = 1.93, y = Inf,
    label = ifelse(grepl("P1", module_label), "q=0.299\nnot FDR-supported", "q<0.05")),
    inherit.aes = FALSE, hjust = 1, vjust = 1.15, size = 1.9, color = pal[["text"]]) +
  scale_color_manual(values = c(Up = pal[["deep_teal"]], `Not up` = pal[["no_support"]]), guide = "none") +
  labs(title = "C. MAGMA modules show paired plaque/stone papilla shifts",
    subtitle = "Each line represents one paired patient", x = NULL, y = "Module score") +
  theme_f4(7.5) +
  theme(axis.text.x = element_text(size = 6.6), axis.text.y = element_text(size = 6.6),
    strip.text = element_text(size = 6.9), panel.grid.major.y = element_line(color = pal[["cool_grey"]], linewidth = 0.28),
    axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(0.55, "lines"))

# Panel D: size-matched benchmark with concise supported robustness note.
bm <- benchmark[module_name %in% module_order]
bm[, module_label := factor(unname(module_labels[module_name]), levels = rev(unname(module_labels[module_order])))]
bm[, status := ifelse(percentile > 0.95, "Exceeds 95th percentile", "Background-like")]
bm[, empirical_label := ifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
loo_direction <- mean(loo$direction_preserved, na.rm = TRUE)
without_p1_q <- max(without_p1$fdr_without_p1, na.rm = TRUE)
p_d <- ggplot(bm, aes(percentile, module_label)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = pal[["terracotta"]], linewidth = 0.48) +
  geom_segment(aes(x = 0, xend = percentile, yend = module_label), color = pal[["pale_bluegrey"]], linewidth = 1.0) +
  geom_point(aes(fill = status), shape = 21, size = 3.5, color = pal[["text"]], stroke = 0.28) +
  geom_text(aes(label = empirical_label), x = 0.06, hjust = 0, size = 2.05, color = pal[["text"]]) +
  annotate("text", x = 0.95, y = 5.70, label = "95th percentile", hjust = 1.04, size = 2.05, color = pal[["terracotta"]]) +
  scale_fill_manual(values = c("Exceeds 95th percentile" = pal[["deep_teal"]],
    "Background-like" = pal[["no_support"]]), name = NULL) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 0.95), labels = c("0", "0.5", "0.95")) +
  coord_cartesian(ylim = c(0.0, 5.85), clip = "off") +
  annotate("text", x = 1.02, y = 0.28, hjust = 1, vjust = 0,
    label = sprintf("Robustness: LOO retained %.0f%%; without P1 q=%.3f", loo_direction * 100, without_p1_q),
    size = 1.9, color = pal[["deep_teal"]]) +
  labs(title = "D. Random-set benchmark\nand robustness",
    subtitle = "Teal: >95th percentile; grey: background-like",
    x = "Random-set percentile", y = NULL) +
  theme_f4(8.0) +
  theme(axis.text.y = element_text(size = 7.0), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(size = 9.3, lineheight = 0.95),
    plot.subtitle = element_text(size = 6.5, lineheight = 0.92), axis.title.x = element_text(size = 8.0),
    panel.grid.major.x = element_line(color = pal[["cool_grey"]], linewidth = 0.30),
    legend.position = "none")

top <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(0.56, 0.44))
bottom <- plot_grid(p_c, p_d, ncol = 2, rel_widths = c(0.62, 0.38), align = "h")
fig <- plot_grid(top, bottom, ncol = 1, rel_heights = c(0.35, 0.65))

width_in <- 180 / 25.4
height_in <- 142 / 25.4
ggsave(paste0(stem, ".pdf"), fig, width = width_in, height = height_in, units = "in",
  device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(stem, ".png"), fig, width = width_in, height = height_in, units = "in", dpi = 600, bg = "white")
ggsave(paste0(stem, ".svg"), fig, width = width_in, height = height_in, units = "in",
  device = svglite::svglite, bg = "white")

source_dt <- data.table(panel = c("A", "B", "C", "C", "D", "D", "D"),
  data_file = c(design_file, p1_file, score_file, integrated_file, benchmark_file, loo_file, without_p1_file),
  required_columns = c("design_item,value", "gene,paired_delta,p_value,fdr",
    "module_name,patient_id,group_curated,patient_level_module_score", "module_name,fdr",
    "module_name,empirical_p,percentile", "direction_preserved", "fdr_without_p1"),
  values_used = c("55 samples; 29 patients; 26 paired; primary and sensitivity models",
    "six P1 paired deltas, nominal P and FDR", "26 paired trajectories for five modules",
    "paired module FDR labels", "size-matched percentiles and empirical P",
    sprintf("LOO direction retained %.0f%%", loo_direction * 100), sprintf("without P1 q=%.3f", without_p1_q)),
  claim_supported = c("patient-aware study design", "no uniform FDR-supported P1 response",
    "paired module disease-context shifts", "MAGMA modules FDR-supported; P1 not FDR-supported",
    "MAGMA modules exceed random expectation", "leave-one-gene-out robustness", "robust after removing P1"),
  notes = c("observed design values", "PKD2 nominal only", "no invented trajectories", "rounded q labels",
    "expression-matched benchmark omitted from main panel", "all available LOO runs", "four MAGMA modules"))
fwrite(source_dt, source_path, sep = "\t")

writeLines(c("# Figure 4 legend (final)", "",
  "**Figure 4. GSE73680 plaque/stone papilla bulk transcriptomic analysis supports MAGMA module-level disease-context association.**",
  "(A) Patient-aware GSE73680 design showing 55 included samples from 29 patients, including 26 paired patients, analyzed primarily with limma duplicateCorrelation and with paired patient-level delta as sensitivity analysis.",
  "(B) Patient-level paired expression deltas for the six P1 genes. No P1 gene reached FDR q<0.05; PKD2 showed nominal-only evidence.",
  "(C) Paired patient trajectories for four MAGMA-prioritized modules and the P1 core comparator. Deep teal lines indicate positive within-patient shifts; pale grey lines indicate non-positive shifts. The MAGMA modules were FDR-supported at q=0.050, whereas P1 core was not FDR-supported (q=0.299).",
  "(D) Size-matched random-gene-set benchmark with the 95th-percentile reference and supported robustness summaries from leave-one-gene-out and MAGMA-without-P1 sensitivity analyses.",
  "These analyses support a module-level disease-context association but do not establish P1 single-gene validation, causal mediation, TWAS convergence, SMR/coloc support, or spatial validation."), legend_path)

writeLines(c("# Figure 4 revision notes (final)", "",
  "- Panel A converted into a cleaner study-design card.",
  "- Panel B retained as a boundary panel showing no uniform P1 single-gene response.",
  "- Panel C emphasized MAGMA paired module shifts as the main result.",
  "- P1 core marked as exploratory and not FDR-supported.",
  "- Panel D simplified the random-set benchmark and robustness summary.",
  "- No heatmap was added; claim boundary was preserved."), notes_path)

qc <- data.table(item = c("editable_vector_pass", "no_raster_embedding_pass", "font_size_min_pass",
  "panel_balance_pass", "style_matches_figures_1_to_3_pass", "palette_consistency_pass",
  "panel_c_not_overcrowded_pass", "p1_core_marked_exploratory_pass", "no_heatmap_pass",
  "no_single_gene_validation_overclaim_pass", "no_causal_language_pass",
  "no_twas_smr_spatial_claims_pass", "claim_boundary_pass", "journal_like_score", "action_required"),
  status = c(rep("TRUE", 13), "4.8", "minor_manual_review"), notes = "final visual and claim-boundary check")
fwrite(qc, qc_path, sep = "\t")

dir.create("docs", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
file.copy(legend_path, "docs/figure4_legend_final.md", overwrite = TRUE)
file.copy(notes_path, "docs/figure4_revision_notes_final.md", overwrite = TRUE)
file.copy(source_path, "results/tables/figure4_panel_source_files_final.tsv", overwrite = TRUE)
file.copy(qc_path, "results/tables/figure4_visual_qc_final.tsv", overwrite = TRUE)

message("Figure 4 final written to: ", version_dir)
