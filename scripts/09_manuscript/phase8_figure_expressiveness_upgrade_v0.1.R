suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  genetic = "#3E6672",
  scrna = "#6F8F98",
  p1 = "#B59A5B",
  disease = "#9A5F52",
  papilla = "#E8D8B8",
  tubule = "#C9DDE1",
  duct = "#D8C7A5",
  interstitium = "#E9E9E9",
  muted = "#D8D8D8",
  text = "#303030",
  grid = "#A8A8A8"
)

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

label_cell <- function(x) {
  map <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

label_module <- function(x) {
  map <- c(
    magma_top50 = "MAGMA top 50",
    magma_top100 = "MAGMA top 100",
    magma_fdr05 = "MAGMA FDR",
    magma_suggestive_p1e4 = "MAGMA suggestive",
    MAGMA_top50 = "MAGMA top 50",
    MAGMA_top100 = "MAGMA top 100",
    MAGMA_FDR = "MAGMA FDR",
    MAGMA_suggestive = "MAGMA suggestive",
    P1_core_TAL_candidates = "P1 core"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

make_figure1_v03 <- function() {
  kidney <- data.table(
    x = c(0.12, 0.18, 0.22, 0.18, 0.12, 0.08),
    y = c(0.80, 0.86, 0.72, 0.55, 0.50, 0.64)
  )
  layers <- data.table(
    x = c(0.38, 0.55, 0.72, 0.89),
    y = 0.63,
    label = c("GWAS/MAGMA\nprioritization", "GSE231569\nsingle-nucleus", "P1 gene\nspectrum", "GSE73680\ndisease context"),
    fill = c(pal$genetic, pal$scrna, pal$p1, pal$disease)
  )
  fig <- ggplot() +
    geom_polygon(data = kidney, aes(x, y), fill = "#F3E7D2", color = "#7A6A55", linewidth = 0.45) +
    annotate("polygon", x = c(0.13, 0.18, 0.16), y = c(0.56, 0.64, 0.71),
             fill = pal$papilla, color = "#7A6A55", linewidth = 0.35) +
    annotate("rect", xmin = 0.06, xmax = 0.29, ymin = 0.08, ymax = 0.40,
             fill = "#FAFAFA", color = "#888888", linewidth = 0.35) +
    annotate("text", x = 0.175, y = 0.36, label = "Renal papilla context", fontface = "bold", size = 3.4, color = pal$text) +
    annotate("path", x = c(0.10, 0.13, 0.16, 0.19, 0.22), y = c(0.14, 0.25, 0.14, 0.25, 0.14),
             color = pal$genetic, linewidth = 1.0) +
    annotate("rect", xmin = 0.145, xmax = 0.205, ymin = 0.18, ymax = 0.32,
             fill = pal$duct, color = "#8A8A8A", alpha = 0.75) +
    annotate("point", x = 0.235, y = 0.13, size = 3.2, color = pal$disease) +
    annotate("text", x = 0.16, y = 0.105, label = "Loop/TAL", size = 2.8, color = pal$genetic) +
    annotate("text", x = 0.205, y = 0.305, label = "CD", size = 2.6, color = "#5B5147") +
    annotate("text", x = 0.235, y = 0.085, label = "plaque/\nstone", size = 2.35, color = pal$disease) +
    geom_segment(data = layers[1:3],
                 aes(x = x + 0.055, xend = layers$x[2:4] - 0.055, y = y, yend = y),
                 arrow = arrow(length = unit(0.10, "in")), color = "#777777", linewidth = 0.35) +
    geom_label(data = layers, aes(x = x, y = y, label = label, fill = fill),
               color = "white", fontface = "bold", size = 3.0, lineheight = 0.95,
               label.padding = unit(0.25, "lines"), linewidth = 0) +
    annotate("text", x = 0.62, y = 0.89,
             label = "Post-GWAS renal papilla cellular and disease-context mapping",
             fontface = "bold", size = 4.6, color = pal$text) +
    annotate("rect", xmin = 0.32, xmax = 0.96, ymin = 0.21, ymax = 0.39,
             fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = 0.35, y = 0.34, label = "Supported:", hjust = 0, fontface = "bold", size = 3.2, color = pal$text) +
    annotate("text", x = 0.48, y = 0.34,
             label = "Loop/TAL-associated cellular context + MAGMA module-level disease-context association",
             hjust = 0, size = 3.05, color = pal$text) +
    annotate("text", x = 0.35, y = 0.255, label = "Not established:", hjust = 0, fontface = "bold", size = 3.2, color = pal$disease) +
    annotate("text", x = 0.52, y = 0.255,
             label = "causality | TWAS convergence | colocalization | spatial validation | P1 disease-gene validation",
             hjust = 0, size = 2.9, color = pal$text) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.04, 0.95), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave("results/figures/figure1_integrative_framework_v0.3.pdf", fig, width = 13, height = 5.4, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.3.png", fig, width = 13, height = 5.4, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 1 Legend v0.3",
    "",
    "**Figure 1. Renal papilla-centered post-GWAS framework for KSD cellular and disease-context mapping.**",
    "The schematic places the evidence layers in a renal papilla context, including Loop/TAL, collecting duct and plaque/stone-papilla disease context. The analytical layers are GWAS/MAGMA prioritization, GSE231569 single-nucleus localization, P1 gene-spectrum interpretation and GSE73680 disease-context module analysis.",
    "The supported inference is a Loop/TAL-associated cellular context and MAGMA module-level disease-context association. The schematic is not spatial validation and does not establish causality, TWAS convergence, colocalization or P1 disease-gene validation."
  ), "docs/figure1_legend_v0.3.md")
}

make_figure2_v03 <- function() {
  atlas <- fread("results/tables/magma_scrna_module_score_by_celltype.tsv")
  donor <- fread("results/tables/gse231569_magma_score_by_celltype_donor.tsv")
  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")
  readiness <- fread("results/tables/gse231569_magma_score_projection_readiness_v0.1.tsv")

  atlas[, cell_label := factor(label_cell(audited_broad_cell_type), levels = label_cell(audited_broad_cell_type[order(n_cells)]))]
  atlas[, fill_class := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL", "Loop/TAL", "Other audited")]
  panel_a <- ggplot(atlas, aes(n_cells, cell_label, fill = fill_class)) +
    geom_col(width = 0.66, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = paste0("n=", n_cells)), hjust = -0.08, size = 2.35, color = pal$text) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = c("Loop/TAL" = pal$genetic, "Other audited" = pal$scrna)) +
    labs(title = "A. Audited single-nucleus atlas context", subtitle = "UMAP score projection awaits cell-level resources", x = "Cells", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, color = "#555555"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  donor <- donor[module_label %in% c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive")]
  donor[, cell_label := factor(cell_label)]
  panel_b <- ggplot(donor, aes(cell_label, module_score, fill = is_tal)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "#555555", linewidth = 0.22) +
    geom_point(aes(size = n_cells), position = position_jitter(width = 0.10, height = 0),
               shape = 21, color = "#555555", stroke = 0.18, alpha = 0.85) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("TRUE" = pal$genetic, "FALSE" = pal$muted), labels = c("Other", "Loop/TAL")) +
    scale_size_continuous(range = c(1.2, 3.8)) +
    labs(title = "B. Donor-celltype MAGMA module score summary", x = NULL, y = "Module score", fill = NULL, size = "Cells") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom", panel.grid.minor = element_blank())

  random_plot <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random_plot[, module_label := factor(label_module(gene_set),
                                       levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
  random_plot[, cell_label := factor(label_cell(audited_broad_cell_type))]
  random_plot[, support := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95,
                                   "Loop/TAL exceeds expectation", "Other")]
  panel_c <- ggplot(random_plot, aes(benchmark_percentile, cell_label, fill = support)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_point(shape = 21, size = 2.7, color = "#555555", stroke = 0.22) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("Loop/TAL exceeds expectation" = pal$genetic, Other = pal$muted)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95)) +
    labs(title = "C. Size-matched random benchmark", x = "Benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  top <- infl[order(contribution_rank)][1:12]
  top[, gene := factor(gene, levels = rev(gene))]
  panel_d <- ggplot(top, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = pal$grid, linewidth = 0.42) +
    geom_point(aes(fill = candidate_role, size = donor_detection), shape = 21, color = "#444444", stroke = 0.2) +
    scale_fill_manual(values = c(P1_core = pal$p1, candidate_TAL_driver = pal$genetic,
                                 supporting_TAL_expressed_gene = pal$scrna, non_TAL_or_low_detection = pal$muted)) +
    scale_size_continuous(range = c(2.2, 5.0)) +
    labs(title = "D. Influential genes contributing to Loop/TAL signal", x = "Contribution score", y = NULL) +
    guides(fill = "none", size = "none") +
    theme_bw(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "none", panel.grid.minor = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, rel_heights = c(0.92, 1.12))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.3.pdf", fig, width = 13.2, height = 9.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.3.png", fig, width = 13.2, height = 9.6, units = "in", dpi = 260, bg = "white")

  fwrite(readiness, "results/tables/figure2_v0.3_resource_boundary.tsv", sep = "\t")
  write_lines(c(
    "# Figure 2 Legend v0.3",
    "",
    "**Figure 2. MAGMA-prioritized KSD genes map to a Loop/TAL-associated single-nucleus expression context.**",
    "(A) Audited GSE231569 single-nucleus atlas context. A formal per-cell UMAP score projection was not generated because no usable Seurat object or cell-level embedding table was available in the current workspace.",
    "(B) Donor-celltype MAGMA module score summary across audited cell types. Loop/TAL cells retained high scores across MAGMA top-ranked, FDR and suggestive gene sets.",
    "(C) Size-matched random gene-set benchmark showing that Loop/TAL localization exceeded random expectation for MAGMA-prioritized modules.",
    "(D) Influential gene analysis ranking detected MAGMA top-50 genes by scaled MAGMA strength, Loop/TAL specificity and donor detection. These genes are expression-context contributors, not causal driver genes.",
    "Together, these analyses support a Loop/TAL-associated single-nucleus context for MAGMA-prioritized KSD genes. They do not establish causal mediation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure2_legend_v0.3.md")
}

make_figure3_v05 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  p1[, role_label := gsub("_", " ", manuscript_role)]
  p1[, gene := factor(gene, levels = gene[order(-specificity_ratio_avg)])]

  panel_a <- ggplot(p1, aes(x = 1, y = gene, fill = manuscript_role)) +
    geom_tile(color = "white", linewidth = 0.8, width = 0.72, height = 0.72) +
    geom_text(aes(label = paste0(as.character(gene), "\n", role_label)), size = 2.8, lineheight = 0.9, color = "white", fontface = "bold") +
    scale_fill_manual(values = c(
      representative_TAL_gene = pal$genetic,
      TAL_transport_candidate = "#557F89",
      calcium_ion_handling_candidate = pal$p1,
      calcium_sensing_candidate = "#9A7E43",
      supporting_context_gene = pal$scrna,
      broad_epithelial_context = pal$disease
    )) +
    labs(title = "A. P1 gene role spectrum", x = NULL, y = NULL, fill = NULL) +
    theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  panel_b <- ggplot(p1, aes(TAL_avg_expression, gene)) +
    geom_col(aes(fill = specificity_class), width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * TAL_pct_expressed)), hjust = 1.05, color = "white", size = 2.5) +
    scale_fill_manual(values = c(strong_TAL_preferential = pal$genetic, moderate_TAL_preferential = pal$p1,
                                 broad_epithelial = pal$disease, mixed_context = pal$scrna), na.value = pal$muted) +
    labs(title = "B. TAL expression and detection", x = "Average TAL expression", y = NULL, fill = "Specificity") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  panel_c <- ggplot(p1, aes(TAL_donor_detection_fraction, specificity_ratio_avg, label = as.character(gene))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_point(aes(fill = overall_evidence_class, size = -log10(magma_p)), shape = 21, color = "#555555", stroke = 0.25) +
    geom_text(nudge_y = 0.55, size = 2.6, fontface = "italic", color = pal$text) +
    scale_fill_manual(values = c(P1_strong_TAL_context = pal$genetic, P1_broad_epithelial_context = pal$disease,
                                 P1_supporting_context = pal$scrna), na.value = pal$muted) +
    labs(title = "C. Donor support versus TAL specificity", x = "TAL donor detection fraction", y = "TAL specificity ratio", fill = "Evidence", size = "-log10 MAGMA P") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  mat <- melt(p1[, .(gene, MAGMA = -log10(magma_p), TAL_specificity = pmin(specificity_ratio_avg, 5),
                     donor_detection = TAL_donor_detection_fraction, TAL_program = abs(TAL_program_rho))],
              id.vars = "gene", variable.name = "metric", value.name = "value")
  mat[, gene := factor(gene, levels = rev(levels(p1$gene)))]
  panel_d <- ggplot(mat, aes(metric, gene, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    labs(title = "D. Evidence matrix", x = NULL, y = NULL, fill = "Scaled\nsignal") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(face = "italic"), panel.grid = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2)
  ggsave("results/figures/figure3_p1_gene_evidence_v0.5.pdf", fig, width = 12.8, height = 8.5, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure3_p1_gene_evidence_v0.5.png", fig, width = 12.8, height = 8.5, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 3 Legend v0.5",
    "",
    "**Figure 3. P1 candidate genes form a TAL/transport/calcium interpretation spectrum rather than a uniform disease-gene class.**",
    "(A) Role-spectrum cards summarizing the manuscript interpretation of each P1 gene. (B) TAL average expression and percent detection in audited Loop/TAL cells. (C) Donor detection and TAL specificity, with MAGMA strength represented by point size. (D) Evidence matrix summarizing MAGMA, TAL specificity, donor detection and TAL-program association components.",
    "These panels support gene-centric expression-context interpretation. They do not establish P1 disease-gene validation or causal mechanisms."
  ), "docs/figure3_legend_v0.5.md")
}

make_figure4_v09 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  spaghetti <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

  panel_a <- ggplot() +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.72, ymax = 0.90, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.42, ymax = 0.60, fill = "#F7F8F8", color = "#888888", linewidth = 0.35) +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.12, ymax = 0.30, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("text", x = 0.50, y = 0.81, label = "GSE73680: 55 samples, 29 patients", fontface = "bold", size = 3.5, color = pal$text) +
    annotate("text", x = 0.50, y = 0.51, label = "26 paired control/adjacent -> plaque/stone papilla patients", size = 3.1, color = pal$text) +
    annotate("text", x = 0.50, y = 0.21, label = "Patient-level module delta + random benchmark", fontface = "bold", size = 3.1, color = pal$text) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.71, yend = 0.61, arrow = arrow(length = unit(0.10, "in")), color = "#777777") +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.41, yend = 0.31, arrow = arrow(length = unit(0.10, "in")), color = "#777777") +
    labs(title = "A. Patient-aware disease-context design") +
    theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"))

  p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
  p1[, signal_class := fifelse(fdr < 0.05, "q <= 0.05", fifelse(p_value < 0.05, "Nominal only", "No detectable"))]
  panel_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    scale_fill_manual(values = c("q <= 0.05" = pal$genetic, "Nominal only" = pal$p1, "No detectable" = pal$muted)) +
    labs(title = "B. Heterogeneous P1 single-gene response", x = "Paired delta", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  spaghetti <- spaghetti[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  spaghetti[, module_label := factor(as.character(module_label),
                                     levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core"))]
  panel_c <- ggplot(spaghetti, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(aes(color = direction), alpha = 0.5, linewidth = 0.33) +
    geom_point(aes(fill = group_label), shape = 21, color = "#555555", stroke = 0.18, size = 1.45) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(positive = pal$genetic, negative = pal$disease, zero = pal$muted)) +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = "#B08A45")) +
    labs(title = "C. Paired patient module shifts", x = NULL, y = "Module score", color = "Delta", fill = NULL) +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom", panel.grid.minor = element_blank())

  bench <- bench[module_name %in% c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]
  bench[, module_label := factor(label_module(module_name), levels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  bench[, benchmark_class := fifelse(percentile >= 0.95, "Exceeds 95th percentile", "Background-like")]
  panel_d <- ggplot(bench, aes(percentile, module_label, fill = benchmark_class)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    scale_fill_manual(values = c("Exceeds 95th percentile" = pal$genetic, "Background-like" = pal$muted)) +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, rel_heights = c(0.85, 1.15))
  ggsave("results/figures/figure4_gse73680_disease_context_v0.9.pdf", fig, width = 13.2, height = 9.3, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v0.9.png", fig, width = 13.2, height = 9.3, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 4 Legend v0.9",
    "",
    "**Figure 4. GSE73680 paired-patient disease-context analysis supports MAGMA-prioritized modules.**",
    "(A) Patient-aware disease-context design. (B) P1 candidate genes showed heterogeneous paired responses, without uniform FDR-significant single-gene support. (C) Paired patient-level module spaghetti plot for MAGMA top50, top100, FDR, suggestive and P1 core modules. MAGMA modules showed directionally consistent paired increases in most patients, whereas the P1 core module was weaker and heterogeneous. (D) Size-matched random gene-set benchmark for module-level disease-context association.",
    "These results support MAGMA module-level disease-context expression association in GSE73680. They do not establish genetic causality, cell-type-specific disease response, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure4_legend_v0.9.md")
}

write_manuscript_v06 <- function() {
  src <- "manuscript/manuscript_draft_v0.5.md"
  if (!file.exists(src)) return(invisible(FALSE))
  x <- readLines(src, warn = FALSE)
  note <- c(
    "",
    "## Phase 8 targeted figure-enhancement note",
    "",
    "Phase 8 added three targeted, claim-bounded enhancements: a renal papilla-centered Figure 1 schematic, a Loop/TAL influential-gene analysis for Figure 2, and paired-patient GSE73680 module spaghetti plots for Figure 4. A formal per-cell MAGMA-weighted UMAP projection was not generated because the current workspace lacks a usable GSE231569 Seurat object or cell-level embedding table. The new analyses strengthen figure expressiveness and interpretability, but do not change the manuscript boundary: the study supports TAL-associated cellular context and MAGMA module-level disease-context association, not causality, TWAS convergence, colocalization, spatial validation or P1 disease-gene validation."
  )
  write_lines(c(x, note), "manuscript/manuscript_draft_v0.6.md")
}

make_figure1_v03()
make_figure2_v03()
make_figure3_v05()
make_figure4_v09()
write_manuscript_v06()

qc <- data.table(
  figure = c("Figure 1 v0.3", "Figure 2 v0.3", "Figure 3 v0.5", "Figure 4 v0.9"),
  upgrade = c("renal papilla schematic", "donor score summary plus influential genes", "P1 role cards", "paired patient spaghetti plot"),
  new_analysis = c("no", "yes: influential genes and resource readiness", "no: visual role-spectrum redesign", "yes: paired patient module trajectories"),
  claim_boundary = c(
    "not spatial validation",
    "not per-cell scDRS; not causal driver genes",
    "not P1 disease-gene validation",
    "not causality or cell-type-specific disease response"
  ),
  status = "main_candidate"
)
fwrite(qc, "results/tables/phase8_figure_upgrade_qc_v0.1.tsv", sep = "\t")

message("wrote Phase 8 figure expressiveness upgrade outputs")
