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
  muted = "#D8D8D8",
  pale = "#EEF3F4",
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
    Perivascular_mural_like = "Perivascular/mural-like",
    Pericyte_smooth_muscle = "Perivascular/mural-like"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

make_figure1_v04 <- function() {
  fig <- ggplot() +
    annotate("polygon", x = c(0.07, 0.13, 0.19, 0.23, 0.17, 0.10),
             y = c(0.67, 0.86, 0.90, 0.74, 0.56, 0.52),
             fill = "#F2E6D2", color = "#7A6A55", linewidth = 0.5) +
    annotate("polygon", x = c(0.13, 0.16, 0.19, 0.15),
             y = c(0.57, 0.72, 0.66, 0.59), fill = "#E4D2AA", color = "#7A6A55", linewidth = 0.35) +
    annotate("text", x = 0.15, y = 0.48, label = "kidney-to-papilla\nzoom", size = 2.8, color = pal$text) +
    annotate("rect", xmin = 0.055, xmax = 0.27, ymin = 0.10, ymax = 0.38, fill = "#FAFAFA", color = "#888888", linewidth = 0.35) +
    annotate("rect", xmin = 0.07, xmax = 0.255, ymin = 0.12, ymax = 0.32, fill = "#F1F1F1", color = NA) +
    annotate("path", x = c(0.09, 0.12, 0.15, 0.18, 0.21), y = c(0.15, 0.26, 0.15, 0.26, 0.15),
             color = pal$genetic, linewidth = 1.0) +
    annotate("rect", xmin = 0.155, xmax = 0.205, ymin = 0.17, ymax = 0.31,
             fill = "#D8C7A5", color = "#8A8A8A", alpha = 0.78) +
    annotate("point", x = 0.232, y = 0.145, size = 3.1, color = pal$disease) +
    annotate("text", x = 0.145, y = 0.105, label = "Loop/TAL", size = 2.55, color = pal$genetic) +
    annotate("text", x = 0.21, y = 0.30, label = "CD", size = 2.45, color = "#5B5147") +
    annotate("text", x = 0.235, y = 0.095, label = "plaque/\nstone", size = 2.15, color = pal$disease) +
    annotate("text", x = 0.62, y = 0.90,
             label = "Renal papilla-centered post-GWAS evidence chain",
             fontface = "bold", size = 4.6, color = pal$text)
  boxes <- data.table(
    x = c(0.36, 0.53, 0.70, 0.88),
    y = 0.65,
    label = c("GWAS/MAGMA\n57 loci\n17,316 genes\n94 Bonferroni",
              "GSE231569\nsnRNA-seq context\naudited cells",
              "P1 candidates\n6 genes\nrole spectrum",
              "GSE73680\n55 samples\n29 patients\n26 paired"),
    fill = c(pal$genetic, pal$scrna, pal$p1, pal$disease)
  )
  fig <- fig +
    geom_segment(data = boxes[1:3], aes(x = x + 0.065, xend = boxes$x[2:4] - 0.065, y = y, yend = y),
                 arrow = arrow(length = unit(0.10, "in")), color = "#777777", linewidth = 0.35) +
    geom_label(data = boxes, aes(x = x, y = y, label = label, fill = fill),
               color = "white", fontface = "bold", size = 2.65, lineheight = 0.92,
               label.padding = unit(0.24, "lines"), linewidth = 0) +
    annotate("rect", xmin = 0.34, xmax = 0.96, ymin = 0.22, ymax = 0.40, fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = 0.365, y = 0.35, label = "Supported:", hjust = 0, fontface = "bold", size = 3.15, color = pal$text) +
    annotate("text", x = 0.485, y = 0.35, label = "Loop/TAL cellular context + MAGMA module-level disease-context association", hjust = 0, size = 2.95, color = pal$text) +
    annotate("text", x = 0.365, y = 0.27, label = "Not established:", hjust = 0, fontface = "bold", size = 3.15, color = pal$disease) +
    annotate("text", x = 0.52, y = 0.27, label = "causality | TWAS convergence | colocalization | spatial validation | P1 disease-gene validation", hjust = 0, size = 2.65, color = pal$text) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.06, 0.95), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave("results/figures/figure1_integrative_framework_v0.4.pdf", fig, width = 13, height = 5.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.4.png", fig, width = 13, height = 5.6, units = "in", dpi = 260, bg = "white")
  write_lines(c("# Figure 1 Legend v0.4", "", "**Figure 1. Renal papilla-centered post-GWAS evidence chain for KSD cellular and disease-context mapping.**",
                "The schematic links primary public KSD GWAS/MAGMA prioritization to audited GSE231569 renal papillary snRNA-seq context, P1 candidate role-spectrum interpretation and GSE73680 plaque/stone papilla disease-context analysis. The supported inference is Loop/TAL-associated cellular context and MAGMA module-level disease-context association. The papilla schematic is conceptual and does not constitute spatial validation."), "docs/figure1_legend_v0.4.md")
}

make_figure2_v04 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := label_cell(audited_broad_cell_type)]
  set.seed(2)
  plot_cells <- if (nrow(top50) > 22000) top50[sample(.N, 22000)] else top50
  plot_cells[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  p_a <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_label, alpha = is_tal), size = 0.17) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.36), guide = "none") +
    labs(title = "A. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right", panel.grid = element_blank())
  qlim <- quantile(plot_cells$celllevel_module_score, c(0.02, 0.98), na.rm = TRUE)
  p_b <- ggplot(plot_cells, aes(UMAP_1, UMAP_2, color = celllevel_module_score)) +
    geom_point(size = 0.17, alpha = 0.76) +
    scale_color_gradient(low = "#E8ECEE", high = pal$genetic, limits = qlim, oob = scales::squish) +
    labs(title = "B. MAGMA top 50 score on UMAP", x = "UMAP 1", y = "UMAP 2", color = "Score") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())

  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  random <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
                                  labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  random[, cell_label := factor(label_cell(audited_broad_cell_type))]
  random[, support := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95, "Loop/TAL exceeded expectation", "Other")]
  p_c <- ggplot(random, aes(benchmark_percentile, cell_label, fill = support)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    annotate("text", x = 0.94, y = Inf, label = "95th percentile", hjust = 1, vjust = 1.2, size = 2.1, color = "#555555") +
    geom_point(shape = 21, size = 2.6, color = "#555555", stroke = 0.22) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("Loop/TAL exceeded expectation" = pal$genetic, Other = pal$muted)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95)) +
    labs(title = "C. Size-matched benchmark percentile", x = "Benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = pal$grid, linewidth = 0.42) +
    geom_point(aes(fill = group), shape = 21, color = "#444444", stroke = 0.22, size = 3.6) +
    scale_fill_manual(values = c("P1 candidate" = pal$p1, "Other MAGMA gene" = pal$genetic)) +
    labs(title = "D. Genes contributing to Loop/TAL-associated signal", x = "Contribution score", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.0) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.05, 1.05))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.4.pdf", fig, width = 13.2, height = 9.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.4.png", fig, width = 13.2, height = 9.6, units = "in", dpi = 260, bg = "white")
  write_lines(c("# Figure 2 Legend v0.4", "", "**Figure 2. MAGMA-prioritized KSD genes project to a Loop/TAL-associated renal papillary single-nucleus context.**",
                "(A) Audited GSE231569 snRNA-seq UMAP colored by harmonized cell type, with Loop/TAL retained as the primary audited epithelial transport context. (B) Per-cell MAGMA top 50 module score projected onto the same UMAP. (C) Size-matched random gene-set benchmark showing that Loop/TAL localization exceeded the 95th percentile for MAGMA-prioritized modules. (D) Genes contributing to the Loop/TAL-associated signal, ranked by MAGMA strength, Loop/TAL specificity and donor detection. These genes are not interpreted as causal drivers."), "docs/figure2_legend_v0.4.md")
}

make_figure3_v06 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  p1[, role_label := gsub("_", " ", manuscript_role)]
  p1[, gene := factor(gene, levels = c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"))]
  role_colors <- c(representative_TAL_gene = pal$genetic, TAL_transport_candidate = "#557F89",
                   calcium_ion_handling_candidate = pal$p1, calcium_sensing_candidate = "#9A7E43",
                   supporting_context_gene = pal$scrna, broad_epithelial_context = pal$disease)
  p_a <- ggplot(p1, aes(1, gene, fill = manuscript_role)) +
    geom_tile(width = 0.82, height = 0.72, color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(as.character(gene), "\n", role_label)), color = "white", fontface = "bold", size = 2.7, lineheight = 0.9) +
    scale_fill_manual(values = role_colors) +
    labs(title = "A. P1 gene role map", x = NULL, y = NULL) +
    theme_void(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(label_cell(cell_type), levels = label_cell(keep_cells))]
  cell[, gene := factor(gene, levels = levels(p1$gene))]
  p_b <- ggplot(cell, aes(cell_label, gene)) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    scale_size_continuous(range = c(0.8, 5.0), labels = scales::percent_format()) +
    labs(title = "B. P1 expression across audited cell types", x = NULL, y = NULL, fill = "Average\nexpression", size = "Detected") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          axis.text.y = element_text(face = "italic"), legend.position = "bottom", panel.grid = element_blank())

  p_c <- ggplot(p1, aes(specificity_ratio_avg, gene)) +
    geom_segment(aes(x = 0, xend = specificity_ratio_avg, yend = gene), color = pal$grid, linewidth = 0.5) +
    geom_point(aes(size = TAL_donor_detection_fraction, fill = specificity_class), shape = 21, color = "#555555", stroke = 0.25) +
    scale_fill_manual(values = c(strong_TAL_preferential = pal$genetic, moderate_TAL_preferential = pal$p1), na.value = pal$muted) +
    scale_size_continuous(range = c(2.2, 5.2), labels = scales::percent_format()) +
    labs(title = "C. TAL specificity and donor support", x = "TAL specificity ratio", y = NULL, fill = "Specificity", size = "Donors") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  mat <- melt(p1[, .(gene, MAGMA = -log10(magma_p), TAL_specificity = pmin(specificity_ratio_avg, 10),
                     donor_detection = TAL_donor_detection_fraction * 10,
                     GSE73680_single_gene = fifelse(gene == "PKD2", 3, 1))],
              id.vars = "gene", variable.name = "metric", value.name = "value")
  mat[, gene := factor(gene, levels = rev(levels(p1$gene)))]
  p_d <- ggplot(mat, aes(metric, gene, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", value)), size = 2.5, color = "#303030") +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    labs(title = "D. Evidence matrix", x = NULL, y = NULL, fill = "Scaled\nsignal") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          axis.text.y = element_text(face = "italic"), panel.grid = element_blank())
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2)
  ggsave("results/figures/figure3_p1_gene_evidence_v0.6.pdf", fig, width = 13.0, height = 8.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure3_p1_gene_evidence_v0.6.png", fig, width = 13.0, height = 8.6, units = "in", dpi = 260, bg = "white")
  write_lines(c("# Figure 3 Legend v0.6", "", "**Figure 3. P1 candidate genes form a TAL, transport and calcium-handling interpretation spectrum.**",
                "(A) Role map for the six P1 genes. (B) Dot plot of P1 expression and detection across audited GSE231569 cell types. (C) TAL specificity and donor support summarized as a compact lollipop plot. (D) Evidence matrix summarizing MAGMA support, TAL specificity, donor detection and GSE73680 single-gene context. The figure supports gene-centric interpretation, not P1 disease-gene validation."), "docs/figure3_legend_v0.6.md")
}

make_figure4_v10 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")
  p_a <- ggplot() +
    annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.65, ymax = 0.88, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.50, y = 0.765, label = "GSE73680: 55 samples | 29 patients | 26 paired", fontface = "bold", size = 3.35, color = pal$text) +
    annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.35, ymax = 0.55, fill = "#F7F8F8", color = "#888888", linewidth = 0.3) +
    annotate("text", x = 0.50, y = 0.45, label = "control/adjacent -> plaque/stone papilla", size = 3.0, color = pal$text) +
    annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.09, ymax = 0.25, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.50, y = 0.17, label = "paired module delta + random benchmark", fontface = "bold", size = 3.0, color = pal$text) +
    labs(title = "A. Compact patient-aware design") + theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"))
  p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No detectable")]
  p1[, p_label := sprintf("P=%.3f", p_value)]
  p_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.35, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$p1, "No detectable" = pal$muted)) +
    labs(title = "B. Heterogeneous P1 single-gene response", x = "Paired delta", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.0) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())
  keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
  long <- long[module_label %in% keep]
  long[, module_label := factor(as.character(module_label), levels = keep)]
  lab <- sumdt[, .(module_label = factor(as.character(module_label), levels = keep),
                   label = sprintf("%d/%d positive", n_positive_delta, n_paired_patients))]
  p_c <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(aes(color = direction), alpha = 0.45, linewidth = 0.32) +
    geom_point(aes(fill = group_label), shape = 21, color = "#555555", stroke = 0.18, size = 1.35) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 0.9, color = "#303030") +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.25, size = 2.45, color = pal$text) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(positive = "#7F9AA3", negative = "#D6A29A", zero = pal$muted)) +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = "#B08A45")) +
    labs(title = "C. Paired patient module shifts", x = NULL, y = "Module score", color = "Delta", fill = NULL) +
    theme_bw(base_size = 8.4) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom", panel.grid.minor = element_blank())
  bench <- bench[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  bench[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                                 labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  bench[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  p_d <- ggplot(bench, aes(percentile, module_label, fill = percentile >= 0.95)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = emp_label), hjust = 1.05, size = 2.45, color = pal$text) +
    scale_fill_manual(values = c("TRUE" = pal$genetic, "FALSE" = pal$muted), labels = c("Background-like", "Exceeds 95th percentile")) +
    coord_cartesian(xlim = c(0, 1.05)) +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.0) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.82, 1.18))
  ggsave("results/figures/figure4_gse73680_disease_context_v1.0.pdf", fig, width = 13.2, height = 9.3, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v1.0.png", fig, width = 13.2, height = 9.3, units = "in", dpi = 260, bg = "white")
  write_lines(c("# Figure 4 Legend v1.0", "", "**Figure 4. GSE73680 supports MAGMA module-level plaque/stone papilla disease-context association.**",
                "(A) Compact patient-aware design. (B) P1 single-gene paired responses were heterogeneous, with no FDR-supported uniform P1 response. (C) Paired patient spaghetti plots show MAGMA module shifts across control/adjacent and plaque/stone papilla samples; labels show positive paired shifts. (D) Size-matched random benchmark with empirical P values. The analysis supports module-level disease-context association, not causal validation or cell-type-specific disease response."), "docs/figure4_legend_v1.0.md")
}

write_manuscript_v08 <- function() {
  x <- readLines("manuscript/manuscript_draft_v0.7.md", warn = FALSE)
  methods_insert <- grep("^### Resource-limited TWAS", x)[1] - 1
  methods_text <- c(
    "",
    "### Cell-level MAGMA score projection and functional-context analyses",
    "",
    "After Phase 9 resource re-audit, an audited GSE231569 Seurat object with UMAP coordinates was available locally. MAGMA top 50, top 100, FDR and suggestive gene sets were scored at the cell level using mean z-scored expression of detected module genes, and donor-cell-type summaries were derived from per-cell scores. GSE73680 injury/remodeling, inflammation/immune, fibrosis/ECM and epithelial-injury programs were scored as mean z-scored marker-set modules. Coupling between MAGMA modules and injury programs was evaluated using Spearman sample-level correlations, paired patient-level delta correlations and patient/group residual correlations. Functional interpretation used local GO Biological Process enrichment through clusterProfiler and org.Hs.eg.db, with curated nephron and functional marker-set enrichment used as an auditable local supplement."
  )
  if (!any(grepl("^### Cell-level MAGMA score projection", x))) x <- append(x, methods_text, after = methods_insert)

  res_insert <- grep("^## Discussion", x)[1] - 1
  res_text <- c(
    "",
    "### Functional-context analyses link prioritized modules to nephron development, transport and papillary injury programs",
    "",
    "Cell-level GSE231569 score projection supported the Figure 2 localization pattern on the audited UMAP, with MAGMA top-ranked module scores visibly enriched in the Loop/TAL-associated region. Functional enrichment added a pathway-level interpretation layer. GO Biological Process analysis of MAGMA top 100 genes identified loop of Henle development, distal tubule development, response to vitamin D, response to metal ion and phosphate ion homeostasis among the leading terms. Curated nephron and functional marker-set analyses further linked prioritized and Loop/TAL-influential genes to TAL transport, calcium-ion handling and epithelial junction contexts. In GSE73680, MAGMA-prioritized modules were strongly coupled with papillary injury/remodeling programs at the module level, including MAGMA top 50 correlation with injury/remodeling scores. These analyses support functional interpretation and disease-context coupling, but do not establish a causal injury mechanism.",
    ""
  )
  if (!any(grepl("^### Functional-context analyses link", x))) x <- append(x, res_text, after = res_insert)

  legend_insert <- length(x)
  fig5_legend <- c(
    "",
    "### Figure 5. Functional interpretation of MAGMA-prioritized TAL-associated KSD genes",
    "",
    "Figure 5 summarizes integrated evidence tiers, GO Biological Process enrichment, curated nephron/functional enrichment and GSE73680 risk-injury module coupling. The figure supports functional interpretation and module-level disease-context coupling of prioritized KSD genes, but does not establish causal mechanism, TWAS convergence, colocalization or spatial validation."
  )
  if (!any(grepl("^### Figure 5", x))) x <- append(x, fig5_legend, after = legend_insert)
  write_lines(x, "manuscript/manuscript_draft_v0.8.md")
}

make_figure1_v04()
make_figure2_v04()
make_figure3_v06()
make_figure4_v10()
write_manuscript_v08()

qc <- data.table(
  figure = c("Figure 1 v0.4", "Figure 2 v0.4", "Figure 3 v0.6", "Figure 4 v1.0", "Figure 5 v0.1"),
  upgrade = c("papilla schematic with key numbers", "true UMAP plus cell-level MAGMA score", "role map plus dotplot and compact evidence", "compact design plus p/empirical labels and median paired lines", "functional enrichment and injury coupling"),
  status = "main_candidate",
  boundary = c("not spatial validation", "not causal cell-type mediation", "not P1 validation", "not causal validation", "functional interpretation only")
)
fwrite(qc, "results/tables/phase9_main_figure_qc_v0.1.tsv", sep = "\t")
message("wrote Phase 9 polished figures and manuscript v0.8")
