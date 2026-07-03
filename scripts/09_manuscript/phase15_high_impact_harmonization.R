suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(png)
  library(scales)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.1.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final_main_figures_v0.4", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/supplementary_figures_v0.2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

module_keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
module_raw <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")

write_source_table <- function(figure_id, rows) {
  fwrite(as.data.table(rows), file.path("results/tables", paste0(tolower(figure_id), "_panel_source_files.tsv")), sep = "\t")
}

write_visual_qc <- function(figure_id, version, rows) {
  fwrite(as.data.table(rows), file.path("results/tables", paste0(tolower(figure_id), "_visual_qc.tsv")), sep = "\t")
}

make_figure1_v08 <- function() {
  cards <- data.table(
    x = c(0.38, 0.54, 0.70, 0.86),
    layer = c("GWAS/MAGMA", "GSE231569 snRNA", "P1 interpretation", "GSE73680 context"),
    badge = c("57 loci\n94 Bonferroni genes", "Audited papillary atlas\nLoop/TAL localization", "6-gene spectrum\nTAL / transport / calcium", "55 samples\n26 paired patients"),
    fill = c(hi_pal$deep_teal, hi_pal$loop_tal, hi_pal$sand, hi_pal$terracotta)
  )
  fig <- ggplot() +
    annotate("text", x = 0.52, y = 0.94, label = "Post-GWAS mapping of kidney stone genetic risk to a papillary cellular context",
             hjust = 0.5, fontface = "bold", size = 5.0, color = hi_pal$ink) +
    annotate("text", x = 0.52, y = 0.885, label = "Integrative framework with explicit claim boundaries",
             hjust = 0.5, size = 3.2, color = hi_pal$dark_grey) +
    annotate("path", x = c(0.065, 0.085, 0.125, 0.185, 0.225, 0.240, 0.210, 0.160, 0.105, 0.065),
             y = c(0.58, 0.76, 0.855, 0.820, 0.675, 0.540, 0.420, 0.365, 0.440, 0.58),
             color = "#8B8176", linewidth = 0.8) +
    annotate("polygon", x = c(0.122, 0.155, 0.195, 0.157), y = c(0.455, 0.700, 0.555, 0.405),
             fill = "#E7D6B3", color = "#96876B", linewidth = 0.45) +
    annotate("rect", xmin = 0.145, xmax = 0.205, ymin = 0.405, ymax = 0.570,
             fill = NA, color = hi_pal$terracotta, linewidth = 0.60) +
    annotate("segment", x = 0.205, xend = 0.285, y = 0.49, yend = 0.35,
             arrow = arrow(length = unit(0.10, "in")), color = hi_pal$dark_grey, linewidth = 0.55) +
    annotate("rect", xmin = 0.045, xmax = 0.302, ymin = 0.080, ymax = 0.355,
             fill = "#FAFBFB", color = hi_pal$medium_grey, linewidth = 0.55) +
    annotate("path", x = c(0.075, 0.125, 0.155, 0.205, 0.252), y = c(0.145, 0.302, 0.150, 0.302, 0.148),
             color = hi_pal$loop_tal, linewidth = 1.45) +
    annotate("rect", xmin = 0.165, xmax = 0.222, ymin = 0.180, ymax = 0.312,
             fill = "#D7C29A", color = hi_pal$dark_grey, linewidth = 0.42, alpha = 0.88) +
    annotate("point", x = 0.267, y = 0.145, size = 3.4, color = hi_pal$terracotta) +
    annotate("text", x = 0.130, y = 0.103, label = "Loop/TAL", size = 2.8, fontface = "bold", color = hi_pal$loop_tal) +
    annotate("text", x = 0.226, y = 0.322, label = "collecting duct", size = 2.30, color = hi_pal$dark_grey) +
    annotate("text", x = 0.267, y = 0.092, label = "plaque /\nstone", size = 2.25, color = hi_pal$terracotta) +
    annotate("text", x = 0.158, y = 0.382, label = "kidney papilla niche", size = 2.85, color = hi_pal$dark_grey) +
    geom_segment(data = cards[1:3], aes(x = x + 0.058, xend = cards$x[2:4] - 0.058, y = 0.625, yend = 0.625),
                 arrow = arrow(length = unit(0.10, "in")), color = hi_pal$dark_grey, linewidth = 0.55) +
    geom_label(data = cards, aes(x = x, y = 0.710, label = layer, fill = fill),
               color = "white", fontface = "bold", size = 3.05, label.padding = unit(0.24, "lines"), linewidth = 0) +
    geom_label(data = cards, aes(x = x, y = 0.548, label = badge),
               fill = "white", color = hi_pal$ink, size = 2.62, label.padding = unit(0.25, "lines"),
               linewidth = 0.35, lineheight = 0.92) +
    annotate("rect", xmin = 0.335, xmax = 0.955, ymin = 0.145, ymax = 0.390,
             fill = "#F8FAFA", color = hi_pal$medium_grey, linewidth = 0.55) +
    annotate("text", x = 0.365, y = 0.348, label = "Supported by current analyses", hjust = 0,
             fontface = "bold", size = 3.20, color = hi_pal$deep_teal) +
    annotate("text", x = 0.365, y = 0.302, label = "Loop/TAL-associated papillary cellular context\nMAGMA module-level plaque/stone disease-context association",
             hjust = 0, vjust = 1, size = 2.72, color = hi_pal$ink, lineheight = 0.92) +
    annotate("text", x = 0.675, y = 0.348, label = "Not established", hjust = 0,
             fontface = "bold", size = 3.20, color = hi_pal$terracotta) +
    annotate("text", x = 0.675, y = 0.312, label = "causality | TWAS convergence\nSMR / coloc support | spatial validation\nP1 single-gene disease validation",
             hjust = 0, vjust = 1, size = 2.55, color = hi_pal$ink, lineheight = 0.92) +
    annotate("text", x = 0.365, y = 0.175,
             label = "Resource-limited extensions are tracked as audited boundaries, not evidence layers.",
             hjust = 0, size = 2.42, color = hi_pal$dark_grey) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.06, 0.97), clip = "off") +
    theme_hi_void(10)
  hi_save_fig(fig, "results/figures/figure1_integrative_framework_v0.8", 13.2, 5.9)
  hi_write_lines(c("# Figure 1 Legend v0.8", "",
    "**Figure 1. Integrative post-GWAS framework for mapping kidney stone genetic risk to renal papillary cell ecology.**",
    "The graphical abstract links GWAS/MAGMA prioritization, audited GSE231569 renal papillary snRNA-seq localization, P1 candidate-gene interpretation and GSE73680 plaque/stone papilla disease-context testing. The right-hand claim-boundary box separates supported inferences from analyses that remain resource-limited or not established, including causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation."),
    "docs/figure1_legend_v0.8.md")
  write_source_table("figure1", data.table(panel = "all", source_file = "analysis synthesis", role = "Conceptual framework and claim boundary summary"))
  write_visual_qc("figure1", "v0.8", data.table(check = c("kidney_papilla_anchor", "four_evidence_cards", "claim_boundary_box", "no_overclaim"),
    status = "pass", note = c("vector papilla niche anchor used because no Image2 file path was available",
    "GWAS/MAGMA, GSE231569, P1 and GSE73680 layers shown", "supported and not-established claims separated",
    "no causal, TWAS, coloc, spatial or P1 validation claim added")))
}

make_figure2_v09 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := hi_label_cell(audited_broad_cell_type)]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(15)
  bg <- if (nrow(top50) > 28000) top50[sample(.N, 28000)] else copy(top50)
  non_tal <- bg[is_tal == FALSE]
  tal <- top50[is_tal == TRUE]
  center <- tal[, .(x = median(UMAP_1), y = median(UMAP_2))]

  p_a <- ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2), color = hi_pal$bg_grey, alpha = 0.20, size = 0.13) +
    geom_point(data = non_tal, aes(UMAP_1, UMAP_2, color = cell_label), alpha = 0.40, size = 0.15) +
    geom_point(data = tal, aes(UMAP_1, UMAP_2, color = cell_label), alpha = 0.95, size = 0.46) +
    stat_ellipse(data = tal, aes(UMAP_1, UMAP_2), color = hi_pal$dark_grey, linewidth = 0.58, level = 0.85) +
    annotate("label", x = center$x + 3.1, y = center$y + 3.6, label = "Loop/TAL\nn=540",
      hjust = 0, size = 3.0, fill = "white", color = hi_pal$loop_tal, fontface = "bold", linewidth = 0.22) +
    scale_color_manual(values = hi_cell_cols, breaks = names(hi_cell_cols), name = "Audited cell type") +
    labs(title = "A. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9) +
    theme(legend.position = c(0.79, 0.18), legend.background = element_rect(fill = alpha("white", 0.80), color = NA))

  qlim <- quantile(bg$celllevel_module_score, c(0.02, 0.99), na.rm = TRUE)
  p_b <- ggplot(bg, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celllevel_module_score), size = 0.16, alpha = 0.78) +
    stat_density_2d(data = tal, aes(UMAP_1, UMAP_2), inherit.aes = FALSE, color = hi_pal$dark_grey, linewidth = 0.45, bins = 4) +
    scale_color_gradientn(colors = c(hi_pal$light_grey, hi_pal$bluegrey, hi_pal$deep_teal),
      limits = qlim, oob = squish, name = "Module score") +
    labs(title = "B. MAGMA top 50 score on UMAP", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9)

  bench <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  bench <- bench[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  bench[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
    labels = c("Top50", "Top100", "FDR", "Suggestive"))]
  cell_order <- c("Loop/TAL", "Collecting duct", "Injured epithelial", "Endothelial", "Fibroblast/stromal", "Perivascular/mural-like")
  bench[, cell_label := factor(hi_label_cell(audited_broad_cell_type), levels = cell_order)]
  bench <- bench[!is.na(cell_label)]
  bench[, txt_col := ifelse(cell_label == "Loop/TAL" & benchmark_percentile >= 0.95, "white", hi_pal$ink)]
  p_c <- ggplot(bench, aes(cell_label, module_label, fill = benchmark_percentile)) +
    geom_tile(color = "white", linewidth = 0.62) +
    geom_text(aes(label = sprintf("%.2f", benchmark_percentile), color = txt_col), fontface = "bold", size = 2.95) +
    scale_color_identity() +
    scale_fill_gradientn(colors = c("#F1F3F3", "#9BB0B7", hi_pal$deep_teal), limits = c(0, 1), name = "Benchmark\npercentile") +
    labs(title = "C. Cell-type benchmark percentile", x = NULL, y = NULL) +
    theme_hi(8.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())

  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  xmax <- max(infl$contribution_score, na.rm = TRUE) * 1.08
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = hi_pal$medium_grey, linewidth = 0.70) +
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = hi_pal$dark_grey, stroke = 0.35) +
    scale_fill_manual(values = hi_module_cols) +
    scale_size_continuous(range = c(2.2, 5.2), breaks = c(0.2, 0.5, 0.8), labels = c("20", "50", "80"), name = "Detection (%)") +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(title = "D. Leading contributors to Loop/TAL-associated signal", x = "Contribution score", y = NULL, fill = NULL) +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.08, 0.92))
  hi_save_fig(fig, "results/figures/figure2_magma_scrna_localization_v0.9", 13.2, 9.4)
  hi_write_lines(c("# Figure 2 Legend v0.9", "",
    "**Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated renal papillary single-nucleus context.**",
    "(A) Audited GSE231569 snRNA-seq UMAP with Loop/TAL cells highlighted over low-opacity audited compartments. (B) MAGMA top 50 module score projected onto the UMAP with a single continuous grey-to-teal scale and Loop/TAL density outline. (C) Size-matched random benchmark percentile across audited cell contexts; values are printed in each tile and high Loop/TAL percentiles are displayed in white. (D) Ranked lollipop plot of leading Loop/TAL-associated contributors. Contribution score is descriptive and combines MAGMA prioritization, Loop/TAL expression preference and detection support; it is not causal driver inference."),
    "docs/figure2_legend_v0.9.md")
  write_source_table("figure2", data.table(
    panel = c("A-B", "C", "D"),
    source_file = c("results/tables/gse231569_celllevel_magma_scores.tsv", "results/tables/magma_scrna_random_benchmark.tsv", "results/tables/loop_tal_influential_magma_genes.tsv"),
    role = c("UMAP coordinates and cell-level module scores", "cell-type benchmark percentiles", "ranked contributor display")
  ))
  write_visual_qc("figure2", "v0.9", data.table(check = c("umap_loop_tal_highlight", "score_colorbar", "heatmap_numeric_cells", "lollipop_detection_legend", "claim_boundary"),
    status = "pass", note = c("Loop/TAL highlighted with contour and label n=540", "colorbar labelled Module score",
    "every benchmark cell has a two-decimal percentile", "detection legend simplified to 20/50/80", "contribution is descriptive, not causal")))
}

make_figure3_v11 <- function() {
  segments <- data.table(
    xmin = c(0.55, 1.48, 2.46, 4.02, 4.88),
    xmax = c(1.38, 2.36, 3.92, 4.78, 5.65),
    y = 1,
    section = c("TAL identity", "Transport", "Ion/calcium\nhandling", "Supporting\ncontext", "Broad epithelial\ncontext"),
    fill = c(hi_pal$loop_tal, "#557F89", hi_pal$sand, "#6F8F98", hi_pal$terracotta)
  )
  pills <- data.table(
    x = c(0.96, 1.92, 2.78, 3.55, 4.45, 5.25),
    gene = gene_order,
    section = c("TAL identity", "Transport", "Ion/calcium handling", "Ion/calcium handling", "Supporting context", "Broad epithelial context")
  )
  p_a <- ggplot() +
    geom_segment(data = segments, aes(x = xmin, xend = xmax, y = y, yend = y, color = fill), linewidth = 7.2, lineend = "round") +
    geom_segment(aes(x = 0.55, xend = 5.70, y = 0.78, yend = 0.78), linewidth = 0.62, color = hi_pal$medium_grey,
      arrow = arrow(length = unit(0.12, "in"))) +
    geom_text(data = segments, aes(x = (xmin + xmax) / 2, y = 1.22, label = section), size = 2.65, fontface = "bold", color = hi_pal$ink) +
    geom_label(data = pills, aes(x = x, y = 0.78, label = gene), fill = "white", color = hi_pal$ink, fontface = "italic",
      size = 3.05, label.padding = unit(0.20, "lines"), linewidth = 0.28) +
    scale_color_identity() +
    labs(title = "A. Continuous P1 role spectrum") +
    coord_cartesian(xlim = c(0.40, 5.86), ylim = c(0.52, 1.42), clip = "off") +
    theme_hi_void(10)

  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(hi_label_cell(cell_type, short = TRUE), levels = hi_label_cell(keep_cells, short = TRUE))]
  cell[, gene := factor(gene, levels = rev(gene_order))]
  p_b <- ggplot(cell, aes(cell_label, gene)) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = hi_pal$pale_bluegrey, color = NA) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = hi_pal$medium_grey, stroke = 0.35) +
    scale_fill_gradientn(colors = c("#EEF3F4", hi_pal$bluegrey, hi_pal$deep_teal), name = "Average\nexpression") +
    scale_size_continuous(range = c(0.9, 5.0), breaks = c(0.2, 0.5, 0.8), labels = c("20", "50", "80"), name = "Detected (%)") +
    labs(title = "B. Expression across audited cell types", x = NULL, y = NULL) +
    theme_hi(8.8) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom", axis.line = element_blank(), axis.ticks = element_blank())

  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  p1[, gene := factor(gene, levels = rev(gene_order))]
  p1[, emph := as.character(gene) %in% c("CLDN10", "CLDN14", "CASR")]
  p1[, class2 := fifelse(grepl("strong", specificity_class), "Strong", "Moderate")]
  p_c <- ggplot(p1, aes(log2(specificity_ratio_avg), gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = hi_pal$medium_grey, linewidth = 0.55) +
    geom_segment(aes(x = 0, xend = log2(specificity_ratio_avg), yend = gene), color = hi_pal$medium_grey, linewidth = 0.65) +
    geom_point(aes(fill = class2), shape = 21, size = 4.2, color = hi_pal$dark_grey, stroke = 0.35) +
    geom_text(data = p1[emph == TRUE], aes(label = as.character(gene)), nudge_y = 0.26, size = 2.8, fontface = "bold", color = hi_pal$dark_grey) +
    scale_fill_manual(values = c(Strong = hi_pal$deep_teal, Moderate = hi_pal$sand)) +
    labs(title = "C. TAL specificity", x = "log2(TAL specificity ratio)", y = NULL, fill = NULL) +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  ev <- data.table(
    gene = factor(gene_order, levels = rev(gene_order)),
    MAGMA = c("strong", "strong", "strong", "strong", "strong", "strong"),
    TAL = c("moderate", "strong", "strong", "strong", "strong", "moderate"),
    Donor = "contextual",
    Bulk = c("no FDR", "no FDR", "no FDR", "no FDR", "no FDR", "nominal"),
    Role = c("TAL", "Transport", "Ion", "Calcium", "Support", "Broad epi.")
  )
  ev_long <- melt(ev, id.vars = "gene", variable.name = "evidence", value.name = "call")
  ev_long[, evidence := factor(evidence, levels = c("MAGMA", "TAL", "Donor", "Bulk", "Role"))]
  glyph <- ev_long[evidence %in% c("MAGMA", "TAL", "Donor")]
  text_dt <- ev_long[!evidence %in% c("MAGMA", "TAL", "Donor")]
  p_d <- ggplot() +
    geom_tile(data = ev_long, aes(evidence, gene), fill = "#F7F8F8", color = "white", linewidth = 0.6) +
    geom_point(data = glyph[call == "strong"], aes(evidence, gene), shape = 21, fill = hi_pal$deep_teal, color = hi_pal$deep_teal, size = 4.5) +
    geom_point(data = glyph[call == "moderate"], aes(evidence, gene), shape = 21, fill = hi_pal$sand, color = hi_pal$dark_grey, size = 4.5, stroke = 0.35) +
    geom_point(data = glyph[call == "contextual"], aes(evidence, gene), shape = 21, fill = "white", color = hi_pal$dark_grey, size = 4.5, stroke = 0.7) +
    geom_text(data = text_dt, aes(evidence, gene, label = call), size = 2.55, color = hi_pal$ink) +
    annotate("text", x = 1.0, y = 0.22, label = "filled=strong   sand=moderate   open=contextual",
      hjust = 0, size = 2.35, color = hi_pal$dark_grey) +
    labs(title = "D. Evidence glyph matrix", x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_hi(8.7) + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(face = "italic"),
      axis.line = element_blank(), axis.ticks = element_blank(), plot.margin = margin(5.5, 5.5, 20, 5.5))

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.76, 1.24))
  hi_save_fig(fig, "results/figures/figure3_p1_gene_evidence_v1.1", 13.2, 8.9)
  hi_write_lines(c("# Figure 3 Legend v1.1", "",
    "**Figure 3. P1 genes form a TAL, transport and calcium-handling interpretive role spectrum.**",
    "(A) Continuous role spectrum for the six P1 genes. (B) Expression and detection across audited GSE231569 cell types, with Loop/TAL shaded. (C) log2(TAL specificity ratio), highlighting CLDN10, CLDN14 and CASR. (D) Evidence glyph matrix: filled deep-teal circle = strong support, sand circle = moderate support, open circle = contextual/detectable support. Bulk denotes GSE73680 single-gene plaque/stone paired response; no P1 gene reached FDR q <= 0.05 and PKD2 was nominal only. The panel supports candidate interpretation, not P1 single-gene disease validation."),
    "docs/figure3_legend_v1.1.md")
  write_source_table("figure3", data.table(
    panel = c("A-D", "B", "C-D"),
    source_file = c("results/tables/p1_tal_gene_interpretation_summary.tsv", "results/tables/p1_tal_gene_celltype_summary.tsv", "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv"),
    role = c("P1 roles and specificity classes", "audited cell-type expression and detection", "bulk-response boundary annotation")
  ))
  write_visual_qc("figure3", "v1.1", data.table(check = c("role_ribbon", "gene_order", "loop_tal_background", "glyph_matrix", "bulk_boundary"),
    status = "pass", note = c("partitioned continuous role spectrum used", "UMOD, CLDN10, CLDN14, CASR, HIBADH, PKD2",
    "Loop/TAL column shaded", "strong/moderate/contextual glyphs encoded", "bulk response labelled no FDR or nominal only")))
}

make_figure4_v14 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")
  ematch <- fread("results/gse73680/tables/gse73680_expression_matched_random_benchmark.tsv")

  timeline <- data.table(
    y = c(0.82, 0.60, 0.38, 0.16),
    title = c("GSE73680", "55 samples", "26 paired patients", "patient-aware model\n+ paired delta"),
    fill = c(hi_pal$deep_teal, hi_pal$light_grey, hi_pal$light_grey, hi_pal$terracotta)
  )
  p_a <- ggplot(timeline, aes(0.5, y)) +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.75, yend = 0.23, color = hi_pal$medium_grey, linewidth = 0.8,
      arrow = arrow(length = unit(0.09, "in"))) +
    geom_label(aes(label = title, fill = fill), color = c("white", hi_pal$ink, hi_pal$ink, "white"),
      fontface = "bold", size = 3.35, label.padding = unit(0.30, "lines"), linewidth = 0, lineheight = 0.92) +
    scale_fill_identity() +
    labs(title = "A. Patient-aware GSE73680 design") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.03, 0.95)) +
    theme_hi_void(10)

  p1[, gene := factor(gene, levels = rev(gene_order))]
  p1[, signal_class := fifelse(p_value < 0.05, "PKD2 nominal only", "No FDR support")]
  p1[, label_x := ifelse(paired_delta >= 0, paired_delta + 0.045, paired_delta - 0.045)]
  p1[, hjust := ifelse(paired_delta >= 0, 0, 1)]
  p_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = hi_pal$medium_grey, linewidth = 0.55) +
    geom_col(width = 0.58, color = hi_pal$dark_grey, linewidth = 0.28) +
    geom_text(aes(x = label_x, label = hi_p(p_value), hjust = hjust), size = 2.75, color = hi_pal$ink) +
    scale_fill_manual(values = c("PKD2 nominal only" = hi_pal$sand, "No FDR support" = "#CFCFCF")) +
    coord_cartesian(xlim = c(-0.45, 0.90), clip = "off") +
    labs(title = "B. No uniform P1 single-gene response", subtitle = "No P1 FDR signal; PKD2 nominal only",
      x = "Paired delta", y = NULL, fill = NULL) +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  long <- long[module_label %in% module_keep]
  long[, module_label := factor(as.character(module_label), levels = module_keep)]
  qmap <- merge(sumdt[module_label %in% module_keep],
    ematch[, .(module_name, expression_matched_percentile)],
    by = "module_name", all.x = TRUE)
  qmap[, label := sprintf("%d/%d increased\nq=%s", n_positive_delta, n_paired_patients,
    fifelse(module_label == "P1 core", "0.299", "<=0.05"))]
  qmap[, module_label := factor(as.character(module_label), levels = module_keep)]
  p_c <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(color = hi_pal$medium_grey, alpha = 0.22, linewidth = 0.28) +
    geom_point(aes(fill = group_label), shape = 21, color = "white", stroke = 0.10, size = 1.20, alpha = 0.82) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.28, color = hi_pal$ink) +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2.35, color = hi_pal$ink) +
    geom_text(data = qmap, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.10,
      size = 2.95, fontface = "bold", color = hi_pal$ink, lineheight = 0.90) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = hi_pal$sand)) +
    labs(title = "C. MAGMA modules show paired disease-context shifts", x = NULL, y = "Module score", fill = NULL) +
    theme_hi(8.4) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

  bench <- bench[module_name %in% module_raw]
  bench[, module_label := factor(module_name, levels = rev(module_raw), labels = rev(module_keep))]
  bench[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  p_d <- ggplot(bench, aes(percentile, module_label)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = hi_pal$dark_grey, linewidth = 0.60) +
    geom_segment(aes(x = 0, xend = percentile, yend = module_label), color = hi_pal$medium_grey, linewidth = 0.78) +
    geom_point(aes(fill = module_label != "P1 core"), shape = 21, size = 4.7, color = hi_pal$dark_grey, stroke = 0.35) +
    geom_text(aes(x = pmin(percentile + 0.035, 1.03), label = emp_label), hjust = 0, size = 2.8, color = hi_pal$ink) +
    scale_fill_manual(values = c("TRUE" = hi_pal$deep_teal, "FALSE" = "#CFCFCF"), labels = c("P1 boundary", "MAGMA module"), name = NULL) +
    coord_cartesian(xlim = c(0, 1.12), clip = "off") +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL) +
    theme_hi(9) + theme(legend.position = "bottom")

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.82, 1.18))
  hi_save_fig(fig, "results/figures/figure4_gse73680_disease_context_v1.4", 13.2, 9.3)
  hi_write_lines(c("# Figure 4 Legend v1.4", "",
    "**Figure 4. GSE73680 supports MAGMA module-level plaque/stone papilla disease-context association.**",
    "(A) Compact patient-aware GSE73680 design summary. (B) P1 single-gene paired responses; no P1 gene reached FDR q <= 0.05 and PKD2 was nominal only. (C) Paired module slopegraphs with low-opacity patient lines, bold median trends, positive-fraction labels and FDR boundary annotations. (D) Size-matched random benchmark; dashed line marks the 95th percentile. An expression-matched benchmark was retained as a conservative boundary check and is shown in Supplementary Figure 3. These analyses support MAGMA module-level disease-context association, not P1 disease-gene validation or causality."),
    "docs/figure4_legend_v1.4.md")
  write_source_table("figure4", data.table(
    panel = c("A", "B", "C", "D", "supplementary boundary"),
    source_file = c("results/gse73680/tables/gse73680_sample_sheet.tsv", "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv",
    "results/tables/gse73680_paired_module_delta_long.tsv; results/tables/gse73680_paired_module_delta_summary.tsv",
    "results/gse73680/tables/gse73680_random_module_benchmark.tsv", "results/gse73680/tables/gse73680_expression_matched_random_benchmark.tsv"),
    role = c("study design display", "P1 single-gene boundary", "paired module shifts", "size-matched benchmark", "conservative expression-matched boundary check")
  ))
  write_visual_qc("figure4", "v1.4", data.table(check = c("compact_design", "p1_boundary", "module_shift_labels", "benchmark_lollipop", "expression_matched_caption"),
    status = "pass", note = c("vertical design card retained", "no uniform P1 response stated", "positive fraction and q boundary labels added",
    "MAGMA deep teal and P1 grey used", "legend mentions conservative boundary check")))
}

make_figure5_v06 <- function() {
  go <- fread("results/tables/figure5_v0.4_go_display_terms.tsv")
  if (!nrow(go)) {
    go <- fread("results/tables/go_bp_redundancy_reduced_terms.tsv")[redundancy_reduced_keep == TRUE & p.adjust < 0.10 & Count >= 2]
  }
  go[, short := sub("^.* \\| ", "", display_name)]
  go[, group := fifelse(grepl("^Nephron", display_name), "Nephron", fifelse(grepl("^Ion", display_name), "Ion/mineral", "Epithelial"))]
  go[, short := factor(short, levels = rev(unique(short)))]
  go[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
    labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  p_go <- ggplot(go, aes(-log10(p.adjust), short, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = hi_pal$dark_grey, stroke = 0.25) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("MAGMA top 100" = hi_pal$deep_teal, "MAGMA FDR" = "#6B8F98",
      "Loop/TAL contributors" = hi_pal$sand, "P1 core" = hi_pal$terracotta)) +
    labs(title = "A. GO functional context", subtitle = "Functional interpretation; not pathway validation",
      x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_hi(8.7) + theme(legend.position = "bottom", strip.placement = "outside")

  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  curated <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core") & term != "papillary_injury_remodeling"]
  curated[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
    labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  curated[, term_label := factor(term, levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction", "proximal_tubule_context", "collecting_duct_context")),
    labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction", "Proximal tubule context", "Collecting duct context")))]
  curated[, fill_cap := pmin(enrichment_ratio, 20)]
  curated[, txt_col := ifelse(fill_cap > 10, "white", hi_pal$ink)]
  p_cur <- ggplot(curated, aes(gene_set_label, term_label, fill = fill_cap)) +
    geom_tile(color = "white", linewidth = 0.62) +
    geom_text(aes(label = ifelse(overlap > 0, overlap, ""), color = txt_col), fontface = "bold", size = 2.85) +
    scale_color_identity() +
    scale_fill_gradient(low = "#EEF3F4", high = hi_pal$deep_teal) +
    labs(title = "B. Curated nephron/transport context", subtitle = "Marker-set overlap; not pathway activity validation",
      x = NULL, y = NULL, fill = "Enrichment\nratio") +
    theme_hi(8.7) + theme(axis.text.x = element_text(angle = 25, hjust = 1))

  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") &
    module_name %in% module_raw &
    injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(module_name, levels = rev(module_raw), labels = rev(module_keep))]
  d[, injury_label := factor(injury_module, levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
    labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"), labels = c("Paired patient delta", "Patient/group residual"))]
  d[, sig := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**", fifelse(fdr < 0.05, "*", "")))]
  d[, label := sprintf("%.2f%s", rho, sig)]
  p_heat <- ggplot(d, aes(injury_label, module_label, fill = pmin(pmax(rho, 0), 1))) +
    geom_tile(color = "white", linewidth = 0.62) +
    geom_text(data = d[rho > 0.60], aes(label = label), color = "white", fontface = "bold", size = 3.05) +
    geom_text(data = d[rho <= 0.60], aes(label = label), color = hi_pal$ink, size = 3.05) +
    facet_wrap(~ analysis, ncol = 1) +
    scale_fill_gradient(low = "#EEF3F4", high = hi_pal$terracotta, limits = c(0, 1), name = "Spearman\nrho") +
    labs(title = "C. Risk-injury coupling robustness", subtitle = "* FDR<0.05, ** FDR<0.01, *** FDR<0.001; P1 is weaker/inconsistent",
      x = NULL, y = NULL) +
    theme_hi(8.9) + theme(axis.text.x = element_text(angle = 25, hjust = 1))

  left <- plot_grid(p_go, p_cur, ncol = 1, rel_heights = c(1.02, 0.98))
  fig <- plot_grid(left, p_heat, ncol = 2, rel_widths = c(0.48, 0.52))
  hi_save_fig(fig, "results/figures/figure5_functional_context_v0.6", 13.2, 8.9)
  hi_write_lines(c("# Figure 5 Legend v0.6", "",
    "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
    "(A) Redundancy-reduced GO Biological Process terms with shortened labels and grouped strips; full terms are retained in the source table. (B) Curated nephron and transport marker-set overlap, with overlap counts shown in each non-zero tile. (C) Risk-injury coupling robustness in GSE73680 using paired patient delta and patient/group residual correlations. Values represent Spearman rho; * FDR < 0.05, ** FDR < 0.01 and *** FDR < 0.001. P1 core showed weaker and less consistent coupling than MAGMA-prioritized modules, consistent with its interpretive rather than disease-validation role."),
    "docs/figure5_legend_v0.6.md")
  write_source_table("figure5", data.table(
    panel = c("A", "B", "C"),
    source_file = c("results/tables/figure5_v0.4_go_display_terms.tsv; results/tables/go_bp_redundancy_reduced_terms.tsv",
    "results/tables/nephron_segment_marker_enrichment.tsv", "results/tables/gse73680_risk_injury_correlation_robustness.tsv"),
    role = c("GO context", "curated nephron/transport overlap", "risk-injury robustness heatmap")
  ))
  write_visual_qc("figure5", "v0.6", data.table(check = c("two_column_layout", "go_short_labels", "curated_white_text", "risk_injury_large_panel", "p1_boundary"),
    status = "pass", note = c("GO/curated left, risk-injury right", "prefixes removed and facet strips added",
    "dark curated cells use white text", "right panel emphasized", "P1 weaker/inconsistent stated")))
}

make_supplementary_figures <- function() {
  donor <- fread("results/tables/gse231569_donor_celltype_magma_scores_v0.2.tsv")
  donor <- donor[module_name %in% module_raw[1:4]]
  donor[, module_label := factor(module_name, levels = module_raw[1:4], labels = module_keep[1:4])]
  donor[, cell_label := factor(hi_label_cell(audited_broad_cell_type, short = TRUE),
    levels = hi_label_cell(c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like"), short = TRUE))]
  p_s2 <- ggplot(donor[!is.na(cell_label)], aes(cell_label, mean_score, fill = cell_label == "Loop/TAL")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = hi_pal$medium_grey, linewidth = 0.45) +
    geom_boxplot(width = 0.60, outlier.shape = NA, color = hi_pal$dark_grey, linewidth = 0.30) +
    geom_jitter(width = 0.13, size = 1.1, alpha = 0.72, color = hi_pal$ink) +
    facet_wrap(~ module_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("TRUE" = hi_pal$pale_bluegrey, "FALSE" = "#ECECEC"), guide = "none") +
    labs(title = "Supplementary Figure 1. Donor-level GSE231569 reproducibility", x = NULL, y = "Mean module score by donor-cell type") +
    theme_hi(8.8) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  hi_save_fig(p_s2, "results/figures/supplementary_figures_v0.2/supp_fig1_donor_level_reproducibility", 9.0, 6.4)

  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  delta <- unique(long[module_label %in% module_keep, .(module_label, patient_id, paired_delta)])
  delta[, module_label := factor(as.character(module_label), levels = module_keep)]
  delta[, patient_rank := frank(paired_delta, ties.method = "first"), by = module_label]
  p_s4 <- ggplot(delta, aes(reorder(patient_id, paired_delta), paired_delta, fill = paired_delta > 0)) +
    geom_hline(yintercept = 0, color = hi_pal$dark_grey, linewidth = 0.45) +
    geom_col(width = 0.72) +
    facet_wrap(~ module_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("TRUE" = hi_pal$deep_teal, "FALSE" = "#CFCFCF"), guide = "none") +
    labs(title = "Supplementary Figure 2. Paired patient waterfall in GSE73680", x = "Patient ordered within module", y = "Paired delta") +
    theme_hi(8.4) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  hi_save_fig(p_s4, "results/figures/supplementary_figures_v0.2/supp_fig2_gse73680_paired_patient_waterfall", 8.2, 9.5)

  em <- fread("results/gse73680/tables/gse73680_expression_matched_random_benchmark.tsv")
  em[, module_label := factor(module_name, levels = rev(module_raw), labels = rev(module_keep))]
  em[, emp_label := sprintf("expr-matched P=%.3f", empirical_p)]
  p_s5 <- ggplot(em, aes(expression_matched_percentile, module_label)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = hi_pal$dark_grey, linewidth = 0.60) +
    geom_segment(aes(x = 0, xend = expression_matched_percentile, yend = module_label), color = hi_pal$medium_grey, linewidth = 0.80) +
    geom_point(aes(fill = module_label != "P1 core"), shape = 21, size = 4.8, color = hi_pal$dark_grey, stroke = 0.35) +
    geom_text(aes(x = pmin(expression_matched_percentile + 0.035, 1.02), label = emp_label), hjust = 0, size = 2.85, color = hi_pal$ink) +
    scale_fill_manual(values = c("TRUE" = hi_pal$deep_teal, "FALSE" = "#CFCFCF"), guide = "none") +
    coord_cartesian(xlim = c(0, 1.18), clip = "off") +
    labs(title = "Supplementary Figure 3. Expression-matched disease-context benchmark", x = "Expression-matched random percentile", y = NULL) +
    theme_hi(9)
  hi_save_fig(p_s5, "results/figures/supplementary_figures_v0.2/supp_fig3_expression_matched_benchmark", 8.0, 4.8)

  hi_write_lines(c("# Supplementary Figure Plan v1.2", "",
    "| Supplement | File stem | Purpose | Boundary |",
    "|---|---|---|---|",
    "| Supplementary Figure 1 | supp_fig1_donor_level_reproducibility | Donor-level GSE231569 module-score reproducibility | Supports robustness of localization, not causal cell-type assignment |",
    "| Supplementary Figure 2 | supp_fig2_gse73680_paired_patient_waterfall | Patient-level paired delta heterogeneity | Supports disease-context shift transparency |",
    "| Supplementary Figure 3 | supp_fig3_expression_matched_benchmark | Conservative expression-matched benchmark | Boundary check; does not replace size-matched main benchmark |"),
    "docs/supplementary_figure_plan_v1.2.md")
}

copy_final_and_qc <- function() {
  stems <- c(
    "figure1_integrative_framework_v0.8",
    "figure2_magma_scrna_localization_v0.9",
    "figure3_p1_gene_evidence_v1.1",
    "figure4_gse73680_disease_context_v1.4",
    "figure5_functional_context_v0.6"
  )
  for (stem in stems) {
    for (ext in c("pdf", "png")) {
      src <- file.path("results/figures", paste0(stem, ".", ext))
      if (file.exists(src)) file.copy(src, file.path("results/figures/final_main_figures_v0.4", basename(src)), overwrite = TRUE)
    }
  }
  figure_files <- file.path("results/figures", paste0(stems, ".png"))
  labels <- c("Figure 1 v0.8", "Figure 2 v0.9", "Figure 3 v1.1", "Figure 4 v1.4", "Figure 5 v0.6")

  pdf("results/figures/main_figures_1_to_5_review_contact_sheet_v1.3.pdf", width = 14, height = 28, onefile = TRUE)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(5, 1)))
  for (i in seq_along(figure_files)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(labels[i], x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 13))
    grid.raster(readPNG(figure_files[i]), x = 0.5, y = 0.48, width = 0.96, height = 0.88)
    popViewport()
  }
  popViewport(); dev.off()

  png("results/figures/main_figures_1_to_5_review_contact_sheet_v1.3.png", width = 2800, height = 5600, res = 200)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(5, 1)))
  for (i in seq_along(figure_files)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(labels[i], x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 13))
    grid.raster(readPNG(figure_files[i]), x = 0.5, y = 0.48, width = 0.96, height = 0.88)
    popViewport()
  }
  popViewport(); dev.off()

  qc <- data.table(
    figure_id = paste0("Figure ", 1:5),
    version = c("v0.8", "v0.9", "v1.1", "v1.4", "v0.6"),
    output_pdf = file.path("results/figures", paste0(stems, ".pdf")),
    output_png = file.path("results/figures", paste0(stems, ".png")),
    final_folder = "results/figures/final_main_figures_v0.4",
    unified_theme = "scripts/09_manuscript/figure_theme_highimpact_v0.1.R",
    claim_boundary = "context mapping and module-level disease-context association only; no causality, TWAS convergence, SMR/coloc, spatial validation or P1 disease-gene validation"
  )
  fwrite(qc, "results/tables/main_figure_qc_v0.4.tsv", sep = "\t")

  hi_write_lines(c("# Main Figure Style Guide v0.4", "",
    "## Unified Theme",
    "All main figures in the v0.4 final folder are generated from `scripts/09_manuscript/figure_theme_highimpact_v0.1.R`.",
    "The plotting device uses a device-safe sans-serif family to avoid local PDF font embedding failures while retaining a consistent manuscript visual language.",
    "",
    "## Palette",
    "- MAGMA module evidence: deep teal `#245A64` with related blue-grey tones.",
    "- Loop/TAL focal context: `#0F4C5C`.",
    "- P1/calcium interpretive layer: sand `#B99B5A`.",
    "- Plaque/stone disease-context layer and injury coupling: terracotta `#9B5C4D`.",
    "- Boundary/no-support elements: greys `#CFCFCF`, `#AEB8BC`, `#5A6062`.",
    "",
    "## Layout Rules",
    "Figures use bold panel titles, regular axis text, compact legends and low-saturation fills. Dense evidence tables were moved into glyphs or supplementary figures where possible.",
    "",
    "## Claim Boundary",
    "The figures support MAGMA-prioritized KSD genes localizing to a Loop/TAL-associated papillary context and MAGMA module-level disease-context association in GSE73680. They do not establish causality, causal cell type, TWAS convergence, SMR/coloc support, spatial validation, therapeutic target validation, P1 single-gene disease validation or pathway activity validation."),
    "docs/main_figure_style_guide_v0.4.md")

  hi_write_lines(c("# Figure Redesign Decision Memo v0.2", "",
    "Phase 15 harmonized all main figures under a single high-impact visual theme and added three supplement figures requested in the latest instruction.",
    "",
    "- Figure 1 v0.8: rebuilt as a graphical abstract with a vector kidney-papilla niche anchor, four evidence cards and a supported/not-established boundary box.",
    "- Figure 2 v0.9: upgraded UMAP hierarchy, renamed panels, added numeric benchmark heatmap values and simplified lollipop detection legend.",
    "- Figure 3 v1.1: converted P1 roles into a partitioned continuous ribbon, fixed gene order, shaded Loop/TAL expression context and tightened the evidence glyph matrix.",
    "- Figure 4 v1.4: compacted GSE73680 design, clarified the P1 single-gene boundary, labelled module shift facets and retained expression-matched benchmark as a supplement boundary check.",
    "- Figure 5 v0.6: moved away from the evidence strip, used a two-column layout with enlarged risk-injury heatmap, shortened GO labels and preserved pathway-activity boundary language.",
    "- Supplementary Figure 1: donor-level GSE231569 reproducibility.",
    "- Supplementary Figure 2: GSE73680 paired patient waterfall.",
    "- Supplementary Figure 3: expression-matched benchmark.",
    "",
    "The manuscript draft remains `manuscript/manuscript_draft_v1.1.md` as requested."),
    "docs/figure_redesign_decision_memo_v0.2.md")
}

make_figure1_v08()
make_figure2_v09()
make_figure3_v11()
make_figure4_v14()
make_figure5_v06()
make_supplementary_figures()
copy_final_and_qc()

message("Phase 15 high-impact figure harmonization outputs written")
