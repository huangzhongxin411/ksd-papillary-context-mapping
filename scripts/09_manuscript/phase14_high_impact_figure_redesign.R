suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(png)
  library(scales)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final_main_figures_v0.3", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  deep_teal = "#245A64",
  loop_tal = "#0F4C5C",
  light_grey = "#DDE4E6",
  bg_grey = "#E6E9EA",
  medium_grey = "#BFC7CA",
  sand = "#B99B5A",
  terracotta = "#9B5C4D",
  sage = "#6F8F72",
  dark_grey = "#4D4D4D",
  text = "#303030"
)

theme_hi <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, color = pal$text),
      plot.subtitle = element_text(size = base_size, color = pal$dark_grey),
      axis.title = element_text(size = base_size + 1, color = pal$text),
      axis.text = element_text(size = base_size, color = pal$text),
      legend.title = element_text(size = base_size + 0.5, face = "bold", color = pal$text),
      legend.text = element_text(size = base_size, color = pal$text),
      strip.text = element_text(size = base_size + 0.5, face = "bold", color = pal$text),
      axis.line = element_line(linewidth = 0.55, color = pal$dark_grey),
      axis.ticks = element_line(linewidth = 0.5, color = pal$dark_grey)
    )
}

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

save_fig <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

label_cell <- function(x, short = FALSE) {
  full <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like",
    Pericyte_smooth_muscle = "Perivascular/mural-like"
  )
  short_map <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fib/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epi.",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivasc.",
    Pericyte_smooth_muscle = "Perivasc."
  )
  map <- if (short) short_map else full
  unname(ifelse(x %in% names(map), map[x], x))
}

make_figure1_v07 <- function() {
  cards <- data.table(
    x = c(0.38, 0.54, 0.70, 0.86),
    layer = c("GWAS/MAGMA", "GSE231569 snRNA", "P1 candidates", "GSE73680"),
    badge = c("57 loci | 94 Bonferroni", "Loop/TAL context", "6-gene role spectrum", "55 samples | 26 paired"),
    fill = c(pal$deep_teal, pal$loop_tal, pal$sand, pal$terracotta)
  )
  fig <- ggplot() +
    annotate("text", x = 0.60, y = 0.92, label = "Post-GWAS framework for KSD papillary context mapping",
             fontface = "bold", size = 5.0, color = pal$text) +
    annotate("path", x = c(0.08, 0.10, 0.14, 0.20, 0.24, 0.25, 0.215, 0.165, 0.115, 0.08),
             y = c(0.62, 0.78, 0.86, 0.82, 0.68, 0.54, 0.42, 0.37, 0.44, 0.62),
             color = "#837B70", linewidth = 0.7) +
    annotate("polygon", x = c(0.128, 0.158, 0.197, 0.158), y = c(0.47, 0.70, 0.55, 0.41),
             fill = "#E8D8B6", color = "#9A8B6D", linewidth = 0.45) +
    annotate("rect", xmin = 0.145, xmax = 0.205, ymin = 0.41, ymax = 0.57, fill = NA, color = pal$terracotta, linewidth = 0.55) +
    annotate("segment", x = 0.205, xend = 0.285, y = 0.49, yend = 0.35,
             arrow = arrow(length = unit(0.10, "in")), color = pal$dark_grey, linewidth = 0.55) +
    annotate("rect", xmin = 0.045, xmax = 0.30, ymin = 0.08, ymax = 0.35, fill = "#FAFBFB", color = pal$medium_grey, linewidth = 0.55) +
    annotate("path", x = c(0.075, 0.125, 0.155, 0.205, 0.252), y = c(0.15, 0.30, 0.15, 0.30, 0.15),
             color = pal$loop_tal, linewidth = 1.3) +
    annotate("rect", xmin = 0.165, xmax = 0.222, ymin = 0.18, ymax = 0.31,
             fill = "#D8C7A5", color = pal$dark_grey, linewidth = 0.45, alpha = 0.85) +
    annotate("point", x = 0.268, y = 0.145, size = 3.3, color = pal$terracotta) +
    annotate("text", x = 0.132, y = 0.105, label = "Loop/TAL", size = 2.8, fontface = "bold", color = pal$loop_tal) +
    annotate("text", x = 0.225, y = 0.315, label = "Collecting duct", size = 2.35, color = pal$dark_grey) +
    annotate("text", x = 0.268, y = 0.092, label = "plaque/\nstone", size = 2.2, color = pal$terracotta) +
    annotate("text", x = 0.157, y = 0.375, label = "kidney -> papilla niche", size = 2.8, color = pal$dark_grey) +
    geom_segment(data = cards[1:3], aes(x = x + 0.055, xend = cards$x[2:4] - 0.055, y = 0.62, yend = 0.62),
                 arrow = arrow(length = unit(0.10, "in")), color = pal$dark_grey, linewidth = 0.55) +
    geom_label(data = cards, aes(x = x, y = 0.69, label = layer, fill = fill), color = "white",
               fontface = "bold", size = 3.0, label.padding = unit(0.22, "lines"), linewidth = 0) +
    geom_label(data = cards, aes(x = x, y = 0.555, label = badge), fill = "white", color = pal$text,
               size = 2.65, label.padding = unit(0.22, "lines"), label.size = 0.35, lineheight = 0.95) +
    annotate("rect", xmin = 0.34, xmax = 0.95, ymin = 0.16, ymax = 0.39, fill = "#F7F8F8", color = pal$medium_grey, linewidth = 0.55) +
    annotate("text", x = 0.37, y = 0.345, label = "Supported inference", hjust = 0, fontface = "bold", size = 3.2, color = pal$deep_teal) +
    annotate("text", x = 0.37, y = 0.295, label = "Loop/TAL-associated cellular context\nMAGMA module-level disease-context association",
             hjust = 0, vjust = 1, size = 2.7, color = pal$text, lineheight = 0.92) +
    annotate("text", x = 0.67, y = 0.345, label = "Not established", hjust = 0, fontface = "bold", size = 3.2, color = pal$terracotta) +
    annotate("text", x = 0.67, y = 0.305, label = "causality | TWAS convergence\ncolocalization | spatial validation\nP1 disease-gene validation",
             hjust = 0, vjust = 1, size = 2.55, color = pal$text, lineheight = 0.92) +
    annotate("text", x = 0.37, y = 0.178, label = "Resource-limited extensions: TWAS / SMR-coloc / spatial audited, not used as evidence layers",
             hjust = 0, size = 2.25, color = pal$dark_grey) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.06, 0.96), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  save_fig(fig, "results/figures/figure1_integrative_framework_v0.7", 13.2, 5.9)
  write_lines(c("# Figure 1 Legend v0.7", "",
                "**Figure 1. Post-GWAS framework for KSD papillary context mapping.**",
                "Four card-based evidence layers link GWAS/MAGMA prioritization, GSE231569 snRNA localization, P1 candidate interpretation and GSE73680 plaque/stone papilla disease-context association. The supported inference and not-established claims are separated in the claim-boundary box; resource-limited TWAS, SMR-coloc and spatial analyses are not used as evidence layers."),
              "docs/figure1_legend_v0.7.md")
}

make_figure2_v08 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := label_cell(audited_broad_cell_type)]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(14)
  bg <- if (nrow(top50) > 26000) top50[sample(.N, 26000)] else copy(top50)
  non_tal <- bg[is_tal == FALSE]
  tal <- top50[is_tal == TRUE]
  cell_cols <- c("Loop/TAL" = pal$loop_tal, "Collecting duct" = "#D9DEE2",
                 "Endothelial" = "#8BA7B0", "Fibroblast/stromal" = "#A9B8A5",
                 "Injured epithelial" = "#C7A39A", "Perivascular/mural-like" = "#B9A9C8")
  center <- tal[, .(x = median(UMAP_1), y = median(UMAP_2))]
  p_a <- ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2), color = pal$bg_grey, alpha = 0.22, size = 0.14) +
    geom_point(data = non_tal, aes(UMAP_1, UMAP_2, color = cell_label), alpha = 0.38, size = 0.15) +
    geom_point(data = tal, aes(UMAP_1, UMAP_2, color = cell_label), alpha = 0.92, size = 0.42) +
    stat_ellipse(data = tal, aes(UMAP_1, UMAP_2), color = pal$dark_grey, linewidth = 0.55, level = 0.85) +
    annotate("label", x = center$x + 3.1, y = center$y + 3.6, label = "Loop/TAL\nn = 540",
             hjust = 0, size = 3.0, fill = "white", color = pal$loop_tal, fontface = "bold") +
    scale_color_manual(values = cell_cols, breaks = names(cell_cols), name = "Audited cell type") +
    labs(title = "A. Highlight UMAP of audited GSE231569 atlas", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9) + theme(legend.position = "right")
  qlim <- quantile(bg$celllevel_module_score, c(0.02, 0.99), na.rm = TRUE)
  p_b <- ggplot(bg, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celllevel_module_score), size = 0.16, alpha = 0.78) +
    stat_density_2d(data = tal, aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                    color = pal$dark_grey, linewidth = 0.55, bins = 4) +
    scale_color_gradientn(colors = c(pal$light_grey, "#8CA7AF", pal$deep_teal),
                          limits = qlim, oob = squish, name = "MAGMA top 50\nscore") +
    labs(title = "B. MAGMA top 50 score", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9)
  bench <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  bench <- bench[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  bench[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
                                 labels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
  cell_order <- c("Loop/TAL", "Collecting duct", "Injured epithelial", "Endothelial", "Fibroblast/stromal", "Perivascular/mural-like")
  bench[, cell_label := factor(label_cell(audited_broad_cell_type), levels = cell_order)]
  bench <- bench[!is.na(cell_label)]
  bench[, label := fifelse(cell_label == "Loop/TAL" & benchmark_percentile >= 0.95, ">95th", "")]
  p_c <- ggplot(bench, aes(cell_label, module_label, fill = benchmark_percentile)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), color = "white", fontface = "bold", size = 3.0) +
    scale_fill_gradientn(colors = c("#F0F2F2", "#9BB0B7", pal$deep_teal), limits = c(0, 1), name = "Benchmark\npercentile") +
    labs(title = "C. Benchmark percentile across audited contexts", x = NULL, y = NULL) +
    theme_hi(8.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = pal$medium_grey, linewidth = 0.7) +
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = pal$dark_grey, stroke = 0.35) +
    scale_fill_manual(values = c("P1 candidate" = pal$sand, "Other MAGMA gene" = pal$deep_teal)) +
    scale_size_continuous(range = c(2.4, 5.0), labels = percent_format()) +
    coord_cartesian(xlim = c(0, 0.78)) +
    labs(title = "D. Leading contributors to Loop/TAL-associated signal", x = "Contribution score", y = NULL, fill = NULL, size = "Detection") +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.08, 0.92))
  save_fig(fig, "results/figures/figure2_magma_scrna_localization_v0.8", 13.2, 9.4)
  write_lines(c("# Figure 2 Legend v0.8", "",
                "**Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated renal papillary single-nucleus context.**",
                "(A) High-impact highlight UMAP of audited GSE231569 cell contexts. (B) MAGMA top 50 score UMAP using a single low-saturation continuous scale; cell-level score is visualization only and donor-cell-type summaries are used for interpretation. (C) Compact benchmark heatmap showing size-matched random benchmark percentiles across audited contexts; >95th marks Loop/TAL benchmark support. (D) Ranked lollipop plot of leading contributors. Contribution score is a descriptive ranking based on MAGMA prioritization, Loop/TAL expression preference and detection support, not causal driver inference."),
              "docs/figure2_legend_v0.8.md")
}

make_figure3_v10 <- function() {
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  ribbon <- data.table(
    x = c(1, 2, 3, 3.55, 4.45, 5.35),
    gene = gene_order,
    role = c("TAL identity", "Transport", "Ion/calcium\nhandling", "Ion/calcium\nhandling", "Supporting\ncontext", "Broad epithelial\ncontext"),
    y = 1
  )
  p_a <- ggplot(ribbon, aes(x, y)) +
    annotate("segment", x = 0.75, xend = 5.65, y = 1, yend = 1, linewidth = 6.5, color = "#EEF3F4", lineend = "round") +
    annotate("segment", x = 0.75, xend = 5.65, y = 1, yend = 1, linewidth = 0.7, color = pal$medium_grey,
             arrow = arrow(length = unit(0.13, "in"))) +
    geom_label(aes(label = paste0(gene, "\n", role), fill = role), color = "white", fontface = "bold",
               size = 2.75, label.padding = unit(0.22, "lines"), linewidth = 0, lineheight = 0.88) +
    scale_fill_manual(values = c("TAL identity" = pal$loop_tal, "Transport" = "#557F89",
                                 "Ion/calcium\nhandling" = pal$sand, "Supporting\ncontext" = "#6F8F98",
                                 "Broad epithelial\ncontext" = pal$terracotta)) +
    annotate("text", x = 0.75, y = 1.42, label = "TAL identity", hjust = 0, size = 3.0, color = pal$dark_grey) +
    annotate("text", x = 5.65, y = 1.42, label = "broader epithelial context", hjust = 1, size = 3.0, color = pal$dark_grey) +
    coord_cartesian(xlim = c(0.55, 5.85), ylim = c(0.55, 1.55), clip = "off") +
    labs(title = "A. P1 role spectrum ribbon") +
    theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 12, color = pal$text), legend.position = "none")
  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(label_cell(cell_type, short = TRUE), levels = label_cell(keep_cells, short = TRUE))]
  cell[, gene := factor(gene, levels = rev(gene_order))]
  bg <- data.table(cell_label = factor("Loop/TAL", levels = levels(cell$cell_label)), xmin = 0.5, xmax = 1.5)
  p_b <- ggplot(cell, aes(cell_label, gene)) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#EFF5F6", color = NA) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = pal$medium_grey, stroke = 0.35) +
    scale_fill_gradientn(colors = c("#EEF3F4", "#83A0A8", pal$deep_teal), name = "Average\nexpression") +
    scale_size_continuous(range = c(0.9, 5.0), labels = percent_format(), name = "Detected") +
    labs(title = "B. P1 expression across audited cell types", x = NULL, y = NULL) +
    theme_hi(8.8) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom", axis.line = element_blank(), axis.ticks = element_blank())
  p1[, gene := factor(gene, levels = rev(gene_order))]
  p1[, emph := as.character(gene) %in% c("CLDN10", "CLDN14", "CASR")]
  p_c <- ggplot(p1, aes(log2(specificity_ratio_avg), gene)) +
    geom_segment(aes(x = 0, xend = log2(specificity_ratio_avg), yend = gene), color = pal$medium_grey, linewidth = 0.65) +
    geom_point(aes(fill = specificity_class, alpha = emph), shape = 21, size = 4.0, color = pal$dark_grey, stroke = 0.35) +
    geom_text(data = p1[emph == TRUE], aes(label = as.character(gene)), nudge_y = 0.26, size = 2.7, fontface = "bold", color = pal$dark_grey) +
    scale_fill_manual(values = c(strong_TAL_preferential = pal$deep_teal, moderate_TAL_preferential = pal$sand), na.value = pal$medium_grey) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.55), guide = "none") +
    labs(title = "C. log2(TAL specificity ratio)", x = "log2(TAL specificity ratio)", y = NULL, fill = "Specificity") +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  ev <- data.table(
    gene = factor(gene_order, levels = rev(gene_order)),
    MAGMA = c("strong", "strong", "strong", "strong", "strong", "strong"),
    `TAL specificity` = c("moderate", "strong", "strong", "strong", "strong", "moderate"),
    `Donor detection` = "contextual",
    `Bulk response` = c("no FDR", "no FDR", "no FDR", "no FDR", "no FDR", "nominal"),
    Role = c("Representative TAL", "Transport", "Ion handling", "Calcium sensing", "Supporting", "Broad epithelial")
  )
  ev_long <- melt(ev, id.vars = "gene", variable.name = "evidence", value.name = "call")
  ev_long[, evidence := factor(evidence, levels = c("MAGMA", "TAL specificity", "Donor detection", "Bulk response", "Role"))]
  glyph <- ev_long[evidence %in% c("MAGMA", "TAL specificity", "Donor detection")]
  text_dt <- ev_long[!evidence %in% c("MAGMA", "TAL specificity", "Donor detection")]
  p_d <- ggplot() +
    geom_tile(data = ev_long, aes(evidence, gene), fill = "#F7F8F8", color = "white", linewidth = 0.6) +
    geom_point(data = glyph[call == "strong"], aes(evidence, gene), shape = 21, fill = pal$deep_teal, color = pal$deep_teal, size = 4.5) +
    geom_point(data = glyph[call == "moderate"], aes(evidence, gene), shape = 21, fill = pal$sand, color = pal$dark_grey, size = 4.5, stroke = 0.35) +
    geom_point(data = glyph[call == "contextual"], aes(evidence, gene), shape = 21, fill = "white", color = pal$dark_grey, size = 4.5, stroke = 0.7) +
    geom_text(data = text_dt, aes(evidence, gene, label = call), size = 2.55, color = pal$text) +
    labs(title = "D. Evidence glyph matrix", x = NULL, y = NULL) +
    theme_hi(8.7) + theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.text.y = element_text(face = "italic"),
                          axis.line = element_blank(), axis.ticks = element_blank())
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.78, 1.22))
  save_fig(fig, "results/figures/figure3_p1_gene_evidence_v1.0", 13.2, 8.9)
  write_lines(c("# Figure 3 Legend v1.0", "",
                "**Figure 3. P1 genes form a TAL, transport and calcium-handling role spectrum.**",
                "(A) Horizontal role-spectrum ribbon. (B) Dotplot of P1 expression and detection across audited cell types, with Loop/TAL highlighted. (C) log2(TAL specificity ratio), emphasizing CLDN10, CLDN14 and CASR. (D) Evidence glyph matrix: filled deep-teal circle = strong support, sand circle = moderate support, open circle = contextual/detectable support. Bulk response denotes GSE73680 single-gene plaque/stone paired response; no P1 gene reached FDR q <= 0.05, and PKD2 was nominal only."),
              "docs/figure3_legend_v1.0.md")
}

make_figure4_v13 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")
  timeline <- data.table(
    y = c(0.82, 0.60, 0.38, 0.16),
    title = c("GSE73680", "55 samples / 29 patients", "26 paired patients", "patient-aware limma + paired delta"),
    fill = c(pal$deep_teal, pal$light_grey, pal$light_grey, pal$terracotta)
  )
  p_a <- ggplot(timeline, aes(0.5, y)) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0.75, yend = 0.23), color = pal$medium_grey, linewidth = 0.8,
                 arrow = arrow(length = unit(0.09, "in"))) +
    geom_label(aes(label = title, fill = fill), color = c("white", pal$text, pal$text, "white"),
               fontface = "bold", size = 3.4, label.padding = unit(0.28, "lines"), linewidth = 0) +
    scale_fill_identity() +
    labs(title = "A. Patient-aware GSE73680 design") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.03, 0.95)) +
    theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 12, color = pal$text))
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1[, gene := factor(gene, levels = rev(gene_order))]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No FDR support")]
  p1[, label_x := ifelse(paired_delta >= 0, paired_delta + 0.035, paired_delta - 0.035)]
  p1[, hjust := ifelse(paired_delta >= 0, 0, 1)]
  p_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$medium_grey, linewidth = 0.6) +
    geom_col(width = 0.58, color = pal$dark_grey, linewidth = 0.3) +
    geom_text(aes(x = label_x, label = sprintf("P=%.3f", p_value), hjust = hjust), size = 2.75, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$sand, "No FDR support" = "#CFCFCF")) +
    coord_cartesian(xlim = c(-0.42, 0.83)) +
    labs(title = "B. P1 single-gene response", subtitle = "No P1 gene reached FDR q <= 0.05; PKD2 nominal only",
         x = "Patient-level paired delta", y = NULL, fill = NULL) +
    theme_hi(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
  long <- long[module_label %in% keep]
  long[, module_label := factor(as.character(module_label), levels = keep)]
  lab <- sumdt[module_label %in% keep, .(module_label = factor(as.character(module_label), levels = keep),
                                        label = sprintf("%d/%d increased", n_positive_delta, n_paired_patients))]
  p_c <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(color = pal$medium_grey, alpha = 0.20, linewidth = 0.28) +
    geom_point(aes(fill = group_label), shape = 21, color = "white", stroke = 0.10, size = 1.08, alpha = 0.78) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.25, color = pal$text) +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2.25, color = pal$text) +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.12,
              size = 3.05, fontface = "bold", color = pal$text) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = pal$sand)) +
    labs(title = "C. Paired patient module shifts", x = NULL, y = "Module score", fill = NULL) +
    theme_hi(8.4) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
  bench <- bench[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  bench[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                                 labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  bench[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  p_d <- ggplot(bench, aes(percentile, module_label)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$dark_grey, linewidth = 0.6) +
    geom_segment(aes(x = 0, xend = percentile, yend = module_label), color = pal$medium_grey, linewidth = 0.8) +
    geom_point(aes(fill = percentile >= 0.95), shape = 21, size = 4.7, color = pal$dark_grey, stroke = 0.35) +
    geom_text(aes(x = pmin(percentile + 0.035, 1.03), label = emp_label), hjust = 0, size = 2.8, color = pal$text) +
    scale_fill_manual(values = c("TRUE" = pal$deep_teal, "FALSE" = "#CFCFCF"),
                      labels = c("Background-like", "Exceeds 95th percentile"), name = NULL) +
    coord_cartesian(xlim = c(0, 1.12)) +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL) +
    theme_hi(9) + theme(legend.position = "bottom")
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.82, 1.18))
  save_fig(fig, "results/figures/figure4_gse73680_disease_context_v1.3", 13.2, 9.3)
  write_lines(c("# Figure 4 Legend v1.3", "",
                "**Figure 4. GSE73680 supports MAGMA module-level plaque/stone papilla disease-context association.**",
                "(A) Vertical patient-aware design timeline. (B) P1 single-gene paired responses with P values shown outside bars; no P1 gene reached FDR q <= 0.05, and PKD2 was nominal only. (C) Paired module slopegraphs with low-opacity patient lines and bold median trends. (D) Lollipop benchmark plot showing size-matched random percentile; dashed line marks the 95th percentile. These analyses support MAGMA module-level disease-context association, not P1 disease-gene validation or causality."),
              "docs/figure4_legend_v1.3.md")
}

make_figure5_v05 <- function() {
  strip <- data.table(
    evidence = rep(c("MAGMA", "Loop/TAL", "GSE73680", "Functional"), each = 4),
    gene_set = rep(c("Top 100", "FDR", "Loop/TAL contrib.", "P1 core"), times = 4),
    call = c("+++", "+++", "++", "+", "++", "++", "+++", "++", "++", "++", "NA", "+", "++", "++", "++", "++")
  )
  strip[, gene_set := factor(gene_set, levels = c("Top 100", "FDR", "Loop/TAL contrib.", "P1 core"))]
  strip[, evidence := factor(evidence, levels = c("MAGMA", "Loop/TAL", "GSE73680", "Functional"))]
  strip[, fill_class := fcase(call == "+++", "strong", call == "++", "moderate", call == "+", "support", default = "na")]
  p_strip <- ggplot(strip, aes(gene_set, evidence, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 2.8, color = pal$text) +
    scale_fill_manual(values = c(strong = pal$deep_teal, moderate = "#7E9AA2", support = pal$sand, na = "#CFCFCF")) +
    labs(title = "Evidence summary", x = NULL, y = NULL) +
    theme_hi(8.2) + theme(axis.text.x = element_text(angle = 0), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none")
  go <- fread("results/tables/figure5_v0.4_go_display_terms.tsv")
  if (!nrow(go)) {
    go <- fread("results/tables/go_bp_redundancy_reduced_terms.tsv")[redundancy_reduced_keep == TRUE & p.adjust < 0.10 & Count >= 2]
  }
  go[, short := sub("^.* \\| ", "", display_name)]
  go[, theme := fifelse(grepl("^Nephron", display_name), "Nephron", fifelse(grepl("^Ion", display_name), "Ion/mineral", "Epithelial"))]
  go[, short := factor(short, levels = rev(unique(short)))]
  go[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  p_go <- ggplot(go, aes(-log10(p.adjust), short, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = pal$dark_grey, stroke = 0.25) +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("MAGMA top 100" = pal$deep_teal, "MAGMA FDR" = "#7E9AA2",
                                 "Loop/TAL contributors" = pal$sand, "P1 core" = pal$terracotta)) +
    labs(title = "A. GO functional context", subtitle = "Functional interpretation; not pathway validation",
         x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_hi(8.6) + theme(legend.position = "bottom", strip.placement = "outside")
  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  curated <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core") & term != "papillary_injury_remodeling"]
  curated[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                     labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  curated[, term_label := factor(term, levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction", "proximal_tubule_context", "collecting_duct_context")),
                                 labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction", "Proximal tubule context", "Collecting duct context")))]
  p_cur <- ggplot(curated, aes(gene_set_label, term_label, fill = pmin(enrichment_ratio, 20))) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(overlap > 0, overlap, "")), size = 2.8, color = pal$text) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$deep_teal) +
    labs(title = "B. Curated nephron/transport context", subtitle = "Curated marker-set overlap; not pathway activity validation",
         x = NULL, y = NULL, fill = "Enrichment\nratio") +
    theme_hi(8.6) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") &
                module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates") &
                injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                             labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  d[, injury_label := factor(injury_module, levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                             labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"),
                         labels = c("Paired patient delta", "Patient/group residual"))]
  d[, sig := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**", fifelse(fdr < 0.05, "*", "")))]
  d[, label := sprintf("%.2f%s", rho, sig)]
  p_heat <- ggplot(d, aes(injury_label, module_label, fill = pmin(pmax(rho, 0), 1))) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(data = d[rho > 0.65], aes(label = label), color = "white", fontface = "bold", size = 3.0) +
    geom_text(data = d[rho <= 0.65], aes(label = label), color = pal$text, size = 3.0) +
    facet_wrap(~ analysis, ncol = 1) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$terracotta, limits = c(0, 1), name = "Spearman\nrho") +
    labs(title = "C. Risk-injury coupling robustness", x = NULL, y = NULL) +
    theme_hi(8.8) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  lower <- plot_grid(plot_grid(p_go, p_cur, ncol = 1, rel_heights = c(1.05, 0.95)), p_heat, ncol = 2, rel_widths = c(0.48, 0.52))
  fig <- plot_grid(p_strip, lower, ncol = 1, rel_heights = c(0.27, 1))
  save_fig(fig, "results/figures/figure5_functional_context_v0.5", 13.2, 9.3)
  write_lines(c("# Figure 5 Legend v0.5", "",
                "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
                "The top strip summarizes evidence levels (+++ strong support, ++ moderate support, + contextual/detectable support, NA not applicable). (A) Redundancy-reduced GO Biological Process terms with shortened display labels; full terms are retained in the source table. (B) Curated nephron/transport marker-set overlap. (C) Risk-injury coupling robustness in GSE73680 using paired patient delta and patient/group residual correlations. Values represent Spearman rho; * FDR < 0.05, ** FDR < 0.01 and *** FDR < 0.001. P1 core showed weaker and less consistent coupling than MAGMA-prioritized modules, consistent with its interpretive rather than disease-validation role."),
              "docs/figure5_legend_v0.5.md")
}

copy_final_and_qc <- function() {
  stems <- c(
    figure1 = "figure1_integrative_framework_v0.7",
    figure2 = "figure2_magma_scrna_localization_v0.8",
    figure3 = "figure3_p1_gene_evidence_v1.0",
    figure4 = "figure4_gse73680_disease_context_v1.3",
    figure5 = "figure5_functional_context_v0.5"
  )
  for (stem in stems) {
    for (ext in c("pdf", "png")) {
      src <- file.path("results/figures", paste0(stem, ".", ext))
      if (file.exists(src)) file.copy(src, file.path("results/figures/final_main_figures_v0.3", basename(src)), overwrite = TRUE)
    }
  }
  figure_files <- file.path("results/figures", paste0(stems, ".png"))
  labels <- c("Figure 1 v0.7", "Figure 2 v0.8", "Figure 3 v1.0", "Figure 4 v1.3", "Figure 5 v0.5")
  pdf("results/figures/main_figures_1_to_5_review_contact_sheet_v1.2.pdf", width = 14, height = 28, onefile = TRUE)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(5, 1)))
  for (i in seq_along(figure_files)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(labels[i], x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 13))
    grid.raster(readPNG(figure_files[i]), x = 0.5, y = 0.48, width = 0.96, height = 0.88)
    popViewport()
  }
  popViewport(); dev.off()
  png("results/figures/main_figures_1_to_5_review_contact_sheet_v1.2.png", width = 2800, height = 5600, res = 200)
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
    version = c("v0.7", "v0.8", "v1.0", "v1.3", "v0.5"),
    typography = "panel title bold; axes and legends regular hierarchy",
    palette = "unified low-saturation deep teal, Loop/TAL teal, sand, terracotta, sage/grey",
    redesign_focus = c("schematic cards and badges", "highlight UMAP plus benchmark heatmap",
                       "role ribbon and evidence glyph matrix", "timeline and lollipop benchmark",
                       "reduced checklist and enlarged coupling heatmap"),
    claim_boundary = "no causality, TWAS convergence, colocalization, spatial validation or P1 disease-gene validation"
  )
  fwrite(qc, "results/tables/figure_style_qc_v0.3.tsv", sep = "\t")
  write_lines(c(
    "# Main Figure Style Guide v0.3", "",
    "## Typography",
    "Panel labels/titles are bold. Axis titles, tick labels and legends remain regular-weight. Figures use a consistent sans-serif R graphics style with larger manuscript-readable labels.",
    "",
    "## Palette",
    "- MAGMA/significant module: #245A64",
    "- Loop/TAL target cells: #0F4C5C",
    "- Other audited/background cells: #DDE4E6, #E6E9EA and low-saturation compartment colors",
    "- P1/calcium-related: #B99B5A",
    "- Disease/plaque-stone: #9B5C4D",
    "- Injury/remodeling context: #6F8F72 where used",
    "- Boundary/no support: #4D4D4D and #CFCFCF",
    "",
    "## Claim Boundary",
    "Figures support context mapping and module-level disease-context association. They do not establish causality, TWAS convergence, colocalization, spatial validation or P1 disease-gene validation."
  ), "docs/main_figure_style_guide_v0.3.md")
  write_lines(c(
    "# Figure Redesign Decision Memo v0.1", "",
    "Phase 14 upgraded the figures from project-report style toward high-impact manuscript style while preserving the claim boundary.",
    "",
    "- Figure 1: redesigned as a card-based visual framework with data badges and a two-column claim boundary.",
    "- Figure 2: changed to a highlight UMAP, single-scale score UMAP, compact benchmark heatmap and ranked lollipop contributor plot.",
    "- Figure 3: changed role cards to a spectrum ribbon and replaced evidence symbols with a glyph matrix.",
    "- Figure 4: changed patient design to a vertical timeline and benchmark bars to lollipop statistics.",
    "- Figure 5: reduced checklist prominence, shortened GO labels and enlarged risk-injury coupling heatmap.",
    "",
    "No new causal, TWAS, colocalization or spatial claims were introduced."
  ), "docs/figure_redesign_decision_memo_v0.1.md")
}

make_figure1_v07()
make_figure2_v08()
make_figure3_v10()
make_figure4_v13()
make_figure5_v05()
copy_final_and_qc()

message("Phase 14 high-impact figure redesign outputs written")
