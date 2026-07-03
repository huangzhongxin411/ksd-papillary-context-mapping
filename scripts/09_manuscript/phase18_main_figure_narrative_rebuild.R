suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(ggalluvial)
  library(scales)
  library(grid)
})

source("scripts/09_manuscript/figure_theme_final_candidate.R")
source("scripts/09_manuscript/figure2a_ksd_gwas_manhattan_publication.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final_main_figures_candidate", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

module_raw <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
module_keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
role_cols <- c(
  "TAL identity" = "#0E5A64", "Transport" = "#6F929B", "Ion/calcium" = "#C69A35",
  "Supporting context" = "#6F929B", "Broad epithelial" = "#A65A49"
)
injury_col <- "#6F8F72"

save_main <- function(plot, stem, width, height) {
  hi_save_fig(plot, file.path("results/figures", stem), width, height)
  for (ext in c("pdf", "png", "svg")) {
    file.copy(file.path("results/figures", paste0(stem, ".", ext)),
      file.path("results/figures/final_main_figures_candidate", paste0(stem, ".", ext)), overwrite = TRUE)
  }
}

write_figure_files <- function(id, legend, sources, qc) {
  hi_write_lines(legend, file.path("docs", paste0("figure", id, "_legend_final_candidate.md")))
  fwrite(sources, file.path("results/tables", paste0("figure", id, "_panel_source_files.tsv")), sep = "\t")
  fwrite(qc, file.path("results/tables", paste0("figure", id, "_visual_qc.tsv")), sep = "\t")
}

make_figure1 <- function() {
  modules <- data.table(
    x = c(42, 57, 72, 87), y = c(62, 67, 62, 67), n = 1:4,
    title = c("GWAS / MAGMA", "GSE231569 snRNA", "P1 candidates", "GSE73680 bulk"),
    metric = c("57 loci | 94 genes", "Loop/TAL context", "6-gene spectrum", "module association"),
    fill = c("#0E5A64", "#6F929B", "#C69A35", "#A65A49")
  )
  connectors <- data.table(x = modules$x[-4] + 4.4, xend = modules$x[-1] - 4.4,
    y = modules$y[-4], yend = modules$y[-1])
  papilla <- data.table(
    x = c(8, 10, 13, 16, 18, 17, 14, 12, 14, 18, 21, 23),
    y = c(31, 48, 42, 56, 75, 84, 86, 78, 67, 55, 42, 24)
  )
  gene_beads <- data.table(x = seq(68.2, 75.8, length.out = 6), y = 61.4,
    gene = gene_order, size = c(2.2, 2.5, 2.5, 3.1, 2.3, 2.1))
  shift <- data.table(y = c(64.5, 66.5, 68.5), x = 84.3, xend = c(88.1, 89.1, 88.6))

  p <- ggplot() +
    annotate("text", x = 50, y = 95, label = "Post-GWAS evidence map for kidney stone papillary cellular context",
      fontface = "bold", size = 7.2, color = hi_pal$ink) +
    annotate("text", x = 50, y = 91, label = "Genetic prioritization to cellular localization and disease-context association",
      size = 3.5, color = hi_pal$bluegrey) +
    annotate("label", x = 5, y = 82.5, label = "RENAL PAPILLARY NICHE", hjust = 0,
      size = 3.1, fontface = "bold", color = hi_pal$deep_teal, fill = alpha("white", 0.92), linewidth = 0) +
    annotate("path", x = c(4, 29, 29, 4, 4), y = c(19, 19, 86, 86, 19),
      color = hi_pal$light_grey, linewidth = 1.0) +
    geom_path(data = papilla, aes(x, y), color = hi_pal$light_grey, linewidth = 8.5, lineend = "round") +
    geom_path(data = papilla, aes(x, y), color = hi_pal$deep_teal, linewidth = 2.4, lineend = "round") +
    annotate("label", x = 5, y = 82.5, label = "RENAL PAPILLARY NICHE", hjust = 0,
      size = 3.1, fontface = "bold", color = hi_pal$deep_teal, fill = alpha("white", 0.96), linewidth = 0) +
    annotate("segment", x = 25, xend = 25, y = 27, yend = 80, color = hi_pal$sand, linewidth = 3.0, lineend = "round") +
    annotate("segment", x = 25, xend = 22.5, y = 62, yend = 66, color = hi_pal$sand, linewidth = 2.0, lineend = "round") +
    annotate("segment", x = 25, xend = 27.5, y = 50, yend = 54, color = hi_pal$sand, linewidth = 2.0, lineend = "round") +
    annotate("point", x = c(23, 26.5), y = c(24, 22.5), shape = 23, size = c(5.0, 4.2),
      fill = hi_pal$terracotta, color = hi_pal$terracotta) +
    annotate("text", x = 6, y = 53, label = "Loop/TAL", hjust = 0, size = 3.4,
      fontface = "bold", color = hi_pal$deep_teal) +
    annotate("text", x = 27, y = 75, label = "Collecting\nduct", hjust = 0, size = 2.8,
      fontface = "bold", color = hi_pal$sand) +
    annotate("text", x = 19, y = 19.8, label = "Plaque / stone", hjust = 0, size = 2.8,
      fontface = "bold", color = hi_pal$terracotta) +
    annotate("text", x = 36, y = 84, label = "FOUR EVIDENCE LAYERS", hjust = 0,
      size = 3.1, fontface = "bold", color = hi_pal$ink) +
    geom_curve(data = connectors, aes(x = x, xend = xend, y = y, yend = yend),
      curvature = 0.12, color = hi_pal$medium_grey, linewidth = 0.8,
      arrow = arrow(length = unit(0.10, "in"))) +
    geom_point(data = modules, aes(x, y, fill = fill), shape = 21, size = 18,
      color = "white", stroke = 1.2) +
    geom_text(data = modules, aes(x, y + 0.2, label = n), color = "white", fontface = "bold", size = 4.0) +
    geom_text(data = modules, aes(x, y + 11, label = title), fontface = "bold", size = 3.2, color = hi_pal$ink) +
    geom_text(data = modules, aes(x, y - 10, label = metric), fontface = "bold", size = 2.8, color = hi_pal$ink) +
    scale_fill_identity() +
    annotate("segment", x = 38.3, xend = 45.7, y = 61.0, yend = 61.0, color = hi_pal$deep_teal, linewidth = 0.5) +
    annotate("segment", x = c(38.8, 40.1, 41.2, 42.0, 43.4, 44.4, 45.0),
      xend = c(38.8, 40.1, 41.2, 42.0, 43.4, 44.4, 45.0), y = 61,
      yend = c(63, 66, 64, 69, 63.5, 67, 64.5), color = hi_pal$deep_teal, linewidth = 0.75) +
    annotate("point", x = c(54.5, 56.2, 57.3, 58.5, 59.0), y = c(66, 69, 64.5, 68, 65.2),
      size = c(2.0, 2.4, 3.3, 2.1, 2.0), color = c(rep(hi_pal$light_grey, 2), hi_pal$deep_teal, hi_pal$light_grey, hi_pal$deep_teal)) +
    geom_segment(data = gene_beads[-6], aes(x = x, xend = gene_beads$x[-1], y = y, yend = y), color = hi_pal$sand, linewidth = 0.6) +
    geom_point(data = gene_beads, aes(x, y, size = size), fill = "white", color = hi_pal$sand, shape = 21, stroke = 1.2) +
    scale_size_identity() +
    geom_segment(data = shift, aes(x = x, xend = xend, y = y, yend = y), color = hi_pal$terracotta,
      linewidth = 1.1, arrow = arrow(length = unit(0.08, "in"))) +
    annotate("rect", xmin = 33, xmax = 97, ymin = 5.8, ymax = 16.8, fill = "#F8FAFA",
      color = hi_pal$medium_grey, linewidth = 0.55) +
    annotate("rect", xmin = 33, xmax = 60, ymin = 5.8, ymax = 16.8, fill = alpha(hi_pal$deep_teal, 0.08), color = NA) +
    annotate("rect", xmin = 60, xmax = 97, ymin = 5.8, ymax = 16.8, fill = alpha(hi_pal$terracotta, 0.06), color = NA) +
    annotate("text", x = 35, y = 14.4, label = "SUPPORTED INFERENCE", hjust = 0, fontface = "bold", size = 3.1, color = hi_pal$deep_teal) +
    annotate("text", x = 35, y = 11.5, label = "Loop/TAL-associated context\nMAGMA module-level disease-context association",
      hjust = 0, vjust = 0.5, size = 2.30, lineheight = 0.95, color = hi_pal$ink) +
    annotate("text", x = 62, y = 14.4, label = "NOT ESTABLISHED", hjust = 0, fontface = "bold", size = 3.1, color = hi_pal$terracotta) +
    annotate("text", x = 62, y = 11.5, label = "causality | TWAS | SMR/coloc\nspatial validation | P1 disease validation",
      hjust = 0, vjust = 0.5, size = 2.25, lineheight = 0.95, color = hi_pal$ink) +
    coord_cartesian(xlim = c(2, 98), ylim = c(3, 98), clip = "off") +
    theme_hi_void(10)

  save_main(p, "figure1_evidence_map_final_candidate", 13.2, 7.6)
  write_figure_files(1,
    c("# Figure 1 Legend Final Candidate", "", "**Figure 1. Graphical evidence map for post-GWAS localization of kidney stone risk to a renal papillary cellular context.**",
      "A renal papillary niche schematic anchors four linked evidence layers: GWAS/MAGMA prioritization, audited GSE231569 single-nucleus localization, six-gene P1 interpretation and patient-aware GSE73680 module-level disease-context association. Mini glyphs are conceptual and do not replace quantitative panels. The bottom ribbon separates supported inference from boundaries that remain unestablished, including causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation."),
    data.table(panel = c("niche", "evidence arc", "boundary ribbon"),
      source_file = c("conceptual renal papillary schematic", "project evidence synthesis", "project claim-boundary audit"),
      role = c("anatomical orientation", "four-layer study framework", "supported versus unestablished claims")),
    data.table(check = c("graphical_abstract_layout", "axis_free_glyphs", "three_color_semantics", "claim_boundary", "editable_svg"),
      status = "pass", note = c("30/70 niche-to-evidence structure", "no quantitative axes in mini glyphs", "teal/sand/terracotta plus grey",
        "supported and not established separated", "all layers generated as vector ggplot elements")))
}

make_figure2 <- function() {
  manhattan <- build_figure2a_manhattan(save_outputs = TRUE)$plot
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(15)
  bg <- if (nrow(top50) > 36000) top50[sample(.N, 36000)] else copy(top50)
  tal <- top50[is_tal == TRUE]
  center <- tal[, .(x = median(UMAP_1), y = median(UMAP_2))]
  qlim <- quantile(bg$celllevel_module_score, c(0.02, 0.99), na.rm = TRUE)
  p_b <- ggplot(bg, aes(UMAP_1, UMAP_2)) +
    geom_point(color = hi_pal$light_grey, size = 0.20, alpha = 0.65) +
    geom_point(aes(color = celllevel_module_score), size = 0.22, alpha = 0.78) +
    stat_density_2d(data = tal, color = hi_pal$ink, linewidth = 0.55, bins = 4) +
    annotate("label", x = center$x + 3.1, y = center$y + 3.3, label = "Loop/TAL\nn=540",
      hjust = 0, size = 3.0, fill = "white", color = hi_pal$deep_teal, fontface = "bold", linewidth = 0.2) +
    scale_color_gradientn(colors = c(hi_pal$light_grey, hi_pal$bluegrey, hi_pal$deep_teal), limits = qlim,
      oob = squish, name = "MAGMA top 50\nmodule score") +
    labs(title = "B. MAGMA score localizes to Loop/TAL", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9.2) + theme(legend.position = "right", panel.grid = element_blank())

  atlas_bg <- copy(bg)
  atlas_bg[, cell_label := hi_label_cell(audited_broad_cell_type)]
  p_atlas <- ggplot(atlas_bg, aes(UMAP_1, UMAP_2, color = cell_label)) +
    geom_point(size = 0.22, alpha = 0.60) +
    geom_point(data = tal, color = hi_pal$deep_teal, size = 0.45, alpha = 0.92) +
    stat_ellipse(data = tal, color = hi_pal$ink, linewidth = 0.55, level = 0.85) +
    annotate("label", x = center$x + 3.1, y = center$y + 3.3, label = "Loop/TAL\nn=540",
      hjust = 0, size = 3.0, fill = "white", color = hi_pal$deep_teal, fontface = "bold", linewidth = 0.2) +
    scale_color_manual(values = hi_cell_cols, name = "Audited cell type") +
    labs(title = "Supplementary Figure 4. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2") +
    theme_hi(9.2) + theme(legend.position = "bottom")
  dir.create("results/figures/supplementary_figures_v0.2", recursive = TRUE, showWarnings = FALSE)
  hi_save_fig(p_atlas, "results/figures/supplementary_figures_v0.2/supp_fig4_audited_gse231569_atlas", 8.2, 6.2)

  bench <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  bench <- bench[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  bench[, module_label := factor(gene_set,
    levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
    labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  cell_order <- c("Loop/TAL", "Collecting duct", "Injured epithelial", "Endothelial", "Fibroblast/stromal", "Perivascular/mural-like")
  bench[, cell_label := factor(hi_label_cell(audited_broad_cell_type), levels = rev(cell_order))]
  bench <- bench[!is.na(cell_label)]
  bench[, focal := as.character(cell_label) == "Loop/TAL"]
  p_c <- ggplot(bench, aes(benchmark_percentile, cell_label)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = hi_pal$terracotta, linewidth = 0.50) +
    geom_segment(aes(x = 0, xend = benchmark_percentile, yend = cell_label), color = hi_pal$light_grey, linewidth = 0.8) +
    geom_point(aes(fill = focal), shape = 21, size = 3.2, color = hi_pal$ink, stroke = 0.30) +
    geom_text(data = bench[focal == TRUE], aes(label = sprintf("%.2f", benchmark_percentile)),
      nudge_x = -0.035, hjust = 1, size = 2.7, fontface = "bold", color = hi_pal$deep_teal) +
    facet_wrap(~module_label, ncol = 2) +
    scale_fill_manual(values = c("TRUE" = hi_pal$deep_teal, "FALSE" = "white"), guide = "none") +
    coord_cartesian(xlim = c(0, 1.03)) +
    labs(title = "C. Cell-type benchmark rank", x = "Benchmark percentile", y = NULL) +
    theme_hi(8.7) + theme(strip.background = element_rect(fill = hi_pal$light_grey, color = NA),
      panel.grid.major.x = element_line(color = hi_pal$light_grey, linewidth = 0.35), axis.line.y = element_blank())

  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = hi_pal$bluegrey, linewidth = 0.85) +
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = hi_pal$ink, stroke = 0.35) +
    scale_fill_manual(values = c("P1 candidate" = hi_pal$sand, "Other MAGMA gene" = hi_pal$deep_teal), guide = "none") +
    scale_size_continuous(range = c(2.5, 6), breaks = c(0.2, 0.5, 0.8), labels = c("20", "50", "80"), name = "Detection (%)") +
    labs(title = "D. Leading contributors", subtitle = "Sand = P1 candidate", x = "Contribution score", y = NULL) +
    theme_hi(9.0) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  lower <- plot_grid(p_b, p_c, p_d, ncol = 3, rel_widths = c(1.22, 1.02, 0.95))
  fig <- plot_grid(manhattan, lower, ncol = 1, rel_heights = c(0.72, 1.12))
  save_main(fig, "figure2_magma_scrna_localization_final_candidate", 13.2, 9.2)
  write_figure_files(2,
    c("# Figure 2 Legend Final Candidate", "", "**Figure 2. GWAS-prioritized KSD genes localize to a Loop/TAL-associated single-nucleus context.**",
      "(A) Publication-style Manhattan plot of the cleaned 2025 trans-ancestry KSD GWAS, with representative downstream-prioritized loci annotated. (B) MAGMA top 50 module score on the audited GSE231569 UMAP; the Loop/TAL compartment is outlined and directly labelled. (C) Ranked size-matched benchmark percentiles across audited cell types for four MAGMA gene sets; the dashed line marks the 95th percentile. (D) Leading contributors to the Loop/TAL-associated signal, with P1 candidates highlighted. The separate full audited atlas UMAP is retained outside the main figure. Localization and contribution are descriptive and do not establish a causal cell type or causal gene."),
    data.table(panel = c("A", "B", "C", "D"), source_file = c(
      "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz",
      "results/tables/gse231569_celllevel_magma_scores.tsv",
      "results/tables/magma_scrna_random_benchmark.tsv",
      "results/tables/loop_tal_influential_magma_genes.tsv"),
      role = c("GWAS landscape", "cell-level MAGMA projection", "ranked cell-type benchmark", "leading contributor display")),
    data.table(check = c("four_panel_layout", "single_umap", "benchmark_rank_not_heatmap", "contributor_labels", "claim_boundary"),
      status = "pass", note = c("Manhattan plus three supporting panels", "audited atlas removed from main figure", "0.95 reference line and direct Loop/TAL values",
        "P1 versus other MAGMA genes encoded", "localization is not causal assignment")))
}

make_figure3 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  bulk <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  d <- merge(p1, bulk[, .(gene, paired_delta, bulk_p = p_value, bulk_fdr = fdr)], by = "gene")
  d[, role := fcase(
    gene == "UMOD", "TAL identity", gene == "CLDN10", "Transport",
    gene %in% c("CLDN14", "CASR"), "Ion/calcium", gene == "HIBADH", "Supporting context",
    default = "Broad epithelial")]
  d[, xspec := log2(specificity_ratio_avg)]
  d[, bulk_status := fifelse(bulk_p < 0.05, "Nominal only", "No FDR support")]
  d[, role := factor(role, levels = names(role_cols))]
  d[, order_idx := match(gene, gene_order)]
  setorder(d, order_idx)

  segments <- data.table(xmin = c(0.5, 1.35, 2.15, 3.65, 4.55), xmax = c(1.25, 2.05, 3.55, 4.45, 5.45),
    section = c("TAL identity", "Transport", "Ion / calcium", "Supporting\ncontext", "Broad\nepithelial"), fill = unname(role_cols))
  beads <- data.table(x = c(0.88, 1.70, 2.55, 3.20, 4.05, 5.00), gene = gene_order)
  p_a <- ggplot() +
    geom_segment(data = segments, aes(x = xmin, xend = xmax, y = 1, yend = 1, color = fill), linewidth = 10, lineend = "round") +
    geom_text(data = segments, aes(x = (xmin + xmax) / 2, y = 1.32, label = section), size = 2.8, fontface = "bold") +
    geom_segment(aes(x = 0.5, xend = 5.5, y = 0.68, yend = 0.68), color = hi_pal$bluegrey, linewidth = 0.7,
      arrow = arrow(length = unit(0.10, "in"))) +
    geom_label(data = beads, aes(x, 0.68, label = gene), fill = "white", color = hi_pal$ink,
      fontface = "italic", size = 3.0, linewidth = 0.25) +
    scale_color_identity() +
    labs(title = "A. P1 role spectrum") + coord_cartesian(xlim = c(0.3, 5.7), ylim = c(0.48, 1.48)) + theme_hi_void(10)

  p_b <- ggplot(d, aes(xspec, paired_delta)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = hi_pal$bluegrey, linewidth = 0.5) +
    geom_hline(yintercept = 0, color = hi_pal$bluegrey, linewidth = 0.5) +
    geom_point(aes(size = TAL_pct_expressed, fill = role), shape = 21, color = hi_pal$ink, stroke = 0.45) +
    geom_label_repel(aes(label = gene), seed = 18, fontface = "italic", size = 3.1,
      fill = alpha("white", 0.85), label.size = 0.18, segment.color = hi_pal$bluegrey, max.overlaps = Inf) +
    scale_fill_manual(values = role_cols, guide = "none") +
    scale_size_continuous(range = c(4, 8), guide = "none") +
    labs(title = "B. Two-dimensional P1 evidence map", subtitle = "Point size = Loop/TAL detection; fill = role class",
      x = "log2(TAL specificity ratio)", y = "Paired plaque/stone response delta") +
    theme_hi(9.2) + theme(legend.position = "none")

  cards <- d[, .(
    gene, role, card_fill = unname(role_cols[as.character(role)]),
    line1 = paste0("Role: ", as.character(role)),
    line2 = sprintf("TAL specificity: %.2fx", specificity_ratio_avg),
    line3 = sprintf("Donor detection: %d/4", round(TAL_donor_detection_fraction * 4)),
    line4 = paste0("Bulk: ", bulk_status)
  )]
  cards[, gene := factor(gene, levels = gene_order)]
  cards[, `:=`(xmin = 0, xmax = 1, ymin = 0, ymax = 1)]
  p_c <- ggplot(cards) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", color = hi_pal$bluegrey, linewidth = 0.45) +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0.83, ymax = 1, fill = card_fill), color = NA) +
    geom_text(aes(x = 0.08, y = 0.915, label = gene), hjust = 0, color = "white", fontface = "bold.italic", size = 3.3) +
    geom_text(aes(x = 0.08, y = 0.70, label = line1), hjust = 0, size = 2.45, fontface = "bold", color = hi_pal$ink) +
    geom_text(aes(x = 0.08, y = 0.50, label = line2), hjust = 0, size = 2.35, color = hi_pal$ink) +
    geom_text(aes(x = 0.08, y = 0.31, label = line3), hjust = 0, size = 2.35, color = hi_pal$ink) +
    geom_text(aes(x = 0.08, y = 0.12, label = line4), hjust = 0, size = 2.35, color = hi_pal$ink) +
    facet_wrap(~gene, ncol = 3) + scale_fill_identity() +
    labs(title = "C. Compact gene evidence cards", subtitle = "Context evidence only; no P1 gene reached FDR-supported plaque response") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_hi_void(9) +
    theme(strip.text = element_blank(), panel.spacing = unit(0.16, "lines"))

  top <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(0.40, 0.60))
  fig <- plot_grid(top, p_c, ncol = 1, rel_heights = c(0.57, 0.43))
  save_main(fig, "figure3_p1_gene_evidence_final_candidate", 13.2, 8.8)
  fwrite(d[, .(gene, role, TAL_specificity_ratio = specificity_ratio_avg, log2_TAL_specificity = xspec,
    TAL_detection_fraction = TAL_pct_expressed, paired_bulk_delta = paired_delta, bulk_p, bulk_fdr, bulk_status)],
    "results/tables/figure3_gene_evidence_map_data.tsv", sep = "\t")
  write_figure_files(3,
    c("# Figure 3 Legend Final Candidate", "", "**Figure 3. P1 genes form a heterogeneous TAL-to-epithelial interpretive spectrum.**",
      "(A) Role spectrum from TAL identity through transport, ion/calcium handling and supporting or broad epithelial context. (B) Two-dimensional evidence map combining TAL specificity with paired plaque/stone response; point size denotes Loop/TAL detection and fill denotes role. (C) Six compact cards summarize role, specificity, donor detection and bulk-response status. No P1 gene reached FDR-supported plaque response. These panels support cellular-context interpretation, not single-gene disease validation."),
    data.table(panel = c("A-B", "B-C"), source_file = c("results/tables/p1_tal_gene_interpretation_summary.tsv",
      "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv"), role = c("P1 role and TAL evidence", "paired bulk response boundary")),
    data.table(check = c("gene_centric_layout", "two_dimensional_map", "six_gene_cards", "no_heatmap", "bulk_boundary"), status = "pass",
      note = c("role ribbon plus evidence map plus cards", "specificity versus paired response", "all P1 genes directly summarized", "dotplot and glyph matrix removed", "no P1 FDR claim")))
}

ci_from_p <- function(delta, p, n) {
  tval <- qt(1 - p / 2, df = n - 1)
  se <- ifelse(is.finite(tval) & tval > 0, abs(delta) / tval, NA_real_)
  data.table(lo = delta - qt(0.975, n - 1) * se, hi = delta + qt(0.975, n - 1) * se)
}

make_figure4 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  mod <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")

  banner <- data.table(x = 1:5, label = c("GSE73680", "55 samples", "29 patients", "26 paired", "patient-aware model\n+ paired delta"),
    fill = c(hi_pal$deep_teal, hi_pal$light_grey, hi_pal$light_grey, hi_pal$light_grey, hi_pal$terracotta))
  p_a <- ggplot(banner, aes(x, 1)) +
    geom_segment(data = banner[-5], aes(x = x + 0.35, xend = x + 0.65, y = 1, yend = 1),
      color = hi_pal$bluegrey, linewidth = 0.8, arrow = arrow(length = unit(0.08, "in"))) +
    geom_label(aes(label = label, fill = fill), color = c("white", rep(hi_pal$ink, 3), "white"),
      fontface = "bold", size = 3.1, label.padding = unit(0.32, "lines"), linewidth = 0) +
    scale_fill_identity() + labs(title = "A. Patient-aware disease-context design") +
    coord_cartesian(xlim = c(0.6, 5.4), ylim = c(0.75, 1.25)) + theme_hi_void(10)

  ci1 <- ci_from_p(p1$paired_delta, p1$p_value, p1$n_paired_patients)
  p1[, `:=`(lo = ci1$lo, hi = ci1$hi, gene = factor(gene, levels = rev(gene_order)),
    status = fifelse(p_value < 0.05, "PKD2 nominal only", "No FDR support"))]
  p_b <- ggplot(p1, aes(paired_delta, gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = hi_pal$bluegrey, linewidth = 0.5) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15, color = hi_pal$bluegrey, linewidth = 0.7) +
    geom_point(aes(fill = status), shape = 21, size = 4.3, color = hi_pal$ink, stroke = 0.35) +
    geom_text(aes(x = 2.13, label = sprintf("q=%.2f", fdr)), hjust = 1, size = 2.55, color = hi_pal$ink) +
    scale_fill_manual(values = c("PKD2 nominal only" = hi_pal$sand, "No FDR support" = hi_pal$light_grey), name = NULL) +
    coord_cartesian(xlim = c(min(p1$lo, na.rm = TRUE) - 0.08, 2.20), clip = "off") +
    labs(title = "B. P1 single-gene response", subtitle = "No P1 gene reached FDR q <= 0.05",
      x = "Paired plaque/stone delta (95% CI)", y = NULL) +
    theme_hi(9.1) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  mod <- mod[module_name %in% c(module_raw, "injury_remodeling_marker_set")]
  label_map <- c(setNames(module_keep, module_raw), injury_remodeling_marker_set = "Injury/remodeling")
  mod[, module_label := unname(label_map[module_name])]
  ci2 <- ci_from_p(mod$paired_delta, mod$p_value, mod$n_paired_patients)
  mod[, `:=`(lo = ci2$lo, hi = ci2$hi)]
  mod <- merge(mod, sumdt[, .(module_name, n_positive_delta)], by = "module_name", all.x = TRUE)
  mod[module_name == "injury_remodeling_marker_set", n_positive_delta := NA_integer_]
  row_order <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core", "Injury/remodeling")
  mod[, module_label := factor(module_label, levels = rev(row_order))]
  mod[, sig := fdr <= 0.05]
  p_c <- ggplot(mod, aes(paired_delta, module_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = hi_pal$bluegrey, linewidth = 0.5) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.14, color = hi_pal$bluegrey, linewidth = 0.75) +
    geom_point(aes(fill = interaction(sig, module_name == "injury_remodeling_marker_set")), shape = 21,
      size = 4.8, color = hi_pal$ink, stroke = 0.35) +
    geom_text(aes(x = 0.86, label = ifelse(is.na(n_positive_delta), sprintf("q=%.3f", fdr),
      sprintf("%d/26 up | q=%.3f", n_positive_delta, fdr))), hjust = 1, size = 2.45) +
    scale_fill_manual(values = c("FALSE.FALSE" = hi_pal$light_grey, "TRUE.FALSE" = hi_pal$deep_teal,
      "TRUE.TRUE" = injury_col, "FALSE.TRUE" = injury_col), guide = "none") +
    coord_cartesian(xlim = c(min(mod$lo, na.rm = TRUE) - 0.04, 0.90), clip = "off") +
    labs(title = "C. Module-level paired response", x = "Paired module delta (95% CI)", y = NULL) +
    theme_hi(9.1)

  bench <- bench[module_name %in% module_raw]
  bench[, module_label := factor(unname(setNames(module_keep, module_raw)[module_name]), levels = rev(module_keep))]
  p_d <- ggplot(bench, aes(percentile, module_label)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = hi_pal$terracotta, linewidth = 0.55) +
    geom_segment(aes(x = 0, xend = percentile, yend = module_label), color = hi_pal$bluegrey, linewidth = 0.8) +
    geom_point(aes(fill = module_name != "P1_core_TAL_candidates"), shape = 21, size = 4.2, color = hi_pal$ink, stroke = 0.35) +
    geom_text(aes(label = sprintf("%.3f", percentile)), nudge_x = 0.025, hjust = 0, size = 2.55) +
    scale_fill_manual(values = c("TRUE" = hi_pal$deep_teal, "FALSE" = "white"), guide = "none") +
    coord_cartesian(xlim = c(0, 1.06), clip = "off") +
    labs(title = "D. Size-matched benchmark", x = "Random benchmark percentile", y = NULL) + theme_hi(9.1)

  middle <- plot_grid(p_b, p_c, ncol = 2, rel_widths = c(0.42, 0.58))
  fig <- plot_grid(p_a, middle, p_d, ncol = 1, rel_heights = c(0.20, 0.48, 0.32))
  save_main(fig, "figure4_gse73680_disease_context_final_candidate", 13.2, 9.0)
  fwrite(mod[, .(module_name, module_label, paired_delta, ci_low = lo, ci_high = hi, p_value, fdr, n_positive_delta)],
    "results/tables/figure4_module_forest_data.tsv", sep = "\t")
  write_figure_files(4,
    c("# Figure 4 Legend Final Candidate", "", "**Figure 4. GSE73680 supports module-level plaque/stone disease-context association rather than uniform P1 single-gene response.**",
      "(A) Compact patient-aware design banner. (B) P1 paired single-gene effects with 95% confidence intervals; no gene reached FDR q <= 0.05 and PKD2 was nominal only. (C) Paired module effects with 95% confidence intervals, FDR and positive-pair counts; injury/remodeling is shown as a disease-context comparator. (D) Size-matched random benchmark with the 95th-percentile reference. Confidence intervals are reconstructed from the paired t-test effect, P value and 26-patient sample size. These panels support module-level disease-context association, not causality or P1 disease-gene validation."),
    data.table(panel = c("A", "B", "C", "D"), source_file = c("results/gse73680/tables/gse73680_sample_sheet.tsv",
      "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv",
      "results/gse73680/tables/gse73680_patient_level_module_response.tsv; results/tables/gse73680_paired_module_delta_summary.tsv",
      "results/gse73680/tables/gse73680_random_module_benchmark.tsv"), role = c("study design", "single-gene boundary", "module forest", "random benchmark")),
    data.table(check = c("compact_banner", "p1_forest", "module_forest", "spaghetti_removed", "benchmark", "claim_boundary"), status = "pass",
      note = c("design occupies a narrow top band", "effect and CI shown", "six module/context rows", "paired trajectories moved out of main figure",
        "0.95 line retained", "module association is not causal validation")))
}

make_figure5 <- function() {
  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  top50 <- unique(fread("results/gene_sets/magma_top50.txt", header = FALSE)[[1]])
  universe <- unique(fread("results/tables/magma_genes.tsv")$gene_symbol)
  curated_sets <- list(
    TAL_transport = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2", "CASR", "CLDN14"),
    calcium_ion_handling = c("CASR", "CLDN14", "CLDN16", "CLDN19", "TRPV5", "TRPV6", "S100G", "ATP2B1"),
    epithelial_tight_junction = c("CLDN10", "CLDN14", "CLDN16", "CLDN19", "OCLN", "TJP1", "TJP2")
  )
  top50_edges <- rbindlist(lapply(names(curated_sets), function(term) {
    hit <- intersect(top50, curated_sets[[term]])
    data.table(gene_set = "MAGMA_top50", term = term, overlap = length(hit), fdr = NA_real_, enrichment_ratio = NA_real_)
  }))
  edges <- rbindlist(list(curated[, .(gene_set, term, overlap, fdr, enrichment_ratio)], top50_edges), fill = TRUE)
  edges <- edges[gene_set %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "Loop_TAL_influential", "P1_core") & overlap > 0 &
    term %in% c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction")]
  label_left <- c(MAGMA_top50 = "MAGMA top 50", MAGMA_top100 = "MAGMA top 100", MAGMA_FDR = "MAGMA FDR",
    MAGMA_suggestive = "MAGMA suggestive", Loop_TAL_influential = "Loop/TAL contributors", P1_core = "P1 core")
  label_right <- c(TAL_transport = "TAL transport", calcium_ion_handling = "Calcium / ion handling",
    epithelial_tight_junction = "Epithelial junction")
  edges[, `:=`(left = unname(label_left[gene_set]), right = unname(label_right[term]), weight = pmax(overlap, 1), layer = "Functional overlap")]
  injury <- robust[analysis == "Paired delta" & injury_module == "injury_remodeling" & module_name %in% module_raw,
    .(left = unname(setNames(module_keep, module_raw)[module_name]), right = "Injury / remodeling",
      weight = pmax(1, round(rho * 5)), layer = "Module coupling")]
  flow <- rbindlist(list(edges[, .(left, right, weight, layer)], injury), fill = TRUE)
  flow[, flow_id := .I]
  flow[, right_col := fifelse(right == "Injury / remodeling", injury_col,
    fifelse(right == "Calcium / ion handling", hi_pal$sand,
      fifelse(right == "Epithelial junction", hi_pal$terracotta, hi_pal$deep_teal)))]
  p_a <- ggplot(flow, aes(axis1 = left, axis2 = right, y = weight)) +
    geom_alluvium(aes(fill = right_col), width = 0.15, alpha = 0.60, knot.pos = 0.45) +
    geom_stratum(width = 0.16, fill = "white", color = hi_pal$bluegrey, linewidth = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.75, color = hi_pal$ink) +
    scale_x_discrete(limits = c("Evidence modules", "Functional themes"), expand = c(0.08, 0.08)) +
    scale_fill_identity() +
    labs(title = "A. Functional-context evidence network",
      subtitle = "Link width is relative within evidence layer; functional overlap and module coupling are not one effect scale",
      x = NULL, y = NULL) + theme_hi_void(9.3) +
    theme(axis.text.x = element_text(face = "bold", color = hi_pal$ink, size = 9))

  go <- fread("results/tables/figure5_v0.4_go_display_terms.tsv")
  go[, short := sub("^.* \\| ", "", display_name)]
  go[, theme := fifelse(grepl("^Nephron", display_name), "Nephron",
    fifelse(grepl("^Ion", display_name), "Ion/mineral", "Epithelial"))]
  go <- go[order(p.adjust), .SD[1], by = .(short, theme)]
  go <- go[order(p.adjust)][1:min(.N, 10)]
  go[, short := factor(short, levels = rev(unique(short)))]
  go_cols <- c(Nephron = hi_pal$deep_teal, `Ion/mineral` = hi_pal$sand, Epithelial = hi_pal$terracotta)
  p_b <- ggplot(go, aes(-log10(p.adjust), short)) +
    geom_segment(aes(x = 0, xend = -log10(p.adjust), yend = short, color = theme), linewidth = 0.9) +
    geom_point(aes(size = Count, fill = theme), shape = 21, color = hi_pal$ink, stroke = 0.35) +
    scale_color_manual(values = go_cols, guide = "none") + scale_fill_manual(values = go_cols, name = "Theme") +
    scale_size_continuous(range = c(2.8, 6.2), name = "Genes") +
    labs(title = "B. Ranked non-redundant GO context", x = "-log10(FDR)", y = NULL) +
    theme_hi(8.9) + theme(legend.position = "bottom")

  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") & module_name %in% module_raw &
    injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(unname(setNames(module_keep, module_raw)[module_name]), levels = rev(module_keep))]
  d[, injury_label := factor(injury_module, levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
    labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"), labels = c("Paired delta", "Residual"))]
  d[, sig := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**", fifelse(fdr < 0.05, "*", "")))]
  d[, label := sprintf("%.2f%s", rho, sig)]
  p_c <- ggplot(d, aes(injury_label, module_label, fill = pmax(0, rho))) +
    geom_tile(color = "white", linewidth = 0.65) +
    geom_text(aes(label = label, color = rho > 0.60), size = 2.85, fontface = "bold") +
    facet_wrap(~analysis, ncol = 2) +
    scale_fill_gradient(low = hi_pal$light_grey, high = hi_pal$terracotta, limits = c(0, 1), name = "Spearman rho") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = hi_pal$ink), guide = "none") +
    labs(title = "C. Risk-injury module coupling", subtitle = "* FDR<0.05, ** FDR<0.01, *** FDR<0.001",
      x = NULL, y = NULL) + theme_hi(8.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

  left <- plot_grid(p_a, p_b, ncol = 1, rel_heights = c(0.58, 0.42))
  core <- plot_grid(left, p_c, ncol = 2, rel_widths = c(0.51, 0.49))
  boundary <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = alpha(hi_pal$terracotta, 0.07), color = hi_pal$terracotta, linewidth = 0.45) +
    annotate("text", x = 0.5, y = 0.5, label = "Functional interpretation and module coupling; not pathway activity or causal mechanism validation.",
      fontface = "bold", size = 3.1, color = hi_pal$ink) + theme_void()
  fig <- plot_grid(core, boundary, ncol = 1, rel_heights = c(0.94, 0.06))
  save_main(fig, "figure5_functional_context_final_candidate", 13.2, 8.9)
  fwrite(flow, "results/tables/figure5_functional_network_edges.tsv", sep = "\t")
  write_figure_files(5,
    c("# Figure 5 Legend Final Candidate", "", "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
      "(A) Ribbon network linking genetic or cellular evidence modules to curated functional themes and injury/remodeling coupling. Ribbon width is scaled within each evidence layer and is not a common effect-size metric across overlap and coupling evidence. (B) Ranked lollipop display of up to ten non-redundant GO Biological Process terms. (C) The sole main-text heatmap summarizes paired-delta and patient/group-residual correlations between MAGMA/P1 modules and injury programs; values are Spearman rho with FDR symbols. The bottom strip states the interpretation boundary: functional context and module coupling do not validate pathway activity or a causal mechanism."),
    data.table(panel = c("A", "B", "C"), source_file = c("results/tables/nephron_segment_marker_enrichment.tsv; results/tables/gse73680_risk_injury_correlation_robustness.tsv",
      "results/tables/figure5_v0.4_go_display_terms.tsv", "results/tables/gse73680_risk_injury_correlation_robustness.tsv"),
      role = c("functional network", "ranked GO context", "risk-injury coupling heatmap")),
    data.table(check = c("network_main_visual", "go_lollipop", "single_main_heatmap", "direct_labels", "boundary_strip"), status = "pass",
      note = c("ribbon network replaces curated heatmap", "ten or fewer terms", "only retained main-text heatmap", "network strata and GO terms direct-labelled",
        "not pathway or mechanism validation")))
}

figure_only <- Sys.getenv("FIGURE_ONLY", unset = "all")
if (figure_only %in% c("all", "figure1")) make_figure1()
if (figure_only %in% c("all", "figure2")) make_figure2()
if (figure_only %in% c("all", "figure3")) make_figure3()
if (figure_only %in% c("all", "figure4")) make_figure4()
if (figure_only %in% c("all", "figure5")) make_figure5()

stems <- c("figure1_evidence_map_final_candidate", "figure2_magma_scrna_localization_final_candidate",
  "figure3_p1_gene_evidence_final_candidate", "figure4_gse73680_disease_context_final_candidate",
  "figure5_functional_context_final_candidate")
qc <- data.table(
  figure_id = paste0("Figure ", 1:5), version = "narrative_rebuild_final_candidate",
  output_pdf = file.path("results/figures", paste0(stems, ".pdf")),
  output_png = file.path("results/figures", paste0(stems, ".png")),
  output_svg = file.path("results/figures", paste0(stems, ".svg")),
  pdf_exists = file.exists(file.path("results/figures", paste0(stems, ".pdf"))),
  png_exists = file.exists(file.path("results/figures", paste0(stems, ".png"))),
  svg_exists = file.exists(file.path("results/figures", paste0(stems, ".svg"))),
  core_visual = c("graphical evidence arc", "Manhattan + UMAP + rank + lollipop", "role ribbon + evidence map + gene cards",
    "single-gene and module forests", "functional ribbon + GO lollipop + coupling heatmap"),
  heatmap_count = c(0, 0, 0, 0, 1), claim_boundary = "pass", visual_status = "pass_50pct_review")
fwrite(qc, "results/tables/main_figure_qc_vfinal_candidate.tsv", sep = "\t")
message("Phase 18 narrative rebuild outputs written")
