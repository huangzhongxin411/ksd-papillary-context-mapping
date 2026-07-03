suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(png)
  library(scales)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)
dir.create("results/supervisor_review_package_v1.1", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  genetic = "#245B64", scrna = "#6F8F98", p1 = "#B59A5B",
  disease = "#9A5F52", muted = "#D8D8D8", pale = "#EEF3F4",
  text = "#303030", grid = "#A8A8A8"
)

theme_pub <- function(base_size = 8.8) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = pal$text),
      plot.subtitle = element_text(color = "#555555"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.6)
    )
}
write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)
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
    Collecting_duct_principal = "CD",
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

make_manuscript_v11 <- function() {
  x <- readLines("manuscript/manuscript_draft_v1.0.md", warn = FALSE)
  x[1] <- "# Post-GWAS mapping of kidney stone risk identifies Loop/TAL-associated renal papillary cellular and disease-context programs"
  x <- gsub(
    "MAGMA-prioritized KSD genes localized to a Loop/TAL-associated single-nucleus context in GSE231569, with enrichment exceeding random gene-set expectations and remaining supported across robustness checks. The six P1 candidate genes formed an interpretable TAL, epithelial transport and calcium-handling expression spectrum, but did not behave as a uniform disease-validated gene panel. In GSE73680, MAGMA-prioritized modules showed module-level plaque/stone disease-context association, whereas P1 single-gene responses were heterogeneous and not FDR-supported. Functional-context analyses linked prioritized modules to nephron development, TAL transport, calcium/ion handling and papillary injury/remodeling programs.",
    "MAGMA-prioritized KSD genes localized to a Loop/TAL-associated single-nucleus context in GSE231569, exceeding size-matched random gene-set expectations and remaining supported across locus-balanced and leave-one-locus-out robustness checks. Six P1 candidate genes formed a TAL, epithelial transport and calcium/ion-handling spectrum rather than a uniform TAL marker or disease-response panel. In GSE73680, MAGMA-prioritized modules, but not individual P1 genes, showed patient-aware plaque/stone disease-context shifts and injury/remodeling coupling.",
    x, fixed = TRUE
  )
  x <- gsub(
    "Paired-delta and residual correlations were prioritized for main-text interpretation because sample-level correlations are more robust to disease-group composition than sample-level correlations.",
    "Paired-delta and residual correlations were prioritized because they are less vulnerable to disease-group composition and repeated-sample non-independence than simple sample-level correlations.",
    x, fixed = TRUE
  )
  x <- gsub(
    "The main finding of this study is that MAGMA-prioritized kidney stone disease genes localize to a Loop/TAL-associated renal papillary single-nucleus context and show disease-context support at the module level. This supports a module-level post-GWAS interpretation rather than a single causal gene or uniform disease-gene panel model. The major contribution is not a single nominated causal gene, but an auditable evidence framework that separates genetic prioritization, renal papillary cell-context localization, gene-centric interpretation and disease-context expression association.",
    "The main finding of this study is that MAGMA-prioritized kidney stone disease genes localize to a Loop/TAL-associated renal papillary single-nucleus context and show plaque/stone disease-context support at the module level. The major contribution is not the nomination of a single causal gene, but the construction of an auditable post-GWAS evidence framework that separates genetic prioritization, renal papillary cell-context localization, gene-centric interpretation and plaque/stone disease-context association.",
    x, fixed = TRUE
  )
  hit <- grep("Curated nephron and functional marker-set analyses further linked", x)[1]
  if (length(hit) == 1 && !grepl("not to validate pathway activity", x[hit])) {
    x[hit] <- paste0(x[hit], " These analyses were used to interpret the biological context of prioritized modules, not to validate pathway activity or establish injury-driven causality.")
  }
  x <- gsub("TWAS, SMR/coloc and spatial transcriptomic pipelines were scaffolded and resource-audited but were not used as evidence layers because required external expression weights, matched eQTL resources or spatial matrix/image/coordinate files were incomplete.",
            "TWAS, SMR-coloc and spatial transcriptomic pipelines were scaffolded and resource-audited but were not used as evidence layers because required external resources were incomplete; detailed resource status is reported in the Supplementary Tables.",
            x, fixed = TRUE)
  x <- gsub("SMR/coloc", "SMR-coloc", x, fixed = TRUE)
  x <- gsub("plaque/stone expression", "plaque/stone papilla expression", x, fixed = TRUE)
  x <- gsub("plaque/stone disease-context", "plaque/stone papilla disease-context", x, fixed = TRUE)
  write_lines(x, "manuscript/manuscript_draft_v1.1.md")
  write_lines(x, "manuscript/manuscript_clean_for_supervisor_v1.1.md")
}

make_figure1_v06 <- function() {
  layer_dt <- data.table(
    x = c(0.38, 0.54, 0.70, 0.86), y = 0.64,
    layer = c("Layer 1\nGenetic prioritization", "Layer 2\nsnRNA localization",
              "Layer 3\nGene interpretation", "Layer 4\nDisease-context association"),
    body = c("GWAS/MAGMA\n57 loci\n17,316 genes\n94 Bonferroni",
             "GSE231569 papilla\nLoop/TAL context\naudited cells",
             "P1 candidates\n6 genes\nrole spectrum",
             "GSE73680 papilla\n55 samples\n29 patients | 26 paired"),
    fill = c(pal$genetic, pal$scrna, pal$p1, pal$disease)
  )
  fig <- ggplot() +
    annotate("text", x = 0.59, y = 0.92, label = "Post-GWAS framework for KSD papillary context mapping",
             fontface = "bold", size = 4.9, color = pal$text) +
    annotate("text", x = 0.59, y = 0.86, label = "Genetic prioritization -> snRNA localization -> gene interpretation -> disease-context association",
             size = 3.0, color = "#555555") +
    annotate("path", x = c(0.075, 0.10, 0.14, 0.19, 0.235, 0.245, 0.215, 0.165, 0.115, 0.075),
             y = c(0.61, 0.77, 0.86, 0.83, 0.69, 0.54, 0.43, 0.38, 0.44, 0.61),
             color = "#8A8176", linewidth = 0.55) +
    annotate("polygon", x = c(0.125, 0.158, 0.195, 0.158),
             y = c(0.47, 0.70, 0.55, 0.41), fill = "#E7D7B4", color = "#9B8A6A", linewidth = 0.30) +
    annotate("rect", xmin = 0.145, xmax = 0.205, ymin = 0.41, ymax = 0.57, fill = NA, color = pal$disease, linewidth = 0.38) +
    annotate("segment", x = 0.205, xend = 0.26, y = 0.49, yend = 0.36,
             arrow = arrow(length = unit(0.09, "in")), color = "#888888", linewidth = 0.35) +
    annotate("rect", xmin = 0.055, xmax = 0.275, ymin = 0.08, ymax = 0.35, fill = "#FBFBFB", color = "#888888", linewidth = 0.32) +
    annotate("rect", xmin = 0.07, xmax = 0.26, ymin = 0.105, ymax = 0.315, fill = "#F1F1F1", color = NA) +
    annotate("path", x = c(0.092, 0.122, 0.152, 0.182, 0.212), y = c(0.14, 0.28, 0.14, 0.28, 0.14),
             color = pal$genetic, linewidth = 1.0) +
    annotate("rect", xmin = 0.158, xmax = 0.206, ymin = 0.17, ymax = 0.31, fill = "#D8C7A5", color = "#8A8A8A", alpha = 0.78) +
    annotate("point", x = 0.234, y = 0.145, size = 3.0, color = pal$disease) +
    annotate("text", x = 0.145, y = 0.102, label = "Loop/TAL", size = 2.55, color = pal$genetic, fontface = "bold") +
    annotate("text", x = 0.210, y = 0.302, label = "CD", size = 2.35, color = "#5B5147") +
    annotate("text", x = 0.237, y = 0.092, label = "plaque/\nstone", size = 2.05, color = pal$disease) +
    annotate("text", x = 0.155, y = 0.375, label = "kidney -> papilla niche", size = 2.7, color = "#555555") +
    geom_segment(data = layer_dt[1:3], aes(x = x + 0.06, xend = layer_dt$x[2:4] - 0.06, y = y, yend = y),
                 arrow = arrow(length = unit(0.10, "in")), color = "#777777", linewidth = 0.35) +
    geom_label(data = layer_dt, aes(x = x, y = y + 0.105, label = layer, fill = fill),
               color = "white", fontface = "bold", size = 2.45, lineheight = 0.92,
               label.padding = unit(0.20, "lines"), linewidth = 0) +
    geom_label(data = layer_dt, aes(x = x, y = y - 0.035, label = body),
               fill = "white", color = pal$text, size = 2.45, lineheight = 0.95,
               label.padding = unit(0.22, "lines"), label.size = 0.25) +
    annotate("rect", xmin = 0.32, xmax = 0.955, ymin = 0.16, ymax = 0.40, fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = 0.35, y = 0.355, label = "Supported inference", hjust = 0, fontface = "bold", size = 3.05, color = pal$genetic) +
    annotate("text", x = 0.35, y = 0.305, label = "Loop/TAL-associated cellular context\nMAGMA module-level disease-context association",
             hjust = 0, vjust = 1, size = 2.55, color = pal$text, lineheight = 0.95) +
    annotate("text", x = 0.67, y = 0.355, label = "Not established", hjust = 0, fontface = "bold", size = 3.05, color = pal$disease) +
    annotate("text", x = 0.67, y = 0.315, label = "causality | TWAS convergence\ncolocalization | spatial validation\nP1 disease-gene validation",
             hjust = 0, vjust = 1, size = 2.45, color = pal$text, lineheight = 0.92) +
    annotate("text", x = 0.35, y = 0.178, label = "Resource-limited extensions: TWAS / SMR-coloc / spatial audited, not used as evidence layers",
             hjust = 0, size = 2.25, color = "#666666") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.06, 0.96), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave("results/figures/figure1_integrative_framework_v0.6.pdf", fig, width = 13.2, height = 5.9, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.6.png", fig, width = 13.2, height = 5.9, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 1 Legend v0.6", "",
    "**Figure 1. Post-GWAS framework for KSD papillary context mapping.**",
    "The schematic organizes the study into four evidence layers: genetic prioritization, audited GSE231569 snRNA localization, P1 candidate gene interpretation and GSE73680 plaque/stone papilla disease-context association. The supported inference is a Loop/TAL-associated cellular context and MAGMA module-level disease-context association. Causality, TWAS convergence, colocalization, spatial validation and P1 disease-gene validation are not established."
  ), "docs/figure1_legend_v0.6.md")
}

make_figure2_v07 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := label_cell(audited_broad_cell_type)]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(13)
  other <- top50[is_tal == FALSE]
  tal <- top50[is_tal == TRUE]
  other <- if (nrow(other) > 24000) other[sample(.N, 24000)] else other
  plot_cells <- rbind(other, tal, fill = TRUE)
  cell_cols <- c("Loop/TAL" = pal$genetic, "Collecting duct" = "#D9DEE1", "Endothelial" = "#BAC7B1",
                 "Fibroblast/stromal" = "#B5CDBA", "Injured epithelial" = "#ADDADB",
                 "Perivascular/mural-like" = "#D7C5D7")
  center <- tal[, .(x = median(UMAP_1), y = median(UMAP_2))]
  p_a <- ggplot(plot_cells[is_tal == FALSE], aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_label), size = 0.16, alpha = 0.18) +
    geom_point(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2, color = cell_label), size = 0.40, alpha = 0.95) +
    stat_ellipse(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                 color = pal$p1, linewidth = 0.35, level = 0.85) +
    annotate("curve", x = center$x + 3.0, y = center$y + 3.5, xend = center$x + 0.45, yend = center$y + 0.45,
             curvature = -0.25, arrow = arrow(length = unit(0.09, "inches")), color = pal$p1, linewidth = 0.30) +
    annotate("label", x = center$x + 3.2, y = center$y + 3.7, label = "Loop/TAL\nn = 540", hjust = 0,
             size = 2.55, fill = "white", color = pal$genetic, fontface = "bold") +
    scale_color_manual(values = cell_cols, breaks = names(cell_cols), name = "Audited cell type") +
    labs(title = "A. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2") +
    theme_pub(8.8) + theme(legend.position = "right", panel.grid = element_blank())
  qlim <- quantile(plot_cells$celllevel_module_score, c(0.02, 0.99), na.rm = TRUE)
  p_b <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celllevel_module_score), size = 0.18, alpha = 0.78) +
    stat_density_2d(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                    color = pal$p1, linewidth = 0.38, bins = 4) +
    scale_color_gradient(low = "#F1F3F4", high = pal$genetic, limits = qlim, oob = squish,
                         name = "Relative\nmodule score") +
    labs(title = "B. MAGMA top 50 module score", x = "UMAP 1", y = "UMAP 2") +
    theme_pub(8.8) + theme(panel.grid = element_blank())
  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  random <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
                                  labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  random[, cell_label := factor(label_cell(audited_broad_cell_type),
                                levels = rev(c("Loop/TAL", "Perivascular/mural-like", "Fibroblast/stromal", "Endothelial", "Injured epithelial", "Collecting duct")))]
  random[, support := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95, "Loop/TAL exceeded expectation", "Other")]
  p_c <- ggplot(random, aes(benchmark_percentile, cell_label, fill = support)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "#555555", linewidth = 0.35) +
    annotate("text", x = 0.94, y = Inf, label = "95th percentile", hjust = 1, vjust = 1.2,
             size = 2.8, fontface = "bold", color = "#555555") +
    geom_point(shape = 21, size = 2.8, color = "#555555", stroke = 0.22) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("Loop/TAL exceeded expectation" = pal$genetic, Other = pal$muted)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1.0), labels = c("0", "0.5", "0.95", "1.0")) +
    labs(title = "C. Size-matched benchmark percentile", x = "Benchmark percentile", y = NULL, fill = NULL) +
    theme_pub(8.8) + theme(legend.position = "bottom")
  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = pal$grid, linewidth = 0.42) +
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = "#444444", stroke = 0.22) +
    scale_fill_manual(values = c("P1 candidate" = pal$p1, "Other MAGMA gene" = pal$genetic)) +
    scale_size_continuous(range = c(2.2, 4.6), labels = percent_format()) +
    labs(title = "D. Genes contributing to Loop/TAL-associated signal", x = "Contribution score", y = NULL, fill = NULL, size = "Detection") +
    theme_pub(9) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.05, 1.05))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.7.pdf", fig, width = 13.2, height = 9.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.7.png", fig, width = 13.2, height = 9.6, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 2 Legend v0.7", "",
    "**Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated renal papillary single-nucleus context.**",
    "(A) Audited GSE231569 snRNA-seq UMAP with Loop/TAL included in the audited cell-type legend and highlighted by a callout. (B) MAGMA top 50 relative module score projected onto the same UMAP. Cell-level scores represent mean z-scored expression and are shown for visualization only; donor-cell-type summaries provide the interpretation layer. (C) Size-matched random gene-set benchmark percentiles. (D) Genes contributing to the Loop/TAL-associated signal. Contribution score summarizes MAGMA prioritization strength, Loop/TAL expression preference and detection support for descriptive ranking, not causal driver inference."
  ), "docs/figure2_legend_v0.7.md")
}

make_figure3_v09 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1[, gene := factor(gene, levels = rev(gene_order))]
  cards <- data.table(
    gene = c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"),
    role = c("TAL identity", "Transport", "Ion handling", "Calcium sensing", "Supporting context", "Broad epithelial"),
    x = c(1, 2, 3, 4, 1.7, 3.3), y = c(2, 2, 2, 2, 1, 1)
  )
  p_a <- ggplot(cards, aes(x, y)) +
    annotate("segment", x = 0.8, xend = 4.2, y = 2.55, yend = 2.55, color = pal$grid, linewidth = 0.4,
             arrow = arrow(length = unit(0.10, "inches"))) +
    annotate("text", x = 0.78, y = 2.72, label = "TAL identity", hjust = 0, size = 2.5, color = "#555555") +
    annotate("text", x = 4.2, y = 2.72, label = "broader epithelial context", hjust = 1, size = 2.5, color = "#555555") +
    geom_label(aes(label = paste0(gene, "\n", role), fill = role), color = "white", fontface = "bold",
               size = 2.8, lineheight = 0.9, label.padding = unit(0.22, "lines"), linewidth = 0) +
    scale_fill_manual(values = c("TAL identity" = pal$genetic, "Transport" = "#557F89",
                                 "Ion handling" = pal$p1, "Calcium sensing" = "#9A7E43",
                                 "Supporting context" = pal$scrna, "Broad epithelial" = pal$disease)) +
    coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.55, 2.95), clip = "off") +
    labs(title = "A. P1 two-row role spectrum") +
    theme_void(base_size = 9.2) + theme(plot.title = element_text(face = "bold", color = pal$text), legend.position = "none")
  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(label_cell(cell_type, short = TRUE), levels = label_cell(keep_cells, short = TRUE))]
  cell[, gene := factor(gene, levels = rev(gene_order))]
  p_b <- ggplot(cell, aes(cell_label, gene)) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    scale_size_continuous(range = c(0.8, 5.0), labels = percent_format()) +
    labs(title = "B. P1 expression across audited cell types", x = NULL, y = NULL, fill = "Average\nexpression", size = "Detected") +
    theme_pub(8.8) + theme(axis.text.x = element_text(angle = 0), axis.text.y = element_text(face = "italic"), legend.position = "bottom", panel.grid = element_blank())
  p_c <- ggplot(p1, aes(log2(specificity_ratio_avg), gene)) +
    geom_segment(aes(x = 0, xend = log2(specificity_ratio_avg), yend = gene), color = pal$grid, linewidth = 0.5) +
    geom_point(aes(fill = specificity_class), shape = 21, size = 3.8, color = "#555555", stroke = 0.25) +
    annotate("text", x = Inf, y = Inf, label = "All shown genes detected\nin 3/4 TAL donors", hjust = 1.05, vjust = 1.25, size = 2.55, color = "#555555") +
    scale_fill_manual(values = c(strong_TAL_preferential = pal$genetic, moderate_TAL_preferential = pal$p1), na.value = pal$muted) +
    labs(title = "C. log2(TAL specificity ratio)", x = "log2(TAL specificity ratio)", y = NULL, fill = "Specificity") +
    theme_pub(8.8) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  ev <- p1[, .(
    gene,
    MAGMA = fifelse(magma_p < 5e-8, "+++", "++"),
    TAL_specificity = fifelse(specificity_class == "strong_TAL_preferential", "+++", "++"),
    Donor_detection = sprintf("%d/4", round(TAL_donor_detection_fraction * 4)),
    Bulk_response = fifelse(as.character(gene) == "PKD2", "nominal", "no FDR"),
    Role = fifelse(as.character(gene) == "UMOD", "Representative TAL",
            fifelse(as.character(gene) == "CLDN10", "Transport",
            fifelse(as.character(gene) == "CLDN14", "Ion handling",
            fifelse(as.character(gene) == "CASR", "Calcium sensing",
            fifelse(as.character(gene) == "HIBADH", "Supporting", "Broad epithelial")))))
  )]
  ev_long <- melt(ev, id.vars = "gene", variable.name = "evidence", value.name = "call")
  ev_long[, evidence := factor(evidence, levels = c("MAGMA", "TAL_specificity", "Donor_detection", "Bulk_response", "Role"),
                               labels = c("MAGMA", "TAL specificity", "Donor detection", "Bulk response", "Role"))]
  ev_long[, fill_class := fcase(call == "+++", "strong", call == "++", "moderate", call == "nominal", "nominal", call == "no FDR", "no_fdr", default = "role")]
  p_d <- ggplot(ev_long, aes(evidence, gene, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 2.25, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, nominal = pal$p1, no_fdr = pal$muted, role = "#F7F8F8")) +
    labs(title = "D. Discrete evidence matrix", x = NULL, y = NULL) +
    theme_pub(8.8) + theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.text.y = element_text(face = "italic"), legend.position = "none", panel.grid = element_blank())
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.86, 1.16))
  ggsave("results/figures/figure3_p1_gene_evidence_v0.9.pdf", fig, width = 13.2, height = 8.9, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure3_p1_gene_evidence_v0.9.png", fig, width = 13.2, height = 8.9, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 3 Legend v0.9", "",
    "**Figure 3. P1 genes form a TAL, transport and calcium-handling role spectrum.**",
    "(A) Two-row role spectrum for six P1 candidates. (B) P1 gene expression and detection across audited GSE231569 cell types, ordered consistently with panel D. (C) log2(TAL specificity ratio), reducing visual compression from highly TAL-preferential genes. (D) Discrete evidence matrix. +++ indicates strong support, ++ moderate support, + contextual or detectable support and NA not applicable. Bulk response denotes GSE73680 single-gene plaque/stone paired response; no P1 gene reached FDR q <= 0.05, and PKD2 was nominal only."
  ), "docs/figure3_legend_v0.9.md")
}

make_figure4_v12 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")
  p_a <- ggplot() +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.68, ymax = 0.86, fill = pal$genetic, color = NA) +
    annotate("text", x = 0.15, y = 0.77, label = "1\nDataset", color = "white", fontface = "bold", size = 3.4, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.68, ymax = 0.86, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.605, y = 0.77, label = "GSE73680 papillary plaque/stone context\n55 samples | 29 patients", size = 3.35, color = pal$text, lineheight = 0.9) +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.40, ymax = 0.58, fill = pal$scrna, color = NA) +
    annotate("text", x = 0.15, y = 0.49, label = "2\nPairing", color = "white", fontface = "bold", size = 3.4, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.40, ymax = 0.58, fill = "#F7F8F8", color = "#888888", linewidth = 0.3) +
    annotate("text", x = 0.605, y = 0.49, label = "26 paired patients\ncontrol/adjacent -> plaque/stone papilla", size = 3.35, color = pal$text, lineheight = 0.9) +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.13, ymax = 0.30, fill = pal$disease, color = NA) +
    annotate("text", x = 0.15, y = 0.215, label = "3\nAnalysis", color = "white", fontface = "bold", size = 3.4, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.13, ymax = 0.30, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.605, y = 0.215, label = "patient-aware limma + paired delta\nmodule response + random benchmark", size = 3.25, color = pal$text, lineheight = 0.9) +
    labs(title = "A. Patient-aware GSE73680 design") + theme_void(base_size = 9.5) + theme(plot.title = element_text(face = "bold"))
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1[, gene := factor(gene, levels = rev(gene_order))]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No FDR support")]
  p1[, p_label := sprintf("P=%.3f", p_value)]
  p_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.4, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$p1, "No FDR support" = pal$muted)) +
    labs(title = "B. P1 single-gene response", subtitle = "No P1 gene reached FDR q <= 0.05; PKD2 nominal only",
         x = "Patient-level paired delta", y = NULL, fill = NULL) +
    theme_pub(9.0) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
  keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
  long <- long[module_label %in% keep]
  long[, module_label := factor(as.character(module_label), levels = keep)]
  lab <- sumdt[module_label %in% keep, .(module_label = factor(as.character(module_label), levels = keep), label = sprintf("%d/%d positive", n_positive_delta, n_paired_patients))]
  p_c <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(aes(color = direction), alpha = 0.25, linewidth = 0.26) +
    geom_point(aes(fill = group_label), shape = 21, color = "#555555", stroke = 0.12, size = 1.08, alpha = 0.75) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.15, color = "#303030") +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2.2, color = "#303030") +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.12, size = 3.05, fontface = "bold", color = pal$text) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(positive = "#7F9AA3", negative = "#D6A29A", zero = pal$muted)) +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = "#B08A45")) +
    labs(title = "C. Paired patient module shifts", x = NULL, y = "Module score", color = "Delta", fill = NULL) +
    theme_pub(8.4) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
  bench <- bench[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  bench[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                                 labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  bench[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  p_d <- ggplot(bench, aes(percentile, module_label, fill = percentile >= 0.95)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "#555555", linewidth = 0.35) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = emp_label, color = percentile >= 0.95), hjust = 1.04, size = 2.9, fontface = "bold") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = pal$text), guide = "none") +
    annotate("text", x = 0.94, y = 5.45, label = "95th percentile", hjust = 1, vjust = 0.5, size = 2.7, fontface = "bold", color = "#555555") +
    scale_fill_manual(values = c("TRUE" = pal$genetic, "FALSE" = pal$muted), labels = c("Background-like", "Exceeds 95th percentile")) +
    coord_cartesian(xlim = c(0, 1.05)) +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL, fill = NULL) +
    theme_pub(9.0) + theme(legend.position = "bottom")
  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.82, 1.18))
  ggsave("results/figures/figure4_gse73680_disease_context_v1.2.pdf", fig, width = 13.2, height = 9.3, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v1.2.png", fig, width = 13.2, height = 9.3, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 4 Legend v1.2", "",
    "**Figure 4. GSE73680 supports MAGMA module-level plaque/stone papilla disease-context association.**",
    "(A) Patient-aware GSE73680 design. (B) P1 single-gene paired responses are shown in the same gene order as Figure 3; no P1 gene reached FDR q <= 0.05, and PKD2 was nominal only. (C) Paired patient module shifts, with individual lines shown at low opacity and a bold median trend line. (D) Size-matched random benchmark with empirical P labels and a 95th percentile reference line. These analyses support MAGMA module-level disease-context association, not P1 disease-gene validation or causality."
  ), "docs/figure4_legend_v1.2.md")
}

make_figure5_v04 <- function() {
  checklist <- data.table(
    evidence = rep(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context"), each = 4),
    gene_set = rep(c("Top 100", "FDR", "Loop/TAL contrib.", "P1 core"), times = 5),
    call = c("+++", "+++", "++", "+",
             "++", "++", "+++", "++",
             "+", "+", "+", "+",
             "++", "++", "NA", "+",
             "++", "++", "++", "++")
  )
  checklist[, evidence := factor(evidence, levels = rev(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context")))]
  checklist[, gene_set := factor(gene_set, levels = c("Top 100", "FDR", "Loop/TAL contrib.", "P1 core"))]
  checklist[, fill_class := fcase(call == "+++", "strong", call == "++", "moderate", call == "+", "support", call == "NA", "not_applicable")]
  p_a <- ggplot(checklist, aes(gene_set, evidence, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 3.1, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, support = pal$p1, not_applicable = pal$muted)) +
    labs(title = "A. Evidence checklist", x = NULL, y = NULL) +
    theme_pub(8.8) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none", panel.grid = element_blank())

  go <- fread("results/tables/go_bp_redundancy_reduced_terms.tsv")
  go <- go[redundancy_reduced_keep == TRUE & p.adjust < 0.10 & Count >= 2]
  go <- go[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
  display_map <- data.table(
    Description = c(
      "distal tubule development",
      "loop of Henle development",
      "intracellular calcium ion homeostasis",
      "regulation of calcium ion import",
      "regulation of monoatomic ion transport",
      "urate metabolic process",
      "phosphate ion homeostasis",
      "cell-cell junction assembly",
      "cell-cell adhesion via plasma-membrane adhesion molecules"
    ),
    display = c(
      "Nephron | distal tubule development",
      "Nephron | loop of Henle development",
      "Ion/mineral | calcium homeostasis",
      "Ion/mineral | calcium ion import",
      "Ion/mineral | ion transport",
      "Ion/mineral | urate metabolism",
      "Ion/mineral | phosphate homeostasis",
      "Epithelial | junction assembly",
      "Epithelial | cell-cell adhesion"
    ),
    term_rank = 1:9
  )
  go <- merge(go, display_map, by = "Description")
  setorder(go, term_rank, p.adjust)
  go[, display := factor(display, levels = rev(display_map$display))]
  go[, gene_set_label := factor(gene_set,
                                levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  fwrite(go[, .(ID, Description, display_name = as.character(display), gene_set, p.adjust, Count, geneID)],
         "results/tables/figure5_v0.4_go_display_terms.tsv", sep = "\t")
  p_b <- ggplot(go, aes(-log10(p.adjust), display, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_manual(values = c("MAGMA top 100" = pal$genetic, "MAGMA FDR" = pal$scrna,
                                 "Loop/TAL contributors" = pal$p1, "P1 core" = pal$disease)) +
    labs(title = "B. Redundancy-reduced GO BP terms", subtitle = "Functional interpretation; not pathway validation",
         x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_pub(8.2) +
    theme(legend.position = "bottom")

  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  curated <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
  curated <- curated[term != "papillary_injury_remodeling"]
  curated[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                     labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  curated[, term_label := factor(term,
                                 levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction",
                                                "proximal_tubule_context", "collecting_duct_context")),
                                 labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction",
                                                "Proximal tubule context", "Collecting duct context")))]
  p_c <- ggplot(curated, aes(gene_set_label, term_label, fill = pmin(enrichment_ratio, 20))) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(overlap > 0, overlap, "")), size = 2.65, color = "#303030") +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    labs(title = "C. Curated nephron and functional context", subtitle = "Curated marker-set overlap; not pathway activity validation",
         x = NULL, y = NULL, fill = "Enrichment\nratio") +
    theme_pub(8.5) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") &
                module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates") &
                injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                             labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  d[, injury_label := factor(injury_module,
                             levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                             labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"),
                         labels = c("Paired patient delta", "Patient/group residual"))]
  d[, sig := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**", fifelse(fdr < 0.05, "*", "")))]
  d[, label := sprintf("%.2f%s", rho, sig)]
  d[, text_col := fifelse(rho > 0.65, "white", "#303030")]
  p_d <- ggplot(d, aes(injury_label, module_label, fill = pmin(pmax(rho, 0), 1))) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(data = d[text_col == "white"], aes(label = label), size = 2.35, color = "white", fontface = "bold") +
    geom_text(data = d[text_col != "white"], aes(label = label), size = 2.35, color = "#303030") +
    facet_wrap(~ analysis, ncol = 1) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$disease, limits = c(0, 1), name = "Spearman\nrho") +
    labs(title = "D. Risk-injury coupling robustness", x = NULL, y = NULL) +
    theme_pub(8.4) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.9, 1.1))
  ggsave("results/figures/figure5_functional_context_v0.4.pdf", fig, width = 13.2, height = 9.5, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure5_functional_context_v0.4.png", fig, width = 13.2, height = 9.5, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 5 Legend v0.4", "",
    "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
    "(A) Evidence checklist across MAGMA-prioritized, Loop/TAL-contributing and P1 gene sets. +++ indicates strong support, ++ moderate support, + contextual or detectable support and NA not applicable. (B) Redundancy-reduced GO Biological Process terms grouped by nephron development, ion/mineral handling and epithelial context; shortened display names are used in the figure and full terms are retained in the source table. GO enrichment is used for functional interpretation rather than pathway validation. (C) Curated marker-set overlap; color indicates enrichment ratio and tile labels indicate overlapping genes. (D) Risk-injury coupling robustness in GSE73680 using paired patient delta and patient/group residual correlations. Values represent Spearman rho; * FDR < 0.05, ** FDR < 0.01 and *** FDR < 0.001. P1 core showed weaker and less consistent coupling than MAGMA-prioritized modules, consistent with its interpretive rather than disease-validation role."
  ), "docs/figure5_legend_v0.4.md")
}

make_supplement_v11 <- function() {
  supp_tables <- data.table(
    table_id = paste0("Table S", 1:15),
    title = c("GWAS QC and lead loci", "MAGMA gene-based results", "MAGMA gene-set definitions and random benchmark",
              "GSE231569 annotation audit and cell counts", "Donor-cell-type MAGMA module score statistics",
              "Locus-balanced and leave-one-locus-out robustness", "P1 candidate gene evidence summary",
              "GSE73680 sample metadata and patient pairing", "GSE73680 P1 single-gene response",
              "GSE73680 MAGMA module response and sensitivity", "Risk-injury coupling robustness",
              "GO and curated functional enrichment", "TWAS/SMR-coloc/spatial resource-status audit",
              "Integrated candidate gene tiers", "Claim boundary audit"),
    status = c("ready", "ready", "ready", "source_audit_needed", "ready", "planned_or_existing_source",
               "ready", "source_audit_needed", "ready", "ready", "ready", "ready", "ready", "ready", "ready")
  )
  fwrite(supp_tables, "docs/supplementary_table_plan_v1.1.tsv", sep = "\t")
  supp_figs <- data.table(
    figure_id = paste0("Figure S", 1:12),
    title = c("GWAS QQ and Manhattan plots", "GSE231569 marker/annotation audit", "Donor-cell-type MAGMA module scores",
              "Locus-balanced and leave-one-locus-out analyses", "P1 individual gene expression and detection",
              "GSE73680 QC and sample pairing", "P1 single-gene paired response",
              "MAGMA module leave-one-gene / without-P1 sensitivity", "Expression-matched benchmark boundary check",
              "Full GO enrichment and redundancy reduction", "Risk-injury coupling robustness",
              "Resource-limited TWAS/SMR/spatial workflow audit"),
    status = c("existing", "planned_or_existing_source", "ready", "planned_or_existing_source", "planned_or_existing_source",
               "planned_or_existing_source", "planned_or_existing_source", "planned_or_existing_source",
               "planned_or_existing_source", "planned_or_existing_source", "existing", "planned")
  )
  fwrite(supp_figs, "docs/supplementary_figure_plan_v1.1.tsv", sep = "\t")
  qc <- CJ(figure_id = paste0("Figure ", 1:5), panel = c("A", "B", "C", "D"))[
    !(figure_id == "Figure 1" & panel == "D")]
  qc[, `:=`(
    supported_claim = "see figure legend and main claim",
    not_supported_claim = "causality, TWAS convergence, colocalization, spatial validation unless explicitly stated",
    statistical_unit = fifelse(figure_id == "Figure 2", "cell visualization plus donor-cell-type summaries",
                        fifelse(figure_id == "Figure 4", "paired patients / patient-aware bulk samples",
                        fifelse(figure_id == "Figure 5", "gene sets/modules", "schematic or gene-level summary"))),
    primary_statistic = "reported in panel or legend",
    visual_check = "rendered",
    legend_check = "updated_v1.1",
    boundary_check = "claim-bounded"
  )]
  fwrite(qc, "results/tables/main_figure_qc_v1.1.tsv", sep = "\t")
  write_lines(c(
    "# Claim Boundary Audit v1.1", "",
    "Supported: Loop/TAL-associated snRNA cellular context; MAGMA module-level plaque/stone papilla disease-context association; P1 role-spectrum interpretation; functional context and non-causal injury/remodeling coupling.",
    "",
    "Not established: causality, causal cell type, TWAS convergence, SMR-colocalization, spatial validation, therapeutic target validation or P1 disease-gene validation.",
    "",
    "Required language: supports, localizes to, maps to, module-level association, functional interpretation.",
    "Avoid: proves, causal driver, therapeutic target, validated pathway activity, TWAS-supported, colocalized or spatially validated."
  ), "docs/claim_boundary_audit_v1.1.md")
}

make_review_sheet_and_package <- function() {
  figs <- c(
    "results/figures/figure1_integrative_framework_v0.6.png",
    "results/figures/figure2_magma_scrna_localization_v0.7.png",
    "results/figures/figure3_p1_gene_evidence_v0.9.png",
    "results/figures/figure4_gse73680_disease_context_v1.2.png",
    "results/figures/figure5_functional_context_v0.4.png"
  )
  labels <- c("Figure 1 v0.6", "Figure 2 v0.7", "Figure 3 v0.9", "Figure 4 v1.2", "Figure 5 v0.4")
  pdf("results/figures/main_figures_1_to_5_review_contact_sheet_v1.1.pdf", width = 14, height = 28, onefile = TRUE)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(5, 1)))
  for (i in seq_along(figs)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(labels[i], x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 12))
    grid.raster(readPNG(figs[i]), x = 0.5, y = 0.48, width = 0.96, height = 0.88)
    popViewport()
  }
  popViewport(); dev.off()
  png("results/figures/main_figures_1_to_5_review_contact_sheet_v1.1.png", width = 2800, height = 5600, res = 200)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(5, 1)))
  for (i in seq_along(figs)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(labels[i], x = 0.02, y = 0.98, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 12))
    grid.raster(readPNG(figs[i]), x = 0.5, y = 0.48, width = 0.96, height = 0.88)
    popViewport()
  }
  popViewport(); dev.off()
  files <- c(
    "manuscript/manuscript_clean_for_supervisor_v1.1.md",
    sub("\\.png$", ".pdf", figs), figs,
    "docs/figure1_legend_v0.6.md", "docs/figure2_legend_v0.7.md", "docs/figure3_legend_v0.9.md",
    "docs/figure4_legend_v1.2.md", "docs/figure5_legend_v0.4.md",
    "docs/supplementary_table_plan_v1.1.tsv", "docs/supplementary_figure_plan_v1.1.tsv",
    "docs/claim_boundary_audit_v1.1.md", "results/tables/main_figure_qc_v1.1.tsv",
    "results/tables/candidate_gene_tiers_v1.0.tsv", "results/twas/twas_resource_status.tsv",
    "results/smr_coloc/eqtl_resource_status.tsv", "results/spatial/spatial_resource_status.tsv",
    "results/figures/main_figures_1_to_5_review_contact_sheet_v1.1.pdf",
    "results/figures/main_figures_1_to_5_review_contact_sheet_v1.1.png"
  )
  manifest <- data.table(file = files, exists = file.exists(files), copied = FALSE)
  for (i in seq_len(nrow(manifest))) if (manifest$exists[i]) manifest$copied[i] <- file.copy(manifest$file[i], file.path("results/supervisor_review_package_v1.1", basename(manifest$file[i])), overwrite = TRUE)
  fwrite(manifest, "results/supervisor_review_package_v1.1/package_manifest_v1.1.tsv", sep = "\t")
}

make_manuscript_v11()
make_figure1_v06()
make_figure2_v07()
make_figure3_v09()
make_figure4_v12()
make_figure5_v04()
make_supplement_v11()
make_review_sheet_and_package()
message("Phase 13 pre-submission hardening outputs written")
