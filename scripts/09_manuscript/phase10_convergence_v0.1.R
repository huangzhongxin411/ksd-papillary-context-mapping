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

make_figure2_v05 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := label_cell(audited_broad_cell_type)]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(3)
  plot_cells <- if (nrow(top50) > 24000) top50[sample(.N, 24000)] else top50
  cell_cols <- c(
    "Collecting duct" = "#D5DCE0",
    "Endothelial" = "#9EAF86",
    "Fibroblast/stromal" = "#86B391",
    "Injured epithelial" = "#78C8C9",
    "Loop/TAL" = "#1F5B65",
    "Perivascular/mural-like" = "#C9A6C9"
  )
  p_a <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_label, alpha = is_tal), size = 0.18) +
    scale_color_manual(values = cell_cols) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.34), guide = "none") +
    annotate("text", x = Inf, y = -Inf, label = "Loop/TAL, n=540", hjust = 1.05, vjust = -0.7,
             size = 2.8, fontface = "bold", color = pal$genetic) +
    labs(title = "A. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right", panel.grid = element_blank())

  qlim <- quantile(plot_cells$celllevel_module_score, c(0.01, 0.99), na.rm = TRUE)
  tal_cells <- plot_cells[is_tal == TRUE]
  p_b <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celllevel_module_score), size = 0.18, alpha = 0.78) +
    stat_density_2d(data = tal_cells, aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                    color = "#B59A5B", linewidth = 0.45, bins = 4) +
    scale_color_gradient(low = "#F1F3F4", high = "#1F5B65", limits = qlim, oob = scales::squish) +
    labs(title = "B. MAGMA top 50 module score on UMAP", x = "UMAP 1", y = "UMAP 2", color = "Score") +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())

  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  random <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
                                  labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  y_order <- c("Loop/TAL", "Perivascular/mural-like", "Fibroblast/stromal", "Endothelial", "Injured epithelial", "Collecting duct")
  random[, cell_label := factor(label_cell(audited_broad_cell_type), levels = rev(y_order))]
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
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = "#444444", stroke = 0.22) +
    scale_fill_manual(values = c("P1 candidate" = pal$p1, "Other MAGMA gene" = pal$genetic)) +
    scale_size_continuous(range = c(2.2, 4.6), labels = scales::percent_format()) +
    labs(title = "D. Genes contributing to Loop/TAL-associated signal", x = "Contribution score", y = NULL, fill = NULL, size = "Detection") +
    theme_bw(base_size = 9.0) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.05, 1.05))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.5.pdf", fig, width = 13.2, height = 9.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.5.png", fig, width = 13.2, height = 9.6, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 2 Legend v0.5",
    "",
    "**Figure 2. MAGMA-prioritized KSD genes project to a Loop/TAL-associated renal papillary single-nucleus context.**",
    "(A) Audited GSE231569 snRNA-seq UMAP with low-saturation cell-type colors and highlighted Loop/TAL cells. (B) Per-cell MAGMA top 50 module score projected onto the same UMAP, with Loop/TAL contour shown for orientation. (C) Size-matched random gene-set benchmark percentiles, with the dashed line marking the 95th percentile. (D) Genes contributing to the Loop/TAL-associated signal, ranked by a contribution score combining MAGMA prioritization strength, Loop/TAL specificity and donor detection. Cell-level UMAP panels are visualization; donor-cell-type summaries provide the analysis unit for context mapping. These analyses do not establish causal cell-type mediation."
  ), "docs/figure2_legend_v0.5.md")
}

make_figure3_v07 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1[, gene := factor(gene, levels = gene_order)]
  axis_dt <- data.table(
    gene = factor(gene_order, levels = gene_order),
    x = seq_along(gene_order),
    role = c("TAL identity", "Transport", "Calcium/ion handling", "Calcium sensing", "Supporting context", "Broad epithelial context")
  )
  role_cols <- c("TAL identity" = pal$genetic, Transport = "#557F89",
                 "Calcium/ion handling" = pal$p1, "Calcium sensing" = "#9A7E43",
                 "Supporting context" = pal$scrna, "Broad epithelial context" = pal$disease)
  p_a <- ggplot(axis_dt, aes(x, 1)) +
    annotate("segment", x = 1, xend = 6, y = 1, yend = 1, color = pal$grid, linewidth = 0.45) +
    geom_label(aes(label = paste0(as.character(gene), "\n", role), fill = role),
               color = "white", fontface = "bold", size = 2.55, lineheight = 0.9,
               label.padding = unit(0.18, "lines"), linewidth = 0) +
    scale_fill_manual(values = role_cols) +
    coord_cartesian(xlim = c(0.55, 6.45), ylim = c(0.65, 1.35), clip = "off") +
    labs(title = "A. P1 functional role axis") +
    theme_void(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(label_cell(cell_type), levels = label_cell(keep_cells))]
  cell[, gene := factor(gene, levels = gene_order)]
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

  ev <- p1[, .(
    gene,
    MAGMA = fifelse(magma_p < 5e-8, "+++", "++"),
    TAL_specificity = fifelse(specificity_class == "strong_TAL_preferential", "+++", "++"),
    Donor_detection = sprintf("%d/4", round(TAL_donor_detection_fraction * 4)),
    GSE73680 = fifelse(as.character(gene) == "PKD2", "nominal", "no FDR"),
    Role = gsub("_", " ", manuscript_role)
  )]
  ev_long <- melt(ev, id.vars = "gene", variable.name = "evidence", value.name = "call")
  ev_long[, gene := factor(gene, levels = rev(gene_order))]
  ev_long[, evidence := factor(evidence, levels = c("MAGMA", "TAL_specificity", "Donor_detection", "GSE73680", "Role"),
                               labels = c("MAGMA", "TAL specificity", "Donor detection", "GSE73680", "Role"))]
  ev_long[, fill_class := fcase(
    call == "+++", "strong",
    call == "++", "moderate",
    call == "nominal", "nominal",
    call == "no FDR", "no_fdr",
    default = "role"
  )]
  p_d <- ggplot(ev_long, aes(evidence, gene, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 2.35, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, nominal = pal$p1,
                                 no_fdr = pal$muted, role = "#F7F8F8")) +
    labs(title = "D. Discrete evidence matrix", x = NULL, y = NULL) +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          axis.text.y = element_text(face = "italic"), legend.position = "none", panel.grid = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.78, 1.22))
  ggsave("results/figures/figure3_p1_gene_evidence_v0.7.pdf", fig, width = 13.2, height = 8.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure3_p1_gene_evidence_v0.7.png", fig, width = 13.2, height = 8.8, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 3 Legend v0.7",
    "",
    "**Figure 3. P1 candidate genes form a TAL, transport and calcium-handling role spectrum.**",
    "(A) Functional role axis from TAL identity to broad epithelial context. (B) P1 gene expression and detection across audited GSE231569 cell types. (C) TAL specificity and donor support. (D) Discrete evidence matrix summarizing MAGMA support, TAL specificity, donor detection, GSE73680 single-gene response and assigned manuscript role. Discrete evidence calls are used to avoid implying that heterogeneous evidence dimensions are directly additive."
  ), "docs/figure3_legend_v0.7.md")
}

make_figure5_v02 <- function() {
  checklist <- data.table(
    evidence = rep(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context"), each = 4),
    gene_set = rep(c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"), times = 5),
    call = c("+++", "+++", "++", "+",
             "++", "++", "+++", "++",
             "+", "+", "+", "+",
             "++", "++", "NA", "+",
             "++", "++", "++", "++")
  )
  checklist[, evidence := factor(evidence, levels = rev(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context")))]
  checklist[, gene_set := factor(gene_set, levels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  checklist[, fill_class := fcase(call == "+++", "strong", call == "++", "moderate", call == "+", "support", call == "NA", "not_applicable")]
  p_a <- ggplot(checklist, aes(gene_set, evidence, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 3.1, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, support = pal$p1, not_applicable = pal$muted)) +
    labs(title = "A. Evidence checklist", x = NULL, y = NULL) +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none", panel.grid = element_blank())

  go <- fread("results/tables/go_bp_redundancy_reduced_terms.tsv")
  go <- go[redundancy_reduced_keep == TRUE & p.adjust < 0.10 & Count >= 2]
  go <- go[grepl("loop|Henle|tubule|water|vitamin D|metal ion|phosphate|ion|transport|epithelial|urate", Description, ignore.case = TRUE)]
  go <- go[order(p.adjust)][1:min(.N, 10)]
  go[, Description_short := ifelse(nchar(Description) > 50, paste0(substr(Description, 1, 47), "..."), Description)]
  go[, Description_short := make.unique(Description_short)]
  go[, Description_short := factor(Description_short, levels = rev(unique(Description_short)))]
  go[, gene_set_label := factor(gene_set,
                                levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential"),
                                labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors"))]
  p_b <- ggplot(go, aes(-log10(p.adjust), Description_short, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_manual(values = c("MAGMA top 100" = pal$genetic, "MAGMA FDR" = pal$scrna, "Loop/TAL contributors" = pal$p1)) +
    labs(title = "B. Redundancy-reduced GO BP terms", subtitle = "Functional interpretation; not pathway validation",
         x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_bw(base_size = 8.6) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, color = "#555555"),
          legend.position = "bottom", panel.grid.minor = element_blank())

  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  curated <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
  curated[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                     labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  curated[, term_label := factor(term,
                                 levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction",
                                                "proximal_tubule_context", "collecting_duct_context", "papillary_injury_remodeling")),
                                 labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction",
                                                "Proximal tubule context", "Collecting duct context", "Papillary injury/remodeling")))]
  p_c <- ggplot(curated, aes(gene_set_label, term_label, fill = pmin(enrichment_ratio, 20))) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(overlap > 0, overlap, "")), size = 2.65, color = "#303030") +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    labs(title = "C. Curated nephron and functional context", subtitle = "Color = enrichment ratio; number = overlapping genes",
         x = NULL, y = NULL, fill = "Enrichment\nratio") +
    theme_bw(base_size = 8.6) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, color = "#555555"),
          axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") &
                module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates") &
                injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                             labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  d[, injury_label := factor(injury_module,
                             levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                             labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"))]
  p_d <- ggplot(d, aes(injury_label, module_label, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.35) +
    facet_wrap(~ analysis, ncol = 1) +
    scale_fill_gradient2(low = "#8AA0A8", mid = "white", high = pal$disease, midpoint = 0, limits = c(-1, 1)) +
    labs(title = "D. Risk-injury coupling robustness", x = NULL, y = NULL, fill = "rho") +
    theme_bw(base_size = 8.4) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.9, 1.1))
  ggsave("results/figures/figure5_functional_context_v0.2.pdf", fig, width = 13.2, height = 9.5, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure5_functional_context_v0.2.png", fig, width = 13.2, height = 9.5, units = "in", dpi = 260, bg = "white")
  write_lines(c(
    "# Figure 5 Legend v0.2",
    "",
    "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
    "(A) Evidence checklist across MAGMA-prioritized, Loop/TAL-contributing and P1 gene sets. (B) Redundancy-reduced GO Biological Process terms filtered for FDR, gene count and theme relevance. GO enrichment is used for functional interpretation rather than pathway-level validation. (C) Curated nephron and functional marker-set enrichment; color indicates enrichment ratio and tile labels indicate overlapping gene counts. (D) Risk-injury coupling robustness in GSE73680 using paired patient delta and patient/group residual correlations. These analyses support module-level disease-context coupling, not causal injury mechanisms."
  ), "docs/figure5_legend_v0.2.md")
}

write_v09 <- function() {
  x <- readLines("manuscript/manuscript_draft_v0.8.md", warn = FALSE)
  x <- gsub("MAGMA-prioritized KSD genes converged on a Loop/TAL-associated single-nucleus context in GSE231569, with enrichment exceeding random gene-set expectations and remaining supported across robustness checks. The six P1 candidate genes formed an interpretable TAL, epithelial transport, and calcium-handling expression spectrum, but did not behave as a uniform disease-validated gene panel. In GSE73680, MAGMA-prioritized modules showed module-level disease-context expression association, whereas P1 single-gene responses were heterogeneous and not FDR-supported.",
            "MAGMA-prioritized KSD genes converged on a Loop/TAL-associated single-nucleus context in GSE231569, with enrichment exceeding random gene-set expectations and remaining supported across robustness checks. The six P1 candidate genes formed an interpretable TAL, epithelial transport, and calcium-handling expression spectrum, but did not behave as a uniform disease-validated gene panel. In GSE73680, MAGMA-prioritized modules showed module-level disease-context expression association, whereas P1 single-gene responses were heterogeneous and not FDR-supported. Functional-context analyses linked prioritized modules to nephron development, TAL transport, calcium/ion handling and papillary injury/remodeling programs.",
            x, fixed = TRUE)

  old_methods <- grep("^### Cell-level MAGMA score projection and functional-context analyses", x)
  if (length(old_methods) == 1) {
    next_methods <- grep("^### ", x)
    end <- next_methods[next_methods > old_methods][1] - 1
    replacement <- c(
      "### Cell-level MAGMA module score projection",
      "",
      "An audited GSE231569 Seurat object with UMAP coordinates was used for cell-level visualization. MAGMA top 50, top 100, FDR and suggestive gene sets were scored at the cell level as mean z-scored expression across detected module genes. Donor-cell-type summaries were then computed from per-cell scores and used as the analysis unit for Loop/TAL module-score support. Cell-level UMAP score maps were treated as visualization rather than independent-cell statistical evidence.",
      "",
      "### Functional enrichment analysis",
      "",
      "GO Biological Process enrichment was performed with clusterProfiler and org.Hs.eg.db using MAGMA-tested genes mapped to Entrez IDs as the background universe. Displayed GO terms were filtered by FDR, gene count and redundancy. Curated nephron and functional marker-set enrichment was analyzed by hypergeometric testing against the MAGMA-tested gene universe. These analyses were used for functional interpretation rather than pathway-level validation.",
      "",
      "### GSE73680 injury/remodeling coupling analysis",
      "",
      "GSE73680 injury/remodeling, inflammation/immune, fibrosis/ECM and epithelial-injury programs were scored as mean z-scored marker-set modules. Coupling between MAGMA modules and injury programs was evaluated using sample-level Spearman correlations, paired patient-level delta correlations and patient/group residual correlations. Paired-delta and residual correlations were prioritized for main-text interpretation because sample-level correlations are more vulnerable to disease-group composition effects."
    )
    x <- c(x[seq_len(old_methods - 1)], replacement, x[(end + 1):length(x)])
  }

  disc_insert <- grep("^GSE73680 provided disease-context support", x)[1]
  add_disc <- "The functional-context analyses further help interpret why the Loop/TAL-associated signal is biologically plausible. GO and curated marker-set enrichment linked prioritized genes to nephron development, TAL transport, calcium/ion handling and epithelial junction contexts. In GSE73680, MAGMA modules showed positive paired-delta and patient/group residual coupling with papillary injury/remodeling programs. These findings connect genetic prioritization to renal epithelial transport and papillary injury contexts, but remain functional annotations and disease-context associations rather than causal mechanism tests."
  if (!any(grepl("The functional-context analyses further help interpret", x, fixed = TRUE)) && length(disc_insert) == 1) {
    x <- append(x, c("", add_disc), after = disc_insert + 1)
  }

  lim_hit <- grep("^TWAS, SMR/coloc and spatial transcriptomic validation", x)[1]
  lim_add <- "Functional enrichment and risk-injury coupling analyses are gene-set- and module-level interpretations and may be influenced by gene-set size, background definition, expression detectability and disease-group composition."
  if (!any(grepl("Functional enrichment and risk-injury coupling analyses", x, fixed = TRUE)) && length(lim_hit) == 1) {
    x[lim_hit] <- paste(x[lim_hit], lim_add)
  }
  x <- x[!grepl("Phase [0-9]+", x)]
  write_lines(x, "manuscript/manuscript_draft_v0.9.md")
}

make_figure2_v05()
make_figure3_v07()
make_figure5_v02()
write_v09()

qc <- data.table(
  item = c("Figure 2 v0.5", "Figure 3 v0.7", "Figure 5 v0.2", "manuscript v0.9"),
  status = "completed",
  main_change = c("UMAP palette and donor-level framing", "functional role axis and discrete evidence matrix", "checklist, filtered GO and robust coupling", "abstract/methods/results/discussion/limitations convergence")
)
fwrite(qc, "results/tables/phase10_convergence_qc_v0.1.tsv", sep = "\t")
message("wrote Phase 10 convergence outputs")
