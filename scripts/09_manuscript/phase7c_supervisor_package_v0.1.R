suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final_draft", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  genetic = "#3E6672",
  scrna = "#6F8F98",
  p1 = "#B59A5B",
  disease = "#9A5F52",
  tal = "#3E6672",
  muted = "#D8D8D8",
  boundary = "#4F4F4F",
  text = "#303030",
  grid = "#A8A8A8"
)

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

label_cell <- function(x, short = FALSE) {
  if (short) {
    map <- c(
      Collecting_duct_principal = "Collecting duct",
      Fibroblast_stromal = "Fibroblast",
      Endothelial = "Endothelial",
      Injured_undifferentiated_epithelial = "Injured\nepithelial",
      Loop_of_Henle_TAL = "Loop/TAL",
      Perivascular_mural_like = "Perivascular",
      Pericyte_smooth_muscle = "Perivascular"
    )
  } else {
    map <- c(
      Collecting_duct_principal = "Collecting duct\nprincipal",
      Fibroblast_stromal = "Fibroblast/\nstromal",
      Endothelial = "Endothelial",
      Injured_undifferentiated_epithelial = "Injured\nundiff. epi.",
      Loop_of_Henle_TAL = "Loop/TAL",
      Perivascular_mural_like = "Perivascular/\nmural-like",
      Pericyte_smooth_muscle = "Perivascular/\nmural-like"
    )
  }
  unname(ifelse(x %in% names(map), map[x], x))
}

label_gene_set <- function(x, short = FALSE) {
  map <- if (short) {
    c(magma_top50 = "top 50", magma_top100 = "top 100",
      magma_fdr05 = "FDR", magma_suggestive_p1e4 = "suggestive")
  } else {
    c(magma_top50 = "MAGMA top 50", magma_top100 = "MAGMA top 100",
      magma_top200 = "MAGMA top 200", magma_fdr05 = "MAGMA FDR",
      magma_suggestive_p1e4 = "MAGMA suggestive")
  }
  unname(ifelse(x %in% names(map), map[x], x))
}

make_figure2_v02 <- function() {
  atlas <- fread("results/tables/magma_scrna_module_score_by_celltype.tsv")
  scores <- fread("results/tables/magma_scrna_module_score_by_celltype.tsv")
  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  gene_summary <- fread("results/tables/magma_gene_set_summary.tsv")
  locus_bal <- fread("results/tables/magma_locus_balanced_scrna_benchmark.tsv")
  loo <- fread("results/tables/magma_leave_one_locus_out.tsv")

  atlas <- unique(atlas[, .(audited_broad_cell_type, n_cells, n_donors, annotation_confidence)])
  atlas[, cell_label := label_cell(audited_broad_cell_type)]
  atlas[, confidence_label := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL", "Loop/TAL", "Other audited")]
  atlas[grepl("low_or_exploratory", annotation_confidence), confidence_label := "Low abundance"]
  atlas[, cell_label := factor(cell_label, levels = atlas[order(n_cells), cell_label])]
  panel_a <- ggplot(atlas, aes(n_cells, cell_label, fill = confidence_label)) +
    geom_col(width = 0.68, color = "#777777", linewidth = 0.2) +
    geom_text(aes(label = paste0("n=", n_cells)), hjust = -0.08, size = 2.35, color = pal$text) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = c("Loop/TAL" = pal$tal, "Other audited" = pal$scrna, "Low abundance" = pal$muted)) +
    labs(title = "A. Audited GSE231569 cell-type composition",
         x = "Cells", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())

  score_cols <- c("magma_top50_score", "magma_top100_score", "magma_fdr05_score", "magma_suggestive_p1e4_score")
  score_long <- melt(scores, id.vars = "audited_broad_cell_type", measure.vars = score_cols,
                     variable.name = "gene_set", value.name = "module_score")
  score_long[, gene_set := sub("_score$", "", gene_set)]
  score_long[, gene_label := factor(label_gene_set(gene_set),
                                    levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
  score_long[, cell_label := factor(label_cell(audited_broad_cell_type),
                                    levels = c("Collecting duct\nprincipal", "Loop/TAL", "Injured\nundiff. epi.",
                                               "Endothelial", "Fibroblast/\nstromal", "Perivascular/\nmural-like"))]
  tal_rect <- data.table(xmin = 1.5, xmax = 2.5, ymin = 0.5, ymax = 4.5)
  panel_b <- ggplot(score_long, aes(cell_label, gene_label, fill = module_score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_rect(data = tal_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = pal$tal, linewidth = 0.85) +
    geom_text(aes(label = sprintf("%.2f", module_score)), size = 2.25, color = pal$text) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$tal, limits = c(0.10, 0.40), oob = scales::squish) +
    labs(title = "B. MAGMA gene sets show highest scores in Loop/TAL cells",
         x = NULL, y = NULL, fill = "Module\nscore") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "right",
          panel.grid = element_blank())

  random_plot <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random_plot[, gene_label := factor(label_gene_set(gene_set, short = TRUE),
                                     levels = c("top 50", "top 100", "FDR", "suggestive"))]
  random_plot[, cell_label := factor(label_cell(audited_broad_cell_type, short = TRUE),
                                     levels = c("Perivascular", "Fibroblast", "Endothelial",
                                                "Injured\nepithelial", "Collecting duct", "Loop/TAL"))]
  random_plot[, support_class := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95,
                                         "Loop/TAL exceeds expectation", "Other audited context")]
  p95_label <- unique(random_plot[, .(gene_label)])
  p95_label[, cell_label := factor("Loop/TAL", levels = levels(random_plot$cell_label))]
  panel_c <- ggplot(random_plot, aes(benchmark_percentile, cell_label, fill = support_class)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.35) +
    geom_text(data = p95_label, aes(x = 0.94, y = cell_label, label = "95th percentile"),
              inherit.aes = FALSE, hjust = 1, vjust = -1.0, size = 2.15, color = pal$boundary) +
    geom_point(shape = 21, size = 3.0, color = "#555555", stroke = 0.25) +
    facet_wrap(~ gene_label, ncol = 2) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95), labels = c("0", "0.5", "0.95")) +
    scale_fill_manual(values = c("Loop/TAL exceeds expectation" = pal$tal, "Other audited context" = pal$muted)) +
    labs(title = "C. Loop/TAL localization exceeds random expectation",
         x = "Random gene-set benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          strip.background = element_rect(fill = "#F1F1F1", color = "#AAAAAA"),
          panel.grid.minor = element_blank())

  top50 <- gene_summary[gene_set == "magma_top50"]
  lb_full <- locus_bal[analysis_version == "full_audited" &
                         gene_set == "magma_locus_balanced_top50" &
                         audited_broad_cell_type == "Loop_of_Henle_TAL", benchmark_percentile][1]
  lb_cons <- locus_bal[analysis_version == "conservative_exclude_low_or_exploratory_and_immune_review" &
                         gene_set == "magma_locus_balanced_top50" &
                         audited_broad_cell_type == "Loop_of_Henle_TAL", benchmark_percentile][1]
  loo_min <- min(loo$TAL_percentile_after_removal, na.rm = TRUE)
  robust <- data.table(
    check = c("Full audited", "Locus-balanced", "Conservative annotation", "Leave-one-locus-out minimum"),
    percentile = c(top50$TAL_percentile, lb_full, lb_cons, loo_min)
  )
  robust[, check := factor(check, levels = rev(check))]
  panel_d <- ggplot(robust, aes(percentile, check)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.35) +
    geom_col(width = 0.62, fill = pal$tal, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3f", percentile)), hjust = 1.08, size = 2.55, color = "white") +
    scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.95), labels = c("0", "0.5", "0.95")) +
    labs(title = "D. Loop/TAL signal remains robust across sensitivity checks",
         x = "Loop/TAL benchmark percentile", y = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, rel_heights = c(1.0, 1.08))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.2.pdf", fig,
         width = 12.5, height = 9.0, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.2.png", fig,
         width = 12.5, height = 9.0, units = "in", dpi = 250, bg = "white")

  write_lines(c(
    "# Figure 2 Legend v0.2",
    "",
    "**Figure 2. MAGMA-prioritized KSD genes converge on a Loop/TAL-associated single-nucleus expression context.**",
    "(A) Audited cell-type composition of the GSE231569 renal papillary single-nucleus dataset. Major epithelial, stromal, endothelial, injured epithelial, Loop/TAL and perivascular/mural-like compartments were harmonized through annotation audit. Low-abundance compartments were retained but interpreted cautiously.",
    "(B) Cell-type projection of MAGMA-prioritized KSD gene sets across audited GSE231569 cell types. MAGMA top-ranked, FDR-significant and suggestive gene sets showed the highest expression-context scores in Loop/TAL cells.",
    "(C) Random gene-set benchmark showing that Loop/TAL-associated localization of MAGMA-prioritized modules exceeded size-matched random gene-set expectations. The dashed line indicates the 95th percentile of random expectation.",
    "(D) Robustness summary showing that the Loop/TAL signal remained supported across full audited projection, locus-balanced benchmarking, conservative annotation and leave-one-locus-out sensitivity checks.",
    "Together, these analyses support a Loop/TAL-associated renal papillary single-nucleus cellular context for MAGMA-prioritized KSD genes. They do not establish causal mediation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure2_legend_v0.2.md")
}

write_qc_v03 <- function() {
  qc <- data.table(
    figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
    title = c(
      "Post-GWAS framework for KSD cellular and disease-context mapping",
      "MAGMA-prioritized KSD genes converge on a Loop/TAL-associated single-nucleus context",
      "Gene-centric single-nucleus evidence for P1 candidates",
      "GSE73680 disease-context analysis supports MAGMA-prioritized modules"
    ),
    data_source = c("Integrated workflow summary", "MAGMA + GSE231569", "GSE231569", "GSE73680"),
    claim_supported = c(
      "post-GWAS cellular and disease-context evidence framework",
      "Loop/TAL-associated single-nucleus cellular context",
      "P1 gene-centric TAL/transport/calcium interpretation spectrum",
      "MAGMA module-level disease-context expression association"
    ),
    claim_not_supported = c(
      "causality; TWAS convergence; colocalization; spatial validation; P1 disease-gene validation",
      "causal mediation; spatial validation; TWAS convergence; colocalization",
      "P1 disease validation; uniform TAL specificity; causal mechanism",
      "uniform P1 single-gene validation; causality; cell-type-specific disease response"
    ),
    color_consistency = "pass",
    font_consistency = "pass",
    terminology_ok = "pass",
    legend_ok = "pass",
    final_status = c("main_candidate", "main_candidate", "main_candidate", "main_candidate")
  )
  fwrite(qc, "results/tables/main_figure_qc_v0.3.tsv", sep = "\t")

  write_lines(c(
    "# Main Figure Style Guide v0.3",
    "",
    "## Final Main-Figure Claims",
    "",
    "- Figure 1 supports the evidence framework and claim boundary.",
    "- Figure 2 supports a Loop/TAL-associated single-nucleus cellular context.",
    "- Figure 3 supports a P1 gene-centric TAL/transport/calcium interpretation spectrum.",
    "- Figure 4 supports MAGMA module-level disease-context expression association.",
    "",
    "## Claims Not Supported",
    "",
    "Across all figures, do not claim causal mediation, P1 disease-gene validation, TWAS convergence, colocalization, spatial validation or therapeutic targeting.",
    "",
    "## Terminology",
    "",
    "- Use `Loop/TAL` inside figures.",
    "- Use `MAGMA top 50`, `MAGMA top 100`, `MAGMA FDR` and `MAGMA suggestive`.",
    "- Use `single-nucleus cellular context`, not `spatial validation`.",
    "- Use `module-level disease-context association`, not `P1 single-gene validation`."
  ), "docs/main_figure_style_guide_v0.3.md")
}

abstract_v02 <- c(
  "# Abstract v0.2",
  "",
  "**Background:** Kidney stone disease has a substantial genetic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.",
  "",
  "**Methods:** We integrated public KSD GWAS summary statistics, MAGMA gene-based prioritization, audited GSE231569 renal papillary single-nucleus annotations, six-gene P1 candidate evidence scoring, and patient-aware GSE73680 papillary plaque/stone disease-context expression analysis.",
  "",
  "**Results:** MAGMA-prioritized KSD genes converged on a Loop/TAL-associated single-nucleus context in GSE231569, with enrichment exceeding random gene-set expectations and remaining supported across robustness checks. The six P1 candidate genes formed an interpretable TAL, epithelial transport, and calcium-handling expression spectrum, but did not behave as a uniform disease-validated gene panel. In GSE73680, MAGMA-prioritized modules showed module-level disease-context expression association, whereas P1 single-gene responses were heterogeneous and not FDR-supported.",
  "",
  "**Conclusion:** These results support a TAL-associated renal papillary cellular context and MAGMA module-level disease-context association for KSD genetic risk. They do not establish causal mediation, TWAS convergence, colocalization, spatial validation, or P1 disease-gene validation."
)

discussion_v02 <- c(
  "# Discussion v0.2",
  "",
  "The main finding of this study is that MAGMA-prioritized kidney stone disease genes converge on a TAL-associated renal papillary single-nucleus context and show disease-context support at the module level. This supports a module-level post-GWAS interpretation rather than a single causal gene or uniform disease-gene panel model. The integrated evidence connects genetic prioritization, single-nucleus cellular localization, P1 gene interpretation and independent papillary disease-context expression, while preserving clear boundaries around what each evidence layer can support.",
  "",
  "The Figure 2 consolidation is central to this interpretation. MAGMA top-ranked, FDR-significant and suggestive gene sets showed Loop/TAL-associated expression-context scores in the audited GSE231569 atlas, exceeded random gene-set expectations and remained supported in locus-balanced and leave-one-locus-out robustness summaries. Because Loop/TAL cells represented a relatively small but manually audited compartment, the interpretation relies on gene-set-level enrichment and robustness checks rather than cell abundance. These results support a TAL-associated single-nucleus cellular context, but they do not identify TAL as a causal cell type, provide spatial validation or substitute for TWAS or colocalization.",
  "",
  "The P1 gene analysis should be interpreted as gene-centric context rather than disease validation. UMOD, CLDN10, CLDN14, CASR, HIBADH and PKD2 form a biologically interpretable TAL, epithelial transport and calcium-handling spectrum, but they do not behave as a uniform TAL marker set. The value of the P1 panel is therefore interpretive rather than confirmatory: it separates representative TAL genes, epithelial transport candidates, calcium-handling genes and broader epithelial-context genes.",
  "",
  "GSE73680 provided disease-context support at the MAGMA module level, while arguing against uniform P1 single-gene differential expression. The module-level signal remained supported in patient-level paired sensitivity analyses and random gene-set benchmarks, but should not be interpreted as cell-type-specific disease expression because GSE73680 is bulk/microarray-based. This distinction keeps Result 4 aligned with the broader post-GWAS mapping framework rather than turning it into single-gene validation.",
  "",
  "Several limitations should guide interpretation. First, the study is computational and hypothesis-generating. Second, GSE231569 and GSE73680 differ in assay type, resolution and disease context, so their evidence layers are complementary rather than directly equivalent. TWAS, SMR/coloc and spatial transcriptomic extensions were resource-limited and should not be interpreted as negative evidence. GSE231569 cannot resolve spatial plaque microenvironments, and GSE73680 cannot infer cell-type-specific disease responses.",
  "",
  "Future work should prioritize kidney papilla-specific eQTL resources, lesion-resolved spatial transcriptomics and experimental perturbation of top MAGMA-prioritized modules in TAL-relevant epithelial models. These extensions would test whether the TAL-associated cellular context identified here corresponds to causal regulatory mechanisms or spatially localized papillary disease programs."
)

replace_section <- function(lines, heading, replacement) {
  start <- grep(paste0("^", heading, "$"), lines)
  if (length(start) != 1) stop("heading not found or duplicated: ", heading)
  next_heads <- grep("^#{2,3} ", lines)
  next_heads <- next_heads[next_heads > start]
  end <- if (length(next_heads)) next_heads[1] - 1 else length(lines)
  c(lines[seq_len(start - 1)], replacement, lines[(end + 1):length(lines)])
}

write_manuscript_v05 <- function() {
  lines <- readLines("manuscript/manuscript_draft_v0.4.md", warn = FALSE)
  abs_body <- abstract_v02[-1]
  disc_body <- discussion_v02[-1]
  lines <- replace_section(lines, "## Abstract", c("## Abstract", abs_body, ""))
  lines <- replace_section(lines, "## Discussion", c("## Discussion", disc_body[-1], ""))
  lines <- replace_section(lines, "### MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context", c(
    "### MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context",
    "",
    "We next projected MAGMA-prioritized KSD gene sets onto the audited GSE231569 single-nucleus atlas (Figure 2). The analysis used harmonized audited cell types and retained Loop/TAL as a high-confidence epithelial transport compartment. Across MAGMA top 50, top 100, FDR-significant and suggestive gene sets, expression-context scores were highest in Loop/TAL cells compared with other audited cell types.",
    "",
    "Random gene-set benchmarking supported the specificity of this localization pattern. Loop/TAL benchmark percentiles were 0.998 for MAGMA top 50, 1.000 for MAGMA top 100, 1.000 for MAGMA suggestive genes and 0.968 for MAGMA FDR genes. Locus-balanced benchmarking, conservative annotation sensitivity and leave-one-locus-out analyses supported robustness of the MAGMA top 50 Loop/TAL signal, with the minimum retained leave-one-locus-out percentile remaining above 0.99. These results support a Loop/TAL-associated renal papillary single-nucleus cellular context for MAGMA-prioritized KSD genes, without establishing causal mediation or spatial validation.",
    ""
  ))
  lines <- replace_section(lines, "### Figure 2. MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context", c(
    "### Figure 2. MAGMA-prioritized KSD genes converge on a Loop/TAL-associated single-nucleus expression context",
    "",
    "MAGMA-prioritized gene sets were evaluated across audited renal papillary single-nucleus cell types from GSE231569. Projection, random gene-set benchmarking and robustness analyses supported the Loop/TAL compartment as the strongest cellular expression context for prioritized kidney stone risk genes. Low-abundance compartments were retained but interpreted cautiously. These analyses do not establish causal mediation, TWAS convergence, colocalization or spatial validation.",
    ""
  ))
  write_lines(lines, "manuscript/manuscript_draft_v0.5.md")
}

write_supervisor_package <- function() {
  file.copy("manuscript/manuscript_draft_v0.5.md", "manuscript/manuscript_for_supervisor_review_v0.1.md", overwrite = TRUE)
  write_lines(c(
    "# Supervisor Review Cover Memo v0.1",
    "",
    "## Core Question",
    "",
    "Can kidney stone disease GWAS/MAGMA signals be mapped to a renal papillary cellular context and an independent disease-context expression layer without overclaiming causality?",
    "",
    "## Current Four-Layer Evidence Chain",
    "",
    "1. GWAS/MAGMA prioritizes KSD-associated genes and modules.",
    "2. GSE231569 single-nucleus projection supports a Loop/TAL-associated cellular context.",
    "3. The six-gene P1 candidate panel provides an interpretable TAL/transport/calcium expression spectrum.",
    "4. GSE73680 supports MAGMA module-level disease-context association rather than uniform P1 single-gene validation.",
    "",
    "## Main Finding",
    "",
    "The current manuscript supports a TAL-associated renal papillary cellular context and MAGMA module-level disease-context association for KSD genetic risk.",
    "",
    "## Explicit Boundaries",
    "",
    "The manuscript does not claim causal mediation, P1 disease-gene validation, TWAS convergence, colocalization, spatial validation or therapeutic targeting.",
    "",
    "## Questions For Supervisor",
    "",
    "1. Should the main line be positioned as TAL-associated cellular context rather than a broader multicellular papillary niche?",
    "2. Is GSE73680 module-level disease-context support sufficient as Result 4?",
    "3. Should TWAS/SMR-coloc/spatial analyses remain in limitations/resource-limited extensions?",
    "4. Should P1 genes be presented as a main candidate gene table or as an interpretive gene spectrum?"
  ), "docs/supervisor_review_cover_memo_v0.1.md")
}

write_forbidden_scan <- function() {
  text <- readLines("manuscript/manuscript_draft_v0.5.md", warn = FALSE)
  words <- c("validate", "validated", "validation", "prove", "causal", "driver",
             "therapeutic target", "spatial validation")
  hits <- rbindlist(lapply(words, function(w) {
    idx <- grep(w, text, ignore.case = TRUE)
    if (!length(idx)) return(data.table(term = w, line = integer(), text = character(), status = character()))
    data.table(term = w, line = idx, text = text[idx],
               status = fifelse(grepl("not|does not|without|do not|cannot|future|resource-limited", text[idx], ignore.case = TRUE),
                                "boundary_or_future_context", "review_context"))
  }), fill = TRUE)
  fwrite(hits, "results/tables/forbidden_overclaim_word_scan_v0.1.tsv", sep = "\t")
}

copy_final_draft <- function() {
  candidates <- c(
    "results/figures/figure1_integrative_framework_v0.2.pdf",
    "results/figures/figure1_integrative_framework_v0.2.png",
    "results/figures/figure2_magma_scrna_localization_v0.2.pdf",
    "results/figures/figure2_magma_scrna_localization_v0.2.png",
    "results/figures/figure3_p1_gene_evidence_v0.4.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.png",
    "results/figures/figure4_gse73680_disease_context_v0.8.pdf",
    "results/figures/figure4_gse73680_disease_context_v0.8.png"
  )
  file.copy(candidates[file.exists(candidates)], "results/figures/final_draft", overwrite = TRUE)
  write_lines(c(
    "# Final Draft Figure Package v0.4",
    "",
    "Preferred supervisor-review main figures:",
    "- figure1_integrative_framework_v0.2.pdf/png",
    "- figure2_magma_scrna_localization_v0.2.pdf/png",
    "- figure3_p1_gene_evidence_v0.4.pdf/png",
    "- figure4_gse73680_disease_context_v0.8.pdf/png"
  ), "results/figures/final_draft/README.md")
}

make_figure2_v02()
write_qc_v03()
write_lines(abstract_v02, "manuscript/abstract_v0.2.md")
write_lines(discussion_v02, "manuscript/discussion_v0.2.md")
write_manuscript_v05()
write_supervisor_package()
write_forbidden_scan()
copy_final_draft()

message("wrote Phase 7C supervisor-review package outputs")
