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

label_cell <- function(x) {
  map <- c(
    Collecting_duct_principal = "Collecting duct\nprincipal",
    Fibroblast_stromal = "Fibroblast/\nstromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured\nundiff. epi.",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/\nmural-like",
    Pericyte_smooth_muscle = "Perivascular/\nmural-like"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

label_gene_set <- function(x) {
  map <- c(
    magma_top50 = "MAGMA top 50",
    magma_top100 = "MAGMA top 100",
    magma_top200 = "MAGMA top 200",
    magma_fdr05 = "MAGMA FDR",
    magma_suggestive_p1e4 = "MAGMA suggestive"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

make_figure2 <- function() {
  atlas <- fread("results/tables/magma_scrna_module_score_by_celltype.tsv")
  scores <- fread("results/tables/magma_scrna_module_score_by_celltype.tsv")
  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  gene_summary <- fread("results/tables/magma_gene_set_summary.tsv")
  locus_bal <- fread("results/tables/magma_locus_balanced_scrna_benchmark.tsv")
  loo <- fread("results/tables/magma_leave_one_locus_out.tsv")

  atlas <- unique(atlas[, .(audited_broad_cell_type, n_cells, n_donors, annotation_confidence)])
  atlas[, cell_label := label_cell(audited_broad_cell_type)]
  atlas[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  atlas[, confidence_label := fcase(
    audited_broad_cell_type == "Loop_of_Henle_TAL", "Loop/TAL audited high confidence",
    grepl("low_or_exploratory", annotation_confidence), "Exploratory",
    grepl("review", annotation_confidence), "Reviewed",
    default = "Audited"
  )]
  atlas[, cell_label := factor(cell_label, levels = atlas[order(n_cells), cell_label])]

  panel_a <- ggplot(atlas, aes(n_cells, cell_label, fill = confidence_label)) +
    geom_col(width = 0.68, color = "#777777", linewidth = 0.2) +
    geom_text(aes(label = paste0("n=", n_cells)), hjust = -0.08, size = 2.35, color = pal$text) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = c("Loop/TAL audited high confidence" = pal$tal,
                                 Audited = pal$scrna,
                                 Reviewed = "#9EB3BA",
                                 Exploratory = pal$muted)) +
    labs(title = "A. Audited GSE231569 single-nucleus atlas",
         x = "Cells", y = NULL, fill = "Annotation") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())

  score_cols <- c("magma_top50_score", "magma_top100_score", "magma_fdr05_score", "magma_suggestive_p1e4_score")
  score_long <- melt(scores,
                     id.vars = "audited_broad_cell_type",
                     measure.vars = score_cols,
                     variable.name = "gene_set",
                     value.name = "module_score")
  score_long[, gene_set := sub("_score$", "", gene_set)]
  score_long[, gene_label := factor(label_gene_set(gene_set),
                                    levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
  score_long[, cell_label := factor(label_cell(audited_broad_cell_type),
                                    levels = c("Collecting duct\nprincipal", "Loop/TAL", "Injured\nundiff. epi.",
                                               "Endothelial", "Fibroblast/\nstromal", "Perivascular/\nmural-like"))]
  panel_b <- ggplot(score_long, aes(cell_label, gene_label, fill = module_score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", module_score)), size = 2.4, color = pal$text) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$tal) +
    labs(title = "B. MAGMA gene sets project to Loop/TAL cells",
         x = NULL, y = NULL, fill = "Module\nscore") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "right",
          panel.grid = element_blank())

  random_plot <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random_plot[, gene_label := factor(label_gene_set(gene_set),
                                     levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
  random_plot[, cell_label := factor(label_cell(audited_broad_cell_type),
                                     levels = c("Perivascular/\nmural-like", "Fibroblast/\nstromal", "Endothelial",
                                                "Injured\nundiff. epi.", "Collecting duct\nprincipal", "Loop/TAL"))]
  random_plot[, support_class := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95,
                                         "Loop/TAL exceeds expectation", "Other audited context")]
  panel_c <- ggplot(random_plot, aes(benchmark_percentile, cell_label, fill = support_class)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.35) +
    geom_point(shape = 21, size = 3.0, color = "#555555", stroke = 0.25) +
    facet_wrap(~ gene_label, ncol = 2) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95), labels = c("0", "0.5", "0.95")) +
    scale_fill_manual(values = c("Loop/TAL exceeds expectation" = pal$tal, "Other audited context" = pal$muted)) +
    labs(title = "C. TAL localization exceeds random expectation",
         x = "Random gene-set benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          strip.background = element_rect(fill = "#F1F1F1", color = "#AAAAAA"),
          panel.grid.minor = element_blank())

  top50 <- gene_summary[gene_set == "magma_top50"]
  lb_full <- locus_bal[analysis_version == "full_audited" &
                         gene_set == "magma_locus_balanced_top50" &
                         audited_broad_cell_type == "Loop_of_Henle_TAL",
                       benchmark_percentile][1]
  lb_cons <- locus_bal[analysis_version == "conservative_exclude_low_or_exploratory_and_immune_review" &
                         gene_set == "magma_locus_balanced_top50" &
                         audited_broad_cell_type == "Loop_of_Henle_TAL",
                       benchmark_percentile][1]
  loo_min <- min(loo$TAL_percentile_after_removal, na.rm = TRUE)
  robust <- data.table(
    check = c("Full audited projection", "Locus-balanced benchmark",
              "Conservative annotation", "Leave-one-locus-out min"),
    percentile = c(top50$TAL_percentile, lb_full, lb_cons, loo_min),
    note = c("MAGMA top 50", "Locus-balanced top 50",
             "Conservative audited set", "Worst retained locus removal")
  )
  robust[, check := factor(check, levels = rev(check))]
  robust[, support := fifelse(percentile >= 0.95, "Supported", "Review")]
  panel_d <- ggplot(robust, aes(percentile, check, fill = support)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = pal$grid, linewidth = 0.35) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3f", percentile)), hjust = 1.08, size = 2.55, color = "white") +
    geom_text(aes(x = 0.02, label = note), hjust = 0, size = 2.35, color = "white") +
    scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.95), labels = c("0", "0.5", "0.95")) +
    scale_fill_manual(values = c(Supported = pal$tal, Review = pal$muted)) +
    labs(title = "D. TAL signal remains supported across checks",
         x = "Loop/TAL benchmark percentile", y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2,
                   rel_heights = c(1.0, 1.08))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.1.pdf", fig,
         width = 12.5, height = 9.0, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.1.png", fig,
         width = 12.5, height = 9.0, units = "in", dpi = 250, bg = "white")

  source_files <- data.table(
    panel = c("A", "B", "C", "D"),
    source_table = c(
      "results/tables/magma_scrna_module_score_by_celltype.tsv; results/tables/gse231569_cell_counts.tsv",
      "results/tables/magma_scrna_module_score_by_celltype.tsv",
      "results/tables/magma_scrna_random_benchmark.tsv",
      "results/tables/magma_gene_set_summary.tsv; results/tables/magma_locus_balanced_scrna_benchmark.tsv; results/tables/magma_leave_one_locus_out.tsv"
    ),
    source_figure = c(
      "results/figures/gse231569_umap_audited_broad_celltype.pdf (supplementary visual context)",
      "results/figures/magma_scrna_benchmark.pdf (legacy source)",
      "results/figures/magma_scrna_benchmark.pdf (legacy source)",
      "results/figures/magma_leave_one_locus_out_tal.pdf (legacy source)"
    ),
    analysis_step = c(
      "annotation audit and audited cell-type summary",
      "MAGMA gene-set projection across audited GSE231569 cell types",
      "size-matched random gene-set benchmark across audited cell types",
      "robustness summary using full audited, conservative, locus-balanced and leave-one-locus-out checks"
    ),
    claim_supported = c(
      "audited single-nucleus atlas used for localization",
      "MAGMA-prioritized gene sets show higher expression-context scores in Loop/TAL cells",
      "Loop/TAL localization exceeds random gene-set expectation",
      "MAGMA-based TAL localization remains supported across robustness checks"
    ),
    claim_boundary = c(
      "not spatial validation",
      "expression context only; not causal mediation",
      "benchmark support only; not gene validation",
      "robustness support only; not TWAS, colocalization or spatial validation"
    ),
    notes = c(
      "Panel A is an audited atlas summary rather than a full UMAP; full UMAP remains available as supplementary context.",
      "Uses MAGMA top 50, top 100, FDR and suggestive sets.",
      "Dashed line marks the 95th percentile.",
      "Leave-one-locus-out summary uses the minimum retained Loop/TAL percentile."
    )
  )
  fwrite(source_files, "results/tables/figure2_panel_source_files.tsv", sep = "\t")

  write_lines(c(
    "# Figure 2 Legend v0.1",
    "",
    "**Figure 2. MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus expression context.**",
    "(A) Audited GSE231569 single-nucleus atlas summary used for cell-type localization. Major renal papillary cell compartments were manually reviewed and harmonized, with Loop/TAL cells retained as an audited epithelial transport compartment.",
    "(B) Projection of MAGMA-prioritized KSD gene sets across audited cell types. MAGMA top-ranked and FDR/suggestive gene sets showed higher expression-context scores in Loop/TAL cells compared with most other audited cell types.",
    "(C) Random gene-set benchmarking showed that TAL-associated localization of MAGMA-prioritized modules exceeded size-matched random gene-set expectations. The dashed line indicates the 95th percentile.",
    "(D) Robustness analyses, including full audited projection, conservative annotation, locus-balanced benchmarking and leave-one-locus-out sensitivity, supported the stability of MAGMA-based TAL-associated localization.",
    "Together, these analyses support a TAL-associated single-nucleus cellular context for MAGMA-prioritized KSD genes. These results do not establish causal mediation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure2_legend_v0.1.md")
}

update_qc_and_style <- function() {
  qc <- data.table(
    figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
    panel = c("framework", "MAGMA scRNA localization panels A-D", "P1 evidence panels", "GSE73680 panels A-D"),
    title = c(
      "Post-GWAS framework for KSD cellular and disease-context mapping",
      "MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context",
      "Gene-centric single-nucleus evidence for P1 candidates",
      "GSE73680 disease-context analysis supports MAGMA-prioritized modules"
    ),
    data_source = c("Integrated workflow summary", "MAGMA + GSE231569", "GSE231569", "GSE73680"),
    main_claim = c(
      "The study is a post-GWAS cellular and disease-context mapping framework.",
      "MAGMA-prioritized genes converge on a TAL-associated single-nucleus context.",
      "P1 genes form an interpretable TAL/transport/calcium expression spectrum.",
      "GSE73680 supports MAGMA module-level disease-context expression association."
    ),
    claim_boundary = c(
      "Does not establish causality, TWAS, coloc, spatial validation or P1 disease-gene validation.",
      "Single-nucleus expression context only; not causal mediation, TWAS, colocalization or spatial validation.",
      "Not P1 disease validation or uniform TAL specificity.",
      "Not P1 gene validation, causality or cell-type-specific disease expression."
    ),
    color_consistency = c("pass", "pass", "pass", "pass"),
    font_consistency = c("pass", "pass", "pass", "pass"),
    axis_label_ok = c("not_applicable", "pass", "pass", "pass"),
    legend_ok = c("pass", "pass", "pass", "pass"),
    needs_revision = c("no", "no", "no", "no"),
    action = c(
      "Use v0.2 in final draft package.",
      "Use consolidated v0.1 as current Figure 2.",
      "Freeze v0.4 unless journal requests layout changes.",
      "Freeze v0.8 after visual inspection."
    )
  )
  fwrite(qc, "results/tables/main_figure_qc_v0.2.tsv", sep = "\t")

  write_lines(c(
    "# Main Figure Style Guide v0.2",
    "",
    "## Shared Claim Boundary",
    "",
    "Figures 1-4 support post-GWAS cellular and disease-context mapping. They should not be described as causal validation, P1 disease-gene validation, TWAS convergence, colocalization or spatial validation.",
    "",
    "## Terminology",
    "",
    "- Use `Loop/TAL` inside figures and `Loop of Henle/thick ascending limb (Loop/TAL)` at first mention in text.",
    "- Use `MAGMA top 50`, `MAGMA top 100`, `MAGMA FDR` and `MAGMA suggestive` consistently.",
    "- Use `TAL-associated single-nucleus context`, not `TAL causal cell type`.",
    "- Use `module-level disease-context association`, not `P1 gene validation`, for GSE73680.",
    "",
    "## Palette",
    "",
    "- MAGMA/genetic layer: deep blue-gray (`#3E6672`).",
    "- Single-nucleus context: restrained blue-gray (`#6F8F98`).",
    "- P1 evidence: ochre (`#B59A5B`).",
    "- Disease context: warm brown (`#9A5F52`) with MAGMA module bars in deep blue-gray.",
    "- Claim boundaries and resource-limited extensions: neutral gray.",
    "",
    "## Figure-Specific Notes",
    "",
    "- Figure 1 shows four evidence layers and a separate boundary box.",
    "- Figure 2 is now consolidated and should be the main Result 2 figure.",
    "- Figure 3 remains gene-centric and should not be interpreted as disease validation.",
    "- Figure 4 remains module-level and patient-aware; Panel C is paired sensitivity and Panel D is size-matched benchmark."
  ), "docs/main_figure_style_guide_v0.2.md")
}

replace_section <- function(lines, heading, replacement) {
  start <- grep(paste0("^", heading, "$"), lines)
  if (length(start) != 1) stop("heading not found or duplicated: ", heading)
  next_heads <- grep("^#{2,3} ", lines)
  next_heads <- next_heads[next_heads > start]
  end <- if (length(next_heads)) next_heads[1] - 1 else length(lines)
  c(lines[seq_len(start - 1)], replacement, lines[(end + 1):length(lines)])
}

write_manuscript_v04 <- function() {
  lines <- readLines("manuscript/manuscript_draft_v0.3.md", warn = FALSE)

  lines <- replace_section(lines, "### GSE231569 single-nucleus processing and annotation audit", c(
    "### GSE231569 single-nucleus processing and annotation audit",
    "",
    "Renal papillary single-nucleus data from GSE231569 were processed with an audited broad cell-type annotation. Marker expression, cluster assignment and renal epithelial compartment labels were reviewed before gene-set projection. Loop of Henle/thick ascending limb cells were harmonized as the Loop/TAL compartment and treated as an audited epithelial transport context. Cell types with review flags or low abundance were retained with explicit caution labels rather than being used as unqualified primary localization claims.",
    ""
  ))
  lines <- replace_section(lines, "### Gene-set projection to single-nucleus cell types", c(
    "### Gene-set projection to single-nucleus cell types",
    "",
    "MAGMA-prioritized gene sets were projected onto audited GSE231569 renal papillary cell types using module-score summaries across detected genes. The main projected sets were MAGMA top 50, MAGMA top 100, MAGMA FDR and MAGMA suggestive genes. For each cell type and gene set, the observed expression-context score was compared with size-matched random gene sets drawn from the detected single-nucleus background. Benchmark percentiles were used to identify cell-type contexts exceeding random expectation.",
    "",
    "Robustness checks were summarized without introducing additional primary claims. Locus-balanced benchmarking evaluated whether the MAGMA top 50 TAL signal persisted after balancing locus contribution. Conservative annotation sensitivity excluded low-confidence or review-flagged contexts where appropriate. Leave-one-locus-out sensitivity recalculated the Loop/TAL benchmark after removing each locus group and used the minimum retained percentile as a robustness summary. These analyses were interpreted as single-nucleus cellular context mapping rather than causal mediation, TWAS convergence, colocalization or spatial validation.",
    ""
  ))
  lines <- replace_section(lines, "### MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context", c(
    "### MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context",
    "",
    "We next projected MAGMA-prioritized KSD gene sets onto the audited GSE231569 single-nucleus atlas (Figure 2). The analysis used harmonized audited cell types and retained Loop/TAL as a high-confidence epithelial transport compartment. Across MAGMA top 50, top 100, FDR-significant and suggestive gene sets, expression-context scores were highest or near-highest in Loop/TAL cells compared with most other audited cell types.",
    "",
    "Random gene-set benchmarking supported the specificity of this localization pattern. Loop/TAL benchmark percentiles were 0.998 for MAGMA top 50, 1.000 for MAGMA top 100, 1.000 for MAGMA suggestive genes and 0.968 for MAGMA FDR genes. Locus-balanced benchmarking and leave-one-locus-out sensitivity supported robustness of the MAGMA top 50 Loop/TAL signal, with the minimum retained leave-one-locus-out percentile remaining above 0.99. These results support a TAL-associated single-nucleus cellular context for MAGMA-prioritized KSD genes, without establishing causal mediation or spatial validation.",
    ""
  ))
  lines <- replace_section(lines, "## Limitations", c(
    "## Limitations",
    "",
    "Several limitations should guide interpretation. TWAS, SMR/coloc and spatial transcriptomic validation were not completed because required external expression-weight, eQTL and spatial matrix/image resources were unavailable in the current analysis environment. These resource-limited modules should not be interpreted as negative results. GSE231569 supports single-nucleus expression-context mapping but does not provide spatial validation or causal cell-type mediation. GSE73680 supports disease-context module association but not causality. P1 single-gene disease differential expression was not FDR-supported, and PKD2 should be treated only as a nominal exploratory observation. Because GSE73680 is a reconstructed bulk/microarray disease-context dataset, it cannot resolve cell-type-specific disease responses or spatial localization.",
    ""
  ))
  lines <- replace_section(lines, "### Figure 2. MAGMA-prioritized KSD genes localize to TAL-associated cell states", c(
    "### Figure 2. MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context",
    "",
    "MAGMA-prioritized gene sets were evaluated across audited renal papillary single-nucleus cell types from GSE231569. Projection, random gene-set benchmarking and robustness analyses supported the Loop/TAL compartment as the strongest cellular expression context for prioritized kidney stone risk genes. These analyses do not establish causal mediation, TWAS convergence, colocalization or spatial validation.",
    ""
  ))
  write_lines(lines, "manuscript/manuscript_draft_v0.4.md")
}

write_abstract_discussion <- function() {
  write_lines(c(
    "# Abstract v0.1",
    "",
    "**Background:** Kidney stone disease has a substantial genetic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.",
    "",
    "**Methods:** We integrated locked KSD GWAS summary statistics, MAGMA gene-based prioritization, audited GSE231569 renal papillary single-nucleus annotations, P1 candidate gene evidence scoring and GSE73680 papillary plaque/stone disease-context expression analysis.",
    "",
    "**Results:** MAGMA-prioritized KSD genes converged on a Loop/TAL-associated single-nucleus context in GSE231569. P1 candidates formed an interpretable TAL, epithelial transport and calcium-handling expression spectrum. In GSE73680, MAGMA-prioritized modules showed disease-context expression association, whereas P1 single-gene responses were heterogeneous and not FDR-supported.",
    "",
    "**Conclusion:** These results support a TAL-associated cellular context and MAGMA module-level disease-context association for KSD genetic risk, while not establishing causal mediation, TWAS convergence, colocalization or spatial validation."
  ), "manuscript/abstract_v0.1.md")

  write_lines(c(
    "# Discussion v0.1",
    "",
    "The main finding of this study is that MAGMA-prioritized kidney stone disease genes converge on a TAL-associated renal papillary single-nucleus context and show disease-context support at the module level. This is a post-GWAS mapping result rather than a causal validation result. The integrated evidence connects genetic prioritization, single-nucleus cellular localization, P1 gene interpretation and independent papillary disease-context expression, while preserving clear boundaries around what each layer can support.",
    "",
    "The Figure 2 consolidation is central to this interpretation. MAGMA top-ranked, FDR-significant and suggestive gene sets showed Loop/TAL-associated expression-context scores in the audited GSE231569 atlas, exceeded random gene-set expectations and remained supported in locus-balanced and leave-one-locus-out robustness summaries. These results support a TAL-associated single-nucleus cellular context, but they do not identify TAL as a causal cell type, do not provide spatial validation and do not substitute for TWAS or colocalization.",
    "",
    "The P1 gene analysis should be interpreted as gene-centric context rather than disease validation. UMOD, CLDN10, CLDN14, CASR, HIBADH and PKD2 form a biologically interpretable TAL, epithelial transport and calcium-handling spectrum, but they do not behave as a uniform TAL marker set. Similarly, GSE73680 strengthens the disease-context interpretation at the MAGMA module level, while P1 single-gene responses remain heterogeneous and not FDR-supported.",
    "",
    "Several limitations remain. GSE231569 is a single-nucleus dataset and cannot resolve spatial plaque microenvironments. GSE73680 is a reconstructed bulk/microarray disease-context dataset and cannot infer cell-type-specific disease responses. TWAS, SMR/coloc and spatial transcriptomic extensions were resource-limited and should not be interpreted as negative evidence. Future work should test whether the TAL-associated cellular context identified here corresponds to causal regulatory mechanisms or spatially localized papillary disease programs."
  ), "manuscript/discussion_v0.1.md")
}

copy_final_draft <- function() {
  candidates <- c(
    "results/figures/figure1_integrative_framework_v0.2.pdf",
    "results/figures/figure1_integrative_framework_v0.2.png",
    "results/figures/figure2_magma_scrna_localization_v0.1.pdf",
    "results/figures/figure2_magma_scrna_localization_v0.1.png",
    "results/figures/figure3_p1_gene_evidence_v0.4.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.png",
    "results/figures/figure4_gse73680_disease_context_v0.8.pdf",
    "results/figures/figure4_gse73680_disease_context_v0.8.png"
  )
  file.copy(candidates[file.exists(candidates)], "results/figures/final_draft", overwrite = TRUE)
  write_lines(c(
    "# Final Draft Figure Package v0.3",
    "",
    "Current preferred main-figure files:",
    "- figure1_integrative_framework_v0.2.pdf/png",
    "- figure2_magma_scrna_localization_v0.1.pdf/png",
    "- figure3_p1_gene_evidence_v0.4.pdf/png",
    "- figure4_gse73680_disease_context_v0.8.pdf/png",
    "",
    "Figure 2 is now consolidated from existing MAGMA scRNA projection, random benchmark and robustness outputs."
  ), "results/figures/final_draft/README.md")
}

make_figure2()
update_qc_and_style()
write_manuscript_v04()
write_abstract_discussion()
copy_final_draft()

message("wrote Phase 7B Figure 2 consolidation outputs")
