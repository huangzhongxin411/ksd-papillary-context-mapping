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
  boundary = "#4F4F4F",
  resource = "#F1F1F1",
  muted = "#D8D8D8",
  text = "#303030",
  grid = "#A9A9A9",
  panel = "#F7F8F8"
)

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

figure1_v02 <- function() {
  layers <- data.table(
    x = c(0.17, 0.39, 0.61, 0.83),
    y = rep(0.62, 4),
    layer = c("Layer 1", "Layer 2", "Layer 3", "Layer 4"),
    heading = c("Genetic prioritization", "Single-nucleus localization",
                "Gene-centric interpretation", "Disease-context association"),
    label = c("GWAS/MAGMA\nprioritization",
              "GSE231569\nsingle-nucleus localization",
              "P1 candidate gene\nevidence",
              "GSE73680 disease-context\nmodule analysis"),
    fill = c(pal$genetic, pal$scrna, pal$p1, pal$disease)
  )

  fig <- ggplot() +
    annotate("text", x = 0.50, y = 0.94,
             label = "Post-GWAS framework for KSD cellular and disease-context mapping",
             fontface = "bold", size = 5.0, color = pal$text) +
    geom_segment(data = layers[1:3],
                 aes(x = x + 0.09, xend = layers$x[2:4] - 0.09, y = y, yend = y),
                 arrow = arrow(length = unit(0.12, "in")), color = "#777777", linewidth = 0.45) +
    geom_label(data = layers, aes(x = x, y = y, label = label, fill = fill),
               color = "white", fontface = "bold", size = 3.45, lineheight = 0.95,
               linewidth = 0, label.padding = unit(0.28, "lines")) +
    geom_text(data = layers, aes(x = x, y = 0.765, label = layer),
              fontface = "bold", size = 3.1, color = pal$boundary) +
    geom_text(data = layers, aes(x = x, y = 0.725, label = heading),
              size = 2.8, color = pal$text) +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.26, ymax = 0.44,
             fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = 0.16, y = 0.395, label = "Supported inference:",
             fontface = "bold", hjust = 0, size = 3.3, color = pal$text) +
    annotate("text", x = 0.36, y = 0.395,
             label = "TAL-associated cellular context + MAGMA module-level disease-context association",
             hjust = 0, size = 3.2, color = pal$text) +
    annotate("text", x = 0.16, y = 0.315, label = "Not established:",
             fontface = "bold", hjust = 0, size = 3.3, color = pal$boundary) +
    annotate("text", x = 0.36, y = 0.315,
             label = "causality | TWAS convergence | colocalization | spatial validation | P1 disease-gene validation",
             hjust = 0, size = 3.1, color = pal$boundary) +
    annotate("rect", xmin = 0.56, xmax = 0.92, ymin = 0.055, ymax = 0.20,
             fill = pal$resource, color = "#A8A8A8", linewidth = 0.35) +
    annotate("text", x = 0.74, y = 0.165, label = "Resource-limited extensions",
             fontface = "bold", size = 3.35, color = pal$text) +
    annotate("text", x = 0.74, y = 0.115, label = "TWAS | SMR/coloc | spatial transcriptomics",
             size = 3.05, color = pal$text) +
    annotate("text", x = 0.74, y = 0.075, label = "prepared but not used for claims",
             size = 2.8, color = pal$boundary) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  ggsave("results/figures/figure1_integrative_framework_v0.2.pdf", fig,
         width = 12.8, height = 5.2, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.2.png", fig,
         width = 12.8, height = 5.2, units = "in", dpi = 260, bg = "white")

  write_lines(c(
    "# Figure 1 Legend v0.2",
    "",
    "**Figure 1. Post-GWAS framework for KSD cellular and disease-context mapping.**",
    "The study is organized as four evidence layers: GWAS/MAGMA gene prioritization, audited GSE231569 single-nucleus localization, P1 candidate gene evidence and GSE73680 disease-context module analysis.",
    "The supported inference is TAL-associated cellular context plus MAGMA module-level disease-context association. TWAS, SMR/coloc and spatial transcriptomic analyses were prepared as resource-limited extensions but were not used for manuscript claims.",
    "The framework does not establish causal mediation, P1 disease-gene validation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure1_legend_v0.2.md")
}

figure4_v08 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  modules <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
  rand_bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

  panel_a <- ggplot() +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.79, ymax = 0.93, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.56, ymax = 0.69, fill = pal$panel, color = "#8A8A8A", linewidth = 0.35) +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.34, ymax = 0.47, fill = pal$panel, color = "#8A8A8A", linewidth = 0.35) +
    annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.10, ymax = 0.24, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("text", x = 0.50, y = 0.86, label = "GSE73680 disease-context dataset", fontface = "bold", size = 3.6, color = pal$text) +
    annotate("text", x = 0.50, y = 0.625, label = "55 included samples: 27 control/adjacent, 28 plaque/stone papilla", size = 3.18, color = pal$text) +
    annotate("text", x = 0.50, y = 0.405, label = "29 patients, including 26 paired patients", size = 3.18, color = pal$text) +
    annotate("text", x = 0.50, y = 0.18, label = "Primary: limma duplicateCorrelation", fontface = "bold", size = 3.12, color = pal$text) +
    annotate("text", x = 0.50, y = 0.125, label = "Sensitivity: paired patient-level delta", fontface = "bold", size = 3.12, color = pal$text) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.78, yend = 0.70, arrow = arrow(length = unit(0.12, "in")), color = "#777777") +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.55, yend = 0.48, arrow = arrow(length = unit(0.12, "in")), color = "#777777") +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.33, yend = 0.25, arrow = arrow(length = unit(0.12, "in")), color = "#777777") +
    labs(title = "A. Patient-aware GSE73680 design") +
    theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 11),
          plot.background = element_rect(fill = "white", color = NA))

  p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No detectable")]
  p1[, p_label := sprintf("P=%.3f", p_value)]
  panel_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
    geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.45, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$p1, "No detectable" = pal$muted)) +
    labs(x = "Patient-level paired delta", y = NULL, fill = "Signal",
         subtitle = "No P1 gene reached FDR q <= 0.05; PKD2 nominal only",
         title = "B. No uniform P1 single-gene response") +
    theme_bw(base_size = 9.5) +
    theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, color = pal$boundary),
          legend.position = "bottom", panel.grid.minor = element_blank())

  module_order <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive",
                    "P1_core_TAL_candidates", "TAL_marker_set", "injury_remodeling_marker_set")
  mod_plot <- modules[module_name %in% module_order]
  mod_plot[, module_name := factor(module_name, levels = rev(module_order))]
  module_label_map <- c(
    MAGMA_top50 = "MAGMA top 50",
    MAGMA_top100 = "MAGMA top 100",
    MAGMA_FDR = "MAGMA FDR",
    MAGMA_suggestive = "MAGMA suggestive",
    P1_core_TAL_candidates = "P1 core candidates",
    TAL_marker_set = "TAL markers",
    injury_remodeling_marker_set = "Injury/remodeling"
  )
  mod_plot[, module_label := factor(module_label_map[as.character(module_name)],
                                    levels = rev(module_label_map[module_order]))]
  mod_plot[, signal_class := fcase(
    fdr <= 0.05 & module_name == "injury_remodeling_marker_set", "Disease-context program",
    fdr <= 0.05, "q <= 0.05",
    fdr <= 0.10, "Borderline",
    default = "No detectable"
  )]
  mod_plot[, q_label := fifelse(fdr <= 0.05, "q <= 0.05", sprintf("q == %.3f", fdr))]
  mod_plot[, label_x := fifelse(abs(paired_delta) < 0.03, 0.03, paired_delta)]
  mod_plot[, label_hjust := fifelse(abs(paired_delta) < 0.03, 0, 1.05)]
  panel_c <- ggplot(mod_plot, aes(paired_delta, module_label, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
    geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
    geom_text(aes(x = label_x, label = q_label, hjust = label_hjust), parse = TRUE, size = 2.35, color = pal$text) +
    scale_fill_manual(
      values = c("q <= 0.05" = pal$genetic, Borderline = pal$p1,
                 "No detectable" = pal$muted, "Disease-context program" = "#6B8E6E"),
      labels = c("q <= 0.05" = expression(q <= 0.05), Borderline = "Borderline",
                 "No detectable" = "No detectable", "Disease-context program" = "Disease-context program")
    ) +
    labs(x = "Paired module delta (plaque/stone papilla - control/adjacent)", y = NULL, fill = "Support",
         title = "C. MAGMA modules show paired disease-context shifts") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  bench_plot <- rand_bench[module_name %in% c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]
  bench_label_map <- c(
    P1_core_TAL_candidates = "P1 core candidates",
    MAGMA_top50 = "MAGMA top 50",
    MAGMA_top100 = "MAGMA top 100",
    MAGMA_FDR = "MAGMA FDR",
    MAGMA_suggestive = "MAGMA suggestive"
  )
  bench_plot[, module_label := factor(bench_label_map[as.character(module_name)],
                                      levels = bench_label_map[as.character(module_name[order(percentile)])])]
  bench_plot[, benchmark_class := fcase(
    benchmark_interpretation == "module response exceeds random expectation", "Exceeds random expectation",
    benchmark_interpretation == "moderate exploratory", "Moderate",
    default = "Background-like"
  )]
  bench_plot[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  panel_d <- ggplot(bench_plot, aes(percentile, module_label, fill = benchmark_class)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "#999999", linewidth = 0.3) +
    geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = emp_label), hjust = 1.05, size = 2.35, color = pal$text) +
    annotate("text", x = 0.94, y = 5.45, label = "95th percentile", hjust = 1,
             vjust = 0.5, size = 2.35, color = "#666666") +
    coord_cartesian(xlim = c(0, 1.05)) +
    scale_fill_manual(values = c("Exceeds random expectation" = pal$genetic, Moderate = pal$p1,
                                 "Background-like" = pal$muted)) +
    labs(x = "Percentile among size-matched random gene sets", y = NULL, fill = "Benchmark",
         title = "D. Size-matched random benchmark") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

  fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, align = "hv")
  ggsave("results/figures/figure4_gse73680_disease_context_v0.8.pdf", fig,
         width = 12.5, height = 8.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v0.8.png", fig,
         width = 12.5, height = 8.8, units = "in", dpi = 240, bg = "white")

  write_lines(c(
    "# Figure 4 Legend v0.8",
    "",
    "**Figure 4. GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression.**",
    "(A) Patient-aware analysis design for GSE73680. A total of 55 included samples, including 27 control/adjacent and 28 plaque/stone papilla samples, were analyzed from 29 patients, including 26 paired patients. Primary analyses used limma duplicateCorrelation, with paired patient-level delta analyses as sensitivity checks.",
    "(B) Paired patient-level responses of six P1 candidate genes. PKD2 showed a nominal exploratory response, but no P1 gene reached q≤0.05 after FDR correction.",
    "(C) Paired patient-level module responses. MAGMA-prioritized modules showed disease-context shifts, whereas the P1 core candidate module and TAL marker set did not show comparable support. The injury/remodeling marker set served as a disease-context reference program, not as a MAGMA-prioritized module.",
    "(D) Size-matched random gene-set benchmark showing that MAGMA module shifts exceeded random expectation, whereas the P1 core candidate module was background-like. The dashed line indicates the 95th percentile of size-matched random gene-set expectations.",
    "Size-matched random gene-set benchmarking was used as the main benchmark; expression-matched benchmarking was retained as a conservative sensitivity analysis and was not used as primary support. Together, these results support disease-context expression association for MAGMA-prioritized modules, rather than uniform P1 single-gene differential expression. These findings do not establish genetic causality, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure4_legend_v0.8.md")
}

figure_qc_outputs <- function() {
  qc <- data.table(
    figure = c("Figure 1", "Figure 2", "Figure 2", "Figure 3", "Figure 4"),
    panel = c("framework", "MAGMA scRNA benchmark", "leave-one-locus-out robustness", "P1 evidence panels", "GSE73680 panels A-D"),
    title = c(
      "Post-GWAS framework for KSD cellular and disease-context mapping",
      "MAGMA-prioritized KSD genes localize to TAL-associated cell states",
      "Robustness of MAGMA TAL localization",
      "Gene-centric single-nucleus evidence for P1 candidates",
      "GSE73680 disease-context analysis supports MAGMA-prioritized modules"
    ),
    data_source = c("Integrated workflow summary", "MAGMA + GSE231569", "MAGMA + GSE231569", "GSE231569", "GSE73680"),
    main_claim = c(
      "The study is a post-GWAS cellular and disease-context mapping framework.",
      "MAGMA-prioritized genes converge on a TAL-associated single-nucleus context.",
      "The TAL-associated MAGMA signal is not explained by a single locus alone.",
      "P1 genes form an interpretable TAL/transport/calcium expression spectrum.",
      "GSE73680 supports MAGMA module-level disease-context expression association."
    ),
    claim_boundary = c(
      "Does not establish causality, TWAS, coloc, spatial validation or P1 disease-gene validation.",
      "Cellular localization only; not causal mediation or spatial validation.",
      "Robustness check only; not causal proof.",
      "Not P1 disease validation or uniform TAL specificity.",
      "Not P1 gene validation, causality or cell-type-specific disease expression."
    ),
    color_consistency = c("pass", "review_needed", "review_needed", "pass", "pass"),
    font_consistency = c("pass", "review_needed", "review_needed", "pass", "pass"),
    axis_label_ok = c("not_applicable", "review_needed", "review_needed", "pass", "pass"),
    legend_ok = c("pass", "review_needed", "review_needed", "pass", "pass"),
    needs_revision = c("no", "yes", "yes", "no", "no"),
    action = c(
      "Use v0.2 in final draft package.",
      "Combine or restyle Figure 2 panels in later layout pass if target journal requires one composite figure.",
      "Pair with MAGMA scRNA benchmark as Figure 2 robustness panel.",
      "Freeze v0.4 unless journal requests layout changes.",
      "Freeze v0.8 after visual inspection."
    )
  )
  fwrite(qc, "results/tables/main_figure_qc_v0.1.tsv", sep = "\t")

  write_lines(c(
    "# Main Figure Style Guide v0.1",
    "",
    "## Shared Claim Boundary",
    "",
    "All main figures should support post-GWAS cellular and disease-context mapping. None should be described as causal validation, P1 disease-gene validation, TWAS convergence, colocalization or spatial validation.",
    "",
    "## Palette",
    "",
    "- GWAS/MAGMA: deep blue-gray (`#3E6672`).",
    "- Single-nucleus/scRNA context: medium blue-gray (`#6F8F98`) or related restrained blue-gray.",
    "- P1 gene evidence: ochre (`#B59A5B`).",
    "- GSE73680 disease context: warm brown (`#9A5F52`) with MAGMA module bars in deep blue-gray.",
    "- Claim boundary and unavailable extensions: neutral gray.",
    "",
    "## Typography and Layout",
    "",
    "- Keep panel titles short and claim-oriented.",
    "- Use consistent figure-internal font scale across Figures 1, 3 and 4.",
    "- Avoid using claim-boundary text as if it were an analysis result.",
    "- Keep resource-limited TWAS/SMR-coloc/spatial analyses visually separated from completed evidence layers.",
    "",
    "## Figure-Specific Notes",
    "",
    "- Figure 1 should show four evidence layers and a separate boundary box.",
    "- Figure 2 should not imply TAL causal mediation; use localization/context language.",
    "- Figure 3 should not imply P1 disease validation.",
    "- Figure 4 should state MAGMA module-level disease-context support and no P1 gene q≤0.05."
  ), "docs/main_figure_style_guide_v0.1.md")
}

supplementary_plans <- function() {
  write_lines(c(
    "# Supplementary Table Plan v0.1",
    "",
    "| Table | Title | Source files | Purpose | Claim boundary |",
    "|---|---|---|---|---|",
    "| Table S1 | GWAS QC and lead loci | `results/tables/phase1_gwas_qc_report.tsv`; `results/tables/phase1_2025_loci.tsv`; `results/tables/phase1_2025_lead_snps.tsv` | Documents GWAS input QC and retained lead loci. | Does not nominate causal variants. |",
    "| Table S2 | MAGMA gene-level results | `results/tables/magma_genes.tsv`; `results/tables/magma_qc_summary.tsv` | Provides full gene-level prioritization output. | Gene-level prioritization only. |",
    "| Table S3 | GSE231569 annotation audit | `results/tables/gse231569_marker_audit.tsv`; `results/tables/gse231569_cluster_marker_assignment_audit.tsv` | Documents cell-type annotation audit. | Annotation support, not causal mediation. |",
    "| Table S4 | MAGMA gene-set localization benchmark | `results/tables/magma_scrna_evidence_summary.tsv`; `results/tables/magma_scrna_random_benchmark.tsv`; `results/tables/magma_locus_balanced_scrna_benchmark.tsv` | Supports TAL-associated single-nucleus localization. | Not spatial validation. |",
    "| Table S5 | P1 candidate gene evidence | `results/tables/p1_tal_gene_interpretation_summary.tsv`; `results/tables/candidate_gene_tiers_v0.4.tsv` | Summarizes P1 evidence spectrum and gene tiers. | Not disease-gene validation. |",
    "| Table S6 | GSE73680 sample metadata and QC | `results/gse73680/tables/gse73680_analysis_design.tsv`; `results/gse73680/tables/gse73680_expression_qc_v2.tsv`; `results/gse73680/tables/gse73680_gene_mapping_qc.tsv` | Documents sample inclusion, mapping and QC. | Bulk disease-context dataset only. |",
    "| Table S7 | GSE73680 P1 single-gene response | `results/gse73680/tables/gse73680_p1_gene_response.tsv`; `results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv` | Shows P1 genes are heterogeneous and not FDR-supported. | PKD2 nominal only. |",
    "| Table S8 | GSE73680 module response and robustness checks | `results/gse73680/tables/gse73680_patient_level_module_response.tsv`; `gse73680_module_leave_one_gene_out.tsv`; `gse73680_magma_without_p1_sensitivity.tsv`; `gse73680_random_module_benchmark.tsv`; `gse73680_expression_matched_random_benchmark.tsv` | Supports MAGMA module-level disease-context association and documents sensitivity checks. | Module-level support only. |",
    "| Table S9 | Claim boundary and resource-limited analyses | `results/tables/claim_boundary_audit_v0.1.tsv`; `results/tables/integrated_evidence_summary_v0.2.tsv`; `docs/twas_resource_acquisition_log.md`; `docs/coloc_run_plan.md`; `docs/phase4_spatial_prep_status.md` | Documents boundaries and pending extensions. | Resource-limited is not negative evidence. |"
  ), "docs/supplementary_table_plan_v0.1.md")

  write_lines(c(
    "# Supplementary Figure Plan v0.1",
    "",
    "| Figure | Title | Candidate files | Purpose | Claim boundary |",
    "|---|---|---|---|---|",
    "| Figure S1 | GWAS QQ and Manhattan plots | `results/figures/phase1_gwas_2025_manhattan_plot.png`; `phase1_gwas_2025_qq_plot.pdf`; `phase1_gwas_qq_no_gws.pdf`; `phase1_gwas_qq_without_known_loci.pdf` | Shows GWAS-level QC visuals. | QC display only. |",
    "| Figure S2 | GSE231569 UMAP and marker audit | `results/figures/gse231569_umap_audited_broad_celltype.pdf`; `gse231569_marker_dotplot_by_cluster.pdf`; `gse231569_transport_marker_audit.pdf` | Supports annotation audit. | Not causal cell-type proof. |",
    "| Figure S3 | Random benchmark for scRNA localization | `results/figures/magma_scrna_benchmark.pdf`; `magma_leave_one_locus_out_tal.pdf`; `audited_locus_scrna_benchmark.pdf` | Shows TAL-associated localization and robustness. | Not spatial validation. |",
    "| Figure S4 | P1 extended donor-level plots | `results/figures/p1_tal_gene_dotplot.pdf`; `p1_tal_gene_by_donor_boxplot.pdf`; `p1_gene_vs_tal_program_correlation.pdf`; `p1_tal_gene_featureplots.pdf` | Provides extended P1 evidence. | Not disease validation. |",
    "| Figure S5 | GSE73680 QC and PCA | `results/gse73680/figures/gse73680_pca.pdf`; `gse73680_density_plot.pdf`; `gse73680_sample_correlation_heatmap.pdf`; `gse73680_boxplot.pdf` | Shows bulk microarray QC. | Bulk dataset only. |",
    "| Figure S6 | GSE73680 robustness checks | `results/gse73680/figures/gse73680_module_leave_one_gene_out.pdf`; `gse73680_random_module_benchmark.pdf`; `gse73680_expression_matched_random_benchmark.pdf`; `gse73680_magma_module_paired_delta_spaghetti.pdf` | Documents module-level robustness and conservative benchmark. | Expression-matched benchmark is a boundary check. |",
    "| Figure S7 | Resource-limited TWAS/coloc/spatial status | Resource manifest tables rendered as a simple status figure if needed | Shows prepared but unavailable extensions. | Not negative TWAS/coloc/spatial evidence. |"
  ), "docs/supplementary_figure_plan_v0.1.md")
}

internal_review_checklist <- function() {
  write_lines(c(
    "# Internal Review Checklist v0.1",
    "",
    "## Claim Boundary",
    "",
    "- [ ] GSE73680 is never described as P1 gene validation.",
    "- [ ] MAGMA is never described as causal proof.",
    "- [ ] TAL localization is described as cellular context, not causal mediation.",
    "- [ ] TWAS/SMR-coloc/spatial modules are described as resource-limited and pending, not negative.",
    "- [ ] GSE73680 is described as bulk/microarray disease context, not cell-type-specific disease validation.",
    "",
    "## Figure and Table Consistency",
    "",
    "- [ ] Figure 1 uses four evidence layers and a separate claim boundary box.",
    "- [ ] Figure 2 uses localization/context language and does not imply causality.",
    "- [ ] Figure 3 states P1 evidence spectrum, not disease validation.",
    "- [ ] Figure 4 states MAGMA module-level disease-context support and no P1 gene q≤0.05.",
    "- [ ] P values and FDR values match source tables.",
    "- [ ] Gene-set names are consistent: MAGMA top 50, MAGMA top 100, MAGMA FDR, MAGMA suggestive, P1 core candidates, TAL markers, injury/remodeling.",
    "",
    "## Methods Reproducibility",
    "",
    "- [ ] GWAS QC and MAGMA inputs are described with reference build and software version.",
    "- [ ] GSE231569 annotation audit and gene-set projection are described separately.",
    "- [ ] GSE73680 reconstruction states 62 TXT.gz files, 44,661 features, 32,055 mapped genes, 55 included samples, 29 patients and 26 paired patients.",
    "- [ ] Primary GSE73680 model is limma duplicateCorrelation.",
    "- [ ] Paired patient-level delta is described as sensitivity analysis.",
    "- [ ] Size-matched benchmark is main Figure 4 benchmark.",
    "- [ ] Expression-matched benchmark is described as conservative sensitivity analysis.",
    "",
    "## Current Status",
    "",
    "Phase 7 package is internally consistent if all boxes above remain checked during final prose polishing."
  ), "docs/internal_review_checklist_v0.1.md")
}

manuscript_v03 <- function() {
  write_lines(c(
    "# Integrative post-GWAS mapping links kidney stone risk genes to TAL-associated renal papillary cellular and disease contexts",
    "",
    "## Abstract",
    "",
    "**Background:** Kidney stone disease has a substantial genetic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.",
    "",
    "**Methods:** We built an integrative post-GWAS framework combining locked kidney stone GWAS summary statistics, MAGMA gene-based prioritization, audited GSE231569 renal papillary single-nucleus annotations, a six-gene P1 candidate evidence spectrum and GSE73680 papillary plaque/stone disease-context expression analysis. GSE73680 analyses used patient-aware limma duplicateCorrelation models, paired patient-level sensitivity analyses and random gene-set benchmarks.",
    "",
    "**Results:** MAGMA-prioritized genes converged on a Loop of Henle/thick ascending limb-associated single-nucleus context in GSE231569. The P1 candidate genes formed an interpretable TAL, epithelial transport and calcium-handling expression spectrum rather than a uniform TAL marker panel. In GSE73680, MAGMA-prioritized modules showed disease-context expression shifts in plaque/stone papilla samples, whereas P1 single-gene responses were heterogeneous and no P1 gene reached FDR significance.",
    "",
    "**Conclusion:** These results support a TAL-associated cellular context and MAGMA module-level disease-context association for KSD genetic risk, while not establishing causal mediation, TWAS convergence, colocalization or spatial validation.",
    "",
    "## Introduction",
    "",
    "Kidney stone disease is a common and recurrent disorder with contributions from urinary chemistry, renal epithelial transport, papillary microenvironments and inherited susceptibility. Genome-wide association studies have identified multiple loci associated with kidney stone risk, but most association signals do not immediately nominate a causal gene, cell type or disease-stage context. Post-GWAS interpretation therefore requires a framework that connects statistical genetic prioritization to renal cell states and disease-relevant expression patterns while preserving clear boundaries around causal inference.",
    "",
    "The renal papilla is a plausible tissue context for stone formation because it contains epithelial, interstitial, vascular and immune compartments involved in urine concentration, mineral handling and papillary plaque biology. Single-cell and single-nucleus transcriptomic resources provide an opportunity to localize genetic signals to renal cell states, but such analyses depend on careful annotation audit and should not be interpreted as causal mediation on their own. Similarly, disease-context expression datasets can provide support for expression shifts in plaque or stone-associated tissue, but they cannot substitute for colocalization, TWAS or spatial validation.",
    "",
    "Here we integrated MAGMA gene-level prioritization from locked kidney stone GWAS summary statistics with audited GSE231569 renal papillary single-nucleus annotations and GSE73680 papillary plaque/stone expression data. The study was designed to ask whether kidney stone genetic risk genes converge on a renal papillary cellular context, whether a core set of candidate genes forms a coherent expression spectrum, and whether an independent disease-context dataset supports the prioritized modules. We explicitly distinguish MAGMA plus scRNA-supported cellular context from genetic causality, TWAS convergence, colocalization and spatial validation.",
    "",
    "## Methods",
    "",
    "### GWAS summary statistics and QC",
    "",
    "Locked kidney stone disease GWAS summary statistics were processed through the fixed Phase 1 workflow. Input columns, sample-size fields, lead-locus records and summary-statistic sanity checks were audited before downstream gene-level analysis. The primary GWAS and lead-locus definitions were not redefined during downstream single-nucleus or disease-context interpretation.",
    "",
    "### MAGMA gene-based prioritization",
    "",
    "MAGMA v1.10 was used for gene-based prioritization with a GRCh37/hg19-compatible reference. Gene-level results were ranked by MAGMA P value, and downstream gene sets were defined from fixed ranking and significance thresholds, including MAGMA top 50, MAGMA top 100, FDR-significant and suggestive gene sets. These MAGMA-derived sets formed the main genetic prioritization layer for all later analyses.",
    "",
    "### GSE231569 single-nucleus processing and annotation audit",
    "",
    "Renal papillary single-nucleus data from GSE231569 were processed with an audited broad cell-type annotation. Marker expression, cluster assignment and renal epithelial compartment labels were reviewed before gene-set projection. The analysis focused on stable audited cell types, including the Loop of Henle/thick ascending limb compartment, and avoided extending annotation-derived observations into causal claims.",
    "",
    "### Gene-set projection to single-nucleus cell types",
    "",
    "MAGMA-prioritized gene sets were projected onto audited GSE231569 renal papillary cell types using module-score summaries. Random gene-set benchmarks were used to compare observed cell-type localization against background expectations. Locus-balanced and leave-one-locus-out analyses evaluated whether TAL-associated localization was robust to locus composition. These analyses were interpreted as single-nucleus cellular context mapping rather than causal mediation or spatial validation.",
    "",
    "### P1 candidate gene evidence scoring",
    "",
    "Six P1 candidate genes, UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2, were evaluated across audited GSE231569 cell types. For each gene, the evidence summary included TAL expression rank, average expression, detection frequency, donor-level detection, TAL specificity ratio, TAL program correlation and manuscript role assignment. The purpose was to define a gene-centric evidence spectrum across TAL, epithelial transport and calcium-handling biology, not to force all P1 genes into a uniform TAL marker class.",
    "",
    "### GSE73680 reconstruction from Agilent feature files",
    "",
    "GSE73680 was reconstructed from the local `GSE73680_RAW.tar` archive. The archive yielded 62 TXT.gz files, all of which passed gzip validation and Agilent text-structure checks. Agilent FEATURES blocks were parsed into a 44,661-feature by 62-sample expression matrix. Direct gene-symbol mapping assigned 32,055 mapped genes, and curated metadata linked expression profiles to patient and tissue-context labels. After inclusion filtering, the disease-context analysis used 55 samples from 29 patients, including 26 paired patients contributing both control/adjacent and plaque/stone papilla samples.",
    "",
    "### Patient-aware disease-context analysis",
    "",
    "Primary GSE73680 disease-context analyses used limma duplicateCorrelation to account for repeated patient-level sampling. Paired patient-level delta analyses were used as sensitivity checks because 26 patients contributed both analysis groups. P1 single-gene responses were analyzed separately from module-level responses. Module scores were computed as mean z-scores across detected genes for P1 core candidates, MAGMA top 50, MAGMA top 100, MAGMA FDR, MAGMA suggestive, TAL marker and injury/remodeling marker sets.",
    "",
    "### Random gene-set benchmarking",
    "",
    "Size-matched random gene-set benchmarking was used as the main Figure 4 benchmark to evaluate whether observed disease-context module shifts exceeded random expectations for gene sets of comparable size. Leave-one-gene-out analyses, MAGMA-without-P1 sensitivity and paired direction-consistency checks were used to evaluate robustness. Expression-matched random benchmarking was retained as a conservative sensitivity analysis; under the current implementation, it did not provide primary support and was not used as the main benchmark.",
    "",
    "### Resource-limited TWAS, SMR/coloc and spatial extensions",
    "",
    "TWAS, SMR/coloc and spatial transcriptomic extensions were prepared through input manifests and resource checks but were not completed as evidence layers because required external expression weights, matched eQTL resources or spatial matrix/image/coordinate files were unavailable in the current analysis environment. These resource-limited extensions were therefore not used to support manuscript claims and should not be interpreted as negative results.",
    "",
    "## Results",
    "",
    "### GWAS and MAGMA prioritization identify KSD-associated genes",
    "",
    "The locked kidney stone GWAS reconstruction retained 57 lead-locus records and was carried forward into MAGMA gene-level prioritization. MAGMA v1.10 tested 17,316 genes using a GRCh37/hg19-compatible reference and identified 94 Bonferroni-significant genes, 369 FDR-significant genes and 187 suggestive genes at P < 1e-4. These outputs defined the genetic prioritization layer for downstream single-nucleus and disease-context analyses.",
    "",
    "### MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context",
    "",
    "MAGMA-prioritized gene sets localized most strongly to the audited Loop of Henle/TAL compartment in GSE231569. TAL benchmark percentiles were 0.998 for the MAGMA top 50 set, 1.000 for the top 100 set, 0.999 for the top 200 set, 1.000 for suggestive genes and 0.968 for FDR-significant genes. Locus-balanced and leave-one-locus-out analyses supported robustness of the top 50 TAL signal. These results support a TAL-associated single-nucleus cellular context for MAGMA-prioritized kidney stone risk genes.",
    "",
    "### P1 candidates form an interpretable TAL/epithelial transport/calcium-handling expression spectrum",
    "",
    "The six P1 candidate genes did not behave as a uniform TAL-specific marker set. UMOD served as a representative TAL-associated gene, CLDN10 supported an epithelial transport pattern, CLDN14 and CASR linked the signal to ion-handling and calcium-sensing biology, HIBADH remained a supporting MAGMA-associated context gene and PKD2 represented a broader renal epithelial context. Figure 3 summarizes this evidence across expression, specificity, donor-level descriptive consistency and gene-role assignment. The supported claim is a MAGMA plus scRNA-based TAL-associated cellular context, not causal mediation or colocalized genetic validation.",
    "",
    "### GSE73680 supports MAGMA module-level disease-context expression association",
    "",
    "We next examined whether the MAGMA- and single-nucleus-prioritized signals were reflected in an independent papillary plaque/stone disease-context dataset. After reconstructing GSE73680 from Agilent supplementary feature files, 55 samples from 29 patients were included, including 26 paired patients with control/adjacent and plaque/stone papilla samples. To account for non-independence among repeated samples, we used a patient-aware limma duplicateCorrelation model and paired patient-level sensitivity analyses.",
    "",
    "At the single-gene level, none of the six P1 candidate genes reached FDR-supported differential expression, although PKD2 showed a nominal exploratory paired response. In contrast, MAGMA-prioritized modules showed consistent paired disease-context shifts. Patient-level MAGMA module responses reached q ≤ 0.05, retained directionality in leave-one-gene-out analyses, remained robust after removing the six P1 genes, and exceeded size-matched random gene-set expectations. Paired direction analyses showed positive shifts in approximately 69% to 73% of paired patients for MAGMA modules. A stricter expression-matched benchmark did not provide primary support under the current implementation and is treated as a conservative boundary check. These findings indicate that GSE73680 supports disease-context expression association at the MAGMA module level rather than uniform single-gene P1 differential expression.",
    "",
    "## Discussion",
    "",
    "The main finding of this study is not that a single KSD gene or a uniform P1 gene panel explains papillary disease biology. Instead, MAGMA-prioritized KSD genes converged on a TAL-associated cellular context and showed module-level disease-context expression association in GSE73680. This distinction is central to the interpretation of the study: the evidence supports a post-GWAS cellular and disease-context mapping model, not causal validation of a single gene or pathway.",
    "",
    "The distinction between single-gene and module-level evidence is important. P1 genes remain MAGMA and scRNA-supported candidates, but GSE73680 did not validate them as uniformly differentially expressed disease genes. Instead, the GSE73680 signal was stronger for MAGMA-prioritized modules. This result is consistent with a polygenic and pathway-level post-GWAS interpretation in which many modestly shifted genes collectively define disease-context expression support.",
    "",
    "These analyses also illustrate the value of explicit claim boundaries. Single-nucleus localization does not prove cell-type-mediated genetic causality, disease-context expression does not prove causal involvement in plaque formation, and random gene-set benchmarking does not substitute for TWAS or colocalization. The current evidence supports a post-GWAS cellular and disease-context expression model that can motivate future causal and spatial validation.",
    "",
    "## Limitations",
    "",
    "Several limitations should guide interpretation. TWAS, SMR/coloc and spatial transcriptomic validation were not completed because required external expression-weight, eQTL and spatial matrix/image resources were unavailable in the current analysis environment. These resource-limited modules should not be interpreted as negative results. GSE73680 supports disease-context module association but not causality. P1 single-gene disease differential expression was not FDR-supported, and PKD2 should be treated only as a nominal exploratory observation. Because GSE73680 is a bulk microarray disease-context dataset, cell-type-specific disease expression cannot be inferred from it.",
    "",
    "## Data availability",
    "",
    "All analyses were performed using public GWAS, GEO single-nucleus and GEO expression resources reconstructed into local analysis-ready files. Intermediate tables, scripts, QC summaries and figure outputs are stored within the project workspace. External resource limitations for TWAS, SMR/coloc and spatial analyses are documented in the corresponding resource status notes.",
    "",
    "## Code availability",
    "",
    "Analysis scripts are organized under `scripts/`, with GSE73680 plaque-context workflows under `scripts/08_plaque_context/` and manuscript figure scripts under `scripts/09_manuscript/`. Reproducible output tables are stored under `results/tables/` and `results/gse73680/tables/`.",
    "",
    "## Figure legends",
    "",
    "### Figure 1. Post-GWAS framework for KSD cellular and disease-context mapping",
    "",
    "The study is organized as four evidence layers: GWAS/MAGMA gene prioritization, audited GSE231569 single-nucleus localization, P1 candidate gene evidence and GSE73680 disease-context module analysis. TWAS, SMR/coloc and spatial transcriptomic analyses were prepared as resource-limited extensions but were not used for manuscript claims.",
    "",
    "### Figure 2. MAGMA-prioritized KSD genes localize to TAL-associated cell states",
    "",
    "MAGMA-prioritized gene sets were evaluated across audited renal papillary single-nucleus cell types from GSE231569. TAL-associated enrichment and robustness analyses supported the Loop of Henle/TAL compartment as the strongest cellular context for prioritized kidney stone risk genes.",
    "",
    "### Figure 3. Gene-centric single-nucleus evidence for P1 core TAL-associated KSD candidates",
    "",
    "P1 candidate genes were evaluated across audited renal papillary cell types, TAL specificity, donor-level descriptive consistency and gene-role assignments. These analyses support a TAL/transport/calcium-handling evidence spectrum but do not establish causal mediation, TWAS convergence, colocalization or spatial validation.",
    "",
    "### Figure 4. GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression",
    "",
    "Patient-aware GSE73680 analyses showed MAGMA module-level disease-context expression association in plaque/stone papilla samples. P1 single-gene responses were heterogeneous and no P1 gene reached q≤0.05 after FDR correction. Size-matched random benchmarking supported MAGMA module shifts, whereas expression-matched benchmarking was retained as a conservative sensitivity analysis rather than primary support. The analysis supports MAGMA-prioritized module-level disease-context association and does not establish genetic causality, TWAS convergence, colocalization or spatial validation."
  ), "manuscript/manuscript_draft_v0.3.md")
}

copy_final_draft <- function() {
  candidates <- c(
    "results/figures/figure1_integrative_framework_v0.2.pdf",
    "results/figures/figure1_integrative_framework_v0.2.png",
    "results/figures/magma_scrna_benchmark.pdf",
    "results/figures/magma_leave_one_locus_out_tal.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.png",
    "results/figures/figure4_gse73680_disease_context_v0.8.pdf",
    "results/figures/figure4_gse73680_disease_context_v0.8.png"
  )
  file.copy(candidates[file.exists(candidates)], "results/figures/final_draft", overwrite = TRUE)
  write_lines(c(
    "# Final Draft Figure Package v0.2",
    "",
    "Included files:",
    "- figure1_integrative_framework_v0.2.pdf/png",
    "- magma_scrna_benchmark.pdf",
    "- magma_leave_one_locus_out_tal.pdf",
    "- figure3_p1_gene_evidence_v0.4.pdf/png",
    "- figure4_gse73680_disease_context_v0.8.pdf/png",
    "",
    "Current status:",
    "- Figure 1 v0.2 and Figure 4 v0.8 are the preferred current versions.",
    "- Figure 3 v0.4 remains frozen.",
    "- Figure 2 is represented by MAGMA scRNA localization and robustness panels; a later layout pass can combine them if needed."
  ), "results/figures/final_draft/README.md")
}

figure1_v02()
figure4_v08()
figure_qc_outputs()
supplementary_plans()
internal_review_checklist()
manuscript_v03()
copy_final_draft()

message("wrote Phase 7 submission package outputs")
