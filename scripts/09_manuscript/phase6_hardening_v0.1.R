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
  scrna = "#6B8E6E",
  p1 = "#B59A5B",
  disease = "#9A5F52",
  boundary = "#555555",
  resource = "#8E8E8E",
  muted = "#D8D8D8",
  text = "#303030",
  panel = "#F7F8F8"
)

write_lines <- function(x, path) {
  writeLines(x, con = path, useBytes = TRUE)
}

figure4_v07 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  modules <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
  rand_bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

  panel_a <- ggplot() +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.78, ymax = 0.94, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.54, ymax = 0.70, fill = pal$panel, color = "#8A8A8A", linewidth = 0.35) +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.30, ymax = 0.46, fill = pal$panel, color = "#8A8A8A", linewidth = 0.35) +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.06, ymax = 0.22, fill = "#EEF3F4", color = pal$genetic, linewidth = 0.4) +
    annotate("text", x = 0.50, y = 0.86, label = "GSE73680 disease-context dataset", fontface = "bold", size = 3.6, color = pal$text) +
    annotate("text", x = 0.50, y = 0.62, label = "55 included samples: 27 control/adjacent, 28 plaque/stone papilla", size = 3.25, color = pal$text) +
    annotate("text", x = 0.50, y = 0.38, label = "29 patients, including 26 paired patients", size = 3.25, color = pal$text) +
    annotate("text", x = 0.50, y = 0.16, label = "Primary: limma duplicateCorrelation", fontface = "bold", size = 3.2, color = pal$text) +
    annotate("text", x = 0.50, y = 0.10, label = "Sensitivity: paired patient-level delta", fontface = "bold", size = 3.2, color = pal$text) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.77, yend = 0.71, arrow = arrow(length = unit(0.13, "in")), color = "#777777") +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.53, yend = 0.47, arrow = arrow(length = unit(0.13, "in")), color = "#777777") +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.29, yend = 0.23, arrow = arrow(length = unit(0.13, "in")), color = "#777777") +
    labs(title = "A. Patient-aware GSE73680 design") +
    theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 11),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No detectable")]
  p1[, p_label := sprintf("P=%.3f", p_value)]
  panel_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
    geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.45, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$p1, "No detectable" = pal$muted)) +
    labs(x = "Patient-level paired delta", y = NULL, fill = "Signal",
         title = "B. No uniform P1 single-gene response") +
    theme_bw(base_size = 9.5) +
    theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(face = "bold"),
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
    scale_fill_manual(values = c("q <= 0.05" = pal$genetic, Borderline = pal$p1,
                                 "No detectable" = pal$muted, "Disease-context program" = pal$scrna)) +
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
  ggsave("results/figures/figure4_gse73680_disease_context_v0.7.pdf", fig, width = 12.5, height = 8.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v0.7.png", fig, width = 12.5, height = 8.8, units = "in", dpi = 240, bg = "white")

  write_lines(c(
    "# Figure 4 Legend v0.7",
    "",
    "**Figure 4. GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression.**",
    "(A) Patient-aware analysis design for GSE73680. A total of 55 included samples, including 27 control/adjacent and 28 plaque/stone papilla samples, were analyzed from 29 patients, including 26 paired patients. Primary analyses used limma duplicateCorrelation, with paired patient-level delta analyses as sensitivity checks.",
    "(B) Paired patient-level responses of six P1 candidate genes. PKD2 showed a nominal exploratory response, but no P1 gene reached FDR significance.",
    "(C) Paired patient-level module responses. MAGMA-prioritized modules showed disease-context shifts, whereas the P1 core candidate module and TAL marker set did not show comparable support. The injury/remodeling marker set served as a disease-context reference program, not as a MAGMA-prioritized module.",
    "(D) Size-matched random gene-set benchmark showing that MAGMA module shifts exceeded random expectation, whereas the P1 core candidate module was background-like. The dashed line indicates the 95th percentile of size-matched random gene-set expectations.",
    "Together, these results support disease-context expression association for MAGMA-prioritized modules, rather than uniform P1 single-gene differential expression. A stricter expression-matched benchmark did not provide primary support under the current implementation and is retained as a conservative boundary check. These findings do not establish genetic causality, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure4_legend_v0.7.md")
}

figure1_v01 <- function() {
  steps <- data.table(
    x = c(0.10, 0.27, 0.44, 0.61, 0.78, 0.93),
    y = rep(0.58, 6),
    label = c(
      "KSD GWAS\nsummary statistics",
      "MAGMA gene\nprioritization",
      "Audited GSE231569\nsingle-nucleus localization",
      "P1 candidate gene\nevidence spectrum",
      "GSE73680 disease-context\nmodule analysis",
      "Claim boundary:\ncellular context +\nmodule-level association"
    ),
    fill = c(pal$genetic, pal$genetic, pal$scrna, pal$p1, pal$disease, pal$boundary)
  )
  fig <- ggplot() +
    annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.18, ymax = 0.86, fill = "white", color = NA) +
    geom_segment(data = steps[1:5], aes(x = x + 0.065, xend = steps$x[2:6] - 0.065, y = y, yend = y),
                 arrow = arrow(length = unit(0.11, "in")), color = "#6F6F6F", linewidth = 0.45) +
    geom_label(data = steps, aes(x = x, y = y, label = label, fill = fill),
               color = "white", fontface = "bold", size = 3.1, linewidth = 0,
               label.padding = unit(0.20, "lines"), lineheight = 0.92) +
    annotate("rect", xmin = 0.70, xmax = 0.98, ymin = 0.02, ymax = 0.25, fill = "#F5F5F5", color = "#9A9A9A", linewidth = 0.35) +
    annotate("text", x = 0.84, y = 0.205, label = "Resource-limited extensions", fontface = "bold", size = 3.6, color = pal$text) +
    annotate("text", x = 0.84, y = 0.135, label = "TWAS   |   SMR/coloc   |   spatial transcriptomics", size = 3.25, color = pal$text) +
    annotate("text", x = 0.50, y = 0.93, label = "Integrative post-GWAS framework for mapping KSD genetic risk to renal papillary cellular and disease contexts",
             fontface = "bold", size = 4.2, color = pal$text) +
    annotate("text", x = 0.50, y = 0.33, label = "Main supported inference: MAGMA-prioritized KSD genes converge on TAL-associated cellular context and show module-level disease-context support.",
             size = 3.25, color = pal$text) +
    annotate("text", x = 0.50, y = 0.285, label = "Not established: causal mediation, P1 disease-gene validation, TWAS convergence, colocalization or spatial validation.",
             size = 3.1, color = pal$boundary) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  ggsave("results/figures/figure1_integrative_framework_v0.1.pdf", fig, width = 13.2, height = 4.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.1.png", fig, width = 13.2, height = 4.8, units = "in", dpi = 260, bg = "white")

  write_lines(c(
    "# Figure 1 Legend v0.1",
    "",
    "**Figure 1. Integrative post-GWAS framework for mapping KSD genetic risk to renal papillary cellular and disease contexts.**",
    "The study links locked kidney stone disease GWAS summary statistics to MAGMA gene-level prioritization, audited GSE231569 renal papillary single-nucleus localization, a P1 candidate gene evidence spectrum and GSE73680 disease-context module analysis.",
    "The supported inference is cellular and disease-context mapping: MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context and show module-level disease-context expression association. TWAS, SMR/coloc and spatial transcriptomic analyses are shown as resource-limited extensions and are not treated as completed evidence layers.",
    "The framework explicitly does not establish causal mediation, P1 disease-gene validation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure1_legend_v0.1.md")
}

claim_boundary_outputs <- function() {
  audit <- data.table(
    claim = c(
      "MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context.",
      "P1 genes form an interpretable TAL/transport/calcium expression spectrum.",
      "GSE73680 supports MAGMA module-level disease-context expression association.",
      "P1 genes are individually validated disease genes.",
      "TAL causally mediates KSD genetic risk.",
      "TWAS or SMR/coloc confirms expression-mediated genetic risk.",
      "Spatial transcriptomics validates TAL localization."
    ),
    supported_by = c(
      "Locked GWAS/MAGMA outputs; audited GSE231569 scRNA projection; random and locus-balanced benchmarks.",
      "GSE231569 P1 dotplot, specificity, donor detection, TAL program correlation and role assignment.",
      "GSE73680 patient-aware module analysis; paired sensitivity; size-matched benchmark; robustness checks.",
      "Not supported by current GSE73680 single-gene FDR results.",
      "Not supported by current observational post-GWAS and scRNA mapping analyses.",
      "Not completed because required external resources were unavailable.",
      "Not completed because required spatial matrix/image/coordinate resources were unavailable."
    ),
    support_strength = c("strong", "moderate", "moderate", "not_supported", "not_supported", "resource_limited", "resource_limited"),
    allowed_wording = c(
      "MAGMA-prioritized KSD genes localize to or converge on a TAL-associated single-nucleus context.",
      "P1 genes define an interpretable TAL/epithelial transport/calcium-handling expression spectrum.",
      "GSE73680 provides MAGMA module-level disease-context expression support.",
      "P1 genes remain MAGMA plus scRNA-supported candidates; GSE73680 single-gene responses were heterogeneous and not FDR-supported.",
      "The results nominate TAL-associated cellular context for future mechanistic testing.",
      "TWAS and SMR/coloc remain resource-limited pending extensions, not negative evidence.",
      "Spatial transcriptomics remains a resource-limited pending extension, not negative evidence."
    ),
    forbidden_wording = c(
      "TAL causally mediates KSD risk; spatially validated TAL mechanism.",
      "All P1 genes are uniform TAL markers or disease-validated genes.",
      "GSE73680 validates P1 genes or proves plaque causality.",
      "P1 genes are validated disease genes; PKD2 is validated by GSE73680.",
      "TAL causally mediates genetic risk; TAL is proven as causal cell type.",
      "TWAS confirms convergence; coloc confirms shared causal variants.",
      "Spatial data validates TAL localization or plaque microenvironment localization."
    ),
    figure = c("Figure 2", "Figure 3", "Figure 4", "Figure 4", "Figures 1-4", "Limitations", "Limitations"),
    result_section = c("Result 2", "Result 3", "Result 4", "Result 4", "Discussion/Limitations", "Limitations", "Limitations"),
    notes = c(
      "Use localization/context language, not mediation language.",
      "Heterogeneity is part of the claim.",
      "This is the main disease-context claim.",
      "PKD2 is nominal exploratory only; no P1 gene reached FDR support.",
      "Requires future causal, perturbational or colocalized regulatory evidence.",
      "Unavailable resources should not be framed as failed analyses.",
      "Unavailable resources should not be framed as failed analyses."
    )
  )
  fwrite(audit, "results/tables/claim_boundary_audit_v0.1.tsv", sep = "\t")
  md <- c(
    "# Claim Boundary Audit v0.1",
    "",
    "## Purpose",
    "",
    "This audit fixes the wording boundary for the manuscript after Phase 6 integration. The manuscript should present the project as post-GWAS cellular and disease-context mapping, not as causal validation.",
    "",
    "## Allowed Core Claims",
    "",
    "- MAGMA-prioritized KSD genes converge on a TAL-associated single-nucleus context.",
    "- P1 genes form an interpretable TAL/epithelial transport/calcium-handling expression spectrum.",
    "- GSE73680 supports MAGMA module-level disease-context expression association.",
    "",
    "## Forbidden Claims",
    "",
    "- TAL causally mediates KSD genetic risk.",
    "- P1 genes are validated disease genes.",
    "- GSE73680 validates P1 genes.",
    "- TWAS confirms expression-mediated risk.",
    "- SMR/coloc confirms shared causal variants.",
    "- Spatial transcriptomics validates TAL localization.",
    "",
    "## Practical Wording Rule",
    "",
    "Use **supports**, **localizes**, **prioritizes**, **is consistent with** and **disease-context expression association**. Avoid **proves**, **validates**, **causes**, **mediates**, **confirms TWAS**, **confirms colocalization** and **spatially validates** unless future analyses directly provide those results.",
    "",
    "## Audit Table",
    "",
    "See `results/tables/claim_boundary_audit_v0.1.tsv`."
  )
  write_lines(md, "docs/claim_boundary_audit_v0.1.md")
}

integration_tables_v02 <- function() {
  integrated <- data.table(
    evidence_layer = c("GWAS/MAGMA", "scRNA localization", "P1 genes", "Disease context", "TWAS", "SMR/coloc", "spatial"),
    dataset = c("2025 KSD GWAS", "GSE231569", "GSE231569", "GSE73680", "External expression-weight resources", "External eQTL resources", "Spatial transcriptomics resources"),
    analysis_method = c(
      "GWAS QC and MAGMA gene-level analysis",
      "Audited annotation and gene-set projection with random/locus-balanced benchmarks",
      "Dotplot, specificity, donor detection, TAL program correlation and role assignment",
      "Patient-aware limma duplicateCorrelation, paired patient-level delta, module benchmarks and robustness checks",
      "Resource manifest and harmonization checks; not run as evidence layer",
      "Resource manifest and pilot-readiness checks; not run as evidence layer",
      "Input inventory/resource checks; not run as evidence layer"
    ),
    main_result = c(
      "KSD-associated genes were prioritized at gene and module level.",
      "MAGMA-prioritized genes localized to a TAL-associated renal papillary single-nucleus context.",
      "Six P1 genes formed an interpretable TAL/transport/calcium-handling expression spectrum.",
      "MAGMA-prioritized modules shifted in plaque/stone papilla disease context; P1 single genes were not uniformly FDR-supported.",
      "Required prediction weights were not sufficiently available for completed TWAS.",
      "Required matched eQTL/summary resource support was not sufficient for completed SMR/coloc.",
      "Required matrix/image/coordinate resources were not sufficient for completed spatial validation."
    ),
    support_strength = c("strong", "strong", "moderate", "moderate", "resource_limited", "resource_limited", "resource_limited"),
    supports_main_claim = c(
      "genetic gene-level prioritization",
      "TAL-associated cellular context",
      "gene-centric TAL/transport/calcium context",
      "MAGMA module-level disease-context association",
      "none; pending extension only",
      "none; pending extension only",
      "none; pending extension only"
    ),
    does_not_support = c(
      "causality not established",
      "not spatial validation or causal mediation",
      "not uniform TAL specificity or P1 disease validation",
      "not uniform P1 single-gene differential expression or causality",
      "not negative TWAS evidence",
      "not negative colocalization evidence",
      "not negative spatial evidence"
    ),
    figure = c("Figure 1/2", "Figure 2", "Figure 3", "Figure 4", "Limitations", "Limitations", "Limitations"),
    table = c(
      "results/tables/magma_qc_summary.tsv; results/tables/magma_genes.tsv",
      "results/tables/magma_scrna_evidence_summary.tsv",
      "results/tables/p1_tal_gene_interpretation_summary.tsv",
      "results/gse73680/tables/gse73680_disease_context_summary_v0.2.tsv",
      "results/tables/twas_resource_manifest_v2.tsv",
      "results/tables/eqtl_resource_manifest_v0.2.tsv",
      "results/tables/spatial_input_inventory_v0.1.tsv"
    ),
    manuscript_use = c("main", "main", "main", "main", "limitations", "limitations", "limitations"),
    claim_boundary = c(
      "Gene prioritization only; do not claim causal genes.",
      "Cellular localization/context only; do not claim causal TAL mediation.",
      "Candidate evidence spectrum only; do not claim P1 gene validation.",
      "Module-level disease-context support only; do not claim P1 gene validation or causality.",
      "Resource-limited pending extension; do not write as negative result.",
      "Resource-limited pending extension; do not write as negative result.",
      "Resource-limited pending extension; do not write as negative result."
    )
  )
  fwrite(integrated, "results/tables/integrated_evidence_summary_v0.2.tsv", sep = "\t")
}

candidate_tiers_v04 <- function() {
  magma <- fread("results/tables/magma_genes.tsv")
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  single <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")

  gene_sets <- magma[, .(
    gene = gene_symbol,
    MAGMA_rank = rank,
    MAGMA_p = p,
    MAGMA_fdr = fdr,
    MAGMA_gene_set = fifelse(rank <= 50, "MAGMA_top50",
                      fifelse(rank <= 100, "MAGMA_top100",
                      fifelse(fdr <= 0.05, "MAGMA_FDR",
                      fifelse(suggestive == TRUE, "MAGMA_suggestive", "MAGMA_ranked"))))
  )]
  gene_sets <- gene_sets[MAGMA_rank <= 100 | MAGMA_fdr <= 0.05 | MAGMA_gene_set == "MAGMA_suggestive" | gene %in% p1$gene]

  out <- merge(gene_sets, p1[, .(gene, scRNA_TAL_context = overall_evidence_class, P1_role = manuscript_role,
                                 P1_interpretation = final_interpretation)], by = "gene", all.x = TRUE)
  out <- merge(out, single[, .(gene, GSE73680_single_gene_response = interpretation,
                               GSE73680_single_gene_p = p_value,
                               GSE73680_single_gene_fdr = fdr)], by = "gene", all.x = TRUE)
  out[is.na(scRNA_TAL_context), scRNA_TAL_context := "not assessed in P1 gene-centric scRNA evidence pack"]
  out[is.na(P1_role), P1_role := "not a P1 core candidate"]
  out[is.na(GSE73680_single_gene_response), GSE73680_single_gene_response := "not tested as P1 single-gene endpoint"]
  out[, GSE73680_module_context := fifelse(
    MAGMA_gene_set %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive"),
    "module-level support only for MAGMA-prioritized gene sets",
    "not part of the supported GSE73680 module claim"
  )]
  out[, TWAS_status := "resource_limited_not_completed"]
  out[, coloc_status := "resource_limited_not_completed"]
  out[, final_tier_v0.4 := fcase(
    gene %in% c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"),
    "Tier 1: MAGMA + scRNA-supported P1 candidates",
    MAGMA_gene_set %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive"),
    "Tier 2: MAGMA module-level disease-context supported genes",
    default = "Tier 3: context candidates"
  )]
  out[, interpretation := fcase(
    final_tier_v0.4 == "Tier 1: MAGMA + scRNA-supported P1 candidates",
    "P1 core candidate supported by MAGMA and scRNA context; GSE73680 support is module-level and does not validate this gene individually.",
    final_tier_v0.4 == "Tier 2: MAGMA module-level disease-context supported genes",
    "MAGMA-prioritized gene within modules supported at disease-context level; gene-specific causality is not established.",
    default = "Context candidate with limited or single-layer support."
  )]
  out[, claim_boundary := "Do not describe as disease-validated, causal, TWAS-supported, colocalized or spatially validated. GSE73680 evidence is module-level unless a single-gene FDR endpoint is explicitly supported."]
  out[, tier_order := fcase(
    grepl("^Tier 1", final_tier_v0.4), 1L,
    grepl("^Tier 2", final_tier_v0.4), 2L,
    default = 3L
  )]
  out[, p1_order := match(gene, c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"))]
  out[is.na(p1_order), p1_order := 9999L]
  setorder(out, tier_order, p1_order, MAGMA_rank)
  out[, c("tier_order", "p1_order") := NULL]
  fwrite(out[, .(gene, MAGMA_rank, MAGMA_p, MAGMA_gene_set, scRNA_TAL_context, P1_role,
                 GSE73680_single_gene_response, GSE73680_module_context, TWAS_status,
                 coloc_status, final_tier_v0.4, interpretation, claim_boundary)],
         "results/tables/candidate_gene_tiers_v0.4.tsv", sep = "\t")
}

manuscript_v02 <- function() {
  write_lines(c(
    "# Integrative post-GWAS mapping links kidney stone risk genes to TAL-associated renal papillary cellular and disease contexts",
    "",
    "## Abstract",
    "",
    "Kidney stone disease has a substantial genetic component, but translating genome-wide association signals into renal papillary cell contexts remains challenging. We built an integrative post-GWAS framework combining locked kidney stone GWAS summary statistics, MAGMA gene-level prioritization, audited renal papillary single-nucleus RNA-seq annotations and an independent papillary plaque/stone disease-context expression dataset. MAGMA-prioritized genes converged on a Loop of Henle/thick ascending limb-associated single-nucleus context in GSE231569. A focused six-gene P1 candidate set showed an interpretable TAL, epithelial transport and calcium-handling expression spectrum rather than behaving as a uniform TAL marker panel. In GSE73680, patient-aware disease-context analyses supported MAGMA-prioritized modules in plaque/stone papilla samples, whereas the six P1 genes did not show uniform FDR-supported single-gene differential expression. These results support a MAGMA and single-nucleus-guided TAL-associated cellular context for kidney stone genetic risk and indicate that disease-context support is strongest at the MAGMA module level. The analyses do not establish genetic causality, TWAS convergence, colocalization or spatial validation.",
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
    "### GWAS QC and MAGMA gene analysis",
    "",
    "Locked kidney stone GWAS summary statistics were processed through a fixed Phase 1 workflow to retain lead-locus records and candidate gene inputs. MAGMA v1.10 was used for gene-level prioritization with a GRCh37/hg19-compatible reference. Gene sets were defined from MAGMA rankings and significance thresholds, including top 50, top 100, FDR-significant and suggestive gene sets. These gene sets were carried forward without redefining the primary GWAS or adding new loci during downstream interpretation.",
    "",
    "### GSE231569 audited annotation",
    "",
    "Renal papillary single-nucleus data from GSE231569 were audited for broad cell-type annotation and used to evaluate MAGMA-prioritized gene-set localization. The analysis focused on audited cell types, including the Loop of Henle/thick ascending limb compartment. Annotation and marker checks were performed before gene-set projection so that downstream localization used a stable cell-type framework.",
    "",
    "### scRNA gene-set projection and random benchmark",
    "",
    "MAGMA-prioritized gene sets were projected onto audited GSE231569 renal papillary cell types using module scores and random gene-set benchmarks. Locus-balanced and leave-one-locus-out analyses were used to evaluate whether the TAL-associated signal was robust to locus composition. Locus-based single-cell observations were frozen as exploratory after sensitivity analyses and were not expanded into causal claims.",
    "",
    "### P1 gene evidence scoring",
    "",
    "Six P1 candidate genes, UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2, were evaluated across audited GSE231569 cell types. For each gene, the analysis summarized TAL expression rank, average expression, detection frequency, donor-level detection, TAL specificity ratio, TAL program correlation and manuscript role assignment. The goal was to define a gene-centric evidence spectrum rather than force all P1 genes into a single mechanistic class.",
    "",
    "### GSE73680 reconstruction from Agilent feature files",
    "",
    "GSE73680 was reconstructed from supplementary Agilent feature files. The RAW tar archive was extracted into 62 TXT.gz files, all of which passed gzip validation. Agilent FEATURES blocks were parsed into a gene-level expression matrix, followed by metadata curation, gene mapping, expression QC and P1 gene availability checks. The final disease-context analysis included 55 samples from 29 patients, including 26 paired patients with control/adjacent and plaque/stone papilla samples.",
    "",
    "### Patient-aware limma duplicateCorrelation and paired sensitivity analyses",
    "",
    "Primary GSE73680 analyses used limma duplicateCorrelation to account for repeated patient-level sampling. Paired patient-level delta analyses were used as sensitivity checks. P1 single-gene responses were tested separately from module-level responses. Module scores were computed as mean z-scores across detected genes for P1 core candidates, MAGMA top 50, MAGMA top 100, MAGMA FDR, MAGMA suggestive, TAL marker and injury/remodeling marker sets.",
    "",
    "### Random gene-set and expression-matched conservative benchmarks",
    "",
    "Size-matched random gene-set benchmarks were used to evaluate whether observed module shifts exceeded random expectations for gene sets of similar size. Leave-one-gene-out analyses, MAGMA-without-P1 sensitivity and paired direction-consistency checks were used to evaluate robustness. A stricter expression-matched random benchmark was retained as a conservative boundary check; under the current implementation, it did not provide primary support and was not used as the main Figure 4 benchmark.",
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
    "### Figure 1. Integrative post-GWAS framework for mapping KSD genetic risk to renal papillary cellular and disease contexts",
    "",
    "The study links locked kidney stone disease GWAS summary statistics to MAGMA gene-level prioritization, audited GSE231569 renal papillary single-nucleus localization, a P1 candidate gene evidence spectrum and GSE73680 disease-context module analysis. TWAS, SMR/coloc and spatial transcriptomic analyses are shown as resource-limited extensions and are not treated as completed evidence layers.",
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
    "Patient-aware GSE73680 analyses showed MAGMA module-level disease-context expression association in plaque/stone papilla samples. P1 single-gene responses were heterogeneous and not FDR-supported. Size-matched random benchmarking supported MAGMA module shifts, whereas a stricter expression-matched benchmark was retained as a conservative boundary check rather than primary support. The analysis supports MAGMA-prioritized module-level disease-context association and does not establish genetic causality, TWAS convergence, colocalization or spatial validation."
  ), "manuscript/manuscript_draft_v0.2.md")
}

copy_final_draft_figures <- function() {
  candidates <- c(
    "results/figures/figure1_integrative_framework_v0.1.pdf",
    "results/figures/figure1_integrative_framework_v0.1.png",
    "results/figures/magma_scrna_benchmark.pdf",
    "results/figures/magma_leave_one_locus_out_tal.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.4.png",
    "results/figures/figure4_gse73680_disease_context_v0.7.pdf",
    "results/figures/figure4_gse73680_disease_context_v0.7.png"
  )
  file.copy(candidates[file.exists(candidates)], "results/figures/final_draft", overwrite = TRUE)
  write_lines(c(
    "# Final Draft Figure Package v0.1",
    "",
    "Included files:",
    "- figure1_integrative_framework_v0.1.pdf/png",
    "- magma_scrna_benchmark.pdf",
    "- magma_leave_one_locus_out_tal.pdf",
    "- figure3_p1_gene_evidence_v0.4.pdf/png",
    "- figure4_gse73680_disease_context_v0.7.pdf/png",
    "",
    "Note: Figure 2 is currently represented by MAGMA scRNA localization and robustness panels. A later layout pass can combine these into a single multi-panel Figure 2 if required by the target journal."
  ), "results/figures/final_draft/README.md")
}

figure4_v07()
figure1_v01()
claim_boundary_outputs()
integration_tables_v02()
candidate_tiers_v04()
manuscript_v02()
copy_final_draft_figures()

message("wrote Phase 6 hardening outputs")
