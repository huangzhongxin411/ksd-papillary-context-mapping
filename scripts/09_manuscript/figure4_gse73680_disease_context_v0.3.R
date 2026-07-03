suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  control = "#8AA0A8",
  plaque = "#B08A45",
  strong = "#3E6672",
  muted = "#C9C9C9",
  warn = "#B59A5B",
  text = "#303030"
)

design <- fread("results/gse73680/tables/gse73680_analysis_design.tsv")
patient_structure <- fread("results/gse73680/tables/gse73680_patient_sample_structure.tsv")
p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
modules <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

get_design <- function(item) design[design_item == item, value][1]

design_counts <- data.table(
  metric = factor(c("Included samples", "Included patients", "Paired patients"),
                  levels = c("Included samples", "Included patients", "Paired patients")),
  value = as.numeric(c(get_design("n_samples_included"), get_design("n_patients_included"),
                       get_design("n_patients_with_both_groups")))
)
group_counts <- data.table(
  group = c("Control/adjacent", "Plaque/stone papilla"),
  value = c(27, 28)
)

panel_a <- ggplot() +
  geom_col(data = design_counts, aes(metric, value), fill = pal$strong, width = 0.62) +
  geom_text(data = design_counts, aes(metric, value, label = value), vjust = -0.35, size = 3) +
  coord_cartesian(ylim = c(0, 60)) +
  labs(x = NULL, y = "Count", title = "A. GSE73680 patient-aware analysis design",
       subtitle = paste0(group_counts$group[1], "=", group_counts$value[1], "; ",
                         group_counts$group[2], "=", group_counts$value[2])) +
  theme_bw(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.3, color = "#555555"),
        panel.grid.minor = element_blank())

p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
p1[, signal_class := fifelse(fdr < 0.05, "FDR < 0.05", fifelse(p_value < 0.05, "Nominal", "No detectable"))]
p1[, p_label := sprintf("P=%.3f", p_value)]
panel_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
  geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.45, color = pal$text) +
  scale_fill_manual(values = c("FDR < 0.05" = pal$strong, Nominal = pal$warn, "No detectable" = pal$muted)) +
  labs(x = "Patient-level paired delta", y = NULL, fill = "Signal",
       title = "B. P1 gene paired response") +
  theme_bw(base_size = 9.5) +
  theme(axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

module_keep <- c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR",
                 "MAGMA_suggestive", "TAL_marker_set", "injury_remodeling_marker_set")
mod_plot <- modules[module_name %in% module_keep]
mod_plot[, module_name := factor(module_name, levels = module_name[order(paired_delta)])]
mod_plot[, signal_class := fifelse(fdr < 0.05, "FDR < 0.05", fifelse(p_value < 0.05, "Nominal", "No detectable"))]
mod_plot[, fdr_label := sprintf("FDR=%.3f", fdr)]
mod_plot[, label_x := fifelse(abs(paired_delta) < 0.03, 0.03, paired_delta)]
mod_plot[, label_hjust := fifelse(abs(paired_delta) < 0.03, 0, 1.05)]
panel_c <- ggplot(mod_plot, aes(paired_delta, module_name, fill = signal_class)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
  geom_text(aes(x = label_x, label = fdr_label, hjust = label_hjust), size = 2.35, color = pal$text) +
  scale_fill_manual(values = c("FDR < 0.05" = pal$strong, Nominal = pal$warn, "No detectable" = pal$muted)) +
  labs(x = "Patient-level paired module delta", y = NULL, fill = "Signal",
       title = "C. Module paired response") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

bench_plot <- bench[module_name %in% c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]
bench_plot[, module_name := factor(module_name, levels = module_name[order(percentile)])]
bench_plot[, benchmark_class := fifelse(benchmark_interpretation == "module response exceeds random expectation", "Exceeds random", "Background-like")]
bench_plot[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
panel_d <- ggplot(bench_plot, aes(percentile, module_name, fill = benchmark_class)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
  geom_text(aes(label = emp_label), hjust = 1.05, size = 2.35, color = pal$text) +
  coord_cartesian(xlim = c(0, 1.08)) +
  scale_fill_manual(values = c("Exceeds random" = pal$strong, "Background-like" = pal$muted)) +
  labs(x = "Observed delta percentile vs random modules", y = NULL, fill = "Benchmark",
       title = "D. MAGMA module shifts exceed random expectation") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, align = "hv")
ggsave("results/figures/figure4_gse73680_disease_context_v0.3.pdf", fig, width = 12.5, height = 8.8, units = "in", device = "pdf")
ggsave("results/figures/figure4_gse73680_disease_context_v0.3.png", fig, width = 12.5, height = 8.8, units = "in", dpi = 240)

writeLines(c(
  "# Figure 4 Legend v0.3",
  "",
  "**Figure 4. GSE73680 disease-context expression support for MAGMA-associated KSD programs.**",
  "(A) Patient-aware GSE73680 analysis design. Included samples were assigned to control/adjacent or plaque/stone papilla contexts using filename-derived curated metadata; repeated patient structure motivated limma duplicateCorrelation and paired patient-level sensitivity analyses.",
  "(B) Patient-level paired response of six P1 core candidate genes. Single-gene responses were heterogeneous and did not provide FDR-supported P1-level differential expression.",
  "(C) Patient-level paired module response for P1, MAGMA, TAL marker and injury-remodeling modules. MAGMA-prioritized modules and the injury-remodeling module showed supportive paired module-level shifts in plaque/stone papilla context.",
  "(D) Random gene-set benchmark for observed module deltas. MAGMA module responses exceeded matched random gene-set expectations, whereas the six-gene P1 module did not.",
  "Together, GSE73680 supports a disease-context expression association for MAGMA-prioritized modules rather than a uniform single-gene P1 differential-expression claim. These analyses do not establish causal mediation, spatial validation, TWAS/coloc convergence or proof of a TAL-mediated mechanism."
), "docs/figure4_legend_v0.3.md")

message("wrote Figure 4 v0.3")
