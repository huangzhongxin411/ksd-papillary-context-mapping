suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(control = "#8AA0A8", plaque = "#B08A45", strong = "#3E6672",
            moderate = "#B59A5B", muted = "#C9C9C9", program = "#7F9AA3",
            text = "#303030")

exact <- fread("results/gse73680/tables/gse73680_exact_statistics_summary.tsv")
p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
modules <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
rand_bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

panel_a <- ggplot() +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.78, ymax = 0.94, fill = "#EEF3F4", color = pal$strong, linewidth = 0.4) +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.54, ymax = 0.70, fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.30, ymax = 0.46, fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.06, ymax = 0.22, fill = "#EEF3F4", color = pal$strong, linewidth = 0.4) +
  annotate("text", x = 0.50, y = 0.86, label = "GSE73680 disease-context dataset", fontface = "bold", size = 3.6, color = pal$text) +
  annotate("text", x = 0.50, y = 0.62, label = "55 included samples: 27 control/adjacent, 28 plaque/stone papilla", size = 3.25, color = pal$text) +
  annotate("text", x = 0.50, y = 0.38, label = "29 patients, including 26 paired patients", size = 3.25, color = pal$text) +
  annotate("text", x = 0.50, y = 0.16, label = "primary: limma duplicateCorrelation", fontface = "bold", size = 3.2, color = pal$text) +
  annotate("text", x = 0.50, y = 0.10, label = "sensitivity: paired patient-level delta", fontface = "bold", size = 3.2, color = pal$text) +
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
  scale_fill_manual(values = c("Nominal only" = pal$moderate, "No detectable" = pal$muted)) +
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
mod_plot[, q_label := fifelse(fdr <= 0.05, "q<=0.05", sprintf("q=%.3f", fdr))]
mod_plot[, label_x := fifelse(abs(paired_delta) < 0.03, 0.03, paired_delta)]
mod_plot[, label_hjust := fifelse(abs(paired_delta) < 0.03, 0, 1.05)]
panel_c <- ggplot(mod_plot, aes(paired_delta, module_label, fill = signal_class)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 0.3) +
  geom_col(width = 0.64, color = "#666666", linewidth = 0.2) +
  geom_text(aes(x = label_x, label = q_label, hjust = label_hjust), size = 2.35, color = pal$text) +
  scale_fill_manual(values = c("q <= 0.05" = pal$strong, Borderline = pal$moderate,
                               "No detectable" = pal$muted, "Disease-context program" = pal$program)) +
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
  scale_fill_manual(values = c("Exceeds random expectation" = pal$strong, Moderate = pal$moderate,
                               "Background-like" = pal$muted)) +
  labs(x = "Percentile among size-matched random gene sets", y = NULL, fill = "Benchmark",
       title = "D. Size-matched random benchmark") +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom", panel.grid.minor = element_blank())

fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, align = "hv")
ggsave("results/figures/figure4_gse73680_disease_context_v0.6.pdf", fig, width = 12.5, height = 8.8, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure4_gse73680_disease_context_v0.6.png", fig, width = 12.5, height = 8.8, units = "in", dpi = 240, bg = "white")

writeLines(c(
  "# Figure 4 Legend v0.6",
  "",
  "**Figure 4. GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression.**",
  "(A) Patient-aware GSE73680 analysis design. A total of 55 included samples, including 27 control/adjacent and 28 plaque/stone papilla samples, were analyzed from 29 patients, including 26 paired patients. Primary analyses used limma duplicateCorrelation, with paired patient-level delta analyses as sensitivity checks.",
  "(B) Paired patient-level responses of six P1 candidate genes. PKD2 showed a nominal exploratory response, whereas no P1 gene reached FDR significance.",
  "(C) Paired patient-level module responses. MAGMA-prioritized modules showed disease-context shifts, while the P1 core candidate module and TAL marker set did not show comparable support. The injury/remodeling marker set served as a disease-context reference program.",
  "(D) Size-matched random gene-set benchmark showing that MAGMA FDR, MAGMA suggestive, MAGMA top 100 and MAGMA top 50 module shifts exceeded random expectation, whereas the P1 core candidate module was background-like. The dashed line indicates the 95th percentile of size-matched random gene-set expectations.",
  "Together, these results support disease-context expression association for MAGMA-prioritized modules, rather than uniform single-gene P1 differential expression. A stricter expression-matched benchmark is treated as a boundary check rather than primary support. These findings do not establish genetic causality, TWAS convergence, colocalization or spatial validation."
), "docs/figure4_legend_v0.6.md")

message("wrote Figure 4 v0.6")
