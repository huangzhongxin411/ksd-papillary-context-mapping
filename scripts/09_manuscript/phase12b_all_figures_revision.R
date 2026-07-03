suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(png)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/supervisor_review_package_v1.0", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  genetic = "#245B64",
  scrna = "#6F8F98",
  p1 = "#B59A5B",
  disease = "#9A5F52",
  muted = "#D8D8D8",
  pale = "#EEF3F4",
  text = "#303030",
  grid = "#A8A8A8"
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

make_figure1_v05 <- function() {
  layer_dt <- data.table(
    x = c(0.38, 0.54, 0.70, 0.86),
    y = 0.64,
    layer = c("Layer 1\nGenetic prioritization", "Layer 2\nsnRNA cellular context",
              "Layer 3\nCandidate gene spectrum", "Layer 4\nDisease-context association"),
    body = c("GWAS/MAGMA\n57 loci\n17,316 genes\n94 Bonferroni",
             "GSE231569 papilla\nLoop/TAL context\naudited cells",
             "P1 candidates\n6 genes\nrole spectrum",
             "GSE73680 papilla\n55 samples\n29 patients | 26 paired"),
    fill = c(pal$genetic, pal$scrna, pal$p1, pal$disease)
  )

  fig <- ggplot() +
    annotate("text", x = 0.59, y = 0.92, label = "Claim-bounded post-GWAS evidence framework",
             fontface = "bold", size = 4.9, color = pal$text) +
    annotate("text", x = 0.59, y = 0.86, label = "KSD GWAS/MAGMA -> renal papillary Loop/TAL context -> gene role spectrum -> plaque/stone disease context",
             size = 3.0, color = "#555555") +
    annotate("path", x = c(0.075, 0.10, 0.14, 0.19, 0.235, 0.245, 0.215, 0.165, 0.115, 0.075),
             y = c(0.61, 0.77, 0.86, 0.83, 0.69, 0.54, 0.43, 0.38, 0.44, 0.61),
             color = "#6D6256", linewidth = 0.6) +
    annotate("polygon", x = c(0.125, 0.158, 0.195, 0.158),
             y = c(0.47, 0.70, 0.55, 0.41), fill = "#E4D2AA", color = "#8A7A62", linewidth = 0.35) +
    annotate("rect", xmin = 0.145, xmax = 0.205, ymin = 0.41, ymax = 0.57, fill = NA, color = pal$disease, linewidth = 0.45) +
    annotate("segment", x = 0.205, xend = 0.26, y = 0.49, yend = 0.36,
             arrow = arrow(length = unit(0.10, "in")), color = "#777777") +
    annotate("rect", xmin = 0.055, xmax = 0.275, ymin = 0.08, ymax = 0.35, fill = "#FBFBFB", color = "#888888", linewidth = 0.35) +
    annotate("rect", xmin = 0.07, xmax = 0.26, ymin = 0.105, ymax = 0.315, fill = "#F1F1F1", color = NA) +
    annotate("path", x = c(0.092, 0.122, 0.152, 0.182, 0.212), y = c(0.14, 0.28, 0.14, 0.28, 0.14),
             color = pal$genetic, linewidth = 1.05) +
    annotate("rect", xmin = 0.158, xmax = 0.206, ymin = 0.17, ymax = 0.31,
             fill = "#D8C7A5", color = "#8A8A8A", alpha = 0.80) +
    annotate("point", x = 0.234, y = 0.145, size = 3.2, color = pal$disease) +
    annotate("text", x = 0.145, y = 0.102, label = "Loop/TAL", size = 2.6, color = pal$genetic, fontface = "bold") +
    annotate("text", x = 0.210, y = 0.302, label = "CD", size = 2.45, color = "#5B5147") +
    annotate("text", x = 0.237, y = 0.092, label = "plaque/\nstone", size = 2.15, color = pal$disease) +
    annotate("text", x = 0.155, y = 0.375, label = "kidney -> papilla niche", size = 2.8, color = "#555555")

  fig <- fig +
    geom_segment(data = layer_dt[1:3], aes(x = x + 0.06, xend = layer_dt$x[2:4] - 0.06, y = y, yend = y),
                 arrow = arrow(length = unit(0.10, "in")), color = "#777777", linewidth = 0.35) +
    geom_label(data = layer_dt, aes(x = x, y = y + 0.105, label = layer, fill = fill),
               color = "white", fontface = "bold", size = 2.45, lineheight = 0.92,
               label.padding = unit(0.20, "lines"), linewidth = 0) +
    geom_label(data = layer_dt, aes(x = x, y = y - 0.035, label = body),
               fill = "white", color = pal$text, size = 2.45, lineheight = 0.95,
               label.padding = unit(0.22, "lines"), label.size = 0.25) +
    annotate("rect", xmin = 0.32, xmax = 0.955, ymin = 0.17, ymax = 0.38, fill = "#F7F8F8", color = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = 0.34, y = 0.335, label = "Claim boundary box", hjust = 0, fontface = "bold", size = 3.15, color = pal$text) +
    annotate("text", x = 0.34, y = 0.285, label = "Supported inference:", hjust = 0, fontface = "bold", size = 2.85, color = pal$genetic) +
    annotate("text", x = 0.49, y = 0.285, label = "Loop/TAL-associated cellular context + MAGMA module-level disease-context association",
             hjust = 0, size = 2.65, color = pal$text) +
    annotate("text", x = 0.34, y = 0.235, label = "Not established:", hjust = 0, fontface = "bold", size = 2.85, color = pal$disease) +
    annotate("text", x = 0.475, y = 0.235, label = "causality | TWAS convergence | colocalization | spatial validation | P1 disease-gene validation",
             hjust = 0, size = 2.45, color = pal$text) +
    annotate("text", x = 0.34, y = 0.188, label = "Resource-limited extensions:", hjust = 0, fontface = "bold", size = 2.6, color = "#555555") +
    annotate("text", x = 0.515, y = 0.188, label = "TWAS / SMR-coloc / spatial prepared and audited, but not used as evidence layers",
             hjust = 0, size = 2.35, color = "#555555") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.98), ylim = c(0.06, 0.96), clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  ggsave("results/figures/figure1_integrative_framework_v0.5.pdf", fig, width = 13.2, height = 5.9, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure1_integrative_framework_v0.5.png", fig, width = 13.2, height = 5.9, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 1 Legend v0.5",
    "",
    "**Figure 1. Claim-bounded post-GWAS framework for mapping KSD genetic risk to renal papillary cellular and disease contexts.**",
    "The schematic organizes the study into four evidence layers: genetic prioritization by GWAS/MAGMA, audited GSE231569 snRNA cellular-context mapping, P1 candidate gene role-spectrum interpretation and GSE73680 plaque/stone papilla disease-context association. The left schematic illustrates a conceptual kidney-to-papilla niche containing Loop/TAL, collecting duct and plaque/stone contexts; it is not spatial validation. The supported inference is a Loop/TAL-associated cellular context plus MAGMA module-level disease-context association. Causality, TWAS convergence, colocalization, spatial validation and P1 disease-gene validation are not established."
  ), "docs/figure1_legend_v0.5.md")
}

make_figure4_v11 <- function() {
  p1 <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
  long <- fread("results/tables/gse73680_paired_module_delta_long.tsv")
  sumdt <- fread("results/tables/gse73680_paired_module_delta_summary.tsv")
  bench <- fread("results/gse73680/tables/gse73680_random_module_benchmark.tsv")

  p_a <- ggplot() +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.68, ymax = 0.86, fill = pal$genetic, color = NA) +
    annotate("text", x = 0.15, y = 0.77, label = "1\nDataset", color = "white", fontface = "bold", size = 3.0, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.68, ymax = 0.86, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.605, y = 0.77, label = "GSE73680 papillary plaque/stone context\n55 samples | 29 patients", size = 3.0, color = pal$text, lineheight = 0.9) +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.40, ymax = 0.58, fill = pal$scrna, color = NA) +
    annotate("text", x = 0.15, y = 0.49, label = "2\nPairing", color = "white", fontface = "bold", size = 3.0, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.40, ymax = 0.58, fill = "#F7F8F8", color = "#888888", linewidth = 0.3) +
    annotate("text", x = 0.605, y = 0.49, label = "26 paired patients\ncontrol/adjacent -> plaque/stone papilla", size = 3.0, color = pal$text, lineheight = 0.9) +
    annotate("rect", xmin = 0.06, xmax = 0.24, ymin = 0.13, ymax = 0.30, fill = pal$disease, color = NA) +
    annotate("text", x = 0.15, y = 0.215, label = "3\nAnalysis", color = "white", fontface = "bold", size = 3.0, lineheight = 0.9) +
    annotate("rect", xmin = 0.27, xmax = 0.94, ymin = 0.13, ymax = 0.30, fill = pal$pale, color = pal$genetic, linewidth = 0.35) +
    annotate("text", x = 0.605, y = 0.215, label = "patient-aware limma + paired delta\nmodule response + random benchmark", size = 2.9, color = pal$text, lineheight = 0.9) +
    annotate("segment", x = 0.605, xend = 0.605, y = 0.66, yend = 0.60, arrow = arrow(length = unit(0.11, "in")), color = "#777777") +
    annotate("segment", x = 0.605, xend = 0.605, y = 0.38, yend = 0.32, arrow = arrow(length = unit(0.11, "in")), color = "#777777") +
    labs(title = "A. Patient-aware GSE73680 design") +
    theme_void(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"))

  p1[, gene := factor(gene, levels = gene[order(paired_delta)])]
  p1[, signal_class := fifelse(p_value < 0.05, "Nominal only", "No FDR support")]
  p1[, p_label := sprintf("P=%.3f", p_value)]
  p_b <- ggplot(p1, aes(paired_delta, gene, fill = signal_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = p_label, hjust = ifelse(paired_delta >= 0, 1.05, -0.05)), size = 2.4, color = pal$text) +
    scale_fill_manual(values = c("Nominal only" = pal$p1, "No FDR support" = pal$muted)) +
    labs(title = "B. P1 single-gene response",
         subtitle = "No P1 gene reached FDR q <= 0.05; PKD2 nominal only",
         x = "Patient-level paired delta", y = NULL, fill = NULL) +
    theme_pub(9.0) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  keep <- c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")
  long <- long[module_label %in% keep]
  long[, module_label := factor(as.character(module_label), levels = keep)]
  lab <- sumdt[module_label %in% keep,
               .(module_label = factor(as.character(module_label), levels = keep),
                 label = sprintf("%d/%d positive", n_positive_delta, n_paired_patients))]
  p_c <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
    geom_line(aes(color = direction), alpha = 0.25, linewidth = 0.28) +
    geom_point(aes(fill = group_label), shape = 21, color = "#555555", stroke = 0.14, size = 1.15, alpha = 0.78) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.05, color = "#303030") +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2.1, color = "#303030") +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.15,
              size = 2.9, fontface = "bold", color = pal$text) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(positive = "#7F9AA3", negative = "#D6A29A", zero = pal$muted)) +
    scale_fill_manual(values = c("Control/adjacent" = "#8AA0A8", "Plaque/stone papilla" = "#B08A45")) +
    labs(title = "C. Paired patient module shifts", x = NULL, y = "Module score", color = "Delta", fill = NULL) +
    theme_pub(8.4) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

  bench <- bench[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  bench[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                                 labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  bench[, emp_label := fifelse(empirical_p == 0, "emp.P<0.001", sprintf("emp.P=%.3f", empirical_p))]
  p_d <- ggplot(bench, aes(percentile, module_label, fill = percentile >= 0.95)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "#555555", linewidth = 0.35) +
    geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
    geom_text(aes(label = emp_label), hjust = 1.04, size = 2.9, fontface = "bold", color = pal$text) +
    annotate("text", x = 0.94, y = 5.45, label = "95th percentile", hjust = 1,
             vjust = 0.5, size = 2.7, fontface = "bold", color = "#555555") +
    scale_fill_manual(values = c("TRUE" = pal$genetic, "FALSE" = pal$muted),
                      labels = c("Background-like", "Exceeds 95th percentile")) +
    coord_cartesian(xlim = c(0, 1.05)) +
    labs(title = "D. Size-matched random benchmark", x = "Random benchmark percentile", y = NULL, fill = NULL) +
    theme_pub(9.0) +
    theme(legend.position = "bottom")

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.82, 1.18))
  ggsave("results/figures/figure4_gse73680_disease_context_v1.1.pdf", fig, width = 13.2, height = 9.3, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure4_gse73680_disease_context_v1.1.png", fig, width = 13.2, height = 9.3, units = "in", dpi = 320, bg = "white")

  loo <- unique(long[, .(module_name, module_label, patient_id, paired_delta)])
  loo_summary <- rbindlist(lapply(unique(loo$module_label), function(m) {
    d <- loo[module_label == m]
    full <- d[, mean(paired_delta, na.rm = TRUE)]
    data.table(
      module_label = m,
      left_out_patient = d$patient_id,
      n_paired_patients = nrow(d) - 1,
      full_mean_delta = full,
      leave_one_mean_delta = vapply(d$patient_id, function(pid) d[patient_id != pid, mean(paired_delta, na.rm = TRUE)], numeric(1))
    )
  }), fill = TRUE)
  loo_summary[, direction_retained := sign(leave_one_mean_delta) == sign(full_mean_delta)]
  fwrite(loo_summary, "results/tables/gse73680_module_leave_one_patient_out.tsv", sep = "\t")
  loo_plot <- ggplot(loo_summary, aes(leave_one_mean_delta, module_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = pal$grid, linewidth = 0.3) +
    geom_point(aes(color = direction_retained), alpha = 0.75, size = 1.8) +
    stat_summary(fun = median, geom = "point", size = 3.0, color = "#303030") +
    scale_color_manual(values = c("TRUE" = pal$genetic, "FALSE" = pal$disease)) +
    labs(title = "GSE73680 leave-one-patient-out module delta stability",
         x = "Leave-one-patient-out mean paired delta", y = NULL, color = "Direction retained") +
    theme_pub(8.8)
  ggsave("results/figures/supp_gse73680_leave_one_patient_out.pdf", loo_plot, width = 8.4, height = 4.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/supp_gse73680_leave_one_patient_out.png", loo_plot, width = 8.4, height = 4.8, units = "in", dpi = 320, bg = "white")

  write_lines(c(
    "# Figure 4 Legend v1.1",
    "",
    "**Figure 4. GSE73680 supports MAGMA module-level plaque/stone papilla disease-context association.**",
    "(A) Patient-aware GSE73680 design showing the dataset, patient pairing and analysis layers. (B) P1 single-gene paired responses. No P1 gene reached FDR q <= 0.05; PKD2 was nominal only. (C) Paired patient module shifts between control/adjacent and plaque/stone papilla samples. Individual patient lines are shown at low opacity, and the bold black trend summarizes the median paired shift. Labels show the number of paired patients with positive shifts. (D) Size-matched random gene-set benchmark with empirical P values and the 95th percentile reference line. Percentile was computed against size-matched random gene sets. These analyses support MAGMA module-level plaque/stone disease-context association, not uniform P1 single-gene disease validation, genetic causality, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure4_legend_v1.1.md")
}

make_combined_review_file <- function() {
  figure_files <- c(
    "results/figures/figure1_integrative_framework_v0.5.png",
    "results/figures/figure2_magma_scrna_localization_v0.6.png",
    "results/figures/figure3_p1_gene_evidence_v0.8.png",
    "results/figures/figure4_gse73680_disease_context_v1.1.png",
    "results/figures/figure5_functional_context_v0.3.png"
  )
  figure_labels <- c("Figure 1 v0.5", "Figure 2 v0.6", "Figure 3 v0.8", "Figure 4 v1.1", "Figure 5 v0.3")
  pdf("results/figures/main_figures_1_to_5_review_contact_sheet_v1.0.pdf", width = 14, height = 28, onefile = TRUE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 1, heights = unit(rep(1, 5), "null"))))
  for (i in seq_along(figure_files)) {
    img <- readPNG(figure_files[i])
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(figure_labels[i], x = 0.02, y = 0.98, just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = 12, col = pal$text))
    grid.raster(img, x = 0.5, y = 0.48, width = 0.96, height = 0.88, interpolate = TRUE)
    popViewport()
  }
  popViewport()
  dev.off()

  png("results/figures/main_figures_1_to_5_review_contact_sheet_v1.0.png", width = 2800, height = 5600, res = 200)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 1, heights = unit(rep(1, 5), "null"))))
  for (i in seq_along(figure_files)) {
    img <- readPNG(figure_files[i])
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(figure_labels[i], x = 0.02, y = 0.98, just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = 12, col = pal$text))
    grid.raster(img, x = 0.5, y = 0.48, width = 0.96, height = 0.88, interpolate = TRUE)
    popViewport()
  }
  popViewport()
  dev.off()
}

update_review_package <- function() {
  file.copy("results/figures/figure1_integrative_framework_v0.5.pdf", "results/supervisor_review_package_v1.0/figure1_integrative_framework_v0.5.pdf", overwrite = TRUE)
  file.copy("results/figures/figure1_integrative_framework_v0.5.png", "results/supervisor_review_package_v1.0/figure1_integrative_framework_v0.5.png", overwrite = TRUE)
  file.copy("results/figures/figure4_gse73680_disease_context_v1.1.pdf", "results/supervisor_review_package_v1.0/figure4_gse73680_disease_context_v1.1.pdf", overwrite = TRUE)
  file.copy("results/figures/figure4_gse73680_disease_context_v1.1.png", "results/supervisor_review_package_v1.0/figure4_gse73680_disease_context_v1.1.png", overwrite = TRUE)
  file.copy("docs/figure1_legend_v0.5.md", "results/supervisor_review_package_v1.0/figure1_legend_v0.5.md", overwrite = TRUE)
  file.copy("docs/figure4_legend_v1.1.md", "results/supervisor_review_package_v1.0/figure4_legend_v1.1.md", overwrite = TRUE)
  file.copy("results/figures/main_figures_1_to_5_review_contact_sheet_v1.0.pdf", "results/supervisor_review_package_v1.0/main_figures_1_to_5_review_contact_sheet_v1.0.pdf", overwrite = TRUE)
  file.copy("results/figures/main_figures_1_to_5_review_contact_sheet_v1.0.png", "results/supervisor_review_package_v1.0/main_figures_1_to_5_review_contact_sheet_v1.0.png", overwrite = TRUE)

  fig_qc <- fread("results/tables/main_figure_qc_v1.0.tsv")
  fig_qc[figure_id == "Figure 1", `:=`(file = "results/figures/figure1_integrative_framework_v0.5.pdf",
                                      visual_status = "updated_v0.5", legend_status = "ready_v0.5")]
  fig_qc[figure_id == "Figure 4", `:=`(file = "results/figures/figure4_gse73680_disease_context_v1.1.pdf",
                                      visual_status = "updated_v1.1", legend_status = "ready_v1.1")]
  fwrite(fig_qc, "results/tables/main_figure_qc_v1.0.tsv", sep = "\t")

  final_figs <- data.table(
    figure_id = paste0("Figure ", 1:5),
    version = c("v0.5", "v0.6", "v0.8", "v1.1", "v0.3"),
    pdf = c(
      "results/figures/figure1_integrative_framework_v0.5.pdf",
      "results/figures/figure2_magma_scrna_localization_v0.6.pdf",
      "results/figures/figure3_p1_gene_evidence_v0.8.pdf",
      "results/figures/figure4_gse73680_disease_context_v1.1.pdf",
      "results/figures/figure5_functional_context_v0.3.pdf"
    ),
    png = c(
      "results/figures/figure1_integrative_framework_v0.5.png",
      "results/figures/figure2_magma_scrna_localization_v0.6.png",
      "results/figures/figure3_p1_gene_evidence_v0.8.png",
      "results/figures/figure4_gse73680_disease_context_v1.1.png",
      "results/figures/figure5_functional_context_v0.3.png"
    )
  )
  fwrite(final_figs, "results/tables/final_main_figure_set_v1.0.tsv", sep = "\t")
  package_manifest <- data.table(
    file = c(
      "figure1_integrative_framework_v0.5.pdf",
      "figure1_integrative_framework_v0.5.png",
      "figure2_magma_scrna_localization_v0.6.pdf",
      "figure2_magma_scrna_localization_v0.6.png",
      "figure3_p1_gene_evidence_v0.8.pdf",
      "figure3_p1_gene_evidence_v0.8.png",
      "figure4_gse73680_disease_context_v1.1.pdf",
      "figure4_gse73680_disease_context_v1.1.png",
      "figure5_functional_context_v0.3.pdf",
      "figure5_functional_context_v0.3.png",
      "main_figures_1_to_5_review_contact_sheet_v1.0.pdf",
      "main_figures_1_to_5_review_contact_sheet_v1.0.png"
    )
  )
  package_manifest[, path := file.path("results/supervisor_review_package_v1.0", file)]
  package_manifest[, exists := file.exists(path)]
  fwrite(package_manifest, "results/supervisor_review_package_v1.0/final_figure_manifest_v1.0.tsv", sep = "\t")
}

make_figure1_v05()
make_figure4_v11()
make_combined_review_file()
update_review_package()

message("Phase 12B all-figure revision outputs written")
