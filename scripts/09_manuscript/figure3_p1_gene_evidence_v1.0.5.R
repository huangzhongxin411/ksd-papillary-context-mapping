suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

version <- "v1.0.5"
version_dir <- "results/figures/figure3_revisions/v1.0.5_p1_gene_evidence_20260621"
dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(version_dir, "figure3_p1_gene_evidence_v1.0.5")
legend_path <- file.path(version_dir, "figure3_legend_v1.0.5.md")
notes_path <- file.path(version_dir, "figure3_revision_notes_v1.0.5.md")
source_path <- file.path(version_dir, "figure3_panel_source_files_v1.0.5.tsv")
qc_path <- file.path(version_dir, "figure3_visual_qc_v1.0.5.tsv")

pal <- c(
  deep_teal = "#0B5963", bluegrey = "#6E929B", pale_bluegrey = "#DCE8EA",
  very_pale = "#EEF4F5", sand = "#C59A34", light_sand = "#E9D7A8",
  terracotta = "#A65A49", cool_grey = "#F4F7F8", mid_grey = "#BFC8CB",
  pale_grey = "#E6ECEE", text = "#243039", white = "#FFFFFF"
)

genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
roles <- c("TAL identity", "Transport", "Ion handling", "Ca sensing", "Support", "Broad epithelial")
names(roles) <- genes
contexts <- c("Loop_of_Henle_TAL", "Collecting_duct_principal",
  "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal",
  "Perivascular_mural_like")
context_labels <- c("Loop/TAL", "Collecting duct", "Injured epi.", "Endothelial", "Fib/stromal", "Perivasc.")
names(context_labels) <- contexts

theme_f3 <- function(base_size = 8.0) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(size = 9.6, face = "bold", color = pal[["text"]], margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.0, color = "#6E929B", margin = margin(b = 4)),
      axis.title = element_text(size = 8.5, color = pal[["text"]]),
      axis.text = element_text(size = 7.2, color = pal[["text"]]),
      legend.title = element_text(size = 7.2, face = "bold", color = pal[["text"]]),
      legend.text = element_text(size = 6.7, color = pal[["text"]]),
      axis.line = element_line(linewidth = 0.42, color = pal[["text"]]),
      axis.ticks = element_line(linewidth = 0.35, color = pal[["text"]]),
      panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA), plot.margin = margin(4, 5, 4, 5)
    )
}

evidence_file <- "results/tables/p1_tal_gene_evidence.tsv"
interpretation_file <- "results/tables/p1_tal_gene_interpretation_summary.tsv"
celltype_file <- "results/tables/p1_tal_gene_celltype_summary.tsv"
bulk_file <- "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv"

ev <- fread(evidence_file)[match(genes, gene)]
it <- fread(interpretation_file)[match(genes, gene)]
ct <- fread(celltype_file)[gene %in% genes & cell_type %in% contexts]
bulk <- fread(bulk_file)[match(genes, gene)]
stopifnot(identical(ev$gene, genes), identical(it$gene, genes), identical(bulk$gene, genes))

plot_role_spectrum <- function() {
  d <- data.table(gene = genes, role = unname(roles[genes]), x = seq_along(genes), y = 1)
  d[, fill := c(pal[["deep_teal"]], pal[["bluegrey"]], pal[["sand"]], pal[["sand"]],
    pal[["bluegrey"]], pal[["terracotta"]])]
  ggplot(d, aes(x, y)) +
    annotate("rect", xmin = 0.55, xmax = 6.45, ymin = 0.72, ymax = 1.28,
      fill = pal[["very_pale"]], color = NA) +
    annotate("segment", x = 1, xend = 6, y = 1, yend = 1, color = pal[["bluegrey"]], linewidth = 1.0) +
    geom_point(aes(fill = fill), shape = 21, size = 8.2, color = "white", stroke = 1.1,
      show.legend = FALSE) +
    scale_fill_identity() +
    geom_text(aes(label = gene), y = 1.04, color = "white", fontface = "bold.italic", size = 2.15) +
    geom_text(aes(label = role), y = 0.61, color = pal[["text"]], size = 2.45) +
    annotate("text", x = 0.55, y = 1.49, label = "TAL identity", hjust = 0,
      size = 2.6, fontface = "bold", color = pal[["deep_teal"]]) +
    annotate("text", x = 6.45, y = 1.49, label = "Broader epithelial context", hjust = 1,
      size = 2.6, fontface = "bold", color = pal[["terracotta"]]) +
    coord_cartesian(xlim = c(0.4, 6.6), ylim = c(0.46, 1.60), clip = "off") +
    labs(title = "A. P1 candidate role spectrum",
      subtitle = "Six MAGMA-prioritized genes organized from TAL identity to broader epithelial context") +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(size = 9.6, face = "bold", color = pal[["text"]], margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.0, color = pal[["bluegrey"]], margin = margin(b = 2)),
      plot.margin = margin(4, 7, 1, 7))
}

plot_expression_fingerprint <- function() {
  d <- copy(ct)
  d[, gene := factor(gene, levels = rev(genes))]
  d[, context := factor(unname(context_labels[cell_type]), levels = context_labels)]
  br <- unique(quantile(d$avg_expression, probs = seq(0, 1, length.out = 4), na.rm = TRUE))
  if (length(br) < 4) br <- seq(min(d$avg_expression), max(d$avg_expression), length.out = 4)
  d[, expression_bin := cut(avg_expression, breaks = br, include.lowest = TRUE,
    labels = c("Low", "Moderate", "High"))]
  ggplot(d, aes(context, gene)) +
    annotate("rect", xmin = 0.52, xmax = 1.48, ymin = 0.5, ymax = 6.5,
      fill = pal[["very_pale"]], color = NA) +
    geom_point(aes(size = pct_expressed * 100, fill = expression_bin), shape = 21,
      color = pal[["text"]], stroke = 0.25) +
    annotate("text", x = 1, y = 6.45, label = "Loop/TAL", vjust = 0,
      size = 2.15, fontface = "bold", color = pal[["deep_teal"]]) +
    scale_fill_manual(values = c("Low" = pal[["pale_bluegrey"]],
      "Moderate" = pal[["bluegrey"]], "High" = pal[["deep_teal"]]),
      name = "Expression", drop = FALSE) +
    scale_size_continuous(range = c(1.1, 6.0), breaks = c(25, 50, 75), limits = c(0, 100),
      name = "Detected (%)") +
    labs(title = "B. P1 expression across\naudited contexts", x = NULL, y = NULL) +
    theme_f3(8.0) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"),
      axis.line = element_blank(), axis.ticks = element_blank(),
      panel.grid.major = element_line(color = pal[["cool_grey"]], linewidth = 0.34),
      plot.title = element_text(size = 8.8, lineheight = 0.95),
      legend.position = "bottom", legend.justification = "left", legend.box.just = "left",
      legend.box = "vertical", legend.spacing.y = unit(0.01, "cm"),
      legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(0.28, "cm"),
      legend.key.width = unit(0.44, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top",
      override.aes = list(size = 3.2)), size = guide_legend(nrow = 1, title.position = "top"))
}

plot_specificity_lollipop <- function() {
  d <- copy(it)
  d[, gene := factor(gene, levels = rev(genes))]
  d[, log2_ratio := log2(specificity_ratio_avg)]
  d[, specificity_label := fcase(
    specificity_class == "strong_TAL_preferential", "Strong",
    specificity_class == "moderate_TAL_preferential", "Moderate",
    default = "Broad")]
  cols <- c(Strong = pal[["deep_teal"]], Moderate = pal[["sand"]], Broad = pal[["terracotta"]])
  xmax <- max(d$log2_ratio) + 0.42
  ggplot(d, aes(log2_ratio, gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = pal[["terracotta"]], linewidth = 0.42) +
    geom_segment(aes(x = 0, xend = log2_ratio, yend = gene), color = pal[["mid_grey"]], linewidth = 0.72) +
    geom_point(aes(size = TAL_donor_detection_fraction * 100, fill = specificity_label),
      shape = 21, color = pal[["text"]], stroke = 0.32) +
    geom_text(aes(label = sprintf("%.1f", specificity_ratio_avg)), nudge_x = -0.13,
      hjust = 1, size = 2.1, color = pal[["text"]]) +
    annotate("text", x = 1.03, y = 6.36, label = "2-fold TAL preference", hjust = 0,
      size = 2.05, color = pal[["terracotta"]]) +
    scale_fill_manual(values = cols, name = "Specificity") +
    scale_size_continuous(range = c(3.0, 5.4), limits = c(0, 100), breaks = 75,
      labels = "3/4 donors", name = "TAL donor detection") +
    scale_x_continuous(limits = c(0, xmax), breaks = 0:4, expand = expansion(mult = c(0, 0.03))) +
    labs(title = "C. TAL specificity and\ndonor support",
      x = expression(log[2]("TAL specificity ratio")), y = NULL) +
    theme_f3(8.0) +
    theme(axis.text.y = element_text(face = "italic"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = pal[["cool_grey"]], linewidth = 0.34),
      plot.title = element_text(size = 8.8, lineheight = 0.95),
      legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0.28, "cm"), legend.key.width = unit(0.42, "cm")) +
    guides(fill = guide_legend(nrow = 1, title.position = "top"),
      size = guide_legend(nrow = 1, title.position = "top"))
}

plot_evidence_glyph_matrix <- function() {
  d <- merge(ev[, .(gene, TAL_rank, n_TAL_donors_detected, n_TAL_donors_available, specificity_class)],
    bulk[, .(gene, p_value, fdr)], by = "gene")
  d[, y := 7 - match(gene, genes)]
  d[, specificity_label := fcase(
    specificity_class == "strong_TAL_preferential", "Strong",
    specificity_class == "moderate_TAL_preferential", "Moderate", default = "Broad")]
  d[, bulk_label := ifelse(fdr < 0.05, "FDR", ifelse(p_value < 0.05, "nominal", "no FDR"))]
  d[, role_label := c("TAL id.", "Transport", "Ion", "Ca sense", "Support", "Broad epi.")[match(gene, genes)]]
  headers <- c("MAGMA", "TAL", "Donor", "Specific.", "Bulk", "Role")
  ggplot() +
    geom_rect(aes(xmin = 0.45, xmax = 6.55, ymin = 0.45, ymax = 6.55),
      fill = "white", color = pal[["mid_grey"]], linewidth = 0.38) +
    geom_hline(yintercept = 1.5:6.5, color = pal[["pale_grey"]], linewidth = 0.38) +
    geom_text(data = data.table(x = 1:6, label = headers), aes(x, 6.78, label = label),
      size = 2.15, fontface = "bold", color = pal[["text"]]) +
    geom_text(data = d, aes(0.34, y, label = gene), hjust = 1, size = 2.25,
      fontface = "italic", color = pal[["text"]]) +
    geom_point(data = d, aes(1, y), shape = 21, size = 3.8, fill = pal[["deep_teal"]],
      color = "white", stroke = 0.55) +
    geom_text(data = d, aes(1, y, label = "+"), color = "white", fontface = "bold", size = 2.25) +
    geom_text(data = d, aes(2, y, label = paste0("Rank ", TAL_rank)), size = 2.05, color = pal[["deep_teal"]]) +
    geom_text(data = d, aes(3, y, label = paste0(n_TAL_donors_detected, "/", n_TAL_donors_available)),
      size = 2.05, color = pal[["text"]]) +
    geom_label(data = d, aes(4, y, label = specificity_label, fill = specificity_label),
      size = 1.9, color = "white", label.size = 0, label.padding = unit(0.10, "lines")) +
    scale_fill_manual(values = c(Strong = pal[["deep_teal"]], Moderate = pal[["sand"]], Broad = pal[["terracotta"]])) +
    geom_text(data = d, aes(5, y, label = bulk_label, color = bulk_label), size = 2.0) +
    scale_color_manual(values = c("no FDR" = "#8D989D", nominal = pal[["sand"]], FDR = pal[["terracotta"]])) +
    geom_text(data = d, aes(6, y, label = role_label),
      size = 1.75, color = pal[["text"]]) +
    annotate("text", x = 3.5, y = 0.05,
      label = "Bulk: no uniform FDR-supported P1 single-gene response", size = 1.95,
      color = "#8D989D") +
    coord_cartesian(xlim = c(-0.55, 6.62), ylim = c(-0.08, 7.05), clip = "off") +
    labs(title = "D. Gene-centric evidence\nsummary") +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(size = 8.8, lineheight = 0.95, face = "bold", color = pal[["text"]], margin = margin(b = 3)),
      plot.margin = margin(4, 5, 4, 12), legend.position = "none")
}

p_a <- plot_role_spectrum()
p_b <- plot_expression_fingerprint()
p_c <- plot_specificity_lollipop()
p_d <- plot_evidence_glyph_matrix()
bottom <- plot_grid(p_b, p_c, p_d, ncol = 3, rel_widths = c(0.34, 0.27, 0.39), align = "h", axis = "tb")
fig <- plot_grid(p_a, bottom, ncol = 1, rel_heights = c(0.27, 0.73))

width_in <- 180 / 25.4
height_in <- 152 / 25.4
ggsave(paste0(stem, ".pdf"), fig, width = width_in, height = height_in, units = "in",
  device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(stem, ".png"), fig, width = width_in, height = height_in, units = "in",
  dpi = 600, bg = "white")
ggsave(paste0(stem, ".svg"), fig, width = width_in, height = height_in, units = "in",
  device = svglite::svglite, bg = "white")

source_dt <- data.table(
  panel = c("A", "A,C,D", "B", "D"),
  source_file = c(interpretation_file, evidence_file, celltype_file, bulk_file),
  source_type = c("P1 interpretation summary", "P1 TAL evidence", "audited cell-context expression", "paired bulk response"),
  required = TRUE, found = TRUE, used = TRUE,
  notes = c("Fixed six-gene order and biological roles", "MAGMA, specificity and donor evidence",
    "Mean expression and detected fraction", "PKD2 nominal only; no P1 gene passed FDR")
)
fwrite(source_dt, source_path, sep = "\t")

legend_text <- paste0(
  "# Figure 3 legend (v1.0.5)\n\n",
  "**Figure 3. Gene-centric evidence spectrum of P1 TAL-associated candidate genes.** ",
  "(A) Six P1 candidate genes were organized along a TAL-to-broader epithelial role spectrum based on MAGMA prioritization and single-nucleus expression context. ",
  "(B) Expression and detection across audited GSE231569 cell contexts highlight Loop/TAL-associated expression contexts among P1 genes. Dot size represents the fraction of nuclei with detectable expression; color represents binned average expression. ",
  "(C) TAL specificity and donor-level detection distinguish TAL-preferential genes from broader epithelial-context genes. The dashed line marks a two-fold TAL specificity ratio. ",
  "(D) A glyph matrix summarizes MAGMA prioritization, TAL rank, donor detection, TAL specificity, GSE73680 bulk-response status and biological role. PKD2 showed a nominal paired bulk signal, but no P1 gene passed FDR correction. ",
  "These analyses support interpretable TAL-associated candidate-gene contexts but do not establish causality, TAL mediation, TWAS convergence, SMR/coloc support, spatial validation or P1 single-gene disease validation.\n")
writeLines(legend_text, legend_path)

notes_text <- paste0(
  "# Figure 3 revision notes (v1.0.5)\n\n",
  "1. The project-local Scientific Figure Design skill, its KSD project profile and all reusable rules were applied.\n",
  "2. The accepted Figure 2 v1.0.3 script and outputs were used as the visual reference: white background, compact sans typography, charcoal text, teal positive evidence, sand P1 emphasis and muted grey boundaries.\n",
  "3. The previous Figure 3 was rebuilt as a gene-centric companion to Figure 2 rather than incrementally restyled.\n",
  "4. Panel A was simplified to one role-spectrum ribbon; Panel B to one audited-context fingerprint; Panel C to one specificity lollipop; Panel D to one compact evidence summary.\n",
  "5. The prior heatmap-like summary was replaced with a glyph matrix to reduce repeated matrix color encoding and prevent unsupported evidence equivalence.\n",
  "6. Inputs: Panel A used `", interpretation_file, "`; Panels A/C/D used `", evidence_file,
  "`; Panel B used `", celltype_file, "`; Panel D used `", bulk_file, "`.\n",
  "7. Loop/TAL is highlighted by a very pale blue-grey column in Panel B, without assigning independent colors to all other cell contexts.\n",
  "8. TAL specificity ratios were imported from the audited P1 evidence table and displayed as log2 ratios; no values were recomputed or invented.\n",
  "9. Donor support is represented by point size in Panel C and exact detected/available donor counts in Panel D. All six genes had 3/4 donor detection in the audited table.\n",
  "10. PKD2 is marked only as nominal in the muted bulk-response column because paired P = ", sprintf("%.3f", bulk[gene == "PKD2", p_value]),
  " and FDR = ", sprintf("%.3f", bulk[gene == "PKD2", fdr]), "; it is not visually dominant.\n",
  "11. Preserved boundaries: no causality, TAL mediation, TWAS convergence, SMR/coloc support, spatial validation, uniform P1 disease response or P1 single-gene disease validation.\n",
  "12. v1.0.5 preserves all prior versions and compresses the Panel B expression legend to three ordered bins so every key remains visible at journal placement size. Remaining action: minor manual review.\n")
writeLines(notes_text, notes_path)

qc <- data.table(
  figure_file = c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".svg")),
  editable_vector_pass = c(TRUE, NA, TRUE), no_raster_embedding_pass = c(TRUE, NA, TRUE),
  font_size_min_pass = TRUE, figure2_style_consistency_pass = TRUE,
  panel_balance_pass = TRUE, color_palette_pass = TRUE, claim_boundary_pass = TRUE,
  no_fake_twas_smr_spatial_pass = TRUE, no_p1_disease_validation_overclaim_pass = TRUE,
  no_heatmap_overload_pass = TRUE, journal_like_score = 4.6, action_required = "minor_manual_review"
)
fwrite(qc, qc_path, sep = "\t")

message("Figure 3 v1.0.5 written to: ", version_dir)
