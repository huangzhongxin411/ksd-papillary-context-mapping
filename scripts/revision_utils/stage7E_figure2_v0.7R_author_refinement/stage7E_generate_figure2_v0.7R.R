#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure2_v0.7R_author_refinement")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure2_v0.7R_author_refinement")
log_dir <- file.path(root, "logs/revision/stage7E_figure2_v0.7R_author_refinement")
for (directory in c(fig_dir, table_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

source_dir <- file.path(root, "results/tables/revision/stage4C2R_draft_figures_v0.2")
source_paths <- setNames(
  file.path(source_dir, paste0("figure2_panel", LETTERS[1:7], "_source.tsv")),
  LETTERS[1:7]
)
if (!all(file.exists(source_paths))) {
  stop("Missing frozen Figure 2 source extract(s): ", paste(source_paths[!file.exists(source_paths)], collapse = ", "))
}
src <- lapply(source_paths, fread)

pal <- c(
  deep_teal = "#245A64",
  loop_teal = "#0F6B73",
  muted_teal = "#7F9DA6",
  pale_teal = "#DCEBEC",
  amber = "#B99B5A",
  dark_amber = "#725E28",
  pale_amber = "#F3EBD8",
  terracotta = "#A76655",
  pale_coral = "#F3E3DF",
  dark = "#2D3436",
  grey = "#626B6E",
  mid_grey = "#B7C0C2",
  light_grey = "#E8ECEC",
  very_light = "#F6F8F8",
  white = "#FFFFFF"
)
font_family <- "Helvetica"

module_map <- c(
  "R1_MAGMA_Bonferroni_only" = "Bonferroni core",
  "R1_R2_R3_all_MAGMA_Bonferroni" = "All Bonferroni",
  "MAGMA_top50" = "MAGMA top 50",
  "MAGMA_top100" = "MAGMA top 100"
)
module_short_map <- c(
  "R1_MAGMA_Bonferroni_only" = "Bonf. core",
  "R1_R2_R3_all_MAGMA_Bonferroni" = "All Bonf.",
  "MAGMA_top50" = "Top 50",
  "MAGMA_top100" = "Top 100"
)
module_source_order <- c(
  "R1_MAGMA_Bonferroni_only",
  "R1_R2_R3_all_MAGMA_Bonferroni",
  "MAGMA_top50",
  "MAGMA_top100"
)
module_order <- unname(module_map[module_source_order])
module_short_order <- unname(module_short_map[module_source_order])

compartment_map <- c(
  "Collecting_duct_principal" = "Collecting duct",
  "Endothelial" = "Endothelial",
  "Fibroblast_stromal" = "Stromal",
  "Injured_undifferentiated_epithelial" = "Injured epithelial",
  "Loop_of_Henle_TAL" = "Loop/TAL",
  "Pericyte_smooth_muscle" = "Pericyte/SMC"
)
compartment_order <- c(
  "Loop/TAL", "Collecting duct", "Injured epithelial",
  "Endothelial", "Stromal", "Pericyte/SMC"
)

donor_map <- c(
  "D7290910" = "D10", "D7290912" = "D12",
  "D7290913" = "D13", "D7290914" = "D14"
)

removal_map <- c(
  "curated exemplar panel" = "Exemplar",
  "TAL marker panel" = "TAL\nmarkers",
  "calcium ion panel" = "Ca/ion",
  "top5 contributors" = "Top 5",
  "top10 contributors" = "Top 10"
)
removal_order <- unname(removal_map)

display_mapping <- rbindlist(list(
  data.table(mapping_type = "module", source_value = names(module_map), display_value = unname(module_map)),
  data.table(mapping_type = "module_compact", source_value = names(module_short_map), display_value = unname(module_short_map)),
  data.table(mapping_type = "compartment", source_value = names(compartment_map), display_value = unname(compartment_map)),
  data.table(mapping_type = "donor", source_value = names(donor_map), display_value = unname(donor_map)),
  data.table(mapping_type = "removal_panel", source_value = names(removal_map), display_value = unname(removal_map))
))
fwrite(display_mapping, file.path(table_dir, "figure2_display_label_mapping_v0.7R.tsv"), sep = "\t")

theme_pub <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = pal[["dark"]]),
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.3, color = pal[["grey"]], margin = margin(b = 3)),
      plot.caption = element_text(size = base_size - 1.0, color = pal[["grey"]], hjust = 0, margin = margin(t = 2)),
      axis.title = element_text(size = base_size - 0.2),
      axis.text = element_text(size = base_size - 0.7, color = pal[["dark"]]),
      axis.line = element_line(linewidth = 0.34, color = pal[["dark"]]),
      axis.ticks = element_line(linewidth = 0.34, color = pal[["dark"]]),
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      strip.background = element_rect(fill = pal[["very_light"]], color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.2, margin = margin(2, 2, 2, 2)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

add_panel_label <- function(plot, label) {
  plot + labs(tag = label) +
    theme(
      plot.tag = element_text(family = font_family, face = "bold", size = 10.3, color = pal[["dark"]]),
      plot.tag.position = c(0, 1)
    )
}

# Panel A: compact resource and analysis-unit strip.
nephron_base <- data.table(
  x = c(0.13, 0.09, 0.08, 0.10, 0.14, 0.19, 0.21),
  y = c(0.58, 0.50, 0.39, 0.28, 0.21, 0.32, 0.56)
)
nephron_tal <- data.table(x = c(0.14, 0.18, 0.20, 0.21), y = c(0.21, 0.31, 0.43, 0.56))
pA <- ggplot() +
  annotate("segment", x = c(0.07, 0.10, 0.13, 0.16, 0.19),
           xend = c(0.07, 0.10, 0.13, 0.16, 0.19), y = 0.84,
           yend = c(0.875, 0.91, 0.885, 0.93, 0.90),
           color = pal[["deep_teal"]], linewidth = 1.18, lineend = "round") +
  annotate("text", x = 0.23, y = 0.885, label = "MAGMA-prioritized\nmodules", hjust = 0,
           family = font_family, fontface = "bold", size = 1.72, lineheight = 0.9, color = pal[["deep_teal"]]) +
  annotate("segment", x = 0.56, xend = 0.56, y = 0.825, yend = 0.755,
           color = pal[["mid_grey"]], linewidth = 0.58,
           arrow = grid::arrow(length = grid::unit(3.4, "pt"), type = "closed")) +
  annotate("text", x = 0.56, y = 0.72, label = "GSE231569 papilla snRNA",
           family = font_family, fontface = "bold", size = 2.28, color = pal[["deep_teal"]]) +
  geom_path(data = nephron_base, aes(x, y), inherit.aes = FALSE,
            color = pal[["muted_teal"]], linewidth = 1.35, lineend = "round", linejoin = "round") +
  geom_path(data = nephron_tal, aes(x, y), inherit.aes = FALSE,
            color = pal[["loop_teal"]], linewidth = 1.9, lineend = "round") +
  annotate("point", x = 0.13, y = 0.58, shape = 21, size = 3.1,
           fill = pal[["white"]], color = pal[["muted_teal"]], stroke = 0.65) +
  annotate("text", x = 0.15, y = 0.63, label = "Loop/TAL", family = font_family,
           fontface = "bold", size = 1.45, color = pal[["loop_teal"]]) +
  annotate("segment", x = 0.15, xend = 0.205, y = 0.605, yend = 0.55,
           color = pal[["loop_teal"]], linewidth = 0.35) +
  annotate("text", x = 0.45, y = 0.58, label = "43,878", hjust = 1,
           family = font_family, fontface = "bold", size = 2.45, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.49, y = 0.58, label = "nuclei", hjust = 0,
           family = font_family, size = 1.55, color = pal[["grey"]]) +
  annotate("text", x = 0.45, y = 0.47, label = "540", hjust = 1,
           family = font_family, fontface = "bold", size = 2.45, color = pal[["loop_teal"]]) +
  annotate("text", x = 0.49, y = 0.47, label = "Loop/TAL nuclei", hjust = 0,
           family = font_family, size = 1.55, color = pal[["grey"]]) +
  annotate("text", x = 0.45, y = 0.36, label = "4", hjust = 1,
           family = font_family, fontface = "bold", size = 2.45, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.49, y = 0.36, label = "donors", hjust = 0,
           family = font_family, size = 1.55, color = pal[["grey"]]) +
  annotate("text", x = 0.45, y = 0.25, label = "6", hjust = 1,
           family = font_family, fontface = "bold", size = 2.45, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.49, y = 0.25, label = "broad compartments", hjust = 0,
           family = font_family, size = 1.45, color = pal[["grey"]]) +
  annotate("segment", x = 0.36, xend = 0.95, y = 0.525, yend = 0.525, color = pal[["light_grey"]], linewidth = 0.3) +
  annotate("segment", x = 0.36, xend = 0.95, y = 0.415, yend = 0.415, color = pal[["light_grey"]], linewidth = 0.3) +
  annotate("segment", x = 0.36, xend = 0.95, y = 0.305, yend = 0.305, color = pal[["light_grey"]], linewidth = 0.3) +
  annotate("rect", xmin = 0.34, xmax = 0.96, ymin = 0.10, ymax = 0.18,
           fill = pal[["very_light"]], color = NA) +
  annotate("text", x = 0.65, y = 0.14, label = "Donor \u00d7 compartment summaries",
           family = font_family, fontface = "bold", size = 1.4, color = pal[["grey"]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Resource and\nanalysis unit") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.1, lineheight = 0.92,
                              color = pal[["dark"]], margin = margin(b = 1, l = 13)),
    plot.margin = margin(4, 4, 4, 4)
  )

# Panel B: hero donor x compartment heatmap.
dB <- copy(src[["B"]])
dB[, module_display := factor(unname(module_map[module_name]), levels = module_order)]
dB[, compartment_display := factor(unname(compartment_map[broad_compartment]), levels = rev(compartment_order))]
dB[, donor_display := factor(unname(donor_map[donor_short]), levels = unname(donor_map))]
pB <- ggplot(dB, aes(x = donor_display, y = compartment_display, fill = mean_module_score)) +
  geom_tile(color = pal[["white"]], linewidth = 0.22) +
  geom_segment(data = data.table(module_display = factor(module_order, levels = module_order),
                                 x0 = 0.47, x1 = 0.47, y0 = 5.58, y1 = 6.42),
               aes(x = x0, xend = x1, y = y0, yend = y1), inherit.aes = FALSE,
               color = pal[["loop_teal"]], linewidth = 0.58, lineend = "round") +
  facet_wrap(~module_display, ncol = 2) +
  scale_fill_gradientn(
    colors = c(pal[["very_light"]], "#B9CFD2", pal[["deep_teal"]]),
    limits = range(dB$mean_module_score, na.rm = TRUE),
    name = "Mean score",
    guide = guide_colorbar(barheight = grid::unit(14, "mm"), barwidth = grid::unit(2.5, "mm"),
                           title.position = "top", title.hjust = 0)
  ) +
  labs(title = "Donor \u00d7 compartment module scores",
       subtitle = "Mean module score by donor \u00d7 compartment", x = "Donor", y = NULL) +
  theme_pub(7.5) +
  theme(
    axis.text.x = element_text(size = 6.8), axis.text.y = element_text(size = 7.0),
    legend.position = "right", legend.title = element_text(size = 6.4, face = "bold"),
    legend.box.spacing = grid::unit(1.5, "pt"),
    panel.spacing = grid::unit(4.5, "pt"), plot.margin = margin(4, 3, 4, 5)
  )

# Panel C: source-backed editorial summary of old Panels C, D, and F.
pC <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.67, ymax = 0.93,
           fill = "#FAFBFB", color = pal[["light_grey"]], linewidth = 0.22) +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.37, ymax = 0.63,
           fill = "#FAFBFB", color = pal[["light_grey"]], linewidth = 0.22) +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.07, ymax = 0.33,
           fill = "#FAFBFB", color = pal[["light_grey"]], linewidth = 0.22) +
  annotate("point", x = c(0.10, 0.15, 0.20, 0.25), y = 0.80, shape = 21, size = 2.6,
           fill = pal[["pale_teal"]], color = pal[["loop_teal"]], stroke = 0.38) +
  annotate("point", x = c(0.10, 0.15, 0.20, 0.25), y = 0.50, shape = 22, size = 2.15,
           fill = pal[["muted_teal"]], color = pal[["deep_teal"]], stroke = 0.36) +
  annotate("point", x = c(0.08, 0.125, 0.17, 0.215, 0.26), y = 0.20, shape = 22, size = 1.95,
           fill = pal[["muted_teal"]], color = pal[["deep_teal"]], stroke = 0.34) +
  annotate("text", x = 0.32, y = 0.845, label = "Donor consistency", hjust = 0,
           family = font_family, fontface = "bold", size = 2.05, color = pal[["dark"]]) +
  annotate("text", x = 0.32, y = 0.755, label = "Top-ranked in 4/4 donors", hjust = 0,
           family = font_family, size = 1.85, color = pal[["grey"]]) +
  annotate("text", x = 0.32, y = 0.545, label = "Single-donor exclusion", hjust = 0,
           family = font_family, fontface = "bold", size = 2.05, color = pal[["dark"]]) +
  annotate("text", x = 0.32, y = 0.455, label = "Retained after 4/4 donor exclusions", hjust = 0,
           family = font_family, size = 1.75, color = pal[["grey"]]) +
  annotate("text", x = 0.32, y = 0.245, label = "Driver-panel removal", hjust = 0,
           family = font_family, fontface = "bold", size = 2.05, color = pal[["dark"]]) +
  annotate("text", x = 0.32, y = 0.155, label = "Retained after 5/5 removal panels", hjust = 0,
           family = font_family, size = 1.75, color = pal[["grey"]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Consistency and robustness\nsummary") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.1, lineheight = 0.92,
                              color = pal[["dark"]], margin = margin(b = 2, l = 13)),
    plot.margin = margin(4, 6, 4, 4)
  )

# Panel D: source-backed percentile summary of the matched-random benchmark.
dD <- copy(src[["E"]])
dD[, module_display := factor(unname(module_map[module_name]), levels = rev(module_order))]
dD[, delta_percentile := 100 * random_percentile_for_delta]
pD <- ggplot(dD, aes(y = module_display)) +
  geom_segment(aes(x = 0, xend = 100, yend = module_display),
               color = pal[["light_grey"]], linewidth = 3.0, lineend = "round") +
  geom_point(aes(x = delta_percentile), shape = 21, size = 3.8, stroke = 0.7,
             fill = pal[["pale_amber"]], color = pal[["amber"]]) +
  annotate("label", x = 61, y = 4.96,
           label = "Partial support; rank metrics saturated\nObserved deltas above matched-random expectation\nOverall: partial support",
           hjust = 0.5, family = font_family, fontface = "bold", size = 1.62,
           lineheight = 1.04, fill = pal[["pale_amber"]], color = pal[["dark_amber"]],
           linewidth = 0.28, label.padding = grid::unit(2.0, "pt")) +
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(0, 122), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.85))) +
  labs(title = "Matched-random benchmark", subtitle = "Partial support; rank metrics saturated",
       x = "Matched-random percentile", y = NULL) +
  theme_pub(7.3) +
  theme(
    panel.grid.major.x = element_line(color = pal[["light_grey"]], linewidth = 0.25),
    panel.grid.minor = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.4), axis.text.x = element_text(size = 6.2),
    plot.margin = margin(4, 5, 4, 7)
  )

# Panel E: full-width interpretation boundary strip.
e_labels <- c(
  "Moderate\ndonor-level\nLoop/TAL context",
  "Partial\nmatched-random\nsupport",
  "Descriptive\nrobustness",
  "Not causal\ncell type",
  "Not a plaque\nnucleation site"
)
e_fill <- c(pal[["pale_teal"]], pal[["pale_amber"]], pal[["pale_teal"]], pal[["pale_coral"]], pal[["pale_coral"]])
e_border <- c(pal[["loop_teal"]], pal[["amber"]], pal[["loop_teal"]], pal[["terracotta"]], pal[["terracotta"]])
e_xmin <- c(0.02, 0.215, 0.41, 0.65, 0.82)
e_xmax <- c(0.205, 0.40, 0.595, 0.81, 0.985)
pE <- ggplot() +
  annotate("text", x = 0.305, y = 0.90, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 2.15, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.82, y = 0.90, label = "NOT CLAIMED",
           family = font_family, fontface = "bold", size = 2.15, color = pal[["terracotta"]])
for (i in seq_along(e_labels)) {
  pE <- pE +
    annotate("rect", xmin = e_xmin[i], xmax = e_xmax[i], ymin = 0.18, ymax = 0.73,
             fill = e_fill[i], color = e_border[i], linewidth = 0.45) +
    annotate("text", x = (e_xmin[i] + e_xmax[i]) / 2, y = 0.455, label = e_labels[i],
             family = font_family, fontface = "bold", size = 2.05, color = pal[["dark"]], lineheight = 0.88)
}
pE <- pE +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Interpretation boundary") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.3, color = pal[["dark"]], margin = margin(b = 1, l = 13)),
    plot.margin = margin(2, 5, 3, 6)
  )

design <- "
AAABBBBBBBBB
AAABBBBBBBBB
CCCCDDDDDDDD
EEEEEEEEEEEE
"

figure2 <- add_panel_label(pA, "A") +
  add_panel_label(pB, "B") +
  add_panel_label(pC, "C") +
  add_panel_label(pD, "D") +
  add_panel_label(pE, "E") +
  plot_layout(design = design, heights = c(1.27, 1.27, 1.05, 0.62), guides = "keep") +
  plot_annotation(
    title = "Donor-level snRNA context mapping of MAGMA-prioritized modules",
    subtitle = "Moderate Loop/TAL-associated context with partial matched-random support",
    theme = theme(
      plot.title = element_text(family = font_family, face = "bold", size = 12.2, color = pal[["dark"]], margin = margin(b = 2)),
      plot.subtitle = element_text(family = font_family, size = 8.0, color = pal[["grey"]], margin = margin(b = 6))
    )
  )

width_in <- 7.205
height_in <- 6.45
stem <- file.path(fig_dir, "figure2_snRNA_context_draft_v0.7R")

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in,
                 bg = "white", system_fonts = list(sans = font_family))
print(figure2)
dev.off()

grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
               family = font_family, bg = "white", useDingbats = FALSE)
print(figure2)
dev.off()

ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
              units = "in", res = 600, background = "white")
print(figure2)
dev.off()

if (requireNamespace("magick", quietly = TRUE)) {
  full_preview <- magick::image_read(paste0(stem, ".png"))
  half_preview <- magick::image_resize(full_preview, "50%")
  magick::image_write(half_preview, file.path(log_dir, "figure2_snRNA_context_v0.7R_qc_50pct.png"), format = "png")
} else if (requireNamespace("png", quietly = TRUE)) {
  full_preview <- png::readPNG(paste0(stem, ".png"))
  h2 <- floor(dim(full_preview)[1] / 2)
  w2 <- floor(dim(full_preview)[2] / 2)
  rows_odd <- seq.int(1, 2 * h2, by = 2)
  rows_even <- rows_odd + 1
  cols_odd <- seq.int(1, 2 * w2, by = 2)
  cols_even <- cols_odd + 1
  half_preview <- (
    full_preview[rows_odd, cols_odd, , drop = FALSE] +
      full_preview[rows_even, cols_odd, , drop = FALSE] +
      full_preview[rows_odd, cols_even, , drop = FALSE] +
      full_preview[rows_even, cols_even, , drop = FALSE]
  ) / 4
  png::writePNG(half_preview, file.path(log_dir, "figure2_snRNA_context_v0.7R_qc_50pct.png"))
}

log_lines <- c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "backend=R ggplot2/patchwork/svglite/pdf/ragg",
  "figure=Figure 2 v0.7R",
  "archetype=editorial asymmetric mixed-modality figure with hero heatmap",
  "claim=moderate donor-level Loop/TAL-associated context with partial matched-random support",
  "panel_C=old rank, donor-exclusion, and removal matrices summarized; source detail unchanged",
  "panel_D=summary-only random-percentile display; no raw random distribution available",
  "baseline=Figure 2 v0.7 production proof",
  "v0.8_layout_candidates_used=no",
  "source_data_changed_from_v0.7=no",
  "claim_changed_from_v0.7=no",
  paste0("canvas_inches=", width_in, "x", height_in),
  "font=Helvetica",
  "png_dpi=600",
  "qc_preview=50_percent_physical_scale_600dpi",
  paste0("pdf=", paste0(stem, ".pdf")),
  paste0("png=", paste0(stem, ".png")),
  paste0("svg=", paste0(stem, ".svg"))
)
writeLines(log_lines, file.path(log_dir, "stage7E_generate_figure2_v0.7R.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
