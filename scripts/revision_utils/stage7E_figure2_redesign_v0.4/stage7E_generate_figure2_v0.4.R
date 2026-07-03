#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure2_redesign_v0.4")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure2_redesign_v0.4")
log_dir <- file.path(root, "logs/revision/stage7E_figure2_redesign_v0.4")
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
fwrite(display_mapping, file.path(table_dir, "figure2_display_label_mapping_v0.4.tsv"), sep = "\t")

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

# Panel A: compact resource overview with an abstract papillary nephron schematic.
nephron_base <- data.table(
  x = c(0.15, 0.11, 0.10, 0.12, 0.18, 0.24, 0.27),
  y = c(0.72, 0.62, 0.42, 0.27, 0.20, 0.34, 0.67)
)
nephron_tal <- data.table(
  x = c(0.18, 0.22, 0.25, 0.27),
  y = c(0.20, 0.30, 0.45, 0.67)
)
pA <- ggplot() +
  annotate("text", x = 0.50, y = 0.91, label = "GSE231569 papilla snRNA",
           family = font_family, fontface = "bold", size = 2.55, color = pal[["deep_teal"]]) +
  geom_path(data = nephron_base, aes(x, y), inherit.aes = FALSE,
            color = pal[["muted_teal"]], linewidth = 1.55, lineend = "round", linejoin = "round") +
  geom_path(data = nephron_tal, aes(x, y), inherit.aes = FALSE,
            color = pal[["loop_teal"]], linewidth = 2.2, lineend = "round") +
  annotate("point", x = 0.15, y = 0.72, shape = 21, size = 4.0,
           fill = pal[["white"]], color = pal[["muted_teal"]], stroke = 0.7) +
  annotate("text", x = 0.20, y = 0.80, label = "papillary nephron",
           family = font_family, size = 1.55, color = pal[["grey"]]) +
  annotate("text", x = 0.285, y = 0.48, label = "Loop/TAL",
           family = font_family, fontface = "bold", size = 1.9, color = pal[["loop_teal"]], angle = 90) +
  annotate("rect", xmin = 0.40, xmax = 0.67, ymin = 0.52, ymax = 0.73,
           fill = pal[["very_light"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("rect", xmin = 0.70, xmax = 0.97, ymin = 0.52, ymax = 0.73,
           fill = pal[["pale_teal"]], color = pal[["loop_teal"]], linewidth = 0.42) +
  annotate("rect", xmin = 0.40, xmax = 0.67, ymin = 0.27, ymax = 0.48,
           fill = pal[["very_light"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("rect", xmin = 0.70, xmax = 0.97, ymin = 0.27, ymax = 0.48,
           fill = pal[["very_light"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("text", x = 0.535, y = 0.64, label = "43,878", family = font_family,
           fontface = "bold", size = 3.25, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.535, y = 0.565, label = "nuclei", family = font_family,
           size = 1.75, color = pal[["grey"]]) +
  annotate("text", x = 0.835, y = 0.64, label = "540", family = font_family,
           fontface = "bold", size = 3.25, color = pal[["loop_teal"]]) +
  annotate("text", x = 0.835, y = 0.565, label = "Loop/TAL", family = font_family,
           size = 1.75, color = pal[["grey"]]) +
  annotate("text", x = 0.535, y = 0.39, label = "4", family = font_family,
           fontface = "bold", size = 3.0, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.535, y = 0.315, label = "donors", family = font_family,
           size = 1.75, color = pal[["grey"]]) +
  annotate("text", x = 0.835, y = 0.39, label = "6", family = font_family,
           fontface = "bold", size = 3.0, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.835, y = 0.335, label = "broad\ncompartments", family = font_family,
           size = 1.45, color = pal[["grey"]], lineheight = 0.88) +
  annotate("segment", x = 0.40, xend = 0.97, y = 0.19, yend = 0.19,
           color = pal[["light_grey"]], linewidth = 0.4) +
  annotate("text", x = 0.685, y = 0.12, label = "Donor x compartment summaries",
           family = font_family, fontface = "bold", size = 1.55, color = pal[["grey"]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Resource overview") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.5, color = pal[["dark"]], margin = margin(b = 1)),
    plot.margin = margin(4, 5, 4, 4)
  )

# Panel B: hero donor x compartment heatmap.
dB <- copy(src[["B"]])
dB[, module_display := factor(unname(module_map[module_name]), levels = module_order)]
dB[, compartment_display := factor(unname(compartment_map[broad_compartment]), levels = rev(compartment_order))]
dB[, donor_display := factor(unname(donor_map[donor_short]), levels = unname(donor_map))]
pB <- ggplot(dB, aes(x = donor_display, y = compartment_display, fill = mean_module_score)) +
  geom_tile(color = pal[["white"]], linewidth = 0.34) +
  facet_wrap(~module_display, ncol = 2) +
  scale_fill_gradientn(
    colors = c(pal[["very_light"]], "#B9CFD2", pal[["deep_teal"]]),
    limits = range(dB$mean_module_score, na.rm = TRUE),
    name = "Mean score",
    guide = guide_colorbar(barheight = grid::unit(22, "mm"), barwidth = grid::unit(3.2, "mm"),
                           title.position = "top", title.hjust = 0)
  ) +
  labs(
    title = "Donor x compartment module scores",
    subtitle = "Mean module score by donor x compartment",
    x = "Donor",
    y = NULL
  ) +
  theme_pub(7.4) +
  theme(
    axis.text.x = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.6),
    legend.position = "right",
    legend.title = element_text(size = 6.5, face = "bold"),
    panel.spacing = grid::unit(5, "pt"),
    plot.margin = margin(4, 3, 4, 5)
  )

# Panel C: compact rank tile matrix.
dC <- copy(src[["C"]])
dC[, module_display := factor(unname(module_short_map[module_name]), levels = rev(module_short_order))]
dC[, donor_display := factor(unname(donor_map[donor_short]), levels = unname(donor_map))]
pC <- ggplot(dC, aes(x = donor_display, y = module_display)) +
  geom_tile(fill = pal[["pale_teal"]], color = pal[["white"]], linewidth = 0.55) +
  geom_text(aes(label = loop_tal_rank), family = font_family, fontface = "bold",
            size = 2.45, color = pal[["deep_teal"]]) +
  labs(title = "Loop/TAL rank", subtitle = "Rank 1 = highest; descriptive only", x = "Donor", y = NULL) +
  theme_pub(7.0) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 6.0), axis.text.y = element_text(size = 6.1))

# Panel D: compact binary single-donor exclusion matrix.
dD <- copy(src[["D"]])
dD[, module_display := factor(unname(module_short_map[module_name]), levels = rev(module_short_order))]
dD[, donor_display := factor(unname(donor_map[donor_short]), levels = unname(donor_map))]
dD[, retained := support_retained == "yes"]
pD <- ggplot(dD, aes(x = donor_display, y = module_display)) +
  geom_tile(fill = pal[["very_light"]], color = pal[["white"]], linewidth = 0.55) +
  geom_point(data = dD[retained == TRUE], shape = 22, size = 4.0, stroke = 0.45,
             fill = pal[["muted_teal"]], color = pal[["deep_teal"]]) +
  labs(title = "Single-donor exclusion",
       subtitle = "Loop/TAL after leave-one-donor-out\nFilled = retained",
       x = "Excluded donor", y = NULL) +
  theme_pub(7.0) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 6.0), axis.text.y = element_text(size = 6.1))

# Panel E: summary-only matched-random benchmark display.
# Raw random-set distributions were not retained; the frozen percentile summary is shown directly.
dE <- copy(src[["E"]])
dE[, module_display := factor(unname(module_short_map[module_name]), levels = rev(module_short_order))]
dE[, delta_percentile := 100 * random_percentile_for_delta]
pE <- ggplot(dE, aes(y = module_display)) +
  geom_segment(aes(x = 0, xend = 100, yend = module_display),
               color = pal[["light_grey"]], linewidth = 2.8, lineend = "round") +
  geom_point(aes(x = delta_percentile), shape = 21, size = 3.5, stroke = 0.65,
             fill = pal[["pale_amber"]], color = pal[["amber"]]) +
  geom_text(aes(x = 103, label = "partial"), hjust = 0, family = font_family,
            fontface = "bold", size = 1.95, color = pal[["dark_amber"]]) +
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(0, 126), expand = c(0, 0)) +
  labs(title = "Matched-random benchmark", subtitle = "Partial matched-random support",
       x = "Observed delta percentile", y = NULL) +
  theme_pub(7.0) +
  theme(panel.grid.major.x = element_line(color = pal[["light_grey"]], linewidth = 0.25),
        panel.grid.minor = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6.1), axis.text.x = element_text(size = 6.0))

# Panel F: compact binary driver-panel removal matrix.
dF <- copy(src[["F"]])
dF[, module_display := factor(unname(module_short_map[base_module]), levels = rev(module_short_order))]
dF[, removal_display := factor(unname(removal_map[as.character(removal_label)]), levels = removal_order)]
dF[, retained := support_change == "unchanged"]
pF <- ggplot(dF, aes(x = removal_display, y = module_display)) +
  geom_tile(fill = pal[["very_light"]], color = pal[["white"]], linewidth = 0.55) +
  geom_point(data = dF[retained == TRUE], shape = 22, size = 4.0, stroke = 0.45,
             fill = pal[["muted_teal"]], color = pal[["deep_teal"]]) +
  labs(title = "Driver-panel removal", subtitle = "Descriptive robustness; filled = retained", x = NULL, y = NULL) +
  theme_pub(7.1) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 5.9, lineheight = 0.88), axis.text.y = element_text(size = 6.2))

# Panel G: concise grouped claim boundary strip.
g_labels <- c(
  "Moderate\ndonor-level\nLoop/TAL context",
  "Partial\nmatched-\nrandom support",
  "Descriptive\ndriver-\nremoval\nrobustness",
  "Not causal\ncell type",
  "Not plaque\nnucleation site"
)
g_fill <- c(pal[["pale_teal"]], pal[["pale_amber"]], pal[["pale_teal"]], pal[["pale_coral"]], pal[["pale_coral"]])
g_border <- c(pal[["loop_teal"]], pal[["amber"]], pal[["loop_teal"]], pal[["terracotta"]], pal[["terracotta"]])
g_xmin <- c(0.02, 0.215, 0.41, 0.64, 0.815)
g_xmax <- c(0.205, 0.40, 0.595, 0.805, 0.98)
pG <- ggplot() +
  annotate("text", x = 0.305, y = 0.90, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 2.0, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.81, y = 0.90, label = "NOT CLAIMED",
           family = font_family, fontface = "bold", size = 2.0, color = pal[["terracotta"]])
for (i in seq_along(g_labels)) {
  pG <- pG +
    annotate("rect", xmin = g_xmin[i], xmax = g_xmax[i], ymin = 0.20, ymax = 0.73,
             fill = g_fill[i], color = g_border[i], linewidth = 0.42) +
    annotate("text", x = (g_xmin[i] + g_xmax[i]) / 2, y = 0.465, label = g_labels[i],
             family = font_family, fontface = "bold", size = 1.62, color = pal[["dark"]], lineheight = 0.86)
}
pG <- pG +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Interpretation boundary") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.5, color = pal[["dark"]], margin = margin(b = 1, l = 14)),
    plot.margin = margin(3, 4, 3, 7)
  )

design <- "
AAAABBBBBBBB
AAAABBBBBBBB
CCCCDDDDEEEE
FFFFGGGGGGGG
"

figure2 <- add_panel_label(pA, "A") +
  add_panel_label(pB, "B") +
  add_panel_label(pC, "C") +
  add_panel_label(pD, "D") +
  add_panel_label(pE, "E") +
  add_panel_label(pF, "F") +
  add_panel_label(pG, "G") +
  plot_layout(design = design, heights = c(1.18, 1.18, 1.02, 0.92), guides = "keep") +
  plot_annotation(
    title = "Donor-level snRNA context mapping of MAGMA-prioritized modules",
    subtitle = "Moderate Loop/TAL-associated context with partial matched-random support",
    theme = theme(
      plot.title = element_text(family = font_family, face = "bold", size = 12.2, color = pal[["dark"]], margin = margin(b = 2)),
      plot.subtitle = element_text(family = font_family, size = 8.0, color = pal[["grey"]], margin = margin(b = 6))
    )
  )

width_in <- 7.205
height_in <- 7.05
stem <- file.path(fig_dir, "figure2_snRNA_context_draft_v0.4")

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
  magick::image_write(half_preview, file.path(log_dir, "figure2_snRNA_context_v0.4_qc_50pct.png"), format = "png")
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
  png::writePNG(half_preview, file.path(log_dir, "figure2_snRNA_context_v0.4_qc_50pct.png"))
}

log_lines <- c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "backend=R ggplot2/patchwork/svglite/pdf/ragg",
  "figure=Figure 2 v0.4",
  "archetype=asymmetric quantitative grid with hero heatmap",
  "claim=moderate donor-level Loop/TAL-associated context with partial matched-random support",
  "panel_E=summary-only random-percentile display; no raw random distribution available",
  "source_data_changed=no",
  "claim_changed=no",
  paste0("canvas_inches=", width_in, "x", height_in),
  "font=Helvetica",
  "png_dpi=600",
  "qc_preview=50_percent_physical_scale_600dpi",
  paste0("pdf=", paste0(stem, ".pdf")),
  paste0("png=", paste0(stem, ".png")),
  paste0("svg=", paste0(stem, ".svg"))
)
writeLines(log_lines, file.path(log_dir, "stage7E_generate_figure2_v0.4.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
