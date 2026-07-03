#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure2_refinement")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure2_refinement")
log_dir <- file.path(root, "logs/revision/stage7E_figure2_refinement")
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
module_order <- unname(module_map[c(
  "R1_MAGMA_Bonferroni_only",
  "R1_R2_R3_all_MAGMA_Bonferroni",
  "MAGMA_top50",
  "MAGMA_top100"
)])

compartment_map <- c(
  "Collecting_duct_principal" = "Collecting duct",
  "Endothelial" = "Endothelial",
  "Fibroblast_stromal" = "Stromal",
  "Injured_undifferentiated_epithelial" = "Injured epithelial",
  "Loop_of_Henle_TAL" = "Loop/TAL",
  "Pericyte_smooth_muscle" = "Pericyte/SMC"
)
compartment_order <- c(
  "Loop/TAL",
  "Collecting duct",
  "Injured epithelial",
  "Endothelial",
  "Stromal",
  "Pericyte/SMC"
)

removal_map <- c(
  "curated exemplar panel" = "Exemplar\npanel",
  "TAL marker panel" = "TAL\nmarkers",
  "calcium ion panel" = "Calcium /\nion",
  "top5 contributors" = "Top 5",
  "top10 contributors" = "Top 10"
)
removal_order <- unname(removal_map)

display_mapping <- rbindlist(list(
  data.table(mapping_type = "module", source_value = names(module_map), display_value = unname(module_map)),
  data.table(mapping_type = "compartment", source_value = names(compartment_map), display_value = unname(compartment_map)),
  data.table(mapping_type = "removal_panel", source_value = names(removal_map), display_value = unname(removal_map))
))
fwrite(display_mapping, file.path(table_dir, "figure2_display_label_mapping_v0.3.tsv"), sep = "\t")

theme_pub <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = pal[["dark"]]),
      plot.title = element_text(face = "bold", size = base_size + 1.3, hjust = 0, margin = margin(b = 2.5)),
      plot.subtitle = element_text(size = base_size - 0.2, color = pal[["grey"]], margin = margin(b = 3.5)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, color = pal[["dark"]]),
      axis.line = element_line(linewidth = 0.35, color = pal[["dark"]]),
      axis.ticks = element_line(linewidth = 0.35, color = pal[["dark"]]),
      legend.title = element_text(face = "bold", size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_rect(fill = pal[["very_light"]], color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.2, margin = margin(2, 2, 2, 2)),
      plot.margin = margin(5, 5, 5, 5)
    )
}

add_panel_label <- function(plot, label) {
  plot + labs(tag = label) +
    theme(
      plot.tag = element_text(family = font_family, face = "bold", size = 10.5, color = pal[["dark"]]),
      plot.tag.position = c(0, 1)
    )
}

# Panel A: audited resource card.
pA <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.08, ymax = 0.92,
           fill = pal[["very_light"]], color = pal[["mid_grey"]], linewidth = 0.45) +
  annotate("text", x = 0.50, y = 0.82, label = "GSE231569 papilla snRNA",
           family = font_family, fontface = "bold", size = 2.25, color = pal[["deep_teal"]]) +
  annotate("rect", xmin = 0.08, xmax = 0.47, ymin = 0.52, ymax = 0.71,
           fill = pal[["white"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("rect", xmin = 0.53, xmax = 0.92, ymin = 0.52, ymax = 0.71,
           fill = pal[["pale_teal"]], color = pal[["loop_teal"]], linewidth = 0.45) +
  annotate("rect", xmin = 0.08, xmax = 0.47, ymin = 0.29, ymax = 0.48,
           fill = pal[["white"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("rect", xmin = 0.53, xmax = 0.92, ymin = 0.29, ymax = 0.48,
           fill = pal[["white"]], color = pal[["light_grey"]], linewidth = 0.35) +
  annotate("text", x = 0.275, y = 0.635, label = "43,878\nnuclei", family = font_family,
           fontface = "bold", size = 3.2, lineheight = 0.90, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.725, y = 0.635, label = "540\nLoop/TAL", family = font_family,
           fontface = "bold", size = 3.2, lineheight = 0.90, color = pal[["loop_teal"]]) +
  annotate("text", x = 0.275, y = 0.405, label = "4\ndonors", family = font_family,
           fontface = "bold", size = 3.0, lineheight = 0.90, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.725, y = 0.405, label = "6 broad\ncompartments", family = font_family,
           fontface = "bold", size = 2.55, lineheight = 0.90, color = pal[["grey"]]) +
  annotate("segment", x = 0.08, xend = 0.92, y = 0.23, yend = 0.23,
           color = pal[["light_grey"]], linewidth = 0.45) +
  annotate("text", x = 0.50, y = 0.15, label = "Donor x compartment summaries",
           family = font_family, fontface = "bold", size = 1.95, color = pal[["grey"]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Resource overview") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.3, color = pal[["dark"]], margin = margin(b = 2)),
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel B: donor x compartment heatmap.
dB <- copy(src[["B"]])
dB[, module_display := factor(unname(module_map[module_name]), levels = module_order)]
dB[, compartment_display := factor(unname(compartment_map[broad_compartment]), levels = rev(compartment_order))]
pB <- ggplot(dB, aes(x = donor_short, y = compartment_display, fill = mean_module_score)) +
  geom_tile(color = pal[["white"]], linewidth = 0.28) +
  facet_wrap(~module_display, ncol = 2) +
  scale_fill_gradientn(
    colors = c(pal[["very_light"]], "#B9CFD2", pal[["deep_teal"]]),
    limits = range(dB$mean_module_score, na.rm = TRUE),
    name = "Mean\nmodule score"
  ) +
  labs(
    title = "Donor x compartment module scores",
    subtitle = "Mean module score by donor x compartment",
    x = "Donor",
    y = NULL
  ) +
  theme_pub(6.8) +
  theme(
    axis.text.x = element_text(size = 5.9),
    axis.text.y = element_text(size = 6.1),
    legend.position = "right",
    panel.spacing = unit(3.5, "pt")
  )

# Panel C: descriptive Loop/TAL ranks.
dC <- copy(src[["C"]])
dC[, module_display := factor(unname(module_map[module_name]), levels = rev(module_order))]
pC <- ggplot(dC, aes(x = donor_short, y = module_display)) +
  geom_point(shape = 21, size = 4.2, stroke = 0.55, fill = pal[["pale_teal"]], color = pal[["loop_teal"]]) +
  geom_text(aes(label = loop_tal_rank), family = font_family, fontface = "bold", size = 2.25, color = pal[["deep_teal"]]) +
  labs(
    title = "Loop/TAL rank by donor",
    subtitle = "Rank 1 = highest; descriptive ranks; no cell-level P values",
    x = "Donor",
    y = NULL
  ) +
  theme_pub(6.8) +
  theme(panel.grid.major = element_line(color = pal[["light_grey"]], linewidth = 0.25),
        axis.line = element_blank(), axis.ticks = element_blank())

# Panel D: single-donor exclusion status.
dD <- copy(src[["D"]])
dD[, module_display := factor(unname(module_map[module_name]), levels = rev(module_order))]
dD[, status_display := fifelse(support_retained == "yes", "retained",
                               fifelse(support_retained == "partial", "partial", "not retained"))]
pD <- ggplot(dD, aes(x = donor_short, y = module_display, fill = status_display)) +
  geom_tile(color = pal[["white"]], linewidth = 0.45) +
  geom_text(aes(label = status_display), family = font_family, size = 2.2, color = pal[["dark"]]) +
  scale_fill_manual(
    values = c("retained" = "#B8D2D4", "partial" = pal[["pale_amber"]], "not retained" = pal[["light_grey"]]),
    limits = c("retained", "partial", "not retained"),
    guide = "none"
  ) +
  labs(title = "Single-donor exclusion", subtitle = "Loop/TAL pattern under leave-one-donor-out summaries", x = "Excluded donor", y = NULL) +
  theme_pub(6.8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

# Panel E: matched-random delta benchmark, explicitly partial.
dE <- copy(src[["E"]])
dE[, module_display := factor(unname(module_map[module_name]), levels = rev(module_order))]
pE <- ggplot(dE, aes(x = observed_loop_tal_delta, y = module_display)) +
  geom_segment(aes(x = 0, xend = observed_loop_tal_delta, yend = module_display),
               color = pal[["amber"]], linewidth = 1.25, lineend = "round") +
  geom_point(shape = 21, size = 3.8, stroke = 0.65, fill = pal[["pale_amber"]], color = pal[["amber"]]) +
  geom_text(aes(label = "Partial"), hjust = -0.20, family = font_family, fontface = "bold",
            size = 2.35, color = "#725E28") +
  geom_vline(xintercept = 0, color = pal[["mid_grey"]], linewidth = 0.35) +
  coord_cartesian(xlim = c(0, max(dE$observed_loop_tal_delta) * 1.30), clip = "off") +
  labs(
    title = "Matched-random benchmark",
    subtitle = "Partial matched-random support; rank metrics saturated",
    x = "Observed Loop/TAL - other compartment delta",
    y = NULL
  ) +
  theme_pub(6.8) +
  theme(panel.grid.major.x = element_line(color = pal[["light_grey"]], linewidth = 0.25),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Panel F: panel-level known-driver removals only.
dF <- copy(src[["F"]])
dF[, module_display := factor(unname(module_map[base_module]), levels = rev(module_order))]
dF[, removal_display := factor(unname(removal_map[as.character(removal_label)]), levels = removal_order)]
dF[, status_display := fifelse(support_change == "unchanged", "retained",
                               fifelse(support_change == "weakened", "partial", "not retained"))]
pF <- ggplot(dF, aes(x = removal_display, y = module_display, fill = status_display)) +
  geom_tile(color = pal[["white"]], linewidth = 0.45) +
  geom_text(aes(label = status_display), family = font_family, size = 1.95, color = pal[["dark"]]) +
  scale_fill_manual(
    values = c("retained" = "#B8D2D4", "partial" = pal[["pale_amber"]], "not retained" = pal[["light_grey"]]),
    limits = c("retained", "partial", "not retained"),
    guide = "none"
  ) +
  labs(title = "Driver-panel removal", subtitle = "Descriptive driver-removal robustness", x = NULL, y = NULL) +
  theme_pub(6.8) +
  theme(
    axis.text.x = element_text(size = 5.7, lineheight = 0.85),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

# Panel G: compact interpretation boundary strip.
g_labels <- c(
  "Moderate\ndonor-level Loop/TAL\ncontext",
  "Partial\nmatched-random\nsupport",
  "Descriptive\ndriver-removal\nrobustness",
  "Not causal\ncell type",
  "Not plaque\nnucleation site"
)
g_fill <- c(pal[["pale_teal"]], pal[["pale_amber"]], pal[["pale_teal"]], pal[["pale_coral"]], pal[["pale_coral"]])
g_border <- c(pal[["loop_teal"]], pal[["amber"]], pal[["loop_teal"]], pal[["terracotta"]], pal[["terracotta"]])
pG <- ggplot() +
  annotate("segment", x = 0.02, xend = 0.98, y = 0.12, yend = 0.12,
           color = pal[["mid_grey"]], linewidth = 0.5) +
  annotate("text", x = 0.31, y = 0.92, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 2.3, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.82, y = 0.92, label = "NOT CLAIMED",
           family = font_family, fontface = "bold", size = 2.3, color = pal[["terracotta"]])
for (i in seq_along(g_labels)) {
  xmin <- 0.02 + (i - 1) * 0.192
  xmax <- xmin + 0.178
  pG <- pG +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = 0.24, ymax = 0.78,
             fill = g_fill[i], color = g_border[i], linewidth = 0.45) +
    annotate("text", x = (xmin + xmax) / 2, y = 0.51, label = g_labels[i],
             family = font_family, fontface = "bold", size = 2.05, color = pal[["dark"]], lineheight = 0.88)
}
pG <- pG +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Interpretation boundary") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(family = font_family, face = "bold", size = 8.3, color = pal[["dark"]], margin = margin(b = 1)),
    plot.margin = margin(3, 5, 2, 5)
  )

design <- "
AABBBB
AABBBB
CCCDDD
EEEFFF
GGGGGG
"

figure2 <- add_panel_label(pA, "A") +
  add_panel_label(pB, "B") +
  add_panel_label(pC, "C") +
  add_panel_label(pD, "D") +
  add_panel_label(pE, "E") +
  add_panel_label(pF, "F") +
  add_panel_label(pG, "G") +
  plot_layout(design = design, heights = c(1.08, 1.08, 0.95, 1.02, 0.62), guides = "keep") +
  plot_annotation(
    title = "Donor-level snRNA context mapping of MAGMA-prioritized modules",
    subtitle = "Moderate Loop/TAL-associated context with partial matched-random support",
    theme = theme(
      plot.title = element_text(family = font_family, face = "bold", size = 12.4, color = pal[["dark"]], margin = margin(b = 2)),
      plot.subtitle = element_text(family = font_family, size = 8.0, color = pal[["grey"]], margin = margin(b = 7))
    )
  )

width_in <- 7.205
height_in <- 7.45
stem <- file.path(fig_dir, "figure2_snRNA_context_draft_v0.3")

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

log_lines <- c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "backend=R ggplot2/patchwork/svglite/pdf/ragg",
  "figure=Figure 2 v0.3",
  "claim=moderate donor-level Loop/TAL-associated context with partial matched-random support",
  "source_data_changed=no",
  "claim_changed=no",
  paste0("canvas_inches=", width_in, "x", height_in),
  "font=Helvetica",
  "png_dpi=600",
  paste0("pdf=", paste0(stem, ".pdf")),
  paste0("png=", paste0(stem, ".png")),
  paste0("svg=", paste0(stem, ".svg"))
)
writeLines(log_lines, file.path(log_dir, "stage7E_generate_figure2_v0.3.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
