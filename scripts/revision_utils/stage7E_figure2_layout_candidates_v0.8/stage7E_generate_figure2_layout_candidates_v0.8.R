#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure2_layout_candidates_v0.8")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure2_layout_candidates_v0.8")
log_dir <- file.path(root, "logs/revision/stage7E_figure2_layout_candidates_v0.8")
for (x in c(fig_dir, table_dir, log_dir)) dir.create(x, recursive = TRUE, showWarnings = FALSE)

source_dir <- file.path(root, "results/tables/revision/stage4C2R_draft_figures_v0.2")
source_paths <- setNames(file.path(source_dir, paste0("figure2_panel", LETTERS[1:7], "_source.tsv")), LETTERS[1:7])
if (!all(file.exists(source_paths))) stop("Missing frozen source tables")
src <- lapply(source_paths, fread)

font_family <- "Helvetica"
pal <- c(dark = "#2D3436", grey = "#667174", line = "#DDE3E4", light = "#F5F7F7",
         teal = "#245A64", teal2 = "#0F6B73", teal_mid = "#7F9DA6", teal_pale = "#DCEBEC",
         gold = "#B99B5A", gold_dark = "#725E28", gold_pale = "#F3EBD8",
         coral = "#A76655", coral_pale = "#F3E3DF", white = "#FFFFFF")

module_map <- c(R1_MAGMA_Bonferroni_only = "Bonferroni core",
                R1_R2_R3_all_MAGMA_Bonferroni = "All Bonferroni",
                MAGMA_top50 = "MAGMA top 50", MAGMA_top100 = "MAGMA top 100")
module_order <- unname(module_map)
compartment_map <- c(Collecting_duct_principal = "Collecting duct", Endothelial = "Endothelial",
                     Fibroblast_stromal = "Stromal", Injured_undifferentiated_epithelial = "Injured epithelial",
                     Loop_of_Henle_TAL = "Loop/TAL", Pericyte_smooth_muscle = "Pericyte/SMC")
compartment_order <- c("Loop/TAL", "Collecting duct", "Injured epithelial", "Endothelial", "Stromal", "Pericyte/SMC")
donor_map <- c(D7290910 = "D10", D7290912 = "D12", D7290913 = "D13", D7290914 = "D14")

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(text = element_text(family = font_family, colour = pal[["dark"]]),
          plot.title = element_text(face = "bold", size = base_size + 0.9, margin = margin(b = 2)),
          plot.subtitle = element_text(size = base_size - 0.2, colour = pal[["grey"]], margin = margin(b = 3)),
          axis.title = element_text(size = base_size - 0.1), axis.text = element_text(size = base_size - 0.5),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_rect(fill = pal[["light"]], colour = NA),
          strip.text = element_text(face = "bold", size = base_size - 0.1, margin = margin(2, 2, 2, 2)),
          legend.title = element_text(face = "bold", size = base_size - 0.2),
          legend.text = element_text(size = base_size - 0.6), plot.margin = margin(3, 3, 3, 3))
}

tag_panel <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag = element_text(family = font_family, face = "bold", size = 9.2), plot.tag.position = c(0, 1))

# Shared horizontal resource strip. The genetic cue is schematic and axis-free.
resource_strip <- function(title = "Resource input") {
  ggplot() +
    annotate("segment", x = c(.035, .05, .065, .08), xend = c(.035, .05, .065, .08), y = .52,
             yend = c(.66, .76, .62, .70), colour = pal[["teal"]], linewidth = 1.0, lineend = "round") +
    annotate("text", x = .10, y = .63, label = "MAGMA-prioritized modules", hjust = 0,
             family = font_family, fontface = "bold", size = 2.0, colour = pal[["teal"]]) +
    annotate("segment", x = .30, xend = .36, y = .63, yend = .63, colour = "#AAB4B6", linewidth = .55,
             arrow = grid::arrow(length = grid::unit(3, "pt"), type = "closed")) +
    annotate("text", x = .38, y = .63, label = "GSE231569 papilla snRNA", hjust = 0,
             family = font_family, fontface = "bold", size = 2.2, colour = pal[["teal2"]]) +
    annotate("text", x = .62, y = .63, label = "43,878 nuclei  |  540 Loop/TAL  |  4 donors  |  6 compartments",
             hjust = 0, family = font_family, size = 1.85, colour = pal[["dark"]]) +
    annotate("rect", xmin = .62, xmax = .97, ymin = .20, ymax = .40, fill = pal[["light"]], colour = NA) +
    annotate("text", x = .795, y = .30, label = "donor × compartment summaries",
             family = font_family, fontface = "bold", size = 1.65, colour = pal[["grey"]]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(title = title) + theme_void(base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = 8.0, margin = margin(b = 0, l = 12)),
          plot.margin = margin(3, 5, 0, 4))
}

dB <- copy(src[["B"]])
dB[, module_display := factor(unname(module_map[module_name]), levels = module_order)]
dB[, compartment_display := factor(unname(compartment_map[broad_compartment]), levels = rev(compartment_order))]
dB[, donor_display := factor(unname(donor_map[donor_short]), levels = unname(donor_map))]

heatmap_panel <- function(with_robustness_line = FALSE) {
  subtitle <- if (with_robustness_line) {
    "4/4 donors · 4/4 donor exclusions · 5/5 removal panels; details in Source Data"
  } else NULL
  ggplot(dB, aes(donor_display, compartment_display, fill = mean_module_score)) +
    geom_tile(colour = pal[["white"]], linewidth = .22) +
    geom_segment(data = data.table(module_display = factor(module_order, levels = module_order),
                                   x0 = .46, x1 = .46, y0 = 5.55, y1 = 6.45),
                 aes(x = x0, xend = x1, y = y0, yend = y1), inherit.aes = FALSE,
                 colour = pal[["teal2"]], linewidth = .75, lineend = "round") +
    facet_wrap(~module_display, ncol = 2) +
    scale_fill_gradientn(colors = c(pal[["light"]], "#B9CFD2", pal[["teal"]]),
                         limits = range(dB$mean_module_score, na.rm = TRUE), name = "Mean score",
                         guide = guide_colorbar(barheight = grid::unit(14, "mm"), barwidth = grid::unit(2.4, "mm"),
                                                title.position = "top", title.hjust = 0)) +
    labs(title = "Donor × compartment module scores", subtitle = subtitle, x = "Donor", y = NULL) +
    theme_pub(7.35) +
    theme(axis.text.x = element_text(size = 6.8), axis.text.y = element_text(size = 6.9),
          legend.position = "right", panel.spacing = grid::unit(4, "pt"), plot.margin = margin(3, 3, 3, 5))
}

robustness_column <- function(title = "Robustness summary", compact = FALSE) {
  ys <- c(.80, .50, .20)
  labels <- c("Donor consistency", "Single-donor exclusion", "Driver-panel removal")
  values <- c("Top-ranked in 4/4 donors", "Retained after 4/4 exclusions", "Retained after 5/5 panels")
  p <- ggplot()
  for (i in 1:3) {
    p <- p + annotate("rect", xmin = .04, xmax = .96, ymin = ys[i] - .11, ymax = ys[i] + .11,
                      fill = pal[["light"]], colour = pal[["line"]], linewidth = .32) +
      annotate("text", x = .10, y = ys[i] + .035, label = labels[i], hjust = 0,
               family = font_family, fontface = "bold", size = if (compact) 1.65 else 1.85) +
      annotate("text", x = .10, y = ys[i] - .045, label = values[i], hjust = 0,
               family = font_family, size = if (compact) 1.48 else 1.65, colour = pal[["grey"]])
  }
  p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + labs(title = title) +
    theme_void(base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = 7.8, margin = margin(b = 2, l = 12)),
          plot.margin = margin(3, 3, 3, 3))
}

dR <- copy(src[["E"]])
dR[, module_display := factor(unname(module_map[module_name]), levels = rev(module_order))]
dR[, pct := 100 * random_percentile_for_delta]
gauge_panel <- function(title = "Matched-random benchmark", mini = FALSE) {
  p <- ggplot(dR, aes(y = module_display)) +
    geom_segment(aes(x = 0, xend = 100, yend = module_display), colour = pal[["line"]], linewidth = 2.4, lineend = "round") +
    geom_point(aes(x = pct), shape = 21, size = if (mini) 2.7 else 3.2, stroke = .6,
               fill = pal[["gold_pale"]], colour = pal[["gold"]]) +
    scale_x_continuous(breaks = c(0, 100), limits = c(0, 103), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(.25, .85))) +
    labs(title = title, subtitle = if (mini) "Partial support\nRank metrics saturated" else NULL,
         x = "Percentile", y = NULL) + theme_pub(if (mini) 6.5 else 7.0) +
    theme(panel.grid.major.x = element_line(colour = pal[["line"]], linewidth = .25),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = if (mini) 5.5 else 6.0),
          plot.margin = margin(3, 4, 3, 4))
  if (!mini) p <- p + annotate("text", x = 0, y = 4.75, label = "Partial support  ·  Rank metrics saturated",
                               hjust = 0, family = font_family, fontface = "bold", size = 1.75,
                               colour = pal[["gold_dark"]])
  p
}

boundary_strip <- function(title = "Interpretation boundary", minimal = FALSE) {
  labels <- c("Moderate donor-level\nLoop/TAL context", "Partial matched-random\nsupport", "Descriptive\nrobustness",
              "Not causal\ncell type", "Not a plaque\nnucleation site")
  xmin <- c(.02, .215, .41, .655, .825); xmax <- c(.205, .40, .595, .815, .985)
  fills <- c(pal[["teal_pale"]], pal[["gold_pale"]], pal[["teal_pale"]], pal[["coral_pale"]], pal[["coral_pale"]])
  borders <- c(pal[["teal2"]], pal[["gold"]], pal[["teal2"]], pal[["coral"]], pal[["coral"]])
  p <- ggplot() +
    annotate("text", x = .305, y = .88, label = "SUPPORTED", family = font_family, fontface = "bold",
             size = 1.9, colour = pal[["teal"]]) +
    annotate("text", x = .82, y = .88, label = "NOT CLAIMED", family = font_family, fontface = "bold",
             size = 1.9, colour = pal[["coral"]])
  for (i in 1:5) p <- p +
    annotate("rect", xmin = xmin[i], xmax = xmax[i], ymin = .16, ymax = .69, fill = fills[i], colour = borders[i], linewidth = .4) +
    annotate("text", x = (xmin[i] + xmax[i]) / 2, y = .425, label = labels[i], family = font_family,
             fontface = "bold", size = if (minimal) 1.7 else 1.82, lineheight = .9)
  p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + labs(title = title) +
    theme_void(base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = 7.8, margin = margin(b = 0, l = 12)),
          plot.margin = margin(1, 4, 2, 4))
}

integrated_evidence <- function() {
  labels <- c("4/4 donors\ntop-ranked", "4/4 exclusions\nretained", "5/5 panels\nretained", "Matched random\npartial support")
  xs <- c(.13, .37, .61, .85)
  p <- ggplot()
  for (i in 1:4) {
    fill <- if (i == 4) pal[["gold_pale"]] else pal[["light"]]
    border <- if (i == 4) pal[["gold"]] else pal[["line"]]
    p <- p + annotate("rect", xmin = xs[i] - .105, xmax = xs[i] + .105, ymin = .25, ymax = .76,
                      fill = fill, colour = border, linewidth = .35) +
      annotate("text", x = xs[i], y = .505, label = labels[i], family = font_family,
               fontface = "bold", size = 1.8, lineheight = .9)
  }
  p + annotate("text", x = .85, y = .14, label = "Rank metrics saturated", family = font_family,
               size = 1.45, colour = pal[["gold_dark"]]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(title = "Integrated evidence summary") + theme_void(base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = 7.8, margin = margin(b = 0, l = 12)),
          plot.margin = margin(2, 4, 2, 4))
}

title_theme <- theme(plot.title = element_text(family = font_family, face = "bold", size = 11.2, colour = pal[["dark"]], margin = margin(b = 5)))

# A: heatmap-first data anchor, with a compact right inset and no robustness cards.
design_A <- "
AAAA
BBBC
BBBC
DDDD
"
fig_A <- tag_panel(resource_strip("Resource and analysis unit"), "A") +
  tag_panel(heatmap_panel(TRUE), "B") + tag_panel(gauge_panel("Matched random", mini = TRUE), "C") +
  tag_panel(boundary_strip(), "D") +
  plot_layout(design = design_A, heights = c(.55, 1.55, 1.55, .67), widths = c(1, 1, 1, .78)) +
  plot_annotation(title = "Donor-level snRNA context mapping", theme = title_theme)

# B: asymmetric evidence column. Heatmap remains larger than all supporting evidence combined.
design_B <- "
AAAA
BBBC
BBBC
BBBD
EEEE
"
fig_B <- tag_panel(resource_strip(NULL), "A") + tag_panel(heatmap_panel(FALSE), "B") +
  tag_panel(robustness_column(NULL, compact = TRUE), "C") + tag_panel(gauge_panel("Matched random", mini = TRUE), "D") +
  tag_panel(boundary_strip(NULL), "E") +
  plot_layout(design = design_B, heights = c(.50, 1.15, 1.15, .92, .64), widths = c(1, 1, 1, .78)) +
  plot_annotation(title = "Donor-level snRNA context mapping", theme = title_theme)

# C: minimum main-text version; all secondary evidence is one quiet strip.
design_C <- "
AAAA
BBBB
BBBB
CCCC
DDDD
"
fig_C <- tag_panel(resource_strip("Input and display unit"), "A") + tag_panel(heatmap_panel(FALSE), "B") +
  tag_panel(integrated_evidence(), "C") + tag_panel(boundary_strip(minimal = TRUE), "D") +
  plot_layout(design = design_C, heights = c(.50, 1.45, 1.45, .58, .62)) +
  plot_annotation(title = "Donor-level snRNA context mapping", theme = title_theme)

save_candidate <- function(plot, stem, width = 7.205, height = 5.45) {
  path <- file.path(fig_dir, stem)
  svglite::svglite(paste0(path, ".svg"), width = width, height = height, bg = "white",
                   system_fonts = list(sans = font_family)); print(plot); dev.off()
  grDevices::pdf(paste0(path, ".pdf"), width = width, height = height, family = font_family,
                 bg = "white", useDingbats = FALSE); print(plot); dev.off()
  ragg::agg_png(paste0(path, ".png"), width = width, height = height, units = "in", res = 600,
                background = "white"); print(plot); dev.off()
  if (requireNamespace("magick", quietly = TRUE)) {
    im <- magick::image_read(paste0(path, ".png")); im <- magick::image_resize(im, "50%")
    magick::image_write(im, file.path(log_dir, paste0(stem, "_qc_50pct.png")), format = "png")
  }
}

save_candidate(fig_A, "figure2_candidate_A_data_anchor")
save_candidate(fig_B, "figure2_candidate_B_evidence_column")
save_candidate(fig_C, "figure2_candidate_C_minimal_maintext")

log_lines <- c(paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
               "backend=R ggplot2/patchwork/svglite/pdf/ragg", "source_data_changed=no", "claim_changed=no",
               "candidate_A=data_anchor", "candidate_B=evidence_column", "candidate_C=minimal_maintext",
               "canvas_inches=7.205x5.45", "png_dpi=600", "figure3_started=no")
writeLines(log_lines, file.path(log_dir, "stage7E_generate_figure2_layout_candidates_v0.8.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
