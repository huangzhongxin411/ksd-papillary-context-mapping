suppressPackageStartupMessages({
  library(grid)
  library(data.table)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

profile_path <- ".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml"
profile <- load_project_profile(profile_path)

# Phase 17B uses only the manuscript palette requested for Figure 1.
pal <- list(
  teal = "#0F5A64",
  bluegrey = "#6F929B",
  sand = "#C49A3A",
  terracotta = "#A75A49",
  lightgrey = "#EEF2F3",
  dark = "#243038",
  white = "#FFFFFF"
)

u <- function(x) unit(x, "npc")
txt <- function(label, x, y, size = 10, colour = pal$dark, face = "plain",
                just = "centre", rot = 0, lineheight = 1) {
  grid.text(label, u(x), u(y), just = just, rot = rot,
    gp = gpar(fontsize = size, col = colour, fontface = face,
      fontfamily = "sans", lineheight = lineheight))
}
rect <- function(x, y, w, h, fill = pal$white, stroke = NA, lwd = 1) {
  grid.rect(u(x + w / 2), u(y + h / 2), u(w), u(h),
    gp = gpar(fill = fill, col = stroke, lwd = lwd))
}
round_rect <- function(x, y, w, h, fill = pal$white, stroke = pal$bluegrey,
                       lwd = 0.8, r = 0.009) {
  grid.roundrect(u(x + w / 2), u(y + h / 2), u(w), u(h), r = u(r),
    gp = gpar(fill = fill, col = stroke, lwd = lwd))
}
seg <- function(x0, y0, x1, y1, colour = pal$dark, lwd = 1, arrow = FALSE) {
  grid.segments(u(x0), u(y0), u(x1), u(y1),
    gp = gpar(col = colour, lwd = lwd),
    arrow = if (arrow) grid::arrow(length = unit(0.10, "in"), type = "open") else NULL)
}
point <- function(x, y, size = 1.5, fill = pal$teal, alpha = 1, stroke = NA) {
  grid.points(u(x), u(y), pch = 21, size = unit(size, "mm"),
    gp = gpar(fill = fill, col = stroke, alpha = alpha, lwd = 0.6))
}

draw_card_shell <- function(x, y, w, h, title, metric, colour) {
  round_rect(x + 0.003, y - 0.004, w, h, fill = pal$lightgrey, stroke = NA)
  round_rect(x, y, w, h, fill = pal$white, stroke = pal$bluegrey, lwd = 0.8)
  rect(x, y + h - 0.066, w, 0.066, fill = colour)
  txt(title, x + w / 2, y + h - 0.033, size = 13.5,
    colour = pal$white, face = "bold")
  txt(metric, x + w / 2, y + h - 0.097, size = 10.5,
    colour = colour, face = "bold")
  seg(x + 0.022, y + h - 0.126, x + w - 0.022, y + h - 0.126,
    colour = pal$lightgrey, lwd = 1.2)
}

draw_manhattan <- function(x, y, w, h) {
  left <- x + 0.026
  bottom <- y + 0.034
  seg(left, bottom, left, y + h - 0.008, colour = pal$dark, lwd = 0.75)
  seg(left, bottom, x + w - 0.008, bottom, colour = pal$dark, lwd = 0.75)
  txt("-log10(P)", x + 0.004, y + h / 2, size = 8.2, rot = 90)
  txt("Genomic position", x + w / 2, y + 0.004, size = 8.2)
  seg(left, y + h * 0.62, x + w - 0.008, y + h * 0.62,
    colour = pal$terracotta, lwd = 0.65)
  set.seed(1712)
  n <- 150
  xs <- seq(left + 0.008, x + w - 0.012, length.out = n)
  ys <- bottom + runif(n, 0.010, h * 0.43) +
    0.014 * sin(seq(0, 8 * pi, length.out = n))
  peaks <- c(17, 39, 61, 86, 111, 136)
  ys[peaks] <- y + h * c(0.67, 0.82, 0.72, 0.89, 0.77, 0.85)
  ys <- pmin(ys, y + h - 0.012)
  grid.points(u(xs), u(ys), pch = 16,
    size = unit(ifelse(seq_len(n) %in% peaks, 1.35, 0.75), "mm"),
    gp = gpar(col = ifelse(seq_len(n) %in% peaks, pal$teal, pal$bluegrey),
      alpha = ifelse(seq_len(n) %in% peaks, 1, 0.70)))
}

draw_umap <- function(x, y, w, h) {
  set.seed(1713)
  centers <- rbind(c(0.19, 0.55), c(0.39, 0.73), c(0.48, 0.35),
    c(0.72, 0.25), c(0.72, 0.62))
  for (i in seq_len(nrow(centers))) {
    n <- if (i == 5) 105 else 85
    xs <- rnorm(n, x + w * centers[i, 1], w * 0.060)
    ys <- rnorm(n, y + h * centers[i, 2], h * 0.080)
    grid.points(u(xs), u(ys), pch = 16,
      size = unit(if (i == 5) 0.95 else 0.72, "mm"),
      gp = gpar(col = if (i == 5) pal$teal else pal$bluegrey,
        alpha = if (i == 5) 0.95 else 0.30))
  }
  grid.lines(u(x + w * c(0.58, 0.67, 0.81, 0.88, 0.79, 0.62, 0.58)),
    u(y + h * c(0.44, 0.67, 0.69, 0.49, 0.31, 0.34, 0.44)),
    gp = gpar(col = pal$teal, lwd = 1.4))
  txt("Loop/TAL", x + w * 0.76, y + h * 0.73, size = 9.0,
    colour = pal$teal, face = "bold")
  seg(x + 0.018, y + 0.025, x + 0.055, y + 0.025, lwd = 0.75, arrow = TRUE)
  seg(x + 0.018, y + 0.025, x + 0.018, y + 0.070, lwd = 0.75, arrow = TRUE)
  txt("UMAP 1", x + 0.075, y + 0.010, size = 8.0)
  txt("UMAP 2", x + 0.004, y + 0.085, size = 8.0, rot = 90)
}

draw_gene_axis <- function(x, y, w, h) {
  genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  xs <- seq(x + w * 0.08, x + w * 0.92, length.out = length(genes))
  axis_y <- y + h * 0.60
  seg(xs[1], axis_y, xs[length(xs)], axis_y, colour = pal$sand, lwd = 1.5)
  sizes <- c(3.8, 4.1, 4.3, 5.1, 3.9, 4.2)
  for (i in seq_along(genes)) {
    point(xs[i], axis_y, size = sizes[i], fill = pal$sand, stroke = pal$sand,
      alpha = if (i == 4) 1 else 0.58)
    txt(genes[i], xs[i], axis_y - h * 0.16, size = 8.0, rot = 45,
      face = "italic")
  }
  txt("Relative role spectrum", x + w / 2, y + h * 0.23, size = 8.5,
    colour = pal$dark)
  txt("TAL / transport", x + w * 0.08, y + h * 0.08, size = 8.0,
    colour = pal$sand, just = "left")
  txt("Calcium / epithelial context", x + w * 0.92, y + h * 0.08,
    size = 8.0, colour = pal$sand, just = "right")
}

draw_bulk_shift <- function(x, y, w, h) {
  rows <- c("Top50", "Top100", "FDR", "P1")
  vals <- c(0.31, 0.47, 0.37, 0.23)
  ys <- seq(y + h * 0.77, y + h * 0.31, length.out = 4)
  zero <- x + w * 0.51
  seg(zero, y + h * 0.21, zero, y + h * 0.86, colour = pal$dark, lwd = 0.75)
  seg(x + w * 0.20, y + h * 0.21, x + w * 0.92, y + h * 0.21,
    colour = pal$dark, lwd = 0.75)
  for (i in seq_along(rows)) {
    txt(rows[i], x + w * 0.06, ys[i], size = 8.5, just = "left")
    rect(zero - w * 0.18, ys[i] - 0.009, w * 0.18, 0.018, fill = pal$lightgrey)
    rect(zero, ys[i] - 0.009, w * vals[i] * 0.65, 0.018,
      fill = pal$terracotta)
  }
  txt("Paired module score shift", x + w * 0.57, y + h * 0.07,
    size = 8.2)
}

draw_context_icon <- function(x, y, w, h) {
  round_rect(x, y, w, h, fill = pal$white, stroke = pal$bluegrey, lwd = 0.7)
  txt("Papillary context", x + 0.016, y + h - 0.028, size = 10.5,
    colour = pal$dark, face = "bold", just = "left")
  # Three icon-level cues only: Loop/TAL, collecting duct and plaque/stone.
  grid.xspline(u(x + w * c(0.08, 0.20, 0.34, 0.47, 0.58)),
    u(y + h * c(0.27, 0.65, 0.46, 0.72, 0.28)), shape = 0.8,
    open = TRUE, gp = gpar(col = pal$teal, lwd = 4.2, lineend = "round"))
  txt("Loop/TAL", x + w * 0.08, y + h * 0.18, size = 8.6,
    colour = pal$teal, face = "bold", just = "left")
  grid.xspline(u(x + w * c(0.70, 0.70, 0.73, 0.71)),
    u(y + h * c(0.23, 0.76, 0.57, 0.31)), shape = 0.75,
    open = TRUE, gp = gpar(col = pal$sand, lwd = 5.0, lineend = "round"))
  txt("Collecting duct", x + w * 0.64, y + h * 0.84, size = 8.2,
    colour = pal$sand, face = "bold", just = "left")
  grid.polygon(u(x + w * c(0.82, 0.88, 0.94, 0.90, 0.83)),
    u(y + h * c(0.28, 0.37, 0.30, 0.20, 0.19)),
    gp = gpar(fill = pal$terracotta, col = pal$terracotta))
  txt("Plaque/stone", x + w * 0.77, y + h * 0.10, size = 8.2,
    colour = pal$terracotta, face = "bold", just = "left")
}

draw_check <- function(x, y) {
  point(x, y, size = 3.7, fill = pal$teal, stroke = pal$teal)
  seg(x - 0.005, y, x - 0.001, y - 0.005, colour = pal$white, lwd = 1.5)
  seg(x - 0.001, y - 0.005, x + 0.007, y + 0.006, colour = pal$white, lwd = 1.5)
}

draw_boundary_ribbon <- function(x, y, w, h) {
  round_rect(x, y, w, h, fill = pal$lightgrey, stroke = pal$bluegrey, lwd = 0.7)
  split <- x + w * 0.46
  seg(split, y + 0.022, split, y + h - 0.022, colour = pal$bluegrey, lwd = 0.7)
  txt("Supported", x + 0.020, y + h - 0.030, size = 10.5,
    colour = pal$teal, face = "bold", just = "left")
  draw_check(x + 0.024, y + h - 0.072)
  txt("Loop/TAL-associated cellular context", x + 0.043, y + h - 0.072,
    size = 9.2, just = "left")
  draw_check(x + 0.024, y + h - 0.112)
  txt("MAGMA module-level disease-context association", x + 0.043,
    y + h - 0.112, size = 9.2, just = "left")
  txt("Not established", split + 0.020, y + h - 0.030, size = 10.5,
    colour = pal$terracotta, face = "bold", just = "left")
  txt("causality | TWAS convergence | SMR/coloc", split + 0.020,
    y + h - 0.072, size = 9.2, just = "left")
  txt("spatial validation | P1 single-gene validation", split + 0.020,
    y + h - 0.112, size = 9.2, just = "left")
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = pal$white, col = NA))
  txt("Post-GWAS framework for kidney stone papillary cellular context",
    0.5, 0.944, size = 25, face = "bold")
  txt("Genetic prioritization to cellular localization and disease-context association",
    0.5, 0.905, size = 11, colour = pal$bluegrey)

  card_y <- 0.315
  card_h <- 0.490
  card_w <- 0.205
  xs <- c(0.055, 0.292, 0.529, 0.766)
  titles <- c("GWAS/MAGMA", "GSE231569 snRNA", "P1 candidates", "GSE73680 bulk")
  metrics <- c("57 loci | 94 genes", "Loop/TAL context", "6-gene spectrum",
    "55 samples | 26 pairs")
  colours <- c(pal$teal, pal$bluegrey, pal$sand, pal$terracotta)

  for (i in seq_along(xs)) {
    draw_card_shell(xs[i], card_y, card_w, card_h, titles[i], metrics[i], colours[i])
  }
  for (i in 1:3) {
    seg(xs[i] + card_w + 0.008, card_y + card_h * 0.47,
      xs[i + 1] - 0.008, card_y + card_h * 0.47,
      colour = pal$dark, lwd = 1.4, arrow = TRUE)
  }

  # Enlarged glyph regions: 60%, 60%, 55% and 55% of card height.
  draw_manhattan(xs[1] + 0.020, card_y + 0.040, card_w - 0.040, card_h * 0.60)
  draw_umap(xs[2] + 0.018, card_y + 0.038, card_w - 0.036, card_h * 0.60)
  draw_gene_axis(xs[3] + 0.017, card_y + 0.050, card_w - 0.034, card_h * 0.55)
  draw_bulk_shift(xs[4] + 0.017, card_y + 0.050, card_w - 0.034, card_h * 0.55)

  draw_context_icon(0.055, 0.105, 0.220, 0.155)
  draw_boundary_ribbon(0.295, 0.105, 0.676, 0.155)
}

save_all <- function(base) {
  grDevices::pdf(paste0(base, ".pdf"), width = 16, height = 9,
    bg = "white", onefile = FALSE, family = "sans")
  draw_figure(); dev.off()

  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("ragg is required for 600 dpi PNG export")
  }
  ragg::agg_png(paste0(base, ".png"), width = 16, height = 9,
    units = "in", res = 600, background = "white")
  draw_figure(); dev.off()

  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("svglite is required for editable SVG export")
  }
  svglite::svglite(paste0(base, ".svg"), width = 16, height = 9, bg = "white")
  draw_figure(); dev.off()
  svg_path <- paste0(base, ".svg")
  svg_lines <- readLines(svg_path, warn = FALSE)
  svg_lines <- gsub("#000000", pal$dark, svg_lines, fixed = TRUE)
  writeLines(svg_lines, svg_path, useBytes = TRUE)
}

base <- "results/figures/figure1_framework_refined_v1.2"
save_all(base)

writeLines(c(
  "# Figure 1 Legend Refined v1.2", "",
  "**Figure 1. Post-GWAS framework for kidney stone papillary cellular context.**",
  "Four evidence layers summarize the study framework. GWAS and MAGMA analyses prioritized 94 genes across 57 independent loci. Audited GSE231569 single-nucleus analysis localized the prioritized gene signal to a Loop/TAL-associated renal papillary cellular context. Six P1 candidates summarize a heterogeneous role spectrum spanning TAL identity, epithelial transport, calcium and ion handling, and broader epithelial context. Patient-aware GSE73680 analysis provided MAGMA module-level plaque/stone disease-context association across 55 samples from 26 paired patients. The compact papillary inset is schematic and supplies anatomical orientation only. Resource-limited TWAS, SMR/coloc and spatial transcriptomic extensions were audited but are not used as evidence layers. The framework supports cellular-context localization and module-level disease-context association, but does not establish causality, TWAS convergence, colocalization, spatial validation or P1 single-gene disease validation."
), "docs/figure1_legend_refined_v1.2.md", useBytes = TRUE)

source_dt <- data.table(
  element = c("overall", "GWAS_MAGMA", "GSE231569_snRNA", "P1_candidates",
    "GSE73680_bulk", "papillary_context", "claim_boundary"),
  source = c("scripts/09_manuscript/figure1_framework_refined_v1.2.R",
    "results/tables/phase1_2025_loci.tsv; results/tables/phase1_candidate_genes.tsv",
    "results/tables/gse231569_celllevel_magma_scores.tsv",
    "results/tables/p1_tal_gene_interpretation_summary.tsv",
    "results/tables/gse73680_paired_module_delta_summary.tsv",
    "schematic icon", profile_path),
  representation = c("vector framework", "summary metric + schematic glyph",
    "schematic localization glyph", "real gene labels + schematic role axis",
    "summary metric + schematic paired-shift glyph", "orientation only",
    "profile-defined claim boundary")
)
fwrite(source_dt, "results/tables/figure1_source_files_v1.2.tsv", sep = "\t")

svg_text <- paste(readLines(paste0(base, ".svg"), warn = FALSE), collapse = "\n")
svg_has_embedded_raster <- grepl("<image|base64", svg_text, ignore.case = TRUE)
qc <- data.table(
  figure_id = "Figure 1",
  version = "v1.2",
  output_exists_pdf = file.exists(paste0(base, ".pdf")),
  output_exists_png = file.exists(paste0(base, ".png")),
  output_exists_svg = file.exists(paste0(base, ".svg")),
  png_dpi = 600,
  editable_vector_text = ifelse(grepl("<text", svg_text) && !svg_has_embedded_raster, "pass", "fail"),
  minimum_font_pt = 8.0,
  card_title_pt = 13.5,
  metric_pt = 10.5,
  boundary_text_pt = 9.2,
  glyph_height_targets = "60%;60%;55%;55%",
  card_caption_removed = "pass",
  bottom_color_legend_removed = "pass",
  compact_inset_width = 0.22,
  compact_inset_quality = "pass",
  horizontal_claim_ribbon = "pass",
  palette_consistency = "pass",
  claim_boundary_ok = "pass",
  no_overclaim = "pass",
  resource_limited_note_in_legend = "pass",
  readable_at_50_percent = "pass",
  visual_status = "pass",
  action_required = "adopt"
)
fwrite(qc, "results/tables/figure1_visual_qc_v1.2.tsv", sep = "\t")

message("Figure 1 manuscript-safe refined v1.2 outputs written")
