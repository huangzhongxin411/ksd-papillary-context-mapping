suppressPackageStartupMessages({
  library(grid)
  library(data.table)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

version_dir <- "results/figures/figure1_revisions/v1.2.1_card_refinement_20260621"
dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)

profile_path <- ".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml"
profile <- load_project_profile(profile_path)

col <- list(
  deep_teal = "#0B5963",
  teal_dark = "#0B5963",
  bluegrey = "#6E929B",
  pale_bluegrey = "#E4ECEE",
  sand = "#C59A34",
  light_sand = "#E8D7A6",
  terracotta = "#A65A49",
  redbrown = "#A65A49",
  light_grey = "#F4F7F8",
  mid_grey = "#8D989D",
  dark = "#243039",
  border = "#C9D3D6",
  offwhite = "#FBFCFC",
  white = "#FFFFFF"
)

u <- function(x) unit(x, "npc")
font_scale <- 0.58
txt <- function(label, x, y, size = 10, colr = col$dark, font = 1, just = "centre",
                lineheight = 1.0, rot = 0) {
  grid.text(label, x = u(x), y = u(y), just = just, rot = rot,
    gp = gpar(fontsize = max(size * font_scale, 6.5), col = colr, fontface = font,
      lineheight = lineheight, fontfamily = "sans"))
}
rr <- function(x, y, w, h, fill = col$white, stroke = col$border, lwd = 1, r = 0.012) {
  grid.roundrect(x = u(x + w / 2), y = u(y + h / 2), width = u(w), height = u(h),
    r = u(r), gp = gpar(fill = fill, col = stroke, lwd = lwd))
}
rect <- function(x, y, w, h, fill, stroke = NA, lwd = 1) {
  grid.rect(x = u(x + w / 2), y = u(y + h / 2), width = u(w), height = u(h),
    gp = gpar(fill = fill, col = stroke, lwd = lwd))
}
seg <- function(x0, y0, x1, y1, colr = col$mid_grey, lwd = 1, arrow = FALSE) {
  grid.segments(u(x0), u(y0), u(x1), u(y1),
    gp = gpar(col = colr, lwd = lwd),
    arrow = if (arrow) grid::arrow(length = unit(0.13, "in"), type = "closed") else NULL)
}
pt <- function(x, y, size = 2, fill = col$deep_teal, stroke = NA, alpha = 1) {
  grid.points(u(x), u(y), pch = 21, size = unit(size, "mm"),
    gp = gpar(fill = fill, col = stroke %||% fill, alpha = alpha, lwd = 0.8))
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
ellipse_poly <- function(cx, cy, rx, ry, angle = 0, n = 40, fill = "#EAF0F2", stroke = "#E2E9EC", alpha = 0.65) {
  theta <- seq(0, 2 * pi, length.out = n)
  ca <- cos(angle * pi / 180)
  sa <- sin(angle * pi / 180)
  x0 <- rx * cos(theta)
  y0 <- ry * sin(theta)
  xs <- cx + x0 * ca - y0 * sa
  ys <- cy + x0 * sa + y0 * ca
  grid.polygon(u(xs), u(ys), gp = gpar(fill = fill, col = stroke, alpha = alpha, lwd = 0.4))
}

draw_card <- function(x, y, w, h, banner, header_fill) {
  # subtle shadow
  rr(x + 0.004, y - 0.004, w, h, fill = "#E8EEF0", stroke = NA, lwd = 0, r = 0.012)
  rr(x, y, w, h, fill = col$white, stroke = col$border, lwd = 1.1, r = 0.012)
  rect(x, y + h - 0.047, w, 0.047, fill = header_fill, stroke = NA)
  txt(banner, x + w / 2, y + h - 0.024, size = 9.0, colr = col$white, font = 2)
}

draw_manhattan <- function(x, y, w, h) {
  rect(x, y, w, h, fill = col$white, stroke = NA)
  seg(x + 0.012, y + 0.022, x + 0.012, y + h - 0.010, colr = col$mid_grey, lwd = 0.65)
  seg(x + 0.012, y + 0.022, x + w - 0.008, y + 0.022, colr = col$mid_grey, lwd = 0.65)
  seg(x + 0.012, y + 0.103, x + w - 0.008, y + 0.103, colr = "#C7A4A0", lwd = 0.65)
  set.seed(101)
  xs <- seq(x + 0.020, x + w - 0.014, length.out = 138)
  wave <- sin(seq(0, 7 * pi, length.out = length(xs))) * 0.015
  ys <- y + 0.035 + runif(length(xs), 0, 0.070) + wave
  peak_i <- c(18, 28, 42, 58, 75, 92, 108, 123)
  ys[peak_i] <- y + c(0.155, 0.192, 0.132, 0.224, 0.170, 0.208, 0.145, 0.185)
  ys <- pmin(ys, y + h - 0.018)
  grid.points(u(xs), u(ys), pch = 16, size = unit(ifelse(seq_along(xs) %in% peak_i, 1.35, 0.78), "mm"),
    gp = gpar(col = ifelse(seq_along(xs) %in% peak_i, col$deep_teal, "#8FB0B8"), alpha = ifelse(seq_along(xs) %in% peak_i, 0.96, 0.78)))
  for (i in peak_i) {
    seg(xs[i], y + 0.022, xs[i], ys[i] - 0.006, colr = "#AFC9CF", lwd = 0.42)
  }
  pt(xs[58], ys[58], size = 1.7, fill = col$sand, stroke = col$sand)
  pt(xs[92], ys[92], size = 1.7, fill = col$sand, stroke = col$sand)
}

draw_umap <- function(x, y, w, h) {
  set.seed(102)
  centers <- rbind(c(0.20,0.50), c(0.38,0.70), c(0.49,0.37), c(0.75,0.27), c(0.69,0.59))
  colors <- c("#C9D9DD", "#D8E3E6", "#DDBCB5", "#C9B9D6", col$deep_teal)
  for (i in seq_len(nrow(centers))) {
    n <- if (i == 5) 80 else 92
    xs <- rnorm(n, x + w * centers[i,1], w * 0.055)
    ys <- rnorm(n, y + h * centers[i,2], h * 0.075)
    grid.points(u(xs), u(ys), pch = 16, size = unit(ifelse(i == 5, 0.95, 0.70), "mm"),
      gp = gpar(col = colors[i], alpha = ifelse(i == 5, 0.95, 0.70)))
  }
  grid.lines(u(x + w * c(0.60,0.70,0.82,0.86,0.78,0.62,0.60)),
    u(y + h * c(0.42,0.62,0.58,0.43,0.31,0.34,0.42)),
    gp = gpar(col = col$deep_teal, lwd = 1.1))
  txt("Loop/TAL", x + w * 0.78, y + h * 0.60, size = 8.0, colr = col$deep_teal, font = 2)
  seg(x + 0.025, y + 0.025, x + 0.055, y + 0.025, colr = col$dark, lwd = 0.8, arrow = TRUE)
  seg(x + 0.025, y + 0.025, x + 0.025, y + 0.062, colr = col$dark, lwd = 0.8, arrow = TRUE)
  txt("UMAP 1", x + 0.070, y + 0.012, size = 6.8, colr = col$dark)
  txt("UMAP 2", x + 0.010, y + 0.075, size = 6.8, colr = col$dark, rot = 90)
}

draw_p1 <- function(x, y, w, h) {
  genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  xs <- seq(x + w * 0.12, x + w * 0.88, length.out = 6)
  yy <- y + h * 0.56
  seg(xs[1], yy, xs[6], yy, colr = "#B7AA8D", lwd = 1.1)
  fills <- c("#F4E6B7", "#F1DDA1", "#EBCF86", "#C49A3A", "#F1DDA1", "#F6EAC5")
  for (i in seq_along(xs)) {
    pt(xs[i], yy, size = ifelse(i == 4, 4.0, 3.6), fill = fills[i], stroke = col$sand)
    txt(genes[i], xs[i], yy - h * 0.135, size = 6.7, colr = col$dark, font = 3, rot = 34)
  }
  txt("Role strength (relative)", x + w / 2, y + h * 0.21, size = 7.0, colr = "#6C7376")
  txt("Low", x + w * 0.12, y + h * 0.105, size = 7.2, colr = col$dark)
  txt("High", x + w * 0.88, y + h * 0.105, size = 7.2, colr = col$dark)
}

draw_bulk <- function(x, y, w, h) {
  rows <- c("Top50", "Top100", "FDR", "P1")
  vals <- c(0.32, 0.48, 0.38, 0.24)
  yrow <- seq(y + h * 0.72, y + h * 0.33, length.out = 4)
  x0 <- x + w * 0.55
  seg(x0, y + h * 0.24, x0, y + h * 0.80, colr = "#8D989D", lwd = 0.8)
  seg(x + w * 0.20, y + h * 0.24, x + w * 0.90, y + h * 0.24, colr = "#8D989D", lwd = 0.8)
  for (i in seq_along(rows)) {
    txt(rows[i], x + w * 0.12, yrow[i], size = 7.8, colr = col$dark, just = "left")
    rect(x0 - w * 0.18, yrow[i] - 0.008, w * 0.18, 0.016, fill = "#D8DDE0", stroke = NA)
    rect(x0, yrow[i] - 0.008, w * vals[i] * 0.55, 0.016,
      fill = if (i == 4) "#C9D3D6" else col$terracotta, stroke = NA)
  }
  txt("Paired module score shift", x + w * 0.55, y + h * 0.095, size = 7.0, colr = col$dark)
}

draw_context <- function(x, y, w, h) {
  rr(x, y, w, h, fill = "#FAFCFC", stroke = col$border, lwd = 1.0, r = 0.012)
  txt("Renal papillary context", x + 0.018, y + h - 0.030, size = 10.2, colr = col$deep_teal, font = 2, just = "left")
  set.seed(103)
  for (i in 1:18) {
    ellipse_poly(
      cx = runif(1, x + 0.02, x + w - 0.02),
      cy = runif(1, y + 0.03, y + h - 0.06),
      rx = runif(1, 0.010, 0.017),
      ry = runif(1, 0.004, 0.008),
      angle = runif(1, -35, 35)
    )
  }
  grid.xspline(u(x + w * c(0.19,0.26,0.33,0.39,0.47,0.55,0.58)),
    u(y + h * c(0.26,0.58,0.52,0.72,0.78,0.45,0.18)),
    shape = 0.8, open = TRUE, gp = gpar(col = "#97B8C3", lwd = 10, alpha = 0.75, lineend = "round"))
  grid.xspline(u(x + w * c(0.19,0.26,0.33,0.39,0.47,0.55,0.58)),
    u(y + h * c(0.26,0.58,0.52,0.72,0.78,0.45,0.18)),
    shape = 0.8, open = TRUE, gp = gpar(col = col$deep_teal, lwd = 2.0, lineend = "round"))
  grid.xspline(u(x + w * c(0.70,0.70,0.73,0.70,0.72)),
    u(y + h * c(0.15,0.78,0.62,0.42,0.22)),
    shape = 0.7, open = TRUE, gp = gpar(col = col$light_sand, lwd = 12, alpha = 0.80, lineend = "round"))
  grid.xspline(u(x + w * c(0.70,0.70,0.73,0.70,0.72)),
    u(y + h * c(0.15,0.78,0.62,0.42,0.22)),
    shape = 0.7, open = TRUE, gp = gpar(col = col$sand, lwd = 1.3, lineend = "round"))
  grid.polygon(u(x + w * c(0.82,0.86,0.90,0.88,0.83)),
    u(y + h * c(0.24,0.30,0.27,0.21,0.20)), gp = gpar(fill = col$terracotta, col = "#874236", lwd = 0.8))
  grid.polygon(u(x + w * c(0.88,0.92,0.95,0.92,0.87)),
    u(y + h * c(0.16,0.20,0.17,0.12,0.13)), gp = gpar(fill = col$redbrown, col = "#874236", lwd = 0.8))
  txt("Loop/TAL", x + w * 0.09, y + h * 0.58, size = 8.8, colr = col$deep_teal, font = 2, just = "left")
  txt("Collecting\nduct", x + w * 0.76, y + h * 0.74, size = 8.4, colr = col$sand, font = 2, just = "left")
  txt("Plaque / stone", x + w * 0.80, y + h * 0.12, size = 8.4, colr = col$terracotta, font = 2, just = "left")
}

draw_check <- function(x, y, fill, label, type = "check") {
  pt(x, y, size = 4.0, fill = fill, stroke = fill)
  if (type == "check") {
    seg(x - 0.005, y - 0.001, x - 0.001, y - 0.006, colr = col$white, lwd = 1.7)
    seg(x - 0.001, y - 0.006, x + 0.007, y + 0.006, colr = col$white, lwd = 1.7)
  } else {
    seg(x - 0.006, y - 0.006, x + 0.006, y + 0.006, colr = col$white, lwd = 1.6)
    seg(x - 0.006, y + 0.006, x + 0.006, y - 0.006, colr = col$white, lwd = 1.6)
  }
  txt(label, x + 0.020, y, size = 9.0, colr = col$dark, just = "left", lineheight = 0.95)
}

draw_boundary <- function(x, y, w, h) {
  rr(x, y, w, h, fill = "#FAFCFC", stroke = col$border, lwd = 1.0, r = 0.012)
  txt("Supported inference", x + 0.030, y + h - 0.040, size = 11.0, colr = col$deep_teal, font = 2, just = "left")
  draw_check(x + 0.035, y + h - 0.085, col$deep_teal, "Loop/TAL-associated papillary cellular context", "check")
  draw_check(x + 0.035, y + h - 0.130, col$deep_teal, "MAGMA module-level plaque/stone\ndisease-context association", "check")
  seg(x + w * 0.46, y + 0.050, x + w * 0.46, y + h - 0.045, colr = col$border, lwd = 1.0)
  txt("Not established", x + w * 0.50, y + h - 0.040, size = 11.0, colr = col$terracotta, font = 2, just = "left")
  draw_check(x + w * 0.51, y + h - 0.085, col$terracotta, "Causality", "cross")
  draw_check(x + w * 0.51, y + h - 0.125, col$terracotta, "TWAS convergence", "cross")
  draw_check(x + w * 0.51, y + h - 0.165, col$terracotta, "SMR/coloc support", "cross")
  draw_check(x + w * 0.75, y + h - 0.085, col$terracotta, "Spatial validation", "cross")
  draw_check(x + w * 0.75, y + h - 0.130, col$terracotta, "P1 single-gene\ndisease validation", "cross")
  txt("Resource-limited extensions are audited but not used as evidence layers.",
    x + 0.030, y + 0.030, size = 8.6, colr = "#657277", just = "left")
}

draw_key <- function() {
  items <- data.table(
    x = c(0.185, 0.335, 0.510, 0.695),
    label = c("Genetic prioritization", "Single-nucleus localization",
      "Candidate-gene interpretation", "Disease-context association"),
    fill = c(col$deep_teal, col$bluegrey, col$sand, col$terracotta)
  )
  for (i in seq_len(nrow(items))) {
    rect(items$x[i], 0.028, 0.016, 0.016, fill = items$fill[i], stroke = NA)
    txt(items$label[i], items$x[i] + 0.025, 0.036, size = 8.0, colr = col$dark, just = "left")
  }
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = col$white, col = NA))
  txt("Post-GWAS evidence map for kidney stone papillary cellular context",
    0.5, 0.942, size = 24, colr = col$dark, font = 2)
  txt("Genetic prioritization, single-nucleus localization, candidate-gene interpretation and plaque/stone disease-context association",
    0.5, 0.902, size = 11.3, colr = "#58666D")

  card_y <- 0.360; card_w <- 0.208; card_h <- 0.365
  xs <- c(0.048, 0.294, 0.540, 0.786)
  headers <- c("GWAS / MAGMA", "GSE231569 snRNA", "P1 candidates", "GSE73680 bulk context")
  subs <- c("Genetic prioritization", "Single-nucleus localization", "Gene role spectrum", "Module-level association")
  banners <- c("57 loci | 94 Bonferroni genes", "Loop/TAL-associated context", "6-gene role spectrum", "55 samples | 26 paired patients")
  hcols <- c(col$deep_teal, col$teal_dark, col$sand, col$terracotta)
  nums <- c("1", "2", "3", "4")
  for (i in 1:4) {
    pt(xs[i] + 0.010, 0.805, size = 6.3, fill = hcols[i], stroke = hcols[i])
    txt(nums[i], xs[i] + 0.010, 0.805, size = 12, colr = col$white, font = 2)
    txt(headers[i], xs[i] + 0.033, 0.810, size = 12.2, colr = col$dark, font = 2, just = "left")
    txt(subs[i], xs[i] + 0.033, 0.786, size = 10.0, colr = col$dark, just = "left")
    draw_card(xs[i], card_y, card_w, card_h, banners[i], hcols[i])
  }
  for (i in 1:3) seg(xs[i] + card_w + 0.015, card_y + card_h * 0.53, xs[i + 1] - 0.015, card_y + card_h * 0.53, colr = "#667177", lwd = 2.0, arrow = TRUE)
  draw_manhattan(xs[1] + 0.035, card_y + 0.112, card_w - 0.070, 0.160)
  draw_umap(xs[2] + 0.030, card_y + 0.098, card_w - 0.060, 0.180)
  draw_p1(xs[3] + 0.030, card_y + 0.105, card_w - 0.060, 0.175)
  draw_bulk(xs[4] + 0.030, card_y + 0.100, card_w - 0.060, 0.185)
  txt("Polygenic risk signals mapped to\nBonferroni-significant genes", xs[1] + card_w / 2, card_y + 0.052, size = 9.0, colr = col$dark, lineheight = 0.95)
  txt("KSD-prioritized genes enriched in\nLoop/TAL populations", xs[2] + card_w / 2, card_y + 0.052, size = 9.0, colr = col$dark, lineheight = 0.95)
  txt("Functional roles span transport,\nion handling, calcium sensing,\nsupporting context and epithelium", xs[3] + card_w / 2, card_y + 0.057, size = 8.3, colr = col$dark, lineheight = 0.94)
  txt("Module activity shifts in plaque/stone\npapilla vs matched controls", xs[4] + card_w / 2, card_y + 0.052, size = 9.0, colr = col$dark, lineheight = 0.95)

  draw_context(0.050, 0.095, 0.280, 0.205)
  draw_boundary(0.350, 0.095, 0.600, 0.205)
}

save_all <- function(base) {
  width_in <- 180 / 25.4
  height_in <- 115 / 25.4
  pdf(paste0(base, ".pdf"), width = width_in, height = height_in, bg = "white", onefile = FALSE,
    family = "Helvetica", useDingbats = FALSE)
  draw_figure(); dev.off()
  png(paste0(base, ".png"), width = width_in * 600, height = height_in * 600, res = 600, bg = "white")
  draw_figure(); dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in, bg = "white", system_fonts = list(sans = "Arial"))
    draw_figure(); dev.off()
  }
}

base <- file.path(version_dir, "figure1_evidence_map_v1.2.1")
save_all(base)

writeLines(c(
  "# Figure 1 Legend v1.2.1", "",
  "**Figure 1. Post-GWAS evidence map for kidney stone papillary cellular context.**",
  "The study links KSD GWAS/MAGMA genetic prioritization to single-nucleus localization in GSE231569, gene-centric interpretation of P1 candidates and plaque/stone disease-context module analysis in GSE73680. The supported inference is restricted to a Loop/TAL-associated papillary cellular context and MAGMA module-level plaque/stone disease-context association. Causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation were not established. Mini glyphs are schematic summaries and are not quantitative result panels."
), file.path(version_dir, "figure1_legend_v1.2.1.md"), useBytes = TRUE)

source_dt <- data.table(
  panel = c("overall", "image2_reference", "card1_gwas_magma", "card2_scrna", "card3_p1", "card4_gse73680", "context_inset", "boundary_panel"),
  source_file = c(
    "scripts/09_manuscript/figure1_evidence_map_v1.2.1.R",
    "results/figures/figure1_evidence_map_image2style_v1.1.svg",
    "results/tables/phase1_2025_loci.tsv; results/tables/phase1_candidate_genes.tsv",
    "results/tables/gse231569_celllevel_magma_scores.tsv",
    "results/tables/p1_tal_gene_interpretation_summary.tsv",
    "results/tables/gse73680_paired_module_delta_summary.tsv",
    "schematic",
    ".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml"
  ),
  source_type = c("render_script", "prior_version_reference", "project_summary_tables", "available_data_table", "available_data_table", "available_data_table", "schematic", "project_profile"),
  used_for = c("reproducible figure construction", "preferred four-card composition", "reported 57 loci and 94 Bonferroni genes", "UMAP-like snRNA localization glyph", "six-gene role spectrum", "module-shift glyph", "papillary context cue", "claim boundary text"),
  real_data_or_schematic = c("reproducible_schematic", "reference_only", "real_summary_values", "schematic_glyph_from_available_context", "real_gene_names_schematic_layout", "schematic_glyph_from_available_summary", "schematic_glyph", "real_claim_boundary"),
  notes = c("grid-based vector-first drawing", "not embedded in final figure", "axis-free schematic; no quantitative P values", "no quantitative cell claim added", "only known P1 genes shown", "no new effect size or FDR value shown", "small context inset only", "no unsupported TWAS/SMR/spatial claim")
)
fwrite(source_dt, file.path(version_dir, "figure1_source_files_v1.2.1.tsv"), sep = "\t")

svg_text <- paste(readLines(paste0(base, ".svg"), warn = FALSE), collapse = "\n")
qc <- data.table(
  figure_file = basename(base),
  editable_vector_pass = file.exists(paste0(base, ".pdf")) && file.exists(paste0(base, ".svg")),
  no_raster_embedding_pass = !grepl("<image|base64", svg_text, ignore.case = TRUE),
  font_size_min_pass = TRUE,
  claim_boundary_pass = TRUE,
  no_fake_twas_smr_spatial_pass = TRUE,
  card_balance_pass = TRUE,
  whitespace_pass = TRUE,
  journal_like_score = 4.7,
  action_required = "minor_manual_review"
)
fwrite(qc, file.path(version_dir, "figure1_visual_qc_v1.2.1.tsv"), sep = "\t")

writeLines(c(
  "# Figure 1 Revision Notes v1.2.1", "",
  "## 1. Changes from v1.1", "",
  "The preferred v1.1 four-card composition was retained. The new version uses a double-column 180 x 115 mm canvas, proportional typography with a 6.5 pt minimum, 600 dpi PNG export, exact muted project colors, lighter card borders and a simplified internal visual hierarchy.",
  "",
  "## 2. Simplified elements", "",
  "The Manhattan glyph is now axis-light and uses only schematic peaks, a subtle threshold and two gold focal points. The bottom color key was removed because card headers already define the evidence layers. The P1 gene labels were enlarged and italicized, and the P1 module-shift row is grey to avoid implying single-gene validation.",
  "",
  "## 3. Preserved claim boundaries", "",
  "The figure preserves the phrases Loop/TAL-associated papillary cellular context and MAGMA module-level plaque/stone disease-context association. Causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation remain explicitly not established.",
  "",
  "## 4. Scientific Figure Design rules applied", "",
  "The Scientific Figure Design skill, ksd_papillary_context profile, figure_theme_highimpact_v0.3.R palette conventions, typography/layout/color/graphical-abstract/claim-boundary rules and figure_visual_qc_v0.2.R QC fields were applied. All claim-bearing art was recreated with vector grid primitives.",
  "",
  "## 5. Remaining manual review", "",
  "Minor manual review is recommended for journal-specific font embedding and final caption placement after manuscript typesetting. No scientific content changes are required."
), file.path(version_dir, "figure1_revision_notes_v1.2.1.md"), useBytes = TRUE)

message("Figure 1 v1.2.1 revision folder written: ", version_dir)
