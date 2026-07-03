#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "results/figures/revision/stage7E_packaging_and_figure_refinement")
table_dir <- file.path(root, "results/tables/revision/stage7E_packaging_and_figure_refinement")
log_dir <- file.path(root, "logs/revision/stage7E_packaging_and_figure_refinement")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stem <- file.path(out_dir, "figure1_evidence_hierarchy_draft_v0.2")
font_family <- "Helvetica"

pal <- c(
  navy = "#245A64",
  teal = "#0F6B73",
  bluegrey = "#7F9DA6",
  gold = "#B99B5A",
  coral = "#A76655",
  pale_blue = "#E8F0F2",
  pale_teal = "#E4F1F0",
  pale_gold = "#F4EEDC",
  pale_coral = "#F4E8E4",
  pale_grey = "#F1F3F3",
  mid_grey = "#B8C0C2",
  dark = "#2D3436",
  grey = "#626B6E",
  white = "#FFFFFF"
)

gp_text <- function(size = 6.1, col = pal[["dark"]], face = "plain") {
  gpar(fontfamily = font_family, fontsize = size, col = col, fontface = face, lineheight = 0.95)
}

add_text <- function(label, x, y, size = 6.1, col = pal[["dark"]], face = "plain",
                     hjust = 0.5, vjust = 0.5, rot = 0) {
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
            just = c(hjust, vjust), rot = rot,
            gp = gp_text(size, col, face))
}

add_rect <- function(x, y, w, h, fill = pal[["white"]], col = pal[["mid_grey"]],
                     lwd = 0.8, lty = 1, radius = 0.008) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(radius, "npc"),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
}

add_line <- function(x0, y0, x1, y1, col = pal[["mid_grey"]], lwd = 0.8,
                     lty = 1, arrow_end = FALSE) {
  grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
             gp = gpar(col = col, lwd = lwd, lty = lty),
             arrow = if (arrow_end) arrow(length = unit(0.07, "inches"), type = "closed") else NULL)
}

panel_box <- function(x, y, w, h, letter, title, evidence, border, fill, lty = 1) {
  add_rect(x + w / 2, y + h / 2, w, h, fill = fill, col = border, lwd = 1.05, lty = lty)
  add_text(letter, x + 0.011, y + h - 0.014, size = 7.8, face = "bold", hjust = 0, vjust = 1)
  add_text(title, x + 0.038, y + h - 0.014, size = 7.0, face = "bold", hjust = 0, vjust = 1)
  add_rect(x + w - 0.078, y + h - 0.055, 0.140, 0.024,
           fill = adjustcolor(border, alpha.f = 0.10), col = border, lwd = 0.6, radius = 0.005)
  add_text(evidence, x + w - 0.078, y + h - 0.055, size = 4.25,
           col = border, face = "bold")
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  add_text("Evidence-stratified renal papillary context mapping of KSD genetic risk",
           0.025, 0.965, size = 10.6, face = "bold", hjust = 0)
  add_text("Genetic prioritization is linked to parallel, graded papillary context evidence - not cumulative proof",
           0.025, 0.928, size = 6.2, col = pal[["grey"]], hjust = 0)

  # Panel A: primary genetic prioritization.
  ax <- 0.022; ay <- 0.595; aw <- 0.245; ah <- 0.285
  panel_box(ax, ay, aw, ah, "A", "GWAS to MAGMA", "PRIMARY PRIORITIZATION",
            pal[["navy"]], pal[["pale_blue"]])
  # chromosome/locus motif
  for (i in 0:4) {
    x <- ax + 0.035 + i * 0.027
    add_line(x, ay + 0.148, x, ay + 0.197, col = pal[["bluegrey"]], lwd = 2.2)
    grid.circle(x = unit(x, "npc"), y = unit(ay + 0.170 + ifelse(i %% 2 == 0, 0.013, -0.008), "npc"),
                r = unit(0.006, "npc"), gp = gpar(fill = pal[["gold"]], col = NA))
  }
  add_text("KSD GWAS", ax + 0.084, ay + 0.123, size = 6.5, face = "bold")
  add_line(ax + 0.145, ay + 0.163, ax + 0.183, ay + 0.163, col = pal[["navy"]], lwd = 1.0, arrow_end = TRUE)
  add_rect(ax + 0.202, ay + 0.163, 0.055, 0.082, fill = "white", col = pal[["navy"]], lwd = 0.8)
  add_text("MAGMA", ax + 0.202, ay + 0.177, size = 6.2, face = "bold")
  add_text("gene sets", ax + 0.202, ay + 0.151, size = 5.2, col = pal[["grey"]])
  add_text("GWAS QC", ax + 0.025, ay + 0.103, size = 5.8, face = "bold", hjust = 0)
  add_text("4,915,033 QC rows", ax + 0.025, ay + 0.077, size = 5.7, hjust = 0)
  add_text("MAGMA output", ax + 0.142, ay + 0.105, size = 5.5, face = "bold", hjust = 0)
  add_text("17,316 genes tested", ax + 0.142, ay + 0.078, size = 5.4, hjust = 0)
  add_text("94 Bonferroni genes", ax + 0.142, ay + 0.052, size = 5.4, hjust = 0)
  add_rect(ax + 0.190, ay + 0.023, 0.095, 0.042, fill = pal[["pale_gold"]], col = pal[["gold"]], lwd = 0.7)
  add_text("1000G EUR LD\nreference", ax + 0.190, ay + 0.023, size = 4.5, col = pal[["dark"]], face = "bold")

  # Panel B: evidence model. Wider hero panel.
  bx <- 0.292; by <- 0.595; bw <- 0.686; bh <- 0.285
  panel_box(bx, by, bw, bh, "B", "Two-axis candidate evidence model", "EVIDENCE SEPARATION",
            pal[["navy"]], "#FAFBFB")
  mx <- bx + 0.058; my <- by + 0.055; mw <- 0.245; mh <- 0.145
  # Matrix 3 x 3 with no rank implication.
  for (row in 0:2) {
    for (col in 0:2) {
      fill <- if (col == 0) adjustcolor(pal[["navy"]], alpha.f = 0.12 + row * 0.06) else
        if (col == 1) adjustcolor(pal[["bluegrey"]], alpha.f = 0.11 + row * 0.04) else pal[["pale_grey"]]
      grid.rect(x = unit(mx + (col + 0.5) * mw / 3, "npc"),
                y = unit(my + (row + 0.5) * mh / 3, "npc"),
                width = unit(mw / 3, "npc"), height = unit(mh / 3, "npc"),
                gp = gpar(fill = fill, col = "white", lwd = 0.7))
    }
  }
  add_text("MAGMA only", mx + mw / 6, my + mh / 2, size = 4.5, face = "bold")
  add_text("MAGMA +\nproxy TWAS", mx + mw / 2, my + mh / 2, size = 4.3, face = "bold")
  add_text("TWAS proxy\nonly", mx + 5 * mw / 6, my + mh / 2, size = 4.3, face = "bold")
  add_text("MAGMA priority", mx + mw / 2, my - 0.020, size = 5.4, face = "bold")
  add_text("Kidney_Cortex\nTWAS proxy", mx - 0.025, my + mh / 2, size = 5.1,
           face = "bold", rot = 90)
  add_text("GENETIC PRIORITY", mx + mw + 0.055, my + mh - 0.015, size = 5.2,
           col = pal[["navy"]], face = "bold", hjust = 0)
  add_text("SUPPLEMENTARY PROXY", mx + mw + 0.055, my + mh - 0.050, size = 5.2,
           col = pal[["bluegrey"]], face = "bold", hjust = 0)
  add_text("42/51 FDR TWAS\nmodels are one-SNP", mx + mw + 0.055, my + mh - 0.091,
           size = 5.0, hjust = 0)
  # Curated exemplar tag as separate side flag.
  ex_x <- bx + 0.557; ex_y <- by + 0.073
  add_rect(ex_x, ex_y, 0.205, 0.092, fill = pal[["pale_gold"]], col = pal[["gold"]], lwd = 0.8)
  add_text("Curated exemplars", ex_x, ex_y + 0.021, size = 6.2, face = "bold")
  add_text("INTERPRETIVE FLAG", ex_x, ex_y - 0.010, size = 5.1, col = pal[["gold"]], face = "bold")
  add_text("Biology-informed; no evidence upgrade", ex_x, ex_y - 0.035, size = 5.0, col = pal[["grey"]])
  # A to B arrow: only serial arrow in the figure.
  add_line(ax + aw + 0.006, ay + ah / 2, bx - 0.008, by + bh / 2,
           col = pal[["navy"]], lwd = 1.0, arrow_end = TRUE)

  # Branching connector into parallel context panels.
  branch_y <- 0.550
  add_line(bx + bw / 2, by - 0.005, bx + bw / 2, branch_y, col = pal[["mid_grey"]], lwd = 0.8)
  add_line(0.170, branch_y, 0.830, branch_y, col = pal[["mid_grey"]], lwd = 0.8)
  add_text("parallel context mapping", 0.500, branch_y + 0.014, size = 5.0, col = pal[["grey"]])

  # Panel C: snRNA.
  cx <- 0.022; cy <- 0.245; cw <- 0.300; ch <- 0.265
  panel_box(cx, cy, cw, ch, "C", "Donor-level snRNA context", "MODERATE CONTEXT SUPPORT",
            pal[["teal"]], pal[["pale_teal"]])
  add_line(0.170, branch_y, cx + cw / 2, cy + ch + 0.004, col = pal[["mid_grey"]], lwd = 0.8, arrow_end = TRUE)
  donor_y <- cy + 0.157
  for (i in 0:3) {
    dx <- cx + 0.046 + i * 0.047
    grid.circle(x = unit(dx, "npc"), y = unit(donor_y + 0.025, "npc"), r = unit(0.010, "npc"),
                gp = gpar(fill = pal[["bluegrey"]], col = "white", lwd = 0.6))
    add_line(dx, donor_y + 0.014, dx, donor_y - 0.020, col = pal[["bluegrey"]], lwd = 1.5)
    add_line(dx - 0.011, donor_y, dx + 0.011, donor_y, col = pal[["bluegrey"]], lwd = 1.2)
  }
  add_text("4 donors", cx + 0.115, cy + 0.105, size = 6.2, face = "bold")
  # Loop/TAL stylized U-shape.
  ux <- cx + 0.237
  grid.lines(x = unit(c(ux - 0.025, ux - 0.025, ux, ux + 0.025, ux + 0.025), "npc"),
             y = unit(c(cy + 0.185, cy + 0.126, cy + 0.100, cy + 0.126, cy + 0.185), "npc"),
             gp = gpar(col = pal[["teal"]], lwd = 3.2, linejoin = "round"))
  add_text("Loop/TAL-associated", cx + 0.150, cy + 0.065, size = 6.2, face = "bold")
  add_text("partial matched-random support", cx + 0.150, cy + 0.037, size = 6.0, col = pal[["grey"]], face = "bold")

  # Panel D: paired bulk with visual attenuation.
  dx <- 0.350; dy <- 0.245; dw <- 0.300; dh <- 0.265
  panel_box(dx, dy, dw, dh, "D", "Paired bulk disease context", "ATTENUATED BULK CONTEXT",
            pal[["coral"]], pal[["pale_coral"]])
  add_line(0.500, branch_y, dx + dw / 2, dy + dh + 0.004, col = pal[["mid_grey"]], lwd = 0.8, arrow_end = TRUE)
  add_rect(dx + 0.078, dy + 0.154, 0.105, 0.065, fill = "white", col = pal[["bluegrey"]], lwd = 0.8)
  add_text("Control / adjacent", dx + 0.078, dy + 0.164, size = 5.5, face = "bold")
  add_text("paired papilla", dx + 0.078, dy + 0.139, size = 4.9, col = pal[["grey"]])
  add_rect(dx + 0.225, dy + 0.154, 0.105, 0.065, fill = "white", col = pal[["gold"]], lwd = 0.8)
  add_text("Plaque / stone", dx + 0.225, dy + 0.164, size = 5.5, face = "bold")
  add_text("paired papilla", dx + 0.225, dy + 0.139, size = 4.9, col = pal[["grey"]])
  add_line(dx + 0.137, dy + 0.154, dx + 0.166, dy + 0.154, col = pal[["coral"]], lwd = 1.0, arrow_end = TRUE)
  add_text("26 paired patients", dx + 0.150, dy + 0.105, size = 6.2, face = "bold")
  # Attenuation bars, decreasing width.
  add_rect(dx + 0.093, dy + 0.065, 0.120, 0.022, fill = adjustcolor(pal[["coral"]], alpha.f = 0.85), col = NA, radius = 0.003)
  add_line(dx + 0.160, dy + 0.065, dx + 0.205, dy + 0.065, col = pal[["coral"]], lwd = 0.8, arrow_end = TRUE)
  add_rect(dx + 0.242, dy + 0.065, 0.065, 0.022, fill = adjustcolor(pal[["coral"]], alpha.f = 0.45), col = NA, radius = 0.003)
  add_text("injury/remodeling-associated", dx + 0.150, dy + 0.034, size = 5.5, face = "bold")
  add_text("attenuated after adjustment", dx + 0.150, dy + 0.014, size = 5.0, col = pal[["grey"]])

  # Panel E: spatial + TWAS supplementary context.
  ex <- 0.678; ey <- 0.245; ew <- 0.300; eh <- 0.265
  panel_box(ex, ey, ew, eh, "E", "Spatial / TWAS context", "SUPPLEMENTARY CONTEXT",
            pal[["bluegrey"]], "#F4F7F7", lty = 2)
  add_line(0.830, branch_y, ex + ew / 2, ey + eh + 0.004, col = pal[["mid_grey"]], lwd = 0.8, arrow_end = TRUE)
  # Five section tiles.
  for (i in 0:4) {
    sx <- ex + 0.045 + i * 0.038
    grid.circle(x = unit(sx, "npc"), y = unit(ey + 0.168, "npc"), r = unit(0.014, "npc"),
                gp = gpar(fill = adjustcolor(pal[["teal"]], alpha.f = 0.14 + i * 0.06), col = pal[["bluegrey"]], lwd = 0.6))
    for (j in 0:2) {
      grid.circle(x = unit(sx - 0.006 + j * 0.006, "npc"), y = unit(ey + 0.164 + (j %% 2) * 0.007, "npc"),
                  r = unit(0.0022, "npc"), gp = gpar(fill = pal[["teal"]], col = NA))
    }
  }
  add_text("5 sections / 7,747 spots", ex + 0.121, ey + 0.117, size = 6.0, face = "bold")
  add_text("Loop/TAL co-distribution", ex + 0.121, ey + 0.090, size = 5.4)
  add_rect(ex + 0.241, ey + 0.158, 0.070, 0.055, fill = pal[["pale_gold"]], col = pal[["gold"]], lwd = 0.7, lty = 2)
  add_text("No lesion ROI", ex + 0.241, ey + 0.158, size = 5.4, face = "bold")
  add_rect(ex + 0.205, ey + 0.055, 0.170, 0.057, fill = "white", col = pal[["bluegrey"]], lwd = 0.7)
  add_text("Kidney_Cortex TWAS proxy", ex + 0.205, ey + 0.064, size = 5.0, face = "bold")
  add_text("supplementary only", ex + 0.205, ey + 0.041, size = 4.6, col = pal[["grey"]])

  # Panel F: full-width boundary strip.
  fx <- 0.022; fy <- 0.035; fw <- 0.956; fh <- 0.155
  add_rect(fx + fw / 2, fy + fh / 2, fw, fh, fill = pal[["pale_grey"]], col = pal[["dark"]], lwd = 0.9)
  add_text("F", fx + 0.011, fy + fh - 0.014, size = 7.8, face = "bold", hjust = 0, vjust = 1)
  add_text("Integrated interpretation boundary", fx + 0.038, fy + fh - 0.014, size = 7.0, face = "bold", hjust = 0, vjust = 1)
  add_text("POST-GWAS RENAL PAPILLARY CONTEXT MAPPING", fx + fw - 0.145, fy + fh - 0.016,
           size = 4.7, col = pal[["navy"]], face = "bold")
  add_text("MAGMA-prioritized genes map to a\nLoop/TAL-associated papillary context and an\ninjury/remodeling-associated bulk background",
           fx + 0.030, fy + 0.067, size = 5.75, face = "bold", hjust = 0)
  bounds <- c("Not causal-gene\nidentification", "Not a causal\ncell type",
              "Not plaque-specific\nlocalization", "TWAS is a\nKidney_Cortex proxy")
  bcols <- c(pal[["coral"]], pal[["coral"]], pal[["coral"]], pal[["gold"]])
  for (i in seq_along(bounds)) {
    xx <- fx + 0.565 + (i - 1) * 0.098
    add_rect(xx, fy + 0.055, 0.095, 0.065,
             fill = adjustcolor(bcols[i], alpha.f = 0.10), col = bcols[i], lwd = 0.7)
    add_text(bounds[i], xx, fy + 0.055, size = 5.0, face = "bold")
  }
}

width_in <- 7.205
height_in <- width_in * 9 / 16

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in,
                  bg = "white", system_fonts = list(sans = font_family))
draw_figure()
dev.off()

grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
               family = font_family, bg = "white", useDingbats = FALSE)
draw_figure()
dev.off()

ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
              units = "in", res = 600, background = "white")
draw_figure()
dev.off()

panel_source <- data.frame(
  figure = "Figure 1 draft v0.2",
  panel = LETTERS[1:6],
  panel_role = c(
    "GWAS QC and EUR-LD-reference MAGMA prioritization",
    "Two-axis genetic-priority and Kidney_Cortex TWAS-proxy model",
    "Donor-level Loop/TAL-associated snRNA context",
    "Paired injury/remodeling-associated bulk context",
    "Supplementary spatial and TWAS context",
    "Integrated interpretation boundary"
  ),
  source_file = c(
    "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv;results/tables/revision/stage2_genetic/magma_output_audit.tsv",
    "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv;results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv;results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
    "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
    "results/tables/revision/stage5C1_gse73680_figure4_draft/gse73680_final_claim_lock_stage5C1.tsv",
    "results/tables/revision/stage6C_spatial_twas_figure5_draft/spatial_twas_final_claim_lock_stage6C.tsv",
    "results/tables/revision/stage7A_manuscript_integration_blueprint/integrated_evidence_hierarchy_v0.1.tsv"
  ),
  evidence_strength = c(
    "primary_prioritization", "mixed_roles", "moderate_context_support",
    "attenuated_bulk_context_support", "moderate_supplementary_context",
    "interpretation_boundary"
  ),
  notes = c(
    "1000G EUR LD limitation displayed.",
    "Curated exemplars displayed as an independent interpretive flag.",
    "Four donors and partial matched-random support displayed.",
    "Visual taper communicates attenuation after adjustment.",
    "No-lesion-ROI boundary and Kidney_Cortex proxy displayed.",
    "No serial validation or causal labels."
  ),
  stringsAsFactors = FALSE
)
write.table(panel_source,
            file.path(table_dir, "figure1_panel_source_manifest_v0.2.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

log_lines <- c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "backend=R grid/svglite/pdf/ragg",
  paste0("canvas_inches=", width_in, "x", round(height_in, 4)),
  "layout=16:9 schematic-led evidence hierarchy",
  "panels=A-F",
  paste0("font=", font_family),
  "version=v0.2 minor visual refinement; scientific meaning unchanged",
  "png_dpi=600",
  paste0("svg=", paste0(stem, ".svg")),
  paste0("pdf=", paste0(stem, ".pdf")),
  paste0("png=", paste0(stem, ".png"))
)
writeLines(log_lines, file.path(log_dir, "stage7E_generate_figure1_v0.2.log"))

cat(paste(log_lines, collapse = "\n"), "\n")
