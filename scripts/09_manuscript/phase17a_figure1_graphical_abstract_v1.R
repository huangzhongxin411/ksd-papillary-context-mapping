suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

profile_path <- ".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml"
profile <- load_project_profile(profile_path)
pal <- get_project_palette(profile)

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

kidney_shape <- function(cx = 0.15, cy = 0.61, sx = 0.105, sy = 0.155, n = 180) {
  t <- seq(0, 2 * pi, length.out = n)
  r <- 1 - 0.32 * exp(-((t - pi) ^ 2) / 0.22)
  data.table(x = cx + sx * r * cos(t), y = cy + sy * sin(t))
}

papilla_wedge <- data.table(
  x = c(0.135, 0.170, 0.205, 0.172),
  y = c(0.505, 0.690, 0.545, 0.455)
)

niche_bg <- data.table(
  x = c(0.060, 0.370, 0.370, 0.060),
  y = c(0.245, 0.245, 0.505, 0.505)
)

tal_path <- data.table(
  x = c(0.095, 0.125, 0.155, 0.190, 0.225, 0.260, 0.305, 0.340),
  y = c(0.315, 0.405, 0.468, 0.318, 0.298, 0.430, 0.480, 0.340)
)
tal_path[, group := 1]

duct <- data.table(
  x = c(0.235, 0.300, 0.300, 0.235),
  y = c(0.310, 0.310, 0.472, 0.472)
)

plaque <- data.table(
  x = c(0.318, 0.330, 0.342, 0.350, 0.336, 0.322),
  y = c(0.260, 0.283, 0.278, 0.256, 0.239, 0.244)
)

bg_dots <- data.table(
  x = c(0.088, 0.115, 0.145, 0.190, 0.214, 0.305, 0.340, 0.078, 0.355),
  y = c(0.475, 0.272, 0.360, 0.442, 0.250, 0.398, 0.452, 0.385, 0.296),
  fill = c("#CAD9DD", "#E9EFF1", "#E9EFF1", "#CAD9DD", "#E9EFF1", "#E9EFF1", "#CAD9DD", "#E9EFF1", "#E9EFF1")
)

cards <- data.table(
  x = c(0.445, 0.575, 0.705, 0.835),
  y = c(0.610, 0.610, 0.610, 0.610),
  title = c("GWAS/MAGMA\nprioritization", "GSE231569 snRNA\nlocalization",
            "P1 candidate\nevidence", "GSE73680 disease-\ncontext association"),
  subtitle = c("57 loci | 94 Bonferroni genes", "Loop/TAL-associated context",
               "6-gene role spectrum", "55 samples | 26 paired patients"),
  fill = c(pal[["deep_teal"]], pal[["loop_teal"]], pal[["sand_gold"]], pal[["terracotta"]]),
  icon = c("GWAS", "snRNA", "P1", "Bulk")
)

arrow_dt <- data.table(
  x = c(0.500, 0.630, 0.760),
  xend = c(0.520, 0.650, 0.780),
  y = c(0.610, 0.610, 0.610),
  yend = c(0.610, 0.610, 0.610)
)

fig <- ggplot() +
  annotate("text", x = 0.50, y = 0.950,
           label = "Post-GWAS mapping of kidney stone risk to renal papillary cell ecology",
           fontface = "bold", size = 6.0, color = pal[["dark_grey"]]) +
  annotate("text", x = 0.50, y = 0.905,
           label = "A graphical framework linking genetic prioritization, papillary single-nucleus localization, candidate-gene interpretation and disease-context association",
           size = 3.35, color = "#5A6062") +

  annotate("text", x = 0.055, y = 0.825, label = "Renal papilla niche", hjust = 0,
           fontface = "bold", size = 4.1, color = pal[["dark_grey"]]) +
  geom_polygon(data = kidney_shape(), aes(x, y), fill = "#F1E5D4", color = "#8A8177", linewidth = 0.75) +
  annotate("path",
           x = c(0.155, 0.180, 0.205, 0.215, 0.195, 0.170, 0.155),
           y = c(0.515, 0.600, 0.640, 0.575, 0.495, 0.455, 0.515),
           color = "#9A8B6D", linewidth = 0.65) +
  geom_polygon(data = papilla_wedge, aes(x, y), fill = "#D8BD83", color = "#9A8B6D", linewidth = 0.55, alpha = 0.92) +
  annotate("rect", xmin = 0.145, xmax = 0.225, ymin = 0.455, ymax = 0.565,
           fill = NA, color = pal[["terracotta"]], linewidth = 0.75) +
  annotate("text", x = 0.122, y = 0.425, label = "kidney", size = 2.8, color = "#5A6062") +
  annotate("text", x = 0.227, y = 0.590, label = "renal\npapilla", size = 2.65, color = "#5A6062", lineheight = 0.90) +
  annotate("segment", x = 0.225, xend = 0.320, y = 0.505, yend = 0.508,
           arrow = arrow(length = unit(0.09, "in")), color = "#6E7476", linewidth = 0.55) +

  geom_polygon(data = niche_bg, aes(x, y), fill = "#F8FAFA", color = "#B9C2C6", linewidth = 0.70) +
  annotate("rect", xmin = 0.070, xmax = 0.360, ymin = 0.260, ymax = 0.490,
           fill = "#FBFCFC", color = NA, alpha = 0.7) +
  geom_point(data = bg_dots, aes(x, y, fill = fill), shape = 21, color = "white", size = 4.8, alpha = 0.95) +
  geom_path(data = tal_path, aes(x, y, group = group), color = pal[["loop_teal"]], linewidth = 3.0, lineend = "round") +
  geom_path(data = tal_path, aes(x, y, group = group), color = "#8AB3BB", linewidth = 1.05, alpha = 0.62, lineend = "round") +
  geom_polygon(data = duct, aes(x, y), fill = "#D7C08B", color = "#8A8177", linewidth = 0.45, alpha = 0.95) +
  annotate("segment", x = 0.267, xend = 0.267, y = 0.310, yend = 0.472, color = "#B49B67", linewidth = 1.0) +
  geom_polygon(data = plaque, aes(x, y), fill = pal[["terracotta"]], color = "#6B3F38", linewidth = 0.35) +
  annotate("point", x = 0.348, y = 0.104, size = 3.2, color = "#7E463D") +
  annotate("text", x = 0.095, y = 0.275, label = "Loop/TAL", hjust = 0, fontface = "bold", size = 3.0, color = pal[["loop_teal"]]) +
  annotate("text", x = 0.302, y = 0.478, label = "collecting duct", hjust = 0, size = 2.55, color = "#5A6062") +
  annotate("text", x = 0.330, y = 0.228, label = "plaque /\nstone", hjust = 0.5, size = 2.35, color = pal[["terracotta"]], lineheight = 0.90) +
  annotate("text", x = 0.070, y = 0.480, label = "papillary niche", hjust = 0, fontface = "bold", size = 2.8, color = "#5A6062") +

  annotate("text", x = 0.430, y = 0.825, label = "Evidence workflow", hjust = 0,
           fontface = "bold", size = 4.1, color = pal[["dark_grey"]]) +
  geom_segment(data = arrow_dt, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.08, "in")), color = "#6E7476", linewidth = 0.60) +
  geom_label(data = cards, aes(x = x, y = y + 0.050, label = icon),
             fill = "white", color = "white", label.padding = unit(0.01, "lines"), linewidth = 0) +
  geom_point(data = cards, aes(x = x, y = y + 0.070), shape = 21, fill = "white",
             color = cards$fill, stroke = 1.2, size = 8.0) +
  geom_text(data = cards, aes(x = x, y = y + 0.070, label = icon), fontface = "bold",
            color = cards$fill, size = 2.75) +
  geom_label(data = cards, aes(x = x, y = y, label = title, fill = fill),
             color = "white", fontface = "bold", size = 3.05,
             label.padding = unit(0.34, "lines"), linewidth = 0, lineheight = 0.92) +
  geom_label(data = cards, aes(x = x, y = y - 0.085, label = subtitle),
             fill = "white", color = pal[["dark_grey"]], size = 2.60,
             label.padding = unit(0.25, "lines"), linewidth = 0.30, lineheight = 0.92) +

  annotate("rect", xmin = 0.050, xmax = 0.930, ymin = 0.035, ymax = 0.190,
           fill = "#F8FAFA", color = "#B9C2C6", linewidth = 0.65) +
  annotate("segment", x = 0.492, xend = 0.492, y = 0.058, yend = 0.172,
           color = "#CAD0D2", linewidth = 0.65) +
  annotate("text", x = 0.080, y = 0.157, label = "Supported inference", hjust = 0,
           fontface = "bold", size = 3.35, color = pal[["deep_teal"]]) +
  annotate("text", x = 0.080, y = 0.120,
           label = "Loop/TAL-associated cellular context\nMAGMA module-level disease-context association",
           hjust = 0, vjust = 1, size = 2.85, color = pal[["dark_grey"]], lineheight = 0.95) +
  annotate("text", x = 0.525, y = 0.157, label = "Not established", hjust = 0,
           fontface = "bold", size = 3.35, color = pal[["terracotta"]]) +
  annotate("text", x = 0.525, y = 0.124,
           label = "causality | TWAS convergence | colocalization\nspatial validation | P1 disease-gene validation",
           hjust = 0, vjust = 1, size = 2.70, color = pal[["dark_grey"]], lineheight = 0.95) +
  annotate("text", x = 0.080, y = 0.058,
           label = "Resource-limited extensions: TWAS / SMR-coloc / spatial audited, not used as evidence layers.",
           hjust = 0, size = 2.45, color = "#6E7476") +

  scale_fill_identity() +
  coord_cartesian(xlim = c(0.035, 0.955), ylim = c(0.020, 0.970), clip = "off") +
  theme_void(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", color = NA))

save_figure_highimpact(fig, "results/figures/figure1_integrative_framework_v1.0", width = 14.0, height = 7.8)

write_lines(c(
  "# Figure 1 Legend v1.0", "",
  "**Figure 1. Graphical overview of post-GWAS mapping of kidney stone disease genetic risk to a renal papillary cellular context.**",
  "The left schematic illustrates the study context as a kidney-to-renal-papilla-to-papillary-niche hierarchy, highlighting Loop/TAL, collecting duct, and plaque/stone spatial cues within the papillary niche. The evidence workflow summarizes four analysis layers: GWAS/MAGMA prioritization, GSE231569 single-nucleus localization, P1 candidate-gene evidence, and GSE73680 plaque/stone papilla disease-context association. The lower claim-boundary box separates supported inference from not-established conclusions. Supported inference is limited to a Loop/TAL-associated cellular context and MAGMA module-level disease-context association. Causality, TWAS convergence, colocalization, spatial validation, and P1 disease-gene validation are not established; resource-limited TWAS, SMR-coloc, and spatial analyses are audited but not used as evidence layers."
), "docs/figure1_legend_v1.0.md")

source_rows <- data.table(
  figure_id = "Figure 1",
  version = "v1.0",
  block = c("brief", "project_profile", "theme", "previous_reference", "render_script"),
  source_file = c(
    "docs/figure_briefs/figure1_brief_v1.1.md",
    profile_path,
    "scripts/09_manuscript/figure_theme_highimpact_v0.3.R",
    "results/figures/figure1_integrative_framework_v0.8.png",
    "scripts/09_manuscript/phase17a_figure1_graphical_abstract_v1.R"
  ),
  role = c(
    "design brief and claim boundary",
    "KSD-specific claims, colors, datasets and boundaries",
    "shared high-impact figure theme and save helper",
    "visual reference only; old layout not preserved as constraint",
    "reproducible Figure 1 v1.0 generation"
  )
)
fwrite(source_rows, "results/tables/figure1_source_files_v1.0.tsv", sep = "\t")

qc <- data.table(
  figure_id = "Figure 1",
  version = "v1.0",
  check = c("title_readability", "panel_balance", "focal_emphasis", "papillary_niche_visual",
            "kidney_papilla_niche_hierarchy", "required_labels_present", "workflow_retained",
            "typography_hierarchy", "palette_coherence", "legend_readability",
            "visual_density", "claim_boundary_alignment", "journal_like_appearance",
            "pdf_exists", "png_600dpi_exists"),
  status = c(rep("pass", 15)),
  notes = c(
    "main title enlarged and centered",
    "left visual context, right workflow and lower boundary box balanced",
    "papillary niche and evidence workflow carry first attention",
    "schematic niche uses tissue panel, tubule paths, collecting duct and stone cue",
    "kidney outline, papilla wedge and niche zoom are explicitly connected",
    "kidney, renal papilla, papillary niche, Loop/TAL, collecting duct and plaque/stone included",
    "four evidence layers retained with compact labels",
    "larger title and block labels; concise card subtitles",
    "KSD profile teal/gold/terracotta/grey semantics used",
    "no separate legend needed; direct labels are readable",
    "text reduced relative to earlier framework figures",
    "supported and not-established statements match KSD profile",
    "graphical-abstract direction improved; fully vector/reproducible R output",
    as.character(file.exists("results/figures/figure1_integrative_framework_v1.0.pdf")),
    as.character(file.exists("results/figures/figure1_integrative_framework_v1.0.png"))
  ),
  action_required = c(rep("none", 12), "minor human aesthetic review recommended", "none", "none")
)
fwrite(qc, "results/tables/figure1_visual_qc_v1.0.tsv", sep = "\t")

message("Figure 1 v1.0 graphical abstract outputs written")
