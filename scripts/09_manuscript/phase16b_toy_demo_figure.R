suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

make_demo <- function(profile_path, stem, figure_id) {
  profile <- load_project_profile(profile_path)
  pal <- get_project_palette(profile)
  claim <- if (!is.null(profile$main_claim)) profile$main_claim else "Project claim from profile"
  toy <- data.table(
    evidence = factor(c("Primary signal", "Supporting evidence", "Boundary"),
      levels = c("Primary signal", "Supporting evidence", "Boundary")),
    score = c(0.92, 0.78, 0.48),
    fill = c("MAGMA top 50", "Audited context", "P1 core")
  )
  p <- ggplot(toy, aes(score, evidence, fill = fill)) +
    geom_col(width = 0.62, color = pal[["dark_grey"]], linewidth = 0.25) +
    geom_text(aes(label = evidence), hjust = 1.05, color = "white", fontface = "bold", size = 3.1) +
    scale_fill_project(profile = profile) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = paste0(figure_id, ": profile-driven figure style"),
      subtitle = substr(claim, 1, 135), x = "Toy evidence score", y = NULL, fill = NULL) +
    theme_highimpact(10) +
    theme(legend.position = "none")
  save_figure_highimpact(p, stem, width = 7.6, height = 3.8)
  qc <- as.data.table(figure_qc_check_v03(
    figure_id = figure_id,
    version = "v0.1",
    pdf_path = paste0(stem, ".pdf"),
    png_path = paste0(stem, ".png"),
    source_table = NA_character_,
    legend_path = NA_character_,
    action_required = "none"
  ))
  qc[, profile_path := profile_path]
  qc
}

qc <- rbindlist(list(
  make_demo(
    ".codex/skills/scientific_figure_design/project_profiles/generic_bioinformatics.profile.yml",
    "results/figures/phase16b_generic_toy_demo",
    "Phase 16B generic toy demo"
  ),
  make_demo(
    ".codex/skills/scientific_figure_design/project_profiles/ksd_papillary_context.profile.yml",
    "results/figures/phase16b_ksd_toy_demo",
    "Phase 16B KSD toy demo"
  )
))

fwrite(qc, "results/tables/phase16b_toy_demo_visual_qc.tsv", sep = "\t")
message("Wrote Phase 16B toy demo figures and QC")
