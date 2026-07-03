suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.2.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

toy <- data.table(
  module = factor(c("MAGMA top 50", "MAGMA top 100", "P1 core"),
    levels = c("MAGMA top 50", "MAGMA top 100", "P1 core")),
  score = c(0.92, 0.84, 0.54),
  message = c("Loop/TAL context", "Module support", "Interpretive boundary")
)

p <- ggplot(toy, aes(score, module, fill = module)) +
  geom_col(width = 0.62, color = "#5A6062", linewidth = 0.25) +
  geom_text(aes(label = message), hjust = 1.05, color = "white", fontface = "bold", size = 3.1) +
  scale_fill_project(values = c(
    "MAGMA top 50" = project_palette[["deep_teal"]],
    "MAGMA top 100" = "#47757E",
    "P1 core" = project_palette[["sand_gold"]]
  )) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Phase 16 toy demo: project figure style", x = "Toy support score", y = NULL, fill = NULL) +
  theme_highimpact(9) +
  theme(legend.position = "none")

save_figure_pdf_png(p, "results/figures/phase16_toy_demo_figure", width = 6.8, height = 3.4)

qc <- figure_qc_check(
  figure_id = "Phase 16 toy demo",
  panel_id = "all",
  pdf_path = "results/figures/phase16_toy_demo_figure.pdf",
  png_path = "results/figures/phase16_toy_demo_figure.png",
  notes = "Toy figure validates theme_highimpact, scale_fill_project and save_figure_pdf_png"
)
fwrite(as.data.table(qc), "results/tables/phase16_toy_demo_visual_qc.tsv", sep = "\t")
message("Wrote toy demo figure and QC table")
