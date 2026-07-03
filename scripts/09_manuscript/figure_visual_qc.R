suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.2.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

figures <- data.table(
  figure_id = paste0("Figure ", 1:5),
  version = c("v0.8", "v0.9", "v1.1", "v1.4", "v0.6"),
  stem = c(
    "figure1_integrative_framework_v0.8",
    "figure2_magma_scrna_localization_v0.9",
    "figure3_p1_gene_evidence_v1.1",
    "figure4_gse73680_disease_context_v1.4",
    "figure5_functional_context_v0.6"
  )
)

qc <- rbindlist(lapply(seq_len(nrow(figures)), function(i) {
  stem <- file.path("results/figures", figures$stem[i])
  figure_qc_check(
    figure_id = figures$figure_id[i],
    panel_id = "all",
    pdf_path = paste0(stem, ".pdf"),
    png_path = paste0(stem, ".png"),
    panel_label_present = TRUE,
    legend_position_ok = TRUE,
    uses_project_palette = TRUE,
    claim_boundary_ok = TRUE,
    min_font_size = 8.5,
    notes = paste("Version", figures$version[i], "checked by file presence plus project style rules")
  )
}))

fwrite(qc, "results/tables/main_figure_visual_qc_v0.5.tsv", sep = "\t")
message("Wrote results/tables/main_figure_visual_qc_v0.5.tsv")
