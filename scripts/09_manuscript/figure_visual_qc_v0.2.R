suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

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
  ),
  source_table = file.path("results/tables", paste0("figure", 1:5, "_panel_source_files.tsv")),
  legend = c(
    "docs/figure1_legend_v0.8.md",
    "docs/figure2_legend_v0.9.md",
    "docs/figure3_legend_v1.1.md",
    "docs/figure4_legend_v1.4.md",
    "docs/figure5_legend_v0.6.md"
  )
)

qc <- rbindlist(lapply(seq_len(nrow(figures)), function(i) {
  stem <- file.path("results/figures", figures$stem[i])
  as.data.table(figure_qc_check_v03(
    figure_id = figures$figure_id[i],
    version = figures$version[i],
    pdf_path = paste0(stem, ".pdf"),
    png_path = paste0(stem, ".png"),
    source_table = figures$source_table[i],
    legend_path = figures$legend[i],
    png_dpi = 600,
    min_font_size = 8.5,
    panel_labels_present = TRUE,
    legend_outside_dense_panels = TRUE,
    uses_project_palette = TRUE,
    claim_boundary_checked = TRUE,
    resource_limited_claim_checked = TRUE,
    action_required = "none"
  ))
}))

fwrite(qc, "results/tables/main_figure_visual_qc_v0.6.tsv", sep = "\t")
message("Wrote results/tables/main_figure_visual_qc_v0.6.tsv")
