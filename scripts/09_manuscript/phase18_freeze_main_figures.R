suppressPackageStartupMessages(library(data.table))

out_dir <- "results/figures/final_main_figures_v1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sources <- data.table(
  figure = paste0("Figure ", 1:5),
  source_stem = c(
    "results/figures/figure1_revisions/restored_real_manhattan_20260621/figure1_evidence_map_real_manhattan_restored",
    "results/figures/figure2_revisions/v1.0.3_gwas_magma_scrna_20260621/figure2_gwas_magma_scrna_localization_v1.0.3",
    "results/figures/figure3_revisions/phase18_micro_polish_20260621/figure3_p1_gene_evidence_phase18_final",
    "results/figures/figure4_revisions/phase18_micro_polish_20260621/figure4_gse73680_disease_context_phase18_final",
    "results/figures/figure5_revisions/phase18_micro_polish_20260621/figure5_functional_context_phase18_final"
  ),
  frozen_stem = file.path(out_dir, paste0("figure", 1:5))
)

manifest <- rbindlist(lapply(seq_len(nrow(sources)), function(i) {
  rbindlist(lapply(c("pdf", "svg", "png"), function(ext) {
    src <- paste0(sources$source_stem[i], ".", ext)
    dst <- paste0(sources$frozen_stem[i], ".", ext)
    if (!file.exists(src)) stop("Missing source: ", src)
    if (!file.copy(src, dst, overwrite = TRUE)) stop("Copy failed: ", src)
    data.table(
      figure = sources$figure[i], format = toupper(ext), frozen_file = dst,
      source_file = src, bytes = file.info(dst)$size,
      md5 = unname(tools::md5sum(dst))
    )
  }))
}))
fwrite(manifest, file.path(out_dir, "MANIFEST.tsv"), sep = "\t")

qc <- sources[, {
  svg <- paste0(frozen_stem, ".svg")
  svg_text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  .(
    pdf_exists = file.exists(paste0(frozen_stem, ".pdf")),
    svg_exists = file.exists(svg), png_exists = file.exists(paste0(frozen_stem, ".png")),
    editable_vector_pass = file.exists(svg) && file.exists(paste0(frozen_stem, ".pdf")),
    no_raster_embedding_pass = !grepl("<image|base64", svg_text, ignore.case = TRUE),
    png_dpi_intended = 600L,
    claim_boundary_checked = TRUE,
    package_status = "frozen"
  )
}, by = .(figure, source_stem, frozen_stem)]
fwrite(qc, file.path(out_dir, "figure_qc_summary.tsv"), sep = "\t")

message("Frozen main figures written to: ", out_dir)
