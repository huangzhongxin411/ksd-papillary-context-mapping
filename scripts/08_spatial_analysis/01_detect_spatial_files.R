suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/spatial", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

roots <- c("data/raw", "data/processed", "results")
paths <- list.files(roots, recursive = TRUE, full.names = TRUE)
sp <- paths[grepl("GSE206306|GSE231630|spatial|Visium|tissue_positions|scalefactors|filtered_feature_bc_matrix|matrix.mtx|barcodes.tsv|features.tsv|\\.tif$|\\.png$", paths, ignore.case = TRUE)]
datasets <- c("GSE206306", "GSE231630")
status <- rbindlist(lapply(datasets, function(ds) {
  h <- sp[grepl(ds, sp, ignore.case = TRUE)]
  data.table(
    dataset = ds,
    matrix = any(grepl("matrix.mtx|filtered_feature_bc_matrix|\\.h5$", h, ignore.case = TRUE)),
    barcodes = any(grepl("barcodes.tsv", h, ignore.case = TRUE)),
    features = any(grepl("features.tsv|genes.tsv", h, ignore.case = TRUE)),
    positions = any(grepl("tissue_positions|positions", h, ignore.case = TRUE)),
    scalefactors = any(grepl("scalefactors", h, ignore.case = TRUE)),
    image = any(grepl("\\.tif$|\\.png$|image", h, ignore.case = TRUE)),
    candidate_files = paste(head(h, 10), collapse = "; ")
  )
}), fill = TRUE)
status[, spatial_ready := matrix & barcodes & features & positions & scalefactors]
fwrite(status, "results/spatial/spatial_resource_status.tsv", sep = "\t")
writeLines(c(
  "# Spatial Resource Blocker Memo v0.1",
  "",
  "Spatial transcriptomics will only be interpreted as spatial context support if matrix/barcodes/features/positions/scalefactors are available.",
  "If image files are missing, overlays may be limited to coordinate plots and cannot be described as H&E-anchored spatial validation."
), "docs/spatial_resource_blocker_memo.md", useBytes = TRUE)
message("wrote spatial resource status")
