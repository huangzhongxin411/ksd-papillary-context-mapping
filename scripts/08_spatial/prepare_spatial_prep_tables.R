suppressPackageStartupMessages(library(data.table))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/gene_sets", recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

spatial_datasets <- data.table(
  dataset = c("GSE206306", "GSE231630"),
  geo_accession = c("GSE206306", "GSE231630"),
  sample_id = c("pending", "pending"),
  data_type = "spatial_transcriptomics",
  platform = "Visium_or_spatial_platform_pending",
  raw_matrix_available = FALSE,
  histology_image_available = FALSE,
  spatial_coordinates_available = FALSE,
  metadata_available = FALSE,
  local_path = c("data/raw/GSE206306", "data/raw/GSE231630"),
  status = "missing",
  notes = "No local spatial matrix/image/coordinate files detected; preprocessing cannot proceed beyond manifest."
)

for (i in seq_len(nrow(spatial_datasets))) {
  path <- spatial_datasets$local_path[[i]]
  files <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else character()
  spatial_datasets[i, raw_matrix_available := any(grepl("matrix\\.mtx|filtered_feature_bc_matrix|\\.h5$|\\.h5ad$|\\.rds$", files, ignore.case = TRUE))]
  spatial_datasets[i, histology_image_available := any(grepl("tissue.*image|\\.png$|\\.jpg$|\\.jpeg$|\\.tif$|\\.tiff$", files, ignore.case = TRUE))]
  spatial_datasets[i, spatial_coordinates_available := any(grepl("tissue_positions|spatial.*coord|coordinates", files, ignore.case = TRUE))]
  spatial_datasets[i, metadata_available := any(grepl("metadata|family\\.soft|series", files, ignore.case = TRUE))]
  spatial_datasets[i, status := ifelse(length(files) > 0, "partial_local_files_detected", "missing")]
}

fwrite(spatial_datasets, "results/tables/spatial_resource_manifest.tsv", sep = "\t")

sample_metadata <- data.table(
  sample_id = c("pending_GSE206306", "pending_GSE231630"),
  dataset = c("GSE206306", "GSE231630"),
  donor_id = "unknown",
  disease_status = "unknown",
  stone_type = "unknown",
  section_id = "unknown",
  region = "renal_papilla_or_spatial_region_pending",
  platform = "unknown",
  has_image = FALSE,
  has_coordinates = FALSE,
  notes = "Template row; replace after spatial raw/processed files and sample metadata are available."
)
fwrite(sample_metadata, "config/spatial_sample_metadata.tsv", sep = "\t")

read_gene_set <- function(path) {
  if (!file.exists(path)) return(character())
  unique(trimws(readLines(path, warn = FALSE)))
}

magma_top50 <- read_gene_set("results/gene_sets/magma_top50.txt")
magma_top100 <- read_gene_set("results/gene_sets/magma_top100.txt")
magma_fdr <- read_gene_set("results/gene_sets/magma_fdr05.txt")
magma_suggestive <- read_gene_set("results/gene_sets/magma_suggestive_p1e4.txt")

tiers <- fread("results/tables/candidate_gene_tiers_v0.2.tsv")
pre_twas_tier2 <- tiers[current_tier == "Tier2_moderate_MAGMA_scRNA", gene]

marker_sets <- list(
  TAL_marker_set = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN14", "CASR"),
  vascular_stromal_marker_set = c("PECAM1", "VWF", "KDR", "COL1A1", "COL1A2", "PDGFRB", "RGS5", "ACTA2"),
  injury_remodeling_marker_set = c("HAVCR1", "LCN2", "SPP1", "MMP7", "MMP9", "GPNMB", "VCAM1")
)

gene_sets <- list(
  MAGMA_top50 = magma_top50,
  MAGMA_top100 = magma_top100,
  MAGMA_FDR = magma_fdr,
  MAGMA_suggestive = magma_suggestive,
  pre_TWAS_Tier2_candidate = pre_twas_tier2,
  TAL_marker_set = marker_sets$TAL_marker_set,
  vascular_stromal_marker_set = marker_sets$vascular_stromal_marker_set,
  injury_remodeling_marker_set = marker_sets$injury_remodeling_marker_set
)

for (nm in names(gene_sets)) {
  fwrite(data.table(gene = gene_sets[[nm]]), file.path("results/gene_sets", paste0("spatial_", nm, ".txt")), col.names = FALSE)
}

availability <- rbindlist(lapply(names(gene_sets), function(nm) {
  genes <- unique(gene_sets[[nm]])
  data.table(
    gene_set = nm,
    n_input_genes = length(genes),
    n_detected_in_spatial = NA_integer_,
    detected_fraction = NA_real_,
    missing_genes = NA_character_,
    notes = "Spatial expression matrix not available locally; gene detection cannot yet be evaluated."
  )
}))
fwrite(availability, "results/tables/spatial_gene_set_availability.tsv", sep = "\t")

qc <- spatial_datasets[, .(
  dataset,
  sample_id,
  status,
  spot_count = NA_integer_,
  median_counts = NA_real_,
  median_features = NA_real_,
  percent_mt = NA_real_,
  tissue_spot_fraction = NA_real_,
  gene_detection_rate = NA_real_,
  notes = "Spatial data unavailable locally; QC not run."
)]
fwrite(qc, "results/tables/spatial_qc_summary.tsv", sep = "\t")

pdf("results/figures/spatial_sample_overview.pdf", width = 8, height = 5)
plot.new()
text(0.5, 0.65, "Spatial sample overview", cex = 1.4)
text(0.5, 0.5, "GSE206306 / GSE231630 spatial files are not available locally.", cex = 1)
text(0.5, 0.38, "Only resource manifest and metadata templates were generated.", cex = 1)
dev.off()

pdf("results/figures/spatial_qc_violin.pdf", width = 8, height = 5)
plot.new()
text(0.5, 0.62, "Spatial QC not run", cex = 1.4)
text(0.5, 0.48, "No local spatial matrix, coordinates, or histology image files were detected.", cex = 1)
dev.off()

message("wrote\tresults/tables/spatial_resource_manifest.tsv")
message("wrote\tconfig/spatial_sample_metadata.tsv")
message("wrote\tresults/tables/spatial_gene_set_availability.tsv")
message("wrote\tresults/tables/spatial_qc_summary.tsv")
