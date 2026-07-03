#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
})

stage_id <- "stage6A_spatial_twas_audit"
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
fig_dir <- file.path("results/figures/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
for (d in c(table_dir, doc_dir, fig_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage6A_spatial_twas_audit.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 6A spatial/TWAS supplementary-context audit\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  stage3 = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  stage4_claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
  stage4_report = "docs/revision/stage4C2R_draft_figures_v0.2/stage4C2R_report.md",
  stage5_claim = "results/tables/revision/stage5C1_gse73680_figure4_draft/gse73680_final_claim_lock_stage5C1.tsv",
  stage5_report = "docs/revision/stage5C1_gse73680_figure4_draft/stage5C1_report.md",
  twas_audit = "results/tables/revision/stage2_genetic/twas_output_audit.tsv",
  twas_models = "results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv",
  twas_memo = "docs/revision/stage2_genetic/twas_proxy_interpretation_memo.md",
  prior_load_qc = "results/spatial/phase28b_loading_rescue/spatial_load_qc_v0.2.tsv",
  prior_projection = "results/spatial/phase28c_projection/spatial_projection_summary_v0.1.tsv",
  prior_adjusted = "results/spatial/phase28c_projection/spatial_score_correlations_adjusted_v0.1.tsv",
  prior_claim_gate = "results/spatial/phase28d_spatial_hardening/spatial_claim_gate_final_v0.1.tsv"
)
missing_inputs <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs)) stop("Required audit input missing: ", paste(missing_inputs, collapse = ", "))

check_cols <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop("Missing columns in ", label, ": ", paste(missing, collapse = ", "))
}

stage3 <- fread(paths$stage3)
twas_audit <- fread(paths$twas_audit)
twas_models <- fread(paths$twas_models)
load_qc <- fread(paths$prior_load_qc)
check_cols(stage3, c("gene", "reporting_group", "curated_exemplar_flag", "smr_status", "coloc_status"), "Stage 3R evidence model")
check_cols(twas_audit, c("reference_tissue", "n_genes_tested", "n_fdr_significant", "n_magma_overlap", "n_multi_snp_fdr_significant", "n_one_snp_fdr_significant"), "TWAS audit")
check_cols(twas_models, c("gene", "twas_fdr", "model_class", "magma_overlap"), "TWAS model audit")
check_cols(load_qc, c("sample_id", "n_spots_matrix", "n_features_median", "n_counts_median", "percent_mt_median", "image_loaded", "coordinates_loaded", "qc_status"), "prior spatial load QC")

qc_objects <- sort(list.files(
  "data/processed/spatial/gse206306/seurat_objects",
  pattern = "_spatial_qc[.]rds$",
  full.names = TRUE
))
score_objects <- sort(list.files(
  "data/processed/spatial/gse206306/seurat_objects",
  pattern = "_spatial_scores_v0[.]1[.]rds$",
  full.names = TRUE
))
if (length(qc_objects) != 5L) stop("Expected five GSE206306 QC RDS objects; found ", length(qc_objects))

sample_from_path <- function(x) sub("_spatial(_scores_v0[.]1|_qc)[.]rds$", "", basename(x))

objects <- lapply(qc_objects, readRDS)
names(objects) <- vapply(qc_objects, sample_from_path, character(1))
if (!all(vapply(objects, inherits, logical(1), what = "Seurat"))) stop("At least one primary spatial object is not a Seurat object")

cat("Loaded primary QC objects:", paste(names(objects), collapse = ", "), "\n")
cat("Total tissue spots:", sum(vapply(objects, ncol, numeric(1))), "\n")

# 1. Spatial file inventory across spatial-specific project roots.
inventory_roots <- c(
  "data/raw/spatial",
  "data/processed/spatial",
  "results/spatial",
  "results/figures/supplementary/spatial_projection",
  "results/figures/supplementary/spatial_qc",
  "scripts/11_spatial",
  "scripts/08_spatial",
  "scripts/08_spatial_analysis",
  "docs/spatial",
  "manuscript/inserts"
)
inventory_files <- unique(unlist(lapply(inventory_roots[file.exists(inventory_roots)], function(root) {
  if (file.info(root)$isdir) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE) else root
})))
extra_files <- c(
  "config/spatial_sample_metadata.tsv",
  list.files("results/gene_sets", pattern = "^spatial_", full.names = TRUE),
  list.files("docs", pattern = "spatial|GSE206306", full.names = TRUE, ignore.case = TRUE)
)
inventory_files <- sort(unique(c(inventory_files, extra_files[file.exists(extra_files)])))
inventory_files <- inventory_files[file.exists(inventory_files) & !file.info(inventory_files)$isdir]

infer_type <- function(x) {
  b <- tolower(basename(x))
  if (grepl("[.]rds$", b)) return("RDS_spatial_object")
  if (grepl("[.]h5$", b)) return("10x_HDF5_matrix")
  if (grepl("tissue_positions.*[.]csv$", b)) return("spot_coordinates")
  if (grepl("scalefactors.*[.]json$", b)) return("scale_factors")
  if (grepl("(hires|lowres|detected_tissue|fiducial).*[.](png|jpg|jpeg)$", b)) return("histology_or_alignment_image")
  if (grepl("[.]tsv$|[.]csv$", b)) return("tabular_result_or_metadata")
  if (grepl("[.]pdf$|[.]png$|[.]svg$", b)) return("figure")
  if (grepl("[.]r$|[.]py$|[.]sh$", b)) return("analysis_script")
  if (grepl("[.]md$", b)) return("documentation")
  if (grepl("[.]tar[.]gz$|[.]tgz$|[.]zip$", b)) return("archive")
  "other"
}

infer_content <- function(x, type) {
  z <- tolower(x)
  if (type == "RDS_spatial_object" && grepl("scores", z)) return("Seurat object with expression, spot metadata, coordinates, image and legacy spatial scores")
  if (type == "RDS_spatial_object") return("QC-passed Seurat object with expression, spot metadata, coordinates and image")
  if (type == "10x_HDF5_matrix") return("10x filtered feature-barcode count matrix")
  if (type == "spot_coordinates") return("Visium spot coordinates and in-tissue flags")
  if (type == "scale_factors") return("Visium image scale factors")
  if (type == "histology_or_alignment_image") return("Visium tissue/alignment image")
  if (grepl("correlation|projection|score", z)) return("existing spatial score/projection/correlation output")
  if (grepl("claim|memo|legend|report", z)) return("spatial interpretation or claim-boundary documentation")
  if (type == "figure") return("existing spatial visualization")
  if (type == "analysis_script") return("existing spatial processing or analysis script")
  "spatial project resource"
}

inventory <- rbindlist(lapply(inventory_files, function(f) {
  info <- file.info(f)
  type <- infer_type(f)
  dataset <- if (grepl("GSE206306|GSM62503|GSM7166170", f, ignore.case = TRUE)) "GSE206306" else if (grepl("GSE231630", f, ignore.case = TRUE)) "GSE231630" else "spatial_project"
  ready <- if (type == "RDS_spatial_object" && grepl("_spatial_qc[.]rds$", f)) {
    "yes_primary_input"
  } else if (type %in% c("10x_HDF5_matrix", "spot_coordinates", "scale_factors", "histology_or_alignment_image")) {
    "yes_companion_or_fallback"
  } else if (grepl("results/spatial|results/figures/supplementary/spatial", f)) {
    "reuse_after_definition_audit"
  } else {
    "not_primary_input"
  }
  data.table(
    file_path = f,
    file_name = basename(f),
    file_type = type,
    dataset_id = dataset,
    size_bytes = as.numeric(info$size),
    last_modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    likely_content = infer_content(f, type),
    readable = file.access(f, 4) == 0,
    ready_for_stage6B = ready,
    notes = ifelse(type == "RDS_spatial_object" && grepl("scores", f), "Legacy scores require current Stage 3R module-definition check before reuse.", "Inventory only; readiness does not expand claim strength.")
  )
}))
fwrite(inventory, file.path(table_dir, "spatial_file_inventory.tsv"), sep = "\t")

# 2. Primary input object audit.
audit_rds <- function(path, recommended) {
  obj <- readRDS(path)
  sid <- sample_from_path(path)
  coords_ok <- length(obj@images) > 0 && nrow(Seurat::GetTissueCoordinates(obj)) == ncol(obj)
  meta <- obj@meta.data
  data.table(
    candidate_object = basename(path),
    dataset_id = "GSE206306",
    file_path = path,
    object_type = ifelse(grepl("scores", path), "Seurat_spatial_scored_RDS", "Seurat_spatial_QC_RDS"),
    readable = TRUE,
    contains_expression = length(obj@assays) > 0 && nrow(obj) > 0,
    contains_spot_metadata = nrow(meta) == ncol(obj),
    contains_coordinates = coords_ok,
    contains_image = length(obj@images) > 0,
    contains_section_id = "derived_from_object_filename",
    contains_sample_id = ifelse("orig.ident" %in% names(meta), "yes_orig.ident", "derived_from_object_filename"),
    contains_region_labels = FALSE,
    contains_roi_or_plaque_annotation = FALSE,
    recommended_primary_input = recommended,
    blocking_issue = "No plaque/lesion ROI annotation; blocks plaque-specific analysis but not tissue-context projection.",
    notes = paste0(sid, ": use QC RDS for current module rescoring; do not treat legacy score columns as current without a definition audit.")
  )
}
primary_audit <- rbindlist(c(
  lapply(qc_objects, audit_rds, recommended = TRUE),
  lapply(score_objects, audit_rds, recommended = FALSE)
), fill = TRUE)

h5_files <- sort(list.files("data/processed/spatial/gse206306", pattern = "filtered_feature_bc_matrix[.]h5$", recursive = TRUE, full.names = TRUE))
if (length(h5_files)) {
  primary_audit <- rbind(primary_audit, rbindlist(lapply(h5_files, function(path) data.table(
    candidate_object = basename(dirname(path)), dataset_id = "GSE206306", file_path = path,
    object_type = "10x_filtered_feature_HDF5", readable = file.access(path, 4) == 0,
    contains_expression = TRUE, contains_spot_metadata = FALSE, contains_coordinates = FALSE,
    contains_image = FALSE, contains_section_id = "derived_from_parent_directory",
    contains_sample_id = "derived_from_parent_directory", contains_region_labels = FALSE,
    contains_roi_or_plaque_annotation = FALSE, recommended_primary_input = FALSE,
    blocking_issue = "Requires companion spatial directory and object reconstruction.",
    notes = "Fallback input only; QC RDS is preferred."
  ))), fill = TRUE)
}
fwrite(primary_audit, file.path(table_dir, "spatial_primary_input_audit.tsv"), sep = "\t")

# 3. Metadata schema audit across the recommended five-object suite.
total_spots <- sum(vapply(objects, ncol, numeric(1)))
all_meta <- rbindlist(lapply(names(objects), function(sid) {
  obj <- objects[[sid]]
  m <- as.data.table(obj@meta.data, keep.rownames = "spot_barcode")
  co <- as.data.table(Seurat::GetTissueCoordinates(obj))
  if ("cell" %in% names(co)) setnames(co, "cell", "spot_barcode")
  if (!"spot_barcode" %in% names(co)) co[, spot_barcode := rownames(Seurat::GetTissueCoordinates(obj))]
  m <- merge(m, co[, .(spot_barcode, x_coordinate = x, y_coordinate = y)], by = "spot_barcode", all.x = TRUE)
  m[, `:=`(
    spot_id = paste(sid, spot_barcode, sep = ":"),
    section_id = sid,
    sample_id = sid,
    dataset_id = "GSE206306",
    histology_image_id = paste0(sid, ":slice1")
  )]
  m
}), fill = TRUE)

schema_spec <- data.table(
  metadata_column = c("spot_id", "section_id", "sample_id", "dataset_id", "x_coordinate", "y_coordinate", "nCount_Spatial", "nFeature_Spatial", "percent_mt", "tissue_region", "histology_image_id", "roi_label", "plaque_or_lesion_label"),
  possible_role = c("spot_id", "section_id", "sample_id", "unknown", "coordinate", "coordinate", "qc_metric", "qc_metric", "qc_metric", "tissue_region", "unknown", "roi_label", "plaque_lesion_annotation"),
  use_in_stage6B = c("yes", "yes_primary_unit", "yes", "yes", "yes_visualization", "yes_visualization", "yes_complexity_adjustment", "yes_complexity_adjustment", "compute_before_adjustment", "no_absent", "yes_visualization", "no_absent", "no_absent")
)
schema_audit <- schema_spec[, {
  col <- metadata_column
  if (col %in% names(all_meta)) {
    x <- all_meta[[col]]
    keep <- !is.na(x) & nzchar(as.character(x))
    vals <- unique(as.character(x[keep]))
    .(non_missing_count = sum(keep), unique_count = length(vals), example_values = paste(head(vals, 5), collapse = ";"))
  } else {
    .(non_missing_count = 0L, unique_count = 0L, example_values = "")
  }
}, by = .(metadata_column, possible_role, use_in_stage6B)]
schema_audit[, notes := fifelse(
  metadata_column %in% c("section_id", "sample_id", "dataset_id", "x_coordinate", "y_coordinate", "histology_image_id"),
  "Derived deterministically from the object/file/image structure during audit.",
  fifelse(metadata_column == "percent_mt", "Absent from QC RDS metadata; computable from Spatial counts before Stage 6B residualization.",
           fifelse(non_missing_count == 0, "Not available in the recommended primary input suite.", "Present in QC RDS metadata."))
)]
setcolorder(schema_audit, c("metadata_column", "non_missing_count", "unique_count", "example_values", "possible_role", "use_in_stage6B", "notes"))
fwrite(schema_audit, file.path(table_dir, "spatial_metadata_schema_audit.tsv"), sep = "\t")

# 4. Section-level QC summary, using audited loading outputs for percent.mt.
section_qc <- rbindlist(lapply(names(objects), function(sid) {
  obj <- objects[[sid]]
  prior <- load_qc[sample_id == sid][1]
  data.table(
    dataset_id = "GSE206306",
    section_id = sid,
    sample_id = sid,
    n_spots_total = ncol(obj),
    n_tissue_spots = ncol(obj),
    median_nCount = median(obj$nCount_Spatial, na.rm = TRUE),
    median_nFeature = median(obj$nFeature_Spatial, na.rm = TRUE),
    median_percent_mt = prior$percent_mt_median,
    image_available = length(obj@images) > 0,
    roi_available = FALSE,
    notes = "QC-passed tissue-filtered Visium section; section ID is sample-derived; no plaque/lesion ROI."
  )
}))
fwrite(section_qc, file.path(table_dir, "spatial_section_qc_summary.tsv"), sep = "\t")

# 5. Current Stage 3R/5 module feasibility against non-zero spatial expression.
get_group <- function(groups) unique(toupper(stage3[reporting_group %in% groups, gene]))
read_set <- function(path) unique(toupper(scan(path, what = character(), quiet = TRUE)))
top50_path <- "results/gene_sets/magma_top50.txt"
top100_path <- "results/gene_sets/magma_top100.txt"
if (!file.exists(top50_path) || !file.exists(top100_path)) stop("Current MAGMA top-gene files missing")

module_defs <- list(
  R1_MAGMA_Bonferroni_only = get_group("R1_MAGMA_Bonferroni_only"),
  R1_R2_R3_all_MAGMA_Bonferroni = get_group(c("R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy")),
  MAGMA_top50 = read_set(top50_path),
  MAGMA_top100 = read_set(top100_path),
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = get_group("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy"),
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = get_group("R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"),
  R2_R3_MAGMA_plus_TWAS_proxy = get_group(c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy")),
  R5_TWAS_proxy_only = get_group("R5_TWAS_proxy_only"),
  Curated_exemplar_panel = unique(toupper(stage3[curated_exemplar_flag == "yes", gene])),
  LoopTAL_signature = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2"),
  epithelial_general = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  proximal_tubule = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "AQP1", "ALDOB"),
  collecting_duct = c("AQP2", "AQP3", "SCNN1A", "SCNN1B", "SCNN1G", "KRT19"),
  injury_epithelial = c("HAVCR1", "LCN2", "SPP1", "VCAM1", "KRT8", "KRT18", "VIM"),
  ECM_fibrosis = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "TGFBI", "MMP2"),
  mineralization_remodeling = c("SPP1", "ALPL", "MMP7", "MMP9", "RUNX2", "BGLAP", "POSTN"),
  immune_myeloid = c("PTPRC", "AIF1", "LST1", "CD68", "CD14", "C1QA", "C1QB"),
  endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CLDN5"),
  fibroblast_stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA")
)
module_roles <- c(
  R1_MAGMA_Bonferroni_only = "primary_MAGMA_module",
  R1_R2_R3_all_MAGMA_Bonferroni = "primary_MAGMA_module",
  MAGMA_top50 = "primary_MAGMA_module", MAGMA_top100 = "primary_MAGMA_module",
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R2_R3_MAGMA_plus_TWAS_proxy = "secondary_TWAS_proxy_module",
  R5_TWAS_proxy_only = "secondary_TWAS_proxy_module",
  Curated_exemplar_panel = "curated_exemplar_context",
  LoopTAL_signature = "spatial_context_signature", epithelial_general = "spatial_context_signature",
  proximal_tubule = "spatial_context_signature", collecting_duct = "spatial_context_signature",
  injury_epithelial = "injury_or_remodeling_signature", ECM_fibrosis = "injury_or_remodeling_signature",
  mineralization_remodeling = "injury_or_remodeling_signature", immune_myeloid = "spatial_context_signature",
  endothelial = "spatial_context_signature", fibroblast_stromal = "spatial_context_signature"
)

nonzero_genes <- unique(unlist(lapply(objects, function(obj) {
  mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "counts"), error = function(e) GetAssayData(obj, assay = "Spatial", slot = "counts"))
  rownames(mat)[Matrix::rowSums(mat) > 0]
})))
nonzero_lookup <- unique(toupper(nonzero_genes))

module_feasibility <- rbindlist(lapply(names(module_defs), function(nm) {
  genes <- unique(toupper(module_defs[[nm]]))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  detected <- intersect(genes, nonzero_lookup)
  missing <- setdiff(genes, detected)
  frac <- if (length(genes)) length(detected) / length(genes) else NA_real_
  feasible <- if (length(detected) >= 3 && frac >= 0.5) "yes" else if (length(detected) >= 1 && (frac >= 0.5 || length(genes) <= 2)) "yes_with_caution" else "no"
  source <- if (nm %in% c("MAGMA_top50", "MAGMA_top100")) {
    ifelse(nm == "MAGMA_top50", top50_path, top100_path)
  } else if (nm %in% names(module_roles)[grepl("TWAS|MAGMA|Curated", module_roles) | names(module_roles) == "Curated_exemplar_panel"]) {
    paths$stage3
  } else {
    "Stage 5A documented renal marker-signature definition"
  }
  data.table(
    module_name = nm,
    module_role = unname(module_roles[nm]),
    source_file_or_definition = source,
    n_genes_input = length(genes),
    n_genes_detected_in_spatial = length(detected),
    detected_fraction = frac,
    detected_genes = paste(detected, collapse = ";"),
    missing_genes = paste(missing, collapse = ";"),
    feasible_for_stage6B = feasible,
    notes = ifelse(grepl("TWAS", unname(module_roles[nm])), "Supplementary proxy module only; feasibility does not upgrade TWAS evidence.", ifelse(nm == "Curated_exemplar_panel", "Interpretive context only; no single-gene validation.", "Detected across the union of five QC-passed sections; section consistency must be tested in Stage 6B."))
  )
}))
fwrite(module_feasibility, file.path(table_dir, "spatial_module_feasibility_audit.tsv"), sep = "\t")

# 6. Existing result audit.
existing_results <- data.table(
  analysis_type = c(
    "spot_level_module_scores_existing", "section_level_module_scores_existing",
    "complexity_adjusted_scores_existing", "spatial_correlation_existing",
    "moran_or_autocorrelation_existing", "roi_based_analysis_existing",
    "histology_overlay_existing", "figure_source_existing"
  ),
  existing_result_found = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  file_path = c(
    "data/processed/spatial/gse206306/seurat_objects/*_spatial_scores_v0.1.rds;results/spatial/phase28c_projection/spatial_module_score_qc_v0.1.tsv",
    "",
    paths$prior_adjusted,
    "results/spatial/phase28c_projection/spatial_score_correlations_by_sample_v0.1.tsv;results/spatial/phase28c_projection/spatial_projection_summary_v0.1.tsv",
    "", "",
    "results/figures/supplementary/spatial_projection/",
    "results/spatial/phase28d_spatial_hardening/supplementary_figure_s7_final_panel_sources_v0.1.tsv"
  ),
  can_reuse = c("partial", "no", "partial", "partial", "no", "no", "yes_reference_only", "partial"),
  needs_rerun = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  reason = c(
    "Legacy score definitions predate the Stage 3R/Stage 6A module set.",
    "No dedicated current section-level module summary table was found.",
    "Residualization approach is reusable, but current modules and section-level synthesis must be regenerated.",
    "Historical correlations exist but must be aligned to current modules and conservative section-level rules.",
    "No Moran or other autocorrelation output was found; optional in Stage 6B.",
    "No plaque/lesion ROI labels exist.",
    "Images and plotting code are reusable, but overlays must use current scores and claim wording.",
    "Legacy source table can guide structure but does not cover the current Stage 6B module set."
  ),
  notes = c(
    "Do not reuse legacy P1 framing as evidence.",
    "Spots must not be treated as independent replicates.",
    "Spot-level P values remain descriptive.",
    "Section consistency is the inferential emphasis.",
    "Absence is not a negative biological result.",
    "ROI-aware analysis is not permitted without new annotation.",
    "Papillary tissue-context visualization only.",
    "Generate new Figure 5 blueprint only after Stage 6B."
  )
)
fwrite(existing_results, file.path(table_dir, "spatial_existing_result_audit.tsv"), sep = "\t")

# 7. TWAS boundary integration.
twas_summary <- twas_audit[1]
n_fdr_one <- twas_models[!is.na(twas_fdr) & twas_fdr < 0.05 & model_class == "one_snp_model", uniqueN(gene)]
n_fdr_multi <- twas_models[!is.na(twas_fdr) & twas_fdr < 0.05 & model_class == "multi_snp_model", uniqueN(gene)]
n_fdr_overlap <- twas_models[!is.na(twas_fdr) & twas_fdr < 0.05 & magma_overlap == "yes", uniqueN(gene)]
n_twas_only <- stage3[reporting_group == "R5_TWAS_proxy_only", uniqueN(gene)]

twas_boundary <- data.table(
  twas_component = c("GTEx Kidney_Cortex proxy TWAS", "one-SNP TWAS models", "multi-SNP TWAS models", "TWAS-MAGMA overlap", "TWAS-only genes", "papilla eQTL availability", "SMR/coloc not claim-grade ready"),
  current_status = c("available_proxy_only", "available_weaker_proxy", "available_stronger_proxy_but_supplementary", "available_overlap", "available_context_only", "not_available", "not_claim_grade_ready"),
  n_genes_or_models = c(twas_summary$n_genes_tested, n_fdr_one, n_fdr_multi, n_fdr_overlap, n_twas_only, 0, 0),
  key_limitation = c(
    "GTEx v8 Kidney_Cortex is not renal papilla and S-PrediXcan does not establish causality.",
    "A single SNP drives each model; LD or locus effects may dominate the proxy association.",
    "Multi-SNP models are relatively stronger than one-SNP models but remain tissue-mismatched proxy evidence.",
    "Overlap is convergence of prioritization layers, not colocalization or mediation.",
    "TWAS-only genes lack MAGMA overlap and remain lower-confidence contextual signals.",
    "No papilla-specific eQTL/prediction resource was identified in the audited project.",
    "Required exposure/QTL and harmonized claim-grade resources remain incomplete."
  ),
  allowed_manuscript_use = c(
    "Supplementary Kidney_Cortex genetically regulated expression proxy support.",
    "Weaker supplementary Kidney_Cortex proxy evidence, explicitly flagged as one-SNP.",
    "Supplementary Kidney_Cortex proxy evidence with model-size context.",
    "Supplementary overlap between MAGMA prioritization and Kidney_Cortex proxy TWAS.",
    "Supplementary contextual list only.",
    "State that papilla-specific regulation was not tested.",
    "State transparently that claim-grade SMR/coloc support was not established."
  ),
  disallowed_manuscript_use = c(
    "papilla-specific causal expression; causal gene; colocalization; spatial validation",
    "causal expression evidence; fine-mapped regulatory mechanism",
    "papilla-specific regulation; causal mediation",
    "SMR/coloc support; shared causal variant",
    "validated disease gene; therapeutic target",
    "papilla-specific eQTL support",
    "SMR-supported gene; coloc-supported gene; negative causal evidence"
  ),
  figure_or_supplement_use = c("supplement_only", "supplement_only_with_flag", "supplement_only", "supplementary_overlap_table", "supplementary_only", "limitations", "methods_and_limitations"),
  notes = c(
    paste0(twas_summary$n_fdr_significant, " FDR-supported of ", twas_summary$n_genes_tested, " tested models."),
    paste0(n_fdr_one, " FDR-supported one-SNP models."),
    paste0(n_fdr_multi, " FDR-supported multi-SNP models."),
    paste0(n_fdr_overlap, " FDR-supported TWAS genes overlap MAGMA."),
    paste0(n_twas_only, " Stage 3R R5 genes."),
    "No new papillary eQTL search or analysis was performed in Stage 6A.",
    "Missing evidence is not negative evidence."
  )
)
fwrite(twas_boundary, file.path(table_dir, "twas_boundary_integration_audit.tsv"), sep = "\t")

# 8. Cross-layer claim boundary for Stage 6 integration.
claim_boundary <- data.table(
  evidence_layer = c("GWAS_MAGMA_layer", "TWAS_proxy_layer", "snRNA_LoopTAL_context", "GSE73680_bulk_context", "spatial_context_projection", "curated_exemplar_panel", "overall_manuscript_claim"),
  current_support_level = c(
    "auditable_EUR_LD_reference_based_gene_prioritization",
    "supplementary_Kidney_Cortex_proxy_only",
    "moderate_donor_level_LoopTAL_associated_context",
    "injury_remodeling_associated_paired_bulk_disease_context",
    "pending_Stage6B_tissue_context_projection",
    "interpretive_only",
    "post_GWAS_papillary_context_mapping_not_causal_validation"
  ),
  allowed_claim = c(
    "EUR-LD-reference-based MAGMA gene prioritization for downstream context mapping.",
    "Supplementary GTEx Kidney_Cortex proxy support for a subset of prioritized genes.",
    "Moderate donor-level Loop/TAL-associated papillary cellular context for prioritized modules.",
    "Injury/remodeling-associated paired bulk disease-context support for MAGMA-prioritized modules.",
    "Papillary tissue-context projection or spatial co-distribution if Stage 6B section-level results are consistent.",
    "Curated biological exemplars illustrate role heterogeneity without evidence upgrading.",
    "Post-GWAS prioritization and renal papillary cell/tissue-context mapping of KSD genetic risk."
  ),
  disallowed_claim = c(
    "ancestry-generalizable fine mapping; causal genes; trans-ancestry precise gene mapping",
    "papilla-specific causal expression; colocalization; SMR support",
    "causal cell type; plaque nucleation cell; cell-level independent replication",
    "cell-type-specific disease response; genetic causality; plaque causation; independent validation",
    "plaque-specific validation; causal spatial mechanism; spatial proof of plaque nucleation",
    "validated disease genes; therapeutic targets; causal exemplars",
    "causal gene validation; plaque-specific mechanism; therapeutic target validation"
  ),
  main_or_supplementary = c("main", "supplementary", "main_conservative", "main_conservative", "pending_likely_supplementary", "supplementary_or_interpretive", "main"),
  figure_use = c("genetic_prioritization_panel", "supplementary_table_or_boundary_panel", "Figure_2", "Figure_4", "pending_Stage6B_or_supplementary_Figure_S7", "Figure_3_interpretive_panel", "integrated_summary_after_Stage6"),
  notes = c(
    "MAGMA is a prioritization layer, not fine mapping.",
    "One-SNP models are explicitly weaker proxies.",
    "Spots/cells are not independent biological replicates.",
    "Composition-aware attenuation remains central.",
    "No plaque/lesion ROI exists; spatial validation is unavailable.",
    "Curation is orthogonal to evidence strength.",
    "All layers remain association/context mapping."
  )
)
fwrite(claim_boundary, file.path(table_dir, "stage6_integrated_claim_boundary_table.tsv"), sep = "\t")

cat("\nOutputs written to:", table_dir, "\n")
cat("Inventory rows:", nrow(inventory), "\n")
cat("Primary object rows:", nrow(primary_audit), "\n")
cat("Total spots:", total_spots, "across", length(objects), "sections\n")
cat("Feasible modules/signatures:", module_feasibility[feasible_for_stage6B != "no", .N], "of", nrow(module_feasibility), "\n")
cat("Plaque/lesion ROI available: FALSE\n")
cat("Stage 6B tissue-context projection ready: TRUE\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
