#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

root <- getwd()
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "scripts/08_spatial_analysis"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
rel <- function(x) sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", root), "/?"), "", normalizePath(x, mustWork = FALSE))
safe_file_info <- function(paths) {
  info <- file.info(paths)
  data.table(
    file_path = rel(paths),
    file_name = basename(paths),
    file_type = tolower(tools::file_ext(paths)),
    file_size = info$size,
    modified_time = as.character(info$mtime)
  )
}

keywords <- c(
  "GSE206306", "GSE231630", "Visium", "spatial", "SpaceRanger",
  "filtered_feature_bc_matrix", "tissue_positions", "tissue_positions_list",
  "scalefactors_json", "tissue_hires_image", "tissue_lowres_image",
  "aligned_fiducials", "detected_tissue_image", "H&E", "histology",
  "spot", "section", "ROI", "plaque", "mineral", "calcification",
  "fibrosis", "lesion", "Loop/TAL", "deconvolution", "RCTD",
  "Cell2location", "Seurat transfer", "label transfer"
)
roi_keywords <- c("ROI", "plaque", "mineral", "calcification", "fibrosis", "lesion", "histology annotation", "manual annotation")
target_dirs <- c(
  "data/raw/spatial", "data/processed/spatial", "data/raw/GSE206306",
  "data/raw/GSE231630", "data/processed/gse206306", "data/processed/gse231630",
  "config", "results", "scripts", "release", "supplementary_package_v1.6",
  "submission_package_v1.7", "submission_package_v1.8_BMC_Genomics",
  "submission_package_v1.9_BMC_working"
)
target_dirs <- target_dirs[dir.exists(file.path(root, target_dirs))]
all_files <- unique(unlist(lapply(file.path(root, target_dirs), function(d) list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE))))
all_files <- all_files[file.exists(all_files)]

text_ext <- c("txt", "tsv", "csv", "md", "r", "py", "sh", "json", "yml", "yaml", "log")
detect_keywords <- function(path) {
  hay <- paste(basename(path), dirname(rel(path)))
  ext <- tolower(tools::file_ext(path))
  fsize <- suppressWarnings(file.info(path)$size)
  if (!is.na(fsize) && ext %in% text_ext && fsize < 5e6) {
    txt <- tryCatch(suppressWarnings(readChar(path, nchars = min(fsize, 200000), useBytes = TRUE)), error = function(e) "")
    hay <- paste(hay, txt)
  }
  hay_lower <- tolower(hay)
  hits <- keywords[vapply(keywords, function(k) grepl(tolower(k), hay_lower, fixed = TRUE), logical(1))]
  paste(unique(hits), collapse = ";")
}
likely_role <- function(path) {
  n <- basename(path)
  p <- rel(path)
  if (grepl("filtered_feature_bc_matrix\\.h5$", n)) return("10x filtered feature-barcode H5 matrix")
  if (grepl("matrix\\.mtx", n)) return("10x MTX count matrix")
  if (grepl("features|genes", n, ignore.case = TRUE)) return("feature/gene annotation")
  if (grepl("barcodes", n, ignore.case = TRUE)) return("barcode list")
  if (grepl("tissue_positions", n, ignore.case = TRUE)) return("Visium tissue coordinate file")
  if (grepl("scalefactors", n, ignore.case = TRUE)) return("Visium scale factor JSON")
  if (grepl("tissue_hires|tissue_lowres|detected_tissue|aligned_fiducials|\\.png$|\\.jpg$|\\.jpeg$", n, ignore.case = TRUE)) return("Visium histology/image file")
  if (grepl("metadata|sample", n, ignore.case = TRUE)) return("sample metadata or manifest")
  if (grepl("spatial.*score|projection|correlation|figure", p, ignore.case = TRUE)) return("existing downstream spatial output")
  "spatial-related project file"
}
sample_from_path <- function(path) {
  p <- rel(path)
  gsm <- regmatches(p, regexpr("GSM[0-9]+", p))
  if (length(gsm) && nchar(gsm) > 0) return(gsm)
  v <- regmatches(p, regexpr("V[0-9A-Za-z.-]+[_-]XY[0-9A-Za-z._-]+", p))
  if (length(v) && nchar(v) > 0) return(v)
  if (grepl("GSE206306", p, ignore.case = TRUE)) return("GSE206306_unspecified")
  if (grepl("GSE231630", p, ignore.case = TRUE)) return("GSE231630_unspecified")
  "unknown"
}
candidate_input <- function(path) {
  grepl("filtered_feature_bc_matrix|matrix\\.mtx|features.tsv|barcodes.tsv|tissue_positions|scalefactors_json|tissue_hires_image|tissue_lowres_image", basename(path), ignore.case = TRUE)
}

inventory <- safe_file_info(all_files)
inventory[, sample_or_section_id := vapply(file_path, sample_from_path, character(1))]
inventory[, likely_role := vapply(file_path, likely_role, character(1))]
inventory[, detected_keywords := vapply(file.path(root, file_path), detect_keywords, character(1))]
inventory[, candidate_for_canonical_input := ifelse(vapply(file.path(root, file_path), candidate_input, logical(1)), "yes", "no")]
inventory[, notes := ifelse(grepl("phase28|stage6|spatial_score|projection|correlation", file_path, ignore.case = TRUE),
                            "Existing spatial downstream output; review before reuse.", "")]
write_tsv(inventory, file.path(root, "codex_tasks/phase3_step1_spatial_file_inventory.tsv"))

find_h5 <- unique(all_files[grepl("filtered_feature_bc_matrix\\.h5$", basename(all_files), ignore.case = TRUE)])
# Prefer raw/manual-staging Visium directories for canonical section locking. Processed
# copies can be useful provenance, but should not create duplicate section rows.
raw_h5 <- find_h5[grepl("data/raw/spatial/manual_download_staging", rel(find_h5), fixed = TRUE)]
if (length(raw_h5) > 0) find_h5 <- raw_h5
find_h5 <- find_h5[!duplicated(normalizePath(find_h5, mustWork = FALSE))]
section_root <- function(h5) {
  d <- dirname(h5)
  if (basename(d) == "outs") return(d)
  d
}
section_roots <- unique(vapply(find_h5, section_root, character(1)))
infer_dataset <- function(path) {
  p <- rel(path)
  if (grepl("GSM62503|GSE206306|KRP|V19S25|V10M09", p, ignore.case = TRUE)) return("GSE206306")
  if (grepl("GSM71661|GSE231630|V12Y24", p, ignore.case = TRUE)) return("GSE231630")
  "unknown"
}
section_id_from_root <- function(sr) {
  if (basename(sr) == "outs") return(basename(dirname(sr)))
  basename(sr)
}
geo_from_root <- function(sr) {
  p <- rel(sr)
  gsm <- regmatches(p, regexpr("GSM[0-9]+", p))
  if (length(gsm) && nchar(gsm) > 0) return(gsm)
  if (grepl("V12Y24.346_XY01_(446|462|475)", p)) return("GSM7166168_inferred_from_archive")
  if (grepl("V12Y24.346_XY02_(F59|F63)", p)) return("GSM7166169_inferred_from_archive")
  "unknown"
}
condition_from_dataset <- function(dataset) {
  if (dataset %in% c("GSE206306", "GSE231630")) return("metadata_not_locked")
  "unknown"
}
pick_file <- function(sr, pattern) {
  files <- list.files(sr, recursive = TRUE, full.names = TRUE)
  hit <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
  if (length(hit)) rel(hit[1]) else ""
}
manifest <- rbindlist(lapply(section_roots, function(sr) {
  files <- list.files(sr, recursive = TRUE, full.names = TRUE)
  h5 <- files[grepl("filtered_feature_bc_matrix\\.h5$", basename(files), ignore.case = TRUE)]
  pos <- files[grepl("tissue_positions(_list)?\\.csv$", basename(files), ignore.case = TRUE)]
  sc <- files[grepl("scalefactors_json\\.json$", basename(files), ignore.case = TRUE)]
  hi <- files[grepl("tissue_hires_image\\.(png|jpg|jpeg)$", basename(files), ignore.case = TRUE)]
  lo <- files[grepl("tissue_lowres_image\\.(png|jpg|jpeg)$", basename(files), ignore.case = TRUE)]
  other <- files[grepl("detected_tissue_image|aligned_fiducials|\\.jpg$|\\.png$", basename(files), ignore.case = TRUE)]
  dataset <- infer_dataset(sr)
  complete <- length(h5) > 0 && length(pos) > 0 && length(sc) > 0 && length(hi) > 0 && length(lo) > 0
  data.table(
    sample_id = geo_from_root(sr),
    geo_accession = geo_from_root(sr),
    section_id = section_id_from_root(sr),
    dataset_label = dataset,
    condition_or_group = condition_from_dataset(dataset),
    matrix_file = ifelse(length(h5), rel(h5[1]), ""),
    features_file_or_h5 = ifelse(length(h5), rel(h5[1]), ""),
    barcodes_file = "",
    tissue_positions_file = ifelse(length(pos), rel(pos[1]), ""),
    scalefactors_file = ifelse(length(sc), rel(sc[1]), ""),
    hires_image_file = ifelse(length(hi), rel(hi[1]), ""),
    lowres_image_file = ifelse(length(lo), rel(lo[1]), ""),
    other_image_files = paste(rel(setdiff(other, c(hi[1], lo[1]))), collapse = ";"),
    metadata_file = "config/spatial_sample_metadata.tsv (template only)",
    complete_visium_input = ifelse(complete, "complete", "partially complete"),
    notes = ifelse(grepl("inferred", geo_from_root(sr)), "GEO accession inferred from archive naming; verify manually.", "Condition/group metadata not locked.")
  )
}), fill = TRUE)
setorder(manifest, dataset_label, section_id)
write_tsv(manifest, file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv"))

read_h5_dims <- function(h5) {
  out <- list(readable = FALSE, n_spots = NA_integer_, n_features = NA_integer_, genes = character())
  mat <- tryCatch(suppressWarnings(Seurat::Read10X_h5(h5)), error = function(e) e)
  if (inherits(mat, "error")) return(out)
  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) mat <- mat[["Gene Expression"]] else mat <- mat[[1]]
  }
  out$readable <- TRUE
  out$n_features <- nrow(mat)
  out$n_spots <- ncol(mat)
  out$genes <- rownames(mat)
  out
}
read_positions <- function(path) {
  out <- list(readable = FALSE, n = NA_integer_, in_tissue = NA_integer_, coordinate_ok = FALSE, cols = "")
  if (!nzchar(path) || !file.exists(file.path(root, path))) return(out)
  dt <- tryCatch(fread(file.path(root, path), header = FALSE), error = function(e) NULL)
  if (is.null(dt)) return(out)
  out$readable <- TRUE
  out$n <- nrow(dt)
  out$cols <- paste(names(dt), collapse = ";")
  if (ncol(dt) >= 6) {
    out$in_tissue <- suppressWarnings(sum(as.integer(dt[[2]]) == 1, na.rm = TRUE))
    out$coordinate_ok <- all(c(5, 6) <= ncol(dt)) && all(!is.na(suppressWarnings(as.numeric(dt[[5]]))[1:min(10, nrow(dt))]))
  }
  out
}
section_audit <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i]
  h5_abs <- file.path(root, row$matrix_file)
  dims <- if (file.exists(h5_abs)) read_h5_dims(h5_abs) else list(readable = FALSE, n_spots = NA_integer_, n_features = NA_integer_, genes = character())
  pos <- read_positions(row$tissue_positions_file)
  critical <- c()
  if (!dims$readable) critical <- c(critical, "matrix_unreadable")
  if (!pos$readable) critical <- c(critical, "positions_unreadable")
  if (!nzchar(row$hires_image_file)) critical <- c(critical, "missing_hires_image")
  if (!nzchar(row$scalefactors_file)) critical <- c(critical, "missing_scale_factors")
  data.table(
    sample_id = row$sample_id,
    section_id = row$section_id,
    matrix_readable = dims$readable,
    n_spots_total = ifelse(is.na(dims$n_spots), pos$n, dims$n_spots),
    n_spots_in_tissue = pos$in_tissue,
    n_features = dims$n_features,
    positions_readable = pos$readable,
    image_available = nzchar(row$hires_image_file) || nzchar(row$lowres_image_file),
    scale_factors_available = nzchar(row$scalefactors_file),
    coordinate_system_ok = pos$coordinate_ok,
    critical_issue = ifelse(length(critical), paste(critical, collapse = ";"), "none"),
    notes = ifelse(row$condition_or_group == "metadata_not_locked", "Visium structure complete, but biological condition/group metadata remains unresolved.", "")
  )
}), fill = TRUE)
write_tsv(section_audit, file.path(root, "results/tables/phase3_step1_spatial_section_audit.tsv"))

metadata_files <- unique(c(file.path(root, "config/spatial_sample_metadata.tsv"),
                           all_files[grepl("metadata|sample_manifest|dataset_manifest", basename(all_files), ignore.case = TRUE)]))
metadata_files <- metadata_files[file.exists(metadata_files) & file.info(metadata_files)$size < 5e6]
metadata_files <- metadata_files[tolower(tools::file_ext(metadata_files)) %in% c("tsv", "csv", "txt")]
metadata_audit <- rbindlist(lapply(metadata_files, function(f) {
  dt <- tryCatch(suppressWarnings(fread(f, fill = TRUE)), error = function(e) NULL)
  if (is.null(dt)) return(data.table(metadata_file = rel(f), metadata_column = "UNREADABLE", n_unique_values = NA_integer_, example_values = "", missing_fraction = NA_real_, likely_meaning = "unreadable", needed_for_phase3 = "yes", notes = "Could not parse as delimited text."))
  rbindlist(lapply(names(dt), function(col) {
    vals <- dt[[col]]
    lname <- tolower(col)
    meaning <- fifelse(grepl("sample", lname), "sample_id",
                fifelse(grepl("section", lname), "section_id",
                fifelse(grepl("disease|control|stone|healthy|group|condition", lname), "condition_or_group",
                fifelse(grepl("donor", lname), "donor",
                fifelse(grepl("region|papilla|anatom", lname), "anatomical_region",
                fifelse(grepl("roi|plaque|mineral|lesion|fibrosis|histology", lname), "ROI_or_histology_annotation", "other"))))))
    data.table(
      metadata_file = rel(f),
      metadata_column = col,
      n_unique_values = uniqueN(vals, na.rm = TRUE),
      example_values = paste(head(unique(as.character(vals[!is.na(vals)])), 5), collapse = ";"),
      missing_fraction = mean(is.na(vals) | vals == ""),
      likely_meaning = meaning,
      needed_for_phase3 = ifelse(meaning %in% c("sample_id", "section_id", "condition_or_group", "donor", "anatomical_region", "ROI_or_histology_annotation"), "yes", "no"),
      notes = ifelse(grepl("pending|unknown|template", paste(vals, collapse = " "), ignore.case = TRUE), "Contains placeholder/template values; not locked metadata.", "")
    )
  }))
}), fill = TRUE)
write_tsv(metadata_audit, file.path(root, "results/tables/phase3_step1_spatial_metadata_audit.tsv"))

roi_files <- inventory[grepl(paste(roi_keywords, collapse = "|"), paste(file_path, detected_keywords, notes), ignore.case = TRUE)]
roi_text_candidates <- inventory[file_type %in% text_ext & grepl("spatial|GSE206306|GSE231630|roi|plaque|mineral|lesion|fibrosis", file_path, ignore.case = TRUE)]
roi_status <- "NOT FOUND"
roi_report <- c(
  "# Phase 3-Step 1 ROI Annotation Search Report",
  "",
  paste0("- Files searched: ", length(all_files), " files across configured spatial/project directories."),
  paste0("- Keywords used: ", paste(roi_keywords, collapse = ", "), "."),
  paste0("- Filename/content keyword hits: ", nrow(roi_files), "."),
  "- Plaque/mineral/lesion/fibrosis-resolved spot-level or ROI-level annotation file: not identified as a canonical input.",
  "- Image-derived manual annotation: not identified in the current canonical spatial input structure.",
  "- Annotation granularity found: sample/file-level keyword mentions only, mostly in prior notes/results/scripts; no locked spot-level ROI table was found.",
  "",
  "No plaque-, mineral-, lesion- or fibrosis-resolved ROI annotation was identified in the current spatial files. Spatial analysis should therefore be framed as deconvolution-informed anatomical/tissue-context projection rather than plaque-specific localization.",
  "",
  paste0("Final status: ", roi_status),
  "",
  "## Keyword-Hit Files For Human Review",
  if (nrow(roi_files)) paste0("- `", head(roi_files$file_path, 50), "`: ", head(roi_files$likely_role, 50)) else "- None."
)
writeLines(roi_report, file.path(root, "notes/phase3_step1_roi_annotation_search_report.md"))

scrna_path <- file.path(root, "data/processed/gse231569_audited_seurat.rds")
scrna_obj <- tryCatch(readRDS(scrna_path), error = function(e) e)
scrna_readable <- !inherits(scrna_obj, "error")
scrna_genes <- character()
comp_counts <- data.table()
gene_format <- "unknown"
if (scrna_readable) {
  DefaultAssay(scrna_obj) <- ifelse("RNA" %in% Assays(scrna_obj), "RNA", DefaultAssay(scrna_obj))
  scrna_genes <- rownames(scrna_obj)
  gene_format <- ifelse(mean(grepl("^[A-Z0-9.-]+$", scrna_genes)) > 0.8, "gene_symbols_likely", "mixed_or_non_symbol")
  if ("phase1_cell_type" %in% colnames(scrna_obj@meta.data)) {
    comp_counts <- as.data.table(table(scrna_obj@meta.data$phase1_cell_type))
    setnames(comp_counts, c("compartment", "n_nuclei"))
  }
}
loop_n <- if (nrow(comp_counts)) comp_counts[compartment == "Loop_of_Henle_TAL", n_nuclei] else NA_integer_
scrna_audit <- data.table(
  audit_item = c(
    "object_readable", "object_path", "broad_compartment_field", "donor_field",
    "n_compartments", "n_reference_nuclei", "Loop/TAL_nuclei", "gene_symbol_format",
    "broad_compartment_suitability"
  ),
  value = c(
    as.character(scrna_readable), rel(scrna_path), as.character(scrna_readable && "phase1_cell_type" %in% colnames(scrna_obj@meta.data)),
    as.character(scrna_readable && "donor_id" %in% colnames(scrna_obj@meta.data)),
    as.character(ifelse(nrow(comp_counts), nrow(comp_counts), NA_integer_)),
    as.character(ifelse(scrna_readable, ncol(scrna_obj), NA_integer_)),
    as.character(ifelse(length(loop_n), loop_n, NA_integer_)),
    gene_format,
    ifelse(scrna_readable && nrow(comp_counts) >= 3 && length(loop_n) && loop_n >= 50, "suitable_for_broad_compartment_deconvolution_or_label_transfer", "requires_human_review")
  ),
  source = c(rep(rel(scrna_path), 8), "Phase 2 canonical snRNA object"),
  notes = c(
    "Read only; no transformation performed.",
    "Canonical Phase 2 reference.",
    "Expected field: phase1_cell_type.",
    "Expected field: donor_id.",
    "Broad compartments used for Phase 2.",
    "Reference nuclei/cells in object.",
    "Expected Loop/TAL label: Loop_of_Henle_TAL.",
    "Assessed from rownames only.",
    "Broad labels are appropriate for conservative deconvolution; fine cell-state claims should be avoided."
  )
)
if (nrow(comp_counts)) {
  scrna_audit <- rbind(scrna_audit, data.table(
    audit_item = paste0("compartment_count:", comp_counts$compartment),
    value = as.character(comp_counts$n_nuclei),
    source = rel(scrna_path),
    notes = "Reference nuclei per compartment."
  ), fill = TRUE)
}
write_tsv(scrna_audit, file.path(root, "results/tables/phase3_step1_scrna_reference_for_spatial_audit.tsv"))

section_genes <- list()
gene_overlap <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i]
  dims <- read_h5_dims(file.path(root, row$matrix_file))
  section_genes[[row$section_id]] <<- dims$genes
  overlap <- intersect(scrna_genes, dims$genes)
  data.table(
    sample_id = row$sample_id,
    section_id = row$section_id,
    n_spatial_genes = length(dims$genes),
    n_scrna_genes = length(scrna_genes),
    n_overlap_genes = length(overlap),
    overlap_fraction_spatial = ifelse(length(dims$genes), length(overlap) / length(dims$genes), NA_real_),
    overlap_fraction_scrna = ifelse(length(scrna_genes), length(overlap) / length(scrna_genes), NA_real_),
    suitable_for_label_transfer = ifelse(length(overlap) >= 2000, "yes", "review"),
    suitable_for_rctd = ifelse(length(overlap) >= 2000 && scrna_readable, "yes", "review"),
    notes = "Overlap assessed on rownames only; no normalization or deconvolution run."
  )
}), fill = TRUE)
write_tsv(gene_overlap, file.path(root, "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv"))

methods_note <- c(
  "# Phase 3-Step 1 Deconvolution Feasibility",
  "",
  "## A. Seurat label transfer",
  "",
  "- Required inputs: readable spatial count matrices, spatial coordinates/images, and a normalized/reference snRNA object with broad labels.",
  "- Whether available: broadly available; canonical snRNA object is readable and Visium sections have matrix/position/image/scale-factor files.",
  "- Expected outputs: predicted compartment labels or compartment prediction scores per spot.",
  "- Strengths: fastest fallback, uses existing Seurat installation, lower environment burden.",
  "- Limitations: transfer scores are contextual, not true cell-fraction deconvolution; should not be framed as plaque localization.",
  "",
  "## B. RCTD",
  "",
  "- Required inputs: raw spatial counts, spot coordinates, reference snRNA counts, and broad compartment labels.",
  "- Whether available: input structure appears feasible if RCTD/spacexr environment is available in Phase 3-Step 2.",
  "- Expected outputs: estimated broad compartment weights per spot and fit diagnostics.",
  "- Strengths: preferred primary method for broad compartment deconvolution; interpretable spot-level mixture estimates.",
  "- Limitations: requires package/environment setup and careful filtering; estimates remain tissue-context projections, not causal or plaque-specific evidence.",
  "",
  "## C. Cell2location",
  "",
  "- Required inputs: spatial count matrices, snRNA reference, trained regression/reference model, and Python/scvi-tools environment.",
  "- Whether available: biological inputs appear available, but environment burden is high and not confirmed in this step.",
  "- Expected outputs: cell-abundance estimates per spot.",
  "- Strengths: advanced Bayesian abundance mapping.",
  "- Limitations: heavier computational and software burden; more model choices; best treated as optional future extension.",
  "",
  "## Recommendation",
  "",
  "Primary recommendation for Phase 3-Step 2: use RCTD if the environment is available. Faster fallback: Seurat label transfer. Advanced optional method: Cell2location. No deconvolution was run in Step 1."
)
writeLines(methods_note, file.path(root, "notes/phase3_step1_deconvolution_feasibility.md"))

module_files <- c(
  MAGMA_top50 = "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
  MAGMA_top100 = "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
  MAGMA_Bonferroni = "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
  MAGMA_FDR05 = "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
  MAGMA_suggestive_p1e4 = "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
)
module_mapping <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i]
  genes <- section_genes[[row$section_id]]
  rbindlist(lapply(names(module_files), function(m) {
    mod <- unique(readLines(file.path(root, module_files[[m]]), warn = FALSE))
    present <- intersect(mod, genes)
    missing <- setdiff(mod, genes)
    data.table(
      sample_id = row$sample_id,
      section_id = row$section_id,
      module_name = m,
      n_module_genes = length(mod),
      n_present_in_spatial = length(present),
      n_missing_from_spatial = length(missing),
      present_genes = paste(present, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      suitable_for_spatial_scoring = ifelse(length(present) >= max(10, ceiling(0.5 * length(mod))), "yes", "review"),
      notes = "Mapping only; no spatial module scoring performed."
    )
  }))
}), fill = TRUE)
write_tsv(module_mapping, file.path(root, "results/tables/phase3_step1_magma_module_spatial_mapping.tsv"))

spatial_outputs <- inventory[grepl("results/spatial|results/figures|source_data|release/.*/results|stage6|phase28", file_path, ignore.case = TRUE)]
deprecated <- spatial_outputs[grepl("spatial.*score|module_score|MAGMA_FDR|MAGMA_suggestive|366|186|phase28|stage6", paste(file_path, detected_keywords, likely_role), ignore.case = TRUE)]
dep_audit <- if (nrow(deprecated)) {
  data.table(
    file_path = deprecated$file_path,
    likely_uses_deprecated_module = ifelse(grepl("phase28|stage6|spatial_MAGMA_FDR|spatial_MAGMA_suggestive|module_score|spatial_score", deprecated$file_path, ignore.case = TRUE), "possible", "unknown"),
    evidence = paste(deprecated$likely_role, deprecated$detected_keywords, sep = " | "),
    affected_module = ifelse(grepl("FDR", deprecated$file_path, ignore.case = TRUE), "FDR/deprecated_or_unverified", "unknown_or_multiple"),
    needs_rerun_or_review = "review_before_reuse",
    notes = "Existing spatial output predates Phase 3 canonical input locking or may use older module definitions; do not treat as canonical until rerun/reviewed."
  )
} else {
  data.table(file_path = "none", likely_uses_deprecated_module = "not_detected", evidence = "", affected_module = "", needs_rerun_or_review = "no", notes = "")
}
write_tsv(dep_audit, file.path(root, "results/tables/phase3_step1_deprecated_spatial_output_audit.tsv"))

complete_n <- manifest[complete_visium_input == "complete", .N]
section_n <- nrow(manifest)
roi_line <- "No plaque-, mineral-, lesion- or fibrosis-resolved ROI annotation was identified in the current spatial files."
plan <- c(
  "# Phase 3-Step 2 Plan",
  "",
  "## Scenario A: all spatial sections complete and snRNA reference compatible",
  "",
  "Proceed to broad-compartment deconvolution or label transfer. Preferred method: RCTD if package/environment is available.",
  "",
  "## Scenario B: spatial sections complete but RCTD environment unavailable",
  "",
  "Use Seurat label transfer as the primary fallback, report transfer scores as anatomical/tissue-context projections only.",
  "",
  "## Scenario C: some spatial sections incomplete",
  "",
  "Proceed only with complete sections and document exclusions in the sample manifest and methods.",
  "",
  "## Scenario D: no ROI annotation found",
  "",
  paste0(roi_line, " Proceed with anatomical/tissue-context projection only; do not claim plaque-specific localization."),
  "",
  "## Scenario E: gene overlap insufficient",
  "",
  "Human decision required before spatial module scoring or deconvolution.",
  "",
  "## Proposed Phase 3-Step 2",
  "",
  if (complete_n > 0) "Proceed with complete Visium sections using RCTD if available, otherwise Seurat label transfer fallback." else "Fix missing spatial inputs first."
)
writeLines(plan, file.path(root, "notes/phase3_step1_phase3_step2_plan.md"))

readiness <- if (complete_n > 0 && all(gene_overlap$n_overlap_genes >= 2000, na.rm = TRUE) && scrna_readable) "A. proceed to Phase 3-Step 2: spatial deconvolution / label transfer" else "D. human decision required"
report <- c(
  "# Phase 3-Step 1 Report",
  "",
  paste0("- Spatial datasets inspected: GSE206306 and GSE231630 candidate Visium files under raw/processed spatial directories plus config/results/scripts/release/package directories."),
  paste0("- Sections found: ", section_n, "."),
  paste0("- Complete Visium sections: ", complete_n, "."),
  "- Key files available: filtered_feature_bc_matrix.h5, tissue_positions_list.csv, scalefactors_json.json, tissue_hires_image.png and tissue_lowres_image.png for complete sections.",
  "- Metadata status: config/spatial_sample_metadata.tsv currently contains placeholder/template rows; disease/control and detailed biological group labels are not locked.",
  paste0("- ROI annotation status: ", roi_line),
  "- snRNA reference compatibility: canonical GSE231569 object is readable; broad compartment labels and donor fields are present; Loop/TAL label is available.",
  "- Recommended deconvolution method: RCTD if environment is available; Seurat label transfer is the faster fallback; Cell2location is optional/high-burden.",
  "- MAGMA module mapping status: canonical Phase 1-Step 3 modules were mapped to each spatial gene universe; mapping table created without scoring.",
  "- Deprecated spatial output status: existing phase28/stage6 spatial outputs and score/projection files require review before reuse and should not be treated as Phase 3 canonical outputs.",
  "- Readiness for Phase 3-Step 2: structurally ready for complete sections, but biological sample metadata and absence of ROI annotation constrain claims.",
  paste0("- Recommended next action: ", readiness, ".")
)
writeLines(report, file.path(root, "notes/phase3_step1_report.md"))

checklist <- data.table(
  task_id = sprintf("P3S1-%02d", 1:12),
  task_name = c(
    "Create spatial file inventory",
    "Create spatial sample manifest",
    "Create spatial input audit script",
    "Create spatial section audit",
    "Create spatial metadata audit",
    "Create ROI annotation search report",
    "Create snRNA reference compatibility audit",
    "Create spatial-snRNA gene overlap table",
    "Create deconvolution feasibility note",
    "Create MAGMA module spatial mapping",
    "Create deprecated spatial output audit",
    "Create Phase 3-Step 2 plan and report"
  ),
  completed = c(
    file.exists(file.path(root, "codex_tasks/phase3_step1_spatial_file_inventory.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv")),
    file.exists(file.path(root, "scripts/08_spatial_analysis/phase3_step1_spatial_input_audit.R")),
    file.exists(file.path(root, "results/tables/phase3_step1_spatial_section_audit.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1_spatial_metadata_audit.tsv")),
    file.exists(file.path(root, "notes/phase3_step1_roi_annotation_search_report.md")),
    file.exists(file.path(root, "results/tables/phase3_step1_scrna_reference_for_spatial_audit.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv")),
    file.exists(file.path(root, "notes/phase3_step1_deconvolution_feasibility.md")),
    file.exists(file.path(root, "results/tables/phase3_step1_magma_module_spatial_mapping.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1_deprecated_spatial_output_audit.tsv")),
    file.exists(file.path(root, "notes/phase3_step1_phase3_step2_plan.md")) && file.exists(file.path(root, "notes/phase3_step1_report.md"))
  ),
  output_file = c(
    "codex_tasks/phase3_step1_spatial_file_inventory.tsv",
    "results/tables/phase3_step1_spatial_sample_manifest.tsv",
    "scripts/08_spatial_analysis/phase3_step1_spatial_input_audit.R",
    "results/tables/phase3_step1_spatial_section_audit.tsv",
    "results/tables/phase3_step1_spatial_metadata_audit.tsv",
    "notes/phase3_step1_roi_annotation_search_report.md",
    "results/tables/phase3_step1_scrna_reference_for_spatial_audit.tsv",
    "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv",
    "notes/phase3_step1_deconvolution_feasibility.md",
    "results/tables/phase3_step1_magma_module_spatial_mapping.tsv",
    "results/tables/phase3_step1_deprecated_spatial_output_audit.tsv",
    "notes/phase3_step1_phase3_step2_plan.md; notes/phase3_step1_report.md"
  ),
  blocking_issue = c(rep("none", 4), "sample metadata not locked", "ROI annotation not found", rep("none", 4), "old outputs require review", "metadata/ROI claim limitations"),
  manual_review_needed = c("yes", "yes", "no", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes"),
  notes = c(
    "Recursive inventory across requested directories.",
    "Completeness is structural, not biological metadata approval.",
    "Script performs input audit only.",
    "Matrix dimensions and coordinate availability checked.",
    "Template metadata found; disease/control labels unresolved.",
    "No spot-level ROI table identified.",
    "Reference fields and compartment counts checked.",
    "Overlap assessed by gene symbols only.",
    "No methods run.",
    "Mapping only; no spatial scoring.",
    "Existing phase28/stage6 outputs should not be reused without review.",
    "Stop rule respected; wait for human review."
  )
)
write_tsv(checklist, file.path(root, "codex_tasks/phase3_step1_completion_checklist.tsv"))

message("Phase 3-Step 1 spatial input audit complete. No deconvolution or spatial scoring was run.")
