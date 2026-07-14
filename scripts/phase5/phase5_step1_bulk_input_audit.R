suppressPackageStartupMessages({
  library(tools)
})

root <- normalizePath(".", mustWork = TRUE)
dir.create(file.path(root, "codex_tasks"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "notes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "results", "tables"), showWarnings = FALSE, recursive = TRUE)

tsv_write <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_exists <- function(path) file.exists(file.path(root, path))
rel <- function(path) sub(paste0("^", gsub("([\\W])", "\\\\\\1", root), "/?"), "", normalizePath(path, mustWork = FALSE))
collapse_vec <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
yn <- function(x) ifelse(isTRUE(x), "yes", "no")

read_tsv_safe <- function(path, n_max = Inf) {
  full <- file.path(root, path)
  if (!file.exists(full)) return(NULL)
  tryCatch(
    read.delim(full, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, nrows = n_max, quote = "", comment.char = ""),
    error = function(e) NULL
  )
}

read_lines_safe <- function(path, n = 60) {
  full <- file.path(root, path)
  if (!file.exists(full)) return(character())
  con <- NULL
  out <- character()
  tryCatch({
    con <- if (grepl("\\.gz$", full, ignore.case = TRUE)) gzfile(full, "rt") else file(full, "rt")
    out <- readLines(con, n = n, warn = FALSE)
  }, error = function(e) character(), finally = if (!is.null(con)) close(con))
  out
}

count_rows_cols <- function(path) {
  full <- file.path(root, path)
  if (!file.exists(full)) return(list(rows = NA_integer_, cols = NA_integer_, header = ""))
  if (!grepl("\\.(tsv|txt|csv|log|md|R|sh)(\\.gz)?$", full, ignore.case = TRUE)) {
    return(list(rows = NA_integer_, cols = NA_integer_, header = ""))
  }
  header <- read_lines_safe(path, 1)
  cols <- if (length(header)) length(strsplit(header[[1]], "\t", fixed = TRUE)[[1]]) else NA_integer_
  rows <- NA_integer_
  if (grepl("\\.gz$", full, ignore.case = TRUE)) {
    con <- NULL
    rows <- tryCatch({
      con <- gzfile(full, "rt")
      length(readLines(con, warn = FALSE)) - 1L
    }, error = function(e) NA_integer_, finally = if (!is.null(con)) close(con))
  } else if (grepl("\\.(tsv|txt|csv)$", full, ignore.case = TRUE)) {
    rows <- tryCatch(length(readLines(full, warn = FALSE)) - 1L, error = function(e) NA_integer_)
  }
  list(rows = rows, cols = cols, header = paste(head(strsplit(ifelse(length(header), header, ""), "\t", fixed = TRUE)[[1]], 12), collapse = ";"))
}

keywords <- c("GSE73680","bulk","microarray","RNA","expression","plaque","stone","papilla",
              "Randall","patient","paired","control","disease","adjacent","metadata","sample",
              "GSM","module_score","MAGMA","injury","remodeling","epithelial","composition",
              "deconvolution","MuSiC","Bisque","CIBERSORT","CIBERSORTx","cell fraction",
              "marker score","patient fixed effect","paired model","Figure4","source data")

scan_roots <- c("data/raw", "data/processed", "results", "scripts", "config", "release",
                Sys.glob(file.path(root, "supplementary_package*")),
                Sys.glob(file.path(root, "submission_package*")))
scan_roots <- unique(scan_roots[file.exists(scan_roots)])
all_files <- unique(unlist(lapply(scan_roots, function(d) list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE))))
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]

detect_keywords <- function(path) {
  path_text <- rel(path)
  text <- path_text
  ext <- tolower(file_ext(path))
  size <- file.info(path)$size
  if (ext %in% c("tsv","txt","csv","md","r","sh","log","gz") && !is.na(size) && size < 5e6) {
    text <- paste(text, paste(read_lines_safe(rel(path), 40), collapse = "\n"))
  }
  hits <- keywords[vapply(keywords, function(k) grepl(k, text, ignore.case = TRUE, fixed = FALSE), logical(1))]
  unique(hits)
}

likely_role <- function(p) {
  b <- basename(p)
  x <- tolower(rel(p))
  if (grepl("gse73680_raw\\.tar$", x)) return("raw supplementary archive for GSE73680")
  if (grepl("gene_expression_matrix", x)) return("gene-level processed expression matrix")
  if (grepl("expression_matrix", x)) return("feature-level processed expression matrix")
  if (grepl("metadata_curated|sample_metadata_curated", x)) return("curated sample metadata")
  if (grepl("sample_metadata", x)) return("sample metadata")
  if (grepl("patient_sample_structure|pairing", x)) return("patient pairing/design table")
  if (grepl("feature_gene_mapping|gene_mapping", x)) return("feature-to-gene mapping/probe annotation")
  if (grepl("module_score|module_response|module_manifest", x)) return("existing bulk module-score/model output")
  if (grepl("injury|remodel", x)) return("injury/remodeling marker-score output")
  if (grepl("figure4|source_data", x)) return("Figure 4 or source-data artifact")
  if (grepl("composition|deconv|cibersort|music|bisque|cell_fraction", x)) return("composition/deconvolution-related artifact")
  if (grepl("\\.r$", tolower(b))) return("analysis script")
  "supporting project file"
}

candidate_input <- function(p) {
  x <- tolower(rel(p))
  grepl("gse73680_raw\\.tar|gse73680_gene_expression_matrix|gse73680_expression_matrix|sample_metadata_curated|sample_metadata|patient_sample_structure|feature_gene_mapping|gene_mapping|gse231569_audited_seurat|magma_", x)
}

inventory_hits <- lapply(all_files, function(p) {
  hits <- detect_keywords(p)
  if (!length(hits) && !candidate_input(p)) return(NULL)
  info <- file.info(p)
  data.frame(
    file_path = rel(p),
    file_name = basename(p),
    file_type = ifelse(grepl("\\.gz$", p), paste0(file_ext(sub("\\.gz$", "", p)), ".gz"), file_ext(p)),
    file_size = info$size,
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    likely_role = likely_role(p),
    detected_keywords = paste(hits, collapse = ";"),
    candidate_for_canonical_input = yn(candidate_input(p)),
    notes = ifelse(info$size == 0, "empty file", ""),
    stringsAsFactors = FALSE
  )
})
inventory <- do.call(rbind, inventory_hits)
if (is.null(inventory)) inventory <- data.frame(file_path=character(), file_name=character(), file_type=character(), file_size=numeric(), modified_time=character(), likely_role=character(), detected_keywords=character(), candidate_for_canonical_input=character(), notes=character())
inventory <- inventory[order(inventory$file_path), ]
tsv_write(inventory, file.path(root, "codex_tasks/phase5_step1_bulk_file_inventory.tsv"))

canonical_paths <- list(
  raw_expression_matrix = "data/raw/GSE73680/GSE73680_RAW.tar",
  processed_expression_matrix = "data/processed/gse73680/gse73680_expression_matrix.tsv.gz",
  normalized_expression_matrix = "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz",
  sample_metadata = "config/gse73680_sample_metadata.tsv",
  curated_sample_metadata = "config/gse73680_sample_metadata_curated.tsv",
  patient_pairing_table = "results/gse73680/tables/gse73680_patient_sample_structure.tsv",
  probe_annotation = "results/gse73680/tables/gse73680_feature_gene_mapping.tsv",
  gene_symbol_mapping = "results/gse73680/tables/gse73680_gene_mapping_qc.tsv",
  existing_bulk_module_score_table = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_sample_module_scores.tsv",
  existing_injury_remodeling_score_table = "results/gse73680/tables/gse73680_injury_program_score_matrix.tsv",
  existing_marker_score_table = "results/gse73680/tables/gse73680_injury_program_score_matrix.tsv",
  existing_bulk_model_results = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_response_model.tsv",
  existing_figure4_source_data = "results/tables/figure4_panel_source_files_final.tsv",
  existing_bulk_deconvolution_results = "",
  snRNA_reference_for_deconvolution = "data/processed/gse231569_audited_seurat.rds",
  canonical_magma_modules = "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt;results/phase1_step3_magma_gene_sets/MAGMA_top100.txt;results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt;results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt;results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
)

alts <- list(
  raw_expression_matrix = "data/raw/GSE73680/GSE73680_series_full.soft.txt;data/raw/GSE73680/GSE73680_suppl_filelist.txt",
  processed_expression_matrix = "data/processed/gse73680/gse73680_patient_level_expression_matrix.tsv.gz",
  normalized_expression_matrix = "data/processed/gse73680/gse73680_patient_level_expression_matrix.tsv.gz",
  patient_pairing_table = "results/tables/revision/stage5B1_gse73680_module_response/gse73680_stage5B1_pairing_summary.tsv",
  existing_bulk_module_score_table = "results/gse73680/tables/gse73680_module_score_matrix.tsv;results/gse73680/tables/gse73680_patient_level_module_score_matrix.tsv",
  existing_bulk_model_results = "results/gse73680/tables/gse73680_module_score_response.tsv;results/gse73680/tables/gse73680_patient_level_module_response.tsv;results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_patient_fixed_effect_adjusted_models.tsv",
  existing_figure4_source_data = "results/tables/revision/stage5C1_gse73680_figure4_draft/figure4_source_data_manifest_v0.1.tsv;results/figures/figure4_revisions/final_gse73680_disease_context_20260621/figure4_panel_source_files_final.tsv",
  existing_bulk_deconvolution_results = "results/tables/revision/stage5B2_gse73680_composition_adjusted/*"
)

manifest_rows <- lapply(names(canonical_paths), function(cat) {
  p <- canonical_paths[[cat]]
  p1 <- strsplit(p, ";", fixed = TRUE)[[1]]
  exists_vec <- nzchar(p1) & file.exists(file.path(root, p1))
  status <- if (cat == "existing_bulk_deconvolution_results") {
    "not_available"
  } else if (all(exists_vec)) {
    "locked"
  } else if (any(exists_vec)) {
    "partial"
  } else {
    "missing"
  }
  primary <- if (length(p1) == 1) p1 else paste(p1, collapse = ";")
  cr <- if (length(p1) == 1 && file.exists(file.path(root, p1))) count_rows_cols(p1) else list(rows = NA_integer_, cols = NA_integer_, header = "")
  info <- if (length(p1) == 1 && file.exists(file.path(root, p1))) file.info(file.path(root, p1)) else NULL
  concern <- switch(cat,
    existing_bulk_deconvolution_results = "No formal deconvolution output was identified; composition-adjusted historical tables are not deconvolution estimates.",
    existing_bulk_module_score_table = "Historical module-score outputs require rerun/review against canonical Phase 1 gene sets before manuscript use.",
    curated_sample_metadata = "Labels are filename-derived and should remain manually reviewable.",
    normalized_expression_matrix = "Matrix appears to contain processed/normalized intensity-like values, not raw counts.",
    snRNA_reference_for_deconvolution = "Requires method-specific object inspection and package availability before formal deconvolution.",
    canonical_magma_modules = "Use these Phase 1 locked modules only; historical frozen modules are deprecated.",
    ""
  )
  action <- if (status == "missing") "locate or regenerate before downstream use" else if (nzchar(concern)) "review before Step 2 use" else "ready for audit use"
  data.frame(
    input_category = cat,
    canonical_status = status,
    canonical_file_path = primary,
    alternative_file_paths = ifelse(!is.null(alts[[cat]]), alts[[cat]], ""),
    selection_reason = switch(cat,
      raw_expression_matrix = "Only local raw GSE73680 supplementary archive found.",
      processed_expression_matrix = "Feature-level matrix preserves original reconstructed feature identifiers.",
      normalized_expression_matrix = "Gene-level matrix is the practical canonical matrix for module mapping/scoring.",
      sample_metadata = "Original local sample metadata table.",
      curated_sample_metadata = "Curated table contains patient IDs, group labels and include flags.",
      patient_pairing_table = "Existing patient-level structure table summarizes paired design.",
      probe_annotation = "Feature-to-gene mapping documents direct symbol mapping and unmapped features.",
      gene_symbol_mapping = "Mapping QC summarizes available gene symbols.",
      existing_bulk_module_score_table = "Most recent revision-stage sample module score table found.",
      existing_injury_remodeling_score_table = "Existing injury/remodeling score matrix found.",
      existing_marker_score_table = "Existing marker-like injury/remodeling score matrix found.",
      existing_bulk_model_results = "Most recent revision-stage module response model table found.",
      existing_figure4_source_data = "Final local Figure 4 source manifest found.",
      existing_bulk_deconvolution_results = "No formal MuSiC/Bisque/CIBERSORTx output found; stage5B2 is adjustment, not deconvolution.",
      snRNA_reference_for_deconvolution = "Phase 2 locked snRNA reference requested by Step 5.1.",
      canonical_magma_modules = "Phase 1 locked module files requested by Step 5.1."
    ),
    file_size = if (!is.null(info)) info$size else NA,
    modified_time = if (!is.null(info)) format(info$mtime, "%Y-%m-%d %H:%M:%S") else "",
    readable = ifelse(status %in% c("locked","partial"), "yes", "no"),
    row_count = cr$rows,
    column_count = cr$cols,
    header_fields = cr$header,
    concerns = concern,
    action_needed = action,
    stringsAsFactors = FALSE
  )
})
manifest <- do.call(rbind, manifest_rows)
tsv_write(manifest, file.path(root, "codex_tasks/phase5_step1_bulk_canonical_input_manifest.tsv"))

expr_path <- "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz"
expr <- read_tsv_safe(expr_path)
meta_path <- "config/gse73680_sample_metadata_curated.tsv"
meta <- read_tsv_safe(meta_path)
mapping <- read_tsv_safe("results/gse73680/tables/gse73680_feature_gene_mapping.tsv")

gene_col <- if (!is.null(expr)) names(expr)[1] else NA_character_
sample_cols <- if (!is.null(expr)) setdiff(names(expr), gene_col) else character()
expr_num <- if (!is.null(expr)) suppressWarnings(as.matrix(data.frame(lapply(expr[sample_cols], as.numeric), check.names = FALSE))) else matrix(numeric(), 0, 0)
expr_values <- as.numeric(expr_num)
expr_values <- expr_values[is.finite(expr_values)]
scale_note <- if (!length(expr_values)) "not_inferable" else if (max(expr_values, na.rm = TRUE) > 1000) "processed_intensity_or_normalized_non_log_scale" else if (max(expr_values, na.rm = TRUE) <= 30 && min(expr_values, na.rm = TRUE) >= -5) "log_like_or_z_like" else "continuous_processed_expression"
meta_ids <- if (!is.null(meta) && "sample_id" %in% names(meta)) meta$sample_id else character()
matched <- intersect(sample_cols, meta_ids)
gene_values <- if (!is.null(expr)) expr[[gene_col]] else character()
nonempty_genes <- gene_values[!is.na(gene_values) & nzchar(gene_values)]
dup_features <- sum(duplicated(nonempty_genes))
missing_frac <- if (!length(expr_num)) NA_real_ else mean(is.na(expr_num))
probe_available <- !is.null(mapping) && all(c("feature_id","gene_symbol") %in% names(mapping))

expr_audit <- data.frame(
  audit_item = c("n_features","n_samples","expression_scale_inferred","gene_symbol_available",
                 "probe_annotation_available","duplicate_features","missing_value_fraction",
                 "sample_id_match_to_metadata","suitable_for_module_scoring","suitable_for_deconvolution"),
  value = c(
    ifelse(is.null(expr), NA, nrow(expr)),
    length(sample_cols),
    scale_note,
    yn(!is.null(expr) && gene_col == "gene"),
    yn(probe_available),
    dup_features,
    signif(missing_frac, 4),
    paste0(length(matched), "/", length(sample_cols)),
    ifelse(length(matched) == length(sample_cols) && length(nonempty_genes) > 1000, "yes_with_gene_filtering_and_z_scoring", "needs_review"),
    "limited_feasibility_normalized_microarray_not_raw_counts"
  ),
  source_file = expr_path,
  notes = c(
    "Rows are interpreted from first column of gene-level matrix.",
    "Sample columns are GSM identifiers.",
    "Large positive values indicate processed intensity-like matrix rather than raw count matrix.",
    "First column is named gene; one blank feature row is present and should be filtered.",
    "Feature-to-gene mapping file is available from prior reconstruction.",
    "Duplicate count excludes empty gene identifiers.",
    "Fraction across numeric expression cells.",
    "Compared expression GSM columns against curated metadata sample_id.",
    "Future scoring should use canonical MAGMA modules and filter missing/empty genes.",
    "Formal bulk deconvolution needs method-specific compatibility checks and sensitivity-only interpretation."
  ),
  stringsAsFactors = FALSE
)
tsv_write(expr_audit, file.path(root, "results/tables/phase5_step1_bulk_expression_audit.tsv"))

if (is.null(meta)) {
  meta_audit <- data.frame(sample_id=character(), geo_accession=character(), patient_id=character(), pair_id=character(), sample_group_raw=character(), sample_group_curated=character(), control_or_disease=character(), papilla_region=character(), stone_or_plaque_status=character(), paired_sample_available=character(), metadata_confidence=character(), metadata_source=character(), notes=character())
  pairing <- data.frame(patient_id=character(), pair_id=character(), n_samples=integer(), control_sample_id=character(), disease_sample_id=character(), pair_complete=character(), group_labels_confident=character(), notes=character())
} else {
  meta$include_bool <- toupper(as.character(meta$include_in_analysis)) == "TRUE"
  control <- meta$group_curated == "control_or_adjacent"
  disease <- meta$group_curated == "plaque_or_stone_papilla"
  included <- meta[meta$include_bool & (control | disease), , drop = FALSE]
  pair_complete_by_patient <- tapply(included$group_curated, included$patient_id, function(g) all(c("control_or_adjacent","plaque_or_stone_papilla") %in% g))
  meta_audit <- data.frame(
    sample_id = meta$sample_id,
    geo_accession = meta$geo_accession,
    patient_id = meta$patient_id,
    pair_id = meta$patient_id,
    sample_group_raw = meta$raw_label,
    sample_group_curated = meta$group_curated,
    control_or_disease = ifelse(meta$group_curated == "control_or_adjacent", "control_or_adjacent",
                                ifelse(meta$group_curated == "plaque_or_stone_papilla", "disease_context", "unclear_or_excluded")),
    papilla_region = meta$tissue_site,
    stone_or_plaque_status = ifelse(nzchar(meta$plaque_status), meta$plaque_status, meta$stone_status),
    paired_sample_available = ifelse(meta$patient_id %in% names(pair_complete_by_patient)[pair_complete_by_patient], "yes", "no"),
    metadata_confidence = ifelse(meta$include_bool & (control | disease), "usable_filename_derived_needs_manual_review", "low_or_excluded"),
    metadata_source = meta$sample_type,
    notes = paste0("include_in_analysis=", meta$include_in_analysis, ifelse(nzchar(meta$exclude_reason), paste0("; exclude_reason=", meta$exclude_reason), "")),
    stringsAsFactors = FALSE
  )
  split_inc <- split(included, included$patient_id)
  pairing <- do.call(rbind, lapply(split_inc, function(d) {
    ctl <- d$sample_id[d$group_curated == "control_or_adjacent"]
    dis <- d$sample_id[d$group_curated == "plaque_or_stone_papilla"]
    data.frame(
      patient_id = d$patient_id[1],
      pair_id = d$patient_id[1],
      n_samples = nrow(d),
      control_sample_id = collapse_vec(ctl),
      disease_sample_id = collapse_vec(dis),
      pair_complete = yn(length(ctl) >= 1 && length(dis) >= 1),
      group_labels_confident = "usable_but_filename_derived",
      notes = ifelse(nrow(d) == 2 && length(ctl) == 1 && length(dis) == 1, "one control/adjacent and one disease-context sample", "incomplete or non-standard pairing"),
      stringsAsFactors = FALSE
    )
  }))
  pairing <- pairing[order(pairing$patient_id), ]
}
tsv_write(meta_audit, file.path(root, "results/tables/phase5_step1_bulk_metadata_audit.tsv"))
tsv_write(pairing, file.path(root, "results/tables/phase5_step1_bulk_pairing_summary.tsv"))

module_files <- c(
  MAGMA_top50 = "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
  MAGMA_top100 = "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
  MAGMA_Bonferroni = "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
  MAGMA_FDR05 = "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
  MAGMA_suggestive_p1e4 = "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
)
bulk_genes <- unique(nonempty_genes)
map_rows <- lapply(names(module_files), function(m) {
  path <- module_files[[m]]
  genes <- if (file.exists(file.path(root, path))) unique(readLines(file.path(root, path), warn = FALSE)) else character()
  genes <- genes[nzchar(genes)]
  present <- intersect(genes, bulk_genes)
  missing <- setdiff(genes, bulk_genes)
  data.frame(
    module_name = m,
    module_file = path,
    n_module_genes = length(genes),
    n_present_in_bulk = length(present),
    n_missing_from_bulk = length(missing),
    present_genes = paste(present, collapse = ";"),
    missing_genes = paste(missing, collapse = ";"),
    suitable_for_bulk_scoring = ifelse(length(present) >= 10 && length(present) / max(1, length(genes)) >= 0.5, "yes", "limited_or_no"),
    notes = "Presence assessed against canonical gene-level GSE73680 matrix only; no scoring was run.",
    stringsAsFactors = FALSE
  )
})
module_mapping <- do.call(rbind, map_rows)
tsv_write(module_mapping, file.path(root, "results/tables/phase5_step1_magma_module_bulk_mapping.tsv"))

bulk_result_candidates <- inventory$file_path[grepl("gse73680|figure4|stage5|bulk|plaque", inventory$file_path, ignore.case = TRUE)]
bulk_result_candidates <- bulk_result_candidates[grepl("\\.(tsv|txt|csv|md|R|log)$", bulk_result_candidates, ignore.case = TRUE)]
bulk_result_candidates <- unique(bulk_result_candidates)
dep_rows <- lapply(bulk_result_candidates, function(p) {
  lines <- read_lines_safe(p, 250)
  txt <- paste(lines, collapse = "\n")
  evidence <- character()
  if (grepl("\\b366\\b", txt)) evidence <- c(evidence, "contains_366")
  if (grepl("\\b186\\b", txt)) evidence <- c(evidence, "contains_186")
  if (grepl("MAGMA_FDR|FDR05|suggestive|p1e-4|p1e4", txt, ignore.case = TRUE)) evidence <- c(evidence, "mentions_FDR_or_suggestive_module")
  if (grepl("P1_core|audited_locus|old|deprecated", txt, ignore.case = TRUE)) evidence <- c(evidence, "mentions_noncanonical_or_deprecated_context")
  uses <- length(evidence) > 0 && grepl("module|figure4|stage5|gse73680|magma|score", p, ignore.case = TRUE)
  data.frame(
    file_path = p,
    likely_uses_deprecated_module = ifelse(uses, "possible", "no_evidence_in_scanned_header"),
    evidence = paste(evidence, collapse = ";"),
    affected_module = ifelse(grepl("FDR", paste(p, txt), ignore.case = TRUE), "FDR_or_suggestive_MAGMA_related", ifelse(uses, "unknown_module", "")),
    needs_rerun_or_review = ifelse(uses, "yes_review_before_use", "no_immediate_action"),
    notes = "Header/early-line scan only; historical bulk/Figure4 outputs should not be reused until canonical Phase 1 modules are confirmed.",
    stringsAsFactors = FALSE
  )
})
dep_audit <- do.call(rbind, dep_rows)
dep_audit <- dep_audit[order(dep_audit$needs_rerun_or_review, dep_audit$file_path), ]
tsv_write(dep_audit, file.path(root, "results/tables/phase5_step1_deprecated_bulk_output_audit.tsv"))

marker_files <- inventory$file_path[grepl("injury|remodel|marker|signature|composition|epithelial|cell_fraction|deconv|cibersort|music|bisque", inventory$file_path, ignore.case = TRUE)]
marker_files <- unique(marker_files[grepl("\\.(tsv|txt|csv|md|R|log)$", marker_files, ignore.case = TRUE)])
marker_rows <- lapply(marker_files, function(p) {
  cr <- count_rows_cols(p)
  data.frame(
    marker_or_score_name = sub("\\.[^.]+$", "", basename(p)),
    file_path = p,
    n_genes_or_samples = cr$rows,
    likely_role = likely_role(p),
    used_in_previous_analysis = ifelse(grepl("results|release|submission|supplementary", p), "yes_or_exported", "script_or_config"),
    concern = ifelse(grepl("composition|adjusted", p, ignore.case = TRUE), "May be adjustment/signature output rather than formal deconvolution.", "Requires canonical-module and label review before reuse."),
    recommended_use = "audit_resource_only_in_step5_1; do_not_use_as_new_evidence_without_step5_2_review",
    notes = paste0("Header fields: ", cr$header),
    stringsAsFactors = FALSE
  )
})
marker_audit <- if (length(marker_rows)) do.call(rbind, marker_rows) else data.frame(marker_or_score_name=character(), file_path=character(), n_genes_or_samples=integer(), likely_role=character(), used_in_previous_analysis=character(), concern=character(), recommended_use=character(), notes=character())
tsv_write(marker_audit, file.path(root, "results/tables/phase5_step1_bulk_marker_resource_audit.tsv"))

n_total_meta <- if (!is.null(meta)) nrow(meta) else NA_integer_
n_included <- if (!is.null(meta)) sum(meta$include_bool & meta$group_curated %in% c("control_or_adjacent","plaque_or_stone_papilla")) else NA_integer_
n_patients <- if (!is.null(meta)) length(unique(meta$patient_id[meta$include_bool & meta$group_curated %in% c("control_or_adjacent","plaque_or_stone_papilla")])) else NA_integer_
n_pairs <- if (exists("pairing") && nrow(pairing)) sum(pairing$pair_complete == "yes") else 0L
metadata_note <- paste0("Curated metadata contains ", n_total_meta, " rows; included disease-context design has ", n_included,
                        " samples, ", n_patients, " patients and ", n_pairs, " complete paired patients.")

writeLines(c(
  "# Phase 5-Step 1 Bulk Deconvolution Feasibility",
  "",
  "No bulk deconvolution was run in this audit.",
  "",
  "## A. MuSiC",
  "- Required inputs: bulk expression matrix, single-cell reference expression, broad cell labels, donor labels.",
  "- Available: gene-level GSE73680 matrix, `data/processed/gse231569_audited_seurat.rds`, requested fields `phase1_cell_type` and `donor_id` for inspection.",
  "- Strength: donor-aware reference design is conceptually aligned with the available snRNA reference.",
  "- Limitation: GSE73680 is microarray/processed intensity-like data rather than RNA-seq counts; gene-platform overlap and Seurat assay/count availability must be verified before use.",
  "",
  "## B. BisqueRNA",
  "- Required inputs: bulk matrix and single-cell reference with sample/donor labels.",
  "- Available: same as MuSiC.",
  "- Strength: can estimate broad cell-state proportions from sc/sn reference.",
  "- Limitation: package availability, microarray compatibility and normalization assumptions require a dry-run feasibility check before formal analysis.",
  "",
  "## C. CIBERSORTx",
  "- Required inputs: signature matrix or external CIBERSORTx setup plus bulk mixture matrix.",
  "- Available: broad snRNA reference may support a signature matrix later, but no local CIBERSORTx output/setup was identified.",
  "- Strength: familiar framework for mixture estimation.",
  "- Limitation: external setup/licensing and cross-platform signature construction make it less suitable as an immediate local primary method.",
  "",
  "## D. Marker-score sensitivity only",
  "- Required inputs: marker gene sets and gene-level bulk expression.",
  "- Available: prior injury/remodeling and signature outputs exist, but should be reviewed against canonical inputs.",
  "- Strength: transparent, local, and less dependent on raw count assumptions.",
  "- Limitation: marker scores are not cell fractions and should be reported only as sensitivity/tissue-state proxies.",
  "",
  "## Recommendation",
  "E. human decision required.",
  "",
  "Reason: MuSiC or BisqueRNA may be feasible after inspecting the snRNA object assays and cross-platform gene overlap, but the current audit does not establish raw-count compatibility. Marker-score sensitivity is immediately feasible as a conservative backup. Any deconvolution should be interpreted as sensitivity analysis, not direct histological cell fractions."
), con = file.path(root, "notes/phase5_step1_bulk_deconvolution_feasibility.md"))

writeLines(c(
  "# Phase 5-Step 1 Bulk Model Plan",
  "",
  "Future analyses only; no model was run in Step 5.1.",
  "",
  "Base paired model:",
  "`module_score ~ group_binary + factor(patient_id)`",
  "",
  "Composition-adjusted sensitivity:",
  "`module_score ~ group_binary + estimated_cell_fraction + factor(patient_id)`",
  "",
  "Tissue-state sensitivity:",
  "`module_score ~ group_binary + injury_score + remodeling_score + factor(patient_id)`",
  "",
  "Compact model:",
  "`module_score ~ group_binary + selected composition proxy + selected injury/remodeling score + factor(patient_id)`",
  "",
  "Interpretation boundary: attenuation after adjustment may mean disease-state embedding rather than simple confounding. Adjustment should be treated as sensitivity analysis, with overfitting avoided because the paired design has limited sample size."
), con = file.path(root, "notes/phase5_step1_bulk_model_plan.md"))

writeLines(c(
  "# Phase 5-Step 1 Bulk Claim Boundary",
  "",
  "Safe wording:",
  "",
  "> Bulk plaque/stone-associated papillary data were used to assess whether MAGMA-prioritized modules were embedded within disease-associated tissue-state programs.",
  "",
  "Unsafe wording:",
  "- Bulk validates causal genes.",
  "- Bulk proves plaque mechanism.",
  "- Deconvolution proves cell fractions.",
  "- Adjusted loss of significance disproves relevance.",
  "",
  "Operational boundary: bulk evidence must not upgrade Figure 3 candidate groups until Phase 5 outputs are reviewed."
), con = file.path(root, "notes/phase5_step1_bulk_claim_boundary.md"))

proceed_rec <- if (!is.na(n_included) && n_included >= 40 && n_pairs >= 20 && all(module_mapping$suitable_for_bulk_scoring == "yes")) {
  "A. proceed to Phase 5-Step 2: bulk module scoring and paired disease-context model using canonical MAGMA modules"
} else {
  "D. human decision required"
}

writeLines(c(
  "# Phase 5-Step 1 Report",
  "",
  "## What was inspected",
  "Local `data/raw/`, `data/processed/`, `results/`, `scripts/`, `config/`, `release/`, `supplementary_package*`, and `submission_package*` files were scanned for GSE73680/bulk/plaque/module/deconvolution/Figure4 resources.",
  "",
  "## Canonical expression and metadata",
  paste0("- Canonical gene-level expression matrix: `", expr_path, "`."),
  paste0("- Expression dimensions: ", ifelse(is.null(expr), "unreadable", paste0(nrow(expr), " features x ", length(sample_cols), " samples")), "."),
  paste0("- Curated metadata: `", meta_path, "`."),
  paste0("- ", metadata_note),
  "- Group labels are usable for planning but remain filename-derived and should be manually reviewable.",
  "",
  "## Module mapping",
  paste0("- Canonical MAGMA module mapping table: `results/tables/phase5_step1_magma_module_bulk_mapping.tsv`. All assessed modules marked: ", paste(unique(module_mapping$suitable_for_bulk_scoring), collapse = ";"), "."),
  "",
  "## Deprecated output audit",
  "- Historical bulk/module/Figure4 outputs exist and should be reviewed before reuse because old frozen FDR/suggestive modules are deprecated.",
  "- The audit table flags files with possible FDR/suggestive/deprecated evidence from path and early-line scans.",
  "",
  "## Deconvolution feasibility",
  "- No formal bulk deconvolution output was identified.",
  "- MuSiC/BisqueRNA may be feasible only after snRNA assay/count and gene-overlap checks; CIBERSORTx requires external setup.",
  "- Marker-score sensitivity is the most immediately conservative backup.",
  "",
  "## Readiness",
  paste0("- Paired base model can proceed after human acceptance of filename-derived labels: ", ifelse(n_pairs >= 20, "yes", "needs review"), "."),
  "- Formal deconvolution should not proceed until method choice and snRNA assay compatibility are confirmed.",
  paste0("- Recommended next: ", proceed_rec, ".")
), con = file.path(root, "notes/phase5_step1_report.md"))

checklist <- data.frame(
  task_id = sprintf("P5S1-%02d", 1:12),
  task_name = c("bulk file inventory","canonical input manifest","expression matrix audit","metadata audit","pairing summary","deprecated bulk output audit","MAGMA module bulk mapping","marker resource audit","deconvolution feasibility note","bulk model plan note","claim boundary note","final report"),
  completed = rep("yes", 12),
  output_file = c(
    "codex_tasks/phase5_step1_bulk_file_inventory.tsv",
    "codex_tasks/phase5_step1_bulk_canonical_input_manifest.tsv",
    "results/tables/phase5_step1_bulk_expression_audit.tsv",
    "results/tables/phase5_step1_bulk_metadata_audit.tsv",
    "results/tables/phase5_step1_bulk_pairing_summary.tsv",
    "results/tables/phase5_step1_deprecated_bulk_output_audit.tsv",
    "results/tables/phase5_step1_magma_module_bulk_mapping.tsv",
    "results/tables/phase5_step1_bulk_marker_resource_audit.tsv",
    "notes/phase5_step1_bulk_deconvolution_feasibility.md",
    "notes/phase5_step1_bulk_model_plan.md",
    "notes/phase5_step1_bulk_claim_boundary.md",
    "notes/phase5_step1_report.md"
  ),
  blocking_issue = c(rep("", 8), "formal deconvolution method not selected", "", "", "human review needed before Step 2"),
  manual_review_needed = c(rep("yes", 12)),
  notes = c(
    "Recursive scan completed without deleting or overwriting existing files.",
    "Canonical status recorded for all required categories.",
    "No scoring, model or deconvolution was run.",
    "All curated metadata rows reported.",
    "Pairing summary uses included control/disease-context samples.",
    "Historical outputs require review before reuse.",
    "Presence only; no module score was computed.",
    "Marker/signature resources inventoried only.",
    "Recommendation is human decision required for formal deconvolution.",
    "Future model formulas documented only.",
    "Safe and unsafe claims documented.",
    "Stop rule satisfied after audit outputs."
  ),
  stringsAsFactors = FALSE
)
tsv_write(checklist, file.path(root, "codex_tasks/phase5_step1_completion_checklist.tsv"))

message("Phase 5-Step 1 bulk input audit completed.")
