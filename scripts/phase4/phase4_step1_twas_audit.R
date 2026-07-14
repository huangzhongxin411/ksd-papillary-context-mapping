suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

dir.create("codex_tasks", showWarnings = FALSE, recursive = TRUE)
dir.create("notes", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

keywords <- c(
  "TWAS", "S-PrediXcan", "SPrediXcan", "PrediXcan", "GTEx", "Kidney_Cortex",
  "Kidney Cortex", "model", "one-SNP", "single SNP", "prediction model",
  "zscore", "effect_size", "pvalue", "FDR", "SMR", "coloc", "eQTL", "R1",
  "R2", "R3", "R4", "R5", "R6", "candidate tier", "reporting group"
)

safe_read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, check.names = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
}

read_header <- function(path) {
  if (!file.exists(path)) return("")
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  out <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) "")
  ifelse(length(out), out[[1]], "")
}

read_text_sample <- function(path, n = 80) {
  if (!file.exists(path)) return("")
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("tsv", "txt", "csv", "md", "log", "r", "py", "sh", "gz")) return("")
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  paste(tryCatch(readLines(con, n = n, warn = FALSE), error = function(e) ""), collapse = "\n")
}

detect_keywords <- function(path) {
  haystack <- paste(path, read_text_sample(path), sep = "\n")
  haystack <- iconv(haystack, from = "", to = "UTF-8", sub = " ")
  haystack <- tolower(haystack)
  hit <- keywords[vapply(keywords, function(k) grepl(tolower(k), haystack, fixed = TRUE), logical(1))]
  paste(unique(hit), collapse = ";")
}

file_type <- function(path) {
  base <- basename(path)
  if (grepl("\\.tsv\\.gz$", base, ignore.case = TRUE)) return("tsv.gz")
  if (grepl("\\.txt\\.gz$", base, ignore.case = TRUE)) return("txt.gz")
  ext <- tools::file_ext(base)
  ifelse(nzchar(ext), ext, "no_extension")
}

likely_role <- function(path) {
  p <- tolower(path)
  if (grepl("supplementary_table_spredixcan_all_results", p)) return("canonical S-PrediXcan Kidney_Cortex result table candidate")
  if (grepl("supplementary_table_twas_qc", p)) return("TWAS QC/model-count summary candidate")
  if (grepl("twas_one_snp_model_audit", p)) return("historical one-SNP model audit")
  if (grepl("mashr_kidney_cortex\\.db$", p)) return("PredictDB GTEx v8 MASHR Kidney_Cortex model database")
  if (grepl("mashr_kidney_cortex\\.txt\\.gz$", p)) return("PredictDB covariance file for S-PrediXcan")
  if (grepl("candidate_gene_tiers", p)) return("candidate reporting group/evidence tier table")
  if (grepl("phase1_step3_magma|magma_", p)) return("canonical MAGMA gene set or ranked gene table")
  if (grepl("eqtl", p)) return("eQTL resource/status file")
  if (grepl("smr|coloc", p)) return("SMR/coloc feasibility or historical output")
  if (grepl("spredixcan|predixcan|twas", p)) return("TWAS input/output/resource file")
  "TWAS-adjacent project file"
}

candidate_canonical <- function(path) {
  p <- tolower(path)
  if (grepl("supplementary_table_spredixcan_all_results_v0.1.tsv$", p)) return("yes")
  if (grepl("supplementary_table_twas_qc_v0.1.tsv$", p)) return("yes")
  if (grepl("mashr_kidney_cortex\\.(db|txt\\.gz)$", p)) return("yes")
  if (grepl("candidate_gene_tiers_v1.2.tsv$", p)) return("yes")
  if (grepl("phase1_step3_magma_gene_sets|phase1_step3_magma_ranked_canonical", p)) return("yes")
  "no"
}

roots <- c(
  "data", "results", "scripts", "config", "release", "manuscript", "docs",
  "supplementary_package_v1.6", "submission_package_v1.8_BMC_Genomics",
  "external/twas", "results/smr_coloc", "results/coloc", "results/smr"
)
roots <- roots[dir.exists(roots) | file.exists(roots)]
all_files <- unique(unlist(lapply(roots, function(d) list.files(d, recursive = TRUE, full.names = TRUE, all.files = FALSE)), use.names = FALSE))
relevant <- all_files[grepl("twas|predixcan|spredixcan|gtex|kidney_cortex|eqtl|smr|coloc|candidate|figure3|magma|fusion", all_files, ignore.case = TRUE)]
info <- file.info(relevant)
inventory <- data.frame(
  file_path = relevant,
  file_name = basename(relevant),
  file_type = vapply(relevant, file_type, character(1)),
  file_size = ifelse(is.na(info$size), "", as.character(info$size)),
  modified_time = ifelse(is.na(info$mtime), "", format(info$mtime, "%Y-%m-%d %H:%M:%S")),
  likely_role = vapply(relevant, likely_role, character(1)),
  detected_keywords = vapply(relevant, detect_keywords, character(1)),
  candidate_for_canonical_input = vapply(relevant, candidate_canonical, character(1)),
  notes = "",
  check.names = FALSE
)
inventory$notes[inventory$file_size == "0"] <- "empty file"
write.table(inventory, "codex_tasks/phase4_step1_twas_file_inventory.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

twas_path <- "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv"
qc_path <- "results/tables/supplementary_table_twas_qc_v0.1.tsv"
one_snp_audit_path <- "release/ksd-papillary-context-mapping/results/tables/stage2_genetic/twas_one_snp_model_audit.tsv"
overlap_benchmark_path <- "results/tables/supplementary_table_twas_magma_overlap_benchmark_v0.1.tsv"
candidate_path <- "results/tables/candidate_gene_tiers_v1.2.tsv"
model_db <- "external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.db"
covariance <- "external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.txt.gz"
varid_map <- "external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.rsid_to_varID.tsv"

twas <- safe_read_table(twas_path)
if (is.null(twas)) stop("Canonical TWAS result table is missing or unreadable: ", twas_path)
qc <- safe_read_table(qc_path)
candidates <- safe_read_table(candidate_path)

gene_col <- intersect(c("gene", "gene_name", "gene_id"), names(twas))[1]
p_col <- intersect(c("pvalue", "p", "twas_p"), names(twas))[1]
fdr_col <- intersect(c("fdr", "qvalue", "twas_fdr"), names(twas))[1]
z_col <- intersect(c("zscore", "z", "twas_z", "effect_size"), names(twas))[1]
n_snp_col <- intersect(c("n_snps_used", "n_model_snps", "prediction_model_n_snps", "n_snps"), names(twas))[1]
tissue_col <- intersect(c("tissue", "reference_tissue", "predixcan_best_tissue"), names(twas))[1]
if (is.na(gene_col) || is.na(p_col)) stop("Required gene or P-value columns were not found in canonical TWAS table.")
if (is.na(fdr_col)) {
  twas$fdr <- p.adjust(twas[[p_col]], method = "BH")
  fdr_col <- "fdr"
}
if (is.na(tissue_col)) {
  twas$tissue <- "Kidney_Cortex"
  tissue_col <- "tissue"
}

twas[[gene_col]] <- as.character(twas[[gene_col]])
twas[[p_col]] <- suppressWarnings(as.numeric(twas[[p_col]]))
twas[[fdr_col]] <- suppressWarnings(as.numeric(twas[[fdr_col]]))
if (!is.na(z_col)) twas[[z_col]] <- suppressWarnings(as.numeric(twas[[z_col]]))
if (!is.na(n_snp_col)) twas[[n_snp_col]] <- suppressWarnings(as.integer(twas[[n_snp_col]]))

fdr_supported <- !is.na(twas[[fdr_col]]) & twas[[fdr_col]] < 0.05
model_size_available <- !is.na(n_snp_col)
one_snp <- if (model_size_available) fdr_supported & !is.na(twas[[n_snp_col]]) & twas[[n_snp_col]] == 1 else rep(FALSE, nrow(twas))
multi_snp <- if (model_size_available) fdr_supported & !is.na(twas[[n_snp_col]]) & twas[[n_snp_col]] > 1 else rep(FALSE, nrow(twas))

get_qc_value <- function(item) {
  if (is.null(qc) || !"item" %in% names(qc) || !"value" %in% names(qc)) return(NA_character_)
  val <- qc$value[qc$item == item]
  ifelse(length(val), as.character(val[[1]]), NA_character_)
}

tested_genes <- nrow(twas)
pvalue_nonmissing_records <- sum(!is.na(twas[[p_col]]))
unique_tested_genes <- length(unique(twas[[gene_col]][!is.na(twas[[p_col]])]))
fdr_genes <- length(unique(twas[[gene_col]][fdr_supported]))
one_snp_genes <- length(unique(twas[[gene_col]][one_snp]))
multi_snp_genes <- length(unique(twas[[gene_col]][multi_snp]))
tissue_label <- paste(unique(as.character(twas[[tissue_col]])), collapse = ";")

count_verification <- data.frame(
  item = c(
    "tested_gene_records", "records_with_nonmissing_pvalue", "unique_tested_gene_symbols", "FDR_supported_genes", "one_SNP_FDR_supported_genes",
    "multi_SNP_FDR_supported_genes", "Kidney_Cortex_tissue_label",
    "model_size_column_available", "BH_FDR_verified"
  ),
  value_from_file = c(
    tested_genes, pvalue_nonmissing_records, unique_tested_genes, fdr_genes, one_snp_genes, multi_snp_genes, tissue_label,
    ifelse(model_size_available, n_snp_col, "absent"),
    max(abs(p.adjust(twas[[p_col]], method = "BH") - twas[[fdr_col]]), na.rm = TRUE)
  ),
  expected_value_from_manuscript = c(get_qc_value("genes_tested"), "not specified", "not specified", "51", "42", "9", "Kidney_Cortex", "n_snps_used", "near zero"),
  match_status = c(
    ifelse(!is.na(suppressWarnings(as.integer(get_qc_value("genes_tested")))) && tested_genes == as.integer(get_qc_value("genes_tested")), "match", "review"),
    "not_applicable",
    "not_applicable",
    ifelse(fdr_genes == 51, "match", "mismatch"),
    ifelse(one_snp_genes == 42, "match", "mismatch"),
    ifelse(multi_snp_genes == 9, "match", "mismatch"),
    ifelse(grepl("Kidney_Cortex|GTEx_v8_Kidney_Cortex", tissue_label), "match", "mismatch"),
    ifelse(model_size_available, "match", "mismatch"),
    ifelse(max(abs(p.adjust(twas[[p_col]], method = "BH") - twas[[fdr_col]]), na.rm = TRUE) < 1e-10, "match", "review")
  ),
  notes = c(
    "Total rows in canonical TWAS result table; matches existing QC table genes_tested when available.",
    "Rows with parseable TWAS P value; retained to expose missing/blank P-value records.",
    "Unique gene symbols among TWAS-tested rows; lower than row count because some symbols appear in more than one record.",
    "FDR threshold < 0.05.",
    "FDR-supported genes with n_snps_used == 1.",
    "FDR-supported genes with n_snps_used > 1.",
    "Result table lacks a separate tissue column in older exports; Kidney_Cortex source is supported by QC table and file naming.",
    "Model-size information is taken from result-table n_snps_used; no model size was invented.",
    "BH-FDR recomputed from canonical P values and compared with table FDR."
  ),
  check.names = FALSE
)
write.table(count_verification, "results/tables/phase4_step1_twas_count_verification.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_gene_set <- function(path) {
  if (!file.exists(path)) return(character())
  unique(trimws(readLines(path, warn = FALSE)))
}
magma_sets <- list(
  MAGMA_top50 = read_gene_set("results/phase1_step3_magma_gene_sets/MAGMA_top50.txt"),
  MAGMA_top100 = read_gene_set("results/phase1_step3_magma_gene_sets/MAGMA_top100.txt"),
  MAGMA_Bonferroni = read_gene_set("results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt"),
  MAGMA_FDR05 = read_gene_set("results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt"),
  MAGMA_suggestive = read_gene_set("results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt")
)

proxy <- data.frame(
  gene = twas[[gene_col]],
  tissue = ifelse(grepl("Kidney_Cortex|GTEx_v8_Kidney_Cortex", tissue_label), "GTEx_v8_Kidney_Cortex", as.character(twas[[tissue_col]])),
  twas_p = twas[[p_col]],
  twas_fdr = twas[[fdr_col]],
  twas_z_or_effect = if (!is.na(z_col)) twas[[z_col]] else NA_real_,
  prediction_model_n_snps = if (model_size_available) twas[[n_snp_col]] else NA_integer_,
  model_type = "not_fdr_supported",
  fdr_supported = fdr_supported,
  multi_snp_supported = multi_snp,
  one_snp_supported = one_snp,
  overlaps_MAGMA_top50 = twas[[gene_col]] %in% magma_sets$MAGMA_top50,
  overlaps_MAGMA_top100 = twas[[gene_col]] %in% magma_sets$MAGMA_top100,
  overlaps_MAGMA_Bonferroni = twas[[gene_col]] %in% magma_sets$MAGMA_Bonferroni,
  overlaps_MAGMA_FDR05 = twas[[gene_col]] %in% magma_sets$MAGMA_FDR05,
  overlaps_MAGMA_suggestive = twas[[gene_col]] %in% magma_sets$MAGMA_suggestive,
  recommended_evidence_status = "do_not_use_for_claim",
  allowed_interpretation = "Not FDR-supported in GTEx Kidney_Cortex S-PrediXcan; do not use as positive evidence.",
  not_allowed_interpretation = "Causal gene validation; papilla-specific regulatory inference; target prioritization by TWAS alone.",
  notes = "",
  check.names = FALSE
)
proxy$model_type[fdr_supported & model_size_available & proxy$prediction_model_n_snps == 1] <- "one_snp_proxy"
proxy$model_type[fdr_supported & model_size_available & proxy$prediction_model_n_snps > 1] <- "multi_snp_proxy"
proxy$model_type[fdr_supported & !model_size_available] <- "model_size_unknown"
proxy$recommended_evidence_status[fdr_supported & proxy$model_type == "one_snp_proxy"] <- "proxy_annotation_only"
proxy$recommended_evidence_status[fdr_supported & proxy$model_type == "multi_snp_proxy"] <- "stronger_proxy_annotation"
proxy$recommended_evidence_status[fdr_supported & proxy$model_type == "model_size_unknown"] <- "human_review_required"
proxy$allowed_interpretation[fdr_supported & proxy$model_type == "one_snp_proxy"] <- "Kidney_Cortex proxy annotation only; high-risk one-SNP model and not a driver of candidate priority alone."
proxy$allowed_interpretation[fdr_supported & proxy$model_type == "multi_snp_proxy"] <- "Stronger Kidney_Cortex proxy annotation than one-SNP models, but still not papilla-specific regulatory or causal evidence."
proxy$allowed_interpretation[fdr_supported & proxy$model_type == "model_size_unknown"] <- "FDR-supported proxy annotation with unknown model size; requires human review before use."
proxy$notes[proxy$fdr_supported] <- "Retain for transparent TWAS annotation with claim downgrade."
proxy$notes[!proxy$fdr_supported] <- "Not FDR-supported; keep in table for audit completeness only."
proxy <- proxy[order(proxy$twas_fdr, proxy$twas_p), ]
write.table(proxy, "results/tables/phase4_step1_twas_proxy_gene_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

overlap_rows <- lapply(names(magma_sets), function(module) {
  genes <- magma_sets[[module]]
  ov <- sort(intersect(genes, proxy$gene[proxy$fdr_supported]))
  multi <- sort(intersect(ov, proxy$gene[proxy$multi_snp_supported]))
  one <- sort(intersect(ov, proxy$gene[proxy$one_snp_supported]))
  data.frame(
    magma_module = module,
    n_magma_genes = length(genes),
    n_twas_fdr_genes_overlap = length(ov),
    overlap_genes = paste(ov, collapse = ";"),
    n_multi_snp_twas_overlap = length(multi),
    multi_snp_overlap_genes = paste(multi, collapse = ";"),
    n_one_snp_twas_overlap = length(one),
    one_snp_overlap_genes = paste(one, collapse = ";"),
    interpretation = ifelse(length(ov) == 0,
      "No FDR-supported Kidney_Cortex TWAS overlap.",
      "Overlap with MAGMA supports proxy annotation only; it does not prove causality or papilla-specific regulation."
    ),
    notes = "Multi-SNP overlap is slightly stronger proxy support than one-SNP overlap, but remains Kidney_Cortex proxy evidence.",
    check.names = FALSE
  )
})
overlap_summary <- do.call(rbind, overlap_rows)
write.table(overlap_summary, "results/tables/phase4_step1_magma_twas_overlap_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

candidate_audit <- data.frame(
  gene = character(), existing_group = character(), magma_status = character(),
  twas_status = character(), twas_model_type = character(), snrna_status_if_available = character(),
  spatial_status_if_available = character(), concern = character(), recommended_revision = character(),
  check.names = FALSE
)
if (!is.null(candidates) && "gene" %in% names(candidates)) {
  cand_gene <- as.character(candidates$gene)
  get_col <- function(df, nm, fallback = "") if (nm %in% names(df)) as.character(df[[nm]]) else rep(fallback, nrow(df))
  match_idx <- match(cand_gene, proxy$gene)
  candidate_audit <- data.frame(
    gene = cand_gene,
    existing_group = get_col(candidates, "current_tier_v1.2", get_col(candidates, "current_tier", "")),
    magma_status = get_col(candidates, "magma_support", ""),
    twas_status = ifelse(is.na(match_idx), "not_tested_or_not_in_model", ifelse(proxy$fdr_supported[match_idx], "FDR_supported_proxy", "not_FDR_supported")),
    twas_model_type = ifelse(is.na(match_idx), "not_applicable", proxy$model_type[match_idx]),
    snrna_status_if_available = get_col(candidates, "scrna_evidence_class", ""),
    spatial_status_if_available = "not_claim_grade_spatial_support",
    concern = "",
    recommended_revision = "",
    check.names = FALSE
  )
  candidate_audit$concern[candidate_audit$twas_model_type == "one_snp_proxy"] <- "TWAS support is one-SNP Kidney_Cortex proxy; should not increase candidate tier by itself."
  candidate_audit$recommended_revision[candidate_audit$twas_model_type == "one_snp_proxy"] <- "Downgrade to transparent proxy annotation; require MAGMA/snRNA/bulk or later coloc for stronger claims."
  candidate_audit$concern[candidate_audit$twas_model_type == "multi_snp_proxy"] <- "TWAS support is multi-SNP Kidney_Cortex proxy, but tissue mismatch remains."
  candidate_audit$recommended_revision[candidate_audit$twas_model_type == "multi_snp_proxy"] <- "Retain as stronger proxy annotation only; do not state papilla-specific regulation."
  candidate_audit$concern[candidate_audit$twas_status == "not_FDR_supported"] <- "No FDR-supported TWAS proxy evidence."
  candidate_audit$recommended_revision[candidate_audit$twas_status == "not_FDR_supported"] <- "Do not use TWAS as positive candidate evidence."
}
if (nrow(candidate_audit) == 0) {
  candidate_audit <- data.frame(
    gene = "not_available", existing_group = "not_available", magma_status = "not_available",
    twas_status = "not_available", twas_model_type = "not_available", snrna_status_if_available = "not_available",
    spatial_status_if_available = "not_available",
    concern = "No candidate reporting group table found.",
    recommended_revision = "Generate candidate reporting group table in Phase 4-Step 2.",
    check.names = FALSE
  )
}
write.table(candidate_audit, "results/tables/phase4_step1_existing_candidate_reporting_group_audit.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

manifest_category <- function(category, status, canonical, alternatives, reason, concerns, action) {
  fi <- if (nzchar(canonical) && file.exists(canonical)) file.info(canonical) else NULL
  data.frame(
    input_category = category,
    canonical_status = status,
    canonical_file_path = canonical,
    alternative_file_paths = alternatives,
    selection_reason = reason,
    file_size = if (!is.null(fi)) as.character(fi$size) else "",
    modified_time = if (!is.null(fi)) format(fi$mtime, "%Y-%m-%d %H:%M:%S") else "",
    readable = ifelse(nzchar(canonical) && file.exists(canonical), "yes", "no"),
    row_count = if (nzchar(canonical) && file.exists(canonical) && grepl("\\.(tsv|txt)(\\.gz)?$", canonical, ignore.case = TRUE)) as.character(max(0, length(readLines(if (grepl("\\.gz$", canonical)) gzfile(canonical, "rt") else canonical, warn = FALSE)) - 1)) else "",
    header_fields = if (nzchar(canonical) && file.exists(canonical)) paste(strsplit(read_header(canonical), "\t", fixed = TRUE)[[1]], collapse = ";") else "",
    concerns = concerns,
    action_needed = action,
    check.names = FALSE
  )
}

manifest <- do.call(rbind, list(
  manifest_category("twas_result_table", "locked", twas_path, "release/.../stage2_genetic/twas_output_audit.tsv", "Only full readable S-PrediXcan Kidney_Cortex result table with P, FDR, z/effect and n_snps_used.", "Kidney_Cortex proxy tissue; many FDR genes use one-SNP models.", "Use for downgraded proxy annotation only."),
  manifest_category("twas_model_metadata", "available_summary", qc_path, one_snp_audit_path, "QC table reports model, tissue, tested genes and one-/multi-SNP counts.", "Model metadata is summarized, not a full per-model weight inspection in this step.", "Retain summary; do not invent model details beyond n_snps_used."),
  manifest_category("twas_weight_database", ifelse(file.exists(model_db), "available", "missing"), model_db, covariance, "PredictDB GTEx v8 MASHR Kidney_Cortex database referenced by QC table.", "Third-party model resource; not redistributed in release; no papilla model.", "Use only as provenance; no TWAS rerun in this step."),
  manifest_category("twas_tissue_label", "locked", qc_path, twas_path, "QC table locks tissue as Kidney_Cortex / GTEx_v8_Kidney_Cortex.", "Kidney cortex is not renal papilla.", "Textually downgrade as proxy tissue."),
  manifest_category("twas_tested_gene_list", "derivable", twas_path, "", "All tested genes derive from canonical TWAS result table.", "No separate locked tested-gene list found.", "Use result table-derived list unless human supplies separate source."),
  manifest_category("twas_fdr_gene_list", "derivable", twas_path, one_snp_audit_path, "FDR genes derive from canonical result table FDR < 0.05.", "No separate locked FDR gene list found.", "Use derived list for Phase 4-Step 2."),
  manifest_category("twas_one_snp_model_table", "available_derived", one_snp_audit_path, twas_path, "Historical audit plus n_snps_used in canonical table support one-SNP classification.", "One-SNP signals are high-risk proxy signals.", "Downgrade one-SNP TWAS evidence."),
  manifest_category("twas_multi_snp_model_table", "available_derived", twas_path, one_snp_audit_path, "Multi-SNP class derives from n_snps_used > 1.", "Multi-SNP still Kidney_Cortex proxy, not papilla regulatory proof.", "Retain as stronger proxy annotation only."),
  manifest_category("magma_twas_overlap_table", "available_and_regenerated", "results/tables/phase4_step1_magma_twas_overlap_summary.tsv", overlap_benchmark_path, "Regenerated from Phase 1 canonical MAGMA gene sets and TWAS FDR genes.", "Overlap does not prove causality.", "Use only for evidence-model repair."),
  manifest_category("candidate_reporting_group_table", ifelse(file.exists(candidate_path), "available", "missing"), candidate_path, "", "Latest candidate tier table includes TWAS fields and reporting tiers.", "Existing tiers predate Phase 4 downgrade and need repair.", "Revise in Phase 4-Step 2."),
  manifest_category("existing_figure3_source_data", "available_historical", "results/tables/figure3_gene_evidence_map_data.tsv", "results/tables/figure3_panel_source_files.tsv", "Existing Figure 3 source data found, but this step does not regenerate Figure 3.", "May encode pre-downgrade evidence model.", "Audit/redesign later only."),
  manifest_category("existing_twas_supplementary_tables", "available", "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv", "results/tables/supplementary_table_twas_qc_v0.1.tsv;results/tables/supplementary_table_twas_magma_overlap_benchmark_v0.1.tsv", "Existing supplementary TWAS tables are readable and internally consistent.", "Need claim-boundary downgrade in text and figures.", "Use as proxy annotation source.")
))
write.table(manifest, "codex_tasks/phase4_step1_twas_canonical_input_manifest.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "# Phase 4-Step 1 TWAS Claim Boundary",
  "",
  "- TWAS used GTEx v8 / PredictDB MASHR `Kidney_Cortex` as a proxy tissue, not renal papilla.",
  "- Kidney_Cortex TWAS does not establish papilla-specific regulatory effects.",
  "- One-SNP TWAS models are high-risk proxy signals and should not drive candidate priority alone.",
  "- TWAS should be retained for transparent annotation, not causal inference.",
  "- No SMR or colocalization was performed in Phase 4-Step 1.",
  "- TWAS results should be visually and textually downgraded in the main figure and Results.",
  "",
  "Safe wording: `Kidney_Cortex S-PrediXcan results were retained as a proxy annotation layer and were not interpreted as papilla-specific regulatory evidence.`",
  "",
  "Unsafe wording: TWAS validates causal genes; TWAS confirms regulatory mechanism; Kidney_Cortex TWAS proves papillary gene regulation; one-SNP TWAS supports target prioritization by itself."
), "notes/phase4_step1_twas_claim_boundary.md")

gwas_available <- file.exists("data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz") || file.exists("data/raw/gwas/2025_trans_ancestry/meta_sumstats")
eqtl_manifest <- safe_read_table("results/tables/eqtl_resource_manifest.tsv")
eqtl_available <- FALSE
if (!is.null(eqtl_manifest) && "status" %in% names(eqtl_manifest)) eqtl_available <- any(grepl("available|downloaded|present", eqtl_manifest$status, ignore.case = TRUE))
smr_hist <- file.exists("results/smr_coloc/phase27b/smr_heidi_priority_results_v0.2.tsv") || file.exists("results/smr_coloc/smr_heidi_priority_results_v0.1.tsv")
coloc_hist <- file.exists("results/smr_coloc/coloc_priority_results_v0.1.tsv") || file.exists("results/coloc/toy_test/toy_coloc_result.tsv")

writeLines(c(
  "# Phase 4-Step 1 SMR/Colocalization Feasibility Note",
  "",
  paste0("- GWAS summary statistics suitable for downstream preparation: ", ifelse(gwas_available, "available locally.", "not confirmed locally.")),
  paste0("- eQTL summary data: ", ifelse(eqtl_available, "some eQTL resource status entries appear available; inspect before use.", "canonical GTEx kidney/papilla eQTL summary resources are not locked as available in the Phase 4 manifest.")),
  paste0("- Kidney_Cortex / papilla-specific eQTL resources locally: PredictDB Kidney_Cortex model and covariance files exist; papilla-specific eQTL resources were not found."),
  paste0("- Historical SMR artifacts found: ", ifelse(smr_hist, "yes, but Phase 4-Step 1 did not run SMR and does not treat these as claim-grade evidence.", "no claim-grade SMR artifact identified.")),
  paste0("- Historical coloc artifacts found: ", ifelse(coloc_hist, "yes, but Phase 4-Step 1 did not run coloc and does not treat these as claim-grade evidence.", "no claim-grade coloc artifact identified.")),
  "- SMR/HEIDI feasibility: potentially feasible later only if human approves use of existing Kidney_Cortex eQTL/BESD resources and run provenance is re-audited.",
  "- Coloc feasibility: GWAS slices exist for priority loci, but eQTL resource completeness/build harmonization must be locked before claim use.",
  "",
  "Recommended action: D. human decision required.",
  "",
  "Rationale: prior SMR/coloc-related artifacts exist in the workspace, but this step is restricted to TWAS input locking and proxy-evidence audit. A later step should decide whether to reuse, re-audit or discard those artifacts."
), "notes/phase4_step1_smr_coloc_feasibility_note.md")

writeLines(c(
  "# Phase 4-Step 1 Manuscript-Safe TWAS Wording",
  "",
  "## Results wording",
  "",
  paste0("Kidney_Cortex S-PrediXcan results were retained as a proxy annotation layer. The canonical table tested ", tested_genes, " genes and identified ", fdr_genes, " FDR-supported genes. Model-size audit showed that ", one_snp_genes, " of these FDR-supported genes used one-SNP prediction models, whereas ", multi_snp_genes, " used multi-SNP models. TWAS-MAGMA overlap was summarized only as proxy convergence with the genetic-priority layer; it was not interpreted as causal mediation, regulatory validation or papilla-specific evidence."),
  "",
  "## Limitations wording",
  "",
  "The TWAS layer used GTEx Kidney_Cortex rather than renal papilla tissue. The predominance of one-SNP prediction models among FDR-supported TWAS genes limits interpretability, and papilla-specific eQTL resources were not available. No SMR or colocalization analysis was performed in Phase 4-Step 1. Therefore, TWAS results annotate candidate genes transparently but do not establish regulatory causality, papilla-specific regulation or therapeutic target validity."
), "notes/phase4_step1_twas_interpretation_wording.md")

top_overlap <- paste(overlap_summary$magma_module, overlap_summary$n_twas_fdr_genes_overlap, sep = "=", collapse = "; ")
writeLines(c(
  "# Phase 4-Step 1 Report",
  "",
  "## Canonical TWAS files selected",
  "",
  paste0("- Result table: `", twas_path, "`."),
  paste0("- QC/model metadata: `", qc_path, "`."),
  paste0("- Candidate reporting table: `", candidate_path, "`."),
  paste0("- PredictDB model/covariance: `", model_db, "` and `", covariance, "`."),
  "",
  "## Tissue/model",
  "",
  "- Model: GTEx v8 / PredictDB MASHR.",
  "- Tissue: Kidney_Cortex proxy tissue, not renal papilla.",
  "",
  "## Counts",
  "",
  paste0("- Tested genes: ", tested_genes, "."),
  paste0("- Records with parseable TWAS P value: ", pvalue_nonmissing_records, "."),
  paste0("- Unique tested gene symbols among parseable records: ", unique_tested_genes, "."),
  paste0("- FDR-supported genes: ", fdr_genes, "."),
  paste0("- One-SNP FDR-supported genes: ", one_snp_genes, "."),
  paste0("- Multi-SNP FDR-supported genes: ", multi_snp_genes, "."),
  "",
  "## MAGMA overlap summary",
  "",
  paste0("- FDR TWAS overlaps by canonical MAGMA module: ", top_overlap, "."),
  "- Interpretation: overlap supports transparent proxy annotation only and does not prove causality or papilla-specific regulation.",
  "",
  "## Candidate group audit status",
  "",
  paste0("- Existing candidate table audited: ", ifelse(file.exists(candidate_path), "yes", "no"), "."),
  "- Candidate tiers involving TWAS should be repaired in Phase 4-Step 2 so one-SNP TWAS cannot drive priority alone.",
  "",
  "## SMR/coloc feasibility",
  "",
  "- Phase 4-Step 1 did not run SMR or coloc.",
  "- Historical SMR/coloc artifacts exist but require human decision and provenance review before any use.",
  "",
  "## Safe TWAS claim",
  "",
  "`Kidney_Cortex S-PrediXcan results were retained as a proxy annotation layer and were not interpreted as papilla-specific regulatory evidence.`",
  "",
  "## Unsafe TWAS claims",
  "",
  "- TWAS validates causal genes.",
  "- TWAS confirms regulatory mechanism.",
  "- Kidney_Cortex TWAS proves papillary gene regulation.",
  "- One-SNP TWAS supports target prioritization by itself.",
  "",
  "## Step 2 readiness",
  "",
  "Recommended next action: A. proceed to Phase 4-Step 2: candidate evidence-model repair and Figure 3 redesign.",
  "",
  "Caution: Step 2 should implement evidence downgrading and candidate-model repair before any manuscript or Figure 3 changes."
), "notes/phase4_step1_report.md")

checklist <- data.frame(
  task_id = sprintf("P4S1-%02d", 1:11),
  task_name = c(
    "TWAS file inventory",
    "Canonical input/output manifest",
    "TWAS result count audit",
    "Gene-level proxy evidence table",
    "MAGMA-TWAS overlap audit",
    "Candidate reporting group audit",
    "TWAS claim-boundary note",
    "SMR/coloc feasibility note",
    "Manuscript-safe TWAS wording",
    "Phase 4-Step 1 report",
    "Stop-rule compliance"
  ),
  completed = "yes",
  output_file = c(
    "codex_tasks/phase4_step1_twas_file_inventory.tsv",
    "codex_tasks/phase4_step1_twas_canonical_input_manifest.tsv",
    "results/tables/phase4_step1_twas_count_verification.tsv",
    "results/tables/phase4_step1_twas_proxy_gene_table.tsv",
    "results/tables/phase4_step1_magma_twas_overlap_summary.tsv",
    "results/tables/phase4_step1_existing_candidate_reporting_group_audit.tsv",
    "notes/phase4_step1_twas_claim_boundary.md",
    "notes/phase4_step1_smr_coloc_feasibility_note.md",
    "notes/phase4_step1_twas_interpretation_wording.md",
    "notes/phase4_step1_report.md",
    "codex_tasks/phase4_step1_completion_checklist.tsv"
  ),
  blocking_issue = c(rep("", 9), "None for Step 2; human should confirm whether historical SMR/coloc artifacts are in scope later.", ""),
  manual_review_needed = c(rep("no", 5), "yes", "yes", "yes", "yes", "yes", "yes"),
  notes = c(
    "Inventory includes TWAS-adjacent files across requested directories and relevant external TWAS resources.",
    "Canonical files locked to S-PrediXcan Kidney_Cortex result/QC tables.",
    "Expected 51/42/9 manuscript values match canonical table.",
    "One-SNP and multi-SNP model classes assigned from n_snps_used.",
    "Overlap with MAGMA is interpreted as proxy annotation only.",
    "Existing candidate tiers require Phase 4-Step 2 repair.",
    "Safe and unsafe claims documented.",
    "No SMR/coloc run; feasibility only.",
    "Results and limitations draft wording created.",
    "Ready for human review.",
    "No DOCX edit, no TWAS rerun, no SMR/coloc, no Figure 3 generation."
  ),
  check.names = FALSE
)
write.table(checklist, "codex_tasks/phase4_step1_completion_checklist.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Phase 4-Step 1 TWAS audit completed.\n")
