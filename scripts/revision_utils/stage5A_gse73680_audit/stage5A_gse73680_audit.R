suppressPackageStartupMessages({
  library(data.table)
})

doc_dir <- "docs/revision/stage5A_gse73680_audit"
table_dir <- "results/tables/revision/stage5A_gse73680_audit"
fig_dir <- "results/figures/revision/stage5A_gse73680_audit"
script_dir <- "scripts/revision_utils/stage5A_gse73680_audit"
log_dir <- "logs/revision/stage5A_gse73680_audit"
for (d in c(doc_dir, table_dir, fig_dir, script_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage5A_gse73680_audit.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 5A GSE73680 resource audit and statistical design\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  stage3_model = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  stage3_counts = "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv",
  exemplar = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
  stage4_claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
  stage4c2r_report = "docs/revision/stage4C2R_draft_figures_v0.2/stage4C2R_report.md",
  gene_expr = "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz",
  feature_expr = "data/processed/gse73680/gse73680_expression_matrix.tsv.gz",
  patient_expr = "data/processed/gse73680/gse73680_patient_level_expression_matrix.tsv.gz",
  metadata = "config/gse73680_sample_metadata_curated.tsv",
  metadata_raw = "config/gse73680_sample_metadata.tsv",
  feature_mapping = "results/gse73680/tables/gse73680_feature_gene_mapping.tsv",
  design = "results/gse73680/tables/gse73680_analysis_design.tsv",
  expression_qc = "results/gse73680/tables/gse73680_expression_qc_v2.tsv",
  module_manifest = "results/tables/revision/stage4B1_scrna_donor_level/stage4B1_gene_module_manifest.tsv",
  gene_availability = "results/gse73680/tables/gse73680_gene_availability.tsv",
  magma_top50 = "results/gene_sets/magma_top50.txt",
  magma_top100 = "results/gene_sets/magma_top100.txt"
)

for (p in paths[1:5]) {
  if (!file.exists(p)) stop("Required upstream revision input missing: ", p)
}

readable_file <- function(p) file.exists(p) && file.access(p, 4) == 0

first_existing <- function(candidates) {
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) hits[1] else NA_character_
}

classify_file <- function(path) {
  low <- tolower(path)
  type <- fifelse(grepl("\\.(r|sh)$", low), "script",
           fifelse(grepl("\\.(tsv|txt|csv)(\\.gz)?$", low), "table_or_text",
           fifelse(grepl("\\.(pdf|png|svg)$", low), "figure",
           fifelse(grepl("\\.(tar|gz)$", low), "raw_or_compressed", "other"))))
  content <- "GSE73680-related file"
  if (grepl("gene_expression|expression_matrix|normalized", low)) content <- "expression matrix"
  if (grepl("metadata|sample_sheet|sample_metadata", low)) content <- "sample metadata"
  if (grepl("feature|probe|mapping|annotation|platform", low)) content <- "feature/probe-to-gene annotation"
  if (grepl("paired|patient_sample|design", low)) content <- "paired sample/design table"
  if (grepl("limma|response|differential|exact_statistics", low)) content <- "existing differential/module response result"
  if (grepl("module|score|camera|roast|gsea|random|benchmark", low)) content <- "existing module or benchmark result"
  if (grepl("composition|signature|injury|stromal|immune|endothelial", low)) content <- "existing context/composition-related result"
  if (grepl("figure4|figures|source", low)) content <- "existing figure/source artifact"
  list(file_type = type, likely_content = content)
}

all_files <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = TRUE)
stage5a_output_patterns <- c(
  "^\\./docs/revision/stage5A_gse73680_audit/",
  "^\\./results/tables/revision/stage5A_gse73680_audit/",
  "^\\./results/figures/revision/stage5A_gse73680_audit/",
  "^\\./scripts/revision_utils/stage5A_gse73680_audit/",
  "^\\./logs/revision/stage5A_gse73680_audit/"
)
for (pat in stage5a_output_patterns) all_files <- all_files[!grepl(pat, all_files)]
gse_files <- all_files[grepl("gse73680|73680|plaque|papilla", all_files, ignore.case = TRUE)]
file_inventory <- rbindlist(lapply(gse_files, function(p) {
  info <- file.info(p)
  cls <- classify_file(p)
  data.table(
    file_path = sub("^\\./", "", p),
    file_name = basename(p),
    file_type = cls$file_type,
    size_bytes = as.numeric(info$size),
    last_modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    likely_content = cls$likely_content,
    readable = readable_file(p),
    ready_for_stage5B = fifelse(readable_file(p) & cls$file_type %in% c("table_or_text", "raw_or_compressed"), "candidate", "supporting_or_not_primary"),
    notes = "Inventory only; readiness for primary analysis determined in dedicated audits."
  )
}), fill = TRUE)
setorder(file_inventory, likely_content, file_path)
fwrite(file_inventory, file.path(table_dir, "gse73680_file_inventory.tsv"), sep = "\t")

primary_expr <- paths$gene_expr
primary_meta <- paths$metadata
primary_mapping <- paths$feature_mapping

if (!readable_file(primary_expr) || !readable_file(primary_meta)) {
  writeLines(c(
    "# GSE73680 resource blocker",
    "",
    paste0("Generated: ", Sys.Date()),
    "",
    "- No usable primary gene expression matrix or curated metadata table was found.",
    paste0("- Candidate expression matrix: `", primary_expr, "` readable = ", readable_file(primary_expr)),
    paste0("- Candidate metadata table: `", primary_meta, "` readable = ", readable_file(primary_meta)),
    "",
    "Stage 5B should not start until these resources are restored."
  ), file.path(doc_dir, "gse73680_resource_blocker.md"))
  stop("Primary GSE73680 expression matrix or metadata missing.")
}

expr <- fread(primary_expr)
feature_expr <- if (readable_file(paths$feature_expr)) fread(paths$feature_expr, nThread = 1) else NULL
meta <- fread(primary_meta)
mapping <- if (readable_file(primary_mapping)) fread(primary_mapping) else data.table()
stage3 <- fread(paths$stage3_model)
exemplar <- fread(paths$exemplar)
module_manifest <- if (readable_file(paths$module_manifest)) fread(paths$module_manifest) else data.table()

id_col <- names(expr)[1]
sample_cols <- setdiff(names(expr), id_col)
expr_genes <- unique(na.omit(as.character(expr[[id_col]])))
expr_genes <- expr_genes[nzchar(expr_genes)]
sample_ids <- sample_cols
expr_num <- as.matrix(expr[, ..sample_cols])
mode(expr_num) <- "numeric"

resource_audit <- rbindlist(list(
  data.table(
    resource = "primary_gene_expression_matrix",
    candidate_file = primary_expr,
    readable = readable_file(primary_expr),
    n_rows = nrow(expr),
    n_columns = ncol(expr),
    contains_gene_symbols_or_probe_ids = id_col,
    contains_sample_ids = paste(head(sample_ids, 6), collapse = ";"),
    recommended_primary_input = "yes",
    blocking_issue = "",
    notes = "Gene-level matrix selected for Stage 5B audit; scale/normalization requires manual confirmation before modeling."
  ),
  data.table(
    resource = "primary_sample_metadata",
    candidate_file = primary_meta,
    readable = readable_file(primary_meta),
    n_rows = nrow(meta),
    n_columns = ncol(meta),
    contains_gene_symbols_or_probe_ids = "not_applicable",
    contains_sample_ids = ifelse("sample_id" %in% names(meta), "sample_id", "missing_sample_id_column"),
    recommended_primary_input = "yes",
    blocking_issue = "",
    notes = "Curated metadata selected; labels are filename-derived and require cautious claim boundaries."
  ),
  data.table(
    resource = "probe_or_feature_gene_mapping",
    candidate_file = primary_mapping,
    readable = readable_file(primary_mapping),
    n_rows = ifelse(nrow(mapping) > 0, nrow(mapping), NA_integer_),
    n_columns = ifelse(ncol(mapping) > 0, ncol(mapping), NA_integer_),
    contains_gene_symbols_or_probe_ids = ifelse(all(c("feature_id", "gene_symbol") %in% names(mapping)), "feature_id;gene_symbol", "unknown"),
    contains_sample_ids = "not_applicable",
    recommended_primary_input = "supporting",
    blocking_issue = ifelse(readable_file(primary_mapping), "", "mapping_missing"),
    notes = "Feature mapping audit supports gene-level matrix provenance; platform annotation is incomplete."
  )
))
fwrite(resource_audit, file.path(table_dir, "gse73680_primary_input_audit.tsv"), sep = "\t")

role_guess <- function(col) {
  low <- tolower(col)
  if (low %in% c("sample_id", "geo_accession")) return("sample_id")
  if (grepl("patient", low)) return("patient_id")
  if (grepl("pair", low)) return("pair_id")
  if (grepl("group|status|context|plaque|stone|site|type", low)) return(ifelse(grepl("stone", low), "stone_type", ifelse(grepl("plaque|site|type", low), "tissue_status", "disease_group")))
  if (grepl("batch|platform", low)) return(ifelse(grepl("platform", low), "batch", "batch"))
  if (grepl("sex|age", low)) return("covariate")
  "unknown"
}

metadata_schema <- rbindlist(lapply(names(meta), function(col) {
  vals <- meta[[col]]
  nonmiss <- vals[!is.na(vals) & vals != "" & vals != "\"\""]
  role <- role_guess(col)
  data.table(
    metadata_column = col,
    non_missing_count = length(nonmiss),
    unique_count = uniqueN(nonmiss),
    example_values = paste(head(unique(as.character(nonmiss)), 6), collapse = ";"),
    possible_role = role,
    use_in_stage5B = fifelse(role %in% c("sample_id", "patient_id", "pair_id", "tissue_status", "disease_group", "stone_type", "batch", "covariate"), "candidate", "no"),
    notes = fifelse(col %in% c("raw_label", "notes"), "Useful for curation trace, not model covariate.", "")
  )
}))
fwrite(metadata_schema, file.path(table_dir, "gse73680_metadata_schema_audit.tsv"), sep = "\t")

include_col <- if ("include_in_analysis" %in% names(meta)) "include_in_analysis" else NA_character_
meta[, include_bool := if (is.na(include_col)) TRUE else as.logical(get(include_col))]
meta[is.na(include_bool), include_bool := FALSE]
if (!"patient_id" %in% names(meta)) meta[, patient_id := sample_id]
if (!"group_curated" %in% names(meta)) meta[, group_curated := disease_status]
meta[, pair_id := patient_id]

included_meta <- meta[include_bool == TRUE & sample_id %in% sample_ids]
pair_counts <- included_meta[, .(
  groups = paste(sort(unique(group_curated)), collapse = ";"),
  n_groups = uniqueN(group_curated),
  sample_ids_all = paste(sample_id, collapse = ";")
), by = patient_id]
included_meta <- merge(included_meta, pair_counts, by = "patient_id", all.x = TRUE)
included_meta[, is_paired := n_groups >= 2]
included_meta[, paired_counterpart_sample := vapply(seq_len(.N), function(i) {
  same <- included_meta[patient_id == included_meta$patient_id[i] & group_curated != included_meta$group_curated[i], sample_id]
  if (length(same)) paste(same, collapse = ";") else ""
}, character(1))]
if (!"plaque_status" %in% names(included_meta)) included_meta[, plaque_status := ""]
if (!"stone_status" %in% names(included_meta)) included_meta[, stone_status := ""]

pairing <- included_meta[, .(
  patient_id,
  pair_id = patient_id,
  sample_id,
  tissue_status = plaque_status,
  disease_group = group_curated,
  stone_type = stone_status,
  is_paired,
  paired_counterpart_sample,
  usable_for_paired_analysis = include_bool & is_paired & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla"),
  blocking_issue = fifelse(!include_bool, "excluded", fifelse(!is_paired, "unpaired_or_single_group_patient", "")),
  notes = "Pairing inferred from patient_id and curated group labels."
)]
setorder(pairing, patient_id, disease_group, sample_id)
fwrite(pairing, file.path(table_dir, "gse73680_sample_pairing_audit.tsv"), sep = "\t")

sample_group_summary <- rbindlist(list(
  included_meta[, .(group_variable = "group_curated", group_label = group_curated, n_samples = .N, n_patients = uniqueN(patient_id), n_paired_samples = sum(is_paired), notes = "Primary Stage 5B group variable."), by = group_curated][, group_curated := NULL],
  included_meta[, .(group_variable = "stone_status", group_label = fifelse(stone_status == "" | is.na(stone_status), "not_available_or_control", stone_status), n_samples = .N, n_patients = uniqueN(patient_id), n_paired_samples = sum(is_paired), notes = "Stone type subgrouping candidate; only reliable if metadata curation is accepted."), by = .(stone_status = fifelse(stone_status == "" | is.na(stone_status), "not_available_or_control", stone_status))][, stone_status := NULL],
  included_meta[, .(group_variable = "plaque_status", group_label = fifelse(plaque_status == "" | is.na(plaque_status), "not_available_or_control", plaque_status), n_samples = .N, n_patients = uniqueN(patient_id), n_paired_samples = sum(is_paired), notes = "Plaque/stone tissue context candidate."), by = .(plaque_status = fifelse(plaque_status == "" | is.na(plaque_status), "not_available_or_control", plaque_status))][, plaque_status := NULL]
), fill = TRUE)
fwrite(sample_group_summary, file.path(table_dir, "gse73680_sample_group_summary.tsv"), sep = "\t")

feature_mapping_status <- if (nrow(mapping) && all(c("feature_id", "gene_symbol") %in% names(mapping))) {
  paste0(mapping[gene_symbol != "" & !is.na(gene_symbol), .N], "/", nrow(mapping), " mapped feature rows")
} else "mapping_file_missing_or_schema_unclear"

gene_counts <- expr[, .N, by = get(id_col)]
setnames(gene_counts, "get", "gene")
dups <- gene_counts[gene != "" & N > 1]
qc <- data.table(
  qc_item = c(
    "n_features_raw",
    "n_features_after_gene_mapping",
    "n_samples",
    "missing_value_fraction",
    "duplicate_gene_handling",
    "log_scale_status",
    "normalization_status",
    "batch_variable_available",
    "platform_annotation_available",
    "gene_symbol_mapping_status",
    "recommended_expression_matrix_for_stage5B"
  ),
  value = c(
    as.character(if (!is.null(feature_expr)) nrow(feature_expr) else NA_integer_),
    as.character(length(expr_genes)),
    as.character(length(sample_cols)),
    sprintf("%.4f", mean(is.na(expr_num))),
    paste0(nrow(dups), " duplicated non-empty gene symbols in primary gene matrix"),
    ifelse(max(expr_num, na.rm = TRUE) > 1000, "not_log2_scale_or_mixed_normalized_intensity", "possibly_log_scale"),
    "normalized_continuous_or_unknown",
    ifelse(any(grepl("batch|platform", names(meta), ignore.case = TRUE)), "candidate_column_available_or_inferable", "not_available"),
    ifelse(readable_file(primary_mapping), "feature mapping file available but platform annotation incomplete", "not_available"),
    feature_mapping_status,
    primary_expr
  ),
  status = c(
    "observed",
    "observed",
    "observed",
    ifelse(mean(is.na(expr_num)) < 0.2, "warning", "needs_fix"),
    ifelse(nrow(dups) == 0, "pass", "requires_documented_collapse"),
    "requires_manual_check",
    "requires_manual_check",
    ifelse(any(grepl("batch|platform", names(meta), ignore.case = TRUE)), "candidate", "limited"),
    ifelse(readable_file(primary_mapping), "candidate", "limited"),
    ifelse(grepl("mapped", feature_mapping_status), "candidate", "limited"),
    "recommended_with_log2_transform_check"
  ),
  source_file = c(paths$feature_expr, primary_expr, primary_expr, primary_expr, primary_expr, primary_expr, primary_expr, primary_meta, primary_mapping, primary_mapping, primary_expr),
  notes = c(
    "Raw/feature-level matrix row count.",
    "Non-empty gene identifiers in selected gene-level matrix.",
    "Sample columns in selected gene-level matrix.",
    "Calculated on selected matrix numeric cells; missingness reflects reconstructed supplementary files.",
    "Stage 5B should document duplicate handling; do not silently average probes.",
    "Large values indicate raw/normalized intensity scale, so Stage 5B should log2-transform or verify prior processing.",
    "Existing scripts call this normalized continuous or unknown; independent verification still needed.",
    "No strong batch variable found beyond platform/sample filename context.",
    "Feature mapping exists, but original platform annotation is not complete.",
    "Mapping relies heavily on original feature IDs treated as gene symbols.",
    "Use this matrix for audit/design; Stage 5B should rerun QC and transform before final analysis."
  )
)
fwrite(qc, file.path(table_dir, "gse73680_expression_qc_audit.tsv"), sep = "\t")

if (nrow(mapping) && all(c("feature_id", "gene_symbol") %in% names(mapping))) {
  map_audit <- copy(mapping)
  map_audit[, gene_symbol_clean := fifelse(is.na(gene_symbol) | gene_symbol == "", NA_character_, gene_symbol)]
  map_audit[, n_probes_per_gene := ifelse(is.na(gene_symbol_clean), NA_integer_, .N), by = gene_symbol_clean]
  map_audit[, recommended_collapse_rule := fifelse(is.na(gene_symbol_clean), "exclude_unmapped", "if multiple probes per gene, use highest mean expression or highest variance; document choice")]
  map_audit[, mapping_status := fifelse(is.na(gene_symbol_clean), "unmapped", mapping_status)]
  map_out <- map_audit[, .(
    probe_id = feature_id,
    gene_symbol = fifelse(is.na(gene_symbol_clean), "", gene_symbol_clean),
    mapping_status,
    n_probes_per_gene,
    recommended_collapse_rule,
    notes
  )]
  fwrite(map_out, file.path(table_dir, "gse73680_probe_gene_mapping_audit.tsv"), sep = "\t")
}

writeLines(c(
  "# GSE73680 gene-level status memo",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "- Primary Stage 5A input is a reconstructed gene-level expression matrix: `data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz`.",
  "- Feature/probe mapping was audited separately in `gse73680_probe_gene_mapping_audit.tsv`.",
  "- The matrix still requires Stage 5B scale verification and an explicit log2 transform decision before final modeling.",
  "- Do not silently average duplicated gene/probe features; use highest mean expression or highest variance if collapsing is needed, and document the chosen rule."
), file.path(doc_dir, "gse73680_gene_level_status_memo.md"))

get_genes <- function(group) stage3[reporting_group == group, unique(gene)]
read_gene_set <- function(path) {
  if (!file.exists(path)) return(character())
  unique(fread(path, header = FALSE)[[1]])
}
module_defs <- list(
  R1_MAGMA_Bonferroni_only = get_genes("R1_MAGMA_Bonferroni_only"),
  R1_R2_R3_all_MAGMA_Bonferroni = stage3[reporting_group %in% c("R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"), unique(gene)],
  MAGMA_top50 = read_gene_set(paths$magma_top50),
  MAGMA_top100 = read_gene_set(paths$magma_top100),
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = get_genes("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy"),
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = get_genes("R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"),
  R2_R3_MAGMA_plus_TWAS_proxy = stage3[reporting_group %in% c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"), unique(gene)],
  R5_TWAS_proxy_only = get_genes("R5_TWAS_proxy_only"),
  Curated_exemplar_panel = exemplar$gene,
  injury_epithelial_signature = c("HAVCR1", "LCN2", "SPP1", "VCAM1", "KRT8", "KRT18", "VIM", "PROM1"),
  immune_signature = c("PTPRC", "AIF1", "LST1", "CD68", "CD14", "C1QA", "C1QB", "MS4A1", "CD3D", "CD3E"),
  stromal_ECM_signature = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "TAGLN", "PDGFRB"),
  endothelial_signature = c("PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CLDN5", "RAMP2"),
  LoopTAL_signature = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2"),
  mineralization_remodeling_signature = c("SPP1", "ALPL", "MMP7", "MMP9", "RUNX2", "BGLAP", "POSTN", "COL1A1", "COL1A2")
)
module_role <- c(
  R1_MAGMA_Bonferroni_only = "primary_MAGMA_module",
  R1_R2_R3_all_MAGMA_Bonferroni = "primary_MAGMA_module",
  MAGMA_top50 = "primary_MAGMA_module",
  MAGMA_top100 = "primary_MAGMA_module",
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R2_R3_MAGMA_plus_TWAS_proxy = "secondary_TWAS_proxy_module",
  R5_TWAS_proxy_only = "secondary_TWAS_proxy_module",
  Curated_exemplar_panel = "curated_exemplar_context",
  injury_epithelial_signature = "injury_or_remodeling_signature",
  immune_signature = "composition_signature",
  stromal_ECM_signature = "composition_signature",
  endothelial_signature = "composition_signature",
  LoopTAL_signature = "composition_signature",
  mineralization_remodeling_signature = "injury_or_remodeling_signature"
)

module_feas <- rbindlist(lapply(names(module_defs), function(nm) {
  genes <- unique(module_defs[[nm]])
  genes <- genes[nzchar(genes)]
  det <- intersect(genes, expr_genes)
  miss <- setdiff(genes, expr_genes)
  data.table(
    module_name = nm,
    module_role = unname(module_role[nm]),
    source_file_or_definition = ifelse(grepl("signature|LoopTAL|immune|stromal|endothelial|mineralization|injury", nm), "Stage 5A documented marker definition", "Stage 3R/Stage 4B module definition"),
    n_genes_input = length(genes),
    n_genes_detected_in_gse73680 = length(det),
    detected_fraction = ifelse(length(genes) > 0, length(det) / length(genes), NA_real_),
    detected_genes = paste(det, collapse = ";"),
    missing_genes = paste(miss, collapse = ";"),
    feasible_for_stage5B = fifelse(length(det) >= 3 | nm %in% c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R5_TWAS_proxy_only"), "yes_with_caution", "no_or_supplement_only"),
    notes = fifelse(length(det) < 5, "Small module; use only as supplementary/descriptive context.", "Feasible for Stage 5B module-score audit.")
  )
}), fill = TRUE)
fwrite(module_feas, file.path(table_dir, "gse73680_module_feasibility_audit.tsv"), sep = "\t")

composition_defs <- list(
  epithelial_general = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  LoopTAL = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2"),
  proximal_tubule = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "AQP1", "ALDOB"),
  collecting_duct = c("AQP2", "AQP3", "SCNN1A", "SCNN1B", "SCNN1G", "KRT19"),
  immune_myeloid = c("PTPRC", "AIF1", "LST1", "CD68", "CD14", "C1QA", "C1QB"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "CD2"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
  fibroblast_stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CLDN5"),
  pericyte_smooth_muscle = c("PDGFRB", "RGS5", "ACTA2", "TAGLN", "MYH11", "MCAM"),
  injury_epithelial = c("HAVCR1", "LCN2", "SPP1", "VCAM1", "KRT8", "KRT18", "VIM"),
  ECM_fibrosis = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "TGFBI", "MMP2"),
  mineralization_remodeling = c("SPP1", "ALPL", "MMP7", "MMP9", "RUNX2", "BGLAP", "POSTN")
)
composition_manifest <- rbindlist(lapply(names(composition_defs), function(nm) {
  markers <- unique(composition_defs[[nm]])
  det <- intersect(markers, expr_genes)
  data.table(
    signature_name = nm,
    marker_genes_input = paste(markers, collapse = ";"),
    n_marker_genes_detected = length(det),
    detected_fraction = length(det) / length(markers),
    source_of_markers = "Stage 5A documented canonical kidney/cell-state marker set; must be rechecked against project marker sources in Stage 5B",
    feasible_for_score = ifelse(length(det) >= 3, "yes", "limited"),
    notes = ifelse(length(det) >= 3, "Feasible for bulk composition/signature sensitivity score.", "Too few detected markers; use cautiously or replace with project-supported marker source.")
  )
}), fill = TRUE)
fwrite(composition_manifest, file.path(table_dir, "gse73680_composition_signature_manifest.tsv"), sep = "\t")

existing_result_specs <- data.table(
  analysis_type = c(
    "paired_limma_existing", "duplicateCorrelation_existing", "module_score_existing",
    "CAMERA_existing", "ROAST_existing", "GSEA_existing",
    "random_set_benchmark_existing", "composition_score_existing",
    "composition_adjusted_model_existing", "injury_marker_correlation_existing",
    "figure_source_existing"
  ),
  pattern = c(
    "paired|patient_level|paired_delta",
    "duplicateCorrelation|module_score_response|p1_gene_response",
    "module_score_matrix|module_score_response|patient_level_module",
    "camera",
    "roast",
    "gsea",
    "random.*benchmark|expression_matched",
    "composition|signature",
    "composition_adjusted|partial_correlations",
    "injury.*correlation|risk_injury",
    "figure4|source"
  )
)
existing_result_audit <- existing_result_specs[, {
  hits <- file_inventory[grepl(pattern, file_path, ignore.case = TRUE) | grepl(pattern, file_name, ignore.case = TRUE)]
  found <- nrow(hits) > 0
    .(
      existing_result_found = found,
      file_path = ifelse(found, paste(head(hits$file_path, 10), collapse = ";"), ""),
      can_reuse = ifelse(found, "audit_reference_only", "no_or_not_found"),
    needs_rerun = ifelse(analysis_type %in% c("composition_adjusted_model_existing", "CAMERA_existing", "ROAST_existing", "GSEA_existing") | found, "yes_for_stage5B_final", "yes_if_required"),
    reason = ifelse(found, "Existing outputs predate Stage 5A claim boundaries; reuse only as audit references unless rerun under Stage 5B design.", "No existing output located in inventory."),
    notes = ifelse(found, paste0("Located ", nrow(hits), " candidate files."), "Not found.")
  )
}, by = .(analysis_type, pattern)][, pattern := NULL]
fwrite(existing_result_audit, file.path(table_dir, "gse73680_existing_result_audit.tsv"), sep = "\t")

n_pairs <- uniqueN(pairing[usable_for_paired_analysis == TRUE, patient_id])
n_included <- nrow(included_meta)
n_group_control <- included_meta[group_curated == "control_or_adjacent", .N]
n_group_plaque <- included_meta[group_curated == "plaque_or_stone_papilla", .N]
expr_qc_status <- qc[qc_item == "log_scale_status", status]
stage5b_ready <- readable_file(primary_expr) && readable_file(primary_meta) && n_pairs >= 3

writeLines(c(
  "# Stage 5B statistical design plan",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## A. Primary Analysis Goal",
  "",
  "Test whether MAGMA-prioritized modules show plaque/stone papilla bulk disease-context association in GSE73680. This is a bulk-tissue context analysis, not cell-type-specific validation.",
  "",
  "## B. Primary Model Options",
  "",
  paste0("- Paired structure is available for ", n_pairs, " patients with both control/adjacent and plaque/stone papilla samples."),
  "- Preferred Stage 5B model: limma on log2-transformed gene/module scores with patient blocking or duplicateCorrelation.",
  "- Paired sensitivity: within-patient plaque/stone minus adjacent/control deltas for primary modules.",
  "- Unpaired model should be used only as a secondary fallback if pairing fails after metadata verification.",
  "",
  "## C. Primary Outputs",
  "",
  "- Sample QC/PCA and sample inclusion table.",
  "- Module score response for R1/R1-R3/MAGMA top50/MAGMA top100.",
  "- Paired delta plots and patient-level module response tables.",
  "- Gene-level response for selected modules as supplementary context only.",
  "- Composition signature scores and composition-adjusted module response.",
  "- Injury/remodeling coupling and sensitivity analysis.",
  "",
  "## D. Primary Claim Allowed",
  "",
  "Bulk plaque/stone papilla disease-context support for MAGMA-prioritized modules, if Stage 5B results remain directionally consistent after paired and composition-aware sensitivity checks.",
  "",
  "## E. Claims Not Allowed",
  "",
  "- Cell-type-specific disease response.",
  "- Genetic causality.",
  "- Plaque causation or plaque nucleation biology.",
  "- Therapeutic target validation.",
  "- Single-gene validation.",
  "",
  "## F. Required Stage 5B Analyses",
  "",
  "- Re-run sample QC/PCA using the selected matrix.",
  "- Confirm scale and apply explicit log2 transform or justify no transform.",
  "- Run paired or patient-blocked module response.",
  "- Score epithelial, Loop/TAL, immune, stromal, endothelial, injury and remodeling signatures.",
  "- Run composition-aware sensitivity models.",
  "- Treat TWAS/exemplar modules as supplementary only.",
  "- Prepare source-data-backed Figure 4 planning after analysis, not before.",
  "",
  "## Recommended Stage 5B Command",
  "",
  "```bash",
  "Rscript scripts/revision_utils/stage5B_gse73680_context/stage5B_gse73680_paired_composition_analysis.R",
  "```"
), file.path(doc_dir, "stage5B_statistical_design_plan.md"))

reviewer_answers <- c(
  "# Stage 5A simulated reviewer check",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "1. Is GSE73680 metadata sufficient for paired analysis?",
  paste0("   Yes, provisionally. The curated metadata supports ", n_pairs, " paired patients, but labels are filename-derived and should remain cautiously described."),
  "",
  "2. Are plaque/stone and adjacent/control labels clear?",
  "   Mostly yes for included samples: `control_or_adjacent` versus `plaque_or_stone_papilla`. The wording should stay plaque/stone papilla context, not plaque nucleation.",
  "",
  "3. Are CaOx/CaP or stone-type labels available?",
  "   Yes for many plaque/stone samples, including CaOx and CaOx_CaP labels, but subgrouping should be secondary because sample sizes shrink.",
  "",
  "4. Is expression matrix gene-level and normalized?",
  "   A gene-level reconstructed matrix exists. Scale/normalization is not independently verified; Stage 5B must explicitly check log2 transformation and normalization.",
  "",
  "5. Can MAGMA modules be mapped to GSE73680 genes?",
  "   Yes. Primary MAGMA modules have high gene detection fractions; small TWAS-proxy modules require caution.",
  "",
  "6. Can composition-aware analysis be done?",
  "   Feasible as bulk marker-signature sensitivity analysis. It cannot deconvolve true cell fractions without stronger reference/deconvolution validation.",
  "",
  "7. What would a reviewer criticize most?",
  "   Metadata curation from filenames, uncertain expression scale/normalization, and possible bulk composition confounding.",
  "",
  "8. What must Stage 5B prove before GSE73680 can support the manuscript?",
  "   MAGMA-prioritized module response must be directionally consistent in patient-blocked/paired analysis and not fully explained by broad injury/composition signatures.",
  "",
  "9. Which claims must remain disallowed?",
  "   Cell-type-specific response, genetic causality, plaque causation/nucleation, therapeutic target validation, and single-gene validation.",
  "",
  "10. Should Stage 5B start or are blockers present?",
  if (stage5b_ready) "   Stage 5B can start after human acceptance, with scale/metadata checks as first operations." else "   Blockers remain; Stage 5B should not start."
)
writeLines(reviewer_answers, file.path(doc_dir, "stage5A_simulated_reviewer_check.md"))

stage5a_report <- c(
  "# Stage 5A report: GSE73680 bulk papillary disease-context audit and design",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Resources Found",
  "",
  paste0("- GSE73680-related files inventoried: ", nrow(file_inventory)),
  paste0("- Existing result artifacts found: ", existing_result_audit[existing_result_found == TRUE, .N], " / ", nrow(existing_result_audit), " audit categories"),
  "",
  "## Primary Expression Matrix and Metadata",
  "",
  paste0("- Primary matrix: `", primary_expr, "`"),
  paste0("- Primary metadata: `", primary_meta, "`"),
  paste0("- Primary feature mapping: `", primary_mapping, "`"),
  "",
  "## Sample and Pairing Status",
  "",
  paste0("- Included samples: ", n_included),
  paste0("- Control/adjacent samples: ", n_group_control),
  paste0("- Plaque/stone papilla samples: ", n_group_plaque),
  paste0("- Paired patients with both groups: ", n_pairs),
  "",
  "## Group Labels",
  "",
  "- Primary comparison: `plaque_or_stone_papilla` versus `control_or_adjacent`.",
  "- Stone-type labels are available for secondary subgroup audit but should not be primary unless power is adequate.",
  "",
  "## Expression QC",
  "",
  paste0("- Non-empty gene identifiers: ", length(expr_genes)),
  paste0("- Sample columns: ", length(sample_cols)),
  paste0("- Missing value fraction: ", sprintf("%.4f", mean(is.na(expr_num)))),
  "- Scale/normalization status: requires manual check before Stage 5B final modeling.",
  "",
  "## Gene Mapping Status",
  "",
  paste0("- Feature mapping status: ", feature_mapping_status),
  "- Gene-level matrix is available; probe/gene mapping audit and gene-level memo generated.",
  "",
  "## Module Feasibility",
  "",
  paste0("- Modules/signatures audited: ", nrow(module_feas)),
  paste0("- Feasible or feasible-with-caution modules: ", module_feas[feasible_for_stage5B == "yes_with_caution", .N]),
  "- Primary MAGMA modules are feasible; very small TWAS-proxy modules should remain supplementary/descriptive.",
  "",
  "## Composition Signature Feasibility",
  "",
  paste0("- Composition/signature rows audited: ", nrow(composition_manifest)),
  paste0("- Feasible marker signatures: ", composition_manifest[feasible_for_score == "yes", .N]),
  "- These are bulk marker-signature scores, not validated cell fractions.",
  "",
  "## Existing Result Audit",
  "",
  "- Existing GSE73680 module, paired, random benchmark, injury-correlation, and figure artifacts are present.",
  "- Stage 5B should rerun final analyses under the new conservative claim boundaries rather than importing old conclusions.",
  "",
  "## Stage 5B Readiness",
  "",
  if (stage5b_ready) "Stage 5B is ready after human acceptance, with mandatory first steps: confirm scale/log2 transform, recheck metadata labels, and run paired/patient-blocked module response." else "Stage 5B is blocked by missing primary resources or insufficient paired samples.",
  "",
  "## Blockers",
  "",
  if (stage5b_ready) "- No hard resource blocker. Caution flags: filename-derived metadata, uncertain expression scale/normalization, and bulk composition confounding." else "- Hard resource blocker present; see resource blocker memo.",
  "",
  "## Exact Recommended Stage 5B Command",
  "",
  "```bash",
  "Rscript scripts/revision_utils/stage5B_gse73680_context/stage5B_gse73680_paired_composition_analysis.R",
  "```"
)
writeLines(stage5a_report, file.path(doc_dir, "stage5A_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 5, `:=`(
    status = "stage5A_completed_audit_and_design",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 5A GSE73680 resource inventory, primary input audit, metadata/pairing audit, expression QC audit, probe/gene mapping audit, module feasibility, composition signature manifest, existing-result audit, Stage 5B design plan, reviewer check, and report generated",
    blocking_issues = "No hard blocker if human accepts filename-derived metadata and requires Stage 5B scale/log2 verification; composition-aware analysis still required before strong disease-context support",
    next_stage_ready = ifelse(stage5b_ready, "stage5B_ready_after_human_acceptance", "no")
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Wrote Stage 5A outputs to", doc_dir, "and", table_dir, "\n")
cat("Stage 5B ready:", stage5b_ready, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
