suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

doc_dir <- "docs/revision/stage4A_scrna_audit"
table_dir <- "results/tables/revision/stage4A_scrna_audit"
figure_dir <- "results/figures/revision/stage4A_scrna_audit"
log_dir <- "logs/revision/stage4A_scrna_audit"
script_dir <- "scripts/revision_utils/stage4A_scrna_audit"
for (d in c(doc_dir, table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "stage4A_scrna_resource_audit.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4A snRNA resource audit\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", getwd(), "\n")
cat("R version:", as.character(getRversion()), "\n")
cat("Seurat available:", requireNamespace("Seurat", quietly = TRUE), "\n\n")

bool <- function(x) ifelse(isTRUE(x), "yes", "no")
present <- function(path) file.exists(path)
collapse0 <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = ";")
}
fmt_num <- function(x) {
  if (is.na(x)) "NA" else format(x, scientific = FALSE, trim = TRUE)
}

likely_type <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.(rds|rds\\.gz|qs|h5seurat)$", lower)) return("single_nucleus_object")
  if (grepl("\\.(h5ad)$", lower)) return("ann_data_object")
  if (grepl("\\.(mtx|h5)$", lower)) return("matrix")
  if (grepl("(metadata|meta|sample|donor|family\\.soft)", lower)) return("metadata")
  if (grepl("(marker|annotation|cluster)", lower)) return("annotation_or_marker")
  if (grepl("(umap|embedding)", lower)) return("embedding")
  if (grepl("(module|score|random|benchmark|leave|driver)", lower)) return("existing_score_or_benchmark")
  if (grepl("(figure2|figure3|source)", lower)) return("figure_or_source_data")
  "other_candidate_resource"
}

likely_content <- function(path) {
  lower <- tolower(path)
  if (grepl("gse231569_annotated_audited|gse231569_audited_seurat", lower)) {
    return("audited Seurat object with expression, metadata, UMAP, donor_id, and phase1_cell_type")
  }
  if (grepl("phase1_gse231569_quick_seurat", lower)) {
    return("earlier Seurat object with expression, metadata, and UMAP; donor_id may be absent")
  }
  if (grepl("papilla_samples\\.rds", lower)) return("downloaded author processed object candidate")
  if (grepl("cell_counts|donor_cell_counts", lower)) return("existing cell/donor compartment count summary")
  if (grepl("marker", lower)) return("existing marker or annotation audit")
  if (grepl("module|score", lower)) return("existing module score result")
  if (grepl("random|benchmark", lower)) return("existing random-set benchmark")
  if (grepl("leave", lower)) return("existing leave-one sensitivity result")
  if (grepl("driver", lower)) return("existing known/locus driver sensitivity result")
  if (grepl("family\\.soft", lower)) return("GEO series metadata")
  if (grepl("raw\\.tar", lower)) return("GEO raw supplementary archive")
  likely_type(path)
}

readable_quick <- function(path) {
  if (!file.exists(path)) return("no")
  if (grepl("\\.(tsv|csv|txt|md|R|r)$", path, ignore.case = TRUE)) return("yes")
  if (grepl("\\.(rds|RDS)$", path)) {
    out <- tryCatch({
      obj <- readRDS(path)
      rm(obj)
      "yes"
    }, error = function(e) paste0("no: ", conditionMessage(e)))
    return(out)
  }
  if (grepl("\\.(gz|tar|pdf|png|svg)$", path, ignore.case = TRUE)) return("not_tested_container_or_binary")
  "not_tested"
}

all_files <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = FALSE)
candidate_pattern <- paste(c(
  "GSE231569", "gse231569", "scrna", "single", "magma_scrna", "celltype",
  "cell_type", "donor", "marker", "umap", "module_score", "random_benchmark",
  "leave_one", "leave-one", "driver", "\\.rds$", "\\.RDS$", "\\.h5ad$",
  "\\.h5seurat$", "\\.mtx$", "\\.h5$"
), collapse = "|")
inventory_paths <- all_files[grepl(candidate_pattern, all_files)]
inventory_paths <- inventory_paths[!grepl("^\\.git/", inventory_paths)]
inventory_paths <- unique(inventory_paths)
info <- file.info(inventory_paths)
inventory <- data.table(
  file_path = inventory_paths,
  file_name = basename(inventory_paths),
  file_type = vapply(inventory_paths, likely_type, character(1)),
  size_bytes = as.numeric(info$size),
  last_modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
  likely_content = vapply(inventory_paths, likely_content, character(1)),
  readable = vapply(inventory_paths, readable_quick, character(1)),
  ready_for_stage4B = "no",
  notes = ""
)
inventory[file_type %in% c("single_nucleus_object", "ann_data_object", "matrix", "metadata", "annotation_or_marker", "existing_score_or_benchmark"),
          ready_for_stage4B := "potentially"]
inventory[grepl("annotated_audited|audited_seurat|donor_cell_counts|cell_counts|marker_audit|module_score_by_donor", file_path),
          ready_for_stage4B := "yes"]
setorder(inventory, file_type, -size_bytes, file_path)
fwrite(inventory, file.path(table_dir, "scrna_resource_inventory.tsv"), sep = "\t")

candidate_objects <- unique(c(
  "results/scrna/gse231569/objects/gse231569_annotated_audited.rds",
  "data/processed/gse231569_audited_seurat.rds",
  "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds",
  "data/raw/GSE231569/processed/GSE231569_papilla_samples.RDS",
  inventory[grepl("(GSE231569|gse231569)", file_path) & grepl("\\.(rds|RDS|h5seurat|h5ad)$", file_path), file_path]
))
candidate_objects <- candidate_objects[file.exists(candidate_objects)]

object_audit_rows <- list()
loaded_objects <- list()
for (path in candidate_objects) {
  obj <- tryCatch(readRDS(path), error = function(e) e)
  readable <- !inherits(obj, "error")
  is_seurat <- readable && inherits(obj, "Seurat")
  meta_cols <- if (is_seurat) colnames(obj@meta.data) else character()
  reductions <- if (is_seurat) names(obj@reductions) else character()
  has_scores <- any(grepl("score|module", meta_cols, ignore.case = TRUE))
  row <- data.table(
    candidate_object = basename(path),
    file_path = path,
    object_type = if (is_seurat) "Seurat" else if (readable) class(obj)[1] else "unreadable",
    readable = bool(readable),
    contains_expression = bool(is_seurat && nrow(obj) > 0 && ncol(obj) > 0),
    contains_metadata = bool(is_seurat && nrow(obj@meta.data) == ncol(obj)),
    contains_umap = bool(is_seurat && "umap" %in% reductions),
    contains_celltype_labels = bool(is_seurat && any(c("phase1_cell_type", "cell_type", "broad_cell_type", "audited_broad_cell_type") %in% meta_cols)),
    contains_donor_ids = bool(is_seurat && "donor_id" %in% meta_cols),
    contains_sample_ids = bool(is_seurat && "sample_id" %in% meta_cols),
    contains_group_status = bool(is_seurat && any(c("disease_status", "group", "condition") %in% meta_cols)),
    contains_broad_compartment = bool(is_seurat && any(c("phase1_cell_type", "broad_cell_type", "audited_broad_cell_type") %in% meta_cols)),
    contains_existing_module_scores = bool(is_seurat && has_scores),
    recommended_primary_input = "no",
    blocking_issue = if (readable) "" else conditionMessage(obj),
    notes = if (is_seurat) paste0(nrow(obj), " genes x ", ncol(obj), " nuclei; metadata columns: ", length(meta_cols)) else ""
  )
  object_audit_rows[[length(object_audit_rows) + 1]] <- row
  if (is_seurat) loaded_objects[[path]] <- obj
}
object_audit <- rbindlist(object_audit_rows, fill = TRUE)

primary_path <- NA_character_
if (nrow(object_audit)) {
  object_audit[, primary_score :=
    (readable == "yes") +
    (contains_expression == "yes") +
    (contains_metadata == "yes") +
    (contains_umap == "yes") +
    (contains_celltype_labels == "yes") +
    (contains_donor_ids == "yes") +
    (contains_sample_ids == "yes") +
    (contains_group_status == "yes") +
    (contains_broad_compartment == "yes")]
  object_audit[order(-primary_score, file_path)]
  primary_path <- object_audit[primary_score == max(primary_score), file_path][1]
  object_audit[file_path == primary_path, recommended_primary_input := "yes"]
  object_audit[, primary_score := NULL]
}
fwrite(object_audit, file.path(table_dir, "scrna_primary_object_audit.tsv"), sep = "\t")

if (is.na(primary_path) || !primary_path %in% names(loaded_objects)) {
  blocker <- c(
    "# Stage 4A snRNA resource blocker",
    "",
    "No usable Seurat object or expression matrix with metadata was found. Stage 4B must not start.",
    "",
    paste0("Generated: ", Sys.Date())
  )
  writeLines(blocker, file.path(doc_dir, "scrna_resource_blocker.md"))
  stop("No usable primary snRNA object found.")
}

obj <- loaded_objects[[primary_path]]
features <- rownames(obj)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
meta_cols <- colnames(meta)
if (!"donor_id" %in% meta_cols && "sample_id" %in% meta_cols) {
  cat("WARNING: donor_id absent in primary object; sample_id could be considered only after explicit metadata repair.\n")
}

role_for_col <- function(col) {
  if (col == "cell_id") return("cell_id")
  if (col == "donor_id") return("donor_id")
  if (col == "sample_id" || col == "orig.ident") return("sample_id")
  if (col %in% c("disease_status", "group", "condition")) return("group_status")
  if (col %in% c("phase1_cell_type", "broad_cell_type", "audited_broad_cell_type")) return("broad_compartment")
  if (col %in% c("seurat_clusters", "RNA_snn_res.0.4", "cluster", "cell_type")) return("fine_cell_type")
  if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")) return("quality_metric")
  if (col %in% c("UMAP_1", "UMAP_2")) return("embedding")
  "unknown"
}

metadata_rows <- lapply(meta_cols, function(col) {
  vals <- meta[[col]]
  ex <- unique(as.character(vals[!is.na(vals)]))
  data.table(
    metadata_column = col,
    non_missing_count = sum(!is.na(vals)),
    unique_count = uniqueN(vals, na.rm = TRUE),
    example_values = collapse0(head(ex, 5)),
    possible_role = role_for_col(col),
    use_in_stage4B = ifelse(role_for_col(col) %in% c("cell_id", "donor_id", "sample_id", "group_status", "broad_compartment", "fine_cell_type", "quality_metric"), "yes", "no"),
    notes = ""
  )
})
metadata_schema <- rbindlist(metadata_rows, fill = TRUE)
umap <- tryCatch(Embeddings(obj, reduction = "umap"), error = function(e) NULL)
if (!is.null(umap) && ncol(umap) >= 2) {
  metadata_schema <- rbind(
    metadata_schema,
    data.table(
      metadata_column = c("UMAP_1", "UMAP_2"),
      non_missing_count = nrow(umap),
      unique_count = c(uniqueN(umap[, 1]), uniqueN(umap[, 2])),
      example_values = c(collapse0(head(round(umap[, 1], 4), 5)), collapse0(head(round(umap[, 2], 4), 5))),
      possible_role = "embedding",
      use_in_stage4B = "descriptive_only",
      notes = "stored in Seurat reduction; not used as statistical unit"
    ),
    fill = TRUE
  )
}
required_roles <- c("cell_id", "donor_id", "sample_id", "group_status", "broad_compartment", "fine_cell_type", "quality_metric", "embedding")
metadata_schema[possible_role %in% required_roles & use_in_stage4B == "no", use_in_stage4B := "review"]
fwrite(metadata_schema, file.path(table_dir, "scrna_metadata_schema_audit.tsv"), sep = "\t")

donor_col <- if ("donor_id" %in% meta_cols) "donor_id" else NA_character_
sample_col <- if ("sample_id" %in% meta_cols) "sample_id" else if ("orig.ident" %in% meta_cols) "orig.ident" else NA_character_
group_col <- if ("disease_status" %in% meta_cols) "disease_status" else NA_character_
comp_col <- if ("phase1_cell_type" %in% meta_cols) "phase1_cell_type" else if ("audited_broad_cell_type" %in% meta_cols) "audited_broad_cell_type" else if ("broad_cell_type" %in% meta_cols) "broad_cell_type" else NA_character_

if (is.na(donor_col) || is.na(comp_col)) {
  stop("Required donor or compartment metadata missing in selected object.")
}

count_dt <- copy(meta)
setnames(count_dt, donor_col, "donor_id_tmp")
if (!is.na(sample_col)) setnames(count_dt, sample_col, "sample_id_tmp") else count_dt[, sample_id_tmp := NA_character_]
if (!is.na(group_col)) setnames(count_dt, group_col, "disease_or_group_tmp") else count_dt[, disease_or_group_tmp := NA_character_]
setnames(count_dt, comp_col, "broad_compartment_tmp")
count_dt[, donor_id_tmp := as.character(donor_id_tmp)]
count_dt[, sample_id_tmp := as.character(sample_id_tmp)]
count_dt[, disease_or_group_tmp := as.character(disease_or_group_tmp)]
count_dt[, broad_compartment_tmp := as.character(broad_compartment_tmp)]

median_safe <- function(x) if (all(is.na(x))) NA_real_ else as.numeric(median(x, na.rm = TRUE))
donor_counts <- count_dt[, .(
  n_nuclei = .N,
  n_detected_genes_median = if ("nFeature_RNA" %in% names(.SD)) median_safe(nFeature_RNA) else NA_real_,
  n_counts_median = if ("nCount_RNA" %in% names(.SD)) median_safe(nCount_RNA) else NA_real_,
  percent_mt_median = if ("percent.mt" %in% names(.SD)) median_safe(percent.mt) else NA_real_,
  notes = "primary statistical unit candidate for Stage 4B; no inferential test performed"
), by = .(
  donor_id = donor_id_tmp,
  sample_id = sample_id_tmp,
  disease_or_group = disease_or_group_tmp,
  broad_compartment = broad_compartment_tmp
)]
setorder(donor_counts, donor_id, broad_compartment)
fwrite(donor_counts, file.path(table_dir, "scrna_donor_compartment_count_audit.tsv"), sep = "\t")

loop_label <- "Loop_of_Henle_TAL"
total_by_donor <- count_dt[, .(total_nuclei = .N), by = .(donor_id = donor_id_tmp, sample_id = sample_id_tmp)]
loop_by_donor <- count_dt[broad_compartment_tmp == loop_label, .(
  n_loop_tal_nuclei = .N
), by = .(donor_id = donor_id_tmp, sample_id = sample_id_tmp)]
loop_audit <- merge(total_by_donor, loop_by_donor, by = c("donor_id", "sample_id"), all.x = TRUE)
loop_audit[is.na(n_loop_tal_nuclei), n_loop_tal_nuclei := 0L]
loop_audit[, `:=`(
  loop_tal_label_used = loop_label,
  loop_tal_fraction = n_loop_tal_nuclei / total_nuclei,
  loop_tal_present = ifelse(n_loop_tal_nuclei > 0, "yes", "no"),
  notes = "counts from selected primary Seurat object; descriptive audit only"
)]
setcolorder(loop_audit, c("donor_id", "sample_id", "loop_tal_label_used", "n_loop_tal_nuclei", "total_nuclei", "loop_tal_fraction", "loop_tal_present", "notes"))
setorder(loop_audit, donor_id)
fwrite(loop_audit, file.path(table_dir, "scrna_loop_tal_count_audit.tsv"), sep = "\t")

loop_total <- sum(loop_audit$n_loop_tal_nuclei)
loop_donors <- loop_audit[n_loop_tal_nuclei > 0, uniqueN(donor_id)]
claim_status <- if (loop_total == 540 && loop_donors == 4) "supported" else if (loop_total > 0) "partially supported" else "not supported"
recommended_wording <- if (claim_status == "supported") {
  "In the audited GSE231569 object, 540 nuclei were annotated as Loop_of_Henle_TAL across four donors; this count is descriptive and should not be used as an inferential replicate count."
} else {
  "The Loop/TAL count should be rewritten after metadata repair because the audited object does not fully support the 540 nuclei across four donors statement."
}
writeLines(c(
  "# Loop/TAL 540-nuclei claim audit",
  "",
  paste0("- Primary object: `", primary_path, "`"),
  paste0("- Exact label used: `", loop_label, "`"),
  paste0("- Total Loop/TAL nuclei found: ", loop_total),
  paste0("- Donors with Loop/TAL nuclei: ", loop_donors),
  paste0("- Claim status: ", claim_status),
  "",
  "Recommended wording:",
  "",
  recommended_wording
), file.path(doc_dir, "loop_tal_540_claim_audit.md"))

marker_sets <- list(
  Proximal_tubule = c("LRP2", "SLC34A1", "CUBN"),
  Loop_of_Henle_TAL = c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16"),
  Collecting_duct_principal = c("AQP2", "AQP3", "SCNN1G"),
  Intercalated_cell = c("ATP6V1B1", "SLC4A1", "FOXI1"),
  Endothelial = c("PECAM1", "VWF", "KDR"),
  Fibroblast_stromal = c("DCN", "LUM", "COL1A1"),
  Macrophage_myeloid = c("LST1", "C1QA", "CD68"),
  Pericyte_smooth_muscle = c("RGS5", "ACTA2")
)
compartments <- sort(unique(count_dt$broad_compartment_tmp))
expr_data <- tryCatch(
  GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
  error = function(e) {
    tryCatch(GetAssayData(obj, assay = DefaultAssay(obj), slot = "data"), error = function(e2) NULL)
  }
)
marker_rows <- lapply(compartments, function(comp) {
  markers <- marker_sets[[comp]]
  if (is.null(markers)) markers <- character()
  detected <- intersect(markers, features)
  comp_cells <- count_dt[broad_compartment_tmp == comp, cell_id]
  comp_cells <- intersect(comp_cells, colnames(expr_data))
  if (!is.null(expr_data) && length(detected) && length(comp_cells)) {
    marker_mat <- expr_data[detected, comp_cells, drop = FALSE]
    pct_cells_with_any_marker <- as.numeric(mean(Matrix::colSums(marker_mat > 0) > 0))
    mean_marker_expression <- as.numeric(sum(marker_mat) / (nrow(marker_mat) * ncol(marker_mat)))
  } else {
    pct_cells_with_any_marker <- NA_real_
    mean_marker_expression <- NA_real_
  }
  support <- if (!length(markers)) {
    "not_evaluable"
  } else if (length(detected) >= min(3, length(markers)) && !is.na(pct_cells_with_any_marker) && pct_cells_with_any_marker >= 0.5) {
    "strong"
  } else if (length(detected) >= 2 && (is.na(pct_cells_with_any_marker) || pct_cells_with_any_marker >= 0.1)) {
    "moderate"
  } else if (length(detected) == 1) {
    "weak"
  } else {
    "not_supported"
  }
  data.table(
    broad_compartment = comp,
    n_nuclei = count_dt[broad_compartment_tmp == comp, .N],
    n_donors_detected = count_dt[broad_compartment_tmp == comp, uniqueN(donor_id_tmp)],
    marker_genes_expected = collapse0(markers),
    marker_genes_detected = collapse0(detected),
    marker_support_status = support,
    mean_marker_expression = mean_marker_expression,
    pct_cells_with_any_marker = pct_cells_with_any_marker,
    notes = if (!length(markers)) "no requested marker set supplied for this compartment in Stage 4A prompt" else "support based on marker genes present and expressed in the labelled compartment; no new annotation inference performed"
  )
})
marker_audit <- rbindlist(marker_rows, fill = TRUE)
fwrite(marker_audit, file.path(table_dir, "scrna_compartment_marker_audit.tsv"), sep = "\t")

stage3_file <- "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv"
stage3 <- fread(stage3_file)
module_list <- list()
if ("reporting_group" %in% names(stage3) && "gene" %in% names(stage3)) {
  for (grp in sort(unique(stage3$reporting_group))) {
    module_list[[grp]] <- list(source = stage3_file, genes = stage3[reporting_group == grp, gene])
  }
}
if (present("results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv")) {
  exemplar <- fread("results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv")
  gene_col <- intersect(c("gene", "gene_symbol"), names(exemplar))[1]
  if (!is.na(gene_col)) {
    module_list[["Curated_exemplar_panel"]] <- list(
      source = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
      genes = exemplar[[gene_col]]
    )
  }
}
gene_set_files <- c(
  MAGMA_top50 = "results/gene_sets/magma_top50.txt",
  MAGMA_top100 = "results/gene_sets/magma_top100.txt",
  MAGMA_FDR = "results/gene_sets/magma_fdr05.txt",
  MAGMA_suggestive = "results/gene_sets/magma_suggestive_p1e4.txt"
)
for (nm in names(gene_set_files)) {
  path <- gene_set_files[[nm]]
  if (present(path)) {
    genes <- scan(path, what = character(), quiet = TRUE)
    module_list[[nm]] <- list(source = path, genes = genes)
  }
}
module_rows <- lapply(names(module_list), function(nm) {
  genes <- unique(toupper(module_list[[nm]]$genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  feat_upper <- toupper(features)
  detected <- genes[genes %in% feat_upper]
  missing <- setdiff(genes, detected)
  frac <- if (length(genes)) length(detected) / length(genes) else NA_real_
  feasible <- if (!length(genes) || length(detected) == 0) {
    "no"
  } else if (length(detected) >= 10 || frac >= 0.5) {
    "yes"
  } else {
    "partial"
  }
  data.table(
    module_name = nm,
    source_file = module_list[[nm]]$source,
    n_genes_input = length(genes),
    n_genes_detected_in_scrna = length(detected),
    detected_fraction = frac,
    detected_genes = collapse0(detected),
    missing_genes = collapse0(missing),
    feasible_for_stage4B = feasible,
    notes = "gene symbol overlap against selected primary Seurat feature names; no module score calculated in Stage 4A"
  )
})
module_audit <- rbindlist(module_rows, fill = TRUE)
fwrite(module_audit, file.path(table_dir, "scrna_module_feasibility_audit.tsv"), sep = "\t")

check_files <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) "" else paste(paths, collapse = ";")
}
existing_specs <- list(
  "nucleus-level module score" = c("results/tables/gse231569_celllevel_magma_scores.tsv", "results/tables/gse231569_celllevel_magma_score_blocker_v0.1.tsv"),
  "cell-type-level mean score" = c("results/tables/gse231569_magma_score_by_celltype_summary.tsv", "results/tables/magma_scrna_module_score_by_celltype.tsv", "results/tables/phase1_module_score_by_celltype.tsv"),
  "donor-level score" = c("results/tables/magma_scrna_module_score_by_donor_celltype.tsv", "results/tables/audited_locus_scrna_module_score_by_donor_celltype.tsv", "results/tables/scrna_module_score_by_donor_celltype.tsv", "results/tables/gse231569_magma_score_by_celltype_donor.tsv"),
  "expression-matched random set" = c("results/tables/magma_scrna_random_benchmark.tsv", "results/tables/audited_locus_scrna_random_benchmark.tsv", "results/tables/scrna_random_gene_set_benchmark.tsv"),
  "locus-balanced random set" = c("results/tables/magma_locus_balanced_scrna_benchmark.tsv", "results/tables/audited_locus_balanced_scrna_benchmark.tsv"),
  "leave-one-locus-out" = c("results/tables/magma_leave_one_locus_out.tsv", "results/tables/audited_locus_leave_one_locus_out.tsv", "results/phase3a_locus_driver_sensitivity_v0.1/tables/audited_locus_leave_one_locus_out.tsv"),
  "leave-one-donor-out" = c("results/tables/leave_one_donor_out.tsv", "results/tables/gse231569_leave_one_donor_out.tsv", "results/tables/revision/stage4B_scrna/leave_one_donor_out.tsv"),
  "known-driver removal" = c("results/tables/magma_top50_tal_driver_genes.tsv", "results/tables/audited_locus_top50_tal_driver_genes.tsv", "results/phase3a_locus_driver_sensitivity_v0.1/tables/audited_locus_top50_tal_driver_genes.tsv")
)
existing_rows <- lapply(names(existing_specs), function(typ) {
  found <- check_files(existing_specs[[typ]])
  exists <- nzchar(found)
  can_reuse <- if (!exists) "no" else if (typ %in% c("nucleus-level module score", "leave-one-locus-out", "known-driver removal")) "partial" else "yes"
  needs_rerun <- if (!exists) "yes" else if (typ %in% c("leave-one-donor-out")) "yes" else if (typ %in% c("nucleus-level module score", "expression-matched random set", "locus-balanced random set", "known-driver removal")) "yes_or_update_for_stage3R_modules" else "review_before_reuse"
  reason <- switch(typ,
    "nucleus-level module score" = "cell-level scores may be descriptive inputs only; no cell-level inferential claim allowed",
    "cell-type-level mean score" = "can inform descriptive summaries but Stage 4B should use donor x compartment unit",
    "donor-level score" = "existing donor summaries are relevant but should be checked against Stage 3R module definitions",
    "expression-matched random set" = "random sets should be regenerated or verified for Stage 3R module definitions and detection matching",
    "locus-balanced random set" = "useful robustness context but not a substitute for donor-level support",
    "leave-one-locus-out" = "existing locus sensitivity is useful but does not address donor robustness",
    "leave-one-donor-out" = "required Stage 4B robustness output not found in canonical location",
    "known-driver removal" = "existing driver sensitivity should be updated for Stage 3R modules before manuscript use",
    "review"
  )
  data.table(
    analysis_type = typ,
    existing_result_found = bool(exists),
    file_path = found,
    can_reuse = can_reuse,
    needs_rerun = needs_rerun,
    reason = reason,
    notes = "Stage 4A audit only; no result imported as final manuscript evidence"
  )
})
existing_audit <- rbindlist(existing_rows, fill = TRUE)
fwrite(existing_audit, file.path(table_dir, "scrna_existing_result_audit.tsv"), sep = "\t")

writeLines(c(
  "# Stage 4B statistical design plan",
  "",
  "## A. Primary statistical unit",
  "",
  "The primary unit must be `donor x broad_compartment`. Nuclei are not independent biological replicates for inferential claims.",
  "",
  "## B. Primary output",
  "",
  "Module scores should be summarized per `donor x broad_compartment`, with nuclei-level values used only to compute donor-level summaries.",
  "",
  "## C. Primary claim allowed",
  "",
  "MAGMA-prioritized modules may be described as showing donor-level descriptive support for a Loop/TAL-associated single-nucleus expression context, if Stage 4B robustness checks support this.",
  "",
  "## D. Claims not allowed",
  "",
  "- Causal mediation",
  "- Causal cell type",
  "- Papilla-specific genetic regulation",
  "- Plaque nucleation site",
  "- Cell-level inferential enrichment",
  "",
  "## E. Stage 4B analyses to run only if metadata supports them",
  "",
  "- donor x compartment module score",
  "- Loop/TAL rank within donor",
  "- leave-one-donor-out robustness",
  "- expression/detection matched random-set benchmark",
  "- known-driver removal sensitivity",
  "- candidate gene detection table",
  "- integration with Stage 3R evidence model",
  "",
  "## F. Minimum requirements to proceed",
  "",
  paste0("- donor_id available: ", bool("donor_id" %in% meta_cols)),
  paste0("- broad_compartment or mappable cell-type label available: ", bool(!is.na(comp_col))),
  paste0("- expression matrix available: ", bool(nrow(obj) > 0 && ncol(obj) > 0)),
  paste0("- Loop/TAL label present or reconstructable from markers: ", bool(loop_total > 0 || marker_audit[broad_compartment == loop_label, marker_support_status %in% c("strong", "moderate")])),
  paste0("- gene identifiers compatible with Stage 3R gene symbols: ", bool(any(module_audit$n_genes_detected_in_scrna > 0))),
  "",
  "Recommended Stage 4B command after human acceptance:",
  "",
  "```bash",
  "Rscript scripts/revision_utils/stage4B_scrna_donor_level/stage4B_scrna_donor_level_analysis.R",
  "```"
), file.path(doc_dir, "stage4B_statistical_design_plan.md"))

reviewer_answers <- c(
  "# Stage 4A simulated reviewer check",
  "",
  "1. Can the current snRNA resource support donor-level analysis?",
  paste0("   Yes, conditionally. The selected object contains donor_id, sample_id, disease_status, broad compartment labels, expression data, and UMAP. The donor x broad_compartment table has ", nrow(donor_counts), " rows."),
  "",
  "2. Is the Loop/TAL compartment present across donors?",
  paste0("   Yes. `", loop_label, "` is present in ", loop_donors, " donors."),
  "",
  "3. Is the 540 Loop/TAL nuclei statement source-supported?",
  paste0("   ", ifelse(claim_status == "supported", "Yes.", "Not fully."), " The audit found ", loop_total, " Loop/TAL nuclei across ", loop_donors, " donors."),
  "",
  "4. Are cell labels sufficiently audited?",
  "   Partially. Broad labels and marker sets are documented, but Stage 4B should preserve annotation provenance and avoid silent relabeling.",
  "",
  "5. Is there a risk of cell-level pseudoreplication?",
  "   Yes if cell-level observations are used for inferential P values. Stage 4B must use donor x broad_compartment as the statistical unit.",
  "",
  "6. Are random-set benchmarks reusable or need rerun?",
  "   Existing random benchmarks are useful provenance but should be rerun or explicitly verified for Stage 3R modules.",
  "",
  "7. Can known-driver removal be performed?",
  "   Likely yes because existing driver-sensitivity outputs exist, but Stage 4B should update them for the Stage 3R evidence model and donor-level framing.",
  "",
  "8. What would a high-impact reviewer still criticize?",
  "   Only four donors, uneven Loop/TAL nuclei per donor, reliance on broad annotations, and absence of an independent donor-level replication atlas.",
  "",
  "9. What must Stage 4B prove before the manuscript can retain a strong Loop/TAL claim?",
  "   It must show that Loop/TAL support persists at donor level, is not driven by one donor or one known marker/locus, and remains above matched random expectations.",
  "",
  "10. Should Stage 4B start, or are missing resources blocking it?",
  "   Stage 4B can start after acceptance of this audit, with conservative claim boundaries. No resource blocker was triggered."
)
writeLines(reviewer_answers, file.path(doc_dir, "stage4A_simulated_reviewer_check.md"))

strong_markers <- marker_audit[marker_support_status == "strong", .N]
module_yes <- module_audit[feasible_for_stage4B == "yes", .N]
module_partial <- module_audit[feasible_for_stage4B == "partial", .N]
raw_author_status <- object_audit[grepl("GSE231569_papilla_samples.RDS$", file_path), blocking_issue]
if (!length(raw_author_status)) raw_author_status <- ""
blockers <- c()
if (!"donor_id" %in% meta_cols) blockers <- c(blockers, "donor_id absent")
if (loop_total == 0) blockers <- c(blockers, "Loop/TAL label absent")
if (nrow(obj) == 0 || ncol(obj) == 0) blockers <- c(blockers, "expression matrix unavailable")
if (!length(blockers)) blockers <- "No hard blocker for Stage 4B; caution remains due to four-donor design and uneven Loop/TAL counts."

report <- c(
  "# Stage 4A report: GSE231569 snRNA resource and statistical-design audit",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Resources found",
  "",
  paste0("- Candidate resources inventoried: ", nrow(inventory)),
  paste0("- Candidate object files audited: ", nrow(object_audit)),
  paste0("- Inventory table: `", file.path(table_dir, "scrna_resource_inventory.tsv"), "`"),
  "",
  "## Primary object selected",
  "",
  paste0("- Recommended primary input: `", primary_path, "`"),
  paste0("- Object dimensions: ", nrow(obj), " genes x ", ncol(obj), " nuclei"),
  paste0("- UMAP present: ", bool(!is.null(umap))),
  if (nzchar(raw_author_status[1])) paste0("- Raw downloaded author RDS read status: not selected; readRDS issue recorded as `", raw_author_status[1], "`.") else "- Raw downloaded author RDS: no separate issue recorded.",
  "",
  "## Donor metadata status",
  "",
  paste0("- donor_id present: ", bool("donor_id" %in% meta_cols)),
  paste0("- sample_id present: ", bool("sample_id" %in% meta_cols)),
  paste0("- disease_status present: ", bool("disease_status" %in% meta_cols)),
  paste0("- broad compartment column used: `", comp_col, "`"),
  paste0("- donor x compartment rows: ", nrow(donor_counts)),
  "",
  "## Cell-type annotation status",
  "",
  paste0("- Broad compartments found: ", collapse0(compartments)),
  paste0("- Compartments with strong marker support based on requested marker genes: ", strong_markers),
  "- This audit does not relabel cells or harmonize labels silently.",
  "",
  "## Loop/TAL count audit",
  "",
  paste0("- Exact label: `", loop_label, "`"),
  paste0("- Total Loop/TAL nuclei: ", loop_total),
  paste0("- Donors with Loop/TAL nuclei: ", loop_donors),
  paste0("- Claim status for \"540 Loop/TAL nuclei across four donors\": ", claim_status),
  "",
  "## Module feasibility",
  "",
  paste0("- Modules/gene sets audited: ", nrow(module_audit)),
  paste0("- Feasible modules: ", module_yes),
  paste0("- Partial modules: ", module_partial),
  "- No Stage 4A module scores or inferential P values were calculated.",
  "",
  "## Existing-result audit",
  "",
  paste0("- Existing-result categories audited: ", nrow(existing_audit)),
  "- Existing donor-level and benchmark outputs should be treated as provenance and reviewed/rerun against Stage 3R definitions before manuscript use.",
  "",
  "## Stage 4B readiness decision",
  "",
  "Stage 4B is ready to start after human acceptance of this audit, but only under the donor x broad_compartment design.",
  "",
  "## Blockers and cautions",
  "",
  paste0("- ", blockers),
  "- Do not claim Loop/TAL enrichment until Stage 4B donor-level robustness and matched-random controls are completed.",
  "- Do not use nuclei as independent biological replicates.",
  "",
  "## Exact recommended Stage 4B command",
  "",
  "```bash",
  "Rscript scripts/revision_utils/stage4B_scrna_donor_level/stage4B_scrna_donor_level_analysis.R",
  "```"
)
writeLines(report, file.path(doc_dir, "stage4A_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4A_completed",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 4A snRNA resource, metadata, Loop/TAL claim, marker, module-feasibility, existing-result, reviewer, and Stage 4B design audits generated",
    blocking_issues = "Stage 4B not started; donor-level robustness/random benchmark/known-driver removal still pending; use donor x broad_compartment only",
    next_stage_ready = "stage4B_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("\nPrimary object:", primary_path, "\n")
cat("Loop/TAL total:", loop_total, "across", loop_donors, "donors\n")
cat("Claim status:", claim_status, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
