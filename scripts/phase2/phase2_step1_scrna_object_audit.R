#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

root <- getwd()
canonical_object <- file.path(root, "data/processed/gse231569_audited_seurat.rds")
out_dir <- file.path(root, "results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

collapse_values <- function(x, max_n = 8) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  ux <- unique(x)
  paste(head(ux, max_n), collapse = ";")
}

classify_metadata <- function(col) {
  low <- tolower(col)
  if (grepl("donor|patient", low)) return("donor/patient identifier")
  if (grepl("sample|orig.ident", low)) return("sample identifier")
  if (grepl("disease|condition|status", low)) return("disease/status field")
  if (grepl("cell.*type|annotation|compartment|broad", low)) return("cell type or compartment annotation")
  if (grepl("cluster|seurat_clusters|res", low)) return("cluster assignment")
  if (grepl("ncount", low)) return("UMI/count QC metric")
  if (grepl("nfeature", low)) return("feature QC metric")
  if (grepl("percent.mt|mito", low)) return("mitochondrial fraction QC metric")
  if (grepl("percent.ribo|ribo", low)) return("ribosomal fraction QC metric")
  "metadata"
}

needed_for_phase2 <- function(col) {
  low <- tolower(col)
  if (grepl("donor|patient|sample|orig.ident|cell.*type|annotation|compartment|broad|cluster|ncount|nfeature|percent.mt|percent.ribo|condition|status", low)) "yes" else "no"
}

pick_col <- function(cols, patterns) {
  for (pat in patterns) {
    hit <- cols[grepl(pat, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  NA_character_
}

obj <- readRDS(canonical_object)
meta <- obj[[]]
features <- rownames(obj)
meta_cols <- colnames(meta)

assays <- names(obj@assays)
default_assay <- DefaultAssay(obj)
reductions <- names(obj@reductions)
idents <- as.character(Idents(obj))

donor_col <- pick_col(meta_cols, c("^donor$", "donor", "patient"))
sample_col <- pick_col(meta_cols, c("^sample_id$", "^sample$", "orig.ident", "sample"))
celltype_col <- pick_col(meta_cols, c("broad_cell_type", "broad.*compartment", "compartment", "cell_type", "celltype", "annotation"))
cluster_col <- pick_col(meta_cols, c("seurat_clusters", "^cluster$", "cluster"))
disease_col <- pick_col(meta_cols, c("disease_status", "condition", "status"))

get_unique_n <- function(col) if (!is.na(col)) length(unique(meta[[col]])) else NA_integer_
get_label <- function(col) if (!is.na(col)) paste(sort(unique(as.character(meta[[col]]))), collapse = ";") else ""

audit <- data.frame(
  audit_item = c(
    "canonical_object",
    "object_class",
    "assays",
    "default_assay",
    "total_nuclei_cells",
    "total_genes_features",
    "reductions_available",
    "metadata_column_count",
    "metadata_columns",
    "identity_classes",
    "donor_field",
    "donors_detected",
    "sample_field",
    "samples_detected",
    "disease_status_field",
    "cell_type_annotation_field",
    "broad_compartment_field",
    "Loop_TAL_compartment_label_present",
    "assay_layer_used_for_expression",
    "donor_level_analysis_possible",
    "donor_compartment_summary_possible"
  ),
  value = c(
    "data/processed/gse231569_audited_seurat.rds",
    paste(class(obj), collapse = ";"),
    paste(assays, collapse = ";"),
    default_assay,
    ncol(obj),
    nrow(obj),
    ifelse(length(reductions) > 0, paste(reductions, collapse = ";"), ""),
    length(meta_cols),
    paste(meta_cols, collapse = ";"),
    collapse_values(idents),
    donor_col,
    get_unique_n(donor_col),
    sample_col,
    get_unique_n(sample_col),
    disease_col,
    celltype_col,
    celltype_col,
    if (!is.na(celltype_col)) any(grepl("Loop|TAL|Henle", as.character(meta[[celltype_col]]), ignore.case = TRUE)) else FALSE,
    default_assay,
    if (!is.na(donor_col) && get_unique_n(donor_col) >= 2) "yes" else "no",
    if (!is.na(donor_col) && !is.na(celltype_col)) "yes" else "no"
  ),
  source_file = "data/processed/gse231569_audited_seurat.rds",
  notes = c(
    "Selected as processed audited Seurat object.",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "First identity values reported only.",
    "",
    "",
    "",
    "",
    "",
    "",
    "Using the same metadata field for broad compartment audit unless a more specific field is introduced later.",
    "",
    "Expression values not inspected or scored in Step 1.",
    "Input audit only; no module scoring performed.",
    "Input audit only; no donor-level score calculation performed."
  )
)
write_tsv(audit, file.path(out_dir, "phase2_step1_scrna_object_audit.tsv"))

field_rows <- lapply(meta_cols, function(col) {
  x <- meta[[col]]
  data.frame(
    metadata_column = col,
    n_unique_values = length(unique(x)),
    example_values = collapse_values(x),
    missing_fraction = mean(is.na(x) | as.character(x) == ""),
    likely_meaning = classify_metadata(col),
    needed_for_phase2 = needed_for_phase2(col),
    notes = "",
    stringsAsFactors = FALSE
  )
})
field_audit <- do.call(rbind, field_rows)
if (length(reductions) > 0) {
  field_audit <- rbind(field_audit, data.frame(
    metadata_column = "reduction_available",
    n_unique_values = length(reductions),
    example_values = paste(reductions, collapse = ";"),
    missing_fraction = 0,
    likely_meaning = "UMAP/PCA or other dimensional reduction availability",
    needed_for_phase2 = "yes",
    notes = "Stored as Seurat reductions, not metadata columns.",
    stringsAsFactors = FALSE
  ))
}
write_tsv(field_audit, file.path(out_dir, "phase2_step1_scrna_metadata_field_audit.tsv"))

if (!is.na(donor_col) && !is.na(celltype_col)) {
  sample_values <- if (!is.na(sample_col)) as.character(meta[[sample_col]]) else "NA"
  combo <- data.frame(
    donor_id = as.character(meta[[donor_col]]),
    sample_id = sample_values,
    compartment = as.character(meta[[celltype_col]]),
    stringsAsFactors = FALSE
  )
  counts <- as.data.frame(table(combo$donor_id, combo$sample_id, combo$compartment), stringsAsFactors = FALSE)
  colnames(counts) <- c("donor_id", "sample_id", "compartment", "n_nuclei")
  counts <- counts[counts$n_nuclei > 0, ]
  donor_totals <- aggregate(n_nuclei ~ donor_id, counts, sum)
  colnames(donor_totals)[2] <- "donor_total"
  counts <- merge(counts, donor_totals, by = "donor_id")
  counts$percent_within_donor <- round(100 * counts$n_nuclei / counts$donor_total, 4)
  counts$notes <- ""
  counts <- counts[, c("donor_id", "sample_id", "compartment", "n_nuclei", "percent_within_donor", "notes")]
  write_tsv(counts, file.path(out_dir, "phase2_step1_donor_compartment_count_audit.tsv"))

  comp_summary <- aggregate(n_nuclei ~ compartment, counts, sum)
  donor_dist <- aggregate(donor_id ~ compartment, counts, function(x) paste(sort(unique(x)), collapse = ";"))
  donor_n <- aggregate(donor_id ~ compartment, counts, function(x) length(unique(x)))
  colnames(donor_n)[2] <- "n_donors_present"
  colnames(donor_dist)[2] <- "donor_distribution"
  comp_summary <- merge(comp_summary, donor_n, by = "compartment")
  comp_summary <- merge(comp_summary, donor_dist, by = "compartment")
  comp_summary$candidate_for_downstream_analysis <- ifelse(comp_summary$n_donors_present >= 2 & comp_summary$n_nuclei >= 50, "yes", "needs_review")
  comp_summary$notes <- ifelse(grepl("Loop|TAL|Henle", comp_summary$compartment, ignore.case = TRUE), "Loop/TAL-like compartment identified.", "")
  colnames(comp_summary)[colnames(comp_summary) == "n_nuclei"] <- "n_nuclei_total"
  write_tsv(comp_summary[, c("compartment", "n_nuclei_total", "n_donors_present", "donor_distribution", "candidate_for_downstream_analysis", "notes")], file.path(out_dir, "phase2_step1_compartment_summary.tsv"))
} else {
  write_tsv(data.frame(donor_id = character(), sample_id = character(), compartment = character(), n_nuclei = integer(), percent_within_donor = numeric(), notes = character()), file.path(out_dir, "phase2_step1_donor_compartment_count_audit.tsv"))
  write_tsv(data.frame(compartment = character(), n_nuclei_total = integer(), n_donors_present = integer(), donor_distribution = character(), candidate_for_downstream_analysis = character(), notes = character()), file.path(out_dir, "phase2_step1_compartment_summary.tsv"))
}

markers <- c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16", "CLDN14", "CASR", "PKD2", "HIBADH")
feature_upper <- toupper(features)
marker_rows <- lapply(markers, function(g) {
  idx <- which(feature_upper == toupper(g))
  data.frame(
    gene = g,
    present_in_scrna = length(idx) > 0,
    feature_id_if_different = if (length(idx) > 0 && features[idx[1]] != g) features[idx[1]] else "",
    notes = "Presence check only; no expression interpretation.",
    stringsAsFactors = FALSE
  )
})
write_tsv(do.call(rbind, marker_rows), file.path(out_dir, "phase2_step1_loop_tal_marker_presence.tsv"))

modules <- data.frame(
  module_name = c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni", "MAGMA_FDR05", "MAGMA_suggestive_p1e4"),
  module_file = c(
    "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
  ),
  stringsAsFactors = FALSE
)
module_rows <- lapply(seq_len(nrow(modules)), function(i) {
  genes <- readLines(file.path(root, modules$module_file[i]), warn = FALSE)
  genes <- genes[nzchar(genes)]
  present <- genes[toupper(genes) %in% feature_upper]
  missing <- setdiff(genes, present)
  non_symbol <- genes[grepl("^[0-9]+$", genes)]
  data.frame(
    module_name = modules$module_name[i],
    module_file = modules$module_file[i],
    n_module_genes = length(genes),
    n_present_in_scrna = length(present),
    n_missing_from_scrna = length(missing),
    present_genes = paste(present, collapse = ";"),
    missing_genes = paste(missing, collapse = ";"),
    non_symbol_gene_ids = paste(non_symbol, collapse = ";"),
    suitable_for_scoring = ifelse(length(present) >= 5, "yes_input_suitable_pending_detected_gene_filtering", "needs_review"),
    notes = "Canonical module definition preserved; no missing genes dropped silently.",
    stringsAsFactors = FALSE
  )
})
write_tsv(do.call(rbind, module_rows), file.path(out_dir, "phase2_step1_magma_module_scrna_mapping.tsv"))
