suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- "config/gse73680_sample_metadata_curated.tsv"
expr_path <- "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz"
stopifnot(file.exists(meta_path), file.exists(expr_path))

meta <- fread(meta_path)
expr_header <- names(fread(expr_path, nrows = 0))
matrix_samples <- setdiff(expr_header, "gene")

meta <- meta[sample_id %in% matrix_samples]
meta[, include_in_analysis := include_in_analysis %in% TRUE]
included <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]

patient_structure <- meta[, .(
  n_samples = .N,
  n_control_or_adjacent = sum(group_curated == "control_or_adjacent" & include_in_analysis),
  n_plaque_or_stone_papilla = sum(group_curated == "plaque_or_stone_papilla" & include_in_analysis),
  has_both_groups = any(group_curated == "control_or_adjacent" & include_in_analysis) &&
    any(group_curated == "plaque_or_stone_papilla" & include_in_analysis),
  sample_ids = paste(sample_id, collapse = ";"),
  group_pattern = paste(sort(unique(group_curated[include_in_analysis])), collapse = ";"),
  n_included = sum(include_in_analysis),
  n_unclear_or_excluded = sum(!include_in_analysis)
), by = patient_id]

patient_structure[, recommended_design := fcase(
  n_included == 0, "exclude_unclear",
  has_both_groups == TRUE, "paired_patient_design",
  n_included > 1, "unpaired_patient_blocked_design",
  n_included == 1, "sample_level_exploratory",
  default = "exclude_unclear"
)]
setorder(patient_structure, -has_both_groups, -n_samples, patient_id)
fwrite(patient_structure[, .(patient_id, n_samples, n_control_or_adjacent, n_plaque_or_stone_papilla,
                             has_both_groups, sample_ids, group_pattern, recommended_design)],
       file.path(table_dir, "gse73680_patient_sample_structure.tsv"), sep = "\t")

n_samples_total <- nrow(meta)
n_samples_included <- nrow(included)
n_patients_total <- uniqueN(meta$patient_id)
n_patients_included <- uniqueN(included$patient_id)
n_groups <- uniqueN(included$group_curated)
group_counts <- included[, .N, by = group_curated][, paste(group_curated, N, sep = "=")] |> paste(collapse = ";")
patient_group_counts <- included[, .(n_samples = .N, n_groups = uniqueN(group_curated)), by = patient_id]
n_patients_with_repeated_samples <- sum(patient_group_counts$n_samples > 1)
n_patients_with_both_groups <- sum(patient_group_counts$n_groups > 1)
included_patient_by_group <- included[, uniqueN(patient_id), by = group_curated]
min_patients_per_group <- min(included_patient_by_group$V1)

recommended_model <- if (n_patients_with_repeated_samples > 0) "limma_duplicateCorrelation" else "limma_sample_level"
paired_sensitivity_required <- n_patients_with_both_groups >= 3
downstream_analysis_grade <- if (min_patients_per_group < 3) "exploratory_only" else "primary_with_patient_blocking"

design <- data.table(
  design_item = c("n_samples_total", "n_samples_included", "n_patients_total", "n_patients_included",
                  "n_groups", "group_counts", "n_patients_with_repeated_samples",
                  "n_patients_with_both_groups", "min_patients_per_group", "recommended_model",
                  "paired_sensitivity_required", "downstream_analysis_grade"),
  value = as.character(c(n_samples_total, n_samples_included, n_patients_total, n_patients_included,
                         n_groups, group_counts, n_patients_with_repeated_samples,
                         n_patients_with_both_groups, min_patients_per_group, recommended_model,
                         paired_sensitivity_required, downstream_analysis_grade)),
  decision = c(rep("observed", 9), "use_for_primary_model",
               ifelse(paired_sensitivity_required, "required", "not_required"),
               "claim_strength_boundary"),
  notes = c(
    "All samples in curated metadata that are present in the expression matrix.",
    "Samples marked include_in_analysis and assigned to the two analysis groups.",
    "Unique patient IDs in curated metadata.",
    "Unique patient IDs among included samples.",
    "Number of included groups.",
    "Included sample counts by group.",
    "Repeated included samples imply non-independent observations.",
    "Patients contributing both control/adjacent and plaque/stone papilla samples.",
    "Minimum included patient count across groups.",
    "Primary model selected from patient structure.",
    "Paired sensitivity is required when at least three patients have both groups.",
    "Disease-context claims must respect sample structure and metadata curation limits."
  )
)
fwrite(design, file.path(table_dir, "gse73680_analysis_design.tsv"), sep = "\t")

expr <- fread(expr_path)
expr <- expr[!is.na(gene) & gene != ""]
patient_groups <- included[, .(sample_ids = list(intersect(sample_id, names(expr)))), by = .(patient_id, group_curated)]
patient_groups <- patient_groups[lengths(sample_ids) > 0]
patient_matrix <- data.table(gene = expr$gene)
for (i in seq_len(nrow(patient_groups))) {
  col_name <- paste(patient_groups$patient_id[i], patient_groups$group_curated[i], sep = "__")
  ids <- patient_groups$sample_ids[[i]]
  vals <- as.matrix(expr[, ..ids])
  mode(vals) <- "numeric"
  patient_matrix[, (col_name) := rowMeans(vals, na.rm = TRUE)]
}
fwrite(patient_matrix, "data/processed/gse73680/gse73680_patient_level_expression_matrix.tsv.gz", sep = "\t")
message("wrote GSE73680 analysis design")
