suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
processed_dir <- "data/processed/gse73680"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

meta_path <- if (file.exists(file.path(table_dir, "gse73680_metadata_audit.tsv"))) file.path(table_dir, "gse73680_metadata_audit.tsv") else "results/tables/gse73680_metadata_audit.tsv"
sample_sheet_path <- file.path(table_dir, "gse73680_sample_sheet.tsv")
meta <- if (file.exists(meta_path)) fread(meta_path) else data.table()
sample_sheet <- if (file.exists(sample_sheet_path)) fread(sample_sheet_path) else data.table(sample_id = character())
matrix_samples <- sample_sheet$sample_id

if (!nrow(meta)) {
  curated <- data.table(sample_id = matrix_samples, patient_id = NA_character_, geo_accession = matrix_samples,
                        raw_label = NA_character_, group_curated = "unclear", disease_status = NA_character_,
                        plaque_status = NA_character_, stone_status = NA_character_, tissue_site = NA_character_,
                        sample_type = NA_character_, batch = NA_character_, include_in_analysis = FALSE,
                        exclude_reason = "metadata_missing", notes = "No metadata audit available.")
} else {
  curated <- copy(meta)
  if (!"raw_label" %in% names(curated)) curated[, raw_label := if ("group_raw" %in% names(curated)) group_raw else sample_id]
  curated[, patient_id := fifelse(grepl("_(ST[0-9]+)_", raw_label),
                                  sub("^.*_(ST[0-9]+)_.*$", "\\1", raw_label),
                                  NA_character_)]
  curated[, group_curated := fcase(
    disease_status == "stone_plaque_context", "plaque_or_stone_papilla",
    disease_status == "normal_context", "control_or_adjacent",
    disease_status == "control_or_other_context", "unclear",
    default = "unclear"
  )]
  curated[, tissue_site := fifelse(group_curated == "plaque_or_stone_papilla", "papillary_plaque_or_stone_context",
                            fifelse(group_curated == "control_or_adjacent", "normal_or_adjacent_papilla", NA_character_))]
  curated[, batch := NA_character_]
  curated[, include_in_analysis := group_curated != "unclear" & sample_id %in% matrix_samples]
  curated[, exclude_reason := fifelse(include_in_analysis, NA_character_,
                               fifelse(!(sample_id %in% matrix_samples), "sample_missing_in_expression_matrix", "group_unclear"))]
  curated[, notes := "Filename-derived curation; manual verification required before differential analysis."]
  needed <- c("sample_id", "patient_id", "geo_accession", "raw_label", "group_curated", "disease_status",
              "plaque_status", "stone_status", "tissue_site", "sample_type", "batch",
              "include_in_analysis", "exclude_reason", "notes")
  for (nm in setdiff(needed, names(curated))) curated[, (nm) := NA_character_]
  curated <- curated[, ..needed]
}
fwrite(curated, "config/gse73680_sample_metadata_curated.tsv", sep = "\t")
n_in_matrix <- length(matrix_samples)
group_counts <- if (nrow(curated)) paste(curated[include_in_analysis == TRUE, .N, by = group_curated][, paste(group_curated, N, sep = "=")], collapse = ";") else ""
n_groups <- if (nrow(curated)) uniqueN(curated[include_in_analysis == TRUE]$group_curated) else 0L
n_included <- sum(curated$include_in_analysis %in% TRUE)
n_unclear <- sum(curated$group_curated == "unclear", na.rm = TRUE)
status <- if (n_included >= 6 && n_groups >= 2 && n_unclear / max(1, n_in_matrix) <= 0.30) "pass" else if (n_included >= 4 && n_groups >= 2) "warning" else "fail"
fwrite(data.table(n_samples_in_matrix = n_in_matrix, n_samples_with_metadata = nrow(curated),
                  n_samples_curated = nrow(curated), n_samples_included = n_included,
                  n_groups = n_groups, group_counts = group_counts, n_unclear = n_unclear,
                  n_excluded = sum(!(curated$include_in_analysis %in% TRUE)),
                  status = status, main_blocker = if (status == "fail") "expression_matrix_or_reliable_group_labels_missing" else NA_character_,
                  notes = "No differential analysis is allowed unless metadata status is pass or explicitly conditional."),
       file.path(table_dir, "gse73680_metadata_curated_audit.tsv"), sep = "\t")
message("wrote GSE73680 curated metadata outputs")
