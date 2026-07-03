suppressPackageStartupMessages({
  library(data.table)
})

audit_path <- "results/tables/gse231569_cluster_marker_assignment_audit.tsv"
cell_counts_path <- "results/tables/gse231569_cell_counts.tsv"
donor_counts_path <- "results/tables/gse231569_donor_cell_counts.tsv"
out_path <- "config/gse231569_cluster_to_celltype_audited.tsv"

dir.create("config", recursive = TRUE, showWarnings = FALSE)

audit <- fread(audit_path)
donor_counts <- fread(donor_counts_path)

donor_summary <- donor_counts[
  ,
  .(n_donors = uniqueN(donor_id), n_cells = sum(n_cells)),
  by = broad_cell_type
]

cfg <- audit[
  ,
  .(
    cluster = as.character(cluster),
    broad_cell_type = fifelse(
      current_phase1_cell_type == "Pericyte_smooth_muscle",
      "Perivascular_mural_like",
      current_phase1_cell_type
    ),
    fine_cell_type = fifelse(
      current_phase1_cell_type == "Pericyte_smooth_muscle",
      "Perivascular/mural-like cells",
      current_phase1_cell_type
    ),
    phase1_cell_type = current_phase1_cell_type,
    annotation_confidence = fifelse(
      current_phase1_cell_type == "Loop_of_Henle_TAL",
      "high",
      fifelse(
        current_phase1_cell_type == "Pericyte_smooth_muscle",
        "low_or_exploratory",
        fifelse(
          immune_review_flag == "immune_marker_signal_requires_manual_review",
          "review",
          confidence
        )
      )
    ),
    main_positive_markers = marker_genes_present,
    immune_review_flag = immune_review_flag == "immune_marker_signal_requires_manual_review",
    mural_review_flag = current_phase1_cell_type == "Pericyte_smooth_muscle",
    top_marker_set = candidate_cell_type,
    top_immune_marker_set,
    top_immune_mean_marker_expression,
    top_immune_pct_cells_with_any_marker,
    n_cluster_cells = n_cells,
    notes = fifelse(
      current_phase1_cell_type == "Pericyte_smooth_muscle",
      "Renamed to Perivascular_mural_like for conservative downstream interpretation; small cell group.",
      fifelse(
        immune_review_flag == "immune_marker_signal_requires_manual_review",
        "Primary label retained, but immune marker signal requires manual review.",
        audit_flag
      )
    )
  )
]

cfg <- merge(
  cfg,
  donor_summary,
  by.x = "phase1_cell_type",
  by.y = "broad_cell_type",
  all.x = TRUE
)
setnames(cfg, c("n_donors", "n_cells"), c("n_celltype_donors", "n_celltype_cells"))

setcolorder(cfg, c(
  "cluster",
  "broad_cell_type",
  "fine_cell_type",
  "annotation_confidence",
  "main_positive_markers",
  "immune_review_flag",
  "mural_review_flag",
  "n_cluster_cells",
  "n_celltype_donors",
  "n_celltype_cells",
  "phase1_cell_type",
  "top_marker_set",
  "top_immune_marker_set",
  "top_immune_mean_marker_expression",
  "top_immune_pct_cells_with_any_marker",
  "notes"
))

cfg[, cluster_order := as.integer(cluster)]
setorder(cfg, cluster_order)
cfg[, cluster_order := NULL]
fwrite(cfg, out_path, sep = "\t")
cat("wrote\t", out_path, "\n", sep = "")
