suppressPackageStartupMessages(library(data.table))

status_path <- "results/coloc/coloc_per_locus_status.tsv"
summary_path <- "results/tables/coloc_summary_v0.1.tsv"

if (!file.exists(status_path)) {
  fwrite(data.table(
    status = "not_run",
    n_targets = 0L,
    n_ready = 0L,
    n_with_coloc_result = 0L,
    n_high_pph4 = 0L,
    notes = "No per-locus coloc status file found. Run 02_run_coloc_per_locus.R after preparing inputs."
  ), summary_path, sep = "\t")
  message("No coloc status file; wrote not_run summary.")
  quit(save = "no", status = 0)
}

status <- fread(status_path)
summary <- data.table(
  status = if (any(status$coloc_status == "completed")) "partial_or_completed" else "not_ready",
  n_targets = nrow(status),
  n_ready = sum(status$coloc_status %in% c("ready", "completed"), na.rm = TRUE),
  n_with_coloc_result = sum(!is.na(status$coloc_pp_h4), na.rm = TRUE),
  n_high_pph4 = sum(status$coloc_pp_h4 >= 0.8, na.rm = TRUE),
  notes = if (any(status$coloc_status == "completed")) {
    "Coloc results available for at least one target; inspect per-locus table before interpretation."
  } else {
    "No coloc statistics available yet because eQTL resources/harmonized per-locus inputs are missing."
  }
)

fwrite(summary, summary_path, sep = "\t")
message("wrote\t", summary_path)
