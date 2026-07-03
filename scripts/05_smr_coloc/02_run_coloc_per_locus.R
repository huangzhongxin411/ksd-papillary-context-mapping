suppressPackageStartupMessages(library(data.table))

pilot_path <- "results/tables/coloc_pilot_targets_v0.1.tsv"
readiness_path <- "results/tables/coloc_input_readiness.tsv"
eqtl_status_path <- "results/tables/eqtl_for_coloc_preparation_status.tsv"
out_dir <- "results/coloc"
out_path <- file.path(out_dir, "coloc_per_locus_status.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pilot_path) || !file.exists(readiness_path)) {
  stop("Missing coloc pilot/readiness tables. Run scripts/05_smr_coloc/prepare_coloc_execution_plan.R first.")
}

pilot <- fread(pilot_path)
readiness <- fread(readiness_path)

if (!file.exists(eqtl_status_path)) {
  fwrite(data.table(
    gene = unique(pilot$gene),
    coloc_status = "not_ready",
    blocking_reason = "eQTL preparation status table missing; run 01_prepare_eqtl_for_coloc.R first.",
    coloc_pp_h4 = NA_real_,
    notes = "No coloc statistics were computed."
  ), out_path, sep = "\t")
  message("Missing eQTL preparation status; wrote not_ready coloc status.")
  quit(save = "no", status = 0)
}

eqtl_status <- fread(eqtl_status_path)
if (!nrow(eqtl_status) || any(eqtl_status$status != "ready")) {
  status_note <- if ("notes" %in% names(eqtl_status)) eqtl_status$notes[1] else "eQTL input is not ready."
  out <- merge(
    pilot[, .(gene, pilot_priority, priority_class, locus_id, chr, locus_start, locus_end)],
    readiness[, .(gene, eqtl_resource_exists, variant_id_matchable, build_matchable, blocking_reason)],
    by = "gene",
    all.x = TRUE
  )
  out[, `:=`(
    coloc_status = "not_ready",
    coloc_pp_h4 = NA_real_,
    coloc_pp_h3 = NA_real_,
    coloc_input_nsnps = NA_integer_,
    coloc_result_file = NA_character_,
    notes = paste("No coloc statistics were computed.", status_note)
  )]
  fwrite(out, out_path, sep = "\t")
  message("eQTL resources are not ready; wrote not_ready coloc status.")
  quit(save = "no", status = 0)
}

stop("Coloc execution is intentionally disabled until a resource-specific parser creates ready per-locus eQTL/GWAS inputs.")
