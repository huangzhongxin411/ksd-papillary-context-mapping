suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

paths <- list.files(c("external/twas", "results/phase2b_twas_ready_v0.1", "data/processed"), recursive = TRUE, full.names = TRUE)
non_test <- paths[!grepl("/tests?/", paths)]
fusion_software <- non_test[grepl("FUSION\\.assoc_test\\.R$|fusion/.+\\.sh$|fusion/.+\\.R$", non_test, ignore.case = TRUE)]
fusion_weights <- non_test[grepl("\\.pos$|\\.wgt\\.RDat$|/weights/", non_test, ignore.case = TRUE)]
ldref <- non_test[grepl("LDREF|1000G|1000G_EUR|LD_?REF", non_test, ignore.case = TRUE)]
predixcan_software <- non_test[grepl("SPrediXcan.py$|S?MultiXcan.py$", non_test)]
predixcan_models <- non_test[grepl("\\.db$", non_test)]
predixcan_cov <- non_test[grepl("covariance|covariances", non_test, ignore.case = TRUE) & grepl("\\.txt\\.gz$|\\.gz$", non_test)]
fusion_input <- c("data/processed/twas_input/ksd_2025_for_fusion.tsv",
                  "results/phase2b_twas_ready_v0.1/data/ksd_2025_for_fusion.tsv")

status <- data.table(
  resource = c("FUSION software", "FUSION weights", "FUSION LD reference", "FUSION input",
               "S-PrediXcan software", "S-PrediXcan model db", "S-PrediXcan covariance"),
  status = c(ifelse(length(fusion_software) > 0, "candidate_found", "missing"),
             ifelse(length(fusion_weights) > 0, "candidate_found", "missing"),
             ifelse(length(ldref) > 0, "candidate_found", "missing"),
             ifelse(any(file.exists(fusion_input)), "available", "missing"),
             ifelse(length(predixcan_software) > 0, "candidate_found", "missing"),
             ifelse(length(predixcan_models) > 0, "candidate_found", "missing"),
             ifelse(length(predixcan_cov) > 0, "candidate_found", "missing")),
  path_or_note = c(paste(head(fusion_software, 5), collapse = "; "),
                   paste(head(fusion_weights, 5), collapse = "; "),
                   paste(head(ldref, 5), collapse = "; "),
                   paste(fusion_input[file.exists(fusion_input)], collapse = "; "),
                   paste(head(predixcan_software, 5), collapse = "; "),
                   paste(head(predixcan_models, 5), collapse = "; "),
                   paste(head(predixcan_cov, 5), collapse = "; "))
)
status[path_or_note == "", path_or_note := "not found in current workspace"]
status[, usable_for_main_analysis := status %in% c("available", "candidate_found")]
status[resource %in% c("FUSION weights", "FUSION LD reference", "S-PrediXcan model db", "S-PrediXcan covariance") &
         status != "candidate_found", usable_for_main_analysis := FALSE]
fwrite(status, "results/twas/twas_resource_status.tsv", sep = "\t")

decision <- if (all(status[resource %in% c("FUSION software", "FUSION weights", "FUSION LD reference", "FUSION input"), status] %in% c("available", "candidate_found"))) {
  "fusion_candidate_resources_found_review_before_run"
} else {
  "twas_not_ready_resource_missing"
}
writeLines(c(
  "# TWAS Resource Blocker Memo v0.1",
  "",
  paste0("Decision: ", decision),
  "",
  "This audit does not run FUSION or S-PrediXcan. It checks whether local software, model weights, LD/covariance files and formatted GWAS inputs are present.",
  "",
  "TWAS results should only be added to the manuscript after real GTEx/PredictDB model resources and matched LD/covariance files are confirmed and smoke tests pass."
), "docs/twas_resource_blocker_memo.md", useBytes = TRUE)

message("wrote TWAS resource status")
