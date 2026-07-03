suppressPackageStartupMessages(library(data.table))

tissues <- c(
  "kidney_cortex",
  "whole_blood",
  "artery_aorta",
  "artery_tibial",
  "adipose_subcutaneous",
  "liver",
  "colon_transverse",
  "small_intestine_terminal_ileum"
)

first_file <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths)) paths[[1]] else NA_character_
}

find_one <- function(path, pattern) {
  if (!dir.exists(path)) return(NA_character_)
  x <- list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(x)) x[[1]] else NA_character_
}

fusion_software <- "external/twas/fusion/software/FUSION.assoc_test.R"
fusion_ld <- "external/twas/fusion/ref_ld/1000G.EUR."
fusion_rows <- rbindlist(lapply(tissues, function(tissue) {
  weight_dir <- file.path("external/twas/fusion/weights", tissue)
  pos <- find_one(weight_dir, "\\.pos$")
  n_weights <- if (dir.exists(weight_dir)) {
    length(list.files(weight_dir, pattern = "\\.RDat$", recursive = TRUE, full.names = TRUE))
  } else 0L
  data.table(
    method = "FUSION",
    resource_type = "expression_weights",
    tissue = tissue,
    source = "pending",
    version = "pending",
    genome_build = "expected_GRCh37_hg19",
    model_file = if (n_weights > 0) weight_dir else NA_character_,
    weight_pos_file = pos,
    covariance_file = NA_character_,
    ld_reference = fusion_ld,
    status = if (file.exists(fusion_software) && !is.na(pos) && n_weights > 0) "ready" else "missing",
    notes = sprintf("FUSION software=%s; n_RDat=%s", ifelse(file.exists(fusion_software), "present", "missing"), n_weights)
  )
}))

predixcan_software <- "external/twas/predixcan/software/SPrediXcan.py"
predixcan_rows <- rbindlist(lapply(tissues, function(tissue) {
  model <- find_one("external/twas/predixcan/models", paste0(tissue, ".*\\.db$"))
  cov <- find_one("external/twas/predixcan/covariance", tissue)
  data.table(
    method = "S-PrediXcan",
    resource_type = "predictdb_model",
    tissue = tissue,
    source = "pending",
    version = "pending",
    genome_build = "expected_GRCh37_hg19",
    model_file = model,
    weight_pos_file = NA_character_,
    covariance_file = cov,
    ld_reference = NA_character_,
    status = if (file.exists(predixcan_software) && !is.na(model) && !is.na(cov)) "ready" else "missing",
    notes = sprintf("SPrediXcan.py=%s", ifelse(file.exists(predixcan_software), "present", "missing"))
  )
}))

manifest <- rbindlist(list(fusion_rows, predixcan_rows), fill = TRUE)
fwrite(manifest, "results/tables/twas_resource_manifest.tsv", sep = "\t")
fwrite(manifest, "external/twas/manifests/twas_resource_manifest.tsv", sep = "\t")

message("wrote\tresults/tables/twas_resource_manifest.tsv")
message("wrote\texternal/twas/manifests/twas_resource_manifest.tsv")
