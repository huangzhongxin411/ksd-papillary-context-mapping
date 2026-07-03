suppressPackageStartupMessages({
  library(data.table)
  library(tools)
})

tissues <- data.table(
  tissue = c(
    "kidney_cortex",
    "whole_blood",
    "artery_aorta",
    "artery_tibial",
    "adipose_subcutaneous",
    "liver",
    "colon_transverse",
    "small_intestine_terminal_ileum"
  ),
  fusion_alias = c(
    "GTExv8.Kidney_Cortex",
    "GTExv8.Whole_Blood",
    "GTExv8.Artery_Aorta",
    "GTExv8.Artery_Tibial",
    "GTExv8.Adipose_Subcutaneous",
    "GTExv8.Liver",
    "GTExv8.Colon_Transverse",
    "GTExv8.Small_Intestine_Terminal_Ileum"
  ),
  predixcan_alias = c(
    "Kidney_Cortex",
    "Whole_Blood",
    "Artery_Aorta",
    "Artery_Tibial",
    "Adipose_Subcutaneous",
    "Liver",
    "Colon_Transverse",
    "Small_Intestine_Terminal_Ileum"
  )
)

file_status <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) "missing" else "downloaded_unchecked"
}

file_size <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  as.character(file.info(path)$size)
}

file_md5 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NA_character_)
  unname(tools::md5sum(path))
}

blocking_level <- function(method, resource_type, tissue) {
  if (method == "FUSION" && resource_type %in% c("software", "ld_reference")) return("critical")
  if (method == "FUSION" && tissue == "kidney_cortex") return("critical")
  if (method == "FUSION" && tissue %in% c("whole_blood", "artery_aorta", "artery_tibial")) return("important")
  if (method == "S-PrediXcan" && resource_type == "software") return("critical")
  if (method == "S-MultiXcan" && resource_type == "software") return("optional")
  if (method == "S-PrediXcan" && tissue == "kidney_cortex") return("critical")
  if (method == "S-PrediXcan" && tissue %in% c("whole_blood", "artery_aorta", "artery_tibial")) return("important")
  "optional"
}

action_decision <- function(status, block) {
  if (status %in% c("checked_ok", "downloaded_unchecked")) return("run_smoke_test")
  if (block == "critical") return("retry")
  if (block == "important") return("manual_download")
  "freeze_missing"
}

decorate <- function(dt) {
  dt[, download_attempted := status != "missing"]
  dt[, download_method := fifelse(download_attempted, "curl_or_git_clone", NA_character_)]
  dt[, download_status := fifelse(download_attempted, status, "not_attempted")]
  dt[, blocking_level := mapply(blocking_level, method, resource_type, tissue)]
  dt[, action_decision := mapply(action_decision, status, blocking_level)]
  dt[]
}

find_one <- function(path, pattern) {
  if (!dir.exists(path)) return(NA_character_)
  x <- list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(x)) x[[1]] else NA_character_
}

source_fusion <- "https://gusevlab.org/projects/fusion/"
source_predictdb <- "https://predictdb.org/post/2019/12/11/gtex-v8-model-release/"
source_metaxcan <- "https://github.com/hakyimlab/MetaXcan"

rows <- list(
  data.table(
    method = "FUSION",
    resource_type = "software",
    tissue = "all",
    tissue_alias = "all",
    expected_file = "FUSION.assoc_test.R",
    local_path = "external/twas/fusion/software/FUSION.assoc_test.R",
    source_page = source_fusion,
    genome_build = "tool",
    version = "pending",
    status = file_status("external/twas/fusion/software/FUSION.assoc_test.R"),
    file_size = file_size("external/twas/fusion/software/FUSION.assoc_test.R"),
    md5 = file_md5("external/twas/fusion/software/FUSION.assoc_test.R"),
    notes = "FUSION association-test script."
  ),
  data.table(
    method = "FUSION",
    resource_type = "ld_reference",
    tissue = "all",
    tissue_alias = "1000G.EUR",
    expected_file = "1000G.EUR.{1..22}.{bed,bim,fam}",
    local_path = "external/twas/fusion/ref_ld/1000G.EUR.",
    source_page = source_fusion,
    genome_build = "expected_GRCh37_hg19",
    version = "pending",
    status = if (file.exists("external/twas/fusion/ref_ld/1000G.EUR.1.bim") ||
      file.exists("external/twas/fusion/ref_ld/1000G.EUR.chr1.bim")) "downloaded_unchecked" else "missing",
    file_size = NA_character_,
    md5 = NA_character_,
    notes = "FUSION LD reference prefix; require chr-level PLINK files."
  ),
  data.table(
    method = "S-PrediXcan",
    resource_type = "software",
    tissue = "all",
    tissue_alias = "all",
    expected_file = "SPrediXcan.py",
    local_path = ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py",
      "external/twas/predixcan/software/SPrediXcan.py"
    ),
    source_page = source_metaxcan,
    genome_build = "tool",
    version = "pending",
    status = file_status(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py",
      "external/twas/predixcan/software/SPrediXcan.py"
    )),
    file_size = file_size(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py",
      "external/twas/predixcan/software/SPrediXcan.py"
    )),
    md5 = file_md5(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py",
      "external/twas/predixcan/software/SPrediXcan.py"
    )),
    notes = "MetaXcan S-PrediXcan script."
  ),
  data.table(
    method = "S-MultiXcan",
    resource_type = "software",
    tissue = "all",
    tissue_alias = "all",
    expected_file = "SMulTiXcan.py",
    local_path = ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py",
      "external/twas/predixcan/software/SMulTiXcan.py"
    ),
    source_page = source_metaxcan,
    genome_build = "tool",
    version = "pending",
    status = file_status(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py",
      "external/twas/predixcan/software/SMulTiXcan.py"
    )),
    file_size = file_size(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py",
      "external/twas/predixcan/software/SMulTiXcan.py"
    )),
    md5 = file_md5(ifelse(
      file.exists("external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py"),
      "external/twas/predixcan/software/MetaXcan/software/SMulTiXcan.py",
      "external/twas/predixcan/software/SMulTiXcan.py"
    )),
    notes = "MetaXcan S-MultiXcan script."
  )
)

for (i in seq_len(nrow(tissues))) {
  tissue <- tissues$tissue[[i]]
  fusion_alias <- tissues$fusion_alias[[i]]
  predixcan_alias <- tissues$predixcan_alias[[i]]
  fusion_dir <- file.path("external/twas/fusion/weights", tissue)
  fusion_pos <- find_one(fusion_dir, "\\.pos$")
  fusion_weight <- find_one(fusion_dir, "\\.RDat$")
  pred_model <- find_one("external/twas/predixcan/models", paste0(predixcan_alias, ".*\\.db$"))
  pred_cov <- find_one("external/twas/predixcan/covariance", predixcan_alias)

  rows[[length(rows) + 1]] <- data.table(
    method = "FUSION",
    resource_type = "weights_pos",
    tissue = tissue,
    tissue_alias = fusion_alias,
    expected_file = paste0(fusion_alias, ".pos"),
    local_path = ifelse(is.na(fusion_pos), file.path(fusion_dir, paste0(fusion_alias, ".pos")), fusion_pos),
    source_page = source_fusion,
    genome_build = "expected_GRCh37_hg19",
    version = "GTEx_v8_expected",
    status = ifelse(is.na(fusion_pos), "missing", file_status(fusion_pos)),
    file_size = ifelse(is.na(fusion_pos), NA_character_, file_size(fusion_pos)),
    md5 = ifelse(is.na(fusion_pos), NA_character_, file_md5(fusion_pos)),
    notes = "FUSION weights position file."
  )
  rows[[length(rows) + 1]] <- data.table(
    method = "FUSION",
    resource_type = "weights_rdat",
    tissue = tissue,
    tissue_alias = fusion_alias,
    expected_file = paste0(fusion_alias, "/*.wgt.RDat"),
    local_path = ifelse(is.na(fusion_weight), fusion_dir, fusion_weight),
    source_page = source_fusion,
    genome_build = "expected_GRCh37_hg19",
    version = "GTEx_v8_expected",
    status = ifelse(is.na(fusion_weight), "missing", file_status(fusion_weight)),
    file_size = ifelse(is.na(fusion_weight), NA_character_, file_size(fusion_weight)),
    md5 = ifelse(is.na(fusion_weight), NA_character_, file_md5(fusion_weight)),
    notes = "At least one RDat file is required; full tissue run requires all referenced weights."
  )
  rows[[length(rows) + 1]] <- data.table(
    method = "S-PrediXcan",
    resource_type = "predictdb_model",
    tissue = tissue,
    tissue_alias = predixcan_alias,
    expected_file = paste0(predixcan_alias, "*.db"),
    local_path = ifelse(is.na(pred_model), "external/twas/predixcan/models/", pred_model),
    source_page = source_predictdb,
    genome_build = "expected_GRCh37_hg19",
    version = "GTEx_v8_expected",
    status = ifelse(is.na(pred_model), "missing", file_status(pred_model)),
    file_size = ifelse(is.na(pred_model), NA_character_, file_size(pred_model)),
    md5 = ifelse(is.na(pred_model), NA_character_, file_md5(pred_model)),
    notes = "PredictDB model database."
  )
  rows[[length(rows) + 1]] <- data.table(
    method = "S-PrediXcan",
    resource_type = "covariance",
    tissue = tissue,
    tissue_alias = predixcan_alias,
    expected_file = paste0(predixcan_alias, "*.txt.gz"),
    local_path = ifelse(is.na(pred_cov), "external/twas/predixcan/covariance/", pred_cov),
    source_page = source_predictdb,
    genome_build = "expected_GRCh37_hg19",
    version = "GTEx_v8_expected",
    status = ifelse(is.na(pred_cov), "missing", file_status(pred_cov)),
    file_size = ifelse(is.na(pred_cov), NA_character_, file_size(pred_cov)),
    md5 = ifelse(is.na(pred_cov), NA_character_, file_md5(pred_cov)),
    notes = "Summary-based S-PrediXcan covariance file."
  )
}

manifest <- rbindlist(rows, fill = TRUE)
manifest <- decorate(manifest)
fwrite(manifest, "results/tables/twas_resource_manifest_v2.tsv", sep = "\t")
fwrite(manifest, "external/twas/manifests/twas_resource_manifest_v2.tsv", sep = "\t")
message("wrote\tresults/tables/twas_resource_manifest_v2.tsv")
message("wrote\texternal/twas/manifests/twas_resource_manifest_v2.tsv")
