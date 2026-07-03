suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
read_status <- function(path, field = "status", default = "fail") {
  if (!file.exists(path)) return(default)
  x <- fread(path)
  if (!nrow(x) || !(field %in% names(x))) return(default)
  as.character(x[[field]][1])
}
read_bool <- function(path, expr_fun, default = FALSE) {
  if (!file.exists(path)) return(default)
  x <- fread(path)
  if (!nrow(x)) return(default)
  isTRUE(expr_fun(x))
}

download_pass <- read_bool(file.path(table_dir, "gse73680_supplementary_download_log.tsv"),
                           function(x) any(x$gzip_valid == TRUE & x$readable_header == TRUE)) ||
  read_bool(file.path(table_dir, "gse73680_extracted_txtgz_validation.tsv"),
            function(x) any(x$gzip_valid == TRUE))
txt_pass <- read_bool(file.path(table_dir, "gse73680_txt_structure_audit.tsv"),
                      function(x) any(x$usable_for_matrix == TRUE))
matrix_status <- read_status(file.path(table_dir, "gse73680_matrix_build_log.tsv"))
mapping_status <- read_status(file.path(table_dir, "gse73680_gene_mapping_qc.tsv"))
metadata_status <- read_status(file.path(table_dir, "gse73680_metadata_curated_audit.tsv"))
qc_status <- read_status(file.path(table_dir, "gse73680_expression_qc_v2.tsv"))
gene_avail_ok <- read_bool(file.path(table_dir, "gse73680_p1_gene_availability.tsv"),
                           function(x) x$n_detected[1] >= 3)
module_ok <- read_bool(file.path(table_dir, "gse73680_gene_availability.tsv"),
                       function(x) any(x$usable_for_module_score %in% c("TRUE", "exploratory_module_only")))

summary <- data.table(
  checkpoint = c("supplementary_download", "txt_structure", "matrix_build", "gene_mapping", "metadata_curation", "expression_qc", "gene_availability", "overall_analysis_ready"),
  status = c(ifelse(download_pass, "pass_or_partial", "fail"),
             ifelse(txt_pass, "pass_or_partial", "fail"),
             matrix_status, mapping_status, metadata_status, qc_status,
             ifelse(gene_avail_ok, "pass", "fail"), NA_character_),
  key_metric = c(
    "at least one valid local supplementary or RAW-tar-extracted TXT.gz",
    "at least one usable matrix TXT structure",
    "matrix build status",
    "gene mapping status",
    "metadata curation status",
    "expression QC status",
    "P1 detected genes >= 3",
    "combined readiness"
  ),
  threshold = c(">=1", ">=1", "pass", "pass or warning", "pass", "pass or warning", ">=3", "see rules"),
  pass = c(download_pass, txt_pass, matrix_status == "pass", mapping_status %in% c("pass", "warning"),
           metadata_status == "pass", qc_status %in% c("pass", "warning"), gene_avail_ok, NA),
  blocker = NA_character_,
  action = NA_character_
)
overall <- if (matrix_status == "pass" && mapping_status %in% c("pass", "warning") &&
               metadata_status == "pass" && qc_status %in% c("pass", "warning") &&
               gene_avail_ok && module_ok) {
  "TRUE"
} else if (matrix_status == "pass" && metadata_status == "warning" && qc_status == "warning" && gene_avail_ok) {
  "CONDITIONAL"
} else {
  "FALSE"
}
summary[checkpoint == "overall_analysis_ready", `:=`(status = overall, pass = overall %in% c("TRUE", "CONDITIONAL"))]
summary[pass == FALSE, blocker := paste0(checkpoint, "_not_ready")]
summary[pass == FALSE, action := fifelse(checkpoint == "supplementary_download", "download_valid_supplementary_txt_gz",
                                  fifelse(checkpoint == "matrix_build", "build_expression_matrix_after_download",
                                          "resolve_checkpoint_before_disease_context_analysis"))]
summary[is.na(action), action := "no_action_or_continue_next_checkpoint"]
fwrite(summary, file.path(table_dir, "gse73680_analysis_readiness_summary.tsv"), sep = "\t")

writeLines(c(
  "# GSE73680 Analysis Readiness Memo",
  "",
  paste0("Overall analysis ready: **", overall, "**"),
  "",
  "## Current interpretation boundary",
  "",
  "GSE73680 resources were curated and evaluated for analysis readiness. No disease-context expression conclusion is generated unless overall_analysis_ready is TRUE or CONDITIONAL.",
  "",
  "## Main blockers",
  "",
  if (overall == "FALSE") paste0("- ", paste(summary[pass == FALSE & checkpoint != "overall_analysis_ready", paste(checkpoint, blocker, sep = ": ")], collapse = "\n- ")) else "- No blocking checkpoint for readiness status.",
  "",
  "## Allowed wording",
  "",
  "GSE73680 was curated for plaque/disease-context evaluation, but downstream disease-context analysis depends on passing expression matrix, metadata, gene mapping and QC readiness criteria.",
  "",
  "## Not allowed wording",
  "",
  "- Do not claim P1 genes are differentially expressed in plaque.",
  "- Do not claim GSE73680 validates TAL localization.",
  "- Do not claim disease-context support until analysis-ready criteria are met and analysis is run."
), "docs/gse73680_analysis_readiness_memo.md")
message("wrote GSE73680 analysis readiness memo")
