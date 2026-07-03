suppressPackageStartupMessages(library(data.table))

dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)
dir.create("results/smr", recursive = TRUE, showWarnings = FALSE)
dir.create("results/coloc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/spatial", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

exists_nonempty <- function(x) file.exists(x) & file.info(x)$size > 0
find_files <- function(dir, pattern) {
  if (!dir.exists(dir)) return(character())
  list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}
first_or_na <- function(x) if (length(x)) x[1] else NA_character_

venv_py <- "external/twas/predixcan/metaxcan_venv/bin/python3"
spredixcan <- "external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"
gwas <- "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
mashr_dir <- "external/twas/predixcan/predictdb_gtex_v8_mashr"
mashr_tar <- file.path(mashr_dir, "mashr_eqtl.tar")
kidney_db <- first_or_na(find_files(mashr_dir, "Kidney.*Cortex.*\\.db$|Kidney_Cortex.*\\.db$|mashr_Kidney_Cortex\\.db$"))
kidney_cov <- first_or_na(find_files(mashr_dir, "Kidney.*Cortex.*\\.txt\\.gz$|Kidney_Cortex.*\\.txt\\.gz$|mashr_Kidney_Cortex\\.txt\\.gz$"))

py_import_ok <- FALSE
py_versions <- "python env not found"
if (exists_nonempty(venv_py)) {
  cmd <- paste(shQuote(venv_py), "- <<'PY'\nimport numpy, scipy, pandas, h5py\nprint('numpy=' + numpy.__version__)\nprint('scipy=' + scipy.__version__)\nprint('pandas=' + pandas.__version__)\nprint('h5py=' + h5py.__version__)\nPY")
  res <- tryCatch(system(cmd, intern = TRUE), error = function(e) conditionMessage(e))
  py_import_ok <- all(grepl("=", res)) && length(res) >= 4
  py_versions <- paste(res, collapse = "; ")
}

twas_resource <- data.table(
  step = c(
    "Step1_directory",
    "Step1_mashr_tar",
    "Step1_kidney_model_db",
    "Step1_kidney_covariance",
    "Step2_python_venv",
    "Step2_python_imports",
    "Step3_gwas_input",
    "Step3_spredixcan_script",
    "Step3_kidney_run_output",
    "Step4_real_summary",
    "Step4_magma_overlap"
  ),
  resource = c(
    mashr_dir,
    mashr_tar,
    ifelse(is.na(kidney_db), "mashr_Kidney_Cortex.db", kidney_db),
    ifelse(is.na(kidney_cov), "mashr_Kidney_Cortex.txt.gz", kidney_cov),
    venv_py,
    "numpy/scipy/pandas/h5py",
    gwas,
    spredixcan,
    "results/twas/spredixcan_Kidney_Cortex_KSD.csv",
    "results/twas/twas_results_real_summary.tsv",
    "results/twas/twas_magma_overlap_real.tsv"
  ),
  status = c(
    if (dir.exists(mashr_dir)) "available" else "missing",
    if (exists_nonempty(mashr_tar)) "available" else "missing_download_failed_or_not_placed",
    if (!is.na(kidney_db) && exists_nonempty(kidney_db)) "available" else "missing",
    if (!is.na(kidney_cov) && exists_nonempty(kidney_cov)) "available" else "missing",
    if (exists_nonempty(venv_py)) "available" else "missing",
    if (py_import_ok) "available" else "missing",
    if (exists_nonempty(gwas)) "available" else "missing",
    if (exists_nonempty(spredixcan)) "available" else "missing",
    if (exists_nonempty("results/twas/spredixcan_Kidney_Cortex_KSD.csv")) "available" else "not_run",
    if (exists_nonempty("results/twas/twas_results_real_summary.tsv")) "available" else "pending_real_run",
    if (exists_nonempty("results/twas/twas_magma_overlap_real.tsv")) "available" else "pending_real_run"
  ),
  note = c(
    "Directory created for PredictDB GTEx v8 MASHR resource.",
    "Download from Zenodo failed in this run; manually place mashr_eqtl.tar here if needed.",
    "Detected actual filename if present; expected Kidney_Cortex model db.",
    "Detected actual filename if present; expected Kidney_Cortex covariance.",
    "Project-local venv; system Python was not modified.",
    py_versions,
    "Header inspected: SNP CHR BP A1 A2 BETA SE Z P N EAF.",
    "MetaXcan SPrediXcan help command succeeded after dependency install.",
    "Run only after model db and covariance are available.",
    "Header-only/blocked summary is allowed until real run completes.",
    "Only FDR-supported TWAS/MAGMA overlap may be promoted."
  )
)
fwrite(twas_resource, "results/twas/phase25b_spredixcan_resource_landing.tsv", sep = "\t")

twas_resource_status <- data.table(
  resource = c(
    "GWAS TWAS input",
    "S-PrediXcan software",
    "S-PrediXcan Python environment",
    "PredictDB MASHR tar",
    "PredictDB Kidney_Cortex model db",
    "PredictDB Kidney_Cortex covariance",
    "Kidney_Cortex S-PrediXcan output",
    "TWAS real summary",
    "TWAS MAGMA overlap"
  ),
  status = c(
    if (exists_nonempty(gwas)) "available" else "missing",
    if (exists_nonempty(spredixcan)) "available" else "missing",
    if (py_import_ok) "available" else "missing",
    if (exists_nonempty(mashr_tar)) "available" else "missing_download_failed_or_not_placed",
    if (!is.na(kidney_db) && exists_nonempty(kidney_db)) "available" else "missing",
    if (!is.na(kidney_cov) && exists_nonempty(kidney_cov)) "available" else "missing",
    if (exists_nonempty("results/twas/spredixcan_Kidney_Cortex_KSD.csv")) "available" else "not_run",
    if (exists_nonempty("results/twas/twas_results_real_summary.tsv")) "available" else "pending_real_run",
    if (exists_nonempty("results/twas/twas_magma_overlap_real.tsv")) "available" else "pending_real_run"
  ),
  path_or_note = c(
    gwas,
    spredixcan,
    paste0(venv_py, "; ", py_versions),
    mashr_tar,
    ifelse(is.na(kidney_db), "not found under external/twas/predixcan/predictdb_gtex_v8_mashr", kidney_db),
    ifelse(is.na(kidney_cov), "not found under external/twas/predixcan/predictdb_gtex_v8_mashr", kidney_cov),
    "results/twas/spredixcan_Kidney_Cortex_KSD.csv",
    "results/twas/twas_results_real_summary.tsv",
    "results/twas/twas_magma_overlap_real.tsv"
  ),
  usable_for_main_analysis = c(
    TRUE,
    TRUE,
    py_import_ok,
    exists_nonempty(mashr_tar),
    !is.na(kidney_db) && exists_nonempty(kidney_db),
    !is.na(kidney_cov) && exists_nonempty(kidney_cov),
    exists_nonempty("results/twas/spredixcan_Kidney_Cortex_KSD.csv"),
    exists_nonempty("results/twas/twas_results_real_summary.tsv") && any(fread("results/twas/twas_results_real_summary.tsv", fill = TRUE)$status == "completed"),
    exists_nonempty("results/twas/twas_magma_overlap_real.tsv") && file.info("results/twas/twas_magma_overlap_real.tsv")$size > 100
  ),
  action_decision = c(
    "ready",
    "ready",
    ifelse(py_import_ok, "ready", "fix_python_dependencies"),
    ifelse(exists_nonempty(mashr_tar), "ready", "manual_place_or_retry_download"),
    ifelse(!is.na(kidney_db) && exists_nonempty(kidney_db), "ready", "extract_mashr_eqtl_tar"),
    ifelse(!is.na(kidney_cov) && exists_nonempty(kidney_cov), "ready", "extract_mashr_eqtl_tar"),
    ifelse(exists_nonempty("results/twas/spredixcan_Kidney_Cortex_KSD.csv"), "completed_check_qc", "run_after_model_and_covariance_available"),
    "do_not_promote_until_completed_and_FDR_checked",
    "only_FDR_supported_overlap_can_enhance_claim"
  )
)
fwrite(twas_resource_status, "results/twas/twas_resource_status.tsv", sep = "\t")

run_exists <- exists_nonempty("results/twas/spredixcan_Kidney_Cortex_KSD.csv")
if (!run_exists) {
  fwrite(data.table(
    status = "not_run_resource_missing",
    tissue = "Kidney_Cortex",
    n_genes_tested = NA_integer_,
    n_genes_with_n_snps_used_gt0 = NA_integer_,
    min_p = NA_real_,
    n_fdr_lt_0_05 = NA_integer_,
    n_magma_overlap_fdr_lt_0_05 = NA_integer_,
    claim_use = "not_used_for_evidence",
    blocking_reason = paste(twas_resource[status %chin% c("missing", "missing_download_failed_or_not_placed"), paste(step, resource, sep = ": ")], collapse = "; ")
  ), "results/twas/twas_results_real_summary.tsv", sep = "\t")
  fwrite(data.table(
    gene = character(),
    twas_p = numeric(),
    twas_fdr = numeric(),
    magma_rank = integer(),
    magma_p = numeric(),
    magma_fdr = numeric(),
    p1_candidate = logical(),
    support_status = character(),
    claim_use = character()
  ), "results/twas/twas_magma_overlap_real.tsv", sep = "\t")
}

smr_files <- find_files("external/smr/gtex_v8", "\\.(besd|esi|epi)$")
smr <- data.table(
  resource = c("GTEx V8 cis-eQTL SMR BESD lite/full", ".besd", ".esi", ".epi"),
  status = c(
    if (length(smr_files)) "candidate_files_found" else "missing",
    if (any(grepl("\\.besd$", smr_files))) "available" else "missing",
    if (any(grepl("\\.esi$", smr_files))) "available" else "missing",
    if (any(grepl("\\.epi$", smr_files))) "available" else "missing"
  ),
  path_or_note = c(
    if (length(smr_files)) paste(head(smr_files, 10), collapse = "; ") else "Place GTEx V8 cis-eQTL SMR BESD lite/full files under external/smr/gtex_v8/",
    first_or_na(smr_files[grepl("\\.besd$", smr_files)]),
    first_or_na(smr_files[grepl("\\.esi$", smr_files)]),
    first_or_na(smr_files[grepl("\\.epi$", smr_files)])
  ),
  decision = c(
    "Do not run SMR until .besd/.esi/.epi are complete.",
    "required",
    "required",
    "required"
  )
)
fwrite(smr, "results/smr/gtex_v8_smr_resource_status.tsv", sep = "\t")

coloc_plan <- "results/coloc/priority_locus_coloc_plan_v0.1.tsv"
if (file.exists(coloc_plan)) {
  coloc <- fread(coloc_plan, fill = TRUE)
} else {
  coloc <- data.table(gene = c("CASR", "CLDN10", "CLDN14", "HIBADH", "PKD2", "UMOD"),
                      chr = c(3, 13, 21, 7, 4, 16),
                      locus_start = NA_integer_,
                      locus_end = NA_integer_)
}
mixqtl_files <- find_files("external/eqtl/mixqtl_kidney_cortex", "\\.(parquet|txt\\.gz|tsv\\.gz)$")
allpairs_files <- find_files("external/eqtl/gtex_kidney_cortex", "allpairs|all_pairs|egenes|signif")
coloc[, required_chromosome := paste0("chr", chr)]
coloc[, kidney_eqtl_resource_status := ifelse(length(mixqtl_files) || length(allpairs_files), "candidate_files_found", "missing")]
coloc[, candidate_eqtl_files := paste(head(c(mixqtl_files, allpairs_files), 20), collapse = "; ")]
coloc[, next_action := ifelse(kidney_eqtl_resource_status == "missing",
                              "Missing variant-level Kidney_Cortex eQTL; do not run coloc.",
                              "Filter eQTL to locus, harmonize genome build/effect allele, then run coloc.")]
fwrite(coloc, "results/coloc/priority_locus_coloc_plan_v0.2.tsv", sep = "\t")

spatial_dir <- "data/raw/spatial/GSE231630"
spatial_files <- find_files(spatial_dir, "filtered_feature_bc_matrix\\.h5|matrix\\.mtx\\.gz|barcodes\\.tsv\\.gz|features\\.tsv\\.gz|tissue_positions|scalefactors_json\\.json|hires|lowres|\\.png$")
spatial_status <- data.table(
  dataset = "GSE231630_or_HuBMAP",
  expression_matrix_or_h5 = any(grepl("filtered_feature_bc_matrix\\.h5|matrix\\.mtx\\.gz", spatial_files)),
  barcodes = any(grepl("barcodes\\.tsv\\.gz", spatial_files)),
  features = any(grepl("features\\.tsv\\.gz", spatial_files)),
  positions = any(grepl("tissue_positions", spatial_files)),
  scalefactors = any(grepl("scalefactors_json\\.json", spatial_files)),
  image = any(grepl("hires|lowres|\\.png$", spatial_files)),
  n_candidate_files = length(spatial_files)
)
spatial_status[, spatial_ready := expression_matrix_or_h5 & positions & scalefactors & image]
spatial_status[, decision := ifelse(spatial_ready, "can_attempt_projection", "resource_missing_do_not_claim_spatial_validation")]
fwrite(spatial_status, "results/spatial/spatial_projection_resource_status_v0.3.tsv", sep = "\t")

memo <- c(
  "# Phase25B resource landing status",
  "",
  "## Step 1. S-PrediXcan resource landing",
  "",
  paste0("- Directory: `", mashr_dir, "`"),
  paste0("- `mashr_eqtl.tar`: ", ifelse(exists_nonempty(mashr_tar), "available", "missing; Zenodo download failed and no local copy was found.")),
  paste0("- Kidney_Cortex model db: ", ifelse(!is.na(kidney_db), kidney_db, "missing")),
  paste0("- Kidney_Cortex covariance: ", ifelse(!is.na(kidney_cov), kidney_cov, "missing")),
  "",
  "## Step 2. Python environment",
  "",
  paste0("- Project-local venv: `", venv_py, "`"),
  paste0("- Import status: ", ifelse(py_import_ok, "PASS", "FAIL")),
  paste0("- Versions/log: ", py_versions),
  "",
  "## Step 3-4. Kidney_Cortex S-PrediXcan and TWAS QC",
  "",
  paste0("- GWAS input header: `SNP CHR BP A1 A2 BETA SE Z P N EAF`."),
  if (run_exists) {
    "- Kidney_Cortex S-PrediXcan was run and produced `results/twas/spredixcan_Kidney_Cortex_KSD.csv`."
  } else {
    "- Kidney_Cortex S-PrediXcan was not run because the PredictDB model db and covariance are missing."
  },
  if (run_exists) {
    "- `twas_results_real_summary.tsv` records the real run QC; use only FDR-supported results according to the claim boundary."
  } else {
    "- `twas_results_real_summary.tsv` records `not_run_resource_missing` and must not be used as evidence."
  },
  "",
  "## Step 5. SMR resource audit",
  "",
  paste0("- SMR resource status written to `results/smr/gtex_v8_smr_resource_status.tsv`. Complete `.besd/.esi/.epi` files are required before running SMR."),
  "",
  "## Step 6. coloc resource audit",
  "",
  "- `results/coloc/priority_locus_coloc_plan_v0.2.tsv` records priority loci and missing Kidney_Cortex variant-level eQTL resources.",
  "",
  "## Step 7. Spatial resource audit",
  "",
  "- `results/spatial/spatial_projection_resource_status_v0.3.tsv` prioritizes GSE231630/HuBMAP Visium resources and remains resource-limited unless complete matrix, positions, scalefactors and images are present.",
  "",
  "## Claim boundary",
  "",
  "No TWAS, SMR, coloc or spatial result was promoted. Only real FDR-supported TWAS results, complete SMR/coloc support, or complete spatial projection outputs may be used in later manuscript revisions."
)
writeLines(memo, "docs/phase25b_resource_landing_status.md")
writeLines(memo, "docs/phase25b_twas_real_run_memo.md")

message("wrote Phase25B status files")
