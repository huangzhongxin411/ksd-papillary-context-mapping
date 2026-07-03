#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
phase_dir <- file.path(project_root, "results", "smr_coloc", "phase27b")
dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "docs"), recursive = TRUE, showWarnings = FALSE)

paths <- list(
  smr_exe = file.path(project_root, "external", "smr", "bin", "smr"),
  smr_zip = file.path(project_root, "external", "smr", "bin", "smr-1.4.1-linux-x86_64.zip"),
  smr_macos_zip = file.path(project_root, "external", "smr", "bin", "smr-1.4.1-macOS-arm64.zip"),
  besd_base = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex"),
  besd = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex.besd"),
  epi = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex.epi"),
  esi = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex.esi"),
  zip = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex.zip"),
  md5 = file.path(project_root, "external", "smr", "gtex_v8_cis_eqtl", "Kidney_Cortex", "Kidney_Cortex.zip.md5sum"),
  md5_check = file.path(phase_dir, "kidney_cortex_zip_md5_check.txt"),
  gwas = file.path(project_root, "data", "processed", "twas_input", "ksd_2025_for_twas.tsv.gz"),
  priority_manifest = file.path(project_root, "results", "smr_coloc", "priority_locus_manifest_v0.1.tsv"),
  twas_magma_overlap = file.path(project_root, "results", "twas", "twas_magma_overlap_real.tsv")
)

out <- list(
  resource_qc = file.path(phase_dir, "kidney_cortex_smr_resource_qc.tsv"),
  gwas_ma = file.path(phase_dir, "ksd_gwas_for_smr.ma"),
  priority_genes = file.path(phase_dir, "priority_genes_for_smr.txt"),
  epi_match = file.path(phase_dir, "priority_gene_epi_match.tsv"),
  all_results = file.path(phase_dir, "smr_heidi_kidney_cortex_all_results.tsv"),
  priority_results = file.path(phase_dir, "smr_heidi_priority_results_v0.2.tsv"),
  priority_summary = file.path(phase_dir, "smr_heidi_priority_summary_v0.2.tsv"),
  run_qc = file.path(phase_dir, "smr_heidi_resource_and_run_qc_v0.2.tsv"),
  memo = file.path(project_root, "docs", "phase27b_smr_heidi_kidney_cortex_memo_v0.1.md")
)

status_row <- function(item, status, detail, path = "") {
  data.frame(item = item, status = status, detail = detail, path = path, stringsAsFactors = FALSE)
}

file_info_row <- function(item, path, min_bytes = 1) {
  if (!file.exists(path)) {
    return(status_row(item, "FAIL", "file_missing", path))
  }
  size <- file.info(path)$size
  status <- if (is.na(size) || size < min_bytes) "FAIL" else "PASS"
  detail <- paste0("size_bytes=", size)
  status_row(item, status, detail, path)
}

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, ...)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

priority_order <- c(
  "UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2",
  "UGT8", "FAM13A", "RGS14", "PTGER1", "GNAZ", "ZCRB1",
  "FBXL20", "CRHR1", "MED1", "ITIH3", "TBX2", "GNL3",
  "TRIOBP", "MAPT"
)
writeLines(priority_order, out$priority_genes)

## Resource QC ---------------------------------------------------------------
qc <- rbind(
  file_info_row("smr_linux_zip_present", paths$smr_zip, 1000000),
  file_info_row("smr_macos_arm64_zip_present", paths$smr_macos_zip, 1000000),
  file_info_row("smr_executable_present", paths$smr_exe, 1000000),
  file_info_row("kidney_cortex_zip_present", paths$zip, 100000000),
  file_info_row("kidney_cortex_md5_present", paths$md5, 20),
  file_info_row("kidney_cortex_besd_present", paths$besd, 100000000),
  file_info_row("kidney_cortex_epi_present", paths$epi, 1000),
  file_info_row("kidney_cortex_esi_present", paths$esi, 1000000),
  file_info_row("gwas_twas_input_present", paths$gwas, 1000000),
  file_info_row("priority_manifest_present", paths$priority_manifest, 1000)
)

md5_status <- "NOT_CHECKED"
md5_detail <- "md5_check_file_missing"
if (file.exists(paths$md5) && file.exists(paths$md5_check)) {
  expected <- strsplit(readLines(paths$md5, warn = FALSE)[1], "[[:space:]]+")[[1]][1]
  observed <- strsplit(readLines(paths$md5_check, warn = FALSE)[1], "[[:space:]]+")[[1]][1]
  md5_status <- if (identical(tolower(expected), tolower(observed))) "PASS" else "FAIL"
  md5_detail <- paste0("expected=", expected, ";observed=", observed)
}
qc <- rbind(qc, status_row("kidney_cortex_zip_md5_match", md5_status, md5_detail, paths$md5_check))

smr_runnable <- FALSE
smr_run_detail <- "not_tested"
if (file.exists(paths$smr_exe)) {
  help_file <- file.path(phase_dir, "smr_runnable_check.txt")
  code <- suppressWarnings(system2(paths$smr_exe, "--help", stdout = help_file, stderr = help_file))
  if (file.exists(help_file)) {
    help_lines <- readLines(help_file, warn = FALSE)
    txt <- paste(head(help_lines, 5), collapse = " | ")
    smr_runnable <- identical(code, 0L) || any(grepl("SMR \\(Summary-data-based Mendelian Randomization\\)|Version 1\\.4\\.1", help_lines))
    smr_run_detail <- paste0("exit_code=", code, ";", txt)
  } else {
    smr_run_detail <- paste0("exit_code=", code, ";no_output")
  }
}
qc <- rbind(qc, status_row(
  "smr_executable_runnable_on_current_host",
  if (smr_runnable) "PASS" else "FAIL",
  smr_run_detail,
  paths$smr_exe
))

write_tsv(qc, out$resource_qc)

## GWAS .ma ------------------------------------------------------------------
ma_summary <- data.frame(metric = character(), value = character())
if (file.exists(paths$gwas)) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required to stream the GWAS TWAS input.")
  }
  gwas <- data.table::fread(paths$gwas, sep = "\t", data.table = FALSE, showProgress = FALSE)
  required <- c("SNP", "A1", "A2", "EAF", "BETA", "SE", "P", "N")
  missing_cols <- setdiff(required, names(gwas))
  if (length(missing_cols) > 0) {
    stop("Missing GWAS columns for SMR .ma: ", paste(missing_cols, collapse = ", "))
  }
  before_n <- nrow(gwas)
  ma <- data.frame(
    SNP = gwas$SNP,
    A1 = toupper(gwas$A1),
    A2 = toupper(gwas$A2),
    freq = suppressWarnings(as.numeric(gwas$EAF)),
    b = suppressWarnings(as.numeric(gwas$BETA)),
    se = suppressWarnings(as.numeric(gwas$SE)),
    p = suppressWarnings(as.numeric(gwas$P)),
    N = suppressWarnings(as.numeric(gwas$N))
  )
  ok <- !is.na(ma$SNP) & ma$SNP != "" &
    ma$A1 %in% c("A", "C", "G", "T") &
    ma$A2 %in% c("A", "C", "G", "T") &
    !is.na(ma$freq) & ma$freq >= 0 & ma$freq <= 1 &
    !is.na(ma$b) &
    !is.na(ma$se) & ma$se > 0 &
    !is.na(ma$p) & ma$p > 0 & ma$p <= 1 &
    !is.na(ma$N) & ma$N > 0
  ma <- ma[ok, ]
  write_tsv(ma, out$gwas_ma)
  ma_summary <- rbind(
    ma_summary,
    data.frame(metric = "gwas_rows_before_filter", value = as.character(before_n)),
    data.frame(metric = "gwas_rows_after_filter", value = as.character(nrow(ma))),
    data.frame(metric = "gwas_rows_removed", value = as.character(before_n - nrow(ma)))
  )
}

## EPI priority matching ------------------------------------------------------
epi_match <- data.frame(
  gene = priority_order,
  probe_id = NA_character_,
  probe_chr = NA_character_,
  probe_bp = NA_character_,
  strand = NA_character_,
  epi_match_status = "resource_missing",
  stringsAsFactors = FALSE
)

if (file.exists(paths$epi)) {
  epi <- read.delim(paths$epi, sep = "\t", header = FALSE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE)
  if (ncol(epi) >= 6) {
    names(epi)[1:6] <- c("probe_chr", "probe_id", "unused", "probe_bp", "gene", "strand")
    for (g in priority_order) {
      hits <- epi[epi$gene == g, ]
      if (nrow(hits) > 0) {
        idx <- match(g, epi_match$gene)
        epi_match$probe_id[idx] <- paste(unique(hits$probe_id), collapse = ";")
        epi_match$probe_chr[idx] <- paste(unique(hits$probe_chr), collapse = ";")
        epi_match$probe_bp[idx] <- paste(unique(hits$probe_bp), collapse = ";")
        epi_match$strand[idx] <- paste(unique(hits$strand), collapse = ";")
        epi_match$epi_match_status[idx] <- "matched"
      } else {
        epi_match$epi_match_status[epi_match$gene == g] <- "not_in_kidney_cortex_epi"
      }
    }
  }
}
write_tsv(epi_match, out$epi_match)

## Placeholder results unless a real SMR output exists ------------------------
priority_manifest <- if (file.exists(paths$priority_manifest)) read_tsv(paths$priority_manifest) else NULL
twas_overlap_genes <- character()
if (file.exists(paths$twas_magma_overlap)) {
  twas_overlap <- read_tsv(paths$twas_magma_overlap)
  if ("gene_name" %in% names(twas_overlap)) twas_overlap_genes <- unique(twas_overlap$gene_name)
  if ("gene" %in% names(twas_overlap)) twas_overlap_genes <- unique(c(twas_overlap_genes, twas_overlap$gene))
}
p1_genes <- character()
if (!is.null(priority_manifest) && all(c("gene", "is_P1") %in% names(priority_manifest))) {
  p1_genes <- priority_manifest$gene[as.logical(priority_manifest$is_P1)]
}

result_cols <- c(
  "ProbeID", "Probe_Chr", "Gene", "Probe_bp", "SNP", "SNP_Chr", "SNP_bp",
  "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL",
  "se_eQTL", "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI",
  "nsnp_HEIDI", "is_priority_gene", "is_P1", "is_TWAS_MAGMA_overlap",
  "smr_fdr", "smr_bonferroni_pass", "heidi_pass", "claim_level"
)

make_resource_limited <- function() {
  data.frame(
    ProbeID = epi_match$probe_id,
    Probe_Chr = epi_match$probe_chr,
    Gene = epi_match$gene,
    Probe_bp = epi_match$probe_bp,
    SNP = NA_character_,
    SNP_Chr = NA_character_,
    SNP_bp = NA_character_,
    A1 = NA_character_,
    A2 = NA_character_,
    Freq = NA_real_,
    b_GWAS = NA_real_,
    se_GWAS = NA_real_,
    p_GWAS = NA_real_,
    b_eQTL = NA_real_,
    se_eQTL = NA_real_,
    p_eQTL = NA_real_,
    b_SMR = NA_real_,
    se_SMR = NA_real_,
    p_SMR = NA_real_,
    p_HEIDI = NA_real_,
    nsnp_HEIDI = NA_real_,
    is_priority_gene = TRUE,
    is_P1 = epi_match$gene %in% p1_genes,
    is_TWAS_MAGMA_overlap = epi_match$gene %in% twas_overlap_genes,
    smr_fdr = NA_real_,
    smr_bonferroni_pass = FALSE,
    heidi_pass = FALSE,
    claim_level = "resource_limited",
    stringsAsFactors = FALSE
  )
}

all_results <- make_resource_limited()
all_results <- all_results[, result_cols]
write_tsv(all_results, out$all_results)
write_tsv(all_results, out$priority_results)

eligible <- all(qc$status[match(c(
  "kidney_cortex_besd_present", "kidney_cortex_epi_present",
  "kidney_cortex_esi_present", "gwas_twas_input_present"
), qc$item)] == "PASS") && smr_runnable && any(epi_match$epi_match_status == "matched")

smr_output <- file.path(phase_dir, "ksd_gtexv8_kidney_cortex_smr_all.smr")
direct_log <- file.path(phase_dir, "ksd_gtexv8_kidney_cortex_smr_all_direct.run.log")
thread1_log <- file.path(phase_dir, "ksd_gtexv8_kidney_cortex_smr_all_thread1.run.log")
query_log <- file.path(phase_dir, "query_umod_kidney_cortex.run.log")
query_nested_log <- file.path(phase_dir, "query_umod_kidney_cortex_nested.run.log")
cmake_log <- file.path(phase_dir, "smr_source_cmake_configure.log")

run_attempted <- file.exists(direct_log) || file.exists(thread1_log) ||
  file.exists(file.path(phase_dir, "ksd_gtexv8_kidney_cortex_smr_all.run.log"))
query_attempted <- file.exists(query_log) || file.exists(query_nested_log)
real_run_completed <- file.exists(smr_output) && file.info(smr_output)$size > 0

run_status <- if (real_run_completed) {
  "REAL_RUN_COMPLETED"
} else if (run_attempted || query_attempted) {
  "NOT_RUN_COMPLETED_SMR_ABORT_OR_BESD_QUERY_SEGFAULT"
} else if (eligible) {
  "READY_TO_RUN_EXTERNAL_SMR"
} else if (!smr_runnable) {
  "NOT_RUN_EXECUTABLE_INCOMPATIBLE_WITH_CURRENT_HOST"
} else {
  "NOT_RUN_RESOURCE_QC_FAILED"
}

summary <- data.frame(
  metric = c(
    "run_completed",
    "run_status",
    "resource_bundle",
    "smr_executable_runnable",
    "priority_genes_requested",
    "priority_genes_matched_epi",
    "strict_supportive_priority_genes",
    "p1_strict_supportive_genes",
    "promotion_rule"
  ),
  value = c(
    as.character(real_run_completed),
    run_status,
    "GTEx_v8_Kidney_Cortex_cis_eQTL_SMR_BESD",
    as.character(smr_runnable),
    as.character(length(priority_order)),
    as.character(sum(epi_match$epi_match_status == "matched")),
    "0",
    "0",
    "Do not promote SMR/HEIDI until real run completes with Bonferroni-significant SMR and HEIDI pass"
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, out$priority_summary)

run_qc <- rbind(
  qc,
  status_row("gwas_ma_created", if (file.exists(out$gwas_ma)) "PASS" else "FAIL",
             if (file.exists(out$gwas_ma)) paste0("size_bytes=", file.info(out$gwas_ma)$size) else "missing", out$gwas_ma),
  status_row("priority_gene_epi_match_created", if (file.exists(out$epi_match)) "PASS" else "FAIL",
             paste0("matched=", sum(epi_match$epi_match_status == "matched"), "/", length(priority_order)), out$epi_match),
  status_row("smr_real_run_completed", if (real_run_completed) "PASS" else "FAIL", run_status, smr_output),
  status_row("smr_all_run_attempt_log", if (file.exists(direct_log) || file.exists(thread1_log)) "INFO" else "WARNING",
             "Direct and/or thread1 SMR run logs record Abort trap 6 on this macOS host.", paste(c(direct_log, thread1_log), collapse = ";")),
  status_row("smr_besd_query_smoke_test", if (query_attempted) "FAIL" else "NOT_RUN",
             "UMOD BESD query produced segmentation fault 11 with the macOS arm64 executable.", paste(c(query_log, query_nested_log), collapse = ";")),
  status_row("smr_source_compile_attempt", if (file.exists(cmake_log)) "FAIL" else "NOT_RUN",
             if (file.exists(cmake_log)) paste(readLines(cmake_log, warn = FALSE), collapse = " ") else "not_attempted", cmake_log),
  status_row("claim_promotion_allowed", "FAIL",
             "No SMR/HEIDI support can be promoted because real SMR did not complete on this host.", "")
)
if (nrow(ma_summary) > 0) {
  for (i in seq_len(nrow(ma_summary))) {
    run_qc <- rbind(run_qc, status_row(ma_summary$metric[i], "INFO", ma_summary$value[i], out$gwas_ma))
  }
}
write_tsv(run_qc, out$run_qc)

memo <- c(
  "# Phase 27B SMR/HEIDI Kidney_Cortex memo",
  "",
  "## Status",
  "",
  paste0("- Resource landing: ", if (all(qc$status[qc$item %in% c("kidney_cortex_besd_present", "kidney_cortex_epi_present", "kidney_cortex_esi_present")] == "PASS")) "PASS" else "INCOMPLETE"),
  paste0("- MD5 check: ", md5_status),
  paste0("- SMR executable runnable on current host: ", smr_runnable),
  paste0("- Run status: ", run_status),
  "",
  "## Interpretation",
  "",
  "GTEx v8 Kidney_Cortex BESD resources were landed and priority genes were matched against the `.epi` file. The official macOS arm64 SMR v1.4.1 executable can start and print the help/version banner, but the real Kidney_Cortex SMR run aborted on this host and a UMOD BESD query smoke test produced a segmentation fault. A source-build attempt was not possible because `cmake` is not installed in the current environment. Therefore, no real SMR/HEIDI result has been generated in this phase.",
  "",
  "## Claim boundary",
  "",
  "SMR/HEIDI remains a resource/execution-limited extension in the current workspace. It must not be promoted to the manuscript, candidate gene tiers, or figure claims until a real run completes and passes strict thresholds: eQTL p < 5e-8, Bonferroni-significant SMR, HEIDI p > 0.01 and nsnp_HEIDI >= 10.",
  "",
  "## Generated files",
  "",
  paste0("- `", sub(project_root, ".", out$resource_qc, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$gwas_ma, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$priority_genes, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$epi_match, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$all_results, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$priority_results, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$priority_summary, fixed = TRUE), "`"),
  paste0("- `", sub(project_root, ".", out$run_qc, fixed = TRUE), "`")
)
writeLines(memo, out$memo)

message("Phase 27B resource audit and boundary outputs complete: ", phase_dir)
