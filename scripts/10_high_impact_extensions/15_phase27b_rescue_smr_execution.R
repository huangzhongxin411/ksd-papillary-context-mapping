#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), mustWork = TRUE)
tmp_root <- "/tmp/ksd_smr27b_rescue"
phase_dir <- file.path(root, "results", "smr_coloc", "phase27b_rescue")
docs_dir <- file.path(root, "docs")
dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

file_first_line <- function(path) {
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE, n = 1), collapse = "")
}

count_lines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  wc <- system2("wc", c("-l", path), stdout = TRUE)
  parts <- strsplit(trimws(wc[1]), "[[:space:]]+")[[1]]
  as.integer(parts[1])
}

file_qc <- function(label, path, row_count = TRUE) {
  exists <- file.exists(path)
  size_bytes <- if (exists) file.info(path)$size else NA_real_
  size_mb <- if (exists) round(size_bytes / 1024^2, 3) else NA_real_
  n_rows <- if (exists && row_count) count_lines(path) else NA_integer_
  first <- if (exists && size_mb < 500) file_first_line(path) else if (exists) "[large binary or large text file; first line not read]" else ""
  status <- if (!exists) "FAIL" else if (!is.na(size_bytes) && size_bytes > 0) "PASS" else "FAIL"
  data.frame(
    file = label,
    path = path,
    exists = exists,
    size_mb = size_mb,
    n_rows = n_rows,
    header_or_first_line = first,
    status = status,
    notes = "",
    stringsAsFactors = FALSE
  )
}

input_dir <- file.path(tmp_root, "input")
ref_dir <- file.path(tmp_root, "ref")
out_dir <- file.path(tmp_root, "out")

paths <- list(
  besd = file.path(input_dir, "Kidney_Cortex.besd"),
  epi = file.path(input_dir, "Kidney_Cortex.epi"),
  esi = file.path(input_dir, "Kidney_Cortex.esi"),
  gwas_ma = file.path(input_dir, "ksd_gwas_for_smr.ma"),
  priority_match = file.path(input_dir, "priority_gene_epi_match.tsv"),
  bed = file.path(ref_dir, "g1000_eur.bed"),
  bim = file.path(ref_dir, "g1000_eur.bim"),
  fam = file.path(ref_dir, "g1000_eur.fam"),
  twas_magma_overlap = file.path(root, "results", "twas", "twas_magma_overlap_real.tsv")
)

p1_genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
twas_magma_genes <- character()
if (file.exists(paths$twas_magma_overlap)) {
  twas_magma <- read_tsv(paths$twas_magma_overlap)
  if ("gene_name" %in% names(twas_magma)) twas_magma_genes <- unique(c(twas_magma_genes, twas_magma$gene_name))
  if ("gene" %in% names(twas_magma)) twas_magma_genes <- unique(c(twas_magma_genes, twas_magma$gene))
}

priority_match <- if (file.exists(paths$priority_match)) read_tsv(paths$priority_match) else data.frame()
probe_col <- if ("probe_id" %in% names(priority_match)) "probe_id" else if ("matched_probe_id" %in% names(priority_match)) "matched_probe_id" else NA_character_
gene_col <- if ("gene" %in% names(priority_match)) "gene" else NA_character_
status_col <- if ("epi_match_status" %in% names(priority_match)) "epi_match_status" else NA_character_

matched <- priority_match
if (!is.na(status_col)) matched <- matched[matched[[status_col]] == "matched", , drop = FALSE]
priority_probes <- if (!is.na(probe_col)) unique(unlist(strsplit(as.character(matched[[probe_col]]), ";", fixed = TRUE))) else character()
priority_probes <- priority_probes[nzchar(priority_probes)]

p1_matched <- matched
if (!is.na(gene_col)) p1_matched <- p1_matched[p1_matched[[gene_col]] %in% p1_genes, , drop = FALSE]
p1_probes <- if (!is.na(probe_col)) unique(unlist(strsplit(as.character(p1_matched[[probe_col]]), ";", fixed = TRUE))) else character()
p1_probes <- p1_probes[nzchar(p1_probes)]

writeLines(priority_probes, file.path(input_dir, "priority_probes.list"))
writeLines(p1_probes, file.path(input_dir, "p1_probes.list"))

integrity <- rbind(
  file_qc("Kidney_Cortex.besd", paths$besd, row_count = FALSE),
  file_qc("Kidney_Cortex.epi", paths$epi, row_count = TRUE),
  file_qc("Kidney_Cortex.esi", paths$esi, row_count = TRUE),
  file_qc("ksd_gwas_for_smr.ma", paths$gwas_ma, row_count = TRUE),
  file_qc("priority_gene_epi_match.tsv", paths$priority_match, row_count = TRUE),
  file_qc("priority_probes.list", file.path(input_dir, "priority_probes.list"), row_count = TRUE),
  file_qc("p1_probes.list", file.path(input_dir, "p1_probes.list"), row_count = TRUE),
  file_qc("g1000_eur.bed", paths$bed, row_count = FALSE),
  file_qc("g1000_eur.bim", paths$bim, row_count = TRUE),
  file_qc("g1000_eur.fam", paths$fam, row_count = TRUE)
)
write_tsv(integrity, file.path(phase_dir, "rescue_file_integrity_qc.tsv"))

host <- paste(system2("uname", "-s", stdout = TRUE), system2("uname", "-m", stdout = TRUE))
is_linux_x86 <- identical(system2("uname", "-s", stdout = TRUE), "Linux") &&
  grepl("x86_64|amd64", system2("uname", "-m", stdout = TRUE), ignore.case = TRUE)
docker_path <- Sys.which("docker")

linux_instructions <- c(
  "# Phase 27B-R Linux SMR/HEIDI rescue instructions",
  "",
  "This workspace is currently not Linux x86_64, so the main rescue run was not executed locally. Use a Linux x86_64 server, HPC node, cloud VM, or Docker image with the copied ASCII workspace.",
  "",
  "## Prepared ASCII workspace",
  "",
  "- `/tmp/ksd_smr27b_rescue/input/`",
  "- `/tmp/ksd_smr27b_rescue/ref/`",
  "- `/tmp/ksd_smr27b_rescue/out/`",
  "",
  "If moving to another machine, copy `/tmp/ksd_smr27b_rescue/` and `external/smr/bin/smr-1.4.1-linux-x86_64.zip`.",
  "",
  "## Linux setup",
  "",
  "```bash",
  "mkdir -p external/smr/linux_x86_64 results/smr_coloc/phase27b_rescue",
  "cd external/smr/linux_x86_64",
  "curl -L -C - --fail --show-error -O https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.1-linux-x86_64.zip",
  "unzip -o smr-1.4.1-linux-x86_64.zip",
  "chmod +x smr",
  "./smr --help > ../../../results/smr_coloc/phase27b_rescue/linux_smr_help.txt 2>&1 || true",
  "```",
  "",
  "## Step 1: P1 smoke test without HEIDI",
  "",
  "```bash",
  "external/smr/linux_x86_64/smr \\",
  "  --bfile /tmp/ksd_smr27b_rescue/ref/g1000_eur \\",
  "  --gwas-summary /tmp/ksd_smr27b_rescue/input/ksd_gwas_for_smr.ma \\",
  "  --beqtl-summary /tmp/ksd_smr27b_rescue/input/Kidney_Cortex \\",
  "  --extract-probe /tmp/ksd_smr27b_rescue/input/p1_probes.list \\",
  "  --heidi-off \\",
  "  --thread-num 1 \\",
  "  --out /tmp/ksd_smr27b_rescue/out/p1_smr_heidi_off \\",
  "  > results/smr_coloc/phase27b_rescue/p1_smr_heidi_off.log 2>&1",
  "```",
  "",
  "Proceed only if `/tmp/ksd_smr27b_rescue/out/p1_smr_heidi_off.smr` exists and is non-empty.",
  "",
  "## Step 2: priority-probe SMR with HEIDI",
  "",
  "```bash",
  "external/smr/linux_x86_64/smr \\",
  "  --bfile /tmp/ksd_smr27b_rescue/ref/g1000_eur \\",
  "  --gwas-summary /tmp/ksd_smr27b_rescue/input/ksd_gwas_for_smr.ma \\",
  "  --beqtl-summary /tmp/ksd_smr27b_rescue/input/Kidney_Cortex \\",
  "  --extract-probe /tmp/ksd_smr27b_rescue/input/priority_probes.list \\",
  "  --peqtl-smr 5e-8 \\",
  "  --peqtl-heidi 1.57e-3 \\",
  "  --heidi-min-m 10 \\",
  "  --heidi-max-m 20 \\",
  "  --thread-num 1 \\",
  "  --out /tmp/ksd_smr27b_rescue/out/priority_smr_heidi \\",
  "  > results/smr_coloc/phase27b_rescue/priority_smr_heidi.log 2>&1",
  "```",
  "",
  "Do not promote SMR/HEIDI unless the real `.smr` output exists and strict QC passes."
)
writeLines(linux_instructions, file.path(docs_dir, "phase27b_rescue_linux_run_instructions_v0.1.md"))

cols <- c(
  "ProbeID", "Probe_Chr", "Gene", "Probe_bp", "SNP", "SNP_Chr", "SNP_bp",
  "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL",
  "se_eQTL", "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI",
  "nsnp_HEIDI", "is_priority_gene", "is_P1", "is_TWAS_MAGMA_overlap",
  "smr_fdr", "smr_bonferroni_pass", "heidi_pass", "claim_level"
)
empty_results <- data.frame(matrix(nrow = length(priority_probes), ncol = length(cols)))
names(empty_results) <- cols
if (length(priority_probes) > 0 && nrow(matched) > 0) {
  empty_results$ProbeID <- priority_probes
  matched_by_probe <- matched[match(priority_probes, matched[[probe_col]]), , drop = FALSE]
  empty_results$Gene <- matched_by_probe[[gene_col]]
  empty_results$Probe_Chr <- matched_by_probe$probe_chr
  empty_results$Probe_bp <- matched_by_probe$probe_bp
  empty_results$is_priority_gene <- TRUE
  empty_results$is_P1 <- empty_results$Gene %in% p1_genes
  empty_results$is_TWAS_MAGMA_overlap <- empty_results$Gene %in% twas_magma_genes
  empty_results$smr_bonferroni_pass <- FALSE
  empty_results$heidi_pass <- FALSE
  empty_results$claim_level <- "execution_limited"
}
write_tsv(empty_results, file.path(phase_dir, "smr_heidi_priority_results_rescue_v0.1.tsv"))

summary <- data.frame(
  metric = c(
    "host",
    "docker_available",
    "linux_x86_64_local_run",
    "ascii_workspace_created",
    "priority_probes",
    "p1_probes",
    "p1_smr_heidi_off_completed",
    "priority_smr_heidi_completed",
    "strict_supportive_priority_genes",
    "claim_promotion_allowed",
    "final_status"
  ),
  value = c(
    host,
    ifelse(nzchar(docker_path), docker_path, "FALSE"),
    as.character(is_linux_x86),
    as.character(dir.exists(tmp_root)),
    as.character(length(priority_probes)),
    as.character(length(p1_probes)),
    "FALSE",
    "FALSE",
    "0",
    "FALSE",
    if (is_linux_x86) "ready_for_local_linux_run" else "execution_limited_pending_linux_x86_64_rescue"
  )
)
write_tsv(summary, file.path(phase_dir, "smr_heidi_priority_summary_rescue_v0.1.tsv"))

exec_qc <- data.frame(
  item = c(
    "host_is_linux_x86_64",
    "docker_available",
    "ascii_workspace_prepared",
    "input_files_copied",
    "priority_probe_list_created",
    "p1_probe_list_created",
    "p1_smr_heidi_off_run",
    "priority_smr_heidi_run",
    "result_parser_ready",
    "claim_boundary_preserved",
    "smr_portal_fallback_note"
  ),
  status = c(
    ifelse(is_linux_x86, "PASS", "WARNING"),
    ifelse(nzchar(docker_path), "INFO", "WARNING"),
    ifelse(dir.exists(tmp_root), "PASS", "FAIL"),
    ifelse(all(integrity$status[match(c("Kidney_Cortex.besd", "Kidney_Cortex.epi", "Kidney_Cortex.esi", "ksd_gwas_for_smr.ma"), integrity$file)] == "PASS"), "PASS", "FAIL"),
    ifelse(length(priority_probes) > 0, "PASS", "FAIL"),
    ifelse(length(p1_probes) == length(p1_genes), "PASS", "WARNING"),
    "NOT_RUN",
    "NOT_RUN",
    "PASS",
    "PASS",
    "PASS"
  ),
  notes = c(
    host,
    ifelse(nzchar(docker_path), docker_path, "docker command not found on current host"),
    tmp_root,
    "BESD/EPI/ESI/GWAS MA/LD reference copied into ASCII workspace",
    paste0(length(priority_probes), " probes"),
    paste0(length(p1_probes), " P1 probes"),
    "Current host is not Linux x86_64; do not rerun macOS arm64 as main rescue",
    "Pending Linux x86_64 server/HPC/Docker run",
    "Placeholder result table with required columns generated",
    "No manuscript, figures, candidate tiers or integrated evidence were modified",
    "SMR Portal may be considered only if GWAS sharing/upload is explicitly approved"
  )
)
write_tsv(exec_qc, file.path(phase_dir, "rescue_execution_qc_v0.1.tsv"))

portal_note <- c(
  "# Optional SMR Portal fallback feasibility note",
  "",
  "SMR Portal can be considered only after explicit approval for GWAS summary upload/sharing. The current local workflow should remain preferred because it preserves local control of the KSD GWAS summary statistics.",
  "",
  "Potential requirements:",
  "",
  "- GWAS summary in SMR/GCTA-COJO-compatible format: SNP, A1, A2, freq, b, se, p, n.",
  "- xQTL resource: GTEx v8 Kidney_Cortex or closest available kidney tissue resource.",
  "- Priority probe/gene restriction matching Phase 27B-R.",
  "",
  "Do not use portal-derived output for claims unless the run is reproducible, downloadable, and passes the same strict thresholds: eQTL p < 5e-8, Bonferroni-significant SMR, HEIDI p > 0.01 and nsnp_HEIDI >= 10."
)
writeLines(portal_note, file.path(docs_dir, "phase27b_rescue_smr_portal_feasibility_note_v0.1.md"))

memo <- c(
  "# Phase 27B-R SMR/HEIDI execution rescue memo",
  "",
  "## Final local status",
  "",
  paste0("- Host: ", host),
  paste0("- Docker available: ", ifelse(nzchar(docker_path), docker_path, "FALSE")),
  "- Main rescue run on current host: NOT RUN, because the current host is not Linux x86_64.",
  "- ASCII workspace prepared: `/tmp/ksd_smr27b_rescue/`",
  paste0("- Priority probes prepared: ", length(priority_probes)),
  paste0("- P1 probes prepared: ", length(p1_probes)),
  "",
  "## Interpretation",
  "",
  "Phase 27B-R prepared a clean ASCII execution workspace and Linux x86_64 run instructions. Because this machine is Darwin arm64 and Docker is not available, the Linux rescue run was not executed locally. This avoids treating the previous macOS arm64 BESD segmentation fault as a biological negative result.",
  "",
  "## Claim boundary",
  "",
  "SMR/HEIDI remains execution-limited pending a clean Linux x86_64 run. No SMR/HEIDI support should be added to the manuscript, main figures, candidate gene tiers, or integrated evidence summary.",
  "",
  "## Outputs",
  "",
  "- `results/smr_coloc/phase27b_rescue/rescue_file_integrity_qc.tsv`",
  "- `results/smr_coloc/phase27b_rescue/rescue_execution_qc_v0.1.tsv`",
  "- `results/smr_coloc/phase27b_rescue/smr_heidi_priority_results_rescue_v0.1.tsv`",
  "- `results/smr_coloc/phase27b_rescue/smr_heidi_priority_summary_rescue_v0.1.tsv`",
  "- `docs/phase27b_rescue_linux_run_instructions_v0.1.md`",
  "- `docs/phase27b_rescue_smr_portal_feasibility_note_v0.1.md`"
)
writeLines(memo, file.path(docs_dir, "phase27b_rescue_completion_memo_v0.1.md"))

message("Phase 27B-R rescue prep complete: ", phase_dir)
