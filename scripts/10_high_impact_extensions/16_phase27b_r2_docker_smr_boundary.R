#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), mustWork = TRUE)
phase_dir <- file.path(root, "results", "smr_coloc", "phase27b_docker")
docs_dir <- file.path(root, "docs")
home_workspace <- path.expand("~/ksd_smr27b_rescue")
dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(phase_dir, "out"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(phase_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

exists_size <- function(path) {
  exists <- file.exists(path)
  size <- if (exists) as.numeric(file.info(path)$size) else NA_real_
  list(exists = exists, size = size)
}

docker_path <- Sys.which("docker")
docker_available <- nzchar(docker_path)
host <- paste(system2("uname", "-s", stdout = TRUE), system2("uname", "-m", stdout = TRUE))

input_paths <- list(
  besd = file.path(home_workspace, "input", "Kidney_Cortex.besd"),
  epi = file.path(home_workspace, "input", "Kidney_Cortex.epi"),
  esi = file.path(home_workspace, "input", "Kidney_Cortex.esi"),
  gwas = file.path(home_workspace, "input", "ksd_gwas_for_smr.ma"),
  p1_probes = file.path(home_workspace, "input", "p1_probes.list"),
  priority_probes = file.path(home_workspace, "input", "priority_probes.list"),
  bed = file.path(home_workspace, "ref", "g1000_eur.bed"),
  bim = file.path(home_workspace, "ref", "g1000_eur.bim"),
  fam = file.path(home_workspace, "ref", "g1000_eur.fam")
)

input_qc <- do.call(rbind, lapply(names(input_paths), function(nm) {
  info <- exists_size(input_paths[[nm]])
  data.frame(
    item = paste0("home_workspace_", nm),
    status = if (isTRUE(info$exists) && !is.na(info$size) && info$size > 0) "PASS" else "FAIL",
    notes = paste0(input_paths[[nm]], "; size_bytes=", ifelse(is.na(info$size), "", as.character(info$size))),
    stringsAsFactors = FALSE
  )
}))

docker_qc <- data.frame(
  item = c(
    "host",
    "home_ascii_workspace_created",
    "docker_cli_available",
    "docker_amd64_uname_test",
    "container_setup",
    "linux_smr_download_inside_container",
    "p1_smr_heidi_off_smoke_test",
    "priority_smr_heidi_real_run",
    "claim_boundary_preserved"
  ),
  status = c(
    "INFO",
    ifelse(dir.exists(home_workspace), "PASS", "FAIL"),
    ifelse(docker_available, "PASS", "FAIL"),
    ifelse(docker_available, "NOT_RUN", "FAIL"),
    "NOT_RUN",
    "NOT_RUN",
    "NOT_RUN",
    "NOT_RUN",
    "PASS"
  ),
  notes = c(
    host,
    home_workspace,
    ifelse(docker_available, docker_path, "docker command not found; Docker Desktop/CLI unavailable from terminal"),
    ifelse(docker_available, "Not executed by this boundary script", "Docker linux/amd64 execution unavailable on current Mac"),
    "Stopped before container setup because docker CLI is unavailable",
    "Stopped before container setup because docker CLI is unavailable",
    "Not run; per instruction, stop if Docker linux/amd64 check fails",
    "Not run; p1 smoke test unavailable",
    "No manuscript, Figure 1-5, candidate_gene_tiers, or integrated evidence files modified"
  ),
  stringsAsFactors = FALSE
)

qc <- rbind(input_qc, docker_qc)
write_tsv(qc, file.path(phase_dir, "docker_smr_execution_qc.tsv"))

rescue_result <- file.path(root, "results", "smr_coloc", "phase27b_rescue", "smr_heidi_priority_results_rescue_v0.1.tsv")
if (file.exists(rescue_result)) {
  res <- read_tsv(rescue_result)
} else {
  cols <- c(
    "ProbeID", "Probe_Chr", "Gene", "Probe_bp", "SNP", "SNP_Chr", "SNP_bp",
    "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL",
    "se_eQTL", "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI",
    "nsnp_HEIDI", "is_priority_gene", "is_P1", "is_TWAS_MAGMA_overlap",
    "smr_fdr", "smr_bonferroni_pass", "heidi_pass", "claim_level"
  )
  res <- data.frame(matrix(nrow = 0, ncol = length(cols)))
  names(res) <- cols
}
res$claim_level <- "execution_limited_after_docker"
res$smr_bonferroni_pass <- FALSE
res$heidi_pass <- FALSE
write_tsv(res, file.path(phase_dir, "smr_heidi_priority_results_docker.tsv"))

summary <- data.frame(
  metric = c(
    "host",
    "home_ascii_workspace",
    "docker_cli_available",
    "docker_linux_amd64_available",
    "p1_smr_heidi_off_completed",
    "priority_smr_heidi_completed",
    "strict_supportive_priority_genes",
    "claim_promotion_allowed",
    "final_status"
  ),
  value = c(
    host,
    home_workspace,
    as.character(docker_available),
    "FALSE",
    "FALSE",
    "FALSE",
    "0",
    "FALSE",
    "execution_limited_after_docker"
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(phase_dir, "docker_smr_execution_summary.tsv"))

log_text <- c(
  "Phase 27B-R2 Docker linux/amd64 check",
  paste0("Host: ", host),
  paste0("Home workspace: ", home_workspace),
  paste0("Docker CLI available: ", docker_available),
  ifelse(docker_available, paste0("Docker path: ", docker_path), "docker command not found"),
  "Decision: Docker linux/amd64 execution unavailable on current Mac; SMR/HEIDI frozen as execution-limited."
)
writeLines(log_text, file.path(phase_dir, "logs", "docker_precheck.log"))

memo <- c(
  "# Phase 27B-R2 Docker SMR completion memo",
  "",
  "## Status",
  "",
  paste0("- Host: ", host),
  paste0("- Home ASCII workspace: `", home_workspace, "`"),
  paste0("- Docker CLI available: ", docker_available),
  "- Docker linux/amd64 container test: NOT RUN / unavailable",
  "- P1 HEIDI-off smoke test: NOT RUN",
  "- Priority-probe SMR/HEIDI: NOT RUN",
  "",
  "## Decision",
  "",
  "Docker linux/amd64 execution is unavailable on the current Mac because the `docker` command is not found. Following the Phase 27B-R2 instruction, no further time was spent on SMR rescue. SMR/HEIDI is frozen as `execution_limited_after_docker`.",
  "",
  "## Claim boundary",
  "",
  "No SMR/HEIDI support is available for manuscript claims. Do not modify the manuscript, main Figures 1-5, candidate gene tiers, or integrated evidence summary based on this phase. The only acceptable statement is that GTEx v8 Kidney_Cortex SMR resources were landed, but execution remained limited after macOS and Docker rescue attempts.",
  "",
  "## Outputs",
  "",
  "- `results/smr_coloc/phase27b_docker/docker_smr_execution_qc.tsv`",
  "- `results/smr_coloc/phase27b_docker/docker_smr_execution_summary.tsv`",
  "- `results/smr_coloc/phase27b_docker/smr_heidi_priority_results_docker.tsv`",
  "- `results/smr_coloc/phase27b_docker/logs/docker_precheck.log`"
)
writeLines(memo, file.path(docs_dir, "phase27b_docker_smr_completion_memo_v0.1.md"))

message("Phase 27B-R2 Docker boundary outputs complete: ", phase_dir)
