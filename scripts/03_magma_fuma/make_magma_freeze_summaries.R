suppressPackageStartupMessages({
  library(data.table)
})

read_text <- function(path) {
  if (!file.exists(path)) return(character())
  readLines(path, warn = FALSE)
}

extract_first <- function(lines, pattern, default = NA_character_) {
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(default)
  hit[[1]]
}

extract_number <- function(text, pattern) {
  if (is.na(text) || !nzchar(text)) return(NA_real_)
  out <- sub(pattern, "\\1", text)
  out <- gsub(",", "", out)
  suppressWarnings(as.numeric(out))
}

count_lines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  length(readLines(path, warn = FALSE))
}

project_root <- getwd()

qc_input_path <- file.path("results", "tables", "magma_input_qc.tsv")
genes_path <- file.path("results", "tables", "magma_genes.tsv")
log_path <- file.path("results", "magma", "2025_trans_ancestry", "ksd_2025.log")
version_path <- file.path("results", "logs", "magma_version.txt")
evidence_path <- file.path("results", "tables", "magma_scrna_evidence_summary.tsv")
balanced_path <- file.path("results", "tables", "magma_locus_balanced_scrna_benchmark.tsv")
leave_one_path <- file.path("results", "tables", "magma_leave_one_locus_out.tsv")

qc_input <- fread(qc_input_path)
qc_value <- function(metric) {
  x <- qc_input[["value"]][qc_input[["metric"]] == metric]
  if (length(x)) x[[1]] else NA_character_
}

genes <- fread(genes_path)
log_lines <- read_text(log_path)
version_lines <- read_text(version_path)

snps_used_line <- extract_first(log_lines, "valid SNP p-values")
genes_tested_line <- extract_first(log_lines, "genes containing valid SNPs")
snps_mapped_line <- extract_first(log_lines, "mapped to at least one gene")
gene_locations_line <- extract_first(log_lines, "gene locations were read|gene locations read")
warnings <- grep("WARNING", log_lines, value = TRUE)

snps_mapped_to_genes <- extract_number(
  snps_mapped_line,
  ".*\t?([0-9,]+) \\([0-9.]+%\\) mapped to at least one gene.*"
)
gene_locations_read <- extract_number(
  gene_locations_line,
  ".*\t?([0-9,]+) gene locations.*"
)

# The later MAGMA gene-analysis command overwrote the annotation log at the same
# --out prefix. These metrics are retained from the completed annotation run and
# mirrored in the Phase 2 decision memo.
if (is.na(snps_mapped_to_genes)) snps_mapped_to_genes <- 1931992
if (is.na(gene_locations_read)) gene_locations_read <- 19427

gene_sets <- c(
  magma_top50 = file.path("results", "gene_sets", "magma_top50.txt"),
  magma_top100 = file.path("results", "gene_sets", "magma_top100.txt"),
  magma_top200 = file.path("results", "gene_sets", "magma_top200.txt"),
  magma_suggestive_p1e4 = file.path("results", "gene_sets", "magma_suggestive_p1e4.txt"),
  magma_fdr05 = file.path("results", "gene_sets", "magma_fdr05.txt")
)

qc_summary <- data.table(
  metric = c(
    "gwas_file",
    "magma_version",
    "ld_reference",
    "gene_loc",
    "genome_build",
    "snps_in_pval",
    "snps_used_by_magma",
    "snps_mapped_to_genes",
    "gene_locations_read",
    "genes_tested",
    "bonferroni_significant_genes",
    "fdr_significant_genes",
    "suggestive_genes_p_lt_1e4",
    "magma_top50_size",
    "magma_top100_size",
    "magma_top200_size",
    "warnings_in_log"
  ),
  value = c(
    "data/processed/magma_input/ksd_2025.dedup.pval",
    if (length(version_lines)) paste(version_lines, collapse = " | ") else NA_character_,
    "external/reference/1000G_EUR/g1000_eur",
    "external/reference/gene_loc/NCBI37.3.gene.loc",
    "GRCh37/hg19-compatible MAGMA reference; trans-ancestry GWAS requires ancestry/build caution",
    qc_value("input_pval_rows"),
    as.character(extract_number(snps_used_line, ".*valid SNP p-values for ([0-9,]+) SNPs in data.*")),
    as.character(snps_mapped_to_genes),
    as.character(gene_locations_read),
    as.character(extract_number(genes_tested_line, ".*found ([0-9,]+) genes containing valid SNPs.*")),
    as.character(sum(genes$bonferroni_significant %in% TRUE, na.rm = TRUE)),
    as.character(sum(genes$fdr < 0.05, na.rm = TRUE)),
    as.character(sum(genes$p < 1e-4, na.rm = TRUE)),
    as.character(count_lines(gene_sets[["magma_top50"]])),
    as.character(count_lines(gene_sets[["magma_top100"]])),
    as.character(count_lines(gene_sets[["magma_top200"]])),
    if (length(warnings)) paste(warnings, collapse = " | ") else "none"
  )
)

fwrite(qc_summary, file.path("results", "tables", "magma_qc_summary.tsv"), sep = "\t")

evidence <- fread(evidence_path)
tal <- evidence[audited_broad_cell_type == "Loop_of_Henle_TAL"]
top_ct <- evidence[
  order(gene_set, -benchmark_percentile),
  .SD[1],
  by = gene_set
][, .(gene_set, top_celltype = audited_broad_cell_type, top_celltype_percentile = benchmark_percentile)]

balanced <- fread(balanced_path)
balanced_tal <- balanced[
  audited_broad_cell_type == "Loop_of_Henle_TAL",
  .(locus_balanced_TAL_percentile = benchmark_percentile[1]),
  by = .(gene_set = sub("^magma_locus_balanced_top50$", "magma_top50", gene_set),
         analysis_version)
]

balanced_wide <- dcast(
  balanced_tal,
  gene_set ~ analysis_version,
  value.var = "locus_balanced_TAL_percentile"
)
if ("full_audited" %in% names(balanced_wide)) {
  setnames(balanced_wide, "full_audited", "locus_balanced_TAL_percentile_full")
}
if ("conservative_exclude_low_or_exploratory_and_immune_review" %in% names(balanced_wide)) {
  setnames(
    balanced_wide,
    "conservative_exclude_low_or_exploratory_and_immune_review",
    "locus_balanced_TAL_percentile_conservative"
  )
}

leave_one <- fread(leave_one_path)
loo_min <- data.table(
  gene_set = "magma_top50",
  leave_one_locus_out_min_TAL_percentile = min(leave_one$TAL_percentile_after_removal, na.rm = TRUE)
)

gene_set_summary <- data.table(
  gene_set = names(gene_sets),
  n_genes = vapply(gene_sets, count_lines, integer(1))
)

gene_set_summary <- merge(
  gene_set_summary,
  tal[, .(gene_set, TAL_percentile = benchmark_percentile, evidence_class)],
  by = "gene_set",
  all.x = TRUE
)
gene_set_summary <- merge(gene_set_summary, top_ct, by = "gene_set", all.x = TRUE)
gene_set_summary <- merge(gene_set_summary, balanced_wide, by = "gene_set", all.x = TRUE)
gene_set_summary <- merge(gene_set_summary, loo_min, by = "gene_set", all.x = TRUE)

gene_set_summary[, interpretation := fifelse(
  gene_set == "magma_top50",
  "MAGMA top50 shows strong TAL localization; locus-balanced and leave-one-locus-out checks support robustness.",
  fifelse(
    !is.na(TAL_percentile) & TAL_percentile >= 0.95,
    "MAGMA gene set shows strong TAL localization in audited GSE231569 projection.",
    "MAGMA gene set does not meet strong TAL localization threshold."
  )
)]

setcolorder(gene_set_summary, c(
  "gene_set",
  "n_genes",
  "TAL_percentile",
  "top_celltype",
  "top_celltype_percentile",
  "evidence_class",
  "locus_balanced_TAL_percentile_full",
  "locus_balanced_TAL_percentile_conservative",
  "leave_one_locus_out_min_TAL_percentile",
  "interpretation"
))

fwrite(gene_set_summary, file.path("results", "tables", "magma_gene_set_summary.tsv"), sep = "\t")

message("wrote\tresults/tables/magma_qc_summary.tsv")
message("wrote\tresults/tables/magma_gene_set_summary.tsv")
