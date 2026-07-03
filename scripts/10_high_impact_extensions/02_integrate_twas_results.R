suppressPackageStartupMessages(library(data.table))

dir.create("results/twas", recursive = TRUE, showWarnings = FALSE)

read_predixcan <- function(file) {
  x <- fread(file, fill = TRUE)
  if (!nrow(x)) return(NULL)
  x[, method := "S-PrediXcan"]
  x[, tissue := sub("\\.spredixcan\\.tsv$", "", basename(file))]
  x[, source_file := file]
  gene_col <- intersect(c("gene", "gene_name", "GENE", "gene_id", "gene_name"), names(x))[1]
  p_col <- intersect(c("pvalue", "p_value", "P", "predixcan_p"), names(x))[1]
  z_col <- intersect(c("zscore", "z_score", "Z", "predixcan_z"), names(x))[1]
  id_col <- intersect(c("gene_id", "gene"), names(x))[1]
  if (is.na(gene_col) || is.na(p_col)) return(NULL)
  data.table(
    method = x$method,
    tissue = x$tissue,
    gene = as.character(x[[gene_col]]),
    gene_id = if (!is.na(id_col)) as.character(x[[id_col]]) else NA_character_,
    z = if (!is.na(z_col)) as.numeric(x[[z_col]]) else NA_real_,
    p = as.numeric(x[[p_col]]),
    source_file = x$source_file
  )
}

files <- list.files("results/twas/predixcan/per_tissue", pattern = "spredixcan\\.tsv$", full.names = TRUE)
res <- rbindlist(lapply(files, read_predixcan), fill = TRUE)

cols <- c("method", "tissue", "gene", "gene_id", "z", "p", "fdr", "source_file", "support_status", "claim_use")
if (!nrow(res)) {
  fwrite(as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols)),
         "results/twas/twas_results.tsv", sep = "\t")
  fwrite(data.table(status = "not_run_or_no_readable_results", n_results = 0L),
         "results/twas/twas_run_status.tsv", sep = "\t")
  quit(save = "no", status = 0)
}

res <- res[!is.na(p)]
res[, fdr := p.adjust(p, method = "BH")]
res[, support_status := fifelse(fdr < 0.05, "FDR_supported", fifelse(p < 1e-4, "suggestive_only", "not_supported"))]
res[, claim_use := fifelse(fdr < 0.05, "eligible_as_TWAS_enhancement", "not_used_for_main_claim")]
setcolorder(res, cols)
fwrite(res[order(fdr, p)], "results/twas/twas_results.tsv", sep = "\t")

magma_path <- "results/phase2_magma_v0.1/tables/magma_genes.tsv"
if (!file.exists(magma_path)) magma_path <- "submission_package_v1.9_BMC_working/supplementary/tables/S3_magma_gene_results_and_gene_sets.tsv"
magma <- fread(magma_path, fill = TRUE)
if ("gene_symbol" %in% names(magma)) setnames(magma, "gene_symbol", "gene")
magma_keep <- intersect(c("gene", "rank", "p", "fdr"), names(magma))
magma <- unique(magma[, ..magma_keep])
if ("p" %in% names(magma)) setnames(magma, "p", "magma_p")
if ("fdr" %in% names(magma)) setnames(magma, "fdr", "magma_fdr")
if ("rank" %in% names(magma)) setnames(magma, "rank", "magma_rank")

best <- res[order(p), .SD[1], by = gene]
ov <- merge(best, magma, by = "gene", all.x = TRUE)
ov[, overlap_class := fifelse(!is.na(magma_rank) & fdr < 0.05, "MAGMA_overlap_FDR_TWAS", fifelse(!is.na(magma_rank), "MAGMA_overlap_not_FDR_TWAS", "TWAS_only"))]
ov[, claim_use := fifelse(fdr < 0.05 & !is.na(magma_rank), "eligible_as_TWAS_enhancement", "not_used_for_main_claim")]
fwrite(ov[order(fdr, p)], "results/twas/twas_magma_overlap.tsv", sep = "\t")

fwrite(data.table(
  status = "completed",
  n_results = nrow(res),
  n_tissues = uniqueN(res$tissue),
  n_fdr_supported = nrow(res[fdr < 0.05]),
  n_magma_overlap_fdr_supported = nrow(ov[fdr < 0.05 & !is.na(magma_rank)]),
  claim_boundary = "Only FDR-supported TWAS/MAGMA overlaps may be described as enhanced genetic evidence."
), "results/twas/twas_run_status.tsv", sep = "\t")
