suppressPackageStartupMessages(library(data.table))

empty_table <- function(path, cols) {
  fwrite(as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols)), path, sep = "\t")
}

files <- list.files(
  "results/twas/predixcan/per_tissue",
  pattern = "spredixcan\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

out_all <- "results/tables/spredixcan_all_tissues.tsv"
out_sig <- "results/tables/spredixcan_significant.tsv"
out_sug <- "results/tables/spredixcan_suggestive.tsv"
out_best <- "results/tables/spredixcan_best_per_gene.tsv"
out_status <- "results/tables/spredixcan_run_status.tsv"

if (!length(files)) {
  cols <- c("gene", "tissue", "predixcan_p", "predixcan_fdr", "predixcan_z", "source_file")
  empty_table(out_all, cols)
  empty_table(out_sig, cols)
  empty_table(out_sug, cols)
  empty_table(out_best, cols)
  fwrite(data.table(status = "not_run", notes = "No S-PrediXcan per-tissue result files found."), out_status, sep = "\t")
  message("No S-PrediXcan result files found; wrote empty summary tables.")
  quit(save = "no", status = 0)
}

res_list <- lapply(files, function(f) {
  x <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x[, tissue := sub("\\.spredixcan\\.tsv$", "", basename(f))]
  x[, source_file := f]
  x
})
res <- rbindlist(res_list, fill = TRUE)
if (!nrow(res)) stop("S-PrediXcan result files exist but no rows could be read.")

gcol <- intersect(c("gene", "gene_name", "GENE", "gene_id"), names(res))[1]
pcol <- intersect(c("pvalue", "p_value", "P", "predixcan_p"), names(res))[1]
zcol <- intersect(c("zscore", "z_score", "Z", "predixcan_z"), names(res))[1]
if (is.na(gcol) || is.na(pcol)) {
  stop("Cannot identify gene or p-value columns in S-PrediXcan output.")
}

setnames(res, gcol, "gene")
setnames(res, pcol, "predixcan_p")
if (!is.na(zcol)) setnames(res, zcol, "predixcan_z")
if (!"predixcan_z" %in% names(res)) res[, predixcan_z := NA_real_]

res[, predixcan_p := as.numeric(predixcan_p)]
res[, predixcan_fdr := p.adjust(predixcan_p, method = "BH")]

fwrite(res, out_all, sep = "\t")
fwrite(res[predixcan_fdr < 0.05], out_sig, sep = "\t")
fwrite(res[predixcan_p < 1e-4], out_sug, sep = "\t")
best <- res[order(predixcan_p), .SD[1], by = gene]
fwrite(best, out_best, sep = "\t")

fwrite(data.table(
  status = "completed",
  n_files = length(files),
  n_rows = nrow(res),
  n_genes = uniqueN(res$gene),
  n_fdr_significant = nrow(res[predixcan_fdr < 0.05]),
  n_suggestive = nrow(res[predixcan_p < 1e-4])
), out_status, sep = "\t")

message("wrote\t", out_all)
