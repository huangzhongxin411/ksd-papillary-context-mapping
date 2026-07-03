suppressPackageStartupMessages(library(data.table))

empty_table <- function(path, cols) {
  fwrite(as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols)), path, sep = "\t")
}

files <- list.files(
  "results/twas/fusion/per_tissue",
  pattern = "fusion\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

out_all <- "results/tables/twas_fusion_all_tissues.tsv"
out_sig <- "results/tables/twas_fusion_significant.tsv"
out_sug <- "results/tables/twas_fusion_suggestive.tsv"
out_best <- "results/tables/twas_fusion_best_per_gene.tsv"
out_overlap <- "results/tables/twas_fusion_magma_overlap.tsv"
out_status <- "results/tables/twas_fusion_run_status.tsv"

if (!length(files)) {
  cols <- c("gene", "tissue", "fusion_p", "fusion_fdr", "fusion_z", "source_file")
  empty_table(out_all, cols)
  empty_table(out_sig, cols)
  empty_table(out_sug, cols)
  empty_table(out_best, cols)
  empty_table(out_overlap, c(cols, "magma_rank", "magma_p", "magma_fdr", "pre_twas_tier", "scrna_top_celltype", "post_fusion_status"))
  fwrite(data.table(status = "not_run", notes = "No FUSION per-tissue result files found."), out_status, sep = "\t")
  message("No FUSION result files found; wrote empty summary tables.")
  quit(save = "no", status = 0)
}

res_list <- lapply(files, function(f) {
  x <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x[, tissue := basename(dirname(f))]
  x[, source_file := f]
  x
})

res <- rbindlist(res_list, fill = TRUE)
if (!nrow(res)) {
  stop("FUSION result files exist but no rows could be read.")
}

pcol <- intersect(c("TWAS.P", "P", "PANEL.P", "twas_p", "fusion_p"), names(res))[1]
gcol <- intersect(c("FILE", "ID", "GENE", "gene", "fusion_gene_raw"), names(res))[1]
zcol <- intersect(c("TWAS.Z", "Z", "z", "fusion_z"), names(res))[1]

if (is.na(pcol) || is.na(gcol)) {
  stop("Cannot identify gene or p-value columns in FUSION output.")
}

setnames(res, pcol, "fusion_p")
setnames(res, gcol, "fusion_gene_raw")
if (!is.na(zcol)) setnames(res, zcol, "fusion_z")
if (!"fusion_z" %in% names(res)) res[, fusion_z := NA_real_]

res[, fusion_p := as.numeric(fusion_p)]
res[, gene := fusion_gene_raw]
res[, gene := sub("\\.wgt\\.RDat$", "", basename(gene))]
res[, gene := sub("\\.RDat$", "", gene)]
res[, fusion_fdr := p.adjust(fusion_p, method = "BH")]

fwrite(res, out_all, sep = "\t")
fwrite(res[fusion_fdr < 0.05], out_sig, sep = "\t")
fwrite(res[fusion_p < 1e-4], out_sug, sep = "\t")

best <- res[order(fusion_p), .SD[1], by = gene]
fwrite(best, out_best, sep = "\t")

magma <- fread("results/tables/magma_genes.tsv")
setnames(magma, "gene_symbol", "gene", skip_absent = TRUE)
magma_small <- unique(magma[, .(
  gene,
  magma_rank = rank,
  magma_p = p,
  magma_fdr = fdr,
  is_magma_top50 = rank <= 50,
  is_magma_top100 = rank <= 100,
  is_magma_fdr = fdr < 0.05
)], by = "gene")

tiers <- fread("results/tables/candidate_gene_tiers_v0.1.tsv")
tiers_small <- unique(tiers[, .(
  gene,
  pre_twas_tier = tier,
  scrna_top_celltype,
  scrna_benchmark_percentile,
  scrna_evidence_class
)], by = "gene")

overlap <- merge(best, magma_small, by = "gene", all.x = TRUE)
overlap <- merge(overlap, tiers_small, by = "gene", all.x = TRUE)
overlap[, post_fusion_status := fifelse(
  fusion_fdr < 0.05,
  "fusion_strong",
  fifelse(fusion_p < 1e-4, "fusion_suggestive", "fusion_not_significant")
)]

fwrite(overlap, out_overlap, sep = "\t")
fwrite(data.table(
  status = "completed",
  n_files = length(files),
  n_rows = nrow(res),
  n_genes = uniqueN(res$gene),
  n_fdr_significant = nrow(res[fusion_fdr < 0.05]),
  n_suggestive = nrow(res[fusion_p < 1e-4])
), out_status, sep = "\t")

message("wrote\t", out_all)
message("wrote\t", out_overlap)
