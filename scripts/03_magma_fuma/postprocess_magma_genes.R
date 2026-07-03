suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

magma_out <- "results/magma/2025_trans_ancestry/ksd_2025.genes.out"
ranked_path <- "results/tables/phase1_candidate_genes_ranked_protein_coding.tsv"
driver_path <- "results/tables/audited_locus_top50_tal_driver_genes.tsv"
out_table <- "results/tables/magma_genes.tsv"
overlap_out <- "results/tables/magma_vs_locus_overlap.tsv"
gene_set_dir <- "results/gene_sets"
dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)

magma <- fread(magma_out)
setnames(magma, old = names(magma), new = tolower(names(magma)))
magma[, gene_id := as.character(gene)]

map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(magma$gene_id),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)
map <- as.data.table(map)
setnames(map, c("ENTREZID", "SYMBOL", "GENENAME"), c("gene_id", "gene_symbol", "gene_name"))
map <- map[!duplicated(gene_id)]

magma <- merge(magma, map, by = "gene_id", all.x = TRUE)
magma[, fdr := p.adjust(p, method = "BH")]
magma[, rank := frank(p, ties.method = "first")]
magma[, bonferroni_significant := p < 0.05 / sum(!is.na(p))]
magma[, suggestive := p < 1e-4]
setorder(magma, p)

out <- magma[, .(
  gene_id,
  gene_symbol,
  gene_name,
  chr,
  start,
  stop,
  n_snps = nsnps,
  nparam,
  n,
  zstat,
  p,
  fdr,
  rank,
  bonferroni_significant,
  suggestive
)]
fwrite(out, out_table, sep = "\t")

write_gene_set <- function(dt, path, n = NULL) {
  x <- dt[!is.na(gene_symbol) & gene_symbol != ""]
  if (!is.null(n)) x <- head(x, n)
  writeLines(unique(x$gene_symbol), path)
}
write_gene_set(out[order(p)], file.path(gene_set_dir, "magma_top50.txt"), 50)
write_gene_set(out[order(p)], file.path(gene_set_dir, "magma_top100.txt"), 100)
write_gene_set(out[order(p)], file.path(gene_set_dir, "magma_top200.txt"), 200)
writeLines(unique(out[p < 1e-4 & !is.na(gene_symbol), gene_symbol]), file.path(gene_set_dir, "magma_suggestive_p1e4.txt"))
writeLines(unique(out[fdr < 0.05 & !is.na(gene_symbol), gene_symbol]), file.path(gene_set_dir, "magma_fdr05.txt"))

top_sets <- rbindlist(list(
  out[rank <= 50, .(gene = gene_symbol, gene_set = "magma_top50", magma_p = p, magma_fdr = fdr, magma_rank = rank)],
  out[rank <= 100, .(gene = gene_symbol, gene_set = "magma_top100", magma_p = p, magma_fdr = fdr, magma_rank = rank)]
))
fwrite(top_sets[!is.na(gene)], "results/tables/magma_top50_top100_gene_sets.tsv", sep = "\t")

ranked <- fread(ranked_path)
drivers <- if (file.exists(driver_path)) fread(driver_path) else data.table(gene = character(), driver_class = character())
locus_top50 <- scan("results/gene_sets/locus_top50.txt", what = character(), quiet = TRUE)
locus_top100 <- scan("results/gene_sets/locus_top100.txt", what = character(), quiet = TRUE)
magma_top50 <- scan(file.path(gene_set_dir, "magma_top50.txt"), what = character(), quiet = TRUE)
magma_top100 <- scan(file.path(gene_set_dir, "magma_top100.txt"), what = character(), quiet = TRUE)

all_genes <- sort(unique(c(locus_top50, locus_top100, magma_top50, magma_top100)))
overlap <- data.table(gene = all_genes)
overlap[, in_locus_top50 := gene %in% locus_top50]
overlap[, in_locus_top100 := gene %in% locus_top100]
overlap[, in_magma_top50 := gene %in% magma_top50]
overlap[, in_magma_top100 := gene %in% magma_top100]
overlap <- merge(
  overlap,
  unique(ranked[, .(gene, locus_id, nearest_lead_snp, lead_snp_p, mapping_method, locus_rank = rank_in_candidate_list)]),
  by = "gene",
  all.x = TRUE
)
overlap <- merge(
  overlap,
  out[, .(gene = gene_symbol, magma_rank = rank, magma_p = p, magma_fdr = fdr)],
  by = "gene",
  all.x = TRUE
)
overlap <- merge(
  overlap,
  unique(drivers[, .(gene, tal_driver_class = driver_class)]),
  by = "gene",
  all.x = TRUE
)
overlap[, is_umod_related_locus := locus_id == "KSD_L018" | gene %in% c("UMOD", "PDILT", "ACSM5")]
overlap[, is_tal_driver_gene := tal_driver_class %in% c("candidate_TAL_driver", "supporting_TAL_expressed_gene")]
setorder(overlap, -in_magma_top50, -in_locus_top50, magma_rank, locus_rank)
fwrite(overlap, overlap_out, sep = "\t")

cat("wrote\t", out_table, "\n", sep = "")
cat("wrote\t", overlap_out, "\n", sep = "")
cat("wrote\t", file.path(gene_set_dir, "magma_top50.txt"), "\n", sep = "")
cat("wrote\t", file.path(gene_set_dir, "magma_top100.txt"), "\n", sep = "")
cat("wrote\t", file.path(gene_set_dir, "magma_top200.txt"), "\n", sep = "")
