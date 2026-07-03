suppressPackageStartupMessages({
  library(data.table)
})

ranked <- fread("results/tables/phase1_candidate_genes_ranked_protein_coding.tsv")
drivers <- fread("results/tables/audited_locus_top50_tal_driver_genes.tsv")
loo <- fread("results/tables/audited_locus_leave_one_locus_out.tsv")
loci <- fread("results/tables/phase1_2025_loci.tsv")

out_path <- "results/tables/umod_locus_context_summary.tsv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

umod_locus_id <- ranked[gene == "UMOD", locus_id][1]
if (is.na(umod_locus_id)) {
  stop("UMOD was not found in ranked candidate genes")
}

locus_genes <- ranked[locus_id == umod_locus_id]
driver_genes <- drivers[locus_id == umod_locus_id]
loo_row <- loo[removed_locus_id == umod_locus_id][1]
locus_row <- loci[locus_id == umod_locus_id][1]

summary <- data.table(
  locus_id = umod_locus_id,
  lead_snp = locus_row$top_lead_snp,
  lead_snp_p = locus_row$min_p,
  locus_chr = locus_row$CHR,
  locus_start = locus_row$start,
  locus_end = locus_row$end,
  genes_in_locus = paste(locus_genes$gene, collapse = ","),
  top50_genes_in_locus = paste(locus_genes[rank_in_candidate_list <= 50, gene], collapse = ","),
  driver_genes = paste(driver_genes[driver_class %in% c("candidate_TAL_driver", "supporting_TAL_expressed_gene"), gene], collapse = ","),
  mapping_methods = paste(unique(locus_genes$mapping_method), collapse = ","),
  distance_to_lead_snp = paste(paste0(locus_genes$gene, ":", locus_genes$distance_to_lead_snp), collapse = ";"),
  TAL_avg_expression = paste(paste0(driver_genes$gene, ":", signif(driver_genes$TAL_avg_expression, 4)), collapse = ";"),
  TAL_pct_expressed = paste(paste0(driver_genes$gene, ":", signif(driver_genes$TAL_pct_expressed, 4)), collapse = ";"),
  TAL_percentile_before_removal = loo_row$TAL_percentile_before_removal,
  TAL_percentile_after_removal = loo_row$TAL_percentile_after_removal,
  leave_one_locus_out_delta = loo_row$delta_percentile,
  interpretation = "UMOD-containing locus is the most influential locus for the conservative TAL top50 signal; retain as exploratory locus-context observation pending MAGMA/TWAS confirmation."
)

fwrite(summary, out_path, sep = "\t")
cat("wrote\t", out_path, "\n", sep = "")
