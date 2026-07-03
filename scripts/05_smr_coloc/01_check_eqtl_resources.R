suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/smr_coloc", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
magma <- fread("results/tables/magma_genes.tsv")
top <- magma[order(p)][1:50]
priority_genes <- unique(c("UMOD", "CLDN14", "CLDN10", "CASR", "PKD2", "HIBADH", "FAM13A", top$gene_symbol[1:20]))
priority <- merge(data.table(gene = priority_genes), magma[, .(gene = gene_symbol, MAGMA_P = p, MAGMA_FDR = fdr, MAGMA_rank = rank)],
                  by = "gene", all.x = TRUE)
priority <- merge(priority, p1[, .(gene, interpretation_role = manuscript_role)], by = "gene", all.x = TRUE)
priority[, window_kb := "500_and_1000"]
priority[, status := "priority_for_resource_landing"]
fwrite(priority, "results/smr_coloc/priority_loci_for_smr_coloc.tsv", sep = "\t")

paths <- list.files(c("external/eqtl", "external/twas", "results/phase2b_twas_ready_v0.1"), recursive = TRUE, full.names = TRUE)
paths <- paths[!grepl("/tests?/", paths)]
eqtl_hits <- paths[grepl("eqtl|GTEx|kidney|Kidney|whole.*blood|Whole_Blood|\\.besd$|\\.esi$|\\.epi$", paths, ignore.case = TRUE)]
coloc_eqtl_hits <- paths[grepl("eqtl|GTEx|eQTL", paths, ignore.case = TRUE) & grepl("\\.txt\\.gz$|\\.tsv\\.gz$", paths, ignore.case = TRUE)]
status <- data.table(
  resource = c("GTEx kidney cortex eQTL", "GTEx whole blood eQTL", "SMR BESD/ESI/EPI", "coloc-ready eQTL summary"),
  status = c(
    ifelse(any(grepl("kidney|Kidney", eqtl_hits, ignore.case = TRUE)), "candidate_found", "missing"),
    ifelse(any(grepl("whole.*blood|Whole_Blood", eqtl_hits, ignore.case = TRUE)), "candidate_found", "missing"),
    ifelse(any(grepl("\\.besd$|\\.esi$|\\.epi$", eqtl_hits, ignore.case = TRUE)), "candidate_found", "missing"),
    ifelse(length(coloc_eqtl_hits) > 0, "candidate_found", "missing")
  ),
  path_or_note = c(
    paste(head(eqtl_hits[grepl("kidney|Kidney", eqtl_hits, ignore.case = TRUE)], 5), collapse = "; "),
    paste(head(eqtl_hits[grepl("whole.*blood|Whole_Blood", eqtl_hits, ignore.case = TRUE)], 5), collapse = "; "),
    paste(head(eqtl_hits[grepl("\\.besd$|\\.esi$|\\.epi$", eqtl_hits, ignore.case = TRUE)], 5), collapse = "; "),
    paste(head(coloc_eqtl_hits, 5), collapse = "; ")
  )
)
status[path_or_note == "", path_or_note := "not found in current workspace"]
fwrite(status, "results/smr_coloc/eqtl_resource_status.tsv", sep = "\t")
writeLines(c(
  "# SMR/coloc Resource Blocker Memo v0.1",
  "",
  "Targeted SMR/coloc is not run unless locus-matched eQTL summary statistics, allele harmonization fields and LD/build compatibility are confirmed.",
  "Current outputs define priority loci and resource status only."
), "docs/smr_coloc_resource_blocker_memo.md", useBytes = TRUE)
message("wrote SMR/coloc resource status")
