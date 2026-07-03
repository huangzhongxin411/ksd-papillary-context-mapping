suppressPackageStartupMessages(library(data.table))

priority <- fread("results/tables/priority_loci_for_smr_coloc.tsv")
tiers <- fread("results/tables/candidate_gene_tiers_v0.2.tsv")
loci <- fread("results/tables/phase1_2025_loci.tsv")
magma <- fread("results/tables/magma_genes.tsv")
phase1_genes <- fread("results/tables/phase1_candidate_genes.tsv")
driver_genes <- fread("results/tables/magma_top50_tal_driver_genes.tsv")
gwas_path <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"

setnames(magma, "gene_symbol", "gene", skip_absent = TRUE)

priority_genes <- unique(c(
  "UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2", "FAM13A",
  magma[rank <= 50, gene],
  tiers[!is.na(magma_fdr) & magma_fdr < 0.05 & scrna_top_celltype == "Loop_of_Henle_TAL", gene],
  tiers[!is.na(magma_p) & magma_p < 1e-4 & scrna_top_celltype == "Loop_of_Henle_TAL", gene]
))

base <- merge(
  tiers[gene %in% priority_genes],
  priority[, .(gene, recommended_for_smr, recommended_for_coloc)],
  by = "gene",
  all.x = TRUE
)
phase1_gene_locus <- phase1_genes[
  gene %in% priority_genes & !is.na(locus_id) & locus_id != "",
  .SD[1],
  by = gene
][, .(
  gene,
  fallback_locus_id = locus_id,
  fallback_lead_snp = nearest_lead_snp,
  fallback_lead_snp_p = lead_snp_p
)]
base <- merge(base, phase1_gene_locus, by = "gene", all.x = TRUE)
driver_locus <- driver_genes[
  gene %in% priority_genes & !is.na(locus_group) & locus_group != "",
  .SD[1],
  by = gene
][, .(
  gene,
  driver_locus_id = locus_group
)]
base <- merge(base, driver_locus, by = "gene", all.x = TRUE)
base[is.na(locus_id) | locus_id == "", locus_id := fallback_locus_id]
base[is.na(locus_id) | locus_id == "", locus_id := driver_locus_id]
base[is.na(lead_snp) | lead_snp == "", lead_snp := fallback_lead_snp]
base[is.na(lead_snp_p), lead_snp_p := fallback_lead_snp_p]
base[is.na(recommended_for_smr), recommended_for_smr := TRUE]
base[is.na(recommended_for_coloc), recommended_for_coloc := TRUE]

base[, locus_id_primary := sub(";.*$", "", locus_id)]
out <- merge(
  base,
  loci[, .(
    locus_id,
    chr = CHR,
    locus_start = start,
    locus_end = end,
    locus_lead_snps = lead_snps,
    locus_min_p = min_p
  )],
  by.x = "locus_id_primary",
  by.y = "locus_id",
  all.x = TRUE
)

out[, required_eqtl_resource := fifelse(
  scrna_top_celltype == "Loop_of_Henle_TAL",
  "kidney_cortex_or_kidney_relevant_eQTL",
  "tissue_matched_or_multi_tissue_eQTL"
)]
out[, status := "resource_preparation_only"]
out[, priority_reason_v02 := fifelse(
  gene %in% c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2", "FAM13A"),
  "predefined_priority_gene_for_SMR_coloc_preparation",
  fifelse(!is.na(magma_rank) & magma_rank <= 50,
    "MAGMA_top50_priority_for_SMR_coloc_preparation",
    "MAGMA_scRNA_priority_for_SMR_coloc_preparation"
  )
)]

priority_v02 <- out[, .(
  gene,
  locus_id,
  lead_snp,
  lead_snp_p,
  chr,
  lead_snp_bp = NA_integer_,
  locus_start,
  locus_end,
  magma_rank,
  magma_p,
  magma_fdr,
  scrna_top_celltype,
  scrna_benchmark_percentile,
  twas_status,
  priority_reason = priority_reason_v02,
  recommended_for_smr,
  recommended_for_coloc,
  required_eqtl_resource,
  resource_status = status
)]
priority_v02 <- unique(priority_v02, by = "gene")
priority_v02[, priority_class := fifelse(
  gene %in% c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2"),
  "P1_core_TAL_candidate",
  fifelse(
    scrna_top_celltype == "Loop_of_Henle_TAL" & !is.na(magma_fdr) & magma_fdr < 0.05 & gene != "FAM13A",
    "P2_Magma_TAL_candidate",
    fifelse(gene == "FAM13A", "P4_exploratory", "P3_context_candidate")
  )
)]
setcolorder(priority_v02, c(
  "gene", "priority_class", "locus_id", "lead_snp", "lead_snp_p", "chr", "lead_snp_bp",
  "locus_start", "locus_end", "magma_rank", "magma_p", "magma_fdr",
  "scrna_top_celltype", "scrna_benchmark_percentile", "twas_status",
  "priority_reason", "recommended_for_smr", "recommended_for_coloc",
  "required_eqtl_resource", "resource_status"
))
setorder(priority_v02, magma_rank, gene)

gwas <- fread(gwas_path)
setnames(gwas, c("CHR", "BP"), c("chr", "bp"), skip_absent = TRUE)

lead_bp_map <- gwas[!is.na(SNP) & !is.na(bp), .(lead_snp = SNP, lead_snp_bp_from_gwas = bp)]
priority_v02 <- merge(priority_v02, lead_bp_map, by = "lead_snp", all.x = TRUE)
priority_v02[is.na(lead_snp_bp), lead_snp_bp := lead_snp_bp_from_gwas]
priority_v02[, lead_snp_bp_from_gwas := NULL]
priority_v02[is.na(lead_snp_bp), lead_snp_bp := as.integer((locus_start + locus_end) / 2)]
setorder(priority_v02, magma_rank, gene)
fwrite(priority_v02, "results/tables/priority_loci_for_smr_coloc_v0.2.tsv", sep = "\t")

windows <- priority_v02[!is.na(chr) & !is.na(lead_snp_bp), .(
  locus_id,
  gene,
  chr = as.integer(chr),
  lead_snp,
  lead_snp_bp = as.integer(lead_snp_bp),
  window_start_500kb = pmax(1L, as.integer(lead_snp_bp) - 500000L),
  window_end_500kb = as.integer(lead_snp_bp) + 500000L,
  window_start_1mb = pmax(1L, as.integer(lead_snp_bp) - 1000000L),
  window_end_1mb = as.integer(lead_snp_bp) + 1000000L,
  recommended_window = fifelse(grepl(";", locus_id), "1mb_sensitivity", "500kb_first_pass"),
  notes = fifelse(grepl(";", locus_id), "complex_locus_or_locus_group", "standard_priority_locus")
)]
windows <- unique(windows, by = c("gene", "locus_id"))
windows[, `:=`(n_gwas_snps_500kb = NA_integer_, n_gwas_snps_1mb = NA_integer_)]
setcolorder(windows, c(
  "locus_id", "gene", "chr", "lead_snp", "lead_snp_bp",
  "window_start_500kb", "window_end_500kb", "window_start_1mb", "window_end_1mb",
  "n_gwas_snps_500kb", "n_gwas_snps_1mb", "recommended_window", "notes"
))
fwrite(windows, "results/tables/priority_locus_windows.tsv", sep = "\t")

dir.create("results/coloc/gwas_slices", recursive = TRUE, showWarnings = FALSE)
readiness <- rbindlist(lapply(seq_len(nrow(windows)), function(i) {
  row <- windows[i]
  safe_gene <- gsub("[^A-Za-z0-9_.-]+", "_", row$gene)
  safe_locus <- gsub("[^A-Za-z0-9_.-]+", "_", row$locus_id)
  slice500 <- gwas[chr == row$chr & bp >= row$window_start_500kb & bp <= row$window_end_500kb]
  slice1mb <- gwas[chr == row$chr & bp >= row$window_start_1mb & bp <= row$window_end_1mb]
  windows[i, `:=`(n_gwas_snps_500kb = nrow(slice500), n_gwas_snps_1mb = nrow(slice1mb))]
  setnames(slice500, c("chr", "bp"), c("CHR", "BP"), skip_absent = TRUE)
  setnames(slice1mb, c("chr", "bp"), c("CHR", "BP"), skip_absent = TRUE)
  out500 <- file.path("results/coloc/gwas_slices", paste0(safe_locus, "_", safe_gene, "_gwas_500kb.tsv.gz"))
  out1mb <- file.path("results/coloc/gwas_slices", paste0(safe_locus, "_", safe_gene, "_gwas_1mb.tsv.gz"))
  fwrite(slice500[, .(SNP, CHR, BP, EA, NEA, BETA, SE, P, EAF, N)], out500, sep = "\t")
  fwrite(slice1mb[, .(SNP, CHR, BP, EA, NEA, BETA, SE, P, EAF, N)], out1mb, sep = "\t")
  data.table(
    locus_id = row$locus_id,
    gene = row$gene,
    gwas_slice_500kb = out500,
    gwas_slice_1mb = out1mb,
    gwas_slice_exists = file.exists(out500) & file.exists(out1mb),
    n_gwas_snps_500kb = nrow(slice500),
    n_gwas_snps_1mb = nrow(slice1mb),
    has_beta_se = all(c("BETA", "SE") %in% names(slice500)),
    has_eaf = "EAF" %in% names(slice500),
    has_n = "N" %in% names(slice500),
    eqtl_resource_available = FALSE,
    ready_for_coloc = FALSE,
    notes = "GWAS slice prepared; eQTL resource missing, so coloc is not ready."
  )
}))
fwrite(windows, "results/tables/priority_locus_windows.tsv", sep = "\t")
fwrite(readiness, "results/tables/coloc_input_readiness.tsv", sep = "\t")

eqtl <- data.table(
  resource_name = c(
    "GTEx_v8_kidney_cortex_eQTL",
    "GTEx_v8_whole_blood_eQTL",
    "GTEx_v8_artery_aorta_eQTL",
    "GTEx_v8_artery_tibial_eQTL",
    "GTEx_v8_adipose_subcutaneous_eQTL",
    "eQTLGen_blood_eQTL",
    "GTEx_v8_multi_tissue_or_eQTL_catalogue"
  ),
  tissue = c(
    "kidney_cortex",
    "whole_blood",
    "artery_aorta",
    "artery_tibial",
    "adipose_subcutaneous",
    "blood",
    "multi_tissue"
  ),
  source = c(
    "GTEx Portal / GTEx v8",
    "GTEx Portal / GTEx v8",
    "GTEx Portal / GTEx v8",
    "GTEx Portal / GTEx v8",
    "GTEx Portal / GTEx v8",
    "eQTLGen",
    "GTEx Portal or eQTL Catalogue"
  ),
  data_type = "cis_eQTL_summary_or_SMR_ready",
  genome_build = c("expected_GRCh38_or_lifted_GRCh37", "expected_GRCh38_or_lifted_GRCh37", "expected_GRCh38_or_lifted_GRCh37", "expected_GRCh38_or_lifted_GRCh37", "expected_GRCh38_or_lifted_GRCh37", "resource_dependent", "resource_dependent"),
  format = c("SMR BESD or coloc-ready summary", "SMR BESD or coloc-ready summary", "SMR BESD or coloc-ready summary", "SMR BESD or coloc-ready summary", "SMR BESD or coloc-ready summary", "SMR BESD or coloc-ready summary", "resource_dependent"),
  local_path = c(
    "external/eqtl/gtex_v8/kidney_cortex/",
    "external/eqtl/gtex_v8/whole_blood/",
    "external/eqtl/gtex_v8/artery_aorta/",
    "external/eqtl/gtex_v8/artery_tibial/",
    "external/eqtl/gtex_v8/adipose_subcutaneous/",
    "external/eqtl/eqtlgen/blood/",
    "external/eqtl/multi_tissue/"
  ),
  required_for = c("P1/P2 kidney and TAL candidates", "systemic context", "vascular context", "vascular context", "metabolic context", "blood eQTL sensitivity", "multi-tissue sensitivity"),
  status = "missing",
  download_attempted = FALSE,
  download_status = "not_attempted",
  blocking_level = c("critical", "important", "important", "important", "optional", "optional", "optional"),
  notes = c(
    "Priority resource for kidney/TAL-related candidate loci; papilla-specific eQTL unavailable.",
    "Systemic immune/inflammatory reference tissue.",
    "Vascular reference tissue.",
    "Vascular reference tissue.",
    "Metabolic context reference tissue.",
    "Large blood eQTL resource if accessible.",
    "Use only after documenting build, variant IDs, and gene ID convention."
  )
)
fwrite(eqtl, "results/tables/eqtl_resource_manifest.tsv", sep = "\t")
fwrite(eqtl, "results/tables/eqtl_resource_manifest_v0.2.tsv", sep = "\t")

message("wrote\tresults/tables/priority_loci_for_smr_coloc_v0.2.tsv")
message("wrote\tresults/tables/priority_locus_windows.tsv")
message("wrote\tresults/tables/coloc_input_readiness.tsv")
message("wrote\tresults/tables/eqtl_resource_manifest_v0.2.tsv")
