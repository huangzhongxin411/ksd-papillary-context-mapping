suppressPackageStartupMessages(library(data.table))

priority <- fread("results/tables/priority_loci_for_smr_coloc_v0.2.tsv")
readiness <- fread("results/tables/coloc_input_readiness.tsv")

p1 <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
pilot <- rbindlist(list(
  priority[gene %in% p1],
  priority[priority_class == "P2_Magma_TAL_candidate"][order(magma_rank)][1:8]
), fill = TRUE)
pilot <- unique(pilot[!is.na(gene)], by = "gene")
pilot[, pilot_priority := fifelse(gene %in% p1, "first_priority_P1_core", "second_priority_P2_Magma_TAL")]
pilot[, pilot_status := "waiting_for_eqtl_resource"]
pilot[, notes := "GWAS slice prepared; run coloc only after tissue-matched eQTL resource is available and variant IDs/build are harmonized."]
setcolorder(pilot, c("gene", "pilot_priority", setdiff(names(pilot), c("gene", "pilot_priority"))))
fwrite(pilot, "results/tables/coloc_pilot_targets_v0.1.tsv", sep = "\t")

upgraded <- merge(
  readiness,
  priority[, .(gene, priority_class, twas_status)],
  by = "gene",
  all.x = TRUE
)
upgraded[, `:=`(
  tissue = "kidney_cortex_first_priority",
  gwas_has_beta = has_beta_se,
  gwas_has_se = has_beta_se,
  gwas_has_eaf = has_eaf,
  gwas_has_n = has_n,
  gwas_case_fraction_ready = FALSE,
  eqtl_resource_exists = FALSE,
  eqtl_has_beta = FALSE,
  eqtl_has_se_or_p = FALSE,
  eqtl_has_maf = FALSE,
  variant_id_matchable = FALSE,
  build_matchable = FALSE,
  blocking_reason = "eQTL resource missing; case fraction not yet documented; variant/build harmonization pending"
)]
fwrite(upgraded, "results/tables/coloc_input_readiness.tsv", sep = "\t")

variant_plan <- data.table(
  resource = c(
    "KSD_2025_GWAS",
    "GTEx_v8_kidney_cortex_eQTL",
    "GTEx_v8_whole_blood_eQTL",
    "GTEx_v8_artery_eQTL",
    "eQTLGen_blood_eQTL"
  ),
  variant_style = c("rsID_with_CHR_BP_EA_NEA", "pending", "pending", "pending", "pending"),
  build = c("GWAS cleaned coordinates; MAGMA run used GRCh37-compatible references", "pending", "pending", "pending", "pending"),
  needs_liftover = c("depends_on_eQTL_resource_build", "pending", "pending", "pending", "pending"),
  needs_rsid_mapping = c("no_for_current_GWAS", "pending", "pending", "pending", "pending"),
  notes = c(
    "GWAS slices include SNP, CHR, BP, EA, NEA, BETA, SE, P, EAF, N.",
    "Check whether eQTL uses rsID or chr_pos_ref_alt_b37/b38 before coloc.",
    "Check whether eQTL uses rsID or chr_pos_ref_alt_b37/b38 before coloc.",
    "Check whether eQTL uses rsID or chr_pos_ref_alt_b37/b38 before coloc.",
    "Check build and allele orientation before using as blood sensitivity."
  )
)
fwrite(variant_plan, "results/tables/variant_id_harmonization_plan.tsv", sep = "\t")

message("wrote\tresults/tables/coloc_pilot_targets_v0.1.tsv")
message("wrote\tresults/tables/coloc_input_readiness.tsv")
message("wrote\tresults/tables/variant_id_harmonization_plan.tsv")
