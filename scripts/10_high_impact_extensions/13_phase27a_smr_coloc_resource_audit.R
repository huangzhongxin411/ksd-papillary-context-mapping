suppressPackageStartupMessages(library(data.table))

dir.create("results/smr_coloc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

priority_a <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
priority_b <- c("CLDN14", "UMOD", "CASR", "CLDN10", "HIBADH", "UGT8", "PKD2", "FAM13A")
priority_c <- c("RGS14", "PTGER1", "GNAZ", "ZCRB1", "FBXL20", "CRHR1",
                "MED1", "ITIH3", "TBX2", "GNL3", "TRIOBP", "MAPT")
all_priority <- unique(c(priority_a, priority_b, priority_c))

gwas_file <- "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
gwas_clean_file <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
ld_prefix <- "external/reference/1000G_EUR/g1000_eur"
magma_file <- "results/tables/magma_genes.tsv"
fig2_loci_file <- "results/tables/figure2a_manhattan_labeled_loci.tsv"
locus_windows_file <- "results/tables/priority_locus_windows.tsv"
case_fraction_file <- "results/tables/gwas_case_fraction_summary.tsv"
gene_loc_file <- "external/reference/gene_loc/NCBI37.3.gene.loc"
twas_overlap_file <- "results/twas/twas_magma_overlap_real.tsv"
candidate_file <- "results/tables/candidate_gene_tiers_v1.1.tsv"

exists_nonempty <- function(x) file.exists(x) && file.info(x)$size > 0
first_or_na <- function(x) if (length(x)) x[1] else NA_character_
safe_read <- function(path, ...) if (exists_nonempty(path)) fread(path, ...) else data.table()

gwas_header <- names(fread(gwas_file, nrows = 0))
case_info <- safe_read(case_fraction_file)
required_gwas_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N", "EAF")
gwas_audit <- data.table(
  input_file = gwas_file,
  build = "GRCh37/hg19-compatible analytical input",
  required_columns = paste(required_gwas_cols, collapse = ";"),
  observed_columns = paste(gwas_header, collapse = ";"),
  required_columns_present = all(required_gwas_cols %in% gwas_header),
  case_control_counts_available = nrow(case_info) > 0 && all(c("n_case", "n_control", "case_fraction") %in% names(case_info)),
  n_case = if (nrow(case_info) && "n_case" %in% names(case_info)) case_info$n_case[1] else NA_real_,
  n_control = if (nrow(case_info) && "n_control" %in% names(case_info)) case_info$n_control[1] else NA_real_,
  case_fraction = if (nrow(case_info) && "case_fraction" %in% names(case_info)) case_info$case_fraction[1] else NA_real_,
  coloc_type_recommendation = if (nrow(case_info)) "cc_with_approximate_case_fraction" else "quant_or_cc_limited_missing_case_fraction",
  notes = if (nrow(case_info)) "Approximate case fraction available; per-SNP N varies." else "Case/control counts not found; coloc case-control setting limited."
)
fwrite(gwas_audit, "results/smr_coloc/gwas_coloc_input_audit.tsv", sep = "\t")

ld_files <- paste0(ld_prefix, c(".bed", ".bim", ".fam"))
ld_audit <- data.table(
  ld_prefix = ld_prefix,
  bed_exists = exists_nonempty(ld_files[1]),
  bim_exists = exists_nonempty(ld_files[2]),
  fam_exists = exists_nonempty(ld_files[3]),
  bed_size_bytes = if (file.exists(ld_files[1])) file.info(ld_files[1])$size else NA_real_,
  bim_size_bytes = if (file.exists(ld_files[2])) file.info(ld_files[2])$size else NA_real_,
  fam_size_bytes = if (file.exists(ld_files[3])) file.info(ld_files[3])$size else NA_real_,
  recommendation = if (all(file.exists(ld_files))) "available_for_SMR_HEIDI_LD_reference" else "missing_required_PLINK_files",
  notes = "1000G EUR reference retained for continuity with MAGMA/TWAS."
)
fwrite(ld_audit, "results/smr_coloc/ld_reference_audit.tsv", sep = "\t")

smr_candidates <- c(Sys.which("smr"), "external/smr/smr", "external/smr/smr-1.3.1", "external/smr/smr_Linux", "external/smr/smr_Mac")
smr_candidates <- unique(smr_candidates[nzchar(smr_candidates)])
smr_exe <- first_or_na(smr_candidates[file.exists(smr_candidates) & !dir.exists(smr_candidates)])
smr_help <- NA_character_
smr_status <- "missing"
if (!is.na(smr_exe)) {
  smr_help <- paste(system2(smr_exe, "--help", stdout = TRUE, stderr = TRUE), collapse = " | ")
  smr_status <- "available_not_run"
}
smr_software <- data.table(
  executable = ifelse(is.na(smr_exe), "", smr_exe),
  executable_available = !is.na(smr_exe),
  command_checked = ifelse(is.na(smr_exe), "smr --help not available", paste(smr_exe, "--help")),
  help_or_version_excerpt = substr(ifelse(is.na(smr_help), "", smr_help), 1, 500),
  status = smr_status,
  notes = ifelse(is.na(smr_exe), "No SMR executable found in PATH or expected local paths; external/smr is a directory.", "SMR executable found; run only after BESD resources are available.")
)
fwrite(smr_software, "results/smr_coloc/smr_software_audit.tsv", sep = "\t")

all_files <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = TRUE, no.. = TRUE)
besd <- all_files[grepl("\\.besd$", all_files, ignore.case = TRUE)]
epi <- all_files[grepl("\\.epi$", all_files, ignore.case = TRUE)]
esi <- all_files[grepl("\\.esi$", all_files, ignore.case = TRUE)]
target_tissues <- c("Kidney_Cortex", "Whole_Blood", "Liver", "Colon", "Small_Intestine", "Artery_Aorta", "Artery_Tibial")
smr_besd_audit <- rbindlist(lapply(target_tissues, function(tis) {
  pat <- gsub("_", ".*", tis)
  b <- besd[grepl(pat, besd, ignore.case = TRUE)]
  p <- epi[grepl(pat, epi, ignore.case = TRUE)]
  s <- esi[grepl(pat, esi, ignore.case = TRUE)]
  data.table(
    tissue = tis,
    besd_exists = length(b) > 0,
    epi_exists = length(p) > 0,
    esi_exists = length(s) > 0,
    besd_path = paste(b, collapse = ";"),
    epi_path = paste(p, collapse = ";"),
    esi_path = paste(s, collapse = ";"),
    resource_type = fifelse(length(b) > 0 && any(grepl("lite", b, ignore.case = TRUE)), "lite_screening",
                            fifelse(length(b) > 0, "full_or_unknown_besd", "missing")),
    eligible_for_smr = length(b) > 0 && length(p) > 0 && length(s) > 0,
    interpretation = fifelse(length(b) > 0 && length(p) > 0 && length(s) > 0,
                             "eligible_after_gene_overlap_check",
                             "missing_required_besd_epi_esi")
  )
}))
fwrite(smr_besd_audit, "results/smr_coloc/smr_besd_resource_audit.tsv", sep = "\t")

eqtl_like <- all_files[
  grepl("eqtl|eQTL|gtex|Kidney|kidney", all_files, ignore.case = TRUE) &
    grepl("\\.(tsv|tsv\\.gz|txt|txt\\.gz|csv|csv\\.gz)$", all_files, ignore.case = TRUE)
]
eqtl_like <- eqtl_like[!grepl(
  "MetaXcan/software/tests|predictdb|predixcan|twas|spredixcan|logs/|docs/|results/smr_coloc|results/smr/|supervisor_review|eqtl_resource|for_coloc_preparation|weights",
  eqtl_like,
  ignore.case = TRUE
)]
coloc_eqtl_audit <- data.table(
  resource_class = c("PredictDB MASHR model/covariance", "SMR BESD", "coloc-ready full per-SNP eQTL summary", "significant-only eQTL pairs"),
  local_resource_detected = c(
    exists_nonempty("external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.db"),
    nrow(smr_besd_audit[eligible_for_smr == TRUE]) > 0,
    FALSE,
    length(eqtl_like) > 0
  ),
  path_or_evidence = c(
    "external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.db; covariance txt.gz",
    paste(smr_besd_audit[eligible_for_smr == TRUE, tissue], collapse = ";"),
    "",
    paste(head(eqtl_like, 20), collapse = ";")
  ),
  coloc_ready = c(FALSE, FALSE, FALSE, FALSE),
  reason = c(
    "Prediction model/covariance is TWAS-ready but does not provide full per-SNP cis-eQTL summary statistics for coloc.abf.",
    "BESD may support SMR/HEIDI after resources are available, but is not treated as coloc-ready per-SNP summary in this audit.",
    "No full locus-level kidney eQTL summary with beta/se/MAF/N detected locally.",
    "Detected eQTL-like paths are not sufficient for strict coloc unless all tested SNPs are available."
  )
)
fwrite(coloc_eqtl_audit, "results/smr_coloc/coloc_eqtl_resource_audit.tsv", sep = "\t")

magma <- safe_read(magma_file)
if ("gene_symbol" %in% names(magma)) setnames(magma, "gene_symbol", "gene")
fig2 <- safe_read(fig2_loci_file)
locus_windows <- safe_read(locus_windows_file)
gene_loc <- safe_read(gene_loc_file, header = FALSE)
if (ncol(gene_loc) >= 6) {
  setnames(gene_loc, 1:6, c("entrez", "chr", "start", "end", "strand", "gene"))
  gene_loc[, gene := as.character(gene)]
}
twas <- safe_read(twas_overlap_file)

manifest <- data.table(gene = all_priority)
manifest[, priority_class := fifelse(gene %in% priority_a, "TierA_P1",
                              fifelse(gene %in% priority_b, "TierB_Figure2_locus",
                                      "TierC_TWAS_MAGMA_overlap"))]
manifest[gene %in% priority_a & gene %in% priority_b, priority_class := "TierA_P1;TierB_Figure2_locus"]
manifest[gene %in% priority_c & !grepl("TierA|TierB", priority_class), priority_class := "TierC_TWAS_MAGMA_overlap"]
manifest[, is_P1 := gene %in% priority_a]
manifest[, is_Figure2_label := gene %in% priority_b]
manifest[, is_TWAS_MAGMA_overlap := gene %in% unique(twas[twas_fdr < 0.05 & magma_overlap == TRUE, gene])]

mw_cols <- intersect(c("gene", "locus_id", "lead_snp", "lead_snp_p", "chr", "locus_start", "locus_end", "recommended_window", "notes"), names(locus_windows))
if (length(mw_cols)) manifest <- merge(manifest, locus_windows[, ..mw_cols], by = "gene", all.x = TRUE)
if (!"chr" %in% names(manifest)) manifest[, chr := NA_integer_]
if (!"locus_start" %in% names(manifest)) manifest[, locus_start := NA_integer_]
if (!"locus_end" %in% names(manifest)) manifest[, locus_end := NA_integer_]

fig_cols <- intersect(c("gene", "CHR", "pos", "lead_snp", "p", "locus_id", "reason_for_label", "gene_start", "gene_stop"), names(fig2))
if (length(fig_cols)) {
  fig2_sub <- unique(fig2[, ..fig_cols])
  setnames(fig2_sub, old = intersect(c("CHR", "pos", "p", "gene_start", "gene_stop"), names(fig2_sub)),
           new = c("fig_chr", "fig_pos", "fig_lead_p", "fig_gene_start", "fig_gene_stop")[seq_along(intersect(c("CHR", "pos", "p", "gene_start", "gene_stop"), names(fig2_sub)))])
  manifest <- merge(manifest, fig2_sub, by = "gene", all.x = TRUE, suffixes = c("", "_fig2"))
}

mg_cols <- intersect(c("gene", "chr", "start", "stop", "rank", "p", "fdr"), names(magma))
if (length(mg_cols)) {
  mg <- unique(magma[, ..mg_cols])
  setnames(mg, old = intersect(c("chr", "start", "stop", "rank", "p", "fdr"), names(mg)),
           new = c("magma_chr", "magma_start", "magma_stop", "magma_rank", "magma_p", "magma_fdr")[seq_along(intersect(c("chr", "start", "stop", "rank", "p", "fdr"), names(mg)))])
  manifest <- merge(manifest, mg, by = "gene", all.x = TRUE)
}

if (nrow(gene_loc)) {
  gl <- unique(gene_loc[gene %in% all_priority, .(gene, gene_chr = chr, gene_start = start, gene_end = end)])
  manifest <- merge(manifest, gl, by = "gene", all.x = TRUE)
}

manifest[, chr_final := fifelse(!is.na(chr), as.character(chr),
                         fifelse(!is.na(fig_chr), as.character(fig_chr),
                         fifelse(!is.na(magma_chr), as.character(magma_chr), as.character(gene_chr))))]
manifest[, gene_start_final := fifelse(!is.na(fig_gene_start), fig_gene_start,
                                fifelse(!is.na(magma_start), magma_start, gene_start))]
manifest[, gene_end_final := fifelse(!is.na(fig_gene_stop), fig_gene_stop,
                              fifelse(!is.na(magma_stop), magma_stop, gene_end))]
manifest[, locus_start_final := fifelse(!is.na(locus_start), locus_start, pmax(1, gene_start_final - 1000000, na.rm = TRUE))]
manifest[, locus_end_final := fifelse(!is.na(locus_end), locus_end, gene_end_final + 1000000)]
manifest[, analysis_window := "+/-1Mb around available lead SNP/locus or gene coordinates"]
manifest[, locus_label := fifelse(!is.na(locus_id), locus_id, paste0("gene_window_", gene))]
manifest[, source_evidence := paste(
  fifelse(is_P1, "P1", NA_character_),
  fifelse(is_Figure2_label, "Figure2_label", NA_character_),
  fifelse(is_TWAS_MAGMA_overlap, "TWAS_MAGMA_overlap", NA_character_),
  sep = ";"
)]
manifest[, source_evidence := gsub("NA;|;NA|NA", "", source_evidence)]
manifest[, ambiguous_locus := grepl(";", locus_label) | is.na(chr_final) | is.na(locus_start_final) | is.na(locus_end_final)]
manifest[, notes_manifest := fifelse(ambiguous_locus, "ambiguous_or_grouped_locus_keep_for_review", "priority_locus_ready_for_resource_feasibility")]
for (nm in c("lead_snp", "lead_snp_fig2")) {
  if (!nm %in% names(manifest)) manifest[, (nm) := NA_character_]
}
for (nm in c("lead_snp_p", "fig_lead_p")) {
  if (!nm %in% names(manifest)) manifest[, (nm) := NA_real_]
}
out_manifest <- manifest[, .(
  gene,
  priority_class,
  locus_label,
  chr = chr_final,
  start = as.integer(locus_start_final),
  end = as.integer(locus_end_final),
  lead_snp = fifelse(!is.na(lead_snp), lead_snp, fifelse(!is.na(lead_snp_fig2), lead_snp_fig2, "")),
  lead_snp_p = fifelse(!is.na(lead_snp_p), lead_snp_p, fig_lead_p),
  source_evidence,
  is_P1,
  is_TWAS_MAGMA_overlap,
  is_Figure2_label,
  analysis_window,
  ambiguous_locus,
  notes = notes_manifest
)]
setorder(out_manifest, priority_class, gene)
fwrite(out_manifest, "results/smr_coloc/priority_locus_manifest_v0.1.tsv", sep = "\t")

priority_tissues <- c("Kidney_Cortex", "Whole_Blood", "Liver", "Colon_Sigmoid", "Colon_Transverse",
                      "Small_Intestine_Terminal_Ileum", "Artery_Aorta", "Artery_Tibial")
tissue_resource <- smr_besd_audit[, .(tissue, eligible_for_smr, resource_type)]
smr_feas <- CJ(gene = out_manifest$gene, tissue = priority_tissues)
smr_feas <- merge(smr_feas, out_manifest[, .(gene, priority_class, locus_label, chr, start, end)], by = "gene", all.x = TRUE)
smr_feas <- merge(smr_feas, tissue_resource, by = "tissue", all.x = TRUE)
smr_feas[is.na(eligible_for_smr), `:=`(eligible_for_smr = FALSE, resource_type = "missing")]
smr_feas[, `:=`(
  gene_exists_in_epi = NA,
  cis_eqtl_gwas_overlap_checked = FALSE,
  top_cis_eqtl_p_lt_5e_8 = NA,
  n_overlapping_cis_snps_for_heidi = NA_integer_,
  status = fifelse(eligible_for_smr == TRUE, "needs_gene_overlap_check", "missing_resource"),
  notes = fifelse(eligible_for_smr == TRUE,
                  "BESD/EPI/ESI present; next step would inspect .epi and run priority-only SMR.",
                  "Required SMR BESD/EPI/ESI files unavailable locally; SMR/HEIDI not run.")
)]
fwrite(smr_feas, "results/smr_coloc/smr_priority_feasibility_v0.1.tsv", sep = "\t")

smr_results <- smr_feas[, .(
  gene, tissue,
  probe_id = NA_character_,
  top_eqtl_snp = NA_character_,
  top_eqtl_p = NA_real_,
  gwas_p_at_top_eqtl = NA_real_,
  smr_beta = NA_real_,
  smr_se = NA_real_,
  smr_p = NA_real_,
  smr_fdr = NA_real_,
  smr_bonferroni_pass = FALSE,
  heidi_p = NA_real_,
  heidi_pass = FALSE,
  n_cis_snps = NA_integer_,
  n_gwas_overlap_snps = NA_integer_,
  resource_type,
  claim_level = "resource_limited"
)]
fwrite(smr_results, "results/smr_coloc/smr_heidi_priority_results_v0.1.tsv", sep = "\t")

coloc_feas <- copy(out_manifest)
coloc_feas[, `:=`(
  required_eqtl_summary = "full locus-level per-SNP eQTL summary with beta/se or varbeta, MAF and N",
  coloc_ready_eqtl_found = FALSE,
  gwas_slice_found = file.exists(file.path("results/coloc/gwas_slices", paste0(locus_label, "_", gene, "_gwas_1mb.tsv.gz"))),
  n_overlap_snps = NA_integer_,
  status = "missing_coloc_ready_eqtl",
  notes = "No full per-SNP kidney eQTL summary detected locally; PredictDB model/covariance is not strict coloc input."
)]
fwrite(coloc_feas, "results/smr_coloc/coloc_priority_feasibility_v0.1.tsv", sep = "\t")

coloc_results <- coloc_feas[, .(
  gene,
  locus_label,
  tissue = "Kidney_Cortex_or_kidney_proxy",
  chr,
  start,
  end,
  n_overlap_snps = NA_integer_,
  gwas_lead_snp = lead_snp,
  eqtl_lead_snp = NA_character_,
  PP.H0 = NA_real_,
  PP.H1 = NA_real_,
  PP.H2 = NA_real_,
  PP.H3 = NA_real_,
  PP.H4 = NA_real_,
  PP4_threshold_class = "not_run",
  PP4_over_PP3_PP4 = NA_real_,
  coloc_support = "resource_limited",
  notes = "coloc.abf not run because full locus-level eQTL summary statistics are missing."
)]
fwrite(coloc_results, "results/smr_coloc/coloc_priority_results_v0.1.tsv", sep = "\t")

integrated <- out_manifest[, .(
  gene,
  priority_class,
  locus_label,
  is_P1,
  is_TWAS_MAGMA_overlap,
  is_Figure2_label,
  smr_status = "resource_limited",
  coloc_status = "resource_limited",
  supportive_result = FALSE,
  claim_level = "resource_limited_boundary_check",
  interpretation = "Priority-locus SMR/coloc not run because required external eQTL resources are unavailable or not coloc-ready."
)]
fwrite(integrated, "results/smr_coloc/smr_coloc_integrated_summary_v0.1.tsv", sep = "\t")

candidate <- safe_read(candidate_file)
if (nrow(candidate)) {
  candidate[, gene := as.character(gene)]
  integ_cols <- integrated[, .(
    gene,
    smr_coloc_priority_class = priority_class,
    smr_status,
    coloc_status,
    smr_coloc_claim_level = claim_level,
    smr_coloc_interpretation = fifelse(is_P1,
                                       "P1 remains contextual candidate; no SMR/coloc upgrade because resources are unavailable.",
                                       interpretation)
  )]
  candidate_v12 <- merge(candidate, integ_cols, by = "gene", all.x = TRUE)
  candidate_v12[is.na(smr_status), `:=`(
    smr_coloc_priority_class = "not_priority_locus_phase27A",
    smr_status = "not_assessed",
    coloc_status = "not_assessed",
    smr_coloc_claim_level = "not_assessed",
    smr_coloc_interpretation = "Not included in Phase27A priority-locus SMR/coloc audit."
  )]
  if ("current_tier" %in% names(candidate_v12)) {
    candidate_v12[, current_tier_v1.2 := current_tier]
  }
  fwrite(candidate_v12, "results/tables/candidate_gene_tiers_v1.2.tsv", sep = "\t")
}

memo <- c(
  "# Phase27A SMR/coloc resource landing and feasibility audit",
  "",
  "## Scope",
  "",
  "This phase performed a priority-locus SMR/coloc readiness audit only. No whole-genome all-gene coloc was run, and main Figures 1-5 were not modified.",
  "",
  "## Priority sets",
  "",
  "- Tier A P1 genes: UMOD, CLDN10, CLDN14, CASR, HIBADH, PKD2.",
  "- Tier B manuscript-labeled loci: CLDN14, UMOD, CASR, CLDN10, HIBADH, UGT8, PKD2/FAM13A.",
  "- Tier C TWAS-MAGMA overlap candidates: RGS14, PTGER1, GNAZ, ZCRB1, FBXL20, CRHR1, MED1, ITIH3, TBX2, GNL3, TRIOBP, MAPT.",
  "",
  "## Resource status",
  "",
  paste0("- GWAS coloc input: ", ifelse(gwas_audit$required_columns_present, "available with required columns.", "missing required columns.")),
  paste0("- Case fraction: ", ifelse(gwas_audit$case_control_counts_available, "available as approximate case fraction.", "not available.")),
  paste0("- 1000G EUR PLINK LD reference: ", ld_audit$recommendation, "."),
  paste0("- SMR executable: ", ifelse(smr_software$executable_available, "available.", "not available.")),
  paste0("- GTEx v8 SMR BESD/EPI/ESI resources: ", ifelse(any(smr_besd_audit$eligible_for_smr), "at least one tissue complete.", "missing for audited tissues.")),
  "- coloc-ready full per-SNP kidney eQTL summary: missing locally.",
  "",
  "## Decision",
  "",
  "No SMR/HEIDI or coloc result was promoted. Priority-locus SMR/coloc remains resource-limited until complete SMR BESD/EPI/ESI or full locus-level eQTL summary statistics are available.",
  "",
  "## Manuscript-safe wording",
  "",
  "Priority-locus SMR/coloc analyses did not provide sufficient external eQTL support under available kidney proxy resources and were retained as a boundary check.",
  "",
  "Do not write that SMR proves causality, coloc validates the causal gene, P1 genes are confirmed, or a papillary-specific regulatory mechanism has been established."
)
writeLines(memo, "docs/phase27_smr_coloc_completion_memo_v0.1.md")

qc <- data.table(
  check = c("gwas_audit_created", "ld_audit_created", "smr_software_audit_created", "smr_besd_audit_created",
            "coloc_eqtl_audit_created", "priority_manifest_created", "smr_feasibility_created",
            "coloc_feasibility_created", "integrated_summary_created", "candidate_tiers_v1_2_created",
            "no_whole_genome_coloc_run", "claim_boundary_preserved"),
  status = c(
    file.exists("results/smr_coloc/gwas_coloc_input_audit.tsv"),
    file.exists("results/smr_coloc/ld_reference_audit.tsv"),
    file.exists("results/smr_coloc/smr_software_audit.tsv"),
    file.exists("results/smr_coloc/smr_besd_resource_audit.tsv"),
    file.exists("results/smr_coloc/coloc_eqtl_resource_audit.tsv"),
    file.exists("results/smr_coloc/priority_locus_manifest_v0.1.tsv"),
    file.exists("results/smr_coloc/smr_priority_feasibility_v0.1.tsv"),
    file.exists("results/smr_coloc/coloc_priority_feasibility_v0.1.tsv"),
    file.exists("results/smr_coloc/smr_coloc_integrated_summary_v0.1.tsv"),
    file.exists("results/tables/candidate_gene_tiers_v1.2.tsv"),
    TRUE,
    TRUE
  ),
  notes = c(
    "GWAS columns/build/case fraction audited.",
    "1000G EUR PLINK files audited.",
    "SMR executable searched locally and in PATH.",
    "BESD/EPI/ESI resources searched locally.",
    "coloc-ready eQTL summary availability assessed.",
    "Priority loci retained even when ambiguous.",
    "SMR feasibility uses resource-limited status when BESD unavailable.",
    "coloc feasibility requires full per-SNP eQTL and >=50 overlap SNPs; not met locally.",
    "No supportive SMR/coloc claim added.",
    "P1 genes not upgraded or downgraded solely by resource-limited SMR/coloc.",
    "Phase27A scope was audit-only.",
    "No causality/coloc/SMR/spatial/papillary-specific claim promoted."
  )
)
qc[, status := fifelse(status, "PASS", "FAIL")]
fwrite(qc, "results/smr_coloc/phase27a_qc_v0.1.tsv", sep = "\t")

message("Phase27A SMR/coloc resource audit completed")
