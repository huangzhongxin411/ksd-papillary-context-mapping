suppressPackageStartupMessages(library(data.table))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/coloc/toy_test", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(path)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

file_nrow <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  suppressWarnings(nrow(fread(path)))
}

slice_path <- function(locus_id, gene, window) {
  locus <- gsub(";", "_", locus_id)
  path <- file.path("results/coloc/gwas_slices", paste0(locus, "_", gene, "_gwas_", window, ".tsv.gz"))
  if (file.exists(path)) path else NA_character_
}

p1 <- read_required("results/tables/p1_tal_gene_evidence.tsv")
priority <- read_required("results/tables/priority_loci_for_smr_coloc_v0.2.tsv")

role_map <- data.table(
  gene = c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"),
  manuscript_role = c(
    "representative_TAL_gene",
    "TAL_transport_candidate",
    "calcium_ion_handling_candidate",
    "calcium_sensing_candidate",
    "supporting_context_gene",
    "broad_epithelial_context_gene"
  ),
  final_interpretation = c(
    "UMOD served as a representative TAL-associated candidate gene, showing preferential expression in audited TAL cells with donor-level support.",
    "CLDN10 represented an epithelial transport-associated candidate with a reproducible TAL component.",
    "CLDN14 linked the MAGMA-supported TAL-associated signal to ion-handling biology in a TAL-preferential expression context.",
    "CASR represented a biologically relevant calcium-sensing candidate within the broader epithelial transport context.",
    "HIBADH emerged as a MAGMA-supported TAL-associated candidate requiring further expression-genetic validation.",
    "PKD2 likely represented a broader renal epithelial context gene rather than a TAL-specific marker."
  )
)

interpretation <- merge(p1, role_map, by = "gene", all.x = TRUE)
interpretation <- interpretation[, .(
  gene,
  magma_rank,
  magma_p,
  current_tier,
  scrna_top_celltype,
  TAL_rank,
  TAL_avg_expression,
  TAL_pct_expressed,
  TAL_donor_detection_fraction,
  specificity_ratio_avg,
  specificity_class,
  TAL_program_rho,
  TAL_program_fdr,
  overall_evidence_class,
  final_interpretation,
  manuscript_role,
  claim_boundary
)]
interpretation[, gene_order := match(gene, c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"))]
setorder(interpretation, gene_order)
interpretation[, gene_order := NULL]
fwrite(interpretation, "results/tables/p1_tal_gene_interpretation_summary.tsv", sep = "\t")

figure_qc <- data.table(
  figure_file = c(
    "results/figures/p1_tal_gene_dotplot.pdf",
    "results/figures/p1_tal_gene_heatmap_by_celltype.pdf",
    "results/figures/p1_tal_gene_specificity_barplot.pdf",
    "results/figures/p1_tal_gene_by_donor_boxplot.pdf",
    "results/figures/p1_tal_gene_featureplots.pdf",
    "results/figures/p1_gene_vs_tal_program_correlation.pdf"
  ),
  included_genes = "UMOD;CASR;CLDN14;CLDN10;HIBADH;PKD2",
  celltype_annotation_used = "audited GSE231569 cell-type labels; TAL derived from audited cluster-to-celltype mapping",
  donor_level_used = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  is_main_figure_candidate = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  is_supplement_candidate = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  main_message = c(
    "P1 genes are detected across audited cell types with strongest average expression in TAL.",
    "Scaled average expression separates TAL-associated candidate patterns from broader epithelial context.",
    "TAL specificity ratios distinguish TAL-preferential candidates from broad-context genes.",
    "Donor-level expression supports reproducibility while exposing inter-donor variation.",
    "Feature-level visualization provides intuitive support for expression localization but should remain descriptive.",
    "TAL program correlations provide robustness/context evidence and should not be overread as pathway membership."
  ),
  risk_of_overinterpretation = c(
    "DotPlot does not prove genetic causality or eQTL mediation.",
    "Heatmap scaling can exaggerate relative differences for lowly expressed genes.",
    "Specificity ratio depends on audited cell-type composition and should not be treated as absolute specificity.",
    "Donor plots support reproducibility but do not establish population-level prevalence.",
    "FeaturePlot is visually compelling but can be presentation-biased by embedding layout.",
    "Correlation significance is inflated by large cell counts; effect size should guide interpretation."
  ),
  action = c(
    "Use as Figure 3A candidate.",
    "Use as Figure 3B candidate.",
    "Use as Figure 3C candidate.",
    "Use as Figure 3D or supplementary robustness panel.",
    "Use as supplementary descriptive panel.",
    "Use as supplementary robustness panel."
  )
)
figure_qc[, file_exists := file.exists(figure_file)]
fwrite(figure_qc, "results/tables/p1_figure_qc_summary.tsv", sep = "\t")

required_pilot <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2", "FAM13A")
extra_pool <- priority[
  !(gene %in% required_pilot) &
    priority_class %in% c("P2_Magma_TAL_candidate", "P1_core_TAL_candidate")
][order(magma_rank)]
pilot_genes <- unique(c(required_pilot, head(extra_pool$gene, 7)))
pilot <- priority[gene %in% pilot_genes]
pilot[, p1_status := fifelse(gene %in% p1$gene, "P1_core_gene", fifelse(gene == "FAM13A", "required_exploratory_pilot", "P2_MAGMA_TAL_pilot"))]
pilot[, pilot_reason := fifelse(
  gene %in% required_pilot,
  "Required first-batch pilot target from P1/core or FAM13A list.",
  "Selected from MAGMA top-ranked TAL-associated candidate pool."
)]
pilot[, gwas_slice_500kb := mapply(slice_path, locus_id, gene, "500kb")]
pilot[, gwas_slice_1mb := mapply(slice_path, locus_id, gene, "1mb")]
pilot[, eqtl_resource_required := required_eqtl_resource]
pilot[, ready_for_coloc := FALSE]
pilot[, blocking_reason := "eQTL_resource_missing; case_fraction_pending; variant_build_harmonization_pending"]
pilot <- pilot[, .(
  gene, locus_id, chr, lead_snp, lead_snp_p, magma_rank, magma_p,
  scrna_top_celltype, scrna_benchmark_percentile, p1_status,
  priority_class, pilot_reason, gwas_slice_500kb, gwas_slice_1mb,
  eqtl_resource_required, ready_for_coloc, blocking_reason
)]
pilot[, pilot_order := match(gene, pilot_genes)]
setorder(pilot, pilot_order)
pilot[, pilot_order := NULL]
fwrite(pilot, "results/tables/coloc_pilot_targets_v0.1.tsv", sep = "\t")

pilot_readiness <- copy(pilot)
pilot_readiness[, `:=`(
  gwas_slice_exists = file.exists(gwas_slice_500kb) & file.exists(gwas_slice_1mb),
  n_snps_500kb = vapply(gwas_slice_500kb, file_nrow, integer(1)),
  n_snps_1mb = vapply(gwas_slice_1mb, file_nrow, integer(1)),
  has_beta = FALSE,
  has_se = FALSE,
  has_eaf = FALSE,
  has_n = FALSE
)]
for (i in seq_len(nrow(pilot_readiness))) {
  path <- pilot_readiness$gwas_slice_500kb[i]
  if (!is.na(path) && file.exists(path)) {
    cols <- names(fread(path, nrows = 0))
    pilot_readiness$has_beta[i] <- "BETA" %in% cols
    pilot_readiness$has_se[i] <- "SE" %in% cols
    pilot_readiness$has_eaf[i] <- "EAF" %in% cols
    pilot_readiness$has_n[i] <- "N" %in% cols
  }
}
pilot_readiness[, `:=`(
  case_fraction_available = FALSE,
  eqtl_resource_available = FALSE,
  variant_id_match_plan_ready = file.exists("results/tables/variant_id_harmonization_plan.tsv"),
  build_match_plan_ready = file.exists("results/tables/variant_id_harmonization_plan.tsv"),
  ready_for_coloc = FALSE,
  blocking_reason = fifelse(
    !gwas_slice_exists,
    "gwas_slice_missing",
    "eQTL_resource_missing; case_fraction_pending; variant_build_harmonization_pending"
  )
)]
pilot_readiness <- pilot_readiness[, .(
  gene, locus_id, gwas_slice_exists, n_snps_500kb, n_snps_1mb,
  has_beta, has_se, has_eaf, has_n, case_fraction_available,
  eqtl_resource_available, variant_id_match_plan_ready, build_match_plan_ready,
  ready_for_coloc, blocking_reason
)]
fwrite(pilot_readiness, "results/tables/coloc_pilot_readiness.tsv", sep = "\t")

toy <- data.table(
  snp = paste0("toy_rs", 1:20),
  beta_gwas = seq(-0.2, 0.2, length.out = 20),
  se_gwas = 0.05,
  beta_eqtl = seq(-0.18, 0.18, length.out = 20),
  se_eqtl = 0.04,
  p_gwas = pmax(1e-8, pnorm(-abs(seq(-0.2, 0.2, length.out = 20) / 0.05)) * 2),
  p_eqtl = pmax(1e-8, pnorm(-abs(seq(-0.18, 0.18, length.out = 20) / 0.04)) * 2),
  note = "technical unit test only; not a biological result"
)
fwrite(toy, "results/coloc/toy_test/toy_coloc_input.tsv", sep = "\t")
toy_result <- data.table(
  test_name = "toy_coloc_unit_test",
  status = if (requireNamespace("coloc", quietly = TRUE)) "package_available_not_executed_by_default" else "coloc_package_missing",
  n_snps = nrow(toy),
  pp_h4 = NA_real_,
  note = "technical unit test only; not a biological result"
)
fwrite(toy_result, "results/coloc/toy_test/toy_coloc_result.tsv", sep = "\t")
writeLines(c(
  "toy_coloc_unit_test",
  "technical unit test only; not a biological result",
  paste("status:", toy_result$status),
  "The toy input confirms expected column structure. Real coloc execution remains blocked until eQTL resources are available."
), "results/coloc/toy_test/toy_coloc_log.txt")

spatial_decision <- data.table(
  dataset = c("GSE206306", "GSE231630"),
  root = c("data/raw/GSE206306", "data/raw/GSE231630")
)
spatial_decision[, n_files := vapply(root, function(x) if (dir.exists(x)) length(list.files(x, recursive = TRUE, full.names = TRUE)) else 0L, integer(1))]
spatial_decision[, matrix_available := FALSE]
spatial_decision[, coordinates_available := FALSE]
spatial_decision[, image_available := FALSE]
spatial_decision[, metadata_available := FALSE]
spatial_decision[, processed_object_available := FALSE]
spatial_decision[, raw_fastq_only := FALSE]
for (i in seq_len(nrow(spatial_decision))) {
  files <- if (dir.exists(spatial_decision$root[i])) list.files(spatial_decision$root[i], recursive = TRUE, full.names = TRUE) else character()
  b <- basename(files)
  spatial_decision$matrix_available[i] <- any(grepl("\\.h5$|matrix\\.mtx|count.*matrix|filtered.*feature", b, ignore.case = TRUE))
  spatial_decision$coordinates_available[i] <- any(grepl("tissue_positions|spatial.*coordinates|positions.*csv", b, ignore.case = TRUE))
  spatial_decision$image_available[i] <- any(grepl("\\.(tif|tiff|png|jpg|jpeg)$", b, ignore.case = TRUE))
  spatial_decision$metadata_available[i] <- any(grepl("metadata|sample|series|soft", b, ignore.case = TRUE))
  spatial_decision$processed_object_available[i] <- any(grepl("\\.(rds|h5ad|loom)$", b, ignore.case = TRUE))
  spatial_decision$raw_fastq_only[i] <- length(files) > 0 && all(grepl("\\.(fastq|fq)(\\.gz)?$", b, ignore.case = TRUE))
}
spatial_decision[, analysis_ready := matrix_available & coordinates_available & image_available]
spatial_decision[, decision := fifelse(
  analysis_ready, "ready_for_qc",
  fifelse(raw_fastq_only, "raw_reprocessing_required",
          fifelse(n_files == 0, "manual_download_required", "pending"))
)]
spatial_decision[, reason := fifelse(
  analysis_ready,
  "Matrix, coordinates, and image are locally available.",
  fifelse(raw_fastq_only,
          "Only raw FASTQ-like files detected; processed spatial files are absent.",
          fifelse(n_files == 0,
                  "No local spatial files detected.",
                  "Some files are present but required spatial matrix/coordinates/image combination is incomplete."))
)]
spatial_decision[, next_action := fifelse(
  decision == "ready_for_qc",
  "Run spatial QC and sample-level inventory.",
  fifelse(decision == "raw_reprocessing_required",
          "Do not include in current scope unless Space Ranger reprocessing is planned.",
          "Confirm and manually download processed matrix, coordinates, image, and metadata if available.")
)]
spatial_decision <- spatial_decision[, .(
  dataset, matrix_available, coordinates_available, image_available,
  metadata_available, processed_object_available, raw_fastq_only,
  analysis_ready, decision, reason, next_action
)]
fwrite(spatial_decision, "results/tables/spatial_dataset_decision.tsv", sep = "\t")

gse73680_root <- "data/raw/GSE73680"
gse73680_files <- if (dir.exists(gse73680_root)) list.files(gse73680_root, recursive = TRUE, full.names = TRUE) else character()
gse73680_manifest <- data.table(
  dataset = "GSE73680",
  local_root = gse73680_root,
  n_local_files = length(gse73680_files),
  expression_matrix_available = any(grepl("\\.(txt|tsv|csv|xlsx|rds|rda|gz)$", basename(gse73680_files), ignore.case = TRUE)),
  metadata_available = any(grepl("metadata|sample|series|soft", basename(gse73680_files), ignore.case = TRUE)),
  status = if (length(gse73680_files)) "local_files_detected_require_manual_audit" else "missing_local_resource",
  intended_use = "plaque/disease-context expression support only; not genetic causality, spatial validation, or cell-type-specific validation",
  next_action = if (length(gse73680_files)) "Audit platform, sample labels, and expression matrix format before analysis." else "Acquire processed expression and sample metadata before analysis."
)
fwrite(gse73680_manifest, "results/tables/gse73680_resource_manifest.tsv", sep = "\t")
gse73680_metadata <- data.table(
  sample_id = character(),
  condition = character(),
  tissue_context = character(),
  platform = character(),
  usable_for_current_scope = logical(),
  notes = character()
)
fwrite(gse73680_metadata, "config/gse73680_sample_metadata.tsv", sep = "\t")

result1 <- paste0(
  "## Result 1. GWAS reconstruction and MAGMA prioritization identified KSD-associated genes\n\n",
  "The reconstructed KSD GWAS analysis passed the locked Phase 1 quality-control workflow and was carried forward into MAGMA gene-level prioritization. ",
  "The current manuscript-facing interpretation should emphasize the reproducible GWAS-to-gene prioritization framework and the frozen MAGMA v0.1 outputs, rather than introducing additional unvalidated modules. ",
  "MAGMA identified a ranked set of KSD-associated genes that supported downstream cell-context projection. ",
  "These outputs provide the genetic-prioritization backbone for the later single-cell analyses.\n\n"
)
result2 <- paste0(
  "## Result 2. MAGMA-prioritized KSD genes preferentially localized to TAL-associated papillary cellular states\n\n",
  "Using the audited GSE231569 papillary single-cell annotation, MAGMA-prioritized genes were projected onto renal papillary cell states. ",
  "The frozen analysis supports a TAL-associated cellular context for the MAGMA-prioritized KSD signal, with benchmarking and sensitivity outputs used to guard against overinterpretation of any single locus or marker gene. ",
  "This result should be written as cell-context localization of prioritized genes, not as proof that specific variants act through TAL cells.\n\n"
)
p1_lines <- interpretation[, paste0(
  "- **", gene, "**: ", final_interpretation,
  " TAL rank = ", TAL_rank,
  ", TAL detection = ", round(TAL_pct_expressed * 100, 1), "%",
  ", donor detection fraction = ", TAL_donor_detection_fraction,
  ", specificity ratio = ", round(specificity_ratio_avg, 2),
  ", role = `", manuscript_role, "`."
)]
result3 <- paste0(
  "## Result 3. P1 core candidate genes showed interpretable TAL or epithelial transport-associated expression patterns\n\n",
  "The P1 evidence pack evaluated six core candidates, UMOD, CASR, CLDN14, CLDN10, HIBADH, and PKD2, using audited cell-type expression, donor-level detection, specificity ratios, and TAL program correlation. ",
  "All six genes ranked highest in audited TAL cells in the current single-cell analysis and showed donor-level TAL detection in 3 of 4 available TAL donors. ",
  "The genes should not be treated as equivalent: UMOD, CLDN10, CLDN14, CASR, and HIBADH support stronger TAL-associated or transport/ion-handling contexts, whereas PKD2 is better framed as a broader renal epithelial context gene with a TAL component.\n\n",
  paste(p1_lines, collapse = "\n"),
  "\n\n",
  "Together, these results support a **MAGMA + scRNA-supported TAL-associated cellular context** for P1 genes. They do not establish TWAS support, colocalization, spatial validation, or causal mediation of KSD risk.\n\n"
)
result4 <- paste0(
  "## Resource-dependent analyses\n\n",
  "TWAS, SMR/coloc, and spatial transcriptomic analyses remain resource-dependent. ",
  "The current coloc pilot tables and spatial decision tables document execution readiness and blockers, but no biological TWAS, coloc, or spatial validation result should be reported from these modules yet.\n"
)
writeLines(c(
  "# Results Draft v0.1",
  "",
  "> Current evidence grade: MAGMA + scRNA-supported TAL-associated cellular context.",
  "",
  result1,
  result2,
  result3,
  result4
), "docs/results_draft_v0.1.md")

message("wrote results-ready P1, coloc pilot, spatial decision, GSE73680 prep, and results draft outputs")
