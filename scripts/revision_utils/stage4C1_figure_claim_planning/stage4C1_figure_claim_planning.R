suppressPackageStartupMessages({
  library(data.table)
})

doc_dir <- "docs/revision/stage4C1_figure_claim_planning"
table_dir <- "results/tables/revision/stage4C1_figure_claim_planning"
figure_dir <- "results/figures/revision/stage4C1_figure_claim_planning"
log_dir <- "logs/revision/stage4C1_figure_claim_planning"
script_dir <- "scripts/revision_utils/stage4C1_figure_claim_planning"
for (d in c(doc_dir, table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "stage4C1_figure_claim_planning.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4C1 conservative figure and claim planning\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

exists_yes <- function(path) ifelse(file.exists(path), "yes", "no")
collapse0 <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = ";")
}

paths <- list(
  stage3_model = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  stage3_counts = "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv",
  stage3_exemplar = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
  stage4a_loop = "docs/revision/stage4A_scrna_audit/loop_tal_540_claim_audit.md",
  stage4a_design = "docs/revision/stage4A_scrna_audit/stage4B_statistical_design_plan.md",
  stage4b1_scores = "results/tables/revision/stage4B1_scrna_donor_level/scrna_donor_compartment_module_scores.tsv",
  stage4b1_ranks = "results/tables/revision/stage4B1_scrna_donor_level/scrna_loop_tal_within_donor_rank.tsv",
  stage4b1_consistency = "results/tables/revision/stage4B1_scrna_donor_level/scrna_loop_tal_donor_consistency_summary.tsv",
  stage4b1_loo = "results/tables/revision/stage4B1_scrna_donor_level/scrna_leave_one_donor_out_module_ranks.tsv",
  stage4b1_low = "results/tables/revision/stage4B1_scrna_donor_level/scrna_low_loop_tal_count_sensitivity.tsv",
  stage4b2_random = "results/tables/revision/stage4B2_scrna_robustness/scrna_random_set_benchmark_summary.tsv",
  stage4b2_driver = "results/tables/revision/stage4B2_scrna_robustness/scrna_known_driver_removal_sensitivity.tsv",
  stage4b2_claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
  stage4b2_evidence = "results/tables/revision/stage4B2_scrna_robustness/evidence_model_with_stage4B2_scrna_robustness_v0.1.tsv"
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) stop("Missing required Stage 4C1 input: ", paths[[nm]])
}

stage3_model <- fread(paths$stage3_model)
stage3_counts <- fread(paths$stage3_counts)
stage3_exemplar <- fread(paths$stage3_exemplar)
stage4b1_scores <- fread(paths$stage4b1_scores)
stage4b1_ranks <- fread(paths$stage4b1_ranks)
stage4b1_consistency <- fread(paths$stage4b1_consistency)
stage4b1_loo <- fread(paths$stage4b1_loo)
stage4b1_low <- fread(paths$stage4b1_low)
stage4b2_random <- fread(paths$stage4b2_random)
stage4b2_driver <- fread(paths$stage4b2_driver)
stage4b2_claim <- fread(paths$stage4b2_claim)
stage4b2_evidence <- fread(paths$stage4b2_evidence)

overview_source <- stage4b1_scores[, .(
  n_donor_compartment_rows = .N,
  n_donors = uniqueN(donor_id),
  n_modules = uniqueN(module_name),
  total_nuclei = sum(n_nuclei),
  low_cell_count_rows = sum(low_cell_count_flag == "yes")
), by = broad_compartment]
fwrite(overview_source, file.path(table_dir, "figure2_panelA_compartment_overview_source.tsv"), sep = "\t")

evidence_status <- rbindlist(list(
  data.table(evidence_component = "MAGMA/GWAS reproducibility", status = "source-backed after Stage 2B-R", support_level = "strong",
             source_file = "results/tables/revision/stage2_genetic/magma_output_audit.tsv",
             allowed_claim = "MAGMA-prioritized genes from reproducible available outputs",
             disallowed_claim = "causal gene mapping; ancestry-generalizable fine mapping",
             figure_use = "Figure 3 evidence model background", manuscript_use = "Methods and Results reproducibility boundary",
             notes = "MAGMA remains EUR-LD-reference-based prioritization."),
  data.table(evidence_component = "LD-reference limitation", status = "EUR LD reference used", support_level = "limitation",
             source_file = "docs/revision/stage2_genetic/ld_reference_mismatch_audit.md",
             allowed_claim = "EUR-LD-reference-based gene prioritization",
             disallowed_claim = "trans-ancestry fine mapping",
             figure_use = "caption/legend boundary note", manuscript_use = "Methods and Limitations",
             notes = "Do not visually imply ancestry-generalizable mapping."),
  data.table(evidence_component = "TWAS proxy status", status = "Kidney_Cortex proxy only; many one-SNP models", support_level = "supplementary_only",
             source_file = "results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv",
             allowed_claim = "proxy support where clearly labelled",
             disallowed_claim = "papilla-specific TWAS mechanism",
             figure_use = "Figure 3 supplementary evidence axis", manuscript_use = "Results boundary and Discussion limitation",
             notes = "R3/R2_R3 snRNA modules remain supplementary."),
  data.table(evidence_component = "SMR/coloc absence", status = "no claim-grade SMR/coloc support", support_level = "limitation",
             source_file = "results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv",
             allowed_claim = "not ready for SMR/coloc-supported claim",
             disallowed_claim = "SMR-supported; coloc-supported",
             figure_use = "Figure 3 evidence strip", manuscript_use = "Results and Limitations",
             notes = "Absence of resource support is not negative biological evidence."),
  data.table(evidence_component = "Stage 3R evidence model", status = paste0(nrow(stage3_model), " genes in two-axis model"), support_level = "strong",
             source_file = paths$stage3_model,
             allowed_claim = "tiered/prioritized evidence model",
             disallowed_claim = "validated or causal candidates",
             figure_use = "Figure 3 Panels A-B", manuscript_use = "Results",
             notes = "Reporting groups are mutually exclusive."),
  data.table(evidence_component = "curated exemplar panel", status = paste0(nrow(stage3_exemplar), " exemplar genes"), support_level = "supplementary_only",
             source_file = paths$stage3_exemplar,
             allowed_claim = "biological role-spectrum exemplars",
             disallowed_claim = "evidence-upgraded genes or validation set",
             figure_use = "Figure 3 Panels C-D", manuscript_use = "Results interpretive context",
             notes = "Exemplar status must not upgrade genetic evidence."),
  data.table(evidence_component = "GSE231569 Loop/TAL label and donor count", status = "540 Loop/TAL nuclei across four donors", support_level = "moderate",
             source_file = paths$stage4a_loop,
             allowed_claim = "source-supported descriptive count",
             disallowed_claim = "inferential biological replicate count",
             figure_use = "Figure 2 legend and Panel A note", manuscript_use = "Results descriptive anchor",
             notes = "Loop/TAL nuclei are imbalanced across donors."),
  data.table(evidence_component = "donor-level Loop/TAL support", status = "primary modules strong_descriptive_support in Stage 4B1", support_level = "moderate",
             source_file = paths$stage4b1_consistency,
             allowed_claim = "donor-level descriptive Loop/TAL-associated pattern",
             disallowed_claim = "cell-level enrichment",
             figure_use = "Figure 2 Panels B-C", manuscript_use = "Results",
             notes = "Primary unit is donor x broad_compartment."),
  data.table(evidence_component = "leave-one-donor-out support", status = "primary modules retained", support_level = "moderate",
             source_file = paths$stage4b1_loo,
             allowed_claim = "not obviously driven by a single donor",
             disallowed_claim = "replicated in independent cohorts",
             figure_use = "Figure 2 Panel D", manuscript_use = "Results robustness sentence",
             notes = "Only four donors remain a limitation."),
  data.table(evidence_component = "low-cell-count sensitivity", status = "low Loop/TAL-count donor exclusion unchanged for primary modules", support_level = "moderate",
             source_file = paths$stage4b1_low,
             allowed_claim = "low-count donor did not change primary module rank",
             disallowed_claim = "low-cell-count issue resolved",
             figure_use = "legend/supplement or Panel D annotation", manuscript_use = "Limitations and Results caveat",
             notes = "GSM7290914 has four Loop/TAL nuclei."),
  data.table(evidence_component = "matched random benchmark", status = "partial_support for 4/4 primary modules", support_level = "partial",
             source_file = paths$stage4b2_random,
             allowed_claim = "partial support beyond expression/detection-matched random expectations",
             disallowed_claim = "robust beyond random expectation",
             figure_use = "Figure 2 Panel E", manuscript_use = "Results qualifying sentence",
             notes = "Rank metrics saturate; delta supports partial benchmark signal."),
  data.table(evidence_component = "known-driver removal", status = "robust for 4/4 primary modules", support_level = "moderate",
             source_file = paths$stage4b2_driver,
             allowed_claim = "robust to known-driver removal",
             disallowed_claim = "not driven by any biological prior in causal sense",
             figure_use = "Figure 2 Panel F", manuscript_use = "Results robustness sentence",
             notes = "Descriptive module-level sensitivity."),
  data.table(evidence_component = "TWAS-proxy module status", status = "weak/inconsistent and supplementary", support_level = "supplementary_only",
             source_file = paths$stage4b2_claim,
             allowed_claim = "supplementary TWAS-proxy context",
             disallowed_claim = "TWAS-supported Loop/TAL mechanism",
             figure_use = "Figure 3/supplement only", manuscript_use = "Results boundary note",
             notes = "Do not use R3/R2_R3 for main Loop/TAL evidence."),
  data.table(evidence_component = "final Loop/TAL claim level", status = "moderate_main_claim_allowed", support_level = "moderate",
             source_file = paths$stage4b2_claim,
             allowed_claim = "MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal.",
             disallowed_claim = "strong enrichment; causal cell type; plaque nucleation site",
             figure_use = "Figure 2 Panel G and legend", manuscript_use = "main Results sentence",
             notes = "Moderate context-mapping evidence, not strong mechanistic proof.")
), fill = TRUE)
fwrite(evidence_status, file.path(table_dir, "stage4C1_integrated_evidence_status.tsv"), sep = "\t")

fig2_panels <- data.table(
  panel_id = LETTERS[1:7],
  panel_title = c(
    "Audited GSE231569 compartment overview",
    "Donor x compartment module scores",
    "Loop/TAL within-donor ranks",
    "Leave-one-donor-out robustness",
    "Matched random-set benchmark",
    "Known-driver removal sensitivity",
    "Claim-decision strip"
  ),
  main_message = c(
    "The snRNA atlas contains source-supported Loop/TAL labels across four donors, with imbalanced Loop/TAL nuclei.",
    "Primary MAGMA modules show donor-level descriptive Loop/TAL-associated scores.",
    "Loop/TAL ranks top across donors for primary modules, but this remains descriptive.",
    "Primary-module Loop/TAL support is not lost when excluding one donor at a time.",
    "Matched random sets provide partial support, mainly through Loop/TAL delta rather than rank metrics.",
    "Primary-module Loop/TAL patterns remain robust after removing known transport/ion-handling drivers.",
    "The final allowable claim is moderate, not strong."
  ),
  source_data_file = c(
    file.path(table_dir, "figure2_panelA_compartment_overview_source.tsv"),
    paths$stage4b1_scores,
    paths$stage4b1_ranks,
    paths$stage4b1_loo,
    paths$stage4b2_random,
    paths$stage4b2_driver,
    paths$stage4b2_claim
  ),
  plot_type = c("compact table or UMAP overview in Stage 4C2", "heatmap", "rank dot/strip plot", "leave-one-donor rank tile plot", "benchmark lollipop/bar with partial-support annotation", "driver-removal tile plot", "claim boundary strip"),
  allowed_interpretation = c(
    "audited descriptive atlas context",
    "donor-level descriptive module localization",
    "within-donor Loop/TAL ranking summary",
    "single-donor exclusion robustness",
    "partial support beyond matched random expectations",
    "robustness to known-driver removal",
    "moderate context-mapping claim"
  ),
  forbidden_interpretation = c(
    "cell-level inferential enrichment",
    "cell-level P value or causal localization",
    "strong enrichment beyond random expectation",
    "independent replication",
    "robust beyond-random enrichment",
    "proof that drivers are irrelevant or non-causal",
    "strong main claim or causal cell type"
  ),
  priority = c("required", "required", "required", "required", "required", "required", "required"),
  notes = c(
    "If UMAP is used in 4C2, export coordinates from the Seurat object as source data.",
    "Anchor panel for Results paragraph.",
    "Use donor-level ranks, not pooled nuclei.",
    "Show all four excluded-donor scenarios.",
    "Visually label partial_support and rank-metric saturation.",
    "Show panel removals more prominently than single-gene removals.",
    "Use wording from claim decision table."
  )
)

fig3_panels <- data.table(
  panel_id = LETTERS[1:5],
  panel_title = c(
    "Two-axis evidence model",
    "Reporting group counts",
    "Curated exemplar role spectrum",
    "Exemplar evidence strip",
    "Allowed versus disallowed claim boundary"
  ),
  main_message = c(
    "Candidate genes are organized by genetic priority and TWAS proxy evidence, not by causal validation.",
    "Reporting groups are mutually exclusive and quantify the evidence model.",
    "Curated exemplars illustrate biological roles but do not upgrade evidence.",
    "Exemplars are MAGMA-prioritized where applicable, with TWAS/SMR/coloc boundaries shown explicitly.",
    "The figure defines what claims are and are not allowed."
  ),
  source_data_file = c(
    paths$stage3_model,
    paths$stage3_counts,
    paths$stage3_exemplar,
    paths$stage3_exemplar,
    file.path(table_dir, "snRNA_claim_wording_decision_table.tsv")
  ),
  plot_type = c("2D evidence matrix schematic", "bar chart", "role-spectrum strip", "evidence strip heatmap", "claim boundary table/strip"),
  allowed_interpretation = c(
    "structured candidate evidence model",
    "distribution of reporting groups",
    "biological interpretation aid",
    "transparent evidence boundary per exemplar",
    "claim-boundary guide"
  ),
  forbidden_interpretation = c(
    "causal ranking",
    "quality score or causal tier",
    "cherry-picked validation",
    "SMR/coloc support or papilla-specific TWAS evidence",
    "therapeutic target or validated gene claim"
  ),
  priority = c("required", "required", "required", "required", "required"),
  notes = c(
    "This is a conceptual/source-data-backed panel, not a mechanistic proof.",
    "Use exact Stage 3R counts.",
    "Keep separate from Figure 2 Loop/TAL proof.",
    "Show absent SMR/coloc as boundary, not as negative biology.",
    "Useful as figure-end boundary panel."
  )
)

write_blueprint <- function(path, title, claim, panels) {
  lines <- c(
    paste0("# ", title),
    "",
    paste0("One-sentence figure claim: ", claim),
    "",
    "| panel_id | panel_title | main_message | source_data_file | plot_type | allowed_interpretation | forbidden_interpretation | priority | notes |",
    "|---|---|---|---|---|---|---|---|---|"
  )
  for (i in seq_len(nrow(panels))) {
    row <- panels[i]
    vals <- vapply(row, as.character, character(1))
    vals <- gsub("\\|", "/", vals)
    lines <- c(lines, paste0("| ", paste(vals, collapse = " | "), " |"))
  }
  writeLines(lines, path)
}

write_blueprint(
  file.path(doc_dir, "figure2_snRNA_context_blueprint.md"),
  "Figure 2 snRNA context blueprint",
  "GSE231569 donor-level snRNA projection supports a moderate Loop/TAL-associated context for primary MAGMA-prioritized modules, with partial random benchmark support and robust driver-removal sensitivity.",
  fig2_panels
)
write_blueprint(
  file.path(doc_dir, "figure3_candidate_evidence_blueprint.md"),
  "Figure 3 candidate evidence blueprint",
  "Stage 3R replaces P1 candidates with a transparent two-axis evidence model and curated exemplars that do not upgrade causal evidence.",
  fig3_panels
)

claim_table <- data.table(
  claim_context = c("main_snRNA_results_claim", "random_benchmark_claim", "driver_removal_claim", "TWAS_proxy_module_claim", "curated_exemplar_claim", "limitations_claim", "figure_caption_claim"),
  recommended_wording = c(
    "MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal.",
    "Expression/detection-matched random benchmarks provided partial support, with rank metrics showing saturation and Loop/TAL-versus-other-compartment delta supporting a conservative signal.",
    "Primary MAGMA module patterns were robust to removal of known transport/ion-handling genes, TAL marker panels, calcium/ion panels, curated exemplars, and top contributors.",
    "TWAS-proxy modules showed weak or inconsistent Loop/TAL patterns and were retained only as supplementary proxy context.",
    "Curated exemplar genes illustrate renal transport and ion-handling biology but do not upgrade evidence strength or validate candidate genes.",
    "The snRNA analysis is donor-level and descriptive, limited by four donors, imbalanced Loop/TAL nuclei, and expression/detection-only random matching.",
    "Figure panels show donor-level descriptive context mapping and robustness boundaries, not causal cell-type mediation."
  ),
  maximum_allowed_strength = c("moderate", "partial", "moderate", "supplementary_only", "supplementary_only", "limitation", "moderate"),
  words_to_avoid = rep("strong enrichment; robust beyond random expectation; causal cell type; causal mediation; papilla-specific regulation; plaque nucleation site; validated candidate genes; SMR-supported; coloc-supported", 7),
  supporting_source = c(paths$stage4b2_claim, paths$stage4b2_random, paths$stage4b2_driver, paths$stage4b2_claim, paths$stage3_exemplar, paths$stage4a_design, paths$stage4b2_claim),
  notes = c(
    "Use as main Results sentence if Stage 4C2 panels preserve visual moderation.",
    "Do not hide partial_support label.",
    "Do not imply driver removal proves independence from known biology.",
    "Keep out of main Loop/TAL proof.",
    "Use role-spectrum framing.",
    "Repeat in Discussion and legend.",
    "Caption must not oversell benchmark."
  )
)
fwrite(claim_table, file.path(table_dir, "snRNA_claim_wording_decision_table.tsv"), sep = "\t")

required_cols <- list(
  "figure2_panelA_compartment_overview_source.tsv" = c("broad_compartment", "n_donor_compartment_rows", "n_donors", "n_modules", "total_nuclei"),
  "scrna_donor_compartment_module_scores.tsv" = c("module_name", "donor_id", "broad_compartment", "mean_module_score", "n_nuclei"),
  "scrna_loop_tal_within_donor_rank.tsv" = c("module_name", "donor_id", "loop_tal_rank", "loop_tal_percentile_rank"),
  "scrna_leave_one_donor_out_module_ranks.tsv" = c("module_name", "excluded_donor", "support_retained", "loop_tal_rank"),
  "scrna_random_set_benchmark_summary.tsv" = c("module_name", "interpretation", "empirical_p_delta", "observed_loop_tal_delta"),
  "scrna_known_driver_removal_sensitivity.tsv" = c("base_module", "removal_type", "support_retained", "support_change", "interpretation"),
  "loop_tal_claim_decision_table.tsv" = c("module_name", "overall_claim_strength", "recommended_manuscript_claim"),
  "candidate_gene_evidence_model_v0.2.tsv" = c("gene", "genetic_priority_level", "twas_proxy_level", "reporting_group"),
  "evidence_model_summary_counts_v0.2.tsv" = c("reporting_group", "n_unique_genes"),
  "curated_exemplar_panel_v0.2.tsv" = c("gene", "biological_role_label", "smr_coloc_status"),
  "snRNA_claim_wording_decision_table.tsv" = c("claim_context", "recommended_wording", "maximum_allowed_strength")
)

panel_manifest_rows <- rbindlist(lapply(list(Figure2 = fig2_panels, Figure3 = fig3_panels), function(panels) {
  fig_name <- names(Filter(function(x) identical(x, panels), list(Figure2 = fig2_panels, Figure3 = fig3_panels)))[1]
  rbindlist(lapply(seq_len(nrow(panels)), function(i) {
    src <- panels$source_data_file[i]
    bname <- basename(src)
    req <- required_cols[[bname]]
    if (is.null(req)) req <- character()
    cols <- if (file.exists(src) && grepl("\\.tsv$", src)) names(fread(src, nrows = 0)) else character()
    missing <- setdiff(req, cols)
    data.table(
      figure = fig_name,
      panel_id = panels$panel_id[i],
      required_source_file = src,
      source_file_exists = exists_yes(src),
      required_columns = collapse0(req),
      columns_present = collapse0(intersect(req, cols)),
      ready_for_plotting = ifelse(file.exists(src) && !length(missing), "yes", "needs_fix"),
      missing_or_problematic_fields = collapse0(missing),
      recommended_fix = ifelse(length(missing), "Add or derive missing source columns before Stage 4C2 plotting.", "Ready for Stage 4C2 draft plotting.")
    )
  }))
}), fill = TRUE)
fwrite(panel_manifest_rows, file.path(table_dir, "figure_panel_source_data_manifest.tsv"), sep = "\t")

risk_table <- data.table(
  risk_id = sprintf("R%02d", 1:9),
  reviewer_criticism = c(
    "Random benchmark only partial support.",
    "Donor number is four.",
    "Loop/TAL nuclei are imbalanced.",
    "Random matching lacks gene length and SNP count.",
    "TWAS proxy modules are weak/supplementary.",
    "Curated exemplars could look cherry-picked.",
    "No SMR/coloc support.",
    "EUR LD used for trans-ancestry GWAS.",
    "57 vs 59 loci reconciliation unresolved."
  ),
  severity = c("high", "high", "high", "medium", "medium", "medium", "high", "high", "medium"),
  current_mitigation = c(
    "Claim downgraded to moderate and benchmark panel labels partial_support.",
    "All snRNA summaries use donor x compartment and no cell-level inference.",
    "Low-count sensitivity and legend caveat included.",
    "Matching limitation explicitly recorded; no gene length/SNP count invented.",
    "TWAS-proxy modules excluded from main Loop/TAL claim.",
    "Exemplar panel separated as interpretive context only.",
    "Figure 3 evidence strip and text state no claim-grade SMR/coloc.",
    "MAGMA described as EUR-LD-reference-based prioritization.",
    "Blocker/reconciliation limitation remains documented from Stage 2."
  ),
  remaining_gap = c(
    "No full beyond-random support for rank metrics.",
    "No independent snRNA replication cohort in Stage 4.",
    "One donor has only four Loop/TAL nuclei.",
    "Potential residual matching bias.",
    "No papilla TWAS/eQTL evidence.",
    "Selection rationale may still need careful captioning.",
    "No causal genetic colocalization layer.",
    "Ancestry-specific LD mismatch remains.",
    "Original 59-locus source unresolved."
  ),
  where_to_address = c("Figure 2E; Results; Limitations", "Figure 2 legend; Limitations", "Figure 2A/D legend; Limitations", "Methods; Limitations", "Figure 3; Results boundary", "Figure 3C-D legend", "Figure 3D; Discussion", "Methods; Limitations", "Supplement/Stage 2 audit"),
  recommended_wording_or_action = c(
    "Use 'partial support beyond matched random expectations', not 'robust enrichment'.",
    "State 'four-donor descriptive analysis'.",
    "Report 29/244/263/4 distribution or cite audit table.",
    "State expression/detection matching only.",
    "Keep R3/R2_R3 supplementary.",
    "Label exemplars as role-spectrum examples, not evidence upgrades.",
    "Use 'no claim-grade SMR/coloc support available'.",
    "Use 'EUR-LD-reference-based prioritization'.",
    "Do not imply 59 original loci were reconciled."
  )
)
fwrite(risk_table, file.path(table_dir, "internal_reviewer_risk_table_stage4C1.tsv"), sep = "\t")

writeLines(c(
  "# Manuscript insertion plan for Stage 4C1",
  "",
  "## A. Methods insertions",
  "",
  "- Candidate-gene evidence modeling: insert after genetic prioritization/MAGMA methods; cite Stage 3R two-axis model.",
  "- Donor-level snRNA module scoring: insert in snRNA analysis methods; state donor x broad_compartment as primary unit.",
  "- Matched random-set benchmark: insert after module scoring; state expression/detection matching and missing gene-length/SNP-count annotations.",
  "- Known-driver removal sensitivity: insert after random benchmark; list single-gene and panel removals.",
  "",
  "## B. Results insertions",
  "",
  "- Evidence model and curated exemplars: introduce before snRNA projection; emphasize P1 removal.",
  "- Donor-level Loop/TAL pattern: report primary MAGMA modules only.",
  "- Random benchmark partial support: explicitly state partial support and rank saturation.",
  "- Known-driver robustness: report robust driver-removal result for primary modules.",
  "- Claim boundary for TWAS-proxy and exemplar modules: keep supplementary/interpretive.",
  "",
  "## C. Discussion insertions",
  "",
  "- Explain why P1 was removed in favor of a two-axis evidence model.",
  "- Explain why Loop/TAL claim is moderate, not strong.",
  "- Explain why random benchmark prevents strong beyond-random claim.",
  "- Explain why driver removal supports robustness but not causality.",
  "- Explain why causal and plaque-nucleation claims remain unsupported.",
  "",
  "## D. Figure legend insertions",
  "",
  "- Figure 2 conservative legend: donor-level descriptive context, partial random support, robust driver removal, no causal cell type.",
  "- Figure 3 conservative legend: two-axis model, Kidney_Cortex TWAS proxy, curated exemplars, no SMR/coloc support, no validated candidate-gene claim."
), file.path(doc_dir, "manuscript_insertion_plan_stage4C1.md"))

writeLines(c(
  "# Draft Figure 2 and Figure 3 legends v0.1",
  "",
  "## Figure 2. Donor-level snRNA context mapping of MAGMA-prioritized modules in GSE231569",
  "",
  "The audited GSE231569 renal papilla single-nucleus object was summarized at the donor x broad-compartment level. Loop/TAL nuclei were present across four donors but were imbalanced across donors, and all module scores are therefore interpreted descriptively. Primary MAGMA-prioritized modules showed donor-level Loop/TAL-associated patterns, supported by within-donor ranking and leave-one-donor-out summaries. Expression/detection-matched random-set benchmarks provided partial support rather than full beyond-random support, while known-driver removal analyses showed robustness to removal of curated exemplars, TAL-marker, calcium/ion, and top-contributor panels. These panels support conservative context mapping only and do not establish causal cell-type mediation, papilla-specific genetic regulation, or a plaque nucleation compartment.",
  "",
  "## Figure 3. Two-axis candidate-gene evidence model and curated biological exemplars",
  "",
  "Candidate genes were organized using a two-axis Stage 3R evidence model separating genetic priority from Kidney_Cortex TWAS proxy support. Reporting groups are mutually exclusive and distinguish MAGMA-prioritized genes, TWAS-proxy-only or TWAS-proxy-supported groups, and lower-priority contextual genes. Curated exemplars illustrate renal transport, calcium/ion handling, and broader epithelial biology but do not upgrade evidence strength. No genes are labeled as SMR-supported, coloc-supported, causal, validated, papilla-specific TWAS-supported, or therapeutic targets."
), file.path(doc_dir, "draft_figure2_figure3_legends_v0.1.md"))

writeLines(c(
  "# Stage 4C1 simulated reviewer check",
  "",
  "1. Would Figure 2 overstate Loop/TAL support?",
  "   No, if the blueprint is followed. Figure 2 must label the claim as moderate and show partial random support explicitly.",
  "",
  "2. Does Figure 2 clearly show random benchmark is partial, not strong?",
  "   Yes. Panel E is dedicated to matched random benchmark and must display partial_support rather than a binary positive mark.",
  "",
  "3. Does Figure 2 separate donor-level evidence from cell-level visualization?",
  "   Yes. Panel A may show atlas context, while Panels B-F use donor x compartment summaries.",
  "",
  "4. Does Figure 3 avoid implying causal gene prioritization?",
  "   Yes. It is a two-axis evidence model and claim-boundary figure, not a causal ranking.",
  "",
  "5. Are TWAS-proxy and curated exemplar modules kept supplementary?",
  "   Yes. They are assigned supplementary/interpretive roles and are not main Loop/TAL proof.",
  "",
  "6. Is the moderate claim wording defensible?",
  "   Yes. It matches donor-level descriptive support, partial matched-random support, and robust driver removal.",
  "",
  "7. What one sentence should be used in Results?",
  "   MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal.",
  "",
  "8. What one sentence should be used in Discussion?",
  "   These snRNA results support a moderate renal papillary Loop/TAL context for genetically prioritized modules, but the partial random benchmark and four-donor design preclude a strong causal or cell-type-mediation claim.",
  "",
  "9. What should not be shown in the main figure?",
  "   Cell-level inferential P values, TWAS-proxy modules as main proof, gene contribution rankings as causal drivers, and curated exemplars as validation.",
  "",
  "10. Is Stage 4C2 draft figure generation ready?",
  "   Yes, after human acceptance of the moderate claim boundary and source-data manifest."
), file.path(doc_dir, "stage4C1_simulated_reviewer_check.md"))

writeLines(c(
  "# Stage 4C1 report: conservative figure planning and claim integration",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Evidence status summary",
  "",
  "- Stage 3R provides the two-axis evidence model and curated exemplar boundaries.",
  "- Stage 4B1 provides donor-level descriptive Loop/TAL support.",
  "- Stage 4B2 provides partial matched-random support and robust known-driver removal.",
  "- The final Loop/TAL claim level is moderate, not strong.",
  "",
  "## Figure 2 recommended panel structure",
  "",
  "- A: audited GSE231569 compartment overview.",
  "- B: donor x compartment module score heatmap.",
  "- C: Loop/TAL within-donor rank summary.",
  "- D: leave-one-donor-out robustness.",
  "- E: matched random-set benchmark with partial_support label.",
  "- F: known-driver removal sensitivity.",
  "- G: claim-decision strip.",
  "",
  "## Figure 3 recommended panel structure",
  "",
  "- A: two-axis evidence model.",
  "- B: reporting group counts.",
  "- C: curated exemplar role spectrum.",
  "- D: exemplar evidence strip.",
  "- E: allowed versus disallowed claim boundary.",
  "",
  "## Source-data readiness",
  "",
  paste0("- Proposed panels checked: ", nrow(panel_manifest_rows)),
  paste0("- Ready for plotting: ", panel_manifest_rows[ready_for_plotting == "yes", .N], " / ", nrow(panel_manifest_rows)),
  "",
  "## Claim wording decision",
  "",
  "Use: MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal.",
  "",
  "## Unresolved risks",
  "",
  "- Random benchmark remains partial, not strong.",
  "- Donor number is four and Loop/TAL nuclei are imbalanced.",
  "- Random matching lacks gene length and SNP-count annotations.",
  "- SMR/coloc and ancestry-matched LD limitations remain.",
  "",
  "## Stage 4C2 readiness",
  "",
  "Stage 4C2 draft figure generation can begin after human acceptance. It must preserve the moderate claim boundary.",
  "",
  "## Manuscript rewrite timing",
  "",
  "Full manuscript rewrite should wait until after Stage 5 and Stage 6, so GSE73680 and spatial/TWAS boundaries can be integrated consistently."
), file.path(doc_dir, "stage4C1_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4C1_completed",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 4A, 4B1, 4B2, and 4C1 completed; conservative Figure 2/3 blueprints, source-data manifest, claim wording table, insertion plan, draft legends, and reviewer risk table generated",
    blocking_issues = "Stage 4C2 draft figure generation not started; full Stage 4 not complete until figures are generated and QC-checked",
    next_stage_ready = "stage4C2_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Completed Stage 4C1 planning\n")
cat("Panels ready:", panel_manifest_rows[ready_for_plotting == "yes", .N], "/", nrow(panel_manifest_rows), "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
