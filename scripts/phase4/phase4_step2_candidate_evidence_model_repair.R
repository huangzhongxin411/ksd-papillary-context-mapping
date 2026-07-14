suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("codex_tasks", recursive = TRUE, showWarnings = FALSE)
dir.create("notes", recursive = TRUE, showWarnings = FALSE)
dir.create("source_data/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing required input: ", path)
    return(NULL)
  }
  read.delim(path, check.names = FALSE, quote = "", comment.char = "")
}

yn <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

first_nonmissing <- function(x, fallback = "") {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x)) as.character(x[[1]]) else fallback
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

magma_path <- "results/tables/phase1_step3_magma_ranked_canonical.tsv"
twas_path <- "results/tables/phase4_step1_twas_proxy_gene_table.tsv"
candidate_path <- "results/tables/candidate_gene_tiers_v1.2.tsv"
scrna_summary_path <- "results/tables/phase2_step5_scrna_integrated_evidence_summary.tsv"
matched_random_path <- "results/tables/phase2_step3_matched_random_benchmark_summary.tsv"
driver_summary_path <- "results/tables/phase2_step4_original_vs_driver_removed_summary.tsv"
driver_random_path <- "results/tables/phase2_step4_driver_removed_matched_random_benchmark_summary.tsv"
spatial_summary_path <- "results/tables/phase3_step4_spatial_integrated_evidence_summary.tsv"

magma <- read_tsv(magma_path)
twas <- read_tsv(twas_path)
cand <- read_tsv(candidate_path)
scrna_summary <- read_tsv(scrna_summary_path, required = FALSE)
matched_random <- read_tsv(matched_random_path, required = FALSE)
driver_summary <- read_tsv(driver_summary_path, required = FALSE)
driver_random <- read_tsv(driver_random_path, required = FALSE)
spatial_summary <- read_tsv(spatial_summary_path, required = FALSE)

for (nm in c("in_top50", "in_top100", "bonferroni_significant", "fdr05_significant", "suggestive_p1e4")) {
  magma[[nm]] <- yn(magma[[nm]])
}

magma_sets <- list(
  MAGMA_top50 = unique(magma$gene_symbol[magma$in_top50]),
  MAGMA_top100 = unique(magma$gene_symbol[magma$in_top100]),
  MAGMA_Bonferroni = unique(magma$gene_symbol[magma$bonferroni_significant]),
  MAGMA_FDR05 = unique(magma$gene_symbol[magma$fdr05_significant]),
  MAGMA_suggestive = unique(magma$gene_symbol[magma$suggestive_p1e4])
)

get_magma_status <- function(g) {
  if (g %in% magma_sets$MAGMA_top50) return("top50")
  if (g %in% magma_sets$MAGMA_top100) return("top100")
  if (g %in% magma_sets$MAGMA_Bonferroni) return("bonferroni")
  if (g %in% magma_sets$MAGMA_FDR05) return("fdr05")
  if (g %in% magma_sets$MAGMA_suggestive) return("suggestive")
  "not_magma_prioritized"
}

gene_universe <- unique(c(as.character(cand$gene), as.character(twas$gene[twas$fdr_supported == TRUE])))
gene_universe <- gene_universe[!is.na(gene_universe) & nzchar(gene_universe)]

magma_idx <- match(gene_universe, magma$gene_symbol)
cand_idx <- match(gene_universe, cand$gene)
twas_idx <- match(gene_universe, twas$gene)

map_scrna_status <- function(x, top_cell = "", percentile = NA_real_) {
  x <- as.character(x)
  top_cell <- as.character(top_cell)
  if (is.na(x) || !nzchar(x)) return("not_assessed")
  if (x == "strong_primary") return("strong_context_support")
  if (x %in% c("review_flagged", "exploratory_low_abundance")) return("partial_context_support")
  if (x == "not_supported") return("no_context_support")
  if (!is.na(percentile) && percentile >= 0.95 && grepl("Loop|TAL", top_cell)) return("moderate_context_support")
  "not_assessed"
}

existing_group <- ifelse(is.na(cand_idx), "", cand$current_tier_v1.2[cand_idx])
scrna_class <- ifelse(is.na(cand_idx), "", cand$scrna_evidence_class[cand_idx])
scrna_top <- ifelse(is.na(cand_idx), "", cand$scrna_top_celltype[cand_idx])
scrna_pct <- ifelse(is.na(cand_idx), NA_real_, safe_num(cand$scrna_benchmark_percentile[cand_idx]))
snrna_status <- mapply(map_scrna_status, scrna_class, scrna_top, scrna_pct, USE.NAMES = FALSE)
snrna_basis <- ifelse(
  snrna_status == "not_assessed",
  "No candidate-level snRNA summary available.",
  paste0("Existing candidate table: scrna_evidence_class=", scrna_class, "; top_celltype=", scrna_top, "; benchmark_percentile=", scrna_pct)
)

twas_model <- ifelse(is.na(twas_idx), "not_tested_or_not_in_model", as.character(twas$model_type[twas_idx]))
twas_model[is.na(twas_model) | !nzchar(twas_model)] <- "not_tested_or_not_in_model"
twas_fdr_supported <- ifelse(is.na(twas_idx), FALSE, twas$fdr_supported[twas_idx] == TRUE)
twas_p <- ifelse(is.na(twas_idx), NA_real_, safe_num(twas$twas_p[twas_idx]))
twas_fdr <- ifelse(is.na(twas_idx), NA_real_, safe_num(twas$twas_fdr[twas_idx]))
twas_proxy_status <- ifelse(
  twas_fdr_supported & twas_model == "multi_snp_proxy", "stronger_proxy_annotation_only",
  ifelse(twas_fdr_supported & twas_model == "one_snp_proxy", "weak_proxy_annotation_only",
    ifelse(twas_fdr_supported & twas_model == "model_size_unknown", "do_not_use_for_claim", "no_positive_twas_evidence")
  )
)

magma_status <- vapply(gene_universe, get_magma_status, character(1))
magma_prioritized <- magma_status != "not_magma_prioritized"
strong_snrna <- snrna_status == "strong_context_support"

assign_group <- function(mg, sn, tw) {
  if (mg && sn && tw == "multi_snp_proxy") return("R2")
  if (mg && sn && tw == "one_snp_proxy") return("R3")
  if (mg && sn) return("R1")
  if (mg) return("R4")
  if (tw %in% c("multi_snp_proxy", "one_snp_proxy", "model_size_unknown")) return("R5")
  "R6"
}
repaired_group <- mapply(assign_group, magma_prioritized, strong_snrna, twas_model, USE.NAMES = FALSE)

group_name <- c(
  R1 = "MAGMA + strong snRNA context",
  R2 = "MAGMA + strong snRNA context + multi-SNP TWAS proxy",
  R3 = "MAGMA + strong snRNA context + one-SNP TWAS proxy",
  R4 = "MAGMA-only or MAGMA + partial snRNA context",
  R5 = "TWAS proxy only",
  R6 = "Contextual or unsupported genes"
)

allowed_interp <- ifelse(
  repaired_group == "R2", "Context-mapped MAGMA-priority candidate with stronger Kidney_Cortex TWAS proxy annotation.",
  ifelse(repaired_group == "R3", "Context-mapped MAGMA-priority candidate with weak one-SNP TWAS proxy annotation.",
    ifelse(repaired_group == "R1", "Context-mapped MAGMA-priority candidate supported by donor-level snRNA context.",
      ifelse(repaired_group == "R4", "Genetic-priority candidate requiring additional context follow-up.",
        ifelse(repaired_group == "R5", "Supplementary Kidney_Cortex TWAS proxy annotation only.", "Contextual background or currently unsupported reporting gene.")
      )
    )
  )
)
not_allowed <- "Causal gene validation; papilla-specific regulatory inference; therapeutic target validity; spatial validation; TWAS-driven priority by one-SNP models."

matrix <- data.frame(
  gene = gene_universe,
  magma_rank = ifelse(is.na(magma_idx), NA_integer_, magma$magma_rank[magma_idx]),
  magma_p = ifelse(is.na(magma_idx), NA_real_, safe_num(magma$magma_p[magma_idx])),
  magma_fdr = ifelse(is.na(magma_idx), NA_real_, safe_num(magma$magma_fdr_bh[magma_idx])),
  magma_status = magma_status,
  in_MAGMA_top50 = gene_universe %in% magma_sets$MAGMA_top50,
  in_MAGMA_top100 = gene_universe %in% magma_sets$MAGMA_top100,
  in_MAGMA_Bonferroni = gene_universe %in% magma_sets$MAGMA_Bonferroni,
  in_MAGMA_FDR05 = gene_universe %in% magma_sets$MAGMA_FDR05,
  in_MAGMA_suggestive = gene_universe %in% magma_sets$MAGMA_suggestive,
  snRNA_support_status = snrna_status,
  snRNA_support_basis = snrna_basis,
  twas_fdr_supported = twas_fdr_supported,
  twas_model_type = twas_model,
  twas_p = twas_p,
  twas_fdr = twas_fdr,
  twas_proxy_status = twas_proxy_status,
  spatial_status = "supplementary_context_only",
  bulk_status_if_available = "pending_not_repaired_in_phase4_step2",
  existing_group = existing_group,
  repaired_reporting_group = repaired_group,
  group_change = ifelse(nzchar(existing_group), paste0(existing_group, "_to_", repaired_group), paste0("new_TWAS_or_context_gene_to_", repaired_group)),
  allowed_interpretation = allowed_interp,
  not_allowed_interpretation = not_allowed,
  notes = "Reporting group is not a causal tier; TWAS and spatial layers cannot upgrade priority beyond MAGMA + snRNA evidence.",
  check.names = FALSE
)
matrix <- matrix[order(matrix$repaired_reporting_group, matrix$magma_rank, matrix$twas_fdr, matrix$gene, na.last = TRUE), ]
write.table(matrix, "results/tables/phase4_step2_candidate_gene_evidence_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

component_rules <- data.frame(
  evidence_component = c(
    "MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR05", "MAGMA_suggestive",
    "snRNA donor-level Loop/TAL-associated context", "matched-random support", "driver-removal sensitivity",
    "TWAS multi-SNP proxy annotation", "TWAS one-SNP proxy annotation",
    "spatial supplementary tissue-context projection", "bulk disease-context status, if present but not yet repaired", "SMR/coloc status"
  ),
  source_file = c(
    magma_path,
    "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt",
    scrna_summary_path,
    matched_random_path,
    paste(driver_summary_path, driver_random_path, sep = ";"),
    twas_path,
    twas_path,
    spatial_summary_path,
    "not_repaired_in_phase4_step2",
    "notes/phase4_step1_smr_coloc_feasibility_note.md"
  ),
  allowed_use = c(
    "Primary genetic-priority layer.",
    "Compact high-priority MAGMA subset for visual emphasis.",
    "Primary/secondary genetic-priority subset for candidate reporting.",
    "FDR-level genetic-priority support.",
    "Suggestive genetic-priority context only.",
    "Donor-level cell-context support; supports Loop/TAL-associated expression context.",
    "Benchmark support that MAGMA module context exceeds matched random expectation.",
    "Sensitivity layer showing contribution of curated drivers without claiming independence.",
    "Stronger Kidney_Cortex TWAS proxy annotation than one-SNP models.",
    "Transparent weak proxy annotation only.",
    "Supplementary anatomical/tissue-context projection only.",
    "Pending layer; do not use to upgrade candidates until repaired.",
    "Resource/provenance status only in this step."
  ),
  not_allowed_use = c(
    "Causal gene evidence.",
    "Causal gene evidence or target validation.",
    "Causal gene evidence or target validation.",
    "Causal gene evidence or target validation.",
    "Primary candidate driver by itself.",
    "Causal cell-type assignment or causal mediation.",
    "Independent validation.",
    "Proof that Loop/TAL biology is dispensable or causal.",
    "Papilla-specific regulatory evidence or causality.",
    "Candidate priority upgrade by itself; papilla-specific regulation.",
    "Candidate priority upgrade; plaque-specific localization; spatial validation.",
    "Candidate priority upgrade in Phase 4-Step 2.",
    "Claim-grade SMR/coloc support unless separately re-audited and approved."
  ),
  strength_label = c(
    "primary_genetic_priority", "primary_genetic_priority_compact", "primary_genetic_priority",
    "primary_genetic_priority", "secondary_genetic_context", "cell_context_support",
    "benchmark_support", "sensitivity_support", "stronger_proxy_annotation", "weak_proxy_annotation",
    "supplementary_context_only", "pending", "not_used_for_claim"
  ),
  notes = c(
    "MAGMA is the genetic-priority layer, not causal validation.",
    "Visual priority can be high because source evidence is MAGMA, not TWAS.",
    "Can populate R1-R4 with snRNA support.",
    "Can populate R1-R4 with snRNA support.",
    "Use cautiously and generally not as a main candidate-driver alone.",
    "Interpret as context mapping across donors.",
    "Use at module level; not gene-level biological validation.",
    "Signals attenuate after driver removal but remain partly retained for some modules.",
    "May annotate R2 but cannot establish regulation.",
    "May annotate R3/R5; should be visually downweighted.",
    "Phase 3 closure retained spatial only as supplementary descriptive context.",
    "Phase 5 should repair this layer before any upgrade use.",
    "Historical artifacts exist but are excluded from Phase 4-Step 2 evidence grouping."
  ),
  check.names = FALSE
)
write.table(component_rules, "results/tables/phase4_step2_evidence_component_rules.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

group_defs <- data.frame(
  reporting_group = paste0("R", 1:6),
  group_name = unname(group_name[paste0("R", 1:6)]),
  required_evidence = c(
    "MAGMA Bonferroni/top100/top50 or FDR05 plus strong donor-level snRNA context.",
    "MAGMA plus strong donor-level snRNA context plus multi-SNP TWAS proxy.",
    "MAGMA plus strong donor-level snRNA context plus one-SNP TWAS proxy.",
    "MAGMA priority without strong snRNA context, or MAGMA plus partial/moderate snRNA context.",
    "TWAS FDR support without MAGMA priority.",
    "Existing table genes without sufficient MAGMA/snRNA support, or contextual background genes."
  ),
  optional_annotation = c(
    "TWAS proxy may be shown but is not required and cannot drive group assignment.",
    "Spatial supplementary context may be listed only as non-upgrading context.",
    "Spatial supplementary context may be listed only as non-upgrading context.",
    "TWAS non-FDR or spatial context may be shown only as non-upgrading annotation.",
    "One- or multi-SNP model type must be displayed; not a priority driver.",
    "None."
  ),
  allowed_claim = c(
    "High-confidence context-mapped genetic-priority candidate.",
    "Context-mapped candidate with stronger proxy TWAS annotation.",
    "Context-mapped candidate with weak TWAS proxy annotation.",
    "Genetic-priority candidate requiring context follow-up.",
    "Proxy TWAS annotation only; supplementary.",
    "Contextual background only."
  ),
  not_allowed_claim = c(
    "Causal gene or validated target.",
    "Papilla-specific regulatory evidence or causal gene.",
    "TWAS-supported regulation or candidate upgrade by one-SNP TWAS alone.",
    "Causal gene, validated target, or context-mapped support if snRNA is weak/absent.",
    "Candidate priority driver, papilla-specific regulation, or causal evidence.",
    "Priority candidate, causal gene, or validated target."
  ),
  visual_priority = c("high", "high_with_proxy_flag", "medium_with_weak_proxy_flag", "medium", "supplementary_low", "background_low"),
  notes = "Reporting groups are not causal tiers.",
  check.names = FALSE
)
write.table(group_defs, "results/tables/phase4_step2_repaired_reporting_group_definitions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

reason <- ifelse(
  matrix$repaired_reporting_group == "R2", "MAGMA priority, strong snRNA context and multi-SNP TWAS proxy annotation.",
  ifelse(matrix$repaired_reporting_group == "R3", "MAGMA priority and strong snRNA context; TWAS is one-SNP proxy only.",
    ifelse(matrix$repaired_reporting_group == "R1", "MAGMA priority and strong snRNA context without FDR TWAS proxy upgrade.",
      ifelse(matrix$repaired_reporting_group == "R4", "MAGMA priority is present but strong snRNA context is absent or partial.",
        ifelse(matrix$repaired_reporting_group == "R5", "FDR TWAS proxy exists without MAGMA priority; supplementary only.", "Insufficient repaired MAGMA/snRNA/TWAS support for priority reporting.")
      )
    )
  )
)
reporting <- data.frame(
  gene = matrix$gene,
  existing_group = matrix$existing_group,
  repaired_reporting_group = matrix$repaired_reporting_group,
  reason_for_group_assignment = reason,
  magma_summary = paste0(matrix$magma_status, ifelse(is.na(matrix$magma_rank), "", paste0("; rank=", matrix$magma_rank))),
  snRNA_summary = paste0(matrix$snRNA_support_status, "; ", matrix$snRNA_support_basis),
  twas_summary = paste0(matrix$twas_proxy_status, "; model=", matrix$twas_model_type, ifelse(is.na(matrix$twas_fdr), "", paste0("; FDR=", signif(matrix$twas_fdr, 3)))),
  spatial_summary = "Supplementary context only; not claim-grade and not used for group assignment.",
  bulk_summary_if_available = matrix$bulk_status_if_available,
  downgrade_reason = "",
  upgrade_reason = "",
  human_review_required = "no",
  notes = "Do not call this a causal tier.",
  check.names = FALSE
)
reporting$downgrade_reason[reporting$repaired_reporting_group %in% c("R3", "R5")] <- "One-SNP or TWAS-only signal is retained only as proxy annotation and cannot drive priority."
reporting$downgrade_reason[reporting$repaired_reporting_group == "R6"] <- "Insufficient repaired priority evidence under conservative rules."
reporting$upgrade_reason[reporting$repaired_reporting_group == "R2"] <- "Multi-SNP TWAS retained as stronger proxy annotation, not regulatory evidence."
reporting$human_review_required[reporting$repaired_reporting_group %in% c("R2", "R3", "R5")] <- "yes"
write.table(reporting, "results/tables/phase4_step2_candidate_reporting_groups_repaired.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

old_candidate_group <- ifelse(is.na(twas_idx), "", existing_group)
twas_audit <- matrix[matrix$twas_fdr_supported | matrix$twas_model_type %in% c("multi_snp_proxy", "one_snp_proxy", "model_size_unknown"), ]
twas_audit$downgrade_or_no_change <- ifelse(
  twas_audit$twas_model_type == "one_snp_proxy" & twas_audit$repaired_reporting_group %in% c("R3", "R5"), "downgraded_or_visually_downweighted",
  ifelse(twas_audit$twas_model_type == "multi_snp_proxy", "retained_as_proxy_annotation_only", "no_change")
)
twas_downgrade <- data.frame(
  gene = twas_audit$gene,
  twas_fdr_supported = twas_audit$twas_fdr_supported,
  twas_model_type = twas_audit$twas_model_type,
  old_candidate_group = twas_audit$existing_group,
  new_reporting_group = twas_audit$repaired_reporting_group,
  downgrade_or_no_change = twas_audit$downgrade_or_no_change,
  reason = ifelse(twas_audit$twas_model_type == "one_snp_proxy",
    "One-SNP Kidney_Cortex TWAS is high-risk proxy annotation and cannot upgrade priority.",
    ifelse(twas_audit$twas_model_type == "multi_snp_proxy",
      "Multi-SNP Kidney_Cortex TWAS is retained only as stronger proxy annotation.",
      "Model-size unknown or unsupported; do not use for claim."
    )
  ),
  safe_wording = "Kidney_Cortex S-PrediXcan was retained as proxy annotation and not interpreted as papilla-specific regulatory evidence.",
  unsafe_wording = "TWAS validates causal genes; TWAS confirms regulation; one-SNP TWAS upgrades candidate priority.",
  check.names = FALSE
)
write.table(twas_downgrade, "results/tables/phase4_step2_twas_downgrade_audit.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

candidate_matrix_source <- matrix[, c(
  "gene", "magma_status", "snRNA_support_status", "twas_model_type", "twas_proxy_status",
  "spatial_status", "bulk_status_if_available", "repaired_reporting_group", "allowed_interpretation"
)]
candidate_matrix_source <- candidate_matrix_source[order(candidate_matrix_source$repaired_reporting_group, candidate_matrix_source$gene), ]
write.table(candidate_matrix_source, "source_data/figures/phase4_step2_Figure3_source_data_candidate_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

twas_burden <- data.frame(
  category = factor(
    c("FDR-supported TWAS genes", "one-SNP proxy models", "multi-SNP proxy models"),
    levels = c("FDR-supported TWAS genes", "one-SNP proxy models", "multi-SNP proxy models")
  ),
  n_genes = c(sum(twas$fdr_supported == TRUE), sum(twas$fdr_supported == TRUE & twas$model_type == "one_snp_proxy"), sum(twas$fdr_supported == TRUE & twas$model_type == "multi_snp_proxy")),
  interpretation = c("Proxy annotation layer", "Weak/high-risk proxy annotation", "Stronger proxy annotation but not regulatory evidence"),
  check.names = FALSE
)
write.table(twas_burden, "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

group_counts <- as.data.frame(table(matrix$repaired_reporting_group), stringsAsFactors = FALSE)
names(group_counts) <- c("reporting_group", "n_genes")
group_counts$group_name <- group_name[group_counts$reporting_group]
all_groups <- data.frame(reporting_group = paste0("R", 1:6), stringsAsFactors = FALSE)
group_counts <- merge(all_groups, group_counts, by = "reporting_group", all.x = TRUE, sort = FALSE)
group_counts$n_genes[is.na(group_counts$n_genes)] <- 0
group_counts$group_name <- group_name[group_counts$reporting_group]
write.table(group_counts, "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

panel_manifest <- data.frame(
  panel = c("A", "B", "C", "D", "E"),
  panel_title = c(
    "Conservative evidence model",
    "Candidate evidence matrix",
    "TWAS model-quality burden",
    "Reporting group counts",
    "Claim boundary"
  ),
  source_data = c(
    "results/tables/phase4_step2_evidence_component_rules.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_candidate_matrix.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv",
    "notes/phase4_step2_figure3_redesign_plan.md"
  ),
  main_message = c(
    "MAGMA is primary; snRNA provides context; TWAS is proxy annotation; spatial is supplementary.",
    "Rows are reporting genes/groups, not causal tiers.",
    "Most FDR-supported TWAS genes use one-SNP models.",
    "Repaired R1-R6 groups prevent TWAS-only or spatial-driven priority.",
    "Supported and unsupported claims are explicit."
  ),
  claim_boundary = c(
    "No causal arrows or regulatory validation.",
    "One-SNP TWAS and spatial context cannot upgrade priority.",
    "TWAS is Kidney_Cortex proxy only.",
    "Reporting groups are not causal tiers.",
    "No causal genes, papilla-specific regulation or therapeutic targets."
  ),
  manual_polishing_needed = c("yes", "yes", "no", "no", "yes"),
  notes = "Draft figure only; final publication-ready Figure 3 requires human design review.",
  check.names = FALSE
)
write.table(panel_manifest, "results/tables/phase4_step2_Figure3_panel_manifest.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "# Phase 4-Step 2 Figure 3 Redesign Plan",
  "",
  "Purpose: show a conservative candidate-reporting framework, not a causal prioritization hierarchy.",
  "",
  "Panel A: Evidence model schematic. MAGMA genetic priority flows into snRNA donor-level context mapping, with TWAS displayed as a Kidney_Cortex proxy annotation layer. Spatial appears as supplementary context only. Do not draw causal arrows.",
  "",
  "Panel B: Candidate evidence matrix. Rows should be selected reporting genes or grouped rows; columns should show MAGMA, snRNA, TWAS model type, spatial context and bulk disease-context reviewed/non-upgrading status. Use distinct marks for primary evidence, proxy annotation, supplementary context and not claim-grade.",
  "",
  "Panel C: TWAS model-quality burden. Show 51 FDR-supported genes, including 42 one-SNP and 9 multi-SNP models. The visual message is downgrade, not support inflation.",
  "",
  "Panel D: Reporting group summary. Display counts for repaired R1-R6 groups and label clearly: `Reporting groups, not causal tiers.`",
  "",
  "Panel E: Claim boundary box. Supported: context-mapped genetic-priority candidates. Not claimed: causal genes, papilla-specific regulation, therapeutic targets, spatial validation or claim-grade SMR/coloc.",
  "",
  "If density becomes excessive, move Panel C to a supplementary TWAS proxy figure and retain only a compact TWAS model-quality flag in Panel B."
), "notes/phase4_step2_figure3_redesign_plan.md")

plot_matrix_all <- matrix[order(matrix$repaired_reporting_group, matrix$magma_rank, matrix$twas_fdr, na.last = TRUE), ]
group_caps <- c(R1 = 10, R2 = 5, R3 = 5, R4 = 10, R5 = 8, R6 = 6)
plot_matrix <- do.call(rbind, lapply(names(group_caps), function(g) {
  rows <- plot_matrix_all[plot_matrix_all$repaired_reporting_group == g, ]
  if (nrow(rows) == 0) return(rows)
  head(rows, group_caps[[g]])
}))
layers <- data.frame(
  gene = rep(plot_matrix$gene, each = 5),
  layer = rep(c("MAGMA", "snRNA context", "TWAS proxy", "Spatial", "Bulk"), times = nrow(plot_matrix)),
  state = NA_character_,
  group = rep(plot_matrix$repaired_reporting_group, each = 5),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(plot_matrix))) {
  idx <- layers$gene == plot_matrix$gene[i]
  layers$state[idx & layers$layer == "MAGMA"] <- ifelse(plot_matrix$magma_status[i] %in% c("top50", "top100", "bonferroni", "fdr05"), "primary", ifelse(plot_matrix$magma_status[i] == "suggestive", "secondary", "none"))
  layers$state[idx & layers$layer == "snRNA context"] <- ifelse(plot_matrix$snRNA_support_status[i] == "strong_context_support", "context", ifelse(plot_matrix$snRNA_support_status[i] %in% c("moderate_context_support", "partial_context_support"), "partial", "none"))
  layers$state[idx & layers$layer == "TWAS proxy"] <- ifelse(plot_matrix$twas_model_type[i] == "multi_snp_proxy", "proxy_stronger", ifelse(plot_matrix$twas_model_type[i] == "one_snp_proxy", "proxy_weak", "none"))
  layers$state[idx & layers$layer == "Spatial"] <- "supplementary"
  layers$state[idx & layers$layer == "Bulk"] <- "pending"
}
layers$gene <- factor(layers$gene, levels = rev(unique(plot_matrix$gene)))
layers$layer <- factor(layers$layer, levels = c("MAGMA", "snRNA context", "TWAS proxy", "Spatial", "Bulk"))

pal <- c(
  primary = "#245A64", secondary = "#7F9DA6", context = "#0F4C5C", partial = "#B99B5A",
  proxy_stronger = "#9B5C4D", proxy_weak = "#D7B8AE", supplementary = "#BFC9CC", pending = "#E6E9EA", none = "#F4F5F5"
)

pdf_path <- "results/figures/phase4_step2_Figure3_candidate_reporting_model_draft.pdf"
png_path <- "results/figures/phase4_step2_Figure3_candidate_reporting_model_draft.png"

draw_figure <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(pdf_path, width = 13, height = 9, family = "Helvetica")
  if (device == "png") {
    png_type <- if (capabilities("aqua")) "quartz" else "Xlib"
    png(png_path, width = 7800, height = 5400, res = 600, type = png_type)
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(0.55, 1.15, 5.2, 1.35), "null"), widths = unit(c(1, 1), "null"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text("Figure 3 draft: repaired candidate reporting model", x = 0.02, y = 0.72, just = c("left", "center"), gp = gpar(fontsize = 17, fontface = "bold", col = "#333333"))
  grid.text("Kidney_Cortex TWAS is proxy annotation only; reporting groups are not causal tiers", x = 0.02, y = 0.26, just = c("left", "center"), gp = gpar(fontsize = 10.5, col = "#333333"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text("A", x = 0.02, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text("Evidence model", x = 0.09, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
  xs <- c(0.18, 0.42, 0.66, 0.86)
  labs <- c("MAGMA\nprimary genetic priority", "snRNA\ncell-context support", "TWAS\nproxy annotation", "Spatial\nsupplement only")
  fills <- c("#245A64", "#0F4C5C", "#D7B8AE", "#E6E9EA")
  for (j in seq_along(xs)) {
    grid.roundrect(x = xs[j], y = 0.45, width = 0.19, height = 0.36, r = unit(0.04, "snpc"), gp = gpar(fill = fills[j], col = "#333333", lwd = 0.6))
    grid.text(labs[j], x = xs[j], y = 0.45, gp = gpar(fontsize = 8.5, col = ifelse(j < 3, "white", "#333333"), fontface = ifelse(j == 1, "bold", "plain")))
  }
  grid.lines(x = unit(c(0.28, 0.32), "npc"), y = unit(c(0.45, 0.45), "npc"), gp = gpar(col = "#333333", lwd = 1.2))
  grid.lines(x = unit(c(0.52, 0.56), "npc"), y = unit(c(0.45, 0.45), "npc"), gp = gpar(col = "#333333", lwd = 1.2))
  grid.text("annotation, not causality", x = 0.66, y = 0.15, gp = gpar(fontsize = 8.5, col = "#9B5C4D"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  burden_plot <- ggplot(twas_burden, aes(x = category, y = n_genes, fill = category)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = n_genes), vjust = -0.25, size = 3.2) +
    scale_fill_manual(values = c("FDR-supported TWAS genes" = "#7F9DA6", "one-SNP proxy models" = "#D7B8AE", "multi-SNP proxy models" = "#9B5C4D")) +
    scale_y_continuous(limits = c(0, max(twas_burden$n_genes) * 1.18), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "C  TWAS model-quality burden", x = NULL, y = "Genes") +
    theme_minimal(base_family = "Helvetica", base_size = 9) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 15, hjust = 1), plot.title = element_text(face = "bold", size = 12))
  print(burden_plot, newpage = FALSE)
  popViewport()

  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  mat_plot <- ggplot(layers, aes(x = layer, y = gene, fill = state)) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = pal, breaks = names(pal), name = "Evidence role") +
    labs(title = "B  Candidate evidence matrix", subtitle = "Selected reporting genes; TWAS and spatial columns are non-upgrading annotations", x = NULL, y = NULL) +
    theme_minimal(base_family = "Helvetica", base_size = 8.5) +
    theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          axis.text.x = element_text(face = "bold"), legend.position = "right",
          plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
  print(mat_plot, newpage = FALSE)
  popViewport()

  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  counts_plot <- ggplot(group_counts, aes(x = reporting_group, y = n_genes, fill = reporting_group)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = n_genes), vjust = -0.25, size = 3.1) +
    scale_fill_manual(values = c(R1 = "#245A64", R2 = "#9B5C4D", R3 = "#D7B8AE", R4 = "#7F9DA6", R5 = "#B99B5A", R6 = "#E6E9EA")) +
    scale_y_continuous(limits = c(0, max(group_counts$n_genes) * 1.18), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "D  Reporting groups, not causal tiers", x = NULL, y = "Genes") +
    theme_minimal(base_family = "Helvetica", base_size = 9) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold", size = 12))
  print(counts_plot, newpage = FALSE)
  popViewport()

  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
  grid.text("E", x = 0.02, y = 0.92, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text("Claim boundary", x = 0.09, y = 0.92, just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
  grid.roundrect(x = 0.5, y = 0.48, width = 0.9, height = 0.62, r = unit(0.03, "snpc"), gp = gpar(fill = "#F7F7F7", col = "#333333", lwd = 0.7))
  grid.text("Supported: context-mapped genetic-priority candidates\n\nNot claimed: causal genes, papilla-specific regulation,\ntherapeutic targets, spatial validation, or claim-grade SMR/coloc", x = 0.09, y = 0.66, just = c("left", "top"), gp = gpar(fontsize = 9.5, col = "#333333"))
  popViewport()

  dev.off()
}

draw_figure("pdf")
draw_figure("png")

legend_text <- c(
  "# Figure 3 Draft Legend: Repaired Candidate Reporting Model",
  "",
  "Draft Figure 3 summarizes the repaired candidate evidence model after downgrading the TWAS layer. Panel A shows the permitted evidence flow: MAGMA defines genetic-priority genes, snRNA provides donor-level Loop/TAL-associated context, TWAS is retained as GTEx Kidney_Cortex proxy annotation, and spatial transcriptomics is supplementary context only. Panel B displays a selected candidate evidence matrix using repaired reporting groups rather than causal tiers. Panel C summarizes TWAS model quality, showing 51 FDR-supported Kidney_Cortex S-PrediXcan genes, of which 42 used one-SNP prediction models and 9 used multi-SNP models. Panel D reports the number of genes assigned to repaired R1-R6 reporting groups. Panel E states the claim boundary. The figure supports context-mapped genetic-priority candidates and does not claim causal genes, papilla-specific regulatory effects, therapeutic target validity, spatial validation, or claim-grade SMR/coloc support."
)
writeLines(legend_text, "notes/phase4_step2_Figure3_draft_legend.md")

figure_qc <- data.frame(
  figure_id = "phase4_step2_Figure3_candidate_reporting_model_draft",
  version = "draft",
  pdf_exists = file.exists(pdf_path),
  png_exists = file.exists(png_path),
  png_dpi_or_intended_dpi = "600",
  minimum_configured_font_size = "8.5 pt",
  panel_label_presence = "A-E present",
  legend_placement_check = "Right-side evidence legend in matrix; no legend inside dense data except tile legend.",
  palette_consistency_check = "Uses project profile colors: deep teal, loop teal, bluegrey, sand gold, terracotta, pale grey.",
  claim_boundary_check = "Boundary text included; TWAS shown as proxy annotation only.",
  resource_limited_claim_check = "No SMR/coloc or spatial validation claim.",
  source_table_existence = all(file.exists(c(
    "source_data/figures/phase4_step2_Figure3_source_data_candidate_matrix.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv"
  ))),
  legend_file_existence = file.exists("notes/phase4_step2_Figure3_draft_legend.md"),
  visual_status = "agent_visual_review_passed_for_draft_requires_human_polish",
  action_required = "Manual polish before any final publication-ready Figure 3; keep draft label until Step 3 approval.",
  check.names = FALSE
)
write.table(figure_qc, "results/tables/phase4_step2_Figure3_visual_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "# Phase 4-Step 2 Results Wording",
  "",
  paste0("The TWAS layer was retained as a transparent proxy annotation rather than as regulatory or causal evidence. In the canonical GTEx v8 Kidney_Cortex S-PrediXcan table, 51 genes reached FDR support; 42 of these used one-SNP prediction models and 9 used multi-SNP models. TWAS-MAGMA overlap was therefore described as proxy convergence with the genetic-priority layer, not as papilla-specific regulation. Candidate reporting groups were repaired so that MAGMA remained the primary genetic-priority layer, snRNA provided donor-level Loop/TAL-associated context, and TWAS model class was displayed only as annotation. One-SNP TWAS no longer drives candidate priority, and R1-R6 labels are reporting groups rather than causal tiers.")
), "notes/phase4_step2_results_wording.md")

writeLines(c(
  "# Phase 4-Step 2 Methods Wording",
  "",
  "Candidate evidence repair used the Phase 1 canonical MAGMA ranked table, the Phase 4-Step 1 TWAS proxy table and the existing candidate reporting table. MAGMA membership was classified hierarchically as top50, top100, Bonferroni, FDR05, suggestive or not MAGMA-prioritized. Candidate-level snRNA context was taken from the existing donor-level evidence fields and harmonized as strong, moderate, partial, absent or not assessed. TWAS model quality was classified from the audited `n_snps_used`-derived model type as multi-SNP proxy, one-SNP proxy, model-size unknown, not FDR-supported or not tested. MAGMA-TWAS overlap was summarized descriptively. Repaired reporting groups prioritized MAGMA plus snRNA context; TWAS multi-SNP signals were retained only as stronger proxy annotation and one-SNP signals only as weak proxy annotation. Spatial evidence was retained as supplementary tissue-context information only. No TWAS rerun, SMR, colocalization or bulk analysis was performed in this step, and no causal inference was made from TWAS."
), "notes/phase4_step2_methods_wording.md")

writeLines(c(
  "# Phase 4-Step 2 Limitations Wording",
  "",
  "The TWAS annotation used GTEx Kidney_Cortex rather than renal papilla, creating a tissue-context mismatch. Most FDR-supported TWAS genes used one-SNP prediction models, which limits interpretability and prevents TWAS from driving candidate priority. Papilla-specific eQTL resources were not available, and no claim-grade SMR or colocalization analysis was used in this step. TWAS results are therefore proxy annotations only. The repaired R1-R6 labels are reporting groups for organizing evidence, not causal tiers, therapeutic target rankings or proof of papilla-specific regulation."
), "notes/phase4_step2_limitations_wording.md")

reviewer_issues <- data.frame(
  reviewer_concern = c(
    "TWAS overinterpretation", "one-SNP TWAS predominance", "Kidney_Cortex tissue mismatch",
    "candidate tier overclaim", "Figure 3 visual overweighting of TWAS", "lack of SMR/coloc",
    "spatial not claim-grade", "reporting groups not causal tiers"
  ),
  action_taken = c(
    "TWAS reclassified as proxy annotation only.",
    "One-SNP and multi-SNP TWAS classes explicitly separated.",
    "Kidney_Cortex proxy tissue named in rules, wording and figure plan.",
    "Old tier language replaced by repaired R1-R6 reporting groups.",
    "Figure 3 draft makes MAGMA/snRNA primary and TWAS visually downweighted as proxy.",
    "SMR/coloc excluded from evidence grouping pending separate provenance audit.",
    "Spatial fixed as supplementary context only.",
    "Definitions and figure labels state reporting groups are not causal tiers."
  ),
  result = c(
    "No TWAS causal/regulatory claim remains in repaired model.",
    "42 one-SNP and 9 multi-SNP FDR TWAS genes reported.",
    "Tissue mismatch preserved as a limitation.",
    "Candidate model repaired with conservative allowed/not-allowed claims.",
    "Draft Figure 3 includes a TWAS burden panel and claim boundary box.",
    "No claim-grade SMR/coloc evidence used.",
    "No spatial upgrade of candidate priority.",
    "R1-R6 used only as reporting labels."
  ),
  remaining_limitation = c(
    "TWAS remains useful only as proxy annotation.",
    "One-SNP burden limits interpretability.",
    "No papilla-specific TWAS/eQTL resource.",
    "Human review needed before manuscript/figure finalization.",
    "Manual design polish needed.",
    "Optional future provenance audit may be needed.",
    "Spatial remains supplementary only.",
    "Reader-facing wording must preserve this distinction."
  ),
  manuscript_or_supplement_location = c(
    "Results; Methods; Supplementary TWAS table.",
    "Figure 3C or Supplementary TWAS figure.",
    "Limitations; Methods.",
    "Figure 3B-D; candidate evidence table.",
    "Figure 3 redesign plan and draft.",
    "Limitations; future work.",
    "Supplementary spatial figure; limitations.",
    "Figure 3D label; Results wording."
  ),
  check.names = FALSE
)
write.table(reviewer_issues, "codex_tasks/phase4_step2_reviewer_response_candidate_TWAS_issue_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

one_snp_downgraded <- sum(twas_downgrade$twas_model_type == "one_snp_proxy" & twas_downgrade$downgrade_or_no_change == "downgraded_or_visually_downweighted")
multi_snp_retained <- sum(twas_downgrade$twas_model_type == "multi_snp_proxy")
count_lines <- paste0("- ", group_counts$reporting_group, " (", group_counts$group_name, "): ", group_counts$n_genes, collapse = "\n")
writeLines(c(
  "# Phase 4-Step 2 Report",
  "",
  "## Candidate Model Status",
  "",
  "Candidate evidence model repaired with conservative R1-R6 reporting groups. Reporting groups are not causal tiers.",
  "",
  "## Repaired Reporting Group Counts",
  "",
  count_lines,
  "",
  "## TWAS Downgrade",
  "",
  paste0("- One-SNP TWAS genes downgraded or visually downweighted: ", one_snp_downgraded, "."),
  paste0("- Multi-SNP TWAS genes retained as stronger proxy annotation only: ", multi_snp_retained, "."),
  "",
  "## Figure 3 Redesign Status",
  "",
  "- Figure 3 redesign plan created.",
  "- Draft Figure 3 generated as a non-final working figure.",
  "- Final publication-ready Figure 3 was not generated.",
  "",
  "## Safe Candidate Claim",
  "",
  "Repaired reporting groups identify context-mapped genetic-priority candidates and transparent proxy annotations.",
  "",
  "## Unsafe Candidate Claims",
  "",
  "- Causal genes.",
  "- Papilla-specific regulatory effects.",
  "- Therapeutic target validity.",
  "- TWAS-driven priority from one-SNP models.",
  "- Spatial validation or claim-grade SMR/coloc support.",
  "",
  "## Recommended Next Step",
  "",
  "A. proceed to Phase 4-Step 3: Figure 3 polishing and Supplementary TWAS figure assembly.",
  "",
  "Human review should confirm the repaired R1-R6 group definitions before manuscript or DOCX edits."
), "notes/phase4_step2_report.md")

checklist <- data.frame(
  task_id = sprintf("P4S2-%02d", 1:15),
  task_name = c(
    "Candidate evidence repair script",
    "Evidence component rules",
    "Gene-level evidence matrix",
    "Reporting group definitions",
    "Candidate groups repaired",
    "TWAS downgrade audit",
    "Figure 3 redesign plan",
    "Figure 3 source data: candidate matrix",
    "Figure 3 source data: TWAS burden",
    "Figure 3 source data: reporting group counts",
    "Draft Figure 3 or panel manifest",
    "Results wording",
    "Methods and limitations wording",
    "Reviewer issue table",
    "Stop-rule compliance"
  ),
  completed = "yes",
  output_file = c(
    "scripts/10_integrated_figures/phase4_step2_candidate_evidence_model_repair.R",
    "results/tables/phase4_step2_evidence_component_rules.tsv",
    "results/tables/phase4_step2_candidate_gene_evidence_matrix.tsv",
    "results/tables/phase4_step2_repaired_reporting_group_definitions.tsv",
    "results/tables/phase4_step2_candidate_reporting_groups_repaired.tsv",
    "results/tables/phase4_step2_twas_downgrade_audit.tsv",
    "notes/phase4_step2_figure3_redesign_plan.md",
    "source_data/figures/phase4_step2_Figure3_source_data_candidate_matrix.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv",
    "results/figures/phase4_step2_Figure3_candidate_reporting_model_draft.pdf;results/figures/phase4_step2_Figure3_candidate_reporting_model_draft.png;results/tables/phase4_step2_Figure3_panel_manifest.tsv",
    "notes/phase4_step2_results_wording.md",
    "notes/phase4_step2_methods_wording.md;notes/phase4_step2_limitations_wording.md",
    "codex_tasks/phase4_step2_reviewer_response_candidate_TWAS_issue_table.tsv",
    "codex_tasks/phase4_step2_completion_checklist.tsv"
  ),
  blocking_issue = c(rep("", 10), "Draft only; manual polish required before final Figure 3.", rep("", 3), ""),
  manual_review_needed = c("no", "yes", "yes", "yes", "yes", "yes", "yes", "no", "no", "no", "yes", "yes", "yes", "yes", "yes"),
  notes = c(
    "Script reads existing inputs only and does not run TWAS/SMR/coloc/bulk.",
    "Rules encode MAGMA primary, snRNA context, TWAS proxy, spatial supplementary status.",
    "Union of existing candidates and FDR TWAS genes.",
    "R1-R6 are reporting groups, not causal tiers.",
    "One-SNP TWAS cannot upgrade priority.",
    "All FDR TWAS genes audited for downgrade/proxy status.",
    "Plan follows conservative figure structure.",
    "Ready for later plotting/polish.",
    "51/42/9 TWAS burden captured.",
    "Group counts include empty groups if any.",
    "Non-final draft generated for review.",
    "Manuscript wording prepared but DOCX not edited.",
    "Methods/limitations wording prepared.",
    "Reviewer-facing response issues prepared.",
    "No manuscript DOCX edit; no SMR/coloc; no bulk analyses; no final Figure 3."
  ),
  check.names = FALSE
)
write.table(checklist, "codex_tasks/phase4_step2_completion_checklist.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Phase 4-Step 2 candidate evidence repair completed.\n")
