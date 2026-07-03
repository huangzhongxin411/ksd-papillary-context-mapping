suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

dir.create("results/phase2d_p1_gene_evidence_v0.1", recursive = TRUE, showWarnings = FALSE)
dir.create("results/phase2d_p1_gene_evidence_v0.1/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/phase2d_p1_gene_evidence_v0.1/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/phase2d_p1_gene_evidence_v0.1/docs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

copy_if_exists <- function(from, to_dir) {
  if (!file.exists(from)) return(FALSE)
  file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
}

freeze_tables <- c(
  "results/tables/p1_tal_gene_evidence.tsv",
  "results/tables/p1_tal_gene_interpretation_summary.tsv",
  "results/tables/p1_figure_qc_summary.tsv",
  "results/tables/p1_tal_gene_celltype_summary.tsv",
  "results/tables/p1_tal_gene_by_donor.tsv",
  "results/tables/p1_tal_gene_specificity.tsv",
  "results/tables/p1_gene_vs_tal_program_correlation.tsv"
)
freeze_figures <- c(
  "results/figures/p1_tal_gene_featureplots.pdf",
  "results/figures/p1_tal_gene_dotplot.pdf",
  "results/figures/p1_tal_gene_heatmap_by_celltype.pdf",
  "results/figures/p1_tal_gene_by_donor_boxplot.pdf",
  "results/figures/p1_tal_gene_specificity_barplot.pdf",
  "results/figures/p1_gene_vs_tal_program_correlation.pdf"
)
freeze_docs <- c(
  "docs/phase2d_p1_gene_evidence_notes.md",
  "docs/results_draft_v0.1.md"
)

copied <- rbindlist(list(
  data.table(source_file = freeze_tables, freeze_subdir = "tables",
             copied = vapply(freeze_tables, copy_if_exists, logical(1), to_dir = "results/phase2d_p1_gene_evidence_v0.1/tables")),
  data.table(source_file = freeze_figures, freeze_subdir = "figures",
             copied = vapply(freeze_figures, copy_if_exists, logical(1), to_dir = "results/phase2d_p1_gene_evidence_v0.1/figures")),
  data.table(source_file = freeze_docs, freeze_subdir = "docs",
             copied = vapply(freeze_docs, copy_if_exists, logical(1), to_dir = "results/phase2d_p1_gene_evidence_v0.1/docs"))
))
copied[, frozen_file := file.path("results/phase2d_p1_gene_evidence_v0.1", freeze_subdir, basename(source_file))]
copied[, notes := ifelse(copied, "included in Phase 2D freeze v0.1", "source file missing at freeze time")]
fwrite(copied, "results/phase2d_p1_gene_evidence_v0.1/MANIFEST.tsv", sep = "\t")

writeLines(c(
  "# Phase 2D P1 Gene Evidence Freeze v0.1",
  "",
  "## Freeze conclusion",
  "",
  "P1 core TAL candidate gene evidence pack supports interpretable TAL-associated or epithelial transport-associated single-cell expression contexts for UMOD, CLDN10, CLDN14, CASR, HIBADH and PKD2.",
  "",
  "Chinese summary: P1 六个核心候选基因已获得可解释的 TAL 相关或上皮转运相关单细胞表达背景证据。该结果支持 MAGMA + scRNA 层面的 cellular context，但不构成 TWAS 收敛、coloc 共定位、空间验证或因果证据。",
  "",
  "## Claim boundary",
  "",
  "- Supported: MAGMA + scRNA-supported TAL-associated cellular context.",
  "- Not supported yet: TWAS convergence, SMR/coloc validation, spatial validation, or causal mediation.",
  "",
  "## Frozen artifacts",
  "",
  "See `MANIFEST.tsv` in this freeze directory."
), "results/phase2d_p1_gene_evidence_v0.1/phase2d_freeze_decision.md")

p1 <- fread("results/tables/p1_tal_gene_evidence.tsv")
interp <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
celltype <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
donor <- fread("results/tables/p1_tal_gene_by_donor.tsv")
specificity <- fread("results/tables/p1_tal_gene_specificity.tsv")

gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
cell_order <- celltype[gene == "UMOD"][order(expression_rank_across_celltypes)]$cell_type
celltype[, gene := factor(gene, levels = rev(gene_order))]
celltype[, cell_type := factor(cell_type, levels = rev(cell_order))]

panel_a <- ggplot(celltype, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = "grey25", stroke = 0.25) +
  scale_size_continuous(range = c(1.4, 7), labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_gradient(low = "#f6f7f1", high = "#1f6f78") +
  labs(x = NULL, y = NULL, size = "Detected", fill = "Avg expr", title = "A. P1 gene expression across audited cell types") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", size = 10),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2))

specificity[, gene := factor(gene, levels = gene_order)]
panel_b <- ggplot(specificity, aes(x = gene, y = specificity_ratio_avg, fill = specificity_class)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  scale_fill_manual(values = c(
    strong_TAL_preferential = "#2c7a7b",
    moderate_TAL_preferential = "#77a65d",
    broad_expression_with_TAL_component = "#b78a3b"
  )) +
  labs(x = NULL, y = "TAL specificity ratio", fill = "Specificity", title = "B. TAL specificity separates gene roles") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold", size = 10),
        legend.position = "bottom")

donor_genes <- c("UMOD", "CLDN10", "CASR", "PKD2")
donor_tal <- donor[is_TAL == TRUE & gene %in% donor_genes]
donor_tal[, gene := factor(gene, levels = donor_genes)]
panel_c <- ggplot(donor_tal, aes(x = gene, y = avg_expression)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "#dce6df", color = "grey25", linewidth = 0.25) +
  geom_point(aes(color = disease_status), position = position_jitter(width = 0.08, height = 0), size = 2.3, alpha = 0.9) +
  scale_color_manual(values = c(healthy_control = "#2f5597", disease = "#b44a3c", stone = "#b44a3c"), na.value = "grey45") +
  labs(x = NULL, y = "Mean expression in TAL donors", color = "Status", title = "C. Donor-level TAL expression") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold", size = 10),
        legend.position = "bottom")

matrix_dt <- copy(interp)
matrix_dt[, magma_support := ifelse(magma_p < 5e-8, "MAGMA+", "MAGMA")]
matrix_dt[, tal_rank_support := ifelse(TAL_rank == 1, "TAL rank 1", "lower")]
matrix_dt[, donor_support := ifelse(TAL_donor_detection_fraction >= 0.75, "donor 3/4+", "lower")]
matrix_dt[, specificity_support := fifelse(grepl("^strong", specificity_class), "strong", fifelse(grepl("^moderate", specificity_class), "moderate", "broad"))]
matrix_dt[, program_support := fifelse(TAL_program_rho >= 0.1, "positive", fifelse(TAL_program_rho > 0, "weak+", "not positive"))]
matrix_dt[, manuscript_role_short := fcase(
  manuscript_role == "representative_TAL_gene", "representative",
  manuscript_role == "TAL_transport_candidate", "transport",
  manuscript_role == "calcium_ion_handling_candidate", "ion handling",
  manuscript_role == "calcium_sensing_candidate", "Ca sensing",
  manuscript_role == "supporting_context_gene", "supporting",
  manuscript_role == "broad_epithelial_context_gene", "broad epithelial",
  default = manuscript_role
)]
long_matrix <- melt(
  matrix_dt[, .(gene, magma_support, tal_rank_support, donor_support, specificity_support, program_support, manuscript_role_short)],
  id.vars = "gene",
  variable.name = "evidence_axis",
  value.name = "evidence_value"
)
long_matrix[, gene := factor(gene, levels = rev(gene_order))]
long_matrix[, evidence_axis := factor(evidence_axis, levels = c(
  "magma_support", "tal_rank_support", "donor_support", "specificity_support", "program_support", "manuscript_role_short"
))]
axis_labels <- c(
  magma_support = "MAGMA",
  tal_rank_support = "TAL rank",
  donor_support = "Donor",
  specificity_support = "Specificity",
  program_support = "TAL program",
  manuscript_role_short = "Role"
)
fill_values <- c(
  "MAGMA+" = "#2c7a7b", "MAGMA" = "#8ab6a4",
  "TAL rank 1" = "#2c7a7b", "lower" = "#d8d8d8",
  "donor 3/4+" = "#2c7a7b",
  "strong" = "#2c7a7b", "moderate" = "#77a65d", "broad" = "#b78a3b",
  "positive" = "#77a65d", "weak+" = "#d0b35f", "not positive" = "#c9c9c9",
  "representative" = "#2c7a7b",
  "transport" = "#3d8f8f",
  "ion handling" = "#7d6fb2",
  "Ca sensing" = "#9c6b43",
  "supporting" = "#6f8f4e",
  "broad epithelial" = "#b78a3b"
)
panel_d <- ggplot(long_matrix, aes(x = evidence_axis, y = gene, fill = evidence_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = evidence_value), size = 2.2, color = "black") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = fill_values, guide = "none") +
  labs(x = NULL, y = NULL, title = "D. Gene-centric evidence matrix") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 10),
        panel.grid = element_blank())

fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, labels = NULL, rel_heights = c(1.05, 1))
caption <- ggdraw() + draw_label(
  "Draft Figure 3. Expression patterns were evaluated in audited GSE231569 single-nucleus data. These analyses support TAL-associated expression context but do not establish causal mediation.",
  x = 0.01, hjust = 0, size = 8
)
fig_all <- plot_grid(fig, caption, ncol = 1, rel_heights = c(1, 0.055))

ggsave("results/figures/figure3_p1_gene_evidence_draft.pdf", fig_all, width = 13, height = 9.5, units = "in", device = "pdf")
ggsave("results/figures/figure3_p1_gene_evidence_draft.png", fig_all, width = 13, height = 9.5, units = "in", dpi = 220)

panel_sources <- data.table(
  panel = c("A", "B", "C", "D"),
  panel_title = c(
    "P1 genes DotPlot",
    "TAL specificity barplot",
    "Donor-level expression",
    "Gene role evidence matrix"
  ),
  source_files = c(
    "results/tables/p1_tal_gene_celltype_summary.tsv",
    "results/tables/p1_tal_gene_specificity.tsv",
    "results/tables/p1_tal_gene_by_donor.tsv",
    "results/tables/p1_tal_gene_interpretation_summary.tsv;results/tables/p1_tal_gene_evidence.tsv"
  ),
  included_genes = c(
    "UMOD;CASR;CLDN14;CLDN10;HIBADH;PKD2",
    "UMOD;CASR;CLDN14;CLDN10;HIBADH;PKD2",
    "UMOD;CLDN10;CASR;PKD2",
    "UMOD;CASR;CLDN14;CLDN10;HIBADH;PKD2"
  ),
  main_message = c(
    "P1 genes show interpretable expression across audited renal papillary cell types with TAL-ranked expression.",
    "Specificity ratios distinguish TAL-preferential candidates from broader epithelial context genes.",
    "Representative gene expression remains observable at donor level in TAL cells.",
    "P1 genes form an interpretable role spectrum rather than a uniform TAL marker set."
  )
)
fwrite(panel_sources, "results/tables/figure3_panel_source_files.tsv", sep = "\t")

writeLines(c(
  "# Figure 3 Legend Draft",
  "",
  "## Figure 3. Gene-centric evidence for P1 TAL-associated candidate genes",
  "",
  "**A.** Dot plot showing average expression and detection fraction for six P1 candidate genes across audited GSE231569 renal papillary cell types.",
  "",
  "**B.** TAL specificity ratio for each P1 gene, comparing average expression in audited TAL cells against the highest non-TAL cell type.",
  "",
  "**C.** Donor-level TAL expression for representative genes UMOD, CLDN10, CASR, and PKD2, used to evaluate whether the signal is observable beyond single-cell-level pseudo-replication.",
  "",
  "**D.** Gene-centric evidence matrix summarizing MAGMA support, TAL rank, donor support, specificity class, TAL program correlation direction, and manuscript role.",
  "",
  "Expression patterns were evaluated in the audited GSE231569 single-nucleus dataset. These analyses support TAL-associated expression context but do not establish causal mediation, TWAS convergence, colocalization, or spatial validation."
), "docs/figure3_legend_draft.md")

case_fraction <- data.table(
  n_case = 31715L,
  n_control = 943655L,
  case_fraction = 31715 / (31715 + 943655),
  source = "reported total case-control sample size from project instruction",
  used_for_coloc = TRUE,
  notes = "A constant case fraction was used as an approximation based on reported total case-control sample size; per-SNP N may vary."
)
fwrite(case_fraction, "results/tables/gwas_case_fraction_summary.tsv", sep = "\t")
cf_value <- case_fraction$case_fraction[1]

harmonization <- data.table(
  resource = c("KSD_2025_GWAS", "GTEx_v8_kidney_cortex_eQTL", "kidney_relevant_eQTL_placeholder"),
  gwas_variant_style = c("rsID plus CHR/BP/EA/NEA in GWAS slices", "rsID plus CHR/BP/EA/NEA in GWAS slices", "rsID plus CHR/BP/EA/NEA in GWAS slices"),
  eqtl_variant_style = c("not_applicable", "missing_resource", "missing_resource"),
  gwas_build = c("GRCh37/hg19-compatible MAGMA reference", "GRCh37/hg19-compatible MAGMA reference", "GRCh37/hg19-compatible MAGMA reference"),
  eqtl_build = c("not_applicable", "unknown_until_resource_available", "unknown_until_resource_available"),
  rsid_available = c(TRUE, NA, NA),
  chr_bp_available = c(TRUE, NA, NA),
  effect_allele_available = c(TRUE, NA, NA),
  other_allele_available = c(TRUE, NA, NA),
  needs_liftover = c(FALSE, NA, NA),
  needs_rsid_mapping = c(FALSE, NA, NA),
  harmonization_status = c("GWAS_ready", "blocked_eqtl_resource_missing", "blocked_eqtl_resource_missing"),
  notes = c(
    "GWAS pilot slices contain SNP, CHR, BP, EA, NEA, BETA, SE, P, EAF and N.",
    "eQTL resource missing, so variant/build harmonization cannot be completed.",
    "Use this row for any manually acquired kidney-relevant eQTL resource until its format/build is audited."
  )
)
fwrite(harmonization, "results/tables/coloc_variant_build_harmonization_plan.tsv", sep = "\t")

software <- data.table(
  software = c("R", "data.table", "coloc"),
  required = c(TRUE, TRUE, TRUE),
  installed = c(TRUE, requireNamespace("data.table", quietly = TRUE), requireNamespace("coloc", quietly = TRUE)),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    if (requireNamespace("data.table", quietly = TRUE)) as.character(packageVersion("data.table")) else "missing",
    if (requireNamespace("coloc", quietly = TRUE)) as.character(packageVersion("coloc")) else "missing"
  ),
  status = c(
    "available",
    if (requireNamespace("data.table", quietly = TRUE)) "available" else "missing",
    if (requireNamespace("coloc", quietly = TRUE)) "available" else "installation_pending"
  ),
  notes = c(
    "Base R is available for local scripts.",
    "Used for table preparation.",
    "Required for real coloc execution; missing package does not affect current MAGMA/scRNA mainline."
  )
)
fwrite(software, "results/tables/coloc_software_status.tsv", sep = "\t")

readiness_path <- "results/tables/coloc_pilot_readiness.tsv"
if (file.exists(readiness_path)) {
  readiness <- fread(readiness_path)
  readiness[, case_fraction_available := TRUE]
  readiness[, case_fraction := cf_value]
  readiness[, eqtl_resource_available := FALSE]
  readiness[, variant_id_match_plan_ready := file.exists("results/tables/coloc_variant_build_harmonization_plan.tsv")]
  readiness[, build_match_plan_ready := file.exists("results/tables/coloc_variant_build_harmonization_plan.tsv")]
  readiness[, ready_for_coloc := FALSE]
  readiness[, blocking_reason := fifelse(
    !gwas_slice_exists,
    "gwas_slice_missing",
    "eQTL_resource_missing; variant_build_harmonization_pending"
  )]
  fwrite(readiness, readiness_path, sep = "\t")
}

spatial <- fread("results/tables/spatial_dataset_decision.tsv")
for (old_col in c("decision_deadline", "fallback_action", "current_status", "missing_files", "manual_download_required", "notes")) {
  if (old_col %in% names(spatial)) spatial[, (old_col) := NULL]
}
spatial[, current_status := decision]
spatial[, missing_files := paste(
  ifelse(matrix_available, NA, "matrix"),
  ifelse(coordinates_available, NA, "coordinates"),
  ifelse(image_available, NA, "image"),
  ifelse(metadata_available, NA, "metadata"),
  sep = ";"
)]
spatial[, missing_files := gsub("NA;|;NA|NA", "", missing_files)]
spatial[, manual_download_required := decision == "manual_download_required"]
spatial[, decision_deadline := as.character(as.Date("2026-06-17") + 5)]
spatial[, fallback_action := fifelse(
  decision == "ready_for_qc",
  "run_spatial_qc",
  fifelse(decision == "raw_reprocessing_required",
          "do_not_use_in_current_scope",
          "move_to_GSE73680_context_if_processed_spatial_files_remain_unavailable")
)]
spatial[, notes := "Spatial datasets were catalogued, but required processed matrices, coordinates and images were not yet available locally."]
spatial_out <- spatial[, .(
  dataset, current_status, missing_files, manual_download_required, decision,
  decision_deadline, fallback_action, notes,
  matrix_available, coordinates_available, image_available, metadata_available,
  processed_object_available, raw_fastq_only, analysis_ready
)]
fwrite(spatial_out, "results/tables/spatial_dataset_decision.tsv", sep = "\t")

gse73680_files <- if (dir.exists("data/raw/GSE73680")) list.files("data/raw/GSE73680", recursive = TRUE, full.names = TRUE) else character()
has_file <- function(pattern) any(grepl(pattern, basename(gse73680_files), ignore.case = TRUE))
file_decision <- data.table(
  resource_type = c("Series matrix", "processed expression matrix", "platform annotation", "sample metadata", "raw CEL or FASTQ"),
  expected_file = c("*series_matrix*.txt.gz", "*.txt/*.tsv/*.csv/*.xlsx/*.rds expression matrix", "GPL annotation table", "sample metadata / SOFT / series matrix header", "*.CEL.gz or FASTQ"),
  local_detected = c(
    has_file("series_matrix"),
    has_file("\\.(txt|tsv|csv|xlsx|rds|rda|gz)$") && !has_file("\\.(cel|fastq|fq)(\\.gz)?$"),
    has_file("GPL|annot|platform"),
    has_file("metadata|sample|soft|series"),
    has_file("\\.(cel|fastq|fq)(\\.gz)?$")
  )
)
file_decision[, download_required := !local_detected]
file_decision[, usable_for_analysis := local_detected & resource_type %in% c("Series matrix", "processed expression matrix", "sample metadata")]
file_decision[, notes := fifelse(
  local_detected,
  "Local file candidate detected; manual format and sample-label audit required before analysis.",
  "Not detected locally; acquire processed expression and metadata before disease-context analysis."
)]
fwrite(file_decision, "results/tables/gse73680_file_decision.tsv", sep = "\t")

magma_qc <- fread("results/tables/magma_qc_summary.tsv")
magma_sets <- fread("results/tables/magma_gene_set_summary.tsv")
n_loci <- nrow(fread("results/tables/phase1_2025_loci.tsv"))
get_metric <- function(metric_name) magma_qc[metric == metric_name, value][1]
top50 <- magma_sets[gene_set == "magma_top50"]
fdr <- magma_sets[gene_set == "magma_fdr05"]
suggestive <- magma_sets[gene_set == "magma_suggestive_p1e4"]
top100 <- magma_sets[gene_set == "magma_top100"]
top200 <- magma_sets[gene_set == "magma_top200"]

p1_lines <- interp[match(gene_order, gene), paste0(
  "- **", gene, "** (`", manuscript_role, "`): TAL rank = ", TAL_rank,
  ", TAL detection = ", round(TAL_pct_expressed * 100, 1), "%",
  ", donor detection fraction = ", TAL_donor_detection_fraction,
  ", specificity ratio = ", round(specificity_ratio_avg, 2),
  ". ", final_interpretation
)]

writeLines(c(
  "# Results Draft v0.2",
  "",
  "> Evidence grade: MAGMA + scRNA-supported TAL-associated cellular context. TWAS/SMR-coloc/spatial modules remain resource-limited and are not treated as biological results.",
  "",
  "## Result 1. GWAS reconstruction and MAGMA prioritization identify KSD-associated genes",
  "",
  paste0("The locked GWAS reconstruction retained ", n_loci, " lead-locus records for downstream prioritization and passed the Phase 1 QC workflow. MAGMA v1.10 was run using a GRCh37/hg19-compatible gene-location reference and 1000G EUR LD reference. The MAGMA input contained ", get_metric("snps_in_pval"), " SNPs, of which ", get_metric("snps_used_by_magma"), " were used by MAGMA and ", get_metric("snps_mapped_to_genes"), " mapped to genes."),
  "",
  paste0("MAGMA tested ", get_metric("genes_tested"), " genes and identified ", get_metric("bonferroni_significant_genes"), " Bonferroni-significant genes, ", get_metric("fdr_significant_genes"), " FDR-significant genes, and ", get_metric("suggestive_genes_p_lt_1e4"), " suggestive genes at P < 1e-4. These outputs form the genetic-prioritization backbone for the single-cell localization analyses and candidate-gene tiering."),
  "",
  "Candidate-gene tiers should be interpreted as prioritization layers rather than causal assignments. The current tiering is appropriate for organizing downstream cell-context and resource-readiness analyses, but not for claiming TWAS, coloc, or spatial support.",
  "",
  "## Result 2. MAGMA-prioritized KSD genes localize to a TAL-associated renal papillary cellular context",
  "",
  paste0("MAGMA-prioritized KSD gene sets were projected onto the audited GSE231569 renal papillary single-nucleus annotation. Across major MAGMA sets, TAL was the top-supported cellular context: MAGMA top50 showed TAL benchmark percentile ", top50$TAL_percentile, ", top100 showed ", top100$TAL_percentile, ", top200 showed ", top200$TAL_percentile, ", suggestive genes showed ", suggestive$TAL_percentile, ", and FDR-significant genes showed ", fdr$TAL_percentile, "."),
  "",
  paste0("The top50 signal remained TAL-associated under locus-balanced checks, with full and conservative TAL percentiles of ", top50$locus_balanced_TAL_percentile_full, " and ", top50$locus_balanced_TAL_percentile_conservative, ", respectively. Leave-one-locus-out analysis retained a minimum TAL percentile of ", top50$leave_one_locus_out_min_TAL_percentile, ", supporting robustness beyond a single dominant locus."),
  "",
  "This MAGMA-based result is distinct from the earlier locus-based UMOD-sensitive observation: it supports a broader MAGMA-prioritized TAL-associated cellular context rather than relying on one locus-mapped marker. It should not be written as evidence that TAL cells causally mediate KSD risk.",
  "",
  "## Result 3. P1 core candidate genes show interpretable TAL or epithelial transport-associated expression contexts",
  "",
  "The P1 candidate genes did not represent a uniform TAL-specific marker set. Instead, they formed an interpretable gene-centric spectrum ranging from representative TAL expression to epithelial transport, calcium/ion handling and broader renal epithelial contexts.",
  "",
  paste(p1_lines, collapse = "\n"),
  "",
  "Together, the six-gene P1 evidence pack supports TAL-associated or epithelial transport-associated expression contexts for prioritized KSD genes in audited single-nucleus data. This supports the MAGMA + scRNA cellular-context model, but it does not establish TWAS convergence, colocalization, spatial validation, or causal mediation.",
  "",
  "## Resource-limited modules",
  "",
  "TWAS, SMR/coloc and spatial transcriptomic analyses remain resource-dependent. Coloc pilot GWAS slices and readiness tables are prepared, but eQTL resources and final variant/build harmonization are still missing. Spatial datasets were catalogued, but required processed matrix, coordinate, image and metadata files are not locally available. These modules should be described as planned or resource-limited, not as biological results."
), "docs/results_draft_v0.2.md")

writeLines(c(
  "# Next Phase Decision Memo",
  "",
  "## Current stable mainline",
  "",
  "The stable mainline is MAGMA + scRNA-supported TAL-associated cellular context. Phase 2D P1 gene evidence is frozen at v0.1 and can support Result 3 and draft Figure 3.",
  "",
  "## Resource-limited branches",
  "",
  "- TWAS: ready/input-prepared but resource-limited.",
  "- SMR/coloc: pilot GWAS slices are ready; eQTL resource and final variant/build harmonization are missing.",
  "- Spatial: GSE206306 and GSE231630 require manual processed-file download before analysis.",
  "- GSE73680: plaque/disease-context backup line is prepared at resource-manifest level; local expression files are missing.",
  "",
  "## Recommended next route",
  "",
  "Proceed with option D: run manuscript draft v0.1 in parallel with resource acquisition for coloc/TWAS/spatial and GSE73680 disease-context backup. Do not delay Results 1-3 on external-resource blockers.",
  "",
  "## Claim boundary",
  "",
  "Allowed: MAGMA-prioritized KSD genes and P1 core genes support a TAL-associated renal papillary single-cell expression context.",
  "",
  "Not allowed yet: TWAS support, colocalization, spatial validation, or causal mediation."
), "docs/next_phase_decision_memo.md")

writeLines(c(
  "# Figure Plan v0.1",
  "",
  "## Figure 1. Workflow and evidence levels",
  "Post-GWAS workflow from GWAS reconstruction to MAGMA prioritization, audited single-nucleus localization, P1 gene evidence, and resource-limited validation branches.",
  "",
  "## Figure 2. MAGMA to GSE231569 TAL localization",
  "MAGMA gene-set projection, TAL benchmark percentiles, locus-balanced robustness, and leave-one-locus-out stability.",
  "",
  "## Figure 3. P1 gene evidence pack",
  "Current draft: `results/figures/figure3_p1_gene_evidence_draft.pdf`. Panels cover DotPlot, TAL specificity, donor-level expression, and gene role evidence matrix.",
  "",
  "## Figure 4. Pending validation branch",
  "Reserved for spatial validation or GSE73680 plaque/disease-context expression if processed resources become available and pass QC."
), "docs/figure_plan_v0.1.md")

message("wrote Phase 2D freeze, Figure 3 draft, Results v0.2, coloc blockers, spatial/GSE73680 decisions, and next-phase docs")
