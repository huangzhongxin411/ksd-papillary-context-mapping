suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
cell_label_map <- c(
  Endothelial = "Endothelial",
  Collecting_duct_principal = "Collecting duct",
  Fibroblast_stromal = "Fibroblast/stromal",
  Perivascular_mural_like = "Perivascular/mural",
  Injured_undifferentiated_epithelial = "Injured epithelial",
  Loop_of_Henle_TAL = "Loop/TAL"
)
cell_order <- names(cell_label_map)

celltype <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
specificity <- fread("results/tables/p1_tal_gene_specificity.tsv")
donor <- fread("results/tables/p1_tal_gene_by_donor.tsv")
interp <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
magma_qc <- fread("results/tables/magma_qc_summary.tsv")
magma_sets <- fread("results/tables/magma_gene_set_summary.tsv")

celltype <- celltype[cell_type %in% names(cell_label_map)]
celltype[, cell_label := cell_label_map[cell_type]]
celltype[, cell_label := factor(cell_label, levels = cell_label_map[cell_order])]
celltype[, gene := factor(gene, levels = rev(gene_order))]
tal_x <- which(levels(celltype$cell_label) == "Loop/TAL")

panel_a <- ggplot(celltype, aes(x = cell_label, y = gene)) +
  annotate("rect", xmin = tal_x - 0.5, xmax = tal_x + 0.5, ymin = -Inf, ymax = Inf,
           fill = "#e9f3ef", alpha = 0.65) +
  geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = "grey25", stroke = 0.25) +
  scale_size_continuous(range = c(1.5, 6.7), labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_gradient(low = "#f6f7f1", high = "#1f6f78") +
  labs(x = NULL, y = NULL, size = "% detected", fill = "Mean expression",
       title = "A. P1 genes across audited cell types") +
  theme_bw(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 10.5),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = "grey91", linewidth = 0.2),
        panel.grid.minor = element_blank())

specificity[, gene := factor(gene, levels = gene_order)]
specificity[, log2_ratio := log2(specificity_ratio_avg)]
specificity[, specificity_short := fcase(
  specificity_class == "strong_TAL_preferential", "Strong",
  specificity_class == "moderate_TAL_preferential", "Moderate",
  specificity_class == "broad_expression_with_TAL_component", "Broad",
  default = specificity_class
)]
panel_b <- ggplot(specificity, aes(x = gene, y = log2_ratio, fill = specificity_short)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", specificity_ratio_avg)), vjust = -0.35, size = 2.6) +
  scale_fill_manual(values = c(Strong = "#2c7a7b", Moderate = "#77a65d", Broad = "#b78a3b")) +
  labs(x = NULL, y = "log2(TAL specificity ratio)", fill = "Specificity",
       title = "B. TAL specificity separates gene roles") +
  theme_bw(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "italic"),
        plot.title = element_text(face = "bold", size = 10.5),
        legend.position = "bottom",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank())

donor_genes <- c("UMOD", "CLDN10", "CASR", "PKD2")
donor_tal <- donor[is_TAL == TRUE & gene %in% donor_genes]
donor_tal[, gene := factor(gene, levels = donor_genes)]
status_values <- sort(unique(donor_tal$disease_status))
status_colors <- c(healthy_control = "#2f5597", stone_disease = "#b44a3c")
panel_c <- ggplot(donor_tal, aes(x = gene, y = avg_expression)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "#dce6df", color = "grey25", linewidth = 0.25) +
  geom_point(aes(color = disease_status), position = position_jitter(width = 0.08, height = 0),
             size = 2.4, alpha = 0.95) +
  scale_color_manual(values = status_colors[status_values], drop = FALSE) +
  labs(x = NULL, y = "Mean expression in TAL donors", color = "Donor status",
       title = "C. Donor-level TAL expression of representative genes",
       subtitle = "Each point = one donor-level average") +
  theme_bw(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "italic"),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.3, color = "grey30"),
        legend.position = "bottom",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank())

matrix_dt <- copy(interp)
matrix_dt[, magma_support := "+"]
matrix_dt[, tal_rank_support := ifelse(TAL_rank == 1, "Rank 1", "Lower")]
matrix_dt[, donor_support := ifelse(TAL_donor_detection_fraction >= 0.75, "3/4 donors", "<3/4")]
matrix_dt[, specificity_support := fcase(
  grepl("^strong", specificity_class), "Strong",
  grepl("^moderate", specificity_class), "Mod",
  default = "Broad"
)]
matrix_dt[, program_support := fifelse(TAL_program_rho >= 0.1, "+", fifelse(TAL_program_rho > 0, "weak", "-"))]
matrix_dt[, manuscript_role_short := fcase(
  manuscript_role == "representative_TAL_gene", "Rep",
  manuscript_role == "TAL_transport_candidate", "Transport",
  manuscript_role == "calcium_ion_handling_candidate", "Ion",
  manuscript_role == "calcium_sensing_candidate", "Ca",
  manuscript_role == "supporting_context_gene", "Support",
  manuscript_role == "broad_epithelial_context_gene", "Broad",
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
  program_support = "Program",
  manuscript_role_short = "Role"
)
fill_values <- c(
  "+" = "#2c7a7b", "Rank 1" = "#2c7a7b", "Lower" = "#d8d8d8",
  "3/4 donors" = "#2c7a7b", "<3/4" = "#d8d8d8",
  "Strong" = "#2c7a7b", "Mod" = "#77a65d", "Broad" = "#b78a3b",
  "weak" = "#d0b35f", "-" = "#c9c9c9",
  "Rep" = "#2c7a7b", "Transport" = "#3d8f8f", "Ion" = "#7d6fb2",
  "Ca" = "#9c6b43", "Support" = "#6f8f4e"
)
panel_d <- ggplot(long_matrix, aes(x = evidence_axis, y = gene, fill = evidence_value)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = evidence_value), size = 2.55, color = "black") +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = fill_values, guide = "none") +
  labs(x = NULL, y = NULL, title = "D. Gene-centric evidence matrix") +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 10.5),
        panel.grid = element_blank())

fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, align = "hv", rel_widths = c(1, 1))
ggsave("results/figures/figure3_p1_gene_evidence_v0.2.pdf", fig, width = 12.5, height = 8.8, units = "in", device = "pdf")
ggsave("results/figures/figure3_p1_gene_evidence_v0.2.png", fig, width = 12.5, height = 8.8, units = "in", dpi = 240)

fwrite(data.table(
  check_item = c(
    "bottom_caption_removed",
    "cell_type_labels_shortened",
    "Loop_TAL_column_highlighted",
    "panel_A_legend_titles_updated",
    "panel_B_log2_specificity_ratio",
    "panel_C_representative_gene_title",
    "panel_C_donor_point_explained",
    "panel_D_text_simplified",
    "gene_names_italicized",
    "claim_boundary_in_legend"
  ),
  status = "done",
  notes = c(
    "No in-figure draft caption retained.",
    "Long audited labels mapped to compact labels.",
    "Loop/TAL column has pale green background.",
    "Mean expression and % detected used.",
    "Dashed line marks ratio = 1 on log2 scale.",
    "Panel title states representative genes.",
    "Subtitle states each point is one donor-level average; donor status legend includes observed statuses.",
    "Matrix labels shortened to +, Rank 1, 3/4 donors, Strong/Mod/Broad, +/weak/-, and compact roles.",
    "Axis gene labels use italic face.",
    "Legend v0.2 states no causal/TWAS/coloc/spatial claim."
  )
), "results/tables/figure3_review_checklist.tsv", sep = "\t")

writeLines(c(
  "# Figure 3 Legend v0.2",
  "",
  "Figure 3. Gene-centric single-nucleus evidence for P1 core TAL-associated KSD candidates.",
  "",
  "(A) Dot plot showing normalized expression and detection frequency of six P1 candidate genes across audited GSE231569 renal papillary cell types. The Loop of Henle/TAL compartment is highlighted as the major epithelial transport context.",
  "",
  "(B) TAL specificity ratios distinguish representative TAL-associated genes, epithelial transport candidates, calcium/ion-handling candidates and broader renal epithelial context genes. Ratios are shown on a log2 scale; the dashed line indicates a TAL/max non-TAL ratio of 1.",
  "",
  "(C) Donor-level TAL expression of representative P1 genes. Representative genes were selected to cover TAL representative, transport, calcium-sensing and broad epithelial context roles. Each point represents one donor-level average within TAL cells.",
  "",
  "(D) Gene-centric evidence matrix summarizing MAGMA support, TAL expression rank, donor-level detection, TAL specificity, TAL program consistency and manuscript role. Matrix labels are abbreviated as follows: Mod, moderate; Rep, representative; Ca, calcium-sensing; Program +, positive TAL-program correlation; Program weak, weak positive correlation.",
  "",
  "These analyses support a MAGMA + scRNA-based TAL-associated cellular context, but do not establish causal mediation, TWAS convergence, colocalization or spatial validation."
), "docs/figure3_legend_v0.2.md")

metric_value <- function(metric_name) magma_qc[metric == metric_name, value][1]
n_loci <- nrow(fread("results/tables/phase1_2025_loci.tsv"))
top50 <- magma_sets[gene_set == "magma_top50"]
top100 <- magma_sets[gene_set == "magma_top100"]
top200 <- magma_sets[gene_set == "magma_top200"]
suggestive <- magma_sets[gene_set == "magma_suggestive_p1e4"]
fdr <- magma_sets[gene_set == "magma_fdr05"]

writeLines(c(
  "# Results Draft v0.3",
  "",
  "> Evidence grade: MAGMA + scRNA-supported TAL-associated cellular context. TWAS/SMR-coloc/spatial modules remain pending/resource-limited and are not biological results.",
  "",
  "## Result 1. GWAS reconstruction and MAGMA prioritization identify KSD-associated genes",
  "",
  paste0("The locked KSD GWAS reconstruction retained ", n_loci, " lead-locus records and was carried forward into MAGMA gene-level prioritization. MAGMA v1.10 was run with a GRCh37/hg19-compatible gene-location reference and 1000G EUR LD reference. The MAGMA input contained ", metric_value("snps_in_pval"), " SNPs, of which ", metric_value("snps_used_by_magma"), " were used and ", metric_value("snps_mapped_to_genes"), " mapped to genes."),
  "",
  paste0("MAGMA tested ", metric_value("genes_tested"), " genes and identified ", metric_value("bonferroni_significant_genes"), " Bonferroni-significant genes, ", metric_value("fdr_significant_genes"), " FDR-significant genes and ", metric_value("suggestive_genes_p_lt_1e4"), " suggestive genes at P < 1e-4. Candidate-gene tiers are therefore used as prioritization layers, not causal assignments."),
  "",
  "## Result 2. MAGMA-prioritized KSD genes localize to a TAL-associated renal papillary cellular context",
  "",
  paste0("MAGMA-prioritized gene sets were projected onto the audited GSE231569 renal papillary single-nucleus annotation. TAL was the top-supported cellular context across MAGMA sets: top50 TAL percentile = ", top50$TAL_percentile, ", top100 = ", top100$TAL_percentile, ", top200 = ", top200$TAL_percentile, ", suggestive = ", suggestive$TAL_percentile, " and FDR-significant = ", fdr$TAL_percentile, "."),
  "",
  paste0("The top50 signal remained TAL-associated under locus-balanced analysis, with full and conservative TAL percentiles of ", top50$locus_balanced_TAL_percentile_full, " and ", top50$locus_balanced_TAL_percentile_conservative, ". Leave-one-locus-out analysis retained a minimum TAL percentile of ", top50$leave_one_locus_out_min_TAL_percentile, ", supporting robustness beyond a single dominant locus."),
  "",
  "This MAGMA-based result should be distinguished from the earlier locus-based UMOD-sensitive observation. It supports a broader MAGMA-prioritized TAL-associated cellular context, but not causal mediation through TAL cells.",
  "",
  "## Result 3. P1 core candidate genes show interpretable TAL, epithelial transport and calcium-handling expression contexts",
  "",
  "To move beyond gene-set-level localization, we next evaluated six P1 core candidate genes in the audited GSE231569 single-nucleus dataset. These genes showed heterogeneous but interpretable expression patterns. UMOD served as a representative TAL-associated gene, whereas CLDN10 showed an epithelial transport-associated pattern with a TAL component. CLDN14 and CASR linked the TAL-associated genetic signal to ion-handling and calcium-sensing biology. HIBADH was retained as a supporting MAGMA-associated context gene, while PKD2 showed a broader renal epithelial expression pattern rather than strict TAL specificity. Together, these analyses indicate that P1 candidates form an interpretable TAL/epithelial transport gene spectrum rather than a uniform TAL-specific marker set.",
  "",
  "The corresponding Figure 3 v0.2 summarizes four complementary evidence layers: audited cell-type expression, TAL specificity ratio, donor-level TAL expression for representative genes and a gene-centric evidence matrix. These results support MAGMA + scRNA-based cellular context and do not establish TWAS convergence, colocalization, spatial validation or causality.",
  "",
  "## Resource-limited modules",
  "",
  "TWAS, SMR/coloc and spatial transcriptomic analyses remain pending/resource-limited. They should be described as resource-readiness or future validation branches until the required prediction models, eQTL resources or spatial processed files are available."
), "docs/results_draft_v0.3.md")

message("wrote Figure 3 v0.2, legend v0.2, review checklist, and Results v0.3")
