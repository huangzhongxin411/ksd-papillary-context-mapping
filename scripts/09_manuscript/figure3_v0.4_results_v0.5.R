suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  strong = "#3E6672",
  moderate = "#7F9AA3",
  light = "#D7E0E3",
  highlight = "#F1F5F6",
  broad = "#B08A45",
  negative = "#C9C9C9",
  donor = "#56616A",
  healthy = "#4C72B0",
  stone = "#B35C44",
  border = "#808080"
)

gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
cell_label_map <- c(
  Endothelial = "Endothelial",
  Collecting_duct_principal = "Collecting duct",
  Fibroblast_stromal = "Fibroblast",
  Perivascular_mural_like = "Perivascular",
  Injured_undifferentiated_epithelial = "Injured epithelial",
  Loop_of_Henle_TAL = "Loop/TAL"
)

celltype <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
specificity <- fread("results/tables/p1_tal_gene_specificity.tsv")
donor <- fread("results/tables/p1_tal_gene_by_donor.tsv")
interp <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
magma_qc <- fread("results/tables/magma_qc_summary.tsv")
magma_sets <- fread("results/tables/magma_gene_set_summary.tsv")

celltype <- celltype[cell_type %in% names(cell_label_map)]
celltype[, cell_label := cell_label_map[cell_type]]
celltype[, cell_label := factor(cell_label, levels = cell_label_map)]
celltype[, gene := factor(gene, levels = rev(gene_order))]
tal_x <- which(levels(celltype$cell_label) == "Loop/TAL")

panel_a <- ggplot(celltype, aes(cell_label, gene)) +
  annotate("rect", xmin = tal_x - 0.5, xmax = tal_x + 0.5, ymin = -Inf, ymax = Inf,
           fill = pal$highlight, alpha = 0.55) +
  geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21,
             color = pal$border, stroke = 0.22) +
  scale_size_continuous(range = c(1.4, 6.4), labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_gradient(low = "#F5F7F8", high = pal$strong) +
  labs(x = NULL, y = NULL, size = "% detected", fill = "Mean expression",
       title = "A. P1 candidate genes across audited cell types") +
  theme_bw(base_size = 9.8) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 10.8),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 8.7),
        legend.text = element_text(size = 8.1),
        legend.key.width = unit(0.28, "cm"),
        legend.key.height = unit(0.38, "cm"))

specificity[, gene := factor(gene, levels = gene_order)]
specificity[, log2_ratio := log2(specificity_ratio_avg)]
specificity[, specificity_short := fcase(
  specificity_class == "strong_TAL_preferential", "Strong",
  specificity_class == "moderate_TAL_preferential", "Moderate",
  default = "Broad"
)]
panel_b <- ggplot(specificity, aes(gene, log2_ratio, fill = specificity_short)) +
  geom_col(width = 0.68, color = "#6F6F6F", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#A8A8A8", linewidth = 0.28) +
  geom_text(aes(label = sprintf("%.1f", specificity_ratio_avg)), vjust = -0.4, size = 2.55, color = "#303030") +
  scale_fill_manual(values = c(Strong = pal$strong, Moderate = pal$moderate, Broad = pal$broad)) +
  labs(x = NULL, y = "log2(TAL specificity ratio)", fill = "Specificity class",
       title = "B. TAL specificity separates gene roles") +
  theme_bw(base_size = 9.8) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, face = "italic"),
        plot.title = element_text(face = "bold", size = 10.8),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

donor_genes <- c("UMOD", "CLDN10", "CASR", "PKD2")
donor_tal <- donor[is_TAL == TRUE & gene %in% donor_genes]
donor_tal[, gene := factor(gene, levels = donor_genes)]
median_dt <- donor_tal[, .(median_expression = median(avg_expression, na.rm = TRUE)), by = gene]
panel_c <- ggplot(donor_tal, aes(gene, avg_expression)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, fill = pal$light, alpha = 0.28,
               color = "#8A8A8A", linewidth = 0.25) +
  geom_crossbar(data = median_dt, aes(x = gene, y = median_expression, ymin = median_expression, ymax = median_expression),
                width = 0.42, inherit.aes = FALSE, color = "#303030", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.09, height = 0),
             size = 2.95, alpha = 0.98, color = pal$donor) +
  labs(x = NULL, y = "Mean expression in TAL donors",
       title = "C. Donor-level TAL expression of representative genes",
       subtitle = "Each point = one donor-level average; descriptive only") +
  theme_bw(base_size = 9.8) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, face = "italic"),
        plot.title = element_text(face = "bold", size = 10.8),
        plot.subtitle = element_text(size = 8.2, color = "#555555"),
        legend.position = "none",
        panel.grid.minor = element_blank())

matrix_dt <- copy(interp)
matrix_dt[, magma_support := "+"]
matrix_dt[, tal_rank_support := ifelse(TAL_rank == 1, "Rank 1", "Lower")]
matrix_dt[, donor_support := ifelse(TAL_donor_detection_fraction >= 0.75, "3/4 donors", "<3/4 donors")]
matrix_dt[, specificity_support := fcase(grepl("^strong", specificity_class), "Strong",
                                         grepl("^moderate", specificity_class), "Mod",
                                         default = "Broad")]
matrix_dt[, program_support := fifelse(TAL_program_rho >= 0.1, "+", fifelse(TAL_program_rho > 0, "weak", "None"))]
matrix_dt[, role_short := fcase(
  manuscript_role == "representative_TAL_gene", "Rep",
  manuscript_role == "TAL_transport_candidate", "Transport",
  manuscript_role == "calcium_ion_handling_candidate", "Ion",
  manuscript_role == "calcium_sensing_candidate", "Ca",
  manuscript_role == "supporting_context_gene", "Support",
  default = "Broad"
)]
long_matrix <- melt(matrix_dt[, .(gene, magma_support, tal_rank_support, donor_support,
                                  specificity_support, program_support, role_short)],
                    id.vars = "gene", variable.name = "axis", value.name = "value")
long_matrix[, gene := factor(gene, levels = rev(gene_order))]
long_matrix[, axis := factor(axis, levels = c("magma_support", "tal_rank_support", "donor_support",
                                              "specificity_support", "program_support", "role_short"))]
axis_labels <- c(magma_support = "MAGMA", tal_rank_support = "TAL rank", donor_support = "Donor",
                 specificity_support = "Specificity", program_support = "Program", role_short = "Role")
long_matrix[, color_group := fifelse(axis %in% c("magma_support", "tal_rank_support", "donor_support"), "core_support", value)]
fill_values <- c(
  core_support = pal$strong,
  Strong = pal$strong, Mod = pal$moderate, Broad = pal$broad,
  "+" = pal$strong, weak = "#B59A5B", None = pal$negative,
  Rep = pal$strong, Transport = pal$moderate, Ion = "#8A809C", Ca = "#8D6F55", Support = "#7A8064"
)
text_values <- c(core_support = "white", Strong = "white", Mod = "black", Broad = "black",
                 "+" = "white", weak = "black", None = "black",
                 Rep = "white", Transport = "black", Ion = "white", Ca = "white", Support = "white")
long_matrix[, text_color := text_values[color_group]]

panel_d <- ggplot(long_matrix, aes(axis, gene, fill = color_group)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = value, color = text_color), size = 2.7) +
  scale_x_discrete(labels = axis_labels) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_identity() +
  labs(x = NULL, y = NULL, title = "D. Gene-centric evidence matrix") +
  theme_minimal(base_size = 9.8) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 10.8),
        panel.grid = element_blank())

fig <- plot_grid(panel_a, panel_b, panel_c, panel_d, ncol = 2, align = "hv")
ggsave("results/figures/figure3_p1_gene_evidence_v0.4.pdf", fig, width = 12.5, height = 8.8, units = "in", device = "pdf")
ggsave("results/figures/figure3_p1_gene_evidence_v0.4.png", fig, width = 12.5, height = 8.8, units = "in", dpi = 240)

writeLines(c(
  "# Figure 3 Legend v0.4",
  "",
  "**Figure 3. Gene-centric single-nucleus evidence for P1 core TAL-associated KSD candidates.**",
  "(A) Normalized expression and detection frequency of six P1 candidate genes across audited GSE231569 renal papillary cell types, with the Loop/TAL compartment highlighted.",
  "(B) TAL specificity ratios separate representative TAL-associated genes, epithelial transport candidates, calcium/ion-handling candidates and broader epithelial context genes. Bar heights show log2-transformed TAL specificity ratios; labels denote raw ratios.",
  "(C) Donor-level TAL expression of representative P1 genes. Each point represents one donor-level average and is shown in a neutral color for descriptive consistency rather than inferential disease comparison.",
  "(D) Gene-centric evidence matrix summarizing MAGMA support, TAL expression rank, donor-level detection, TAL specificity, TAL program consistency and manuscript role assignment. Rep denotes representative TAL gene, Ca denotes calcium-sensing candidate and Ion denotes ion-handling candidate.",
  "Together, these analyses support a MAGMA + scRNA-based TAL-associated cellular context, but do not establish causal mediation, TWAS convergence, colocalization or spatial validation."
), "docs/figure3_legend_v0.4.md")

metric_value <- function(metric_name) magma_qc[metric == metric_name, value][1]
n_loci <- nrow(fread("results/tables/phase1_2025_loci.tsv"))
top50 <- magma_sets[gene_set == "magma_top50"]
top100 <- magma_sets[gene_set == "magma_top100"]
top200 <- magma_sets[gene_set == "magma_top200"]
suggestive <- magma_sets[gene_set == "magma_suggestive_p1e4"]
fdr <- magma_sets[gene_set == "magma_fdr05"]

writeLines(c(
  "# Results Draft v0.5",
  "",
  "> Current evidence grade: MAGMA + scRNA-supported TAL-associated cellular context. GSE73680 has passed analysis-ready checks, but disease-context conclusions remain pending formal plaque/control analysis.",
  "",
  "## Result 1. KSD GWAS reconstruction and MAGMA prioritization",
  paste0("The locked KSD GWAS reconstruction retained ", n_loci, " lead-locus records and was carried forward into MAGMA gene-level prioritization. MAGMA v1.10 tested ", metric_value("genes_tested"), " genes using a GRCh37/hg19-compatible reference and identified ", metric_value("bonferroni_significant_genes"), " Bonferroni-significant genes, ", metric_value("fdr_significant_genes"), " FDR-significant genes and ", metric_value("suggestive_genes_p_lt_1e4"), " suggestive genes at P < 1e-4."),
  "",
  "## Result 2. MAGMA-prioritized genes converge on a TAL-associated cellular context",
  paste0("MAGMA-prioritized gene sets localized most strongly to the audited Loop of Henle/TAL compartment in GSE231569, with TAL benchmark percentiles of ", top50$TAL_percentile, " for top50, ", top100$TAL_percentile, " for top100, ", top200$TAL_percentile, " for top200, ", suggestive$TAL_percentile, " for suggestive genes and ", fdr$TAL_percentile, " for FDR-significant genes. Locus-balanced and leave-one-locus-out analyses supported robustness of the top50 TAL signal."),
  "",
  "## Result 3. P1 core candidate genes define an interpretable TAL/epithelial transport evidence spectrum",
  "Six P1 core candidate genes were evaluated in the audited GSE231569 single-nucleus dataset. The genes did not behave as a uniform TAL-specific marker set. Instead, UMOD represented a TAL-associated gene, CLDN10 supported an epithelial transport pattern, CLDN14 and CASR linked the signal to ion-handling and calcium-sensing biology, HIBADH remained a supporting MAGMA-associated context gene and PKD2 represented a broader renal epithelial context. Figure 3 v0.4 summarizes this evidence across expression, specificity, donor-level descriptive consistency and gene-role assignment.",
  "",
  "## Result 4. Disease-context evaluation in GSE73680",
  "GSE73680 has now passed analysis-ready criteria after RAW tar extraction, Agilent FEATURES-block parsing, expression matrix construction, metadata curation, gene mapping, expression QC and P1 gene availability checks. This supports moving GSE73680 into a formal plaque/control disease-context analysis module, but it must not yet be written as a disease-context result until that analysis is run and interpreted.",
  "",
  "## Boundary",
  "The current manuscript-ready claim remains MAGMA + scRNA-supported TAL-associated cellular context. TWAS, SMR/coloc and spatial analyses remain pending/resource-limited. GSE73680 is analysis-ready, but disease-context conclusions remain pending formal downstream testing."
), "docs/results_draft_v0.5.md")

message("wrote Figure 3 v0.4 and Results v0.5")
