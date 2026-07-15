#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

root <- normalizePath(getwd())
out_main <- file.path(root, "figures_v3.3")
out_supp <- file.path(root, "supplementary_figures_v3.3")
dir.create(out_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_supp, recursive = TRUE, showWarnings = FALSE)

COL <- c(
  teal = "#205B66", deep = "#123F4A", bluegrey = "#7897A1",
  gold = "#B5964A", terra = "#A55846", rose = "#D5B2A7",
  olive = "#697A5B", ink = "#273238", grey = "#68757B",
  pale = "#E8EDEE", light = "#F5F7F7", white = "#FFFFFF"
)
module_cols <- c(
  "Top 50" = COL["bluegrey"], "Top 100" = COL["terra"],
  "Bonferroni" = COL["teal"], "FDR 0.05" = COL["deep"],
  "Suggestive" = COL["gold"]
)
module_cols <- setNames(unname(module_cols), c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive"))

theme_journal <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(color = COL["ink"]),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL["grey"], hjust = 0),
      plot.caption = element_text(size = base_size - 0.5, color = COL["grey"], hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, color = COL["ink"]),
      strip.background = element_rect(fill = COL["pale"], color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.2),
      legend.title = element_text(face = "bold", size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.size = unit(0.32, "cm"),
      panel.grid.major = element_line(color = "#E8ECEC", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_pair <- function(plot, pdf_path, png_path, width, height) {
  pdf_device <- function(filename, ...) grDevices::pdf(filename, ..., family = "Helvetica", useDingbats = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = pdf_device, bg = "white")
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

panel_title <- function(letter, title) paste0(letter, "  ", title)
read_tsv <- function(path) fread(path, sep = "\t", na.strings = c("NA", ""), data.table = TRUE)
read_gz <- function(path) fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", na.strings = c("NA", ""))

# Figure 1 ------------------------------------------------------------------
gwas_path <- file.path(root, "source_data/figures/phase1_step4_gwas_plot_source_data.tsv.gz")
gwas <- read_gz(gwas_path)
gwas <- gwas[is.finite(P) & P > 0 & P <= 1 & is.finite(CHR) & is.finite(BP)]
gwas[, `:=`(CHR = as.integer(CHR), BP = as.numeric(BP), neglog10p = -log10(P))]
chr_sizes <- gwas[, .(chr_max = max(BP)), by = CHR][order(CHR)]
chr_sizes[, offset := shift(cumsum(as.numeric(chr_max)), fill = 0)]
gwas <- merge(gwas, chr_sizes[, .(CHR, offset)], by = "CHR", all.x = TRUE)
gwas[, cum_pos := BP + offset]
axis_pos <- gwas[, .(axis = mean(range(cum_pos))), by = CHR]
gwas[, chromosome_group := factor(CHR %% 2)]
threshold <- -log10(5e-8)

p1 <- ggplot(gwas, aes(cum_pos, neglog10p, color = chromosome_group)) +
  geom_point(size = 0.28, alpha = 0.62, stroke = 0) +
  geom_hline(yintercept = threshold, color = COL["terra"], linetype = "dashed", linewidth = 0.55) +
  annotate("label", x = max(gwas$cum_pos) * 0.995, y = threshold + 1.4,
           label = "P = 5 x 10^-8", hjust = 1, size = 2.8,
           family = "Helvetica", color = COL["terra"], fill = "white", linewidth = 0) +
  scale_color_manual(values = c("0" = unname(COL["teal"]), "1" = "#34495E"), guide = "none") +
  scale_x_continuous(breaks = axis_pos$axis, labels = axis_pos$CHR, expand = expansion(mult = c(0.005, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "GWAS quality-control Manhattan plot",
    subtitle = "Downsampled cleaned summary statistics with all sampled genome-wide significant variants retained",
    x = "Chromosome", y = expression(-log[10](italic(P))),
    caption = "Input visualization only; no locus fine mapping or causal-gene inference."
  ) +
  theme_journal(10) +
  theme(panel.grid.major.x = element_blank(), plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 9.5), plot.caption = element_text(size = 8.5),
        axis.text.x = element_text(size = 8.2))

save_pair(
  p1,
  file.path(out_main, "Figure_1_GWAS_QC_Manhattan_polished.pdf"),
  file.path(out_main, "Figure_1_GWAS_QC_Manhattan_polished.png"),
  7.2, 4.2
)

# Figure 2 ------------------------------------------------------------------
umap <- read_gz(file.path(root, "source_data/figures/phase2_step2_umap_source_data.tsv.gz"))
markers <- read_tsv(file.path(root, "source_data/figures/phase2_step2_loop_tal_marker_dotplot_source_data.tsv"))
donor_scores <- read_tsv(file.path(root, "results/tables/phase2_step2_donor_compartment_module_scores.tsv"))
random_dist <- read_gz(file.path(root, "source_data/figures/phase2_step3_matched_random_distribution_source_data.tsv.gz"))
driver_summary <- read_tsv(file.path(root, "results/tables/phase2_step4_original_vs_driver_removed_summary.tsv"))

comp_labels <- c(
  Collecting_duct_principal = "Collecting duct", Endothelial = "Endothelial",
  Fibroblast_stromal = "Fibroblast/stromal",
  Injured_undifferentiated_epithelial = "Injured epithelial",
  Loop_of_Henle_TAL = "Loop/TAL", Pericyte_smooth_muscle = "Pericyte/smooth muscle"
)
comp_order <- unname(comp_labels)
comp_cols <- c(
  "Collecting duct" = COL["bluegrey"], "Endothelial" = "#587985",
  "Fibroblast/stromal" = COL["gold"], "Injured epithelial" = COL["terra"],
  "Loop/TAL" = COL["deep"], "Pericyte/smooth muscle" = COL["olive"]
)
comp_cols <- setNames(unname(comp_cols), c("Collecting duct", "Endothelial", "Fibroblast/stromal", "Injured epithelial", "Loop/TAL", "Pericyte/smooth muscle"))
mod_labels <- c(
  MAGMA_top50 = "Top 50", MAGMA_top100 = "Top 100",
  MAGMA_Bonferroni = "Bonferroni", MAGMA_FDR05 = "FDR 0.05",
  MAGMA_suggestive_p1e4 = "Suggestive"
)
as_comp <- function(x) factor(comp_labels[x], levels = comp_order)
as_mod <- function(x) factor(mod_labels[x], levels = unname(mod_labels))

resource <- data.table(
  metric = c("Nuclei", "Donors", "Broad\ncompartments", "Loop/TAL\nnuclei"),
  value = c("43,878", "4", "6", "540"), x = c(1, 2, 1, 2), y = c(2, 2, 1, 1)
)
p2a <- ggplot(resource, aes(x, y)) +
  geom_tile(width = 0.92, height = 0.76, fill = "white", color = COL["pale"], linewidth = 0.55) +
  geom_text(aes(y = y + 0.12, label = value), size = 5.1, fontface = "bold", color = COL["deep"]) +
  geom_text(aes(y = y - 0.12, label = metric), size = 2.7, color = COL["ink"]) +
  annotate("text", x = 1.5, y = 0.36, label = "Biological unit: donor x compartment", size = 2.75, fontface = "bold") +
  coord_cartesian(xlim = c(0.48, 2.52), ylim = c(0.18, 2.45), clip = "off") +
  labs(title = panel_title("A", "snRNA resource")) + theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold", size = 9.5), plot.margin = margin(3, 3, 3, 3))

set.seed(2505)
umap[, compartment_display := as_comp(compartment)]
umap_plot <- umap[sample(.N, min(.N, 25000))]
p2b <- ggplot(umap_plot, aes(UMAP_1, UMAP_2, color = compartment_display)) +
  geom_point(size = 0.12, alpha = 0.46, stroke = 0) +
  geom_point(data = umap_plot[compartment_display == "Loop/TAL"], size = 0.22, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = comp_cols, drop = FALSE, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(title = panel_title("B", "broad-compartment UMAP"), x = "UMAP 1", y = "UMAP 2") +
  theme_journal(8.2) + theme(legend.position = "right", aspect.ratio = 0.92)

marker_order <- c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16", "CASR", "CLDN14", "PKD2")
marker_plot <- markers[gene %in% marker_order]
marker_plot[, `:=`(gene = factor(gene, levels = marker_order), compartment_display = as_comp(compartment))]
marker_plot[, avg_scaled := as.numeric(scale(average_expression)), by = gene]
marker_plot[is.na(avg_scaled), avg_scaled := 0]
p2c <- ggplot(marker_plot, aes(gene, compartment_display)) +
  geom_point(aes(size = percent_expressing, fill = avg_scaled), shape = 21, color = "white", stroke = 0.25) +
  scale_size_area(max_size = 4.7, name = "% expressing", breaks = c(5, 25, 50, 75)) +
  scale_fill_gradient2(low = "#D8E0E2", mid = "white", high = COL["deep"], midpoint = 0, name = "Mean expression\n(scaled)") +
  labs(title = panel_title("C", "Loop/TAL marker context"), x = NULL, y = NULL) +
  theme_journal(8.1) + theme(axis.text.x = element_text(angle = 42, hjust = 1, face = "italic"), legend.position = "right")

heat <- donor_scores[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni")]
heat[, `:=`(module_display = as_mod(module_name), compartment_display = as_comp(compartment))]
donors <- sort(unique(heat$donor_id))
heat[, donor_display := factor(paste0("D", match(donor_id, donors)), levels = paste0("D", seq_along(donors)))]
heat[, score_z := as.numeric(scale(mean_module_score)), by = module_name]
heat[is.na(score_z), score_z := 0]
p2d <- ggplot(heat, aes(donor_display, compartment_display, fill = score_z)) +
  geom_tile(color = "white", linewidth = 0.35) + facet_wrap(~module_display, nrow = 1) +
  scale_fill_gradient2(low = "#D8E1E3", mid = "white", high = COL["deep"], midpoint = 0, name = "Module score\n(z)") +
  labs(title = panel_title("D", "donor x compartment module scores"), x = "Donor", y = NULL) +
  theme_journal(8.1) + theme(panel.grid = element_blank(), legend.position = "right")

bench <- random_dist[
  module_name %in% c("MAGMA_Bonferroni", "MAGMA_top100", "MAGMA_suggestive_p1e4") &
    donor_subset == "full_4_donors" & statistic_name == "median_loop_tal_minus_other"
]
bench[, module_display := factor(mod_labels[module_name], levels = c("Bonferroni", "Top 100", "Suggestive"))]
obs <- unique(bench[, .(module_display, observed_value)])
p2e <- ggplot(bench, aes(random_value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = COL["pale"], color = "white", linewidth = 0.15) +
  geom_density(color = COL["bluegrey"], linewidth = 0.45) +
  geom_vline(data = obs, aes(xintercept = observed_value), color = COL["terra"], linewidth = 0.75) +
  facet_wrap(~module_display, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels = number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(0, 0.18)) +
  labs(title = panel_title("E", "matched-random benchmark"), x = "Median Loop/TAL - other score", y = "Density") +
  theme_journal(8.1) + theme(legend.position = "none")

driver <- driver_summary[
  module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni", "MAGMA_suggestive_p1e4") &
    statistic_name == "median LoopTAL-minus-other across all 4 donors",
  .(module_name, original_value, driver_removed_value)
]
driver_long <- melt(driver, id.vars = "module_name", measure.vars = c("original_value", "driver_removed_value"),
                    variable.name = "score_set", value.name = "value")
driver_long[, `:=`(
  module_display = factor(mod_labels[module_name], levels = c("Top 50", "Top 100", "Bonferroni", "Suggestive")),
  score_set = factor(score_set, levels = c("original_value", "driver_removed_value"), labels = c("Original", "Driver removed"))
)]
p2f <- ggplot(driver_long, aes(score_set, value, group = module_display, color = module_display)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_color_manual(values = module_cols[c("Top 50", "Top 100", "Bonferroni", "Suggestive")], name = NULL) +
  labs(title = panel_title("F", "known-driver sensitivity"), x = NULL, y = "Median Loop/TAL - other score") +
  theme_journal(8.1) + theme(legend.position = "bottom", panel.grid.major.x = element_blank())

fig2_top <- p2a + p2b + p2c + plot_layout(widths = c(0.78, 1.28, 1.24))
fig2_bottom <- p2d + p2e + p2f + plot_layout(widths = c(1.35, 1.42, 1.08))
fig2 <- (fig2_top / fig2_bottom) +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(
    title = "Donor-level snRNA context of genetically prioritized modules",
    subtitle = "GSE231569: 43,878 nuclei from four papillary donors; Loop/TAL context assessed at donor level",
    caption = "Context mapping only; the figure does not assign a causal cell type.",
    theme = theme(
      plot.title = element_text(family = "Helvetica", face = "bold", size = 15, color = COL["ink"]),
      plot.subtitle = element_text(family = "Helvetica", size = 10, color = COL["grey"]),
      plot.caption = element_text(family = "Helvetica", size = 9, color = COL["grey"], hjust = 0)
    )
  )
save_pair(fig2,
          file.path(out_main, "Figure_2_snRNA_context_polished.pdf"),
          file.path(out_main, "Figure_2_snRNA_context_polished.png"), 11.5, 7.5)

# Figure 3 ------------------------------------------------------------------
matrix <- read_tsv(file.path(root, "supplementary_tables/Supplementary_Table_5_candidate_reporting_matrix.tsv"))
selection <- read_tsv(file.path(root, "results/tables/phase6_step4_figure3_candidate_matrix_row_selection.tsv"))
burden <- read_tsv(file.path(root, "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv"))
group_counts <- read_tsv(file.path(root, "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv"))
selected_genes <- selection[toupper(as.character(included_in_main_figure)) == "TRUE", gene]
main_rows <- matrix[match(selected_genes, gene), nomatch = NULL]

state_cols <- c(
  genetic_priority = COL["teal"], context = COL["deep"], partial = COL["gold"],
  proxy_multi = COL["terra"], proxy_one = COL["rose"], supplementary = "#B9C6CA",
  reviewed = COL["pale"], none = "#F4F6F6"
)
state_cols <- setNames(unname(state_cols), c("genetic_priority", "context", "partial", "proxy_multi", "proxy_one", "supplementary", "reviewed", "none"))
make_layers <- function(rows) {
  z <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i]
    data.table(
      gene = r$gene,
      group = r$repaired_reporting_group,
      layer = c("MAGMA", "snRNA", "TWAS proxy", "Spatial", "Bulk review"),
      state = c(
        ifelse(r$magma_status %in% c("top50", "top100", "bonferroni", "fdr05"), "genetic_priority", ifelse(r$magma_status == "suggestive", "partial", "none")),
        ifelse(r$snRNA_support_status == "strong_context_support", "context", ifelse(r$snRNA_support_status %in% c("moderate_context_support", "partial_context_support"), "partial", "none")),
        ifelse(r$twas_model_type == "multi_snp_proxy", "proxy_multi", ifelse(r$twas_model_type == "one_snp_proxy", "proxy_one", "none")),
        "supplementary",
        "reviewed"
      )
    )
  }))
  z[, `:=`(
    gene = factor(gene, levels = rev(rows$gene)),
    group = factor(group, levels = paste0("R", 1:6)),
    layer = factor(layer, levels = c("MAGMA", "snRNA", "TWAS proxy", "Spatial", "Bulk review"))
  )]
  z
}
main_layers <- make_layers(main_rows)

evidence_roles <- data.table(
  role = factor(c("MAGMA", "snRNA", "TWAS", "Spatial", "Bulk"), levels = c("MAGMA", "snRNA", "TWAS", "Spatial", "Bulk")),
  label = c("Genetic\npriority", "Cell\ncontext", "Cortex\nproxy", "Tissue\ncontext", "Disease-state\nreview"),
  state = c("genetic_priority", "context", "proxy_one", "supplementary", "reviewed")
)
p3a <- ggplot(evidence_roles, aes(role, 1, fill = state)) +
  geom_tile(width = 0.82, height = 0.55, color = COL["grey"], linewidth = 0.35) +
  geom_text(aes(label = label, color = role %in% c("MAGMA", "snRNA")), size = 2.7, lineheight = 0.95) +
  scale_fill_manual(values = state_cols, guide = "none") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = unname(COL["ink"])), guide = "none") +
  coord_cartesian(ylim = c(0.65, 1.35), clip = "off") +
  labs(title = panel_title("A", "parallel evidence roles"), x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica") + theme(plot.title = element_text(face = "bold", size = 9.5))

burden[, category_short := factor(
  category,
  levels = c("FDR-supported TWAS genes", "one-SNP proxy models", "multi-SNP proxy models"),
  labels = c("FDR supported", "One-SNP", "Multi-SNP")
)]
p3c <- ggplot(burden, aes(category_short, n_genes, fill = category_short)) +
  geom_col(width = 0.65) + geom_text(aes(label = n_genes), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("FDR supported" = unname(COL["bluegrey"]), "One-SNP" = unname(COL["rose"]), "Multi-SNP" = unname(COL["terra"])), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = panel_title("C", "Kidney_Cortex proxy models"), x = NULL, y = "Genes") +
  theme_journal(8.2) + theme(panel.grid.major.x = element_blank())

p3b <- ggplot(main_layers, aes(layer, gene, fill = state)) +
  geom_tile(color = "white", linewidth = 0.28) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = state_cols, breaks = names(state_cols), name = "Evidence role",
                    labels = c("MAGMA priority", "Strong context", "Partial context", "Multi-SNP proxy", "One-SNP proxy", "Supplementary", "Reviewed, non-upgrading", "None")) +
  labs(title = panel_title("B", "candidate evidence matrix"), subtitle = "Representative rows; the complete 235-gene matrix is paginated in Supplementary Figure S3", x = NULL, y = NULL) +
  theme_journal(8.1) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.2, face = "italic"),
        axis.text.x = element_text(size = 7.5, face = "bold"), strip.text.y = element_text(angle = 0, size = 7.5),
        legend.position = "right", legend.text = element_text(size = 6.3), plot.subtitle = element_text(size = 7.2))

group_counts[, reporting_group := factor(reporting_group, levels = paste0("R", 1:6))]
group_cols <- c(R1 = COL["teal"], R2 = COL["terra"], R3 = COL["rose"], R4 = COL["bluegrey"], R5 = COL["gold"], R6 = COL["pale"])
group_cols <- setNames(unname(group_cols), paste0("R", 1:6))
p3d <- ggplot(group_counts, aes(reporting_group, n_genes, fill = reporting_group)) +
  geom_col(width = 0.68) + geom_text(aes(label = n_genes), vjust = -0.3, size = 2.7) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = panel_title("D", "reporting-group counts"), x = NULL, y = "Genes") +
  theme_journal(8.2) + theme(panel.grid.major.x = element_blank())

p3e <- ggplot() +
  annotate("text", x = 0, y = 0.82, hjust = 0, label = "E  Interpretation", fontface = "bold", size = 3.5, color = COL["ink"]) +
  annotate("text", x = 0, y = 0.52, hjust = 0, label = "Supported", fontface = "bold", size = 2.9, color = COL["teal"]) +
  annotate("text", x = 0.34, y = 0.52, hjust = 0, label = "context-mapped genetic-priority candidates", size = 2.7, color = COL["ink"]) +
  annotate("text", x = 0, y = 0.24, hjust = 0, label = "Not inferred", fontface = "bold", size = 2.9, color = COL["terra"]) +
  annotate("text", x = 0.34, y = 0.24, hjust = 0, label = "causal genes or papilla-specific regulation", size = 2.7, color = COL["ink"]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") + theme_void(base_family = "Helvetica")

fig3 <- ((p3a | p3c) / p3b / (p3d | p3e)) +
  plot_layout(heights = c(0.9, 3.6, 1.05), widths = c(1.18, 0.82)) +
  plot_annotation(
    title = "Candidate reporting model",
    subtitle = "R1-R6 are reporting groups, not causal tiers; Kidney_Cortex TWAS is proxy annotation",
    theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 13.5, color = COL["ink"]),
                  plot.subtitle = element_text(family = "Helvetica", size = 9, color = COL["grey"]))
  )
save_pair(fig3,
          file.path(out_main, "Figure_3_candidate_reporting_model_polished.pdf"),
          file.path(out_main, "Figure_3_candidate_reporting_model_polished.png"), 7.2, 7.1)

# Figure 4 ------------------------------------------------------------------
f4a <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelA_source_data.tsv"))
f4b <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelB_source_data.tsv"))
f4c <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelC_source_data.tsv"))
f4d <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelD_source_data.tsv"))
f4e <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelE_source_data.tsv"))
f4f <- read_tsv(file.path(root, "source_data/figures/phase5_step4_Figure4_panelF_source_data.tsv"))

p4a <- ggplot() +
  annotate("rect", xmin = 0.1, xmax = 1.9, ymin = 1.35, ymax = 2.2, fill = COL["pale"], color = COL["bluegrey"], linewidth = 0.6) +
  annotate("text", x = 1, y = 1.93, label = "GSE73680", fontface = "bold", size = 4.0, color = COL["ink"]) +
  annotate("text", x = 1, y = 1.62, label = "55 samples | 29 patients", size = 3.2) +
  annotate("segment", x = 1, xend = 1, y = 1.31, yend = 1.05, arrow = arrow(length = unit(0.08, "in")), color = COL["teal"], linewidth = 0.55) +
  annotate("rect", xmin = 0.05, xmax = 0.93, ymin = 0.25, ymax = 0.95, fill = "#EDF2F3", color = COL["bluegrey"], linewidth = 0.5) +
  annotate("rect", xmin = 1.07, xmax = 1.95, ymin = 0.25, ymax = 0.95, fill = "#F4EEDF", color = COL["gold"], linewidth = 0.5) +
  annotate("text", x = 0.49, y = 0.70, label = "Control", fontface = "bold", size = 2.8) +
  annotate("text", x = 0.49, y = 0.45, label = "27 samples", size = 2.7) +
  annotate("text", x = 1.51, y = 0.70, label = "Disease context", fontface = "bold", size = 2.8) +
  annotate("text", x = 1.51, y = 0.45, label = "28 samples", size = 2.7) +
  annotate("text", x = 1, y = 0.05, label = "26 complete patient pairs", fontface = "bold", size = 2.9, color = COL["terra"]) +
  coord_cartesian(xlim = c(0, 2), ylim = c(-0.05, 2.25), clip = "off") +
  labs(title = panel_title("A", "paired bulk design")) + theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold", size = 10))

f4b <- f4b[module_name %in% c("Top 50", "FDR 0.05", "Suggestive")]
f4b[, `:=`(
  module_name = factor(module_name, levels = c("Top 50", "FDR 0.05", "Suggestive")),
  group_curated = factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))
)]
p4b <- ggplot(f4b, aes(group_curated, module_score, group = patient_id)) +
  geom_line(color = COL["bluegrey"], alpha = 0.46, linewidth = 0.3) +
  geom_point(aes(fill = group_curated), shape = 21, color = COL["ink"], stroke = 0.2, size = 1.2) +
  facet_wrap(~module_name, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(control_or_adjacent = unname(COL["bluegrey"]), plaque_or_stone_papilla = unname(COL["gold"])), guide = "none") +
  scale_x_discrete(labels = c(control_or_adjacent = "Control", plaque_or_stone_papilla = "Disease")) +
  labs(title = panel_title("B", "selected paired module scores"), x = NULL, y = "Module score") +
  theme_journal(8.5) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

f4c[, module_label := factor(module_label, levels = rev(c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive")))]
p4c <- ggplot(f4c, aes(group_coefficient, module_label)) +
  geom_vline(xintercept = 0, color = COL["ink"], linewidth = 0.35) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.16, color = COL["bluegrey"], linewidth = 0.6) +
  geom_point(aes(fill = nominal), shape = 21, color = COL["ink"], size = 2.2, stroke = 0.3) +
  scale_fill_manual(values = c("TRUE" = unname(COL["gold"]), "FALSE" = unname(COL["pale"])), guide = "none") +
  labs(title = panel_title("C", "base disease-context coefficients"), subtitle = "All module FDR values > 0.05", x = "Coefficient (95% CI)", y = NULL) +
  theme_journal(8.8)

f4d <- f4d[signature_name %in% c("injury_KIM1_LCN2_like", "epithelial_transport", "mineralization_remodeling_proxy")]
f4d[, `:=`(
  signature_label = factor(signature_label, levels = c("Injury", "Epithelial state", "Mineralization/remodeling"), labels = c("Injury", "Epithelial", "Mineralization")),
  group_curated = factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))
)]
p4d <- ggplot(f4d, aes(group_curated, signature_score, group = patient_id)) +
  geom_line(color = COL["bluegrey"], alpha = 0.46, linewidth = 0.3) +
  geom_point(aes(fill = group_curated), shape = 21, color = COL["ink"], stroke = 0.2, size = 1.2) +
  facet_wrap(~signature_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(control_or_adjacent = unname(COL["bluegrey"]), plaque_or_stone_papilla = unname(COL["gold"])), guide = "none") +
  scale_x_discrete(labels = c(control_or_adjacent = "Control", plaque_or_stone_papilla = "Disease")) +
  labs(title = panel_title("D", "paired tissue-state proxy shifts"), subtitle = "Expression scores, not cell fractions", x = NULL, y = "Proxy score") +
  theme_journal(8.5) + theme(axis.text.x = element_text(angle = 20, hjust = 1), strip.text = element_text(size = 7.6))

f4e[, `:=`(
  module_label = factor(module_label, levels = rev(c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive"))),
  signature_label = factor(signature_label, levels = c("Injury", "Epithelial state", "Stromal/ECM", "Endothelial", "Fibrosis/remodeling"))
)]
p4e <- ggplot(f4e, aes(signature_label, module_label, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.4) + geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.65) +
  scale_fill_gradient2(low = COL["bluegrey"], mid = "white", high = COL["terra"], midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(title = panel_title("E", "module-tissue-state coupling"), subtitle = "Paired disease-minus-control differences", x = NULL, y = NULL) +
  theme_journal(8.4) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "right")

f4f[, `:=`(
  model_label = factor(model_label, levels = c("Base", "Injury", "Injury + fibrosis")),
  module_label = factor(module_label, levels = c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive"))
)]
p4f <- ggplot(f4f, aes(model_label, adjusted_group_coefficient, group = module_label, color = module_label)) +
  geom_hline(yintercept = 0, color = COL["ink"], linewidth = 0.35) + geom_line(linewidth = 0.65, alpha = 0.82) + geom_point(size = 1.9) +
  scale_color_manual(values = module_cols, name = NULL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = panel_title("F", "coefficient attenuation"), subtitle = "Tissue-state sensitivity", x = NULL, y = "Disease-context coefficient") +
  theme_journal(8.8) + theme(legend.position = "bottom", legend.text = element_text(size = 7.2), legend.key.width = unit(0.38, "cm"))

design4 <- c(
  area(t = 1, l = 1, b = 1, r = 3), area(t = 1, l = 4, b = 1, r = 8), area(t = 1, l = 9, b = 1, r = 12),
  area(t = 2, l = 1, b = 2, r = 6), area(t = 2, l = 7, b = 2, r = 12),
  area(t = 3, l = 1, b = 3, r = 12)
)
fig4 <- wrap_plots(p4a, p4c, p4f, p4b, p4d, p4e, design = design4) +
  plot_annotation(
    title = "Paired bulk disease context and tissue-state sensitivity",
    subtitle = "GSE73680: 55 samples from 29 patients, including 26 complete pairs",
    caption = "Positive module coefficients did not remain FDR significant and attenuated after tissue-state adjustment.",
    theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 15, color = COL["ink"]),
                  plot.subtitle = element_text(family = "Helvetica", size = 9.5, color = COL["grey"]),
                  plot.caption = element_text(family = "Helvetica", size = 8.5, color = COL["grey"], hjust = 0))
  )
save_pair(fig4,
          file.path(out_main, "Figure_4_bulk_disease_context_polished.pdf"),
          file.path(out_main, "Figure_4_bulk_disease_context_polished.png"), 10.5, 7.8)

# Supplementary Figure S1 ---------------------------------------------------
magma_pack <- read_tsv(file.path(root, "supplementary_tables/Supplementary_Table_2_MAGMA_results_and_modules.tsv"))
ranked <- magma_pack[record_type == "ranked_MAGMA_gene"][order(magma_rank)][1:20]
sets <- magma_pack[record_type == "gene_set_summary"]
if (nrow(sets) == 0) {
  sets <- unique(magma_pack[!is.na(gene_set_name), .(gene_set_name, n_genes)])
}
qq_p <- sort(gwas$P)
qq <- data.table(expected = -log10(ppoints(length(qq_p))), observed = -log10(qq_p))
ps1a <- ggplot(qq, aes(expected, observed)) + geom_point(size = 0.35, alpha = 0.55, color = COL["teal"]) +
  geom_abline(slope = 1, intercept = 0, color = COL["terra"], linewidth = 0.55) +
  labs(title = panel_title("A", "GWAS quantile-quantile diagnostic"), x = expression(Expected~~-log[10](italic(P))), y = expression(Observed~~-log[10](italic(P)))) +
  theme_journal(9)
ps1b <- ggplot(gwas, aes(P)) + geom_histogram(bins = 70, fill = COL["bluegrey"], color = "white", linewidth = 0.15) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(title = panel_title("B", "GWAS P-value distribution"), x = "P value", y = "Variants") + theme_journal(9)
ranked[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
ps1c <- ggplot(ranked, aes(-log10(magma_p), gene_symbol)) + geom_col(fill = COL["teal"], width = 0.72) +
  labs(title = panel_title("C", "top 20 MAGMA genes by rank"), x = expression(-log[10](MAGMA~~italic(P))), y = NULL) + theme_journal(8.7) +
  theme(axis.text.y = element_text(face = "italic", size = 7.3))
set_label <- function(x) {
  x <- gsub("MAGMA_", "", x); x <- gsub("_p1e4", "", x); x <- gsub("FDR05", "FDR 0.05", x)
  x <- gsub("top50", "Top 50", x); x <- gsub("top100", "Top 100", x); tools::toTitleCase(x)
}
sets[, label := factor(set_label(gene_set_name), levels = set_label(gene_set_name))]
ps1d <- ggplot(sets, aes(label, n_genes, fill = label)) + geom_col(width = 0.68) + geom_text(aes(label = n_genes), vjust = -0.3, size = 2.7) +
  scale_fill_manual(values = unname(rep(c(COL["teal"], COL["bluegrey"], COL["terra"], COL["gold"], COL["grey"]), length.out = nrow(sets))), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = panel_title("D", "canonical MAGMA module sizes"), x = NULL, y = "Genes") + theme_journal(8.7) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

supp1 <- ((ps1a | ps1b) / (ps1c | ps1d)) +
  plot_annotation(
    title = "Supplementary Figure S1 | GWAS and MAGMA diagnostics",
    subtitle = "Input behavior and gene-priority summaries",
    caption = "Diagnostic and prioritization displays only; no fine mapping or causal-gene inference.",
    theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 14, color = COL["ink"]),
                  plot.subtitle = element_text(family = "Helvetica", size = 9, color = COL["grey"]),
                  plot.caption = element_text(family = "Helvetica", size = 8, color = COL["grey"], hjust = 0))
  )
save_pair(supp1,
          file.path(out_supp, "Supplementary_Figure_S1_GWAS_MAGMA_diagnostics_polished.pdf"),
          file.path(out_supp, "Supplementary_Figure_S1_GWAS_MAGMA_diagnostics_polished.png"), 8.5, 7.0)

message("Phase 7-Step 2 main figures and Supplementary Figure S1 created from locked source tables.")
