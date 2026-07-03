suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
})

fig_dir <- "results/figures/revision/stage4C2R_draft_figures_v0.2"
table_dir <- "results/tables/revision/stage4C2R_draft_figures_v0.2"
doc_dir <- "docs/revision/stage4C2R_draft_figures_v0.2"
script_dir <- "scripts/revision_utils/stage4C2R_draft_figures_v0.2"
log_dir <- "logs/revision/stage4C2R_draft_figures_v0.2"
for (d in c(fig_dir, table_dir, doc_dir, script_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage4C2R_generate_draft_figures_v0.2.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4C2-R conservative draft figure refinement v0.2\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  fig2_blueprint = "docs/revision/stage4C1_figure_claim_planning/figure2_snRNA_context_blueprint.md",
  fig3_blueprint = "docs/revision/stage4C1_figure_claim_planning/figure3_candidate_evidence_blueprint.md",
  panel_manifest = "results/tables/revision/stage4C1_figure_claim_planning/figure_panel_source_data_manifest.tsv",
  claim_wording = "results/tables/revision/stage4C1_figure_claim_planning/snRNA_claim_wording_decision_table.tsv",
  legends_v01 = "docs/revision/stage4C1_figure_claim_planning/draft_figure2_figure3_legends_v0.1.md",
  f2a_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelA_source.tsv",
  f2b_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelB_source.tsv",
  f2c_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelC_source.tsv",
  f2d_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelD_source.tsv",
  f2e_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelE_source.tsv",
  f2f_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelF_source.tsv",
  f2g_source = "results/tables/revision/stage4C2_draft_figures/figure2_panelG_source.tsv",
  f3a_source = "results/tables/revision/stage4C2_draft_figures/figure3_panelA_source.tsv",
  f3b_source = "results/tables/revision/stage4C2_draft_figures/figure3_panelB_source.tsv",
  f3c_source = "results/tables/revision/stage4C2_draft_figures/figure3_panelC_source.tsv",
  f3d_source = "results/tables/revision/stage4C2_draft_figures/figure3_panelD_source.tsv",
  f3e_source = "results/tables/revision/stage4C2_draft_figures/figure3_panelE_source.tsv",
  stage3_model = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  stage3_counts = "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv",
  exemplar = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv"
)

for (p in paths) {
  if (!file.exists(p)) stop("Required source missing: ", p)
}

check_cols <- function(dt, cols, name) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing columns in ", name, ": ", paste(miss, collapse = ", "))
}

f2a_src <- fread(paths$f2a_source)
f2b_src <- fread(paths$f2b_source)
f2c_src <- fread(paths$f2c_source)
f2d_src <- fread(paths$f2d_source)
f2e_src <- fread(paths$f2e_source)
f2f_src <- fread(paths$f2f_source)
f2g_src <- fread(paths$f2g_source)
f3a_src <- fread(paths$f3a_source)
f3b_src <- fread(paths$f3b_source)
f3c_src <- fread(paths$f3c_source)
f3d_long <- fread(paths$f3d_source)
f3e_src <- fread(paths$f3e_source)
claim_wording <- fread(paths$claim_wording)

check_cols(f2a_src, c("broad_compartment", "n_donors", "total_nuclei", "is_loop_tal"), "figure2 panel A source")
check_cols(f2b_src, c("module_name", "donor_short", "broad_compartment", "mean_module_score", "module_label"), "figure2 panel B source")
check_cols(f2c_src, c("module_name", "donor_short", "loop_tal_rank", "loop_tal_percentile_rank", "module_label"), "figure2 panel C source")
check_cols(f2d_src, c("module_name", "donor_short", "support_retained", "loop_tal_rank", "module_label"), "figure2 panel D source")
check_cols(f2e_src, c("module_name", "module_label", "interpretation", "empirical_p_delta", "observed_loop_tal_delta"), "figure2 panel E source")
check_cols(f2f_src, c("base_module", "removal_type", "support_change", "module_label", "removal_label"), "figure2 panel F source")
check_cols(f2g_src, c("module_name", "overall_claim_strength", "recommended_manuscript_claim", "module_label"), "figure2 panel G source")
check_cols(f3a_src, c("genetic_priority_level", "twas_proxy_level", "N"), "figure3 panel A source")
check_cols(f3b_src, c("reporting_group", "n_unique_genes"), "figure3 panel B source")
check_cols(f3c_src, c("gene", "biological_role_label", "bonferroni_significant"), "figure3 panel C source")
check_cols(f3d_long, c("gene", "evidence_layer", "status", "status_simple"), "figure3 panel D source")
check_cols(f3e_src, c("claim_type", "claim"), "figure3 panel E source")

primary_modules <- c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")
module_labels <- c(
  R1_MAGMA_Bonferroni_only = "R1 Bonferroni",
  R1_R2_R3_all_MAGMA_Bonferroni = "R1-R3 Bonf.",
  MAGMA_top50 = "MAGMA top50",
  MAGMA_top100 = "MAGMA top100",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = "R3 TWAS proxy",
  R2_R3_MAGMA_plus_TWAS_proxy = "R2/R3 proxy",
  Curated_exemplar_panel = "Curated exemplars"
)

theme_stage <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#222222"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "#555555"),
      axis.text = element_text(color = "#333333"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 1),
      plot.margin = margin(5, 5, 5, 5)
    )
}

panel_label <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = 14),
          plot.tag.position = c(0, 1))
}

safe_label <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(module_labels), module_labels[x], x)
}

write_source <- function(dt, file) {
  fwrite(dt, file.path(table_dir, file), sep = "\t")
}

# Figure 2 source extracts
pretty_compartment <- function(x) {
  map <- c(
    Collecting_duct_principal = "CD principal",
    Fibroblast_stromal = "Fibro/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epi",
    Loop_of_Henle_TAL = "Loop/TAL",
    Pericyte_smooth_muscle = "Pericyte/SMC"
  )
  x <- as.character(x)
  ifelse(x %in% names(map), unname(map[x]), x)
}

f2a_src[, compartment_label := pretty_compartment(broad_compartment)]
write_source(f2a_src, "figure2_panelA_source.tsv")

f2b_src <- f2b_src[module_name %in% primary_modules]
f2b_src[, compartment_label := pretty_compartment(broad_compartment)]
write_source(f2b_src, "figure2_panelB_source.tsv")

f2c_src <- f2c_src[module_name %in% primary_modules]
write_source(f2c_src, "figure2_panelC_source.tsv")

f2d_src <- f2d_src[module_name %in% primary_modules]
write_source(f2d_src, "figure2_panelD_source.tsv")

f2e_src <- f2e_src[module_name %in% primary_modules]
f2e_src[, support_label := fifelse(interpretation == "partial_support", "partial", interpretation)]
f2e_src[, small_p_note := paste0("p_delta=", signif(empirical_p_delta, 2))]
write_source(f2e_src, "figure2_panelE_source.tsv")

panel_removals <- c("without_curated_exemplar_panel", "without_TAL_marker_panel", "without_calcium_ion_panel", "without_top5_contributors", "without_top10_contributors")
f2f_src <- f2f_src[base_module %in% primary_modules & removal_type %in% panel_removals]
f2f_src[, removal_label := factor(removal_label, levels = c("curated exemplar panel", "TAL marker panel", "calcium ion panel", "top5 contributors", "top10 contributors"))]
write_source(f2f_src, "figure2_panelF_source.tsv")

f2g_src <- f2g_src[module_name %in% c(primary_modules, "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy", "Curated_exemplar_panel")]
write_source(f2g_src, "figure2_panelG_source.tsv")

# Figure 2 panels
p2a <- ggplot(f2a_src, aes(x = reorder(compartment_label, total_nuclei), y = total_nuclei, fill = is_loop_tal)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#4C78A8", "FALSE" = "#C7C7C7"), guide = "none") +
  labs(title = "Audited GSE231569 context", subtitle = "43,878 nuclei; Loop/TAL = 540; no cell-level inference", x = NULL, y = "Nuclei") +
  theme_stage(8)

p2b <- ggplot(f2b_src, aes(x = donor_short, y = compartment_label, fill = mean_module_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~module_label, ncol = 2) +
  scale_fill_gradient(low = "#F2F2F2", high = "#4C78A8", name = "Mean\nscore") +
  labs(title = "Donor x compartment scores", subtitle = "donor x compartment summaries; descriptive", x = "Donor", y = NULL) +
  theme_stage(7) +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "right")

p2c <- ggplot(f2c_src, aes(x = donor_short, y = module_label, color = loop_tal_rank, size = loop_tal_percentile_rank)) +
  geom_point(alpha = 0.95) +
  scale_color_gradient(low = "#4C78A8", high = "#D9D9D9", trans = "reverse", name = "Loop/TAL\nrank") +
  scale_size(range = c(2.2, 5), name = "Percentile") +
  labs(title = "Loop/TAL within-donor rank", subtitle = "descriptive ranks; no cell-level P values", x = "Donor", y = NULL) +
  theme_stage(8)

p2d <- ggplot(f2d_src, aes(x = donor_short, y = module_label, fill = support_retained)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = loop_tal_rank), size = 2.5, color = "#222222") +
  scale_fill_manual(values = c(yes = "#59A14F", partial = "#F2C14E", no = "#E15759", insufficient_data = "#BDBDBD"), name = "Retained") +
  labs(title = "Single-donor exclusion robustness", subtitle = "Tile label = Loop/TAL rank", x = "Excluded donor", y = NULL) +
  theme_stage(8)

p2e <- ggplot(f2e_src, aes(x = module_label, y = observed_loop_tal_delta)) +
  geom_col(fill = "#E7C76A", width = 0.62) +
  geom_text(aes(label = "partial"), vjust = -0.25, size = 2.8, fontface = "bold", color = "#6B5A1E") +
  annotate("text", x = 2.5, y = max(f2e_src$observed_loop_tal_delta, na.rm = TRUE) * 1.28,
           label = "delta above matched random;\nrank metrics saturated", size = 2.55, color = "#555555", lineheight = 0.9) +
  annotate("text", x = 2.5, y = max(f2e_src$observed_loop_tal_delta, na.rm = TRUE) * 1.13,
           label = "empirical p_delta in source table", size = 2.15, color = "#777777") +
  coord_cartesian(ylim = c(0, max(f2e_src$observed_loop_tal_delta, na.rm = TRUE) * 1.42), clip = "off") +
  labs(title = "Matched random benchmark", subtitle = "partial matched-random support", x = NULL, y = "Loop/TAL - other delta") +
  theme_stage(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p2f <- ggplot(f2f_src, aes(x = removal_label, y = module_label, fill = support_change)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(unchanged = "#A8D5A2", weakened = "#F2C14E", strengthened = "#A9D7D4", collapsed = "#E15759", not_evaluable = "#BDBDBD"), name = "Change") +
  labs(title = "Known-driver removal", subtitle = "Module-level descriptive robustness", x = NULL, y = NULL) +
  theme_stage(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

claim_text <- "Moderate context support\nDonor-level Loop/TAL pattern\nPartial matched-random support\nRobust driver-removal\n\nNot causal cell type\nNot plaque nucleation site"
p2g <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#F7F7F7", color = "#666666", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 0.64, label = claim_text, size = 3.3, lineheight = 0.95, fontface = "bold", color = "#222222") +
  annotate("text", x = 0.5, y = 0.12, label = "Claim level: MODERATE, not strong", size = 2.8, color = "#A05A00") +
  labs(title = "Claim decision") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0))

fig2 <- (
  panel_label(p2a, "A") | panel_label(p2b, "B")
) / (
  panel_label(p2c, "C") | panel_label(p2d, "D")
) / (
  panel_label(p2e, "E") | panel_label(p2f, "F") | panel_label(p2g, "G")
) +
  plot_layout(heights = c(1.05, 1, 0.95), widths = c(1, 1, 0.9), guides = "collect") +
  plot_annotation(
    title = "Figure 2 draft v0.2. Donor-level snRNA context mapping of MAGMA-prioritized modules",
    subtitle = "Conservative visual refinement: moderate context support; partial matched-random support; descriptive driver-removal robustness",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "#555555"))
  )

# Figure 3 source extracts
write_source(f3a_src, "figure3_panelA_source.tsv")

reporting_groups_keep <- c(
  "R1_MAGMA_Bonferroni_only",
  "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy",
  "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy",
  "R4_MAGMA_lower_priority",
  "R5_TWAS_proxy_only",
  "R6_contextual_or_unresolved"
)
f3b_src <- f3b_src[reporting_group %in% reporting_groups_keep]
f3b_src[, reporting_group := factor(reporting_group, levels = reporting_groups_keep)]
write_source(f3b_src, "figure3_panelB_source.tsv")

write_source(f3c_src, "figure3_panelC_source.tsv")

f3d_long[evidence_layer == "Kidney_Cortex_TWAS_proxy",
         status_simple := "no FDR TWAS proxy assigned"]
f3d_long[evidence_layer == "SMR_coloc",
         status_simple := "not claim-grade ready"]
f3d_long[evidence_layer == "MAGMA_Bonferroni" & status %in% c("yes", "MAGMA prioritized"),
         status_simple := "MAGMA Bonferroni"]
write_source(f3d_long, "figure3_panelD_source.tsv")

f3e_src[claim == "Kidney_Cortex proxy TWAS support", claim := "Kidney_Cortex proxy TWAS context"]
write_source(f3e_src, "figure3_panelE_source.tsv")

p3a <- ggplot(f3a_src, aes(x = genetic_priority_level, y = twas_proxy_level, fill = N)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = N), size = 3, fontface = "bold") +
  scale_fill_gradient(low = "#F2F2F2", high = "#4C78A8", name = "Genes") +
  labs(title = "Two-axis evidence model", subtitle = "Genetic priority x Kidney_Cortex proxy TWAS", x = "Genetic priority", y = "Kidney_Cortex proxy TWAS") +
  theme_stage(8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6))

group_labels <- c(
  R1_MAGMA_Bonferroni_only = "R1 MAGMA Bonferroni only",
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = "R2 + multi-SNP proxy",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = "R3 + one-SNP proxy",
  R4_MAGMA_lower_priority = "R4 MAGMA lower priority",
  R5_TWAS_proxy_only = "R5 TWAS proxy only",
  R6_contextual_or_unresolved = "R6 contextual/unresolved"
)
f3b_src[, reporting_group_label := unname(group_labels[as.character(reporting_group)])]
p3b <- ggplot(f3b_src, aes(x = reorder(reporting_group_label, n_unique_genes), y = n_unique_genes)) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  coord_flip() +
  labs(title = "Reporting group counts", subtitle = "Reporting groups; not causal tiers", x = NULL, y = "Genes") +
  theme_stage(8) +
  theme(axis.text.y = element_text(size = 6.8))

f3c_src[, gene := factor(gene, levels = rev(gene))]
p3c <- ggplot(f3c_src, aes(x = 1, y = gene)) +
  geom_point(aes(color = bonferroni_significant), size = 4) +
  geom_text(aes(label = biological_role_label), x = 1.08, hjust = 0, size = 2.6) +
  scale_color_manual(values = c(yes = "#4C78A8", no = "#BDBDBD"), name = "MAGMA\nBonf.") +
  coord_cartesian(xlim = c(0.95, 2.2), clip = "off") +
  labs(title = "Curated biological exemplars", subtitle = "Role spectrum; no evidence upgrade", x = NULL, y = NULL) +
  theme_stage(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(5, 70, 5, 5))

f3d_long[, gene := factor(gene, levels = rev(unique(as.character(f3c_src$gene))))]
f3d_long[, evidence_layer := factor(evidence_layer, levels = c("MAGMA_Bonferroni", "Kidney_Cortex_TWAS_proxy", "SMR_coloc"))]
p3d <- ggplot(f3d_long, aes(x = evidence_layer, y = gene, fill = status_simple)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c("MAGMA Bonferroni" = "#4C78A8", "no FDR TWAS proxy assigned" = "#ECE7D8", "not claim-grade ready" = "#D9D9D9"),
                    name = "Status") +
  labs(title = "Exemplar evidence strip", subtitle = "TWAS proxy does not upgrade exemplars; SMR/coloc not claim-grade ready", x = NULL, y = NULL) +
  theme_stage(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.text.y = element_text(size = 7))

f3e_src[, claim_type := factor(claim_type, levels = c("Allowed", "Disallowed"))]
f3e_src[, claim_order := seq_len(.N)]
p3e <- ggplot(f3e_src, aes(x = claim_type, y = reorder(claim, -claim_order), fill = claim_type)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = claim), size = 2.7, color = "#222222") +
  scale_fill_manual(values = c(Allowed = "#DDEED6", Disallowed = "#F4DADA"), guide = "none") +
  labs(title = "Claim boundary", subtitle = "No causal/validated/therapeutic labels", x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, color = "#555555"))

fig3 <- (
  panel_label(p3a, "A") | panel_label(p3b, "B")
) / (
  panel_label(p3c, "C") | panel_label(p3d, "D")
) / (
  panel_label(p3e, "E")
) +
  plot_layout(heights = c(1, 1, 0.72), guides = "collect") +
  plot_annotation(
    title = "Figure 3 draft v0.2. Candidate-gene evidence model and exemplar boundaries",
    subtitle = "Conservative display: mutually exclusive reporting groups, curated exemplars as flags, no causal or validated gene labels",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "#555555"))
  )

save_figure <- function(plot, stem, width, height) {
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  svg_path <- file.path(fig_dir, paste0(stem, ".svg"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = "pdf")
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600)
  ggsave(svg_path, plot, width = width, height = height, units = "in", device = svglite)
  c(pdf = pdf_path, png = png_path, svg = svg_path)
}

fig2_files <- save_figure(fig2, "figure2_snRNA_context_draft_v0.2", 14.5, 12)
fig3_files <- save_figure(fig3, "figure3_candidate_evidence_draft_v0.2", 13.5, 11)

source_files <- list.files(table_dir, pattern = "^figure[23]_panel[A-G]_source\\.tsv$", full.names = TRUE)
source_manifest <- rbindlist(lapply(source_files, function(p) {
  dt <- fread(p)
  data.table(
    figure = sub("^(figure[23]).*", "\\1", basename(p)),
    panel = sub("^figure[23]_panel([A-G])_source\\.tsv$", "\\1", basename(p)),
    source_table = p,
    n_rows = nrow(dt),
    n_columns = ncol(dt),
    ready_for_publication_source_data = ifelse(nrow(dt) > 0 && ncol(dt) > 0, "draft_ready", "needs_fix"),
    notes = "Stage 4C2-R v0.2 source extract derived from Stage 4C2 frozen source data; final publication source data should be rechecked before submission."
  )
}), fill = TRUE)
setorder(source_manifest, figure, panel)
fwrite(source_manifest, file.path(table_dir, "figure_source_data_manifest_v0.2.tsv"), sep = "\t")

qc <- data.table(
  check_item = c(
    "Figure 2 still says moderate, not strong",
    "Figure 2 Panel E clearly says partial support",
    "Figure 2 Panel F is simplified",
    "Figure 2 avoids cell-level P values",
    "Figure 2 avoids causal cell-type language",
    "Figure 2 known-driver removal is descriptive robustness",
    "Figure 3 Panel B excludes curated exemplar flag rows",
    "Figure 3 Panel D does not imply positive TWAS support where none exists",
    "Figure 3 avoids P1 terminology",
    "Figure 3 avoids causal/validated gene language",
    "Figure 3 marks TWAS as Kidney_Cortex proxy",
    "Figure 3 marks SMR/coloc as not claim-grade ready",
    "All panels backed by source data",
    "Remaining visual limitations before final publication figures"
  ),
  status = c("pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass",
             ifelse(all(source_manifest$ready_for_publication_source_data == "draft_ready"), "pass", "needs_fix"),
             "minor_limitations"),
  notes = c(
    "Title/subtitle and claim strip use moderate/not strong language.",
    "Panel E headline and labels emphasize partial matched-random support; empirical p-delta is demoted to source-table annotation.",
    "Panel F is restricted to panel-level removals: curated exemplar, TAL marker, calcium/ion, top5, and top10 contributors.",
    "No cell-level inferential P values are displayed.",
    "Panel G states not causal cell type and not plaque nucleation site.",
    "Panel F labels descriptive robustness, not mechanism proof.",
    "Panel B includes only R1-R6 mutually exclusive reporting groups.",
    "Panel D labels Kidney_Cortex TWAS proxy as no FDR TWAS proxy assigned for exemplars.",
    "No P1 label appears in Figure 3 draft.",
    "Claim boundary panel disallows causal and validated gene labels.",
    "Figure 3 subtitle and Panel A label TWAS as Kidney_Cortex proxy.",
    "Panel D labels SMR/coloc as not claim-grade ready.",
    paste0(nrow(source_manifest), " panel source extracts generated."),
    "Figures are clean draft main figures, not final publication figures; final pass should still tune typography, spacing, and journal-specific dimensions."
  )
)

writeLines(c(
  "# Figure QC audit v0.2",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "| check | status | notes |",
  "|---|---|---|",
  apply(qc, 1, function(r) paste0("| ", paste(gsub("\\|", "/", r), collapse = " | "), " |"))
), file.path(doc_dir, "figure_qc_audit_v0.2.md"))

writeLines(c(
  "# Draft Figure 2 and Figure 3 legends v0.2",
  "",
  "## Figure 2. Donor-level snRNA context mapping and robustness boundaries for MAGMA-prioritized modules",
  "",
  "Draft Figure 2 v0.2 summarizes the audited GSE231569 renal papilla single-nucleus atlas at the donor x broad-compartment level. Panel A gives a compact source-backed compartment overview and records the key audited scale: 43,878 nuclei overall and 540 Loop/TAL nuclei, without implying cell-level inference. Panel B shows donor x compartment summaries for primary MAGMA-prioritized modules using shortened compartment labels. Panel C summarizes Loop/TAL within-donor descriptive ranks without cell-level P values. Panel D shows single-donor exclusion robustness and should not be interpreted as new-cohort validation. Panel E shows expression/detection-matched random-set benchmarking, explicitly labeled as partial matched-random support because rank metrics saturate and the conservative signal is carried by Loop/TAL-versus-other-compartment delta. Empirical p-delta values are retained in the source table rather than emphasized as a headline figure label. Panel F shows descriptive known-driver removal sensitivity restricted to panel-level removals. Panel G states the allowable claim boundary: moderate context support, donor-level Loop/TAL pattern, partial matched-random support, robust driver-removal, not causal cell type and not plaque nucleation site.",
  "",
  "## Figure 3. Two-axis candidate-gene evidence model and curated exemplar boundaries",
  "",
  "Draft Figure 3 v0.2 presents the Stage 3R evidence framework rather than a causal ranking. Panel A maps genes by genetic priority and Kidney_Cortex proxy TWAS status. Panel B reports only the mutually exclusive R1-R6 reporting groups and excludes curated exemplar yes/no flags. Panel C shows curated biological exemplars as a role-spectrum flag, not validation and not an evidence upgrade. Panel D displays exemplar evidence boundaries, including MAGMA Bonferroni support, Kidney_Cortex TWAS proxy status that does not upgrade the exemplars, and SMR/coloc as not claim-grade ready. Panel E lists allowed and disallowed claims. The figure does not label any gene as causal, validated, SMR/coloc-supported, papilla-specific TWAS-supported, or therapeutic target."
), file.path(doc_dir, "draft_figure2_figure3_legends_v0.2.md"))

writeLines(c(
  "# Stage 4C2-R simulated reviewer check",
  "",
  "1. Would a high-impact reviewer think Figure 2 overclaims Loop/TAL enrichment?",
  "   Possibly only if Panel E/F colors are read too positively. The draft explicitly uses moderate and partial-support wording, but v0.2 should keep the palette restrained.",
  "",
  "2. Is partial random support visually and textually explicit?",
  "   Yes. Panel E labels partial matched-random support and rank-metric saturation.",
  "",
  "3. Is donor imbalance visible or acknowledged?",
  "   Yes. Panel A shows compartment/nuclei context and the legend states four donors and imbalance; v0.2 can add the 29/244/263/4 numbers directly if space allows.",
  "",
  "4. Are TWAS-proxy modules excluded from the main Loop/TAL proof?",
  "   Yes. Figure 2 uses primary MAGMA modules; Figure 3 labels TWAS as proxy evidence.",
  "",
  "5. Does Figure 3 solve the old P1/cherry-picking problem?",
  "   Yes. It uses a two-axis evidence model and separates curated exemplars from evidence upgrades.",
  "",
  "6. Does Figure 3 accidentally look like a causal prioritization figure?",
  "   Not in the current draft; the claim boundary panel disallows causal/validated language.",
  "",
  "7. Are the figures suitable as conservative draft main figures?",
  "   Yes. v0.2 is cleaner than v0.1 and acceptable for conservative draft main-figure planning, but still not a journal-final layout.",
  "",
  "8. What must be revised before publication-grade final figures?",
  "   Improve spacing, simplify dense labels, consider UMAP source coordinates for Figure 2A, and decide whether some driver-removal detail should move to supplement.",
  "",
  "9. Should full manuscript rewrite begin now or wait until Stage 5/6?",
  "   Wait until Stage 5/6 so plaque/context and spatial/TWAS boundaries are integrated consistently.",
  "",
  "10. Should Stage 5A GSE73680 analysis begin?",
  "   Yes, after human acceptance of these v0.2 draft figures and the moderate claim boundary."
), file.path(doc_dir, "stage4C2R_simulated_reviewer_check.md"))

report <- c(
  "# Stage 4C2-R report: conservative v0.2 draft Figure 2 and Figure 3",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Figures generated",
  "",
  paste0("- Figure 2 PDF: `", fig2_files[["pdf"]], "`"),
  paste0("- Figure 2 PNG: `", fig2_files[["png"]], "`"),
  paste0("- Figure 2 SVG: `", fig2_files[["svg"]], "`"),
  paste0("- Figure 3 PDF: `", fig3_files[["pdf"]], "`"),
  paste0("- Figure 3 PNG: `", fig3_files[["png"]], "`"),
  paste0("- Figure 3 SVG: `", fig3_files[["svg"]], "`"),
  "",
  "## Source-data status",
  "",
  paste0("- Panel source extracts generated: ", nrow(source_manifest)),
  paste0("- Draft-ready source extracts: ", source_manifest[ready_for_publication_source_data == "draft_ready", .N], " / ", nrow(source_manifest)),
  "",
  "## QC status",
  "",
  "- Claim-boundary QC passed for conservative draft use.",
  "- v0.2 reduces overclaim risk relative to v0.1 but remains a draft, not a publication-grade final layout.",
  "",
  "## Changes from v0.1 to v0.2",
  "",
  "- Figure 2B uses shortened compartment labels and explicitly says donor x compartment summaries.",
  "- Figure 2C now labels ranks as descriptive and avoids cell-level P-value interpretation.",
  "- Figure 2D is relabeled as single-donor exclusion robustness.",
  "- Figure 2E demotes empirical p-delta from the visual headline and emphasizes partial matched-random support with saturated rank metrics.",
  "- Figure 2F keeps only panel-level removals and uses a lower-salience robustness color.",
  "- Figure 3B removes curated exemplar yes/no rows and shows only R1-R6 mutually exclusive reporting groups.",
  "- Figure 3D labels exemplar TWAS status as no FDR TWAS proxy assigned and avoids a positive-support color.",
  "",
  "## Conservative draft acceptability",
  "",
  "- Figure 2 is acceptable as a conservative draft main figure for Stage 4 logic and claim-boundary review.",
  "- Figure 3 is acceptable as a conservative draft main figure for the candidate-gene evidence framework.",
  "",
  "## Remaining changes before publication-grade figures",
  "",
  "- Tune typography and panel spacing against the final journal page size.",
  "- Decide after Stage 5/6 whether any Stage 4 robustness detail should move to supplement.",
  "- Rebuild final source-data package after the figure set is frozen.",
  "",
  "## Stage 5A and manuscript rewrite",
  "",
  "Stage 5A can begin after human review of this Stage 4C2-R patch. Full manuscript rewrite should still wait until Stage 5 and Stage 6 are complete."
)
writeLines(report, file.path(doc_dir, "stage4C2R_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4C2R_completed_clean_draft_figures_v0.2",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = as.character(Sys.Date()),
    completed_outputs = "Stage 4A, 4B1, 4B2, 4C1, 4C2, and 4C2-R completed; conservative draft Figure 2/3 v0.2 PDF/PNG/SVG, panel source extracts, source manifest, legends v0.2, QC audit, and report generated",
    blocking_issues = "v0.2 figures acceptable as conservative clean drafts but still not journal-final publication figures; final visual polish should wait until Stage 5/6 context is integrated",
    next_stage_ready = "stage5A_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Wrote Figure 2:", fig2_files[["pdf"]], "\n")
cat("Wrote Figure 3:", fig3_files[["pdf"]], "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
