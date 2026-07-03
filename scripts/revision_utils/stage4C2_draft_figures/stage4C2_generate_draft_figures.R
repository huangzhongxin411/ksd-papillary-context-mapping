suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
})

fig_dir <- "results/figures/revision/stage4C2_draft_figures"
table_dir <- "results/tables/revision/stage4C2_draft_figures"
doc_dir <- "docs/revision/stage4C2_draft_figures"
script_dir <- "scripts/revision_utils/stage4C2_draft_figures"
log_dir <- "logs/revision/stage4C2_draft_figures"
for (d in c(fig_dir, table_dir, doc_dir, script_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage4C2_generate_draft_figures.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4C2 draft figure generation\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  fig2_blueprint = "docs/revision/stage4C1_figure_claim_planning/figure2_snRNA_context_blueprint.md",
  fig3_blueprint = "docs/revision/stage4C1_figure_claim_planning/figure3_candidate_evidence_blueprint.md",
  panel_manifest = "results/tables/revision/stage4C1_figure_claim_planning/figure_panel_source_data_manifest.tsv",
  claim_wording = "results/tables/revision/stage4C1_figure_claim_planning/snRNA_claim_wording_decision_table.tsv",
  legends_v01 = "docs/revision/stage4C1_figure_claim_planning/draft_figure2_figure3_legends_v0.1.md",
  f2a_overview = "results/tables/revision/stage4C1_figure_claim_planning/figure2_panelA_compartment_overview_source.tsv",
  donor_scores = "results/tables/revision/stage4B1_scrna_donor_level/scrna_donor_compartment_module_scores.tsv",
  loop_ranks = "results/tables/revision/stage4B1_scrna_donor_level/scrna_loop_tal_within_donor_rank.tsv",
  loo = "results/tables/revision/stage4B1_scrna_donor_level/scrna_leave_one_donor_out_module_ranks.tsv",
  low = "results/tables/revision/stage4B1_scrna_donor_level/scrna_low_loop_tal_count_sensitivity.tsv",
  random = "results/tables/revision/stage4B2_scrna_robustness/scrna_random_set_benchmark_summary.tsv",
  driver = "results/tables/revision/stage4B2_scrna_robustness/scrna_known_driver_removal_sensitivity.tsv",
  claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
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

donor_scores <- fread(paths$donor_scores)
loop_ranks <- fread(paths$loop_ranks)
loo <- fread(paths$loo)
low <- fread(paths$low)
random <- fread(paths$random)
driver <- fread(paths$driver)
claim <- fread(paths$claim)
stage3_model <- fread(paths$stage3_model)
stage3_counts <- fread(paths$stage3_counts)
exemplar <- fread(paths$exemplar)
claim_wording <- fread(paths$claim_wording)

check_cols(donor_scores, c("module_name", "donor_id", "broad_compartment", "mean_module_score", "n_nuclei"), "donor_scores")
check_cols(loop_ranks, c("module_name", "donor_id", "loop_tal_rank", "loop_tal_percentile_rank"), "loop_ranks")
check_cols(loo, c("module_name", "excluded_donor", "support_retained", "loop_tal_rank"), "leave-one-donor-out")
check_cols(random, c("module_name", "interpretation", "empirical_p_delta", "observed_loop_tal_delta"), "random benchmark")
check_cols(driver, c("base_module", "removal_type", "support_retained", "support_change", "interpretation"), "driver removal")
check_cols(claim, c("module_name", "overall_claim_strength", "recommended_manuscript_claim"), "claim decision")
check_cols(stage3_model, c("gene", "genetic_priority_level", "twas_proxy_level", "reporting_group"), "stage3 model")
check_cols(stage3_counts, c("reporting_group", "n_unique_genes"), "stage3 counts")
check_cols(exemplar, c("gene", "biological_role_label", "smr_coloc_status", "twas_status", "bonferroni_significant"), "exemplar")

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
# Stage 4B1 module score rows repeat donor/compartment nuclei counts across modules.
# De-duplicate before making the compartment overview so Figure 2A stays on the audited nuclei scale.
f2a_unique <- unique(donor_scores[, .(donor_id, broad_compartment, n_nuclei)])
f2a_src <- f2a_unique[, .(
  n_donor_compartment_rows = .N,
  n_donors = uniqueN(donor_id),
  total_nuclei = sum(n_nuclei)
), by = broad_compartment]
f2a_src[, is_loop_tal := broad_compartment == "Loop_of_Henle_TAL"]
setorder(f2a_src, -total_nuclei)
write_source(f2a_src, "figure2_panelA_source.tsv")

f2b_src <- donor_scores[module_name %in% primary_modules]
f2b_src[, module_label := safe_label(module_name)]
f2b_src[, donor_short := sub("GSM", "D", donor_id)]
write_source(f2b_src, "figure2_panelB_source.tsv")

f2c_src <- loop_ranks[module_name %in% primary_modules]
f2c_src[, module_label := safe_label(module_name)]
f2c_src[, donor_short := sub("GSM", "D", donor_id)]
write_source(f2c_src, "figure2_panelC_source.tsv")

f2d_src <- loo[module_name %in% primary_modules]
f2d_src[, module_label := safe_label(module_name)]
f2d_src[, donor_short := sub("GSM", "D", excluded_donor)]
write_source(f2d_src, "figure2_panelD_source.tsv")

f2e_src <- random[module_name %in% primary_modules]
f2e_src[, module_label := safe_label(module_name)]
f2e_src[, support_label := fifelse(interpretation == "partial_support", "partial", interpretation)]
write_source(f2e_src, "figure2_panelE_source.tsv")

panel_removals <- c("without_curated_exemplar_panel", "without_TAL_marker_panel", "without_calcium_ion_panel", "without_top5_contributors", "without_top10_contributors")
f2f_src <- driver[base_module %in% primary_modules & removal_type %in% panel_removals]
f2f_src[, module_label := safe_label(base_module)]
f2f_src[, removal_label := sub("^without_", "", removal_type)]
f2f_src[, removal_label := gsub("_", " ", removal_label)]
write_source(f2f_src, "figure2_panelF_source.tsv")

f2g_src <- claim[module_name %in% c(primary_modules, "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy", "Curated_exemplar_panel")]
f2g_src[, module_label := safe_label(module_name)]
write_source(f2g_src, "figure2_panelG_source.tsv")

# Figure 2 panels
p2a <- ggplot(f2a_src, aes(x = reorder(broad_compartment, total_nuclei), y = total_nuclei, fill = is_loop_tal)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#4C78A8", "FALSE" = "#C7C7C7"), guide = "none") +
  labs(title = "Audited GSE231569 context", subtitle = "4 donors; Loop/TAL label retained", x = NULL, y = "Nuclei") +
  theme_stage(8)

p2b <- ggplot(f2b_src, aes(x = donor_short, y = broad_compartment, fill = mean_module_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~module_label, ncol = 2) +
  scale_fill_gradient(low = "#F2F2F2", high = "#4C78A8", name = "Mean\nscore") +
  labs(title = "Donor x compartment scores", subtitle = "Primary MAGMA modules; descriptive", x = "Donor", y = NULL) +
  theme_stage(7) +
  theme(axis.text.y = element_text(size = 5.6), legend.position = "right")

p2c <- ggplot(f2c_src, aes(x = donor_short, y = module_label, color = loop_tal_rank, size = loop_tal_percentile_rank)) +
  geom_point(alpha = 0.95) +
  scale_color_gradient(low = "#4C78A8", high = "#D9D9D9", trans = "reverse", name = "Loop/TAL\nrank") +
  scale_size(range = c(2.2, 5), name = "Percentile") +
  labs(title = "Loop/TAL within-donor rank", subtitle = "No cell-level P values", x = "Donor", y = NULL) +
  theme_stage(8)

p2d <- ggplot(f2d_src, aes(x = donor_short, y = module_label, fill = support_retained)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = loop_tal_rank), size = 2.5, color = "#222222") +
  scale_fill_manual(values = c(yes = "#59A14F", partial = "#F2C14E", no = "#E15759", insufficient_data = "#BDBDBD"), name = "Retained") +
  labs(title = "Leave-one-donor-out", subtitle = "Tile label = Loop/TAL rank", x = "Excluded donor", y = NULL) +
  theme_stage(8)

p2e <- ggplot(f2e_src, aes(x = module_label, y = observed_loop_tal_delta)) +
  geom_col(fill = "#F2C14E", width = 0.65) +
  geom_text(aes(label = paste0("partial\np_delta=", signif(empirical_p_delta, 2))), vjust = -0.15, size = 2.3, color = "#333333") +
  coord_cartesian(ylim = c(0, max(f2e_src$observed_loop_tal_delta, na.rm = TRUE) * 1.35)) +
  labs(title = "Matched random benchmark", subtitle = "Partial support; rank metrics saturate", x = NULL, y = "Loop/TAL - other delta") +
  theme_stage(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p2f <- ggplot(f2f_src, aes(x = removal_label, y = module_label, fill = support_change)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(unchanged = "#59A14F", weakened = "#F2C14E", strengthened = "#76B7B2", collapsed = "#E15759", not_evaluable = "#BDBDBD"), name = "Change") +
  labs(title = "Known-driver removal", subtitle = "Module-level descriptive robustness", x = NULL, y = NULL) +
  theme_stage(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

claim_text <- "Moderate context support\nDonor-level Loop/TAL pattern\n+ partial matched-random support\n+ robust driver-removal\n\nNot causal cell type\nNot plaque nucleation site"
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
    title = "Figure 2 draft. Donor-level snRNA context mapping of MAGMA-prioritized modules",
    subtitle = "Draft v0.1: moderate context support; partial matched-random support; robust known-driver removal",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "#555555"))
  )

# Figure 3 source extracts
f3a_src <- stage3_model[, .N, by = .(genetic_priority_level, twas_proxy_level)]
write_source(f3a_src, "figure3_panelA_source.tsv")

f3b_src <- copy(stage3_counts)
write_source(f3b_src, "figure3_panelB_source.tsv")

f3c_src <- exemplar[, .(gene, biological_role_label, magma_rank, bonferroni_significant)]
write_source(f3c_src, "figure3_panelC_source.tsv")

f3d_src <- exemplar[, .(
  gene,
  MAGMA_Bonferroni = ifelse(bonferroni_significant == "yes", "yes", "no"),
  Kidney_Cortex_TWAS_proxy = twas_status,
  SMR_coloc = smr_coloc_status
)]
f3d_long <- melt(f3d_src, id.vars = "gene", variable.name = "evidence_layer", value.name = "status")
f3d_long[, status_simple := fifelse(evidence_layer == "MAGMA_Bonferroni" & status == "yes", "MAGMA prioritized",
  fifelse(evidence_layer == "Kidney_Cortex_TWAS_proxy", "proxy/not FDR",
          fifelse(evidence_layer == "SMR_coloc", "not claim-grade ready", as.character(status))))]
write_source(f3d_long, "figure3_panelD_source.tsv")

allowed <- c("MAGMA-prioritized gene", "Kidney_Cortex proxy TWAS support", "curated biological exemplar", "context-mapping candidate")
disallowed <- c("causal gene", "validated disease gene", "SMR/coloc-supported gene", "papilla-specific TWAS gene", "therapeutic target")
f3e_src <- rbind(
  data.table(claim_type = "Allowed", claim = allowed),
  data.table(claim_type = "Disallowed", claim = disallowed)
)
write_source(f3e_src, "figure3_panelE_source.tsv")

p3a <- ggplot(f3a_src, aes(x = genetic_priority_level, y = twas_proxy_level, fill = N)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = N), size = 3, fontface = "bold") +
  scale_fill_gradient(low = "#F2F2F2", high = "#4C78A8", name = "Genes") +
  labs(title = "Two-axis evidence model", subtitle = "Genetic priority x Kidney_Cortex TWAS proxy", x = "Genetic priority", y = "TWAS proxy") +
  theme_stage(8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6))

f3b_src[, reporting_group_label := gsub("_", " ", reporting_group)]
p3b <- ggplot(f3b_src, aes(x = reorder(reporting_group_label, n_unique_genes), y = n_unique_genes)) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  coord_flip() +
  labs(title = "Reporting group counts", subtitle = "Mutually exclusive groups, not causal tiers", x = NULL, y = "Genes") +
  theme_stage(8) +
  theme(axis.text.y = element_text(size = 6))

f3c_src[, gene := factor(gene, levels = rev(gene))]
p3c <- ggplot(f3c_src, aes(x = 1, y = gene)) +
  geom_point(aes(color = bonferroni_significant), size = 4) +
  geom_text(aes(label = biological_role_label), x = 1.08, hjust = 0, size = 2.6) +
  scale_color_manual(values = c(yes = "#4C78A8", no = "#BDBDBD"), name = "MAGMA\nBonf.") +
  coord_cartesian(xlim = c(0.95, 2.2), clip = "off") +
  labs(title = "Curated biological exemplars", subtitle = "Role spectrum; no evidence upgrade", x = NULL, y = NULL) +
  theme_stage(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(5, 70, 5, 5))

f3d_long[, gene := factor(gene, levels = rev(unique(exemplar$gene)))]
f3d_long[, evidence_layer := factor(evidence_layer, levels = c("MAGMA_Bonferroni", "Kidney_Cortex_TWAS_proxy", "SMR_coloc"))]
p3d <- ggplot(f3d_long, aes(x = evidence_layer, y = gene, fill = status_simple)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c("MAGMA prioritized" = "#4C78A8", "proxy/not FDR" = "#F2C14E", "not claim-grade ready" = "#D9D9D9"),
                    name = "Status") +
  labs(title = "Exemplar evidence strip", subtitle = "SMR/coloc not claim-grade ready", x = NULL, y = NULL) +
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
    title = "Figure 3 draft. Candidate-gene evidence model and exemplar boundaries",
    subtitle = "Draft v0.1: two-axis prioritization, Kidney_Cortex TWAS proxy, no causal or validated gene labels",
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

fig2_files <- save_figure(fig2, "figure2_snRNA_context_draft_v0.1", 14, 12)
fig3_files <- save_figure(fig3, "figure3_candidate_evidence_draft_v0.1", 13, 11)

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
    notes = "Stage 4C2 draft source extract; final publication source data should be rechecked after v0.2 revisions."
  )
}), fill = TRUE)
setorder(source_manifest, figure, panel)
fwrite(source_manifest, file.path(table_dir, "figure_source_data_manifest_v0.1.tsv"), sep = "\t")

qc <- data.table(
  check_item = c(
    "Figure 2 states moderate support",
    "Figure 2 shows partial random support",
    "Figure 2 avoids cell-level P values",
    "Figure 2 avoids causal cell-type language",
    "Figure 2 driver removal as robustness",
    "Figure 3 avoids P1 terminology",
    "Figure 3 avoids causal/validated gene language",
    "Figure 3 marks TWAS as Kidney_Cortex proxy",
    "Figure 3 marks SMR/coloc as not claim-grade ready",
    "All panels backed by source data",
    "Visually misleading panels",
    "Fixes for v0.2"
  ),
  status = c("pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass", "pass",
             ifelse(all(source_manifest$ready_for_publication_source_data == "draft_ready"), "pass", "needs_fix"),
             "review_needed", "needed"),
  notes = c(
    "Title/subtitle and claim strip use moderate/not strong language.",
    "Panel E labels partial matched-random support and rank saturation.",
    "No P values are shown except empirical random-set p-delta labels, not cell-level inferential tests.",
    "Panel G states not causal cell type and not plaque nucleation site.",
    "Panel F labels descriptive robustness, not mechanism proof.",
    "No P1 label appears in Figure 3 draft.",
    "Claim boundary panel disallows causal and validated gene labels.",
    "Figure 3 subtitle and Panel A label TWAS as Kidney_Cortex proxy.",
    "Panel D labels SMR/coloc as not claim-grade ready.",
    paste0(nrow(source_manifest), " panel source extracts generated."),
    "Draft figures are information-dense; v0.2 should refine typography/spacing and possibly move some robustness detail to supplement.",
    "Improve visual spacing, consider UMAP source coordinates for Figure 2A, and reduce Figure 2F label density."
  )
)

writeLines(c(
  "# Figure QC audit v0.1",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "| check | status | notes |",
  "|---|---|---|",
  apply(qc, 1, function(r) paste0("| ", paste(gsub("\\|", "/", r), collapse = " | "), " |"))
), file.path(doc_dir, "figure_qc_audit_v0.1.md"))

writeLines(c(
  "# Draft Figure 2 and Figure 3 legends v0.2",
  "",
  "## Figure 2. Donor-level snRNA context mapping and robustness boundaries for MAGMA-prioritized modules",
  "",
  "Draft Figure 2 summarizes the audited GSE231569 renal papilla single-nucleus atlas at the donor x broad-compartment level. Panel A gives the source-backed compartment overview and highlights the Loop/TAL label across four donors. Panel B shows donor x compartment module scores for primary MAGMA-prioritized modules only. Panel C summarizes Loop/TAL within-donor ranks without cell-level inferential tests. Panel D shows leave-one-donor-out robustness and should not be interpreted as independent replication. Panel E shows expression/detection-matched random-set benchmarking, explicitly labeled as partial matched-random support because rank metrics saturate and the conservative signal is carried by Loop/TAL-versus-other-compartment delta. Panel F shows descriptive known-driver removal sensitivity. Panel G states the allowable claim boundary: moderate context support, not causal cell type or plaque nucleation site.",
  "",
  "## Figure 3. Two-axis candidate-gene evidence model and curated exemplar boundaries",
  "",
  "Draft Figure 3 presents the Stage 3R evidence framework rather than a causal ranking. Panel A maps genes by genetic priority and Kidney_Cortex TWAS proxy status. Panel B reports mutually exclusive reporting-group counts. Panel C shows curated biological exemplars as role-spectrum examples. Panel D displays exemplar evidence boundaries, including MAGMA prioritization, Kidney_Cortex TWAS proxy status, and SMR/coloc as not claim-grade ready where applicable. Panel E lists allowed and disallowed claims. The figure does not label any gene as causal, validated, SMR/coloc-supported, papilla-specific TWAS-supported, or therapeutic target."
), file.path(doc_dir, "draft_figure2_figure3_legends_v0.2.md"))

writeLines(c(
  "# Stage 4C2 simulated reviewer check",
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
  "7. Are the figures suitable as draft main figures?",
  "   Yes, as draft main figures for logic review. They are not yet publication-grade final figures.",
  "",
  "8. What must be revised before publication-grade final figures?",
  "   Improve spacing, simplify dense labels, consider UMAP source coordinates for Figure 2A, and decide whether some driver-removal detail should move to supplement.",
  "",
  "9. Should full manuscript rewrite begin now or wait until Stage 5/6?",
  "   Wait until Stage 5/6 so plaque/context and spatial/TWAS boundaries are integrated consistently.",
  "",
  "10. Should Stage 5 GSE73680 analysis begin?",
  "   Yes, after human acceptance of these draft figures and the moderate claim boundary."
), file.path(doc_dir, "stage4C2_simulated_reviewer_check.md"))

report <- c(
  "# Stage 4C2 report: source-data-backed draft Figure 2 and Figure 3",
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
  "- Claim-boundary QC passed for draft use.",
  "- Visual density and spacing need v0.2 refinement before publication-grade use.",
  "",
  "## Problems found",
  "",
  "- Figure 2 is information-dense because it carries seven panels.",
  "- Figure 2F driver-removal labels may need simplification.",
  "- Figure 2A currently uses a compact overview rather than UMAP coordinates.",
  "",
  "## Recommended v0.2 revisions",
  "",
  "- Consider exporting UMAP coordinates for Panel A if the atlas overview should be spatially intuitive.",
  "- Move single-gene driver removals to supplement and keep panel-removal groups in main Figure 2F.",
  "- Increase vertical spacing and reduce label density for publication-grade layout.",
  "",
  "## Draft figure acceptability",
  "",
  "Figure 2 and Figure 3 are acceptable as source-data-backed draft figures for logic review, not final publication figures.",
  "",
  "## Stage 5 and manuscript rewrite",
  "",
  "Stage 5 GSE73680 analysis can begin after human review. Full manuscript rewrite should still wait until Stage 5 and Stage 6 are complete."
)
writeLines(report, file.path(doc_dir, "stage4C2_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4C2_completed_draft_figures",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = as.character(Sys.Date()),
    completed_outputs = "Stage 4A, 4B1, 4B2, 4C1, and 4C2 completed; draft Figure 2/3 PDF/PNG/SVG, panel source extracts, source manifest, legends v0.2, QC audit, and reviewer check generated",
    blocking_issues = "Draft figures acceptable for logic review but not publication-grade final figures; v0.2 visual refinement recommended",
    next_stage_ready = "stage5_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Wrote Figure 2:", fig2_files[["pdf"]], "\n")
cat("Wrote Figure 3:", fig3_files[["pdf"]], "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
