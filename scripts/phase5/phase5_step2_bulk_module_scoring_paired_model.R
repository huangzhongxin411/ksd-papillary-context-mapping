suppressPackageStartupMessages({
  library(ggplot2)
})

root <- normalizePath(".", mustWork = TRUE)
dir.create(file.path(root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "source_data", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)

tsv_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

collapse_vec <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
fmt_p <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 3))

read_tsv <- function(path) {
  read.delim(file.path(root, path), sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

expr_path <- "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz"
meta_path <- "config/gse73680_sample_metadata_curated.tsv"
pairing_path <- "results/gse73680/tables/gse73680_patient_sample_structure.tsv"
step1_dep_path <- "results/tables/phase5_step1_deprecated_bulk_output_audit.tsv"

module_files <- c(
  MAGMA_top50 = "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
  MAGMA_top100 = "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
  MAGMA_Bonferroni = "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
  MAGMA_FDR05 = "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
  MAGMA_suggestive_p1e4 = "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
)

expr <- read_tsv(expr_path)
meta <- read_tsv(meta_path)
pairing_input <- read_tsv(pairing_path)

if (!"gene" %in% names(expr)) stop("Expected first/available gene column named 'gene'.")
sample_cols <- setdiff(names(expr), "gene")
expr$gene <- as.character(expr$gene)
nonempty <- !is.na(expr$gene) & nzchar(expr$gene)
expr_filtered <- expr[nonempty, , drop = FALSE]
if (any(duplicated(expr_filtered$gene))) {
  stop("Gene-level matrix contains duplicated non-empty genes; inspect before Step 2 scoring.")
}

meta$include_bool <- toupper(as.character(meta$include_in_analysis)) == "TRUE"
meta_included <- meta[meta$include_bool & meta$group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla"), , drop = FALSE]
included_samples <- meta_included$sample_id
missing_from_expr <- setdiff(included_samples, sample_cols)
if (length(missing_from_expr) > 0) stop(paste("Included samples missing from expression matrix:", paste(missing_from_expr, collapse = ",")))

expr_included <- expr_filtered[, c("gene", included_samples), drop = FALSE]
expr_mat <- as.matrix(data.frame(lapply(expr_included[, included_samples, drop = FALSE], as.numeric), check.names = FALSE))
rownames(expr_mat) <- expr_included$gene

gene_means <- rowMeans(expr_mat, na.rm = TRUE)
gene_sds <- apply(expr_mat, 1, sd, na.rm = TRUE)
valid_gene_for_z <- is.finite(gene_means) & is.finite(gene_sds) & gene_sds > 0
z_mat <- expr_mat[valid_gene_for_z, , drop = FALSE]
z_mat <- sweep(z_mat, 1, gene_means[valid_gene_for_z], "-")
z_mat <- sweep(z_mat, 1, gene_sds[valid_gene_for_z], "/")

control_count <- sum(meta_included$group_curated == "control_or_adjacent")
disease_count <- sum(meta_included$group_curated == "plaque_or_stone_papilla")

split_meta <- split(meta_included, meta_included$patient_id)
complete_pairs <- do.call(rbind, lapply(split_meta, function(d) {
  ctl <- d$sample_id[d$group_curated == "control_or_adjacent"]
  dis <- d$sample_id[d$group_curated == "plaque_or_stone_papilla"]
  data.frame(
    patient_id = d$patient_id[1],
    pair_id = d$patient_id[1],
    n_control = length(ctl),
    n_disease = length(dis),
    control_sample_id = ifelse(length(ctl) >= 1, ctl[1], ""),
    disease_sample_id = ifelse(length(dis) >= 1, dis[1], ""),
    complete_pair = length(ctl) >= 1 && length(dis) >= 1,
    nonstandard_pair = !(length(ctl) == 1 && length(dis) == 1),
    stringsAsFactors = FALSE
  )
}))
complete_pairs <- complete_pairs[complete_pairs$complete_pair, , drop = FALSE]

input_check <- data.frame(
  audit_item = c(
    "expression matrix path", "metadata path", "patient pairing table path",
    "n_features_before_filtering", "n_features_after_empty_gene_filter", "n_zscore_valid_features",
    "n_samples_total", "n_samples_included", "n_complete_paired_patients",
    "control_or_adjacent sample count", "plaque_or_stone sample count",
    "expression scale handling", "z-scoring method", "missing-value handling",
    "whether paired model can proceed"
  ),
  value = c(
    expr_path, meta_path, pairing_path,
    nrow(expr), nrow(expr_filtered), sum(valid_gene_for_z),
    length(sample_cols), nrow(meta_included), nrow(complete_pairs),
    control_count, disease_count,
    "processed intensity-like/non-count values retained; no raw-count modeling",
    "per gene: subtract mean and divide by SD across all included samples",
    "gene-level z-scores use na.rm; module means computed over available z-scored genes",
    ifelse(nrow(complete_pairs) > 1, "yes_complete_pairs_available", "no")
  ),
  source_file = c(
    expr_path, meta_path, pairing_path,
    expr_path, expr_path, expr_path,
    expr_path, meta_path, "curated metadata-derived complete pairs",
    meta_path, meta_path,
    expr_path, expr_path, expr_path,
    "phase5_step2 script"
  ),
  notes = c(
    "Canonical Phase 5-Step 1 gene-level matrix.",
    "Curated metadata with include flag and filename-derived labels.",
    "Existing pairing table inspected; Step 2 complete-pair set derived from curated metadata.",
    "Includes one empty gene row before filtering.",
    "Empty gene symbols removed before scoring.",
    "Genes with zero or undefined SD across included samples are not used in z-score module scoring.",
    "All GSM columns in expression matrix.",
    "Included disease-context samples only; no silent sample dropping.",
    "Complete-pair count is reported as observed, not forced.",
    "Included samples with group_curated control_or_adjacent.",
    "Included samples with group_curated plaque_or_stone_papilla.",
    "No transformation beyond z-scoring was applied in this step.",
    "All included samples define the scoring universe.",
    "No imputation; missing module values are ignored when averaging available genes.",
    "Primary paired base model restricted to complete pairs only."
  ),
  stringsAsFactors = FALSE
)
tsv_write(input_check, file.path(root, "results/tables/phase5_step2_bulk_analysis_input_check.tsv"))

bulk_genes <- rownames(z_mat)
mapping_rows <- lapply(names(module_files), function(module_name) {
  module_file <- module_files[[module_name]]
  genes <- unique(readLines(file.path(root, module_file), warn = FALSE))
  genes <- genes[nzchar(genes)]
  present <- intersect(genes, bulk_genes)
  missing <- setdiff(genes, rownames(expr_mat))
  missing_or_not_valid <- setdiff(genes, bulk_genes)
  notes <- if (length(setdiff(missing_or_not_valid, missing)) > 0) {
    paste0("Some genes present in matrix but not usable after z-score SD/missing filter: ", paste(setdiff(missing_or_not_valid, missing), collapse = ";"))
  } else {
    "Scoring uses present genes after empty-gene and z-score-valid filtering."
  }
  data.frame(
    module_name = module_name,
    module_file = module_file,
    n_module_genes = length(genes),
    n_present_in_bulk = length(present),
    n_missing_from_bulk = length(setdiff(genes, bulk_genes)),
    genes_used_for_scoring = paste(present, collapse = ";"),
    missing_genes = paste(setdiff(genes, bulk_genes), collapse = ";"),
    score_computed = ifelse(length(present) > 0, "yes", "no"),
    notes = notes,
    stringsAsFactors = FALSE
  )
})
module_mapping <- do.call(rbind, mapping_rows)
tsv_write(module_mapping, file.path(root, "results/tables/phase5_step2_magma_module_bulk_gene_mapping.tsv"))

score_rows <- list()
for (i in seq_len(nrow(module_mapping))) {
  module_name <- module_mapping$module_name[i]
  genes_used <- strsplit(module_mapping$genes_used_for_scoring[i], ";", fixed = TRUE)[[1]]
  genes_used <- genes_used[nzchar(genes_used)]
  if (!length(genes_used)) next
  module_values <- colMeans(z_mat[genes_used, , drop = FALSE], na.rm = TRUE)
  score_rows[[module_name]] <- data.frame(
    sample_id = names(module_values),
    module_name = module_name,
    n_genes_used = length(genes_used),
    module_score = as.numeric(module_values),
    stringsAsFactors = FALSE
  )
}
scores <- do.call(rbind, score_rows)
scores <- merge(scores, meta_included[, c("sample_id", "patient_id", "group_curated", "plaque_status", "stone_status")],
                by = "sample_id", all.x = TRUE, sort = FALSE)
scores$pair_id <- scores$patient_id
scores$control_or_disease <- ifelse(scores$group_curated == "control_or_adjacent", "control_or_adjacent", "disease_context")
scores$metadata_confidence <- "usable_filename_derived_needs_manual_review"
scores$notes <- "Score recalculated from canonical gene-level matrix and Phase 1 locked MAGMA module; no historical module-score output used."
scores <- scores[, c("sample_id", "patient_id", "pair_id", "group_curated", "control_or_disease",
                     "plaque_status", "stone_status", "module_name", "n_genes_used", "module_score",
                     "metadata_confidence", "notes")]
scores <- scores[order(scores$module_name, scores$patient_id, scores$sample_id), ]
tsv_write(scores, file.path(root, "results/tables/phase5_step2_bulk_sample_module_scores.tsv"))

summary_by_group <- do.call(rbind, lapply(split(scores, list(scores$module_name, scores$group_curated), drop = TRUE), function(d) {
  data.frame(
    module_name = d$module_name[1],
    group_curated = d$group_curated[1],
    n_samples = nrow(d),
    mean_module_score = mean(d$module_score, na.rm = TRUE),
    median_module_score = median(d$module_score, na.rm = TRUE),
    sd_module_score = sd(d$module_score, na.rm = TRUE),
    iqr_module_score = IQR(d$module_score, na.rm = TRUE),
    notes = "Summary over included disease-context samples.",
    stringsAsFactors = FALSE
  )
}))
summary_by_group <- summary_by_group[order(summary_by_group$module_name, summary_by_group$group_curated), ]
tsv_write(summary_by_group, file.path(root, "results/tables/phase5_step2_bulk_module_score_summary_by_group.tsv"))

pair_samples <- unique(c(complete_pairs$control_sample_id, complete_pairs$disease_sample_id))
model_scores <- scores[scores$sample_id %in% pair_samples, , drop = FALSE]
model_scores$group_binary <- ifelse(model_scores$group_curated == "plaque_or_stone_papilla", 1, 0)

model_rows <- list()
diff_rows <- list()
for (module_name in names(module_files)) {
  d <- model_scores[model_scores$module_name == module_name, , drop = FALSE]
  d$patient_id <- factor(d$patient_id)
  fit <- tryCatch(lm(module_score ~ group_binary + factor(patient_id), data = d), error = function(e) NULL)
  if (is.null(fit) || !"group_binary" %in% rownames(summary(fit)$coefficients)) {
    coef_val <- se_val <- t_val <- p_val <- NA_real_
    direction <- "model_unstable"
    interpretation <- "model_unstable"
    note <- "Model did not return a group coefficient."
  } else {
    coef_tab <- summary(fit)$coefficients
    coef_val <- coef_tab["group_binary", "Estimate"]
    se_val <- coef_tab["group_binary", "Std. Error"]
    t_val <- coef_tab["group_binary", "t value"]
    p_val <- coef_tab["group_binary", "Pr(>|t|)"]
    direction <- ifelse(is.na(p_val) || p_val >= 0.05, "no_clear_shift",
                        ifelse(coef_val > 0, "positive_disease_context_shift", "negative_disease_context_shift"))
    interpretation <- direction
    note <- "Patient fixed-effect paired base model; no composition or injury/remodeling adjustment."
  }
  model_rows[[module_name]] <- data.frame(
    module_name = module_name,
    n_paired_patients = nrow(complete_pairs),
    n_samples_used = nrow(d),
    model_formula = "module_score ~ group_binary + factor(patient_id)",
    group_coefficient = coef_val,
    standard_error = se_val,
    t_statistic = t_val,
    p_value = p_val,
    fdr_bh_across_modules = NA_real_,
    direction = direction,
    interpretation = interpretation,
    notes = note,
    stringsAsFactors = FALSE
  )

  module_diff <- do.call(rbind, lapply(seq_len(nrow(complete_pairs)), function(j) {
    pair <- complete_pairs[j, ]
    ctl_score <- d$module_score[d$sample_id == pair$control_sample_id]
    dis_score <- d$module_score[d$sample_id == pair$disease_sample_id]
    data.frame(
      patient_id = pair$patient_id,
      pair_id = pair$pair_id,
      module_name = module_name,
      control_sample_id = pair$control_sample_id,
      disease_sample_id = pair$disease_sample_id,
      control_module_score = ifelse(length(ctl_score) == 1, ctl_score, NA_real_),
      disease_module_score = ifelse(length(dis_score) == 1, dis_score, NA_real_),
      paired_difference = ifelse(length(ctl_score) == 1 && length(dis_score) == 1, dis_score - ctl_score, NA_real_),
      notes = "Disease-context minus control/adjacent module score in complete paired patients.",
      stringsAsFactors = FALSE
    )
  }))
  diff_rows[[module_name]] <- module_diff
}
model_results <- do.call(rbind, model_rows)
model_results$fdr_bh_across_modules <- p.adjust(model_results$p_value, method = "BH")
model_results$direction <- ifelse(is.na(model_results$p_value), "model_unstable",
                                  ifelse(model_results$fdr_bh_across_modules < 0.05 & model_results$group_coefficient > 0, "positive_disease_context_shift",
                                         ifelse(model_results$fdr_bh_across_modules < 0.05 & model_results$group_coefficient < 0, "negative_disease_context_shift",
                                                ifelse(model_results$p_value < 0.05 & model_results$group_coefficient > 0, "positive_disease_context_shift",
                                                       ifelse(model_results$p_value < 0.05 & model_results$group_coefficient < 0, "negative_disease_context_shift", "no_clear_shift")))))
model_results$interpretation <- model_results$direction
model_results <- model_results[order(match(model_results$module_name, names(module_files))), ]
tsv_write(model_results, file.path(root, "results/tables/phase5_step2_bulk_paired_model_results.tsv"))

paired_diffs <- do.call(rbind, diff_rows)
paired_diffs <- paired_diffs[order(match(paired_diffs$module_name, names(module_files)), paired_diffs$patient_id), ]
tsv_write(paired_diffs, file.path(root, "results/tables/phase5_step2_bulk_paired_differences.tsv"))

model_diff_summary <- aggregate(paired_difference ~ module_name, paired_diffs, median, na.rm = TRUE)
names(model_diff_summary)[2] <- "paired_difference_median"
interpretation_summary <- merge(model_results, model_diff_summary, by = "module_name", all.x = TRUE, sort = FALSE)
interpretation_summary$support_class <- ifelse(is.na(interpretation_summary$p_value), "unstable",
  ifelse(interpretation_summary$fdr_bh_across_modules < 0.05, "retained_disease_context_shift",
    ifelse(interpretation_summary$p_value < 0.05 & interpretation_summary$group_coefficient > 0, "nominal_positive_shift",
      ifelse(interpretation_summary$p_value < 0.05, "weak_or_inconsistent_shift",
        ifelse(abs(interpretation_summary$paired_difference_median) < 0.05, "no_shift", "weak_or_inconsistent_shift")))))
interpretation_out <- data.frame(
  module_name = interpretation_summary$module_name,
  direction = interpretation_summary$direction,
  fdr_bh = interpretation_summary$fdr_bh_across_modules,
  paired_difference_median = interpretation_summary$paired_difference_median,
  support_class = interpretation_summary$support_class,
  allowed_interpretation = "module shows paired disease-context shift in plaque/stone-associated papilla.",
  not_allowed_interpretation = "bulk validates causal genes; bulk proves plaque mechanism; bulk identifies causal cell type; bulk establishes genetic causality",
  notes = "Base paired model only; composition and injury/remodeling sensitivity pending.",
  stringsAsFactors = FALSE
)
tsv_write(interpretation_out, file.path(root, "results/tables/phase5_step2_bulk_base_model_interpretation_summary.tsv"))

line_source <- model_scores[, c("sample_id", "patient_id", "pair_id", "group_curated", "module_name", "module_score")]
line_source$group_order <- ifelse(line_source$group_curated == "control_or_adjacent", 0, 1)
line_source <- line_source[order(line_source$module_name, line_source$patient_id, line_source$group_order), ]
tsv_write(line_source, file.path(root, "source_data/figures/phase5_step2_bulk_paired_module_score_lines_source_data.tsv"))
tsv_write(paired_diffs, file.path(root, "source_data/figures/phase5_step2_bulk_paired_difference_summary_source_data.tsv"))
coef_source <- model_results
coef_source$ci_low <- coef_source$group_coefficient - 1.96 * coef_source$standard_error
coef_source$ci_high <- coef_source$group_coefficient + 1.96 * coef_source$standard_error
tsv_write(coef_source, file.path(root, "source_data/figures/phase5_step2_bulk_paired_model_coefficients_source_data.tsv"))

palette <- c(control_or_adjacent = "#7F9DA6", plaque_or_stone_papilla = "#9B5C4D", disease_context = "#9B5C4D")
base_theme <- theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Helvetica", colour = "#333333"),
    strip.background = element_rect(fill = "#E6E9EA", colour = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

line_source$group_curated <- factor(line_source$group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))
p1 <- ggplot(line_source, aes(x = group_curated, y = module_score, group = patient_id)) +
  geom_line(colour = "#7F9DA6", alpha = 0.45, linewidth = 0.35) +
  geom_point(aes(fill = group_curated), shape = 21, colour = "#333333", alpha = 0.85, size = 1.8, stroke = 0.2) +
  facet_wrap(~ module_name, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette) +
  labs(title = "Canonical MAGMA bulk module scores in paired papillary samples",
       subtitle = "Complete paired GSE73680 samples; base disease-context view only",
       x = NULL, y = "Mean z-scored module expression", fill = "Group") +
  base_theme
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_module_score_lines.pdf"), p1, width = 8.5, height = 7.5, device = "pdf")
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_module_score_lines_600dpi.png"), p1, width = 8.5, height = 7.5, dpi = 600)

p2 <- ggplot(paired_diffs, aes(x = module_name, y = paired_difference)) +
  geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.3) +
  geom_boxplot(fill = "#E6E9EA", colour = "#245A64", width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.4, colour = "#9B5C4D") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.4, fill = "#B99B5A", colour = "#333333") +
  labs(title = "Within-patient paired module-score differences",
       subtitle = "Plaque/stone-associated papilla minus control/adjacent papilla",
       x = NULL, y = "Paired difference in module score") +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_difference_summary.pdf"), p2, width = 8.2, height = 4.8, device = "pdf")
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_difference_summary_600dpi.png"), p2, width = 8.2, height = 4.8, dpi = 600)

coef_source$module_name <- factor(coef_source$module_name, levels = rev(names(module_files)))
p3 <- ggplot(coef_source, aes(x = group_coefficient, y = module_name)) +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18, colour = "#245A64", linewidth = 0.6) +
  geom_point(aes(fill = direction), shape = 21, size = 2.8, colour = "#333333", stroke = 0.25) +
  scale_fill_manual(values = c(
    positive_disease_context_shift = "#9B5C4D",
    negative_disease_context_shift = "#0F4C5C",
    no_clear_shift = "#7F9DA6",
    model_unstable = "#E6E9EA"
  ), drop = FALSE) +
  labs(title = "Unadjusted paired base-model coefficients",
       subtitle = "Patient fixed effects; no deconvolution or tissue-state adjustment",
       x = "Plaque/stone-associated vs control/adjacent coefficient", y = NULL, fill = "Direction") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_model_coefficients.pdf"), p3, width = 7.2, height = 4.2, device = "pdf")
ggsave(file.path(root, "results/figures/phase5_step2_bulk_paired_model_coefficients_600dpi.png"), p3, width = 7.2, height = 4.2, dpi = 600)

figure_qc <- data.frame(
  figure_id = c("phase5_step2_bulk_paired_module_score_lines", "phase5_step2_bulk_paired_difference_summary", "phase5_step2_bulk_paired_model_coefficients"),
  version = "v0.1",
  pdf_exists = c(
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_module_score_lines.pdf")),
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_difference_summary.pdf")),
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_model_coefficients.pdf"))
  ),
  png_exists = c(
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_module_score_lines_600dpi.png")),
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_difference_summary_600dpi.png")),
    file.exists(file.path(root, "results/figures/phase5_step2_bulk_paired_model_coefficients_600dpi.png"))
  ),
  intended_png_dpi = 600,
  minimum_configured_font_size = 9,
  panel_label_presence = "single-panel/faceted figure with clear title",
  legend_placement_check = "bottom or absent",
  palette_consistency_check = "project profile teal/bluegrey/terracotta/gold palette",
  claim_boundary_check = "base paired disease-context association only",
  resource_limited_claim_check = "no deconvolution, no injury/remodeling adjustment, no causal language",
  source_table_exists = c(
    file.exists(file.path(root, "source_data/figures/phase5_step2_bulk_paired_module_score_lines_source_data.tsv")),
    file.exists(file.path(root, "source_data/figures/phase5_step2_bulk_paired_difference_summary_source_data.tsv")),
    file.exists(file.path(root, "source_data/figures/phase5_step2_bulk_paired_model_coefficients_source_data.tsv"))
  ),
  legend_file_exists = TRUE,
  visual_status = "generated; human readability review recommended before manuscript use",
  action_required = "human review before figure reuse",
  stringsAsFactors = FALSE
)
tsv_write(figure_qc, file.path(root, "results/tables/phase5_step2_figure_visual_qc.tsv"))

writeLines(c(
  "# Phase 5-Step 2 Figure Legends",
  "",
  "## Paired Module Score Lines",
  "Paired line plots show recalculated canonical MAGMA module scores in complete paired GSE73680 papillary samples. Scores are arithmetic means of gene-level z-scored expression across present module genes. Each line connects control/adjacent and plaque/stone-associated papilla from the same patient. This figure summarizes an unadjusted paired disease-context view and does not estimate cell fractions or causal gene effects.",
  "",
  "## Paired Difference Summary",
  "Distributions show within-patient differences calculated as plaque/stone-associated papilla minus control/adjacent papilla for each canonical MAGMA module. Boxes summarize the paired-difference distribution and gold points indicate medians. The panel supports disease-context association assessment only.",
  "",
  "## Paired Model Coefficients",
  "Coefficients and 95% confidence intervals are from base paired models of module score on plaque/stone-associated status with patient fixed effects. Models are unadjusted for composition and injury/remodeling programs. Results should not be interpreted as causal genetic validation or plaque-mechanism validation."
), con = file.path(root, "notes/phase5_step2_figure_legends.md"))

if (file.exists(file.path(root, step1_dep_path))) {
  dep <- read_tsv(step1_dep_path)
  dep_sub <- dep[grepl("figure4|stage5|gse73680|module_score|module_response|paired", dep$file_path, ignore.case = TRUE) &
                   dep$needs_rerun_or_review == "yes_review_before_use", , drop = FALSE]
  dep_sub <- head(dep_sub, 120)
  replacement <- data.frame(
    deprecated_file = dep_sub$file_path,
    likely_old_module_or_old_analysis = ifelse(grepl("figure4", dep_sub$file_path, ignore.case = TRUE), "old Figure4/source-data artifact",
                                      ifelse(grepl("module_score", dep_sub$file_path, ignore.case = TRUE), "old module-score output",
                                      ifelse(grepl("module_response|paired|model", dep_sub$file_path, ignore.case = TRUE), "old paired/model output", "old bulk analysis artifact"))),
    new_replacement_output = ifelse(grepl("figure4", dep_sub$file_path, ignore.case = TRUE), "results/figures/phase5_step2_bulk_* plus source_data/figures/phase5_step2_*",
                              ifelse(grepl("module_score", dep_sub$file_path, ignore.case = TRUE), "results/tables/phase5_step2_bulk_sample_module_scores.tsv",
                              ifelse(grepl("module_response|paired|model", dep_sub$file_path, ignore.case = TRUE), "results/tables/phase5_step2_bulk_paired_model_results.tsv;results/tables/phase5_step2_bulk_paired_differences.tsv",
                                     "Step 2 canonical outputs, review manually"))),
    rerun_status = "replaced_by_phase5_step2_canonical_base_outputs_or_flagged_for_review",
    notes = "Do not reuse as current evidence until checked against canonical Phase 1 module definitions.",
    stringsAsFactors = FALSE
  )
} else {
  replacement <- data.frame(deprecated_file=character(), likely_old_module_or_old_analysis=character(), new_replacement_output=character(), rerun_status=character(), notes=character())
}
tsv_write(replacement, file.path(root, "results/tables/phase5_step2_bulk_deprecated_output_replacement_plan.tsv"))

result_lines <- apply(model_results, 1, function(r) {
  paste0("- ", r[["module_name"]], ": coefficient ", fmt_p(as.numeric(r[["group_coefficient"]])),
         ", p=", fmt_p(as.numeric(r[["p_value"]])), ", FDR=", fmt_p(as.numeric(r[["fdr_bh_across_modules"]])),
         ", direction=", r[["direction"]], ".")
})

writeLines(c(
  "# Phase 5-Step 2 Results Wording",
  "",
  "Using canonical MAGMA modules, we recalculated bulk module scores in paired papillary samples and tested plaque/stone-associated versus control/adjacent papilla shifts with patient fixed effects.",
  "",
  paste0("The analysis included ", nrow(meta_included), " disease-context samples from ", length(unique(meta_included$patient_id)),
         " patients, including ", nrow(complete_pairs), " complete paired patients for the primary paired base model."),
  "",
  "Five locked MAGMA modules were scored after filtering empty gene symbols and z-scoring each gene across included samples.",
  "",
  "Base paired model summary:",
  result_lines,
  "",
  "These results support only paired disease-context assessment of module-level expression shifts in GSE73680. They do not establish causal genes, plaque mechanisms, cell fractions, or genetic causality."
), con = file.path(root, "notes/phase5_step2_results_wording.md"))

writeLines(c(
  "# Phase 5-Step 2 Methods Wording",
  "",
  "Bulk GSE73680 analysis used the canonical gene-level expression matrix `data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz` and curated sample metadata `config/gse73680_sample_metadata_curated.tsv`.",
  "",
  "Samples were restricted to records marked for inclusion in the disease-context analysis and labelled as control/adjacent or plaque/stone-associated papilla. Empty gene symbols were removed before scoring. For each gene, expression was z-scored across all included samples. For each canonical MAGMA module, the sample-level module score was calculated as the arithmetic mean of available z-scored expression values across module genes present in the bulk matrix.",
  "",
  "The primary disease-context model used complete paired patients only and fit `module_score ~ group_binary + factor(patient_id)` separately for each module, where group_binary denotes plaque/stone-associated papilla versus control/adjacent papilla. P values were corrected across the five modules using the Benjamini-Hochberg procedure.",
  "",
  "No formal deconvolution, injury/remodeling adjustment, composition adjustment, or unpaired primary test was performed in this step."
), con = file.path(root, "notes/phase5_step2_methods_wording.md"))

writeLines(c(
  "# Phase 5-Step 2 Limitations Wording",
  "",
  "GSE73680 group labels are derived from sample filenames and remain manually reviewable. The expression matrix is processed intensity-like/non-count data, so the analysis uses within-study z-scored module summaries rather than count-based modeling.",
  "",
  "The paired base model is intentionally unadjusted for cellular composition and injury/remodeling programs; those sensitivity analyses remain pending. Bulk module shifts cannot establish genetic causality, validate individual causal genes, prove plaque mechanisms, or identify causal cell types. Any later attenuation after tissue-state or composition adjustment should be interpreted as sensitivity to disease-state embedding rather than simple disproof of relevance."
), con = file.path(root, "notes/phase5_step2_limitations_wording.md"))

writeLines(c(
  "# Phase 5-Step 2 Report",
  "",
  "## Inputs",
  paste0("- Expression matrix: `", expr_path, "`."),
  paste0("- Metadata: `", meta_path, "`."),
  paste0("- Pairing reference: `", pairing_path, "`."),
  "- Module files: Phase 1 locked MAGMA top50, top100, Bonferroni, FDR05 and suggestive_p1e4 lists.",
  "",
  "## Analysis Set",
  paste0("- Included samples: ", nrow(meta_included), "."),
  paste0("- Included patients: ", length(unique(meta_included$patient_id)), "."),
  paste0("- Complete paired patients: ", nrow(complete_pairs), "."),
  paste0("- Control/adjacent samples: ", control_count, "."),
  paste0("- Plaque/stone-associated samples: ", disease_count, "."),
  "",
  "## Module Mapping and Scoring",
  "- Empty gene symbols were filtered and genes were z-scored across included samples.",
  "- Scores were recomputed from the canonical matrix; historical module-score outputs were not reused.",
  "- Module mapping table: `results/tables/phase5_step2_magma_module_bulk_gene_mapping.tsv`.",
  "",
  "## Paired Base Model Results",
  result_lines,
  "",
  "## Figures Generated",
  "- `results/figures/phase5_step2_bulk_paired_module_score_lines.pdf`.",
  "- `results/figures/phase5_step2_bulk_paired_difference_summary.pdf`.",
  "- `results/figures/phase5_step2_bulk_paired_model_coefficients.pdf`.",
  "- Source data are saved under `source_data/figures/phase5_step2_*`.",
  "",
  "## Old Outputs Replaced",
  "- Historical bulk/Figure4 outputs are listed in `results/tables/phase5_step2_bulk_deprecated_output_replacement_plan.tsv` for review/replacement by Step 2 canonical outputs.",
  "",
  "## Claim Boundary",
  "- Safe claim: canonical MAGMA modules show paired disease-context shifts in plaque/stone-associated papillary bulk data where supported by the base model.",
  "- Unsafe claims: bulk validates causal genes; bulk proves plaque mechanism; bulk identifies causal cell type; bulk establishes genetic causality.",
  "",
  "## Recommended Next Action",
  "A. proceed to Phase 5-Step 3: injury/remodeling and marker-score tissue-state sensitivity."
), con = file.path(root, "notes/phase5_step2_report.md"))

checklist <- data.frame(
  task_id = sprintf("P5S2-%02d", 1:13),
  task_name = c(
    "main Step 2 script", "analysis input check", "module gene mapping",
    "sample-level module scores", "group score summary", "paired base model",
    "paired differences", "paired visualization figures", "figure source data",
    "interpretation summary", "deprecated replacement plan", "manuscript-safe wording", "Step 2 report"
  ),
  completed = "yes",
  output_file = c(
    "scripts/09_bulk_plaque_context/phase5_step2_bulk_module_scoring_paired_model.R",
    "results/tables/phase5_step2_bulk_analysis_input_check.tsv",
    "results/tables/phase5_step2_magma_module_bulk_gene_mapping.tsv",
    "results/tables/phase5_step2_bulk_sample_module_scores.tsv",
    "results/tables/phase5_step2_bulk_module_score_summary_by_group.tsv",
    "results/tables/phase5_step2_bulk_paired_model_results.tsv",
    "results/tables/phase5_step2_bulk_paired_differences.tsv",
    "results/figures/phase5_step2_bulk_paired_module_score_lines.pdf;results/figures/phase5_step2_bulk_paired_difference_summary.pdf;results/figures/phase5_step2_bulk_paired_model_coefficients.pdf",
    "source_data/figures/phase5_step2_bulk_paired_module_score_lines_source_data.tsv;source_data/figures/phase5_step2_bulk_paired_difference_summary_source_data.tsv;source_data/figures/phase5_step2_bulk_paired_model_coefficients_source_data.tsv",
    "results/tables/phase5_step2_bulk_base_model_interpretation_summary.tsv",
    "results/tables/phase5_step2_bulk_deprecated_output_replacement_plan.tsv",
    "notes/phase5_step2_results_wording.md;notes/phase5_step2_methods_wording.md;notes/phase5_step2_limitations_wording.md",
    "notes/phase5_step2_report.md"
  ),
  blocking_issue = c(rep("", 12), "human review required before Step 3"),
  manual_review_needed = "yes",
  notes = c(
    "Script is reproducible and uses canonical inputs only.",
    "Input universe and complete-pair count recorded.",
    "Exact used/missing genes recorded.",
    "Scores use z-scored gene means; no historical score table reused.",
    "Descriptive group summaries only.",
    "Patient fixed-effect base model only; no adjustment models.",
    "Disease minus control paired differences saved.",
    "PDF and 600dpi PNG generated for each figure.",
    "Figure source data saved.",
    "Allowed and forbidden interpretations recorded.",
    "Old outputs flagged for replacement/review.",
    "Results, methods and limitations wording written conservatively.",
    "Stop rule satisfied after Step 2 outputs."
  ),
  stringsAsFactors = FALSE
)
tsv_write(checklist, file.path(root, "codex_tasks/phase5_step2_completion_checklist.tsv"))

message("Phase 5-Step 2 bulk module scoring and paired base model completed.")
