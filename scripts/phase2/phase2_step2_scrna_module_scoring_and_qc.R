#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

root <- getwd()
object_path <- file.path(root, "data/processed/gse231569_audited_seurat.rds")
table_dir <- file.path(root, "results/tables")
figure_dir <- file.path(root, "results/figures")
source_dir <- file.path(root, "source_data/figures")
note_dir <- file.path(root, "notes")
task_dir <- file.path(root, "codex_tasks")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(note_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

collapse_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  paste(unique(x), collapse = ";")
}

pick_col <- function(cols, patterns) {
  for (pat in patterns) {
    hit <- cols[grepl(pat, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  NA_character_
}

safe_var <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_se <- function(x) {
  x <- as.numeric(x)
  n <- sum(is.finite(x))
  if (n <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(n)
}

pick_expression_layer <- function(obj) {
  assay <- DefaultAssay(obj)
  assay_obj <- obj[[assay]]
  layers <- if ("layers" %in% slotNames(assay_obj)) names(assay_obj@layers) else character()
  if ("data" %in% layers) {
    mat <- LayerData(obj, assay = assay, layer = "data")
    return(list(assay = assay, slot_or_layer = "data", matrix = mat,
                notes = "Selected log-normalized RNA data layer by preference."))
  }
  if ("data" %in% slotNames(assay_obj)) {
    mat <- GetAssayData(obj, assay = assay, slot = "data")
    return(list(assay = assay, slot_or_layer = "data", matrix = mat,
                notes = "Selected log-normalized RNA data slot by preference."))
  }
  if ("counts" %in% layers) {
    mat <- LayerData(obj, assay = assay, layer = "counts")
    return(list(assay = assay, slot_or_layer = "counts", matrix = mat,
                notes = "No data layer detected; counts layer selected as fallback and flagged for review."))
  }
  mat <- GetAssayData(obj, assay = assay, slot = "counts")
  list(assay = assay, slot_or_layer = "counts", matrix = mat,
       notes = "No data slot/layer detected; counts selected as fallback and flagged for review.")
}

module_files <- data.table(
  module_name = c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni", "MAGMA_FDR05", "MAGMA_suggestive_p1e4"),
  module_file = c(
    "results/phase1_step3_magma_gene_sets/MAGMA_top50.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_top100.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_Bonferroni.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_FDR05.txt",
    "results/phase1_step3_magma_gene_sets/MAGMA_suggestive_p1e4.txt"
  )
)

loop_label <- "Loop_of_Henle_TAL"
marker_genes <- c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16", "CLDN14", "CASR", "PKD2", "HIBADH")
known_driver_path <- file.path(root, "data/reference/phase1_step3_known_driver_gene_list_draft.tsv")

obj <- readRDS(object_path)
meta <- as.data.table(obj[[]], keep.rownames = "cell_id")
meta_cols <- setdiff(colnames(meta), "cell_id")
donor_col <- pick_col(meta_cols, c("^donor_id$", "^donor$", "donor", "patient"))
sample_col <- pick_col(meta_cols, c("^sample_id$", "^sample$", "orig.ident", "sample"))
compartment_col <- pick_col(meta_cols, c("^phase1_cell_type$", "broad_cell_type", "compartment", "cell_type", "celltype", "annotation"))
if (is.na(donor_col) || is.na(compartment_col)) {
  blocking <- data.table(
    audit_item = c("Seurat object path", "blocking_issue", "module_scoring_can_proceed"),
    value = c("data/processed/gse231569_audited_seurat.rds", "Required donor or compartment metadata column missing.", "no"),
    source = "Seurat metadata",
    notes = "Step 2 stopped before scoring; metadata must be repaired manually."
  )
  write_tsv(blocking, file.path(table_dir, "phase2_step2_scrna_scoring_input_check.tsv"))
  stop("Blocking metadata issue: donor or compartment column missing.")
}
if (is.na(sample_col)) {
  meta[, sample_id_fallback := get(donor_col)]
  sample_col <- "sample_id_fallback"
}

expr_pick <- pick_expression_layer(obj)
expr <- expr_pick$matrix
features <- rownames(expr)
feature_upper <- toupper(features)
names(features) <- feature_upper
reductions <- names(obj@reductions)
umap_available <- "umap" %in% reductions
donor_comp_possible <- length(unique(meta[[donor_col]])) >= 2 && length(unique(meta[[compartment_col]])) >= 2
loop_present <- any(as.character(meta[[compartment_col]]) == loop_label)

input_check <- data.table(
  audit_item = c(
    "Seurat object path", "object_class", "assays_available", "assay used", "expression slot/layer used",
    "expression_layer_selection_reason", "number of nuclei/cells", "number of features", "donor metadata column",
    "sample metadata column", "compartment metadata column", "Loop/TAL label", "reductions available",
    "UMAP availability", "known-driver draft list", "known-driver handling",
    "donor x compartment summary is possible", "module scoring can proceed"
  ),
  value = c(
    "data/processed/gse231569_audited_seurat.rds", paste(class(obj), collapse = ";"),
    paste(names(obj@assays), collapse = ";"), expr_pick$assay, expr_pick$slot_or_layer,
    expr_pick$notes, ncol(obj), nrow(expr), donor_col, sample_col, compartment_col,
    ifelse(loop_present, loop_label, "not_detected"), paste(reductions, collapse = ";"),
    ifelse(umap_available, "yes", "no"),
    ifelse(file.exists(known_driver_path), "available", "not_found"),
    "annotation_only_not_used_for_removal",
    ifelse(donor_comp_possible, "yes", "no"),
    ifelse(donor_comp_possible && nrow(expr) > 0, "yes", "no")
  ),
  source = c(rep("Seurat object", 8), rep("Seurat metadata", 4), rep("Seurat reductions", 2),
             rep("data/reference/phase1_step3_known_driver_gene_list_draft.tsv", 2), rep("Seurat metadata", 2)),
  notes = c(
    "Canonical Phase 2-Step 1 object.", "", "", "Default assay retained; no annotation altered.",
    "Preference order was log-normalized data before counts.", "", "", "", "", "",
    "Existing annotation field used without relabeling.", "Expected accepted label from Step 1.",
    "", "", "Draft list availability recorded only.", "Driver-removal sensitivity is prohibited in Step 2.",
    "Donors, not nuclei, remain the biological support units.", "Transparent mean-expression scoring."
  )
)
write_tsv(input_check, file.path(table_dir, "phase2_step2_scrna_scoring_input_check.tsv"))

if (!donor_comp_possible) stop("Donor x compartment summary is not possible; stopping before scoring.")

meta_min <- meta[, .(
  cell_id,
  donor_id = as.character(get(donor_col)),
  sample_id = as.character(get(sample_col)),
  compartment = as.character(get(compartment_col))
)]

donor_totals <- meta_min[, .(donor_total = .N), by = donor_id]
donor_counts <- meta_min[, .(n_nuclei = .N), by = .(donor_id, sample_id, compartment)]
donor_counts <- merge(donor_counts, donor_totals, by = "donor_id", all.x = TRUE)
donor_counts[, percent_within_donor := round(100 * n_nuclei / donor_total, 4)]
donor_counts[, notes := ifelse(compartment == loop_label, "Loop/TAL compartment.", "")]
setorder(donor_counts, donor_id, compartment, sample_id)
write_tsv(donor_counts[, .(donor_id, sample_id, compartment, n_nuclei, percent_within_donor, notes)],
          file.path(table_dir, "phase2_step2_donor_compartment_counts.tsv"))

comp_summary <- donor_counts[, .(
  n_nuclei_total = sum(n_nuclei),
  n_donors_present = uniqueN(donor_id),
  min_nuclei_per_donor = min(n_nuclei),
  max_nuclei_per_donor = max(n_nuclei)
), by = compartment]
comp_summary[, candidate_for_interpretation := fifelse(n_donors_present >= 3 & n_nuclei_total >= 50, "yes_descriptive", "needs_review")]
comp_summary[, notes := fifelse(compartment == loop_label & n_nuclei_total == 540 & n_donors_present == 4,
                                "Loop/TAL matches Step 1: 540 nuclei across 4 donors.",
                                fifelse(compartment == loop_label, "Loop/TAL count differs from Step 1; review before interpretation.", ""))]
setorder(comp_summary, compartment)
write_tsv(comp_summary, file.path(table_dir, "phase2_step2_compartment_summary.tsv"))

if (umap_available) {
  umap <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "cell_id")
  setnames(umap, old = colnames(umap)[2:3], new = c("UMAP_1", "UMAP_2"))
  umap_src <- merge(umap, meta_min, by = "cell_id", all.x = TRUE)
  write_tsv_gz(umap_src, file.path(source_dir, "phase2_step2_umap_source_data.tsv.gz"))

  p_comp <- ggplot(umap_src, aes(UMAP_1, UMAP_2, color = compartment)) +
    geom_point(size = 0.08, alpha = 0.55) +
    theme_classic(base_size = 8) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Compartment")
  ggsave(file.path(figure_dir, "phase2_step2_umap_by_compartment.pdf"), p_comp, width = 7.2, height = 5.6)

  p_donor <- ggplot(umap_src, aes(UMAP_1, UMAP_2, color = donor_id)) +
    geom_point(size = 0.08, alpha = 0.55) +
    theme_classic(base_size = 8) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Donor")
  ggsave(file.path(figure_dir, "phase2_step2_umap_by_donor.pdf"), p_donor, width = 7.2, height = 5.6)
} else {
  write_tsv_gz(data.table(), file.path(source_dir, "phase2_step2_umap_source_data.tsv.gz"))
}

write_tsv(donor_counts[, .(donor_id, sample_id, compartment, n_nuclei, percent_within_donor)],
          file.path(source_dir, "phase2_step2_donor_compartment_barplot_source_data.tsv"))
p_bar <- ggplot(donor_counts, aes(donor_id, n_nuclei, fill = compartment)) +
  geom_col(width = 0.72) +
  theme_classic(base_size = 8) +
  labs(x = "Donor", y = "Nuclei", fill = "Compartment")
ggsave(file.path(figure_dir, "phase2_step2_donor_compartment_barplot.pdf"), p_bar, width = 7.0, height = 4.8)

marker_rows <- rbindlist(lapply(marker_genes, function(g) {
  hit <- which(feature_upper == toupper(g))
  data.table(
    gene = g,
    present_in_object = length(hit) > 0,
    feature_id_if_different = ifelse(length(hit) > 0 && features[hit[1]] != g, features[hit[1]], ""),
    notes = "Presence check before descriptive marker summary."
  )
}))
write_tsv(marker_rows, file.path(table_dir, "phase2_step2_loop_tal_marker_presence.tsv"))

marker_present <- marker_rows[present_in_object == TRUE, gene]
marker_summary <- rbindlist(lapply(marker_present, function(g) {
  feature <- features[which(feature_upper == toupper(g))[1]]
  vals <- as.numeric(expr[feature, ])
  dt <- copy(meta_min)
  dt[, expression := vals]
  dt[, gene := g]
  dt[, .(
    present_in_object = TRUE,
    percent_expressing = round(100 * mean(expression > 0, na.rm = TRUE), 4),
    average_expression = mean(expression, na.rm = TRUE),
    n_nuclei = .N
  ), by = .(gene, compartment)]
}), fill = TRUE)
if (nrow(marker_summary) > 0) {
  marker_summary[, mean_rank_desc := frank(-average_expression, ties.method = "min"), by = gene]
  marker_summary[, highest_or_enriched_expression_in_loop_tal := compartment == loop_label & mean_rank_desc == 1]
  marker_summary[, notes := fifelse(highest_or_enriched_expression_in_loop_tal,
                                    "Loop/TAL has the highest average expression descriptively.",
                                    "Descriptive marker summary only.")]
}
write_tsv(marker_summary, file.path(table_dir, "phase2_step2_loop_tal_marker_summary_by_compartment.tsv"))
write_tsv(marker_summary, file.path(source_dir, "phase2_step2_loop_tal_marker_dotplot_source_data.tsv"))
p_marker <- ggplot(marker_summary, aes(gene, compartment)) +
  geom_point(aes(size = percent_expressing, color = average_expression)) +
  scale_color_viridis_c(option = "C") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Marker", y = "Compartment", size = "% expressing", color = "Average expression")
ggsave(file.path(figure_dir, "phase2_step2_loop_tal_marker_dotplot.pdf"), p_marker, width = 7.2, height = 4.8)

module_mapping <- rbindlist(lapply(seq_len(nrow(module_files)), function(i) {
  genes <- readLines(file.path(root, module_files$module_file[i]), warn = FALSE)
  genes <- genes[nzchar(genes)]
  present <- genes[toupper(genes) %in% feature_upper]
  missing <- setdiff(genes, present)
  non_symbol <- genes[grepl("^[0-9]+$", genes)]
  data.table(
    module_name = module_files$module_name[i],
    module_file = module_files$module_file[i],
    n_module_genes = length(genes),
    n_present_in_scrna = length(present),
    n_missing_from_scrna = length(missing),
    present_genes = paste(present, collapse = ";"),
    missing_genes = paste(missing, collapse = ";"),
    non_symbol_gene_ids = paste(non_symbol, collapse = ";"),
    used_for_scoring = ifelse(length(present) >= 5, "yes_present_genes_only", "no_insufficient_present_genes"),
    notes = "Canonical module preserved; missing genes reported and not silently dropped."
  )
}))
write_tsv(module_mapping, file.path(table_dir, "phase2_step2_magma_module_gene_mapping.tsv"))

score_list <- list()
score_summary <- list()
for (i in seq_len(nrow(module_mapping))) {
  module_name <- module_mapping$module_name[i]
  used_genes <- unlist(strsplit(module_mapping$present_genes[i], ";", fixed = TRUE))
  used_genes <- used_genes[nzchar(used_genes)]
  used_features <- features[toupper(used_genes)]
  used_features <- used_features[!is.na(used_features)]
  if (length(used_features) < 1) next
  sub <- expr[used_features, , drop = FALSE]
  score <- Matrix::colMeans(sub)
  detect_frac <- Matrix::colMeans(sub > 0)
  score_dt <- data.table(
    cell_id = colnames(expr),
    module_name = module_name,
    module_score = as.numeric(score),
    detection_fraction = as.numeric(detect_frac)
  )
  score_list[[module_name]] <- score_dt
  score_summary[[module_name]] <- data.table(
    module_name = module_name,
    n_module_genes = module_mapping$n_module_genes[i],
    n_genes_used = length(used_features),
    n_genes_missing = module_mapping$n_missing_from_scrna[i],
    score_method = "arithmetic_mean_log_normalized_expression_across_present_module_genes",
    assay = expr_pick$assay,
    slot_or_layer = expr_pick$slot_or_layer,
    min_score = min(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    notes = "No Seurat AddModuleScore control-gene subtraction used."
  )
}
scores <- rbindlist(score_list)
score_summary_dt <- rbindlist(score_summary)
scores_meta <- merge(scores, meta_min, by = "cell_id", all.x = TRUE)
write_tsv(score_summary_dt, file.path(table_dir, "phase2_step2_nucleus_module_scores_summary.tsv"))
write_tsv_gz(scores_meta[, .(cell_id, donor_id, sample_id, compartment, module_name, module_score)],
             file.path(table_dir, "phase2_step2_nucleus_module_scores.tsv.gz"))

donor_module <- scores_meta[, .(
  n_nuclei = .N,
  mean_module_score = mean(module_score, na.rm = TRUE),
  median_module_score = median(module_score, na.rm = TRUE),
  sd_module_score = safe_var(module_score),
  se_module_score = safe_se(module_score),
  detection_fraction_mean = mean(detection_fraction, na.rm = TRUE)
), by = .(module_name, donor_id, sample_id, compartment)]
donor_module <- merge(donor_module, score_summary_dt[, .(module_name, n_genes_used)], by = "module_name", all.x = TRUE)
donor_module[, notes := "Descriptive donor x compartment summary; no cell-level significance test."]
setcolorder(donor_module, c("module_name", "donor_id", "sample_id", "compartment", "n_nuclei", "n_genes_used",
                            "mean_module_score", "median_module_score", "sd_module_score", "se_module_score",
                            "detection_fraction_mean", "notes"))
setorder(donor_module, module_name, donor_id, compartment)
write_tsv(donor_module, file.path(table_dir, "phase2_step2_donor_compartment_module_scores.tsv"))

loop_other <- rbindlist(lapply(split(donor_module, donor_module$module_name), function(dt) {
  rbindlist(lapply(split(dt, dt$donor_id), function(dd) {
    loop <- dd[compartment == loop_label]
    other <- dd[compartment != loop_label]
    if (nrow(loop) == 0 || nrow(other) == 0) return(NULL)
    ranks <- dd[, .(compartment, mean_module_score)]
    ranks[, rank_desc := frank(-mean_module_score, ties.method = "min")]
    data.table(
      module_name = loop$module_name[1],
      donor_id = loop$donor_id[1],
      loop_tal_mean_score = loop$mean_module_score[1],
      other_compartments_mean_score = stats::weighted.mean(other$mean_module_score, other$n_nuclei, na.rm = TRUE),
      loop_tal_minus_other = loop$mean_module_score[1] - stats::weighted.mean(other$mean_module_score, other$n_nuclei, na.rm = TRUE),
      loop_tal_rank_within_donor = ranks[compartment == loop_label, rank_desc][1],
      n_loop_tal_nuclei = loop$n_nuclei[1],
      n_other_nuclei = sum(other$n_nuclei),
      notes = "Descriptive donor-level contrast; no significance test."
    )
  }))
}))
write_tsv(loop_other, file.path(table_dir, "phase2_step2_loop_tal_vs_other_summary.tsv"))

heat_src <- copy(donor_module)
write_tsv(heat_src, file.path(source_dir, "phase2_step2_donor_compartment_module_heatmap_source_data.tsv"))
p_heat <- ggplot(heat_src, aes(donor_id, compartment, fill = mean_module_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ module_name, ncol = 2, scales = "free_x") +
  scale_fill_viridis_c(option = "C") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) +
  labs(x = "Donor", y = "Compartment", fill = "Mean score")
ggsave(file.path(figure_dir, "phase2_step2_donor_compartment_module_heatmap.pdf"), p_heat, width = 8.2, height = 7.8)

if (umap_available) {
  module_umap_src <- merge(scores_meta[, .(cell_id, donor_id, sample_id, compartment, module_name, module_score)],
                           umap_src[, .(cell_id, UMAP_1, UMAP_2)], by = "cell_id", all.x = TRUE)
  write_tsv_gz(module_umap_src, file.path(source_dir, "phase2_step2_module_score_umap_source_data.tsv.gz"))
  p_score_umap <- ggplot(module_umap_src, aes(UMAP_1, UMAP_2, color = module_score)) +
    geom_point(size = 0.05, alpha = 0.55) +
    facet_wrap(~ module_name, ncol = 2) +
    scale_color_viridis_c(option = "C") +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank()) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Module score")
  ggsave(file.path(figure_dir, "phase2_step2_module_score_umap.pdf"), p_score_umap, width = 8.4, height = 8.0)
} else {
  write_tsv_gz(data.table(), file.path(source_dir, "phase2_step2_module_score_umap_source_data.tsv.gz"))
}

write_tsv_gz(scores_meta[, .(cell_id, donor_id, sample_id, compartment, module_name, module_score)],
             file.path(source_dir, "phase2_step2_module_score_violin_source_data.tsv.gz"))
p_violin <- ggplot(scores_meta, aes(compartment, module_score, fill = compartment)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75) +
  facet_wrap(~ module_name, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", strip.background = element_blank()) +
  labs(x = "Compartment", y = "Module score")
ggsave(file.path(figure_dir, "phase2_step2_module_score_violin_by_compartment.pdf"), p_violin, width = 8.0, height = 10.5)

write_tsv(donor_module, file.path(source_dir, "phase2_step2_donor_level_plot_source_data.tsv"))
p_donor_level <- ggplot(donor_module, aes(compartment, mean_module_score, color = donor_id, group = donor_id)) +
  geom_point(size = 1.7) +
  geom_line(alpha = 0.5, linewidth = 0.25) +
  facet_wrap(~ module_name, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) +
  labs(x = "Compartment", y = "Donor x compartment mean score", color = "Donor")
ggsave(file.path(figure_dir, "phase2_step2_donor_level_boxplot_or_pointplot.pdf"), p_donor_level, width = 8.0, height = 10.2)

replacement_src_path <- file.path(table_dir, "phase2_step1_deprecated_module_usage_audit.tsv")
if (file.exists(replacement_src_path)) {
  dep <- fread(replacement_src_path)
  dep <- dep[needs_rerun_or_review == "yes" & grepl("figure2|scrna|module|source", file_path, ignore.case = TRUE)]
  dep[, new_replacement_output := fifelse(
    grepl("umap", file_path, ignore.case = TRUE), "results/figures/phase2_step2_umap_by_compartment.pdf; source_data/figures/phase2_step2_umap_source_data.tsv.gz",
    fifelse(grepl("donor|celltype|module_score|scores", file_path, ignore.case = TRUE),
            "results/tables/phase2_step2_donor_compartment_module_scores.tsv; results/tables/phase2_step2_loop_tal_vs_other_summary.tsv",
            fifelse(grepl("Figure2|figure2", file_path),
                    "Step 2 provides component replacements only; final Figure 2 assembly remains prohibited in this step.",
                    "Review against Phase 2-Step 2 tables/figures before reuse"))
  )]
  dep[, rerun_status := "replaced_or_component_ready_in_step2"]
  dep[, notes := "Old output should not be reused without checking canonical Step 3 module definitions and Step 2 source data."]
  plan <- dep[, .(
    deprecated_file = file_path,
    likely_old_module = affected_module,
    new_replacement_output,
    rerun_status,
    notes
  )]
  write_tsv(plan, file.path(table_dir, "phase2_step2_deprecated_output_replacement_plan.tsv"))
} else {
  write_tsv(data.table(deprecated_file = character(), likely_old_module = character(), new_replacement_output = character(),
                       rerun_status = character(), notes = character()),
            file.path(table_dir, "phase2_step2_deprecated_output_replacement_plan.tsv"))
}

mapping_line <- paste(module_mapping[, paste0(module_name, " ", n_present_in_scrna, "/", n_module_genes)], collapse = "; ")
loop_status <- comp_summary[compartment == loop_label]
loop_rank_summary <- loop_other[, .(
  donors_with_loop_top_rank = sum(loop_tal_rank_within_donor == 1, na.rm = TRUE),
  donors_total = .N,
  median_loop_minus_other = median(loop_tal_minus_other, na.rm = TRUE)
), by = module_name]
loop_pattern <- paste(loop_rank_summary[, paste0(module_name, ": top-rank donors ", donors_with_loop_top_rank, "/", donors_total,
                                                ", median Loop/TAL-minus-other ", signif(median_loop_minus_other, 4))], collapse = "; ")

figures <- c(
  "results/figures/phase2_step2_umap_by_compartment.pdf",
  "results/figures/phase2_step2_umap_by_donor.pdf",
  "results/figures/phase2_step2_donor_compartment_barplot.pdf",
  "results/figures/phase2_step2_loop_tal_marker_dotplot.pdf",
  "results/figures/phase2_step2_donor_compartment_module_heatmap.pdf",
  "results/figures/phase2_step2_module_score_umap.pdf",
  "results/figures/phase2_step2_module_score_violin_by_compartment.pdf",
  "results/figures/phase2_step2_donor_level_boxplot_or_pointplot.pdf"
)

report <- c(
  "# Phase 2-Step 2 Report",
  "",
  "## Object and Scoring Input",
  "",
  "- Object used: `data/processed/gse231569_audited_seurat.rds`.",
  paste0("- Assay/layer used: `", expr_pick$assay, "` / `", expr_pick$slot_or_layer, "`; ", expr_pick$notes),
  paste0("- Total nuclei/cells: ", ncol(obj), "."),
  paste0("- Donors detected: ", uniqueN(meta_min$donor_id), "."),
  paste0("- Compartment field: `", compartment_col, "`; donor field: `", donor_col, "`."),
  "",
  "## Loop/TAL Status",
  "",
  paste0("- Loop/TAL nuclei and donor coverage: ", loop_status$n_nuclei_total, " nuclei across ", loop_status$n_donors_present, " donors."),
  "- This remains a descriptive cell-context mapping result, not causal cell-type evidence.",
  "",
  "## MAGMA Module Mapping",
  "",
  paste0("- Mapping rates: ", mapping_line, "."),
  "- Missing genes are explicitly reported in `results/tables/phase2_step2_magma_module_gene_mapping.tsv` and were not silently dropped.",
  "",
  "## Score Method and Donor-Level Summary",
  "",
  "- Module score method: arithmetic mean of log-normalized expression across present module genes; no Seurat `AddModuleScore` control-gene subtraction.",
  "- Known-driver draft list status: recorded as annotation-only input; no known-driver removal sensitivity was performed.",
  "- Primary interpretation unit: donor x compartment summaries.",
  paste0("- Descriptive donor-level Loop/TAL-associated pattern: ", loop_pattern, "."),
  "",
  "## Figures Generated",
  "",
  paste0("- `", figures, "`"),
  "",
  "## Deprecated Outputs",
  "",
  "- Existing old Figure 2/snRNA/source-data outputs are superseded for Step 2 components by the new canonical tables and figures listed above.",
  "- Replacement targets are listed in `results/tables/phase2_step2_deprecated_output_replacement_plan.tsv`.",
  "",
  "## Readiness",
  "",
  "- Step 3 can proceed after human review: matched-random benchmark remains intentionally not run in Step 2."
)
writeLines(report, file.path(note_dir, "phase2_step2_report.md"))

limitations <- c(
  "# Phase 2-Step 2 Limitations and Next Steps",
  "",
  "- Donor number limitation: the accepted GSE231569 object contains 4 donors, so donor-level summaries are descriptive and underpowered for broad inference.",
  paste0("- Loop/TAL nuclei limitation: Loop/TAL contains ", loop_status$n_nuclei_total, " nuclei across ", loop_status$n_donors_present, " donors; donor-specific Loop/TAL counts should be reviewed before figure claims."),
  "- Missing module genes: module-specific missingness is reported; scores use only genes present in the snRNA expression matrix.",
  "- No matched-random benchmark was run in this step.",
  "- No known-driver removal sensitivity was run in this step.",
  "- UMAP and violin plots are visualization only and should not be treated as biological significance evidence.",
  "",
  "Recommended next step: A. proceed to Phase 2-Step 3: matched-random benchmark."
)
writeLines(limitations, file.path(note_dir, "phase2_step2_limitations_and_next_steps.md"))

checklist <- data.table(
  task_id = sprintf("P2S2-%02d", 1:14),
  task_name = c(
    "Create main Step 2 R script",
    "Verify metadata and expression layer before scoring",
    "Reconfirm donor x compartment counts",
    "Generate UMAP and donor-compartment QC figures",
    "Generate Loop/TAL marker presence table and dot plot",
    "Map canonical MAGMA modules to snRNA gene universe",
    "Compute nucleus-level module scores",
    "Aggregate donor x compartment module scores",
    "Generate donor x compartment module heatmap",
    "Generate module score UMAPs and violin/point plots",
    "Create deprecated-output replacement plan",
    "Create Step 2 report",
    "Create limitations and next-steps note",
    "Stop before prohibited downstream analyses"
  ),
  completed = "yes",
  output_file = c(
    "scripts/06_scrna_processing/phase2_step2_scrna_module_scoring_and_qc.R",
    "results/tables/phase2_step2_scrna_scoring_input_check.tsv",
    "results/tables/phase2_step2_donor_compartment_counts.tsv; results/tables/phase2_step2_compartment_summary.tsv",
    "results/figures/phase2_step2_umap_by_compartment.pdf; results/figures/phase2_step2_umap_by_donor.pdf; results/figures/phase2_step2_donor_compartment_barplot.pdf",
    "results/tables/phase2_step2_loop_tal_marker_presence.tsv; results/tables/phase2_step2_loop_tal_marker_summary_by_compartment.tsv; results/figures/phase2_step2_loop_tal_marker_dotplot.pdf",
    "results/tables/phase2_step2_magma_module_gene_mapping.tsv",
    "results/tables/phase2_step2_nucleus_module_scores_summary.tsv; results/tables/phase2_step2_nucleus_module_scores.tsv.gz",
    "results/tables/phase2_step2_donor_compartment_module_scores.tsv; results/tables/phase2_step2_loop_tal_vs_other_summary.tsv",
    "results/figures/phase2_step2_donor_compartment_module_heatmap.pdf",
    "results/figures/phase2_step2_module_score_umap.pdf; results/figures/phase2_step2_module_score_violin_by_compartment.pdf; results/figures/phase2_step2_donor_level_boxplot_or_pointplot.pdf",
    "results/tables/phase2_step2_deprecated_output_replacement_plan.tsv",
    "notes/phase2_step2_report.md",
    "notes/phase2_step2_limitations_and_next_steps.md",
    "codex_tasks/phase2_step2_completion_checklist.tsv"
  ),
  blocking_issue = "",
  manual_review_needed = "yes",
  notes = c(rep("Generated within Phase 2-Step 2 boundary.", 13),
            "Matched-random benchmark, driver-removal sensitivity, final Figure 2, spatial/TWAS/bulk and manuscript edits were not performed.")
)
write_tsv(checklist, file.path(task_dir, "phase2_step2_completion_checklist.tsv"))
