#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

root <- getwd()
object_path <- file.path(root, "data/processed/gse231569_audited_seurat.rds")
table_dir <- file.path(root, "results/tables")
figure_dir <- file.path(root, "results/figures")
source_dir <- file.path(root, "source_data/figures")
note_dir <- file.path(root, "notes")
task_dir <- file.path(root, "codex_tasks")
for (d in c(table_dir, figure_dir, source_dir, note_dir, task_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

random_seed <- 20260709L
n_random_requested <- 1000L
set.seed(random_seed)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(as.data.frame(x), con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
collapse_values <- function(x) paste(sort(unique(as.character(x))), collapse = ";")
pick_col <- function(cols, patterns) {
  for (pat in patterns) {
    hit <- cols[grepl(pat, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  NA_character_
}
assign_bins <- function(x, n = 10L) {
  if (length(unique(x[is.finite(x)])) <= 1) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = unique(stats::quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE)),
                 include.lowest = TRUE, labels = FALSE))
}
support_class <- function(percentile, pval, observed, random_median, subset_label) {
  if (!is.finite(observed) || !is.finite(random_median)) return("not_evaluable")
  if (observed <= random_median) return("not_supported")
  if (grepl("exclude_GSM7290914|min_loop_tal_ge20|min_loop_tal_ge50", subset_label) && percentile < 90) {
    return("unstable_due_to_low_cell_count")
  }
  if (percentile >= 95 && pval <= 0.05) return("strong_matched_random_support")
  if (percentile >= 90 && pval <= 0.10) return("moderate_matched_random_support")
  "partial_support"
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
low_count_donor <- "GSM7290914"
known_driver_path <- file.path(root, "data/reference/phase1_step3_known_driver_gene_list_draft.tsv")

obj <- readRDS(object_path)
DefaultAssay(obj) <- "RNA"
expr <- LayerData(obj, assay = "RNA", layer = "data")
features <- rownames(expr)
feature_upper <- toupper(features)
feature_lookup <- features
names(feature_lookup) <- feature_upper

meta <- as.data.table(obj[[]], keep.rownames = "cell_id")
meta_cols <- setdiff(colnames(meta), "cell_id")
donor_col <- pick_col(meta_cols, c("^donor_id$", "^donor$", "donor", "patient"))
sample_col <- pick_col(meta_cols, c("^sample_id$", "^sample$", "orig.ident", "sample"))
compartment_col <- pick_col(meta_cols, c("^phase1_cell_type$", "broad_cell_type", "compartment", "cell_type", "celltype", "annotation"))
if (is.na(donor_col) || is.na(compartment_col)) stop("Required donor or compartment metadata field missing.")
if (is.na(sample_col)) {
  meta[, sample_id_fallback := get(donor_col)]
  sample_col <- "sample_id_fallback"
}
meta_min <- meta[, .(
  cell_id,
  donor_id = as.character(get(donor_col)),
  sample_id = as.character(get(sample_col)),
  compartment = as.character(get(compartment_col))
)]
if (!identical(meta_min$cell_id, colnames(expr))) {
  meta_min <- meta_min[match(colnames(expr), cell_id)]
  if (!identical(meta_min$cell_id, colnames(expr))) stop("Metadata cells do not align with expression matrix columns.")
}

modules <- rbindlist(lapply(seq_len(nrow(module_files)), function(i) {
  genes <- readLines(file.path(root, module_files$module_file[i]), warn = FALSE)
  genes <- genes[nzchar(genes)]
  data.table(
    module_name = module_files$module_name[i],
    module_file = module_files$module_file[i],
    module_gene = genes,
    present_in_scrna = toupper(genes) %in% feature_upper,
    feature = unname(feature_lookup[toupper(genes)])
  )
}))
used_modules <- modules[present_in_scrna == TRUE]
any_module_features <- unique(used_modules$feature)

known_driver_genes <- character()
if (file.exists(known_driver_path)) {
  kd <- fread(known_driver_path)
  gene_cols <- colnames(kd)[grepl("gene|symbol", colnames(kd), ignore.case = TRUE)]
  if (length(gene_cols) > 0) {
    known_driver_genes <- toupper(unique(as.character(unlist(kd[, ..gene_cols], use.names = FALSE))))
    known_driver_genes <- known_driver_genes[!is.na(known_driver_genes) & known_driver_genes != ""]
  }
}

global_mean <- Matrix::rowMeans(expr)
global_detect <- Matrix::rowMeans(expr > 0)
n_cells_detected <- Matrix::rowSums(expr > 0)
gene_universe <- data.table(
  gene = features,
  global_mean_expression = as.numeric(global_mean),
  global_detection_fraction = as.numeric(global_detect),
  n_cells_detected = as.integer(n_cells_detected)
)
gene_universe[, mean_expression_bin := assign_bins(global_mean_expression, 10L)]
gene_universe[, detection_fraction_bin := assign_bins(global_detection_fraction, 10L)]
gene_universe[, in_any_canonical_magma_module := gene %in% any_module_features]
gene_universe[, in_known_driver_draft := toupper(gene) %in% known_driver_genes]
gene_universe[, eligible_for_random_sampling := n_cells_detected > 0]
gene_universe[, notes := fifelse(eligible_for_random_sampling,
                                 "Eligible: detected in at least one nucleus; target module genes excluded only module-wise.",
                                 "Excluded from random sampling: zero detection.")]
write_tsv(gene_universe, file.path(table_dir, "phase2_step3_scrna_gene_universe_metrics.tsv"))

settings <- data.table(
  setting = c("random_seed", "n_random_sets_requested_per_module", "assay", "layer", "donor_field", "compartment_field",
              "loop_tal_label", "score_method", "random_sampling_universe", "target_module_exclusion",
              "known_driver_handling", "matched_bins"),
  value = c(random_seed, n_random_requested, "RNA", "data", donor_col, compartment_col, loop_label,
            "arithmetic_mean_log_normalized_expression_across_present_module_genes",
            "genes with nonzero global detection in the snRNA expression matrix",
            "target observed module genes excluded for that module when sampling",
            "known-driver draft genes annotated but not globally excluded",
            "global mean-expression decile and global detection-fraction decile")
)
write_tsv(settings, file.path(table_dir, "phase2_step3_random_benchmark_settings.tsv"))

eligible <- gene_universe[eligible_for_random_sampling == TRUE]
make_random_sets <- function(observed_features, module_name) {
  target_exclude <- unique(observed_features)
  obs_info <- gene_universe[match(observed_features, gene)]
  candidate_lists <- lapply(seq_along(observed_features), function(i) {
    mb <- obs_info$mean_expression_bin[i]
    db <- obs_info$detection_fraction_bin[i]
    exact <- eligible[mean_expression_bin == mb & detection_fraction_bin == db & !(gene %in% target_exclude), gene]
    near <- eligible[abs(mean_expression_bin - mb) <= 1 & abs(detection_fraction_bin - db) <= 1 & !(gene %in% target_exclude), gene]
    global <- eligible[!(gene %in% target_exclude), gene]
    if (length(exact) > 0) {
      pool <- exact
      status <- "exact_bin"
    } else if (length(near) > 0) {
      pool <- near
      status <- "nearest_bin"
    } else {
      pool <- global
      status <- "unmatched_fallback"
    }
    list(observed_gene = observed_features[i], observed_mean_bin = mb, observed_detection_bin = db,
         pool = pool, global = global, status = status)
  })
  sets <- matrix(NA_character_, nrow = length(observed_features), ncol = n_random_requested)
  status_mat <- matrix(NA_character_, nrow = length(observed_features), ncol = n_random_requested)
  for (i in seq_along(candidate_lists)) {
    cl <- candidate_lists[[i]]
    sets[i, ] <- sample(cl$pool, n_random_requested, replace = length(cl$pool) < n_random_requested)
    status_mat[i, ] <- cl$status
  }
  for (rid in seq_len(n_random_requested)) {
    dup_idx <- which(duplicated(sets[, rid]))
    if (length(dup_idx) > 0) {
      for (i in dup_idx) {
        cl <- candidate_lists[[i]]
        cand <- setdiff(cl$global, sets[, rid])
        if (length(cand) == 0) cand <- cl$global
        sets[i, rid] <- sample(cand, 1)
        status_mat[i, rid] <- paste0(status_mat[i, rid], "_duplicate_repaired")
      }
    }
  }
  manifest <- data.table(
    module_name = rep(module_name, length(observed_features) * n_random_requested),
    random_set_id = rep(sprintf("R%04d", seq_len(n_random_requested)), each = length(observed_features)),
    random_gene = as.vector(sets),
    matched_to_observed_gene = rep(observed_features, times = n_random_requested),
    observed_gene_expression_bin = rep(obs_info$mean_expression_bin, times = n_random_requested),
    observed_gene_detection_bin = rep(obs_info$detection_fraction_bin, times = n_random_requested),
    matching_status = as.vector(status_mat)
  )
  rg_info <- gene_universe[match(manifest$random_gene, gene)]
  manifest[, random_gene_expression_bin := rg_info$mean_expression_bin]
  manifest[, random_gene_detection_bin := rg_info$detection_fraction_bin]
  manifest[, notes := "Matched by global expression/detection bins; duplicate genes repaired within random sets where needed; known-driver draft genes annotated but not excluded."]
  setcolorder(manifest, c("module_name", "random_set_id", "random_gene", "matched_to_observed_gene",
                          "observed_gene_expression_bin", "observed_gene_detection_bin",
                          "random_gene_expression_bin", "random_gene_detection_bin", "matching_status", "notes"))
  list(sets = sets, manifest = manifest)
}

random_set_matrices <- list()
manifest_list <- list()
sampling_summary <- list()
for (mn in module_files$module_name) {
  observed <- used_modules[module_name == mn, feature]
  generated <- make_random_sets(observed, mn)
  random_set_matrices[[mn]] <- generated$sets
  manifest_list[[mn]] <- generated$manifest
  sampling_summary[[mn]] <- data.table(
    module_name = mn,
    n_observed_genes_used = length(observed),
    n_random_sets_requested = n_random_requested,
    n_random_sets_successful = ncol(generated$sets),
    n_failed_matches = generated$manifest[grepl("unmatched|reuse", matching_status), .N],
    fallback_matching_used = any(generated$manifest$matching_status != "exact_bin_unique"),
    random_seed = random_seed,
    notes = "Expression/detection-bin matched random gene sets generated from detected snRNA genes."
  )
}
random_manifest <- rbindlist(manifest_list)
sampling_summary_dt <- rbindlist(sampling_summary)
write_tsv_gz(random_manifest, file.path(table_dir, "phase2_step3_random_gene_set_manifest.tsv.gz"))
write_tsv(sampling_summary_dt, file.path(table_dir, "phase2_step3_random_sampling_summary.tsv"))

group_dt <- unique(meta_min[, .(donor_id, sample_id, compartment)])
setorder(group_dt, donor_id, compartment, sample_id)
group_dt[, group_id := paste(donor_id, sample_id, compartment, sep = "|")]
group_indices <- split(seq_len(nrow(meta_min)), paste(meta_min$donor_id, meta_min$sample_id, meta_min$compartment, sep = "|"))
gene_group_mean <- matrix(NA_real_, nrow = nrow(expr), ncol = length(group_indices),
                          dimnames = list(features, names(group_indices)))
for (gid in names(group_indices)) {
  gene_group_mean[, gid] <- Matrix::rowMeans(expr[, group_indices[[gid]], drop = FALSE])
}
group_counts <- meta_min[, .(n_nuclei = .N), by = .(donor_id, sample_id, compartment)]
group_counts[, group_id := paste(donor_id, sample_id, compartment, sep = "|")]
group_dt <- merge(group_dt, group_counts[, .(group_id, n_nuclei)], by = "group_id", all.x = TRUE)

module_group_scores <- function(gene_sets) {
  if (is.vector(gene_sets)) {
    out <- matrix(NA_real_, nrow = 1, ncol = ncol(gene_group_mean), dimnames = list("observed", colnames(gene_group_mean)))
    genes <- unique(gene_sets[gene_sets %in% rownames(gene_group_mean)])
    for (g in colnames(gene_group_mean)) out[1, g] <- mean(gene_group_mean[genes, g], na.rm = TRUE)
    return(out)
  }
  out <- matrix(NA_real_, nrow = ncol(gene_sets), ncol = ncol(gene_group_mean),
                dimnames = list(sprintf("R%04d", seq_len(ncol(gene_sets))), colnames(gene_group_mean)))
  idx_mat <- matrix(match(gene_sets, rownames(gene_group_mean)), nrow = nrow(gene_sets), ncol = ncol(gene_sets))
  for (g in colnames(gene_group_mean)) {
    vals <- gene_group_mean[, g]
    out[, g] <- colMeans(matrix(vals[idx_mat], nrow = nrow(gene_sets), ncol = ncol(gene_sets)), na.rm = TRUE)
  }
  out
}

calc_stats_from_score_matrix <- function(score_mat, module_name) {
  donors <- sort(unique(group_dt$donor_id))
  contrast_mat <- matrix(NA_real_, nrow = nrow(score_mat), ncol = length(donors), dimnames = list(rownames(score_mat), donors))
  top_mat <- matrix(FALSE, nrow = nrow(score_mat), ncol = length(donors), dimnames = list(rownames(score_mat), donors))
  loop_counts <- integer(length(donors))
  names(loop_counts) <- donors
  for (d in donors) {
    loop_gid <- group_dt[donor_id == d & compartment == loop_label, group_id]
    other <- group_dt[donor_id == d & compartment != loop_label]
    if (length(loop_gid) != 1 || nrow(other) == 0) next
    weights <- other$n_nuclei / sum(other$n_nuclei)
    loop_score <- score_mat[, loop_gid]
    other_score <- as.numeric(score_mat[, other$group_id, drop = FALSE] %*% weights)
    contrast_mat[, d] <- loop_score - other_score
    top_mat[, d] <- loop_score >= apply(score_mat[, c(loop_gid, other$group_id), drop = FALSE], 1, max, na.rm = TRUE)
    loop_counts[d] <- group_dt[group_id == loop_gid, n_nuclei]
  }
  subset_list <- list(
    full_4_donors = donors,
    exclude_GSM7290914 = setdiff(donors, low_count_donor),
    min_loop_tal_ge20 = names(loop_counts)[loop_counts >= 20],
    min_loop_tal_ge50 = names(loop_counts)[loop_counts >= 50]
  )
  for (d in donors) subset_list[[paste0("leave_one_out_exclude_", d)]] <- setdiff(donors, d)
  rbindlist(lapply(names(subset_list), function(subset_name) {
    inc <- subset_list[[subset_name]]
    inc <- inc[inc %in% colnames(contrast_mat)]
    if (length(inc) == 0) return(NULL)
    med <- apply(contrast_mat[, inc, drop = FALSE], 1, median, na.rm = TRUE)
    avg <- rowMeans(contrast_mat[, inc, drop = FALSE], na.rm = TRUE)
    topn <- rowSums(top_mat[, inc, drop = FALSE], na.rm = TRUE)
    data.table(
      module_name = module_name,
      random_set_id = rownames(score_mat),
      donor_subset = subset_name,
      n_donors = length(inc),
      included_donors = collapse_values(inc),
      excluded_donors = collapse_values(setdiff(donors, inc)),
      statistic_name = rep(c("median_loop_tal_minus_other", "mean_loop_tal_minus_other", "n_donors_loop_tal_top_rank"),
                           each = nrow(score_mat)),
      value = c(med, avg, topn),
      notes = "Donor-level descriptive statistic; no cell-level significance test."
    )
  }))
}

observed_stats <- rbindlist(lapply(module_files$module_name, function(mn) {
  score_mat <- module_group_scores(used_modules[module_name == mn, feature])
  calc_stats_from_score_matrix(score_mat, mn)
}))
observed_out <- observed_stats[, .(
  module_name, donor_subset, n_donors, included_donors, excluded_donors,
  statistic_name, observed_value = value, notes
)]
write_tsv(observed_out, file.path(table_dir, "phase2_step3_observed_loop_tal_statistics.tsv"))

random_stats <- rbindlist(lapply(module_files$module_name, function(mn) {
  score_mat <- module_group_scores(random_set_matrices[[mn]])
  calc_stats_from_score_matrix(score_mat, mn)
}))
random_out <- random_stats[, .(
  module_name, random_set_id, donor_subset, n_donors, included_donors,
  statistic_name, random_value = value
)]
write_tsv_gz(random_out, file.path(table_dir, "phase2_step3_random_loop_tal_statistics.tsv.gz"))

obs_key <- observed_out[, .(module_name, donor_subset, statistic_name, observed_value)]
random_summary <- merge(random_out, obs_key, by = c("module_name", "donor_subset", "statistic_name"), all.x = TRUE)
benchmark <- random_summary[, .(
  n_random_sets_successful = uniqueN(random_set_id),
  observed_value = observed_value[1],
  random_median = median(random_value, na.rm = TRUE),
  random_p025 = as.numeric(quantile(random_value, 0.025, na.rm = TRUE)),
  random_p975 = as.numeric(quantile(random_value, 0.975, na.rm = TRUE)),
  observed_minus_random_median = observed_value[1] - median(random_value, na.rm = TRUE),
  observed_percentile = mean(random_value <= observed_value[1], na.rm = TRUE) * 100,
  empirical_p_one_sided = (sum(random_value >= observed_value[1], na.rm = TRUE) + 1) / (.N + 1)
), by = .(module_name, donor_subset, statistic_name)]
benchmark <- merge(benchmark, sampling_summary_dt[, .(module_name, n_observed_genes_used)], by = "module_name", all.x = TRUE)
benchmark[, benchmark_support_class := mapply(support_class, observed_percentile, empirical_p_one_sided,
                                             observed_value, random_median, donor_subset)]
benchmark[, notes := "Matched random benchmark; interpret donor-level statistics conservatively."]
setcolorder(benchmark, c("module_name", "donor_subset", "statistic_name", "n_observed_genes_used",
                         "n_random_sets_successful", "observed_value", "random_median", "random_p025",
                         "random_p975", "observed_minus_random_median", "observed_percentile",
                         "empirical_p_one_sided", "benchmark_support_class", "notes"))
write_tsv(benchmark, file.path(table_dir, "phase2_step3_matched_random_benchmark_summary.tsv"))

sens_stat <- "median_loop_tal_minus_other"
sensitivity <- dcast(observed_out[statistic_name == sens_stat],
                     module_name + statistic_name ~ donor_subset, value.var = "observed_value")
rename_map <- c(full_4_donors = "full_4_donor_value", exclude_GSM7290914 = "exclude_GSM7290914_value",
                min_loop_tal_ge20 = "donors_loop_tal_ge20_value", min_loop_tal_ge50 = "donors_loop_tal_ge50_value")
for (old in names(rename_map)) if (old %in% names(sensitivity)) setnames(sensitivity, old, rename_map[[old]])
for (nm in unname(rename_map)) if (!nm %in% names(sensitivity)) sensitivity[, (nm) := NA_real_]
sensitivity[, interpretation := fifelse(is.na(donors_loop_tal_ge50_value), "insufficient_donors_for_strong_claim",
                                 fifelse(exclude_GSM7290914_value >= 0.75 * full_4_donor_value &
                                           donors_loop_tal_ge20_value >= 0.75 * full_4_donor_value,
                                         "robust_to_low_cell_count_donor",
                                 fifelse(exclude_GSM7290914_value > 0,
                                         "attenuated_after_excluding_low_cell_count_donor",
                                         "dependent_on_low_cell_count_donor")))]
sensitivity[, notes := "GSM7290914 has 4 Loop/TAL nuclei; sensitivity uses median Loop/TAL-minus-other."]
write_tsv(sensitivity[, .(module_name, statistic_name, full_4_donor_value, exclude_GSM7290914_value,
                          donors_loop_tal_ge20_value, donors_loop_tal_ge50_value, interpretation, notes)],
          file.path(table_dir, "phase2_step3_low_cell_count_sensitivity.tsv"))

dist_src <- random_out[statistic_name == sens_stat & donor_subset %in% c("full_4_donors", "min_loop_tal_ge20")]
obs_dist <- observed_out[statistic_name == sens_stat & donor_subset %in% c("full_4_donors", "min_loop_tal_ge20"),
                         .(module_name, donor_subset, observed_value)]
write_tsv_gz(merge(dist_src, obs_dist, by = c("module_name", "donor_subset"), all.x = TRUE),
             file.path(source_dir, "phase2_step3_matched_random_distribution_source_data.tsv.gz"))
p_dist <- ggplot(dist_src, aes(random_value)) +
  geom_histogram(bins = 40, fill = "grey75", color = "white", linewidth = 0.1) +
  geom_vline(data = obs_dist, aes(xintercept = observed_value), color = "#C44E52", linewidth = 0.5) +
  facet_grid(donor_subset ~ module_name, scales = "free_x") +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Random median Loop/TAL-minus-other", y = "Random sets")
ggsave(file.path(figure_dir, "phase2_step3_matched_random_distribution.pdf"), p_dist, width = 10.5, height = 5.2)

sens_long <- melt(sensitivity[, .(module_name, statistic_name, full_4_donor_value, exclude_GSM7290914_value,
                                  donors_loop_tal_ge20_value, donors_loop_tal_ge50_value, interpretation)],
                  id.vars = c("module_name", "statistic_name", "interpretation"),
                  variable.name = "donor_subset", value.name = "observed_value")
write_tsv(sens_long, file.path(source_dir, "phase2_step3_low_cell_count_sensitivity_source_data.tsv"))
p_sens <- ggplot(sens_long, aes(donor_subset, observed_value, group = module_name, color = module_name)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.35) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Donor subset", y = "Observed median Loop/TAL-minus-other", color = "Module")
ggsave(file.path(figure_dir, "phase2_step3_low_cell_count_sensitivity_plot.pdf"), p_sens, width = 7.5, height = 4.8)

rank_src <- random_out[statistic_name == "n_donors_loop_tal_top_rank" & donor_subset %in% c("full_4_donors", "min_loop_tal_ge20")]
rank_obs <- observed_out[statistic_name == "n_donors_loop_tal_top_rank" & donor_subset %in% c("full_4_donors", "min_loop_tal_ge20"),
                         .(module_name, donor_subset, observed_value)]
write_tsv_gz(merge(rank_src, rank_obs, by = c("module_name", "donor_subset"), all.x = TRUE),
             file.path(source_dir, "phase2_step3_rank_metric_saturation_source_data.tsv.gz"))
p_rank <- ggplot(rank_src, aes(factor(random_value))) +
  geom_bar(fill = "grey70") +
  geom_vline(data = rank_obs, aes(xintercept = observed_value + 1), color = "#C44E52", linewidth = 0.5) +
  facet_grid(donor_subset ~ module_name) +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank()) +
  labs(x = "Random donors with Loop/TAL top rank", y = "Random sets")
ggsave(file.path(figure_dir, "phase2_step3_rank_metric_saturation_plot.pdf"), p_rank, width = 10.5, height = 5.2)

primary_summary <- benchmark[donor_subset %in% c("full_4_donors", "min_loop_tal_ge20", "min_loop_tal_ge50") &
                               statistic_name == sens_stat]
rank_summary <- benchmark[donor_subset == "full_4_donors" & statistic_name == "n_donors_loop_tal_top_rank"]
support_lines <- primary_summary[, paste0(
  "- ", module_name, " / ", donor_subset, ": observed=", signif(observed_value, 4),
  ", random median=", signif(random_median, 4), ", percentile=", signif(observed_percentile, 4),
  ", empirical P=", signif(empirical_p_one_sided, 4), ", class=", benchmark_support_class, "."
)]
rank_lines <- rank_summary[, paste0(
  "- ", module_name, ": observed top-rank donors=", observed_value,
  ", random median=", signif(random_median, 4), "; rank metric may show ceiling/saturation."
)]

writeLines(c(
  "# Phase 2-Step 3 Manuscript-Safe Interpretation Wording",
  "",
  "## Results Wording",
  "",
  "Expression- and detection-matched random gene-set benchmarks were used to contextualize the donor-level Loop/TAL-versus-other compartment score differences for canonical MAGMA genetic-priority modules. These results should be described as matched-benchmark support for renal papillary cell-context mapping, not as causal cell-type assignment.",
  "",
  "## Limitations Wording",
  "",
  "The benchmark remains limited by four donors and uneven Loop/TAL representation, including one donor with only four Loop/TAL nuclei. Rank-based summaries can be affected by ceiling effects and low-cell-count sensitivity; donor-level Loop/TAL-minus-other contrasts should be interpreted descriptively.",
  "",
  "## Figure Legend Wording",
  "",
  "Matched random distributions show expression- and detection-matched random gene-set statistics; vertical lines denote observed canonical MAGMA module statistics. UMAP and rank panels are contextual visualizations, whereas donor-level Loop/TAL-minus-other contrasts provide the primary descriptive benchmark statistic.",
  "",
  "## Prohibited Wording",
  "",
  "- Do not write: validated Loop/TAL enrichment.",
  "- Do not write: causal Loop/TAL cell type.",
  "- Do not write: independent validation.",
  "- Do not write: strong cell-type specificity unless support remains strong across low-cell-count sensitivity subsets."
), file.path(note_dir, "phase2_step3_interpretation_wording.md"))

writeLines(c(
  "# Phase 2-Step 3 Report",
  "",
  "## Inputs",
  "",
  "- Object used: `data/processed/gse231569_audited_seurat.rds`.",
  "- Assay/layer used: `RNA` / `data`.",
  paste0("- Modules used: ", paste(module_files$module_name, collapse = ", "), "."),
  paste0("- Random seed: ", random_seed, "; requested random sets per module: ", n_random_requested, "."),
  "",
  "## Random Matching Strategy",
  "",
  "- Gene universe: detected genes in the snRNA expression matrix.",
  "- Matching: global mean-expression decile plus global detection-fraction decile.",
  "- Target observed module genes were excluded while sampling that module; known-driver draft genes were annotated but not removed.",
  "- Scoring method matched Step 2: arithmetic mean of log-normalized expression across present module genes.",
  "",
  "## Random Sets Generated",
  "",
  paste0("- ", sampling_summary_dt$module_name, ": ", sampling_summary_dt$n_random_sets_successful, " random sets; observed genes used=", sampling_summary_dt$n_observed_genes_used, "."),
  "",
  "## Benchmark Summary",
  "",
  support_lines,
  "",
  "## Rank Metric Saturation",
  "",
  rank_lines,
  "",
  "## Low-Cell-Count Sensitivity",
  "",
  paste0("- ", sensitivity$module_name, ": ", sensitivity$interpretation, "; full=", signif(sensitivity$full_4_donor_value, 4),
         ", excluding GSM7290914=", signif(sensitivity$exclude_GSM7290914_value, 4),
         ", Loop/TAL >=20=", signif(sensitivity$donors_loop_tal_ge20_value, 4),
         ", Loop/TAL >=50=", signif(sensitivity$donors_loop_tal_ge50_value, 4), "."),
  "",
  "## Claim Boundary",
  "",
  "- Safe claim: expression- and detection-matched random benchmarks contextualize a descriptive donor-level Loop/TAL-associated pattern.",
  "- Unsafe claim: causal Loop/TAL cell-type assignment, validated enrichment, independent validation, or plaque-specific localization.",
  "",
  "## Readiness",
  "",
  "- Step 4 can proceed after human review: known-driver removal sensitivity remains intentionally not run in Step 3."
), file.path(note_dir, "phase2_step3_report.md"))

writeLines(c(
  "# Phase 2-Step 3 Limitations and Next Steps",
  "",
  "- Donor number limitation: only 4 donors are available.",
  "- Low-cell-count limitation: GSM7290914 has only 4 Loop/TAL nuclei.",
  "- Matching limitations: random sets are matched by binned global expression and detection; nearest-bin fallback may occur where exact bins are sparse.",
  "- Known-driver removal sensitivity was not run in this step.",
  "- Spatial validation was not run in this step.",
  "- Recommended next step: A. proceed to Phase 2-Step 4: known-driver removal sensitivity."
), file.path(note_dir, "phase2_step3_limitations_and_next_steps.md"))

checklist <- data.table(
  task_id = sprintf("P2S3-%02d", 1:11),
  task_name = c(
    "Create matched-random benchmark script",
    "Define eligible snRNA gene universe",
    "Build expression/detection-matched random gene sets",
    "Compute observed and random Loop/TAL donor-level statistics",
    "Calculate empirical benchmark percentiles and P values",
    "Create low-cell-count sensitivity summary",
    "Generate benchmark figures and source data",
    "Create manuscript-safe interpretation wording",
    "Create Step 3 report",
    "Create limitations and next-steps note",
    "Stop before prohibited downstream analyses"
  ),
  completed = "yes",
  output_file = c(
    "scripts/06_scrna_processing/phase2_step3_scrna_matched_random_benchmark.R",
    "results/tables/phase2_step3_scrna_gene_universe_metrics.tsv",
    "results/tables/phase2_step3_random_gene_set_manifest.tsv.gz; results/tables/phase2_step3_random_sampling_summary.tsv",
    "results/tables/phase2_step3_observed_loop_tal_statistics.tsv; results/tables/phase2_step3_random_loop_tal_statistics.tsv.gz",
    "results/tables/phase2_step3_matched_random_benchmark_summary.tsv",
    "results/tables/phase2_step3_low_cell_count_sensitivity.tsv",
    "results/figures/phase2_step3_matched_random_distribution.pdf; results/figures/phase2_step3_low_cell_count_sensitivity_plot.pdf; results/figures/phase2_step3_rank_metric_saturation_plot.pdf",
    "notes/phase2_step3_interpretation_wording.md",
    "notes/phase2_step3_report.md",
    "notes/phase2_step3_limitations_and_next_steps.md",
    "codex_tasks/phase2_step3_completion_checklist.tsv"
  ),
  blocking_issue = "",
  manual_review_needed = "yes",
  notes = c(rep("Generated within Phase 2-Step 3 boundary.", 10),
            "Known-driver removal, final Figure 2, spatial/TWAS/bulk and manuscript edits were not performed.")
)
write_tsv(checklist, file.path(task_dir, "phase2_step3_completion_checklist.tsv"))
