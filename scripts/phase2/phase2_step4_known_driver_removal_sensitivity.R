#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

root <- getwd()
table_dir <- file.path(root, "results/tables")
figure_dir <- file.path(root, "results/figures")
source_dir <- file.path(root, "source_data/figures")
note_dir <- file.path(root, "notes")
task_dir <- file.path(root, "codex_tasks")
gene_set_dir <- file.path(root, "results/phase2_step4_driver_removed_gene_sets")
for (d in c(table_dir, figure_dir, source_dir, note_dir, task_dir, gene_set_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(root, "data/processed/gse231569_audited_seurat.rds")
driver_path <- file.path(root, "data/reference/phase1_step3_known_driver_gene_list_draft.tsv")
driver_used_path <- file.path(root, "data/reference/phase2_step4_known_driver_gene_list_used.tsv")
random_seed <- 20260709L + 4L
n_random_requested <- 1000L
set.seed(random_seed)
loop_label <- "Loop_of_Henle_TAL"
low_count_donor <- "GSM7290914"

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
safe_sd <- function(x) if (sum(is.finite(x)) > 1) sd(x, na.rm = TRUE) else NA_real_
safe_se <- function(x) if (sum(is.finite(x)) > 1) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))) else NA_real_
assign_bins <- function(x, n = 10L) {
  if (length(unique(x[is.finite(x)])) <= 1) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = unique(stats::quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE)),
                 include.lowest = TRUE, labels = FALSE))
}
support_class <- function(percentile, pval, observed, random_median, n_genes, subset_label) {
  if (n_genes < 10) return("unstable_due_to_small_gene_set")
  if (grepl("min_loop_tal_ge|exclude_GSM7290914", subset_label) && percentile < 90) return("unstable_due_to_low_cell_count")
  if (!is.finite(observed) || !is.finite(random_median)) return("not_supported")
  if (observed <= random_median) return("not_supported")
  if (percentile >= 95 && pval <= 0.05) return("strong_matched_random_support")
  if (percentile >= 90 && pval <= 0.10) return("moderate_matched_random_support")
  "partial_support"
}

initial_driver <- data.table(
  gene = c("UMOD", "SLC12A1", "CLDN10", "CLDN14", "KCNJ1", "CLDN16", "CASR", "PKD2", "HIBADH"),
  category = c("known KSD gene; Loop/TAL marker", "Loop/TAL marker; renal transport", "Loop/TAL marker; renal transport",
               "known KSD gene; calcium handling", "Loop/TAL marker; renal transport", "renal transport; calcium handling",
               "known KSD gene; calcium handling", "known KSD gene; calcium handling", "curated exemplar"),
  reason_for_inclusion = "Required initial Phase 2-Step 4 sensitivity driver gene.",
  source_or_basis = "Phase 2-Step 4 instruction",
  recommended_for_driver_removal_sensitivity = "yes",
  human_review_required = "yes"
)
if (file.exists(driver_path)) {
  driver <- fread(driver_path)
} else {
  driver <- copy(initial_driver)
}
if (!"gene" %in% names(driver)) setnames(driver, names(driver)[1], "gene")
if (!"category" %in% names(driver)) driver[, category := "curated exemplar"]
driver[, gene := toupper(as.character(gene))]
driver <- unique(rbindlist(list(driver, initial_driver[, lapply(.SD, as.character)]), fill = TRUE), by = "gene")
driver[is.na(recommended_for_driver_removal_sensitivity) | recommended_for_driver_removal_sensitivity == "", recommended_for_driver_removal_sensitivity := "yes"]
write_tsv(driver, driver_used_path)

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
if (is.na(donor_col) || is.na(compartment_col)) stop("Required metadata fields are missing.")
if (is.na(sample_col)) {
  meta[, sample_id_fallback := get(donor_col)]
  sample_col <- "sample_id_fallback"
}
meta_min <- meta[, .(cell_id, donor_id = as.character(get(donor_col)), sample_id = as.character(get(sample_col)), compartment = as.character(get(compartment_col)))]
if (!identical(meta_min$cell_id, colnames(expr))) {
  meta_min <- meta_min[match(colnames(expr), cell_id)]
  if (!identical(meta_min$cell_id, colnames(expr))) stop("Metadata and expression columns do not align.")
}

modules <- rbindlist(lapply(seq_len(nrow(module_files)), function(i) {
  genes <- readLines(file.path(root, module_files$module_file[i]), warn = FALSE)
  genes <- genes[nzchar(genes)]
  data.table(module_name = module_files$module_name[i], module_file = module_files$module_file[i],
             module_gene = genes, gene_upper = toupper(genes),
             present_in_scrna = toupper(genes) %in% feature_upper,
             feature = unname(feature_lookup[toupper(genes)]))
}))
driver_genes <- unique(driver[recommended_for_driver_removal_sensitivity == "yes", gene])
present_driver_scrna <- driver_genes[driver_genes %in% feature_upper]

audit <- copy(driver[, .(gene, category)])
audit[, present_in_driver_list := TRUE]
audit[, present_in_scrna := gene %in% feature_upper]
for (mn in module_files$module_name) {
  audit[, (paste0("present_in_", mn)) := gene %in% modules[module_name == mn, gene_upper]]
}
audit[, included_for_removal := present_in_scrna & Reduce(`|`, lapply(module_files$module_name, function(mn) get(paste0("present_in_", mn))))]
audit[, reason := fifelse(included_for_removal, "Present in driver list, snRNA expression matrix and at least one canonical MAGMA module.",
                   fifelse(!present_in_scrna, "Driver gene not present in snRNA expression matrix.",
                           "Driver gene not present in canonical MAGMA modules."))]
audit[, notes := "Sensitivity only; not interpreted as causal gene evidence."]
write_tsv(audit, file.path(table_dir, "phase2_step4_known_driver_gene_audit.tsv"))

driver_removed_modules <- list()
module_summary <- rbindlist(lapply(module_files$module_name, function(mn) {
  dt <- modules[module_name == mn]
  removable <- dt[present_in_scrna == TRUE & gene_upper %in% present_driver_scrna, gene_upper]
  canonical_after <- dt[!(gene_upper %in% removable), module_gene]
  present_after <- dt[present_in_scrna == TRUE & !(gene_upper %in% removable), feature]
  driver_removed_modules[[mn]] <<- present_after
  out_path <- file.path(gene_set_dir, paste0(mn, "_driver_removed.txt"))
  writeLines(canonical_after, out_path)
  n_before <- dt[present_in_scrna == TRUE, .N]
  n_removed <- length(removable)
  n_after <- length(present_after)
  data.table(
    module_name = mn,
    n_canonical_genes = nrow(dt),
    n_present_in_scrna_before_removal = n_before,
    driver_genes_removed = paste(removable, collapse = ";"),
    n_driver_genes_removed = n_removed,
    n_present_in_scrna_after_removal = n_after,
    percent_genes_removed = ifelse(n_before > 0, 100 * n_removed / n_before, NA_real_),
    sufficient_genes_remaining = ifelse(n_after >= 10 && (n_removed / max(n_before, 1)) <= 0.30, "yes", "unstable_review"),
    notes = ifelse(n_after == 0, "No genes remain after driver removal; downstream scores unavailable.",
                   "Driver removal sensitivity only; canonical module remains primary.")
  )
}))
write_tsv(module_summary, file.path(table_dir, "phase2_step4_driver_removed_module_summary.tsv"))

score_cells <- rbindlist(lapply(names(driver_removed_modules), function(mn) {
  genes <- driver_removed_modules[[mn]]
  if (length(genes) == 0) return(NULL)
  score <- Matrix::colMeans(expr[genes, , drop = FALSE])
  data.table(cell_id = colnames(expr), module_name = mn, module_score = as.numeric(score))
}))
score_cells <- merge(score_cells, meta_min, by = "cell_id", all.x = TRUE)
donor_scores <- score_cells[, .(
  n_nuclei = .N,
  n_genes_used_after_driver_removal = length(driver_removed_modules[[module_name[1]]]),
  mean_module_score = mean(module_score, na.rm = TRUE),
  median_module_score = median(module_score, na.rm = TRUE),
  sd_module_score = safe_sd(module_score),
  se_module_score = safe_se(module_score)
), by = .(module_name, donor_id, sample_id, compartment)]
donor_scores[, notes := "Driver-removed module score; donor x compartment descriptive summary."]
setcolorder(donor_scores, c("module_name", "donor_id", "sample_id", "compartment", "n_nuclei",
                            "n_genes_used_after_driver_removal", "mean_module_score", "median_module_score",
                            "sd_module_score", "se_module_score", "notes"))
write_tsv(donor_scores, file.path(table_dir, "phase2_step4_driver_removed_donor_compartment_module_scores.tsv"))

loop_summary <- rbindlist(lapply(split(donor_scores, donor_scores$module_name), function(dt) {
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
      other_compartments_mean_score = weighted.mean(other$mean_module_score, other$n_nuclei, na.rm = TRUE),
      loop_tal_minus_other = loop$mean_module_score[1] - weighted.mean(other$mean_module_score, other$n_nuclei, na.rm = TRUE),
      loop_tal_rank_within_donor = ranks[compartment == loop_label, rank_desc][1],
      n_loop_tal_nuclei = loop$n_nuclei[1],
      n_other_nuclei = sum(other$n_nuclei),
      notes = "Driver-removed descriptive donor-level contrast; no cell-level significance test."
    )
  }))
}))
write_tsv(loop_summary, file.path(table_dir, "phase2_step4_driver_removed_loop_tal_vs_other_summary.tsv"))

calc_observed_stats <- function(loop_dt) {
  subsets <- list(
    full_4_donors = loop_dt$donor_id,
    exclude_GSM7290914 = setdiff(loop_dt$donor_id, low_count_donor),
    min_loop_tal_ge20 = loop_dt[n_loop_tal_nuclei >= 20, donor_id],
    min_loop_tal_ge50 = loop_dt[n_loop_tal_nuclei >= 50, donor_id]
  )
  for (d in loop_dt$donor_id) subsets[[paste0("leave_one_out_exclude_", d)]] <- setdiff(loop_dt$donor_id, d)
  rbindlist(lapply(names(subsets), function(sn) {
    inc <- subsets[[sn]]
    dd <- loop_dt[donor_id %in% inc]
    if (nrow(dd) == 0) return(NULL)
    data.table(
      donor_subset = sn,
      n_donors = nrow(dd),
      included_donors = collapse_values(dd$donor_id),
      excluded_donors = collapse_values(setdiff(loop_dt$donor_id, dd$donor_id)),
      statistic_name = c("median_loop_tal_minus_other", "mean_loop_tal_minus_other", "n_donors_loop_tal_top_rank"),
      value = c(median(dd$loop_tal_minus_other, na.rm = TRUE),
                mean(dd$loop_tal_minus_other, na.rm = TRUE),
                sum(dd$loop_tal_rank_within_donor == 1, na.rm = TRUE))
    )
  }))
}
driver_removed_stats <- rbindlist(lapply(split(loop_summary, loop_summary$module_name), function(dt) {
  out <- calc_observed_stats(dt)
  out[, module_name := dt$module_name[1]]
  out
}))
setcolorder(driver_removed_stats, c("module_name", "donor_subset", "n_donors", "included_donors", "excluded_donors", "statistic_name", "value"))

orig_loop <- fread(file.path(table_dir, "phase2_step2_loop_tal_vs_other_summary.tsv"))
orig_stats <- rbindlist(lapply(split(orig_loop, orig_loop$module_name), function(dt) {
  out <- calc_observed_stats(dt)
  out[, module_name := dt$module_name[1]]
  out
}))
stat_map <- data.table(
  donor_subset = c("full_4_donors", "full_4_donors", "full_4_donors", "exclude_GSM7290914", "min_loop_tal_ge20", "min_loop_tal_ge50"),
  statistic_name = c("median_loop_tal_minus_other", "mean_loop_tal_minus_other", "n_donors_loop_tal_top_rank",
                     "median_loop_tal_minus_other", "median_loop_tal_minus_other", "median_loop_tal_minus_other"),
  statistic_label = c("median LoopTAL-minus-other across all 4 donors", "mean LoopTAL-minus-other across all 4 donors",
                      "number of donors where Loop/TAL is top-ranked", "median LoopTAL-minus-other excluding GSM7290914",
                      "median LoopTAL-minus-other among donors with Loop/TAL nuclei >=20",
                      "median LoopTAL-minus-other among donors with Loop/TAL nuclei >=50")
)
comparison <- merge(
  orig_stats[, .(module_name, donor_subset, statistic_name, original_value = value)],
  driver_removed_stats[, .(module_name, donor_subset, statistic_name, driver_removed_value = value)],
  by = c("module_name", "donor_subset", "statistic_name"), all = FALSE
)
comparison <- merge(comparison, stat_map, by = c("donor_subset", "statistic_name"), all.y = TRUE, allow.cartesian = TRUE)
comparison[, absolute_change := driver_removed_value - original_value]
comparison[, percent_change := ifelse(original_value != 0, 100 * absolute_change / abs(original_value), NA_real_)]
comparison <- merge(comparison, module_summary[, .(module_name, sufficient_genes_remaining, n_present_in_scrna_after_removal)], by = "module_name", all.x = TRUE)
comparison[, interpretation := fifelse(sufficient_genes_remaining != "yes", "unstable_due_to_small_remaining_gene_set",
                                fifelse(driver_removed_value <= 0, "largely_driver_dependent",
                                fifelse(driver_removed_value >= 0.70 * original_value, "retained_after_driver_removal",
                                fifelse(driver_removed_value > 0, "attenuated_but_retained", "insufficient_data"))))]
comparison[, notes := "Original canonical module remains primary; this is sensitivity analysis only."]
write_tsv(comparison[, .(module_name, statistic_name = statistic_label, original_value, driver_removed_value,
                         absolute_change, percent_change, interpretation, notes)],
          file.path(table_dir, "phase2_step4_original_vs_driver_removed_summary.tsv"))

group_dt <- unique(meta_min[, .(donor_id, sample_id, compartment)])
setorder(group_dt, donor_id, compartment, sample_id)
group_dt[, group_id := paste(donor_id, sample_id, compartment, sep = "|")]
group_indices <- split(seq_len(nrow(meta_min)), paste(meta_min$donor_id, meta_min$sample_id, meta_min$compartment, sep = "|"))
gene_group_mean <- matrix(NA_real_, nrow = nrow(expr), ncol = length(group_indices), dimnames = list(features, names(group_indices)))
for (gid in names(group_indices)) gene_group_mean[, gid] <- Matrix::rowMeans(expr[, group_indices[[gid]], drop = FALSE])
group_counts <- meta_min[, .(n_nuclei = .N), by = .(donor_id, sample_id, compartment)]
group_counts[, group_id := paste(donor_id, sample_id, compartment, sep = "|")]
group_dt <- merge(group_dt, group_counts[, .(group_id, n_nuclei)], by = "group_id", all.x = TRUE)

global_mean <- Matrix::rowMeans(expr)
global_detect <- Matrix::rowMeans(expr > 0)
eligible <- data.table(gene = features, mean_expression_bin = assign_bins(global_mean), detection_fraction_bin = assign_bins(global_detect),
                       detected = Matrix::rowSums(expr > 0) > 0)
eligible <- eligible[detected == TRUE]

make_random_sets <- function(observed_features) {
  target_exclude <- unique(observed_features)
  obs_info <- eligible[match(observed_features, gene)]
  sets <- matrix(NA_character_, nrow = length(observed_features), ncol = n_random_requested)
  status <- matrix(NA_character_, nrow = length(observed_features), ncol = n_random_requested)
  for (i in seq_along(observed_features)) {
    mb <- obs_info$mean_expression_bin[i]
    db <- obs_info$detection_fraction_bin[i]
    exact <- eligible[mean_expression_bin == mb & detection_fraction_bin == db & !(gene %in% target_exclude), gene]
    near <- eligible[abs(mean_expression_bin - mb) <= 1 & abs(detection_fraction_bin - db) <= 1 & !(gene %in% target_exclude), gene]
    if (length(exact) > 0) {
      pool <- exact; st <- "exact_bin"
    } else if (length(near) > 0) {
      pool <- near; st <- "nearest_bin"
    } else {
      pool <- eligible[!(gene %in% target_exclude), gene]; st <- "unmatched_fallback"
    }
    sets[i, ] <- sample(pool, n_random_requested, replace = length(pool) < n_random_requested)
    status[i, ] <- st
  }
  for (rid in seq_len(n_random_requested)) {
    dup <- which(duplicated(sets[, rid]))
    if (length(dup) > 0) {
      for (i in dup) {
        cand <- setdiff(eligible[!(gene %in% target_exclude), gene], sets[, rid])
        if (length(cand) == 0) cand <- eligible[!(gene %in% target_exclude), gene]
        sets[i, rid] <- sample(cand, 1)
        status[i, rid] <- paste0(status[i, rid], "_duplicate_repaired")
      }
    }
  }
  list(sets = sets, status = status)
}
module_group_scores <- function(gene_sets) {
  out <- matrix(NA_real_, nrow = ncol(gene_sets), ncol = ncol(gene_group_mean),
                dimnames = list(sprintf("R%04d", seq_len(ncol(gene_sets))), colnames(gene_group_mean)))
  idx_mat <- matrix(match(gene_sets, rownames(gene_group_mean)), nrow = nrow(gene_sets), ncol = ncol(gene_sets))
  for (g in colnames(gene_group_mean)) {
    vals <- gene_group_mean[, g]
    out[, g] <- colMeans(matrix(vals[idx_mat], nrow = nrow(gene_sets), ncol = ncol(gene_sets)), na.rm = TRUE)
  }
  out
}
calc_stats_from_score_matrix <- function(score_mat, mn) {
  donors <- sort(unique(group_dt$donor_id))
  contrast <- matrix(NA_real_, nrow = nrow(score_mat), ncol = length(donors), dimnames = list(rownames(score_mat), donors))
  top <- matrix(FALSE, nrow = nrow(score_mat), ncol = length(donors), dimnames = list(rownames(score_mat), donors))
  loop_counts <- integer(length(donors)); names(loop_counts) <- donors
  for (d in donors) {
    loop_gid <- group_dt[donor_id == d & compartment == loop_label, group_id]
    other <- group_dt[donor_id == d & compartment != loop_label]
    weights <- other$n_nuclei / sum(other$n_nuclei)
    loop_score <- score_mat[, loop_gid]
    other_score <- as.numeric(score_mat[, other$group_id, drop = FALSE] %*% weights)
    contrast[, d] <- loop_score - other_score
    top[, d] <- loop_score >= apply(score_mat[, c(loop_gid, other$group_id), drop = FALSE], 1, max, na.rm = TRUE)
    loop_counts[d] <- group_dt[group_id == loop_gid, n_nuclei]
  }
  subsets <- list(full_4_donors = donors, min_loop_tal_ge20 = names(loop_counts)[loop_counts >= 20],
                  min_loop_tal_ge50 = names(loop_counts)[loop_counts >= 50])
  for (d in donors) subsets[[paste0("leave_one_out_exclude_", d)]] <- setdiff(donors, d)
  rbindlist(lapply(names(subsets), function(sn) {
    inc <- subsets[[sn]]
    if (length(inc) == 0) return(NULL)
    data.table(
      module_name = mn,
      random_set_id = rownames(score_mat),
      donor_subset = sn,
      n_donors = length(inc),
      included_donors = collapse_values(inc),
      statistic_name = rep(c("median_loop_tal_minus_other", "mean_loop_tal_minus_other", "n_donors_loop_tal_top_rank"), each = nrow(score_mat)),
      random_value = c(apply(contrast[, inc, drop = FALSE], 1, median, na.rm = TRUE),
                       rowMeans(contrast[, inc, drop = FALSE], na.rm = TRUE),
                       rowSums(top[, inc, drop = FALSE], na.rm = TRUE))
    )
  }))
}

random_summary_list <- list()
random_stats_list <- list()
for (mn in names(driver_removed_modules)) {
  genes <- driver_removed_modules[[mn]]
  rs <- make_random_sets(genes)
  random_summary_list[[mn]] <- data.table(
    module_name = mn,
    n_driver_removed_genes_used = length(genes),
    n_random_sets_requested = n_random_requested,
    n_random_sets_successful = ncol(rs$sets),
    n_failed_matches = sum(grepl("unmatched|reuse", as.vector(rs$status))),
    fallback_matching_used = any(rs$status != "exact_bin"),
    random_seed = random_seed,
    notes = "Driver-removed expression/detection matched random sets."
  )
  random_stats_list[[mn]] <- calc_stats_from_score_matrix(module_group_scores(rs$sets), mn)
}
random_sampling_summary <- rbindlist(random_summary_list)
write_tsv(random_sampling_summary, file.path(table_dir, "phase2_step4_driver_removed_random_sampling_summary.tsv"))
random_stats <- rbindlist(random_stats_list)

observed_for_bench <- driver_removed_stats[donor_subset %in% unique(random_stats$donor_subset),
                                           .(module_name, donor_subset, statistic_name, observed_value = value)]
bench_input <- merge(random_stats, observed_for_bench, by = c("module_name", "donor_subset", "statistic_name"), all.x = TRUE)
bench <- bench_input[, .(
  n_random_sets_successful = uniqueN(random_set_id),
  observed_value = observed_value[1],
  random_median = median(random_value, na.rm = TRUE),
  random_p025 = as.numeric(quantile(random_value, 0.025, na.rm = TRUE)),
  random_p975 = as.numeric(quantile(random_value, 0.975, na.rm = TRUE)),
  observed_minus_random_median = observed_value[1] - median(random_value, na.rm = TRUE),
  observed_percentile = mean(random_value <= observed_value[1], na.rm = TRUE) * 100,
  empirical_p_one_sided = (sum(random_value >= observed_value[1], na.rm = TRUE) + 1) / (.N + 1)
), by = .(module_name, donor_subset, statistic_name)]
bench <- merge(bench, random_sampling_summary[, .(module_name, n_driver_removed_genes_used)], by = "module_name", all.x = TRUE)
bench[, benchmark_support_class := mapply(support_class, observed_percentile, empirical_p_one_sided,
                                          observed_value, random_median, n_driver_removed_genes_used, donor_subset)]
bench[, notes := "Driver-removed matched-random benchmark; sensitivity analysis only."]
setcolorder(bench, c("module_name", "donor_subset", "statistic_name", "n_driver_removed_genes_used",
                     "n_random_sets_successful", "observed_value", "random_median", "random_p025",
                     "random_p975", "observed_minus_random_median", "observed_percentile",
                     "empirical_p_one_sided", "benchmark_support_class", "notes"))
write_tsv(bench, file.path(table_dir, "phase2_step4_driver_removed_matched_random_benchmark_summary.tsv"))

fig_comp <- comparison[statistic_label %in% c("median LoopTAL-minus-other across all 4 donors",
                                              "median LoopTAL-minus-other among donors with Loop/TAL nuclei >=20",
                                              "median LoopTAL-minus-other among donors with Loop/TAL nuclei >=50")]
fig_long <- melt(fig_comp[, .(module_name, statistic_name = statistic_label, original_value, driver_removed_value)],
                 id.vars = c("module_name", "statistic_name"), variable.name = "score_set", value.name = "value")
write_tsv(fig_long, file.path(source_dir, "phase2_step4_original_vs_driver_removed_source_data.tsv"))
p_comp <- ggplot(fig_long, aes(score_set, value, group = module_name, color = module_name)) +
  geom_point(size = 1.8) + geom_line(linewidth = 0.35) +
  facet_wrap(~ statistic_name, scales = "free_y") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), strip.background = element_blank()) +
  labs(x = "", y = "Loop/TAL-minus-other", color = "Module")
ggsave(file.path(figure_dir, "phase2_step4_original_vs_driver_removed_loop_tal_contrast.pdf"), p_comp, width = 8.2, height = 4.8)

dist_src <- bench_input[statistic_name == "median_loop_tal_minus_other" & donor_subset %in% c("full_4_donors", "min_loop_tal_ge20")]
write_tsv_gz(dist_src, file.path(source_dir, "phase2_step4_driver_removed_matched_random_distribution_source_data.tsv.gz"))
p_dist <- ggplot(dist_src, aes(random_value)) +
  geom_histogram(bins = 40, fill = "grey75", color = "white", linewidth = 0.1) +
  geom_vline(aes(xintercept = observed_value), color = "#C44E52", linewidth = 0.5) +
  facet_grid(donor_subset ~ module_name, scales = "free_x") +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Driver-removed random median Loop/TAL-minus-other", y = "Random sets")
ggsave(file.path(figure_dir, "phase2_step4_driver_removed_matched_random_distribution.pdf"), p_dist, width = 10.5, height = 5.2)

write_tsv(donor_scores, file.path(source_dir, "phase2_step4_driver_removed_heatmap_source_data.tsv"))
p_heat <- ggplot(donor_scores, aes(donor_id, compartment, fill = mean_module_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_wrap(~ module_name, ncol = 2, scales = "free_x") +
  scale_fill_viridis_c(option = "C") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank()) +
  labs(x = "Donor", y = "Compartment", fill = "Mean score")
ggsave(file.path(figure_dir, "phase2_step4_driver_removed_heatmap.pdf"), p_heat, width = 8.2, height = 7.8)

primary <- bench[donor_subset %in% c("full_4_donors", "min_loop_tal_ge20", "min_loop_tal_ge50") &
                   statistic_name == "median_loop_tal_minus_other"]
primary_lines <- primary[, paste0("- ", module_name, " / ", donor_subset, ": observed=", signif(observed_value, 4),
                                  ", random median=", signif(random_median, 4),
                                  ", percentile=", signif(observed_percentile, 4),
                                  ", empirical P=", signif(empirical_p_one_sided, 4),
                                  ", class=", benchmark_support_class, ".")]
removed_lines <- module_summary[, paste0("- ", module_name, ": removed ", n_driver_genes_removed, " genes (",
                                         driver_genes_removed, "); remaining present genes=", n_present_in_scrna_after_removal, ".")]
comparison_lines <- comparison[statistic_label == "median LoopTAL-minus-other across all 4 donors",
                               paste0("- ", module_name, ": original=", signif(original_value, 4),
                                      ", driver-removed=", signif(driver_removed_value, 4),
                                      ", change=", signif(percent_change, 4), "%, interpretation=", interpretation, ".")]

writeLines(c(
  "# Phase 2-Step 4 Manuscript-Safe Interpretation Wording",
  "",
  "## Results Wording",
  "",
  "Driver-removal sensitivity indicated that the donor-level Loop/TAL-associated pattern was retained or attenuated after removing curated Loop/TAL, renal transport, calcium-handling and KSD exemplar genes. This sensitivity analysis suggests the signal was not solely explained by the curated driver set, while preserving the primary canonical MAGMA module results.",
  "",
  "Matched-random benchmarks after driver removal were used to contextualize the reduced modules against expression- and detection-matched random gene sets. These analyses are sensitivity checks, not mediation analyses, causal tests or independent validation.",
  "",
  "## Limitations Wording",
  "",
  "The driver list is curated and not exhaustive. Removal is conservative because some removed genes are biologically relevant members of the genetic-priority modules. Retention after driver removal does not prove independence from Loop/TAL biology, and donor number remains limited.",
  "",
  "## Prohibited Wording",
  "",
  "- Do not write: not driven by any known biology.",
  "- Do not write: independent of Loop/TAL genes.",
  "- Do not write: mechanistically validated.",
  "- Do not write: causal after driver removal."
), file.path(note_dir, "phase2_step4_interpretation_wording.md"))

writeLines(c(
  "# Phase 2-Step 4 Report",
  "",
  "## Driver List Used",
  "",
  paste0("- Driver list used: `", sub(paste0(root, "/"), "", driver_used_path, fixed = TRUE), "`."),
  "- Driver genes are sensitivity-analysis inputs only and are not labeled causal.",
  "",
  "## Genes Removed Per Module",
  "",
  removed_lines,
  "",
  "## Original Versus Driver-Removed Statistics",
  "",
  comparison_lines,
  "",
  "## Driver-Removed Matched-Random Benchmark",
  "",
  primary_lines,
  "",
  "## Low-Cell-Count Sensitivity After Driver Removal",
  "",
  "- Low-cell-count sensitivity is represented by full_4_donors, min_loop_tal_ge20, min_loop_tal_ge50 and leave-one-donor-out donor subsets in the benchmark summary.",
  "",
  "## Claim Boundary",
  "",
  "- Safe claim: driver-removal sensitivity tests whether the Loop/TAL-associated pattern is solely explained by a curated driver set.",
  "- Unsafe claim: not driven by any known biology, independent of Loop/TAL genes, mechanistically validated, causal after driver removal.",
  "",
  "## Readiness",
  "",
  "- Step 5 can proceed after human review: snRNA evidence integration and Figure 2 assembly."
), file.path(note_dir, "phase2_step4_report.md"))

writeLines(c(
  "# Phase 2-Step 4 Limitations and Next Steps",
  "",
  "- The curated driver list is useful but not exhaustive.",
  "- Donor number remains small, with 4 total donors.",
  "- GSM7290914 has only 4 Loop/TAL nuclei.",
  "- Driver removal is conservative because removed genes may be biologically relevant MAGMA-prioritized module members.",
  "- Figure 2 component evidence can now be assembled after human review, but this step did not generate final Figure 2.",
  "- Recommended next step: A. proceed to Phase 2-Step 5: snRNA evidence integration and Figure 2 assembly."
), file.path(note_dir, "phase2_step4_limitations_and_next_steps.md"))

checklist <- data.table(
  task_id = sprintf("P2S4-%02d", 1:11),
  task_name = c("Create main driver-removal sensitivity script", "Audit known-driver genes",
                "Create driver-removed module definitions", "Recompute donor-level module scores after driver removal",
                "Compare original versus driver-removed donor-level statistics",
                "Run matched-random benchmark for driver-removed modules",
                "Create driver-removal sensitivity figures", "Create manuscript-safe interpretation wording",
                "Create Step 4 report", "Create limitations and next-steps note",
                "Stop before prohibited downstream analyses"),
  completed = "yes",
  output_file = c(
    "scripts/06_scrna_processing/phase2_step4_known_driver_removal_sensitivity.R",
    "results/tables/phase2_step4_known_driver_gene_audit.tsv",
    "results/phase2_step4_driver_removed_gene_sets/; results/tables/phase2_step4_driver_removed_module_summary.tsv",
    "results/tables/phase2_step4_driver_removed_donor_compartment_module_scores.tsv; results/tables/phase2_step4_driver_removed_loop_tal_vs_other_summary.tsv",
    "results/tables/phase2_step4_original_vs_driver_removed_summary.tsv",
    "results/tables/phase2_step4_driver_removed_random_sampling_summary.tsv; results/tables/phase2_step4_driver_removed_matched_random_benchmark_summary.tsv",
    "results/figures/phase2_step4_original_vs_driver_removed_loop_tal_contrast.pdf; results/figures/phase2_step4_driver_removed_matched_random_distribution.pdf; results/figures/phase2_step4_driver_removed_heatmap.pdf",
    "notes/phase2_step4_interpretation_wording.md",
    "notes/phase2_step4_report.md",
    "notes/phase2_step4_limitations_and_next_steps.md",
    "codex_tasks/phase2_step4_completion_checklist.tsv"
  ),
  blocking_issue = "",
  manual_review_needed = "yes",
  notes = c(rep("Generated within Phase 2-Step 4 boundary.", 10),
            "Final Figure 2, spatial/TWAS/bulk and manuscript edits were not performed.")
)
write_tsv(checklist, file.path(task_dir, "phase2_step4_completion_checklist.tsv"))
