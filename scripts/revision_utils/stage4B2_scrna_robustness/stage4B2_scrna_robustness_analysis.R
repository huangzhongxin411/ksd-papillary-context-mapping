suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

doc_dir <- "docs/revision/stage4B2_scrna_robustness"
table_dir <- "results/tables/revision/stage4B2_scrna_robustness"
figure_dir <- "results/figures/revision/stage4B2_scrna_robustness"
log_dir <- "logs/revision/stage4B2_scrna_robustness"
script_dir <- "scripts/revision_utils/stage4B2_scrna_robustness"
for (d in c(doc_dir, table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "stage4B2_scrna_robustness_analysis.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4B2 matched random benchmark and driver-removal robustness\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", getwd(), "\n")
cat("R version:", as.character(getRversion()), "\n\n")

set.seed(20260626)
random_seed <- 20260626L
n_random_default <- 1000L
loop_label <- "Loop_of_Henle_TAL"

obj_path <- "results/scrna/gse231569/objects/gse231569_annotated_audited.rds"
stage3_path <- "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv"
exemplar_path <- "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv"
stage4b1_score_path <- "results/tables/revision/stage4B1_scrna_donor_level/scrna_donor_compartment_module_scores.tsv"
stage4b1_consistency_path <- "results/tables/revision/stage4B1_scrna_donor_level/scrna_loop_tal_donor_consistency_summary.tsv"
stage4b1_manifest_path <- "results/tables/revision/stage4B1_scrna_donor_level/stage4B1_gene_module_manifest.tsv"
stage4b1_evidence_path <- "results/tables/revision/stage4B1_scrna_donor_level/evidence_model_with_stage4B1_scrna_context_v0.1.tsv"

collapse0 <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = ";")
}
bool <- function(x) ifelse(isTRUE(x), "yes", "no")

check_item <- function(input_item, expected, found, path, notes = "") {
  data.table(
    input_item = input_item,
    expected = expected,
    found = found,
    status = ifelse(found %in% c("yes", "available", "readable"), "pass", "fail"),
    file_path = path,
    notes = notes
  )
}

if (!file.exists(obj_path)) stop("Missing Seurat object: ", obj_path)
obj <- readRDS(obj_path)
if (!inherits(obj, "Seurat")) stop("Object is not Seurat: ", obj_path)
assay <- DefaultAssay(obj)
mat <- tryCatch(
  GetAssayData(obj, assay = assay, layer = "data"),
  error = function(e) tryCatch(GetAssayData(obj, assay = assay, slot = "data"), error = function(e2) NULL)
)
if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) stop("Expression matrix unavailable")
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
req_cols <- c("donor_id", "sample_id", "phase1_cell_type", "disease_status")
missing_cols <- setdiff(req_cols, names(meta))
if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
if (!loop_label %in% meta$phase1_cell_type) stop("Loop/TAL label missing: ", loop_label)

for (p in c(stage3_path, exemplar_path, stage4b1_score_path, stage4b1_consistency_path, stage4b1_manifest_path, stage4b1_evidence_path)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

stage3 <- fread(stage3_path)
exemplar <- fread(exemplar_path)
stage4b1_scores <- fread(stage4b1_score_path)
stage4b1_consistency <- fread(stage4b1_consistency_path)
stage4b1_manifest <- fread(stage4b1_manifest_path)
stage4b1_evidence <- fread(stage4b1_evidence_path)

features <- rownames(mat)
feature_lookup <- setNames(features, toupper(features))

input_check <- rbindlist(list(
  check_item("Seurat object", "exists and readable", "readable", obj_path, paste0(nrow(mat), " genes x ", ncol(mat), " nuclei")),
  check_item("donor_id", "metadata column exists", bool("donor_id" %in% names(meta)), obj_path, paste0(uniqueN(meta$donor_id), " donors")),
  check_item("phase1_cell_type", "metadata column exists", bool("phase1_cell_type" %in% names(meta)), obj_path, paste0(uniqueN(meta$phase1_cell_type), " compartments")),
  check_item("Loop_of_Henle_TAL", "label exists", bool(loop_label %in% meta$phase1_cell_type), obj_path, paste0(sum(meta$phase1_cell_type == loop_label), " nuclei")),
  check_item("Stage 4B1 module score table", "exists", bool(file.exists(stage4b1_score_path)), stage4b1_score_path),
  check_item("Stage 4B1 donor consistency table", "exists", bool(file.exists(stage4b1_consistency_path)), stage4b1_consistency_path),
  check_item("Stage 3R evidence model", "exists", bool(file.exists(stage3_path)), stage3_path),
  check_item("Module definitions", "reconstructable from Stage 4B1 manifest", bool(nrow(stage4b1_manifest) > 0), stage4b1_manifest_path, paste0(nrow(stage4b1_manifest), " modules")),
  check_item("Gene symbol overlap", "Stage 3R genes match object rownames", bool(length(intersect(toupper(stage3$gene), toupper(features))) > 0), obj_path, paste0(length(intersect(toupper(stage3$gene), toupper(features))), " overlapping genes"))
), fill = TRUE)
fwrite(input_check, file.path(table_dir, "stage4B2_input_consistency_check.tsv"), sep = "\t")

manifest_genes <- setNames(
  lapply(stage4b1_manifest$detected_genes, function(x) {
    y <- strsplit(as.character(x), ";", fixed = TRUE)[[1]]
    unique(toupper(y[nzchar(y)]))
  }),
  stage4b1_manifest$module_name
)

module_role <- function(nm, n_detected) {
  if (nm %in% c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")) return("primary_main_claim")
  if (nm %in% c("R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy")) return("proxy_TWAS_caution")
  if (nm %in% c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R5_TWAS_proxy_only") || n_detected < 3) return("small_module_not_for_main_claim")
  if (nm == "Curated_exemplar_panel") return("secondary_supporting")
  if (grepl("^MAGMA_", nm)) return("legacy_comparison")
  "secondary_supporting"
}

module_priority <- merge(
  stage4b1_manifest[, .(module_name, n_genes_input, n_genes_detected, module_status)],
  stage4b1_consistency[, .(module_name, stage4B1_support_summary = support_summary)],
  by = "module_name",
  all.x = TRUE
)
module_priority[, module_role := mapply(module_role, module_name, n_genes_detected)]
module_priority[, stage4B2_required_tests := fifelse(
  module_role == "small_module_not_for_main_claim",
  "flag_only_no_main_benchmark",
  "matched_random_benchmark;known_driver_removal_if_primary_or_claim_relevant"
)]
module_priority[, interpretation_priority := fifelse(module_role == "primary_main_claim", "main",
  fifelse(module_role %in% c("secondary_supporting", "legacy_comparison"), "supporting",
          fifelse(module_role == "proxy_TWAS_caution", "caution", "not_for_main_claim")))]
module_priority[, notes := fifelse(module_role == "small_module_not_for_main_claim", "Small module; do not use for strong Loop/TAL claim.", "Stage 4B2 robustness module priority.")]
setcolorder(module_priority, c("module_name", "module_role", "n_genes_input", "n_genes_detected", "stage4B1_support_summary", "stage4B2_required_tests", "interpretation_priority", "notes", "module_status"))
fwrite(module_priority, file.path(table_dir, "stage4B2_module_priority_table.tsv"), sep = "\t")

cell_group <- paste(meta$donor_id, meta$sample_id, meta$disease_status, meta$phase1_cell_type, sep = "||")
group_levels <- unique(cell_group)
group_meta <- tstrsplit(group_levels, "\\|\\|")
group_dt <- data.table(
  group_id = group_levels,
  donor_id = group_meta[[1]],
  sample_id = group_meta[[2]],
  disease_status = group_meta[[3]],
  broad_compartment = group_meta[[4]]
)
group_factor <- factor(cell_group, levels = group_levels)
group_mat <- sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), length(group_levels)),
  dimnames = list(colnames(mat), group_levels)
)
group_counts <- as.numeric(Matrix::colSums(group_mat))
names(group_counts) <- group_levels
gene_by_group_sum <- mat %*% group_mat
gene_by_group_mean <- t(t(gene_by_group_sum) / group_counts)
gene_by_group_detect <- (mat > 0) %*% group_mat
gene_by_group_detect <- t(t(gene_by_group_detect) / group_counts)

global_mean <- Matrix::rowMeans(mat)
global_detect <- Matrix::rowMeans(mat > 0)
loop_cells <- meta$phase1_cell_type == loop_label
loop_mean <- Matrix::rowMeans(mat[, loop_cells, drop = FALSE])
loop_detect <- Matrix::rowMeans(mat[, loop_cells, drop = FALSE] > 0)
module_genes_all <- unique(unlist(manifest_genes, use.names = FALSE))
universe <- data.table(
  gene = toupper(features),
  feature = features,
  mean_expression_global = as.numeric(global_mean),
  detection_fraction_global = as.numeric(global_detect),
  n_cells_detected = as.integer(global_detect * ncol(mat)),
  mean_expression_loop_tal = as.numeric(loop_mean),
  detection_fraction_loop_tal = as.numeric(loop_detect),
  gene_length_if_available = NA_real_,
  magma_or_module_gene = ifelse(toupper(features) %in% module_genes_all, "yes", "no"),
  eligible_for_random_sampling = ifelse(global_detect > 0 & is.finite(global_mean), "yes", "no"),
  notes = "Gene length and SNP-count annotations unavailable in Stage 4B2; expression/detection matching used."
)
fwrite(universe, file.path(table_dir, "scrna_gene_matching_universe.tsv"), sep = "\t")

eligible_genes <- universe[eligible_for_random_sampling == "yes", gene]
make_bins <- function(x, n = 10) {
  qs <- unique(quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  if (length(qs) < 3) return(rep("all", length(x)))
  as.character(cut(x, breaks = qs, include.lowest = TRUE))
}
universe[, expr_bin := make_bins(mean_expression_global, 10)]
universe[, detect_bin := make_bins(detection_fraction_global, 10)]
bin_key <- setNames(paste(universe$expr_bin, universe$detect_bin, sep = "||"), universe$gene)
expr_key <- setNames(universe$expr_bin, universe$gene)
detect_key <- setNames(universe$detect_bin, universe$gene)
combined_pools <- split(eligible_genes, bin_key[eligible_genes])
expr_pools <- split(eligible_genes, expr_key[eligible_genes])
detect_pools <- split(eligible_genes, detect_key[eligible_genes])

sample_matched_genes <- function(real_genes) {
  real_genes <- unique(real_genes[real_genes %in% eligible_genes])
  sampled <- character()
  fallback <- character()
  need_by_bin <- table(bin_key[real_genes])
  for (bk in names(need_by_bin)) {
    need <- as.integer(need_by_bin[[bk]])
    pool <- setdiff(combined_pools[[bk]], c(real_genes, sampled))
    take <- min(length(pool), need)
    if (take > 0) sampled <- c(sampled, sample(pool, take))
    missing_n <- need - take
    if (missing_n > 0) {
      expr_bin <- strsplit(bk, "\\|\\|")[[1]][1]
      pool <- setdiff(expr_pools[[expr_bin]], c(real_genes, sampled))
      take <- min(length(pool), missing_n)
      if (take > 0) sampled <- c(sampled, sample(pool, take))
      missing_n <- missing_n - take
      fallback <- c(fallback, "expression_only")
    }
    if (missing_n > 0) {
      detect_bin <- strsplit(bk, "\\|\\|")[[1]][2]
      pool <- setdiff(detect_pools[[detect_bin]], c(real_genes, sampled))
      take <- min(length(pool), missing_n)
      if (take > 0) sampled <- c(sampled, sample(pool, take))
      missing_n <- missing_n - take
      fallback <- c(fallback, "detection_only")
    }
    if (missing_n > 0) {
      pool <- setdiff(eligible_genes, c(real_genes, sampled))
      take <- min(length(pool), missing_n)
      if (take > 0) sampled <- c(sampled, sample(pool, take))
      fallback <- c(fallback, "global_pool")
    }
  }
  list(genes = unique(sampled), fallback = collapse0(fallback))
}

score_module_group <- function(genes) {
  genes <- unique(toupper(genes))
  genes <- genes[genes %in% rownames(gene_by_group_mean)]
  if (!length(genes)) return(NULL)
  score <- Matrix::colMeans(gene_by_group_mean[genes, , drop = FALSE])
  detect <- Matrix::colMeans(gene_by_group_detect[genes, , drop = FALSE])
  out <- copy(group_dt)
  out[, `:=`(
    mean_module_score = as.numeric(score[group_id]),
    detection_fraction_mean = as.numeric(detect[group_id]),
    n_genes_detected = length(genes)
  )]
  out
}

rank_loop <- function(score_dt) {
  score_dt <- copy(score_dt)
  score_dt[, rank_tmp := frank(-mean_module_score, ties.method = "min"), by = donor_id]
  score_dt[, n_comp := .N, by = donor_id]
  score_dt[, pct_tmp := fifelse(n_comp == 1, 100, 100 * (n_comp - as.numeric(rank_tmp)) / (n_comp - 1))]
  loop <- score_dt[broad_compartment == loop_label]
  top <- score_dt[order(donor_id, rank_tmp, -mean_module_score), .SD[1], by = donor_id]
  merge(
    loop[, .(donor_id, loop_tal_rank = as.integer(rank_tmp), loop_tal_percentile = pct_tmp, loop_tal_score = mean_module_score)],
    top[, .(donor_id, top_compartment = broad_compartment, top_score = mean_module_score)],
    by = "donor_id",
    all.x = TRUE
  )
}

summary_metrics <- function(score_dt) {
  ranks <- rank_loop(score_dt)
  loop_scores <- score_dt[broad_compartment == loop_label, .(donor_id, loop_score = mean_module_score)]
  other_scores <- score_dt[broad_compartment != loop_label, .(other_mean = mean(mean_module_score, na.rm = TRUE)), by = donor_id]
  delta <- merge(loop_scores, other_scores, by = "donor_id", all.x = TRUE)
  list(
    ranks = ranks,
    mean_rank = mean(ranks$loop_tal_rank, na.rm = TRUE),
    n_top_rank = sum(ranks$loop_tal_rank == 1, na.rm = TRUE),
    n_top_two = sum(ranks$loop_tal_rank <= 2, na.rm = TRUE),
    median_percentile = median(ranks$loop_tal_percentile, na.rm = TRUE),
    delta = mean(delta$loop_score - delta$other_mean, na.rm = TRUE)
  )
}

benchmark_modules <- module_priority[module_role != "small_module_not_for_main_claim", module_name]
design_rows <- list()
summary_rows <- list()
donor_detail_rows <- list()
random_metric_store <- list()

for (nm in benchmark_modules) {
  real_genes <- manifest_genes[[nm]]
  real_genes <- real_genes[real_genes %in% rownames(gene_by_group_mean)]
  n_req <- if (module_priority[module_name == nm, module_role] == "primary_main_claim") n_random_default else n_random_default
  observed_scores <- score_module_group(real_genes)
  obs <- summary_metrics(observed_scores)
  random_metrics <- vector("list", n_req)
  donor_random <- list()
  fallbacks <- character()
  generated <- 0L
  for (i in seq_len(n_req)) {
    smp <- sample_matched_genes(real_genes)
    if (length(smp$genes) < max(3L, floor(0.8 * length(real_genes)))) next
    fallbacks <- c(fallbacks, smp$fallback)
    rscores <- score_module_group(smp$genes)
    if (is.null(rscores)) next
    rm <- summary_metrics(rscores)
    generated <- generated + 1L
    random_metrics[[generated]] <- data.table(
      iter = generated,
      mean_rank = rm$mean_rank,
      n_top_rank = rm$n_top_rank,
      n_top_two = rm$n_top_two,
      median_percentile = rm$median_percentile,
      delta = rm$delta
    )
    donor_random[[generated]] <- rm$ranks[, .(iter = generated, donor_id, loop_tal_rank, loop_tal_percentile)]
  }
  random_metrics <- rbindlist(random_metrics[seq_len(generated)], fill = TRUE)
  donor_random <- rbindlist(donor_random[seq_len(generated)], fill = TRUE)
  status <- if (!generated) "failed_insufficient_pool" else if (any(nzchar(fallbacks))) "complete_expression_detection_matched" else "complete_expression_detection_matched"
  design_rows[[nm]] <- data.table(
    module_name = nm,
    n_module_genes_detected = length(real_genes),
    n_random_sets_requested = n_req,
    n_random_sets_generated = generated,
    matching_variables = "global_mean_expression;global_detection_fraction",
    binning_strategy = "10 quantile bins for global mean expression and global detection fraction",
    fallback_used = ifelse(any(nzchar(fallbacks)), collapse0(fallbacks), "none"),
    random_seed = random_seed,
    benchmark_status = status,
    notes = "Random benchmark computed from gene x donor-compartment means; no cell-level inferential replicate used."
  )
  if (generated) {
    pct_top <- mean(random_metrics$n_top_rank <= obs$n_top_rank, na.rm = TRUE)
    pct_med <- mean(random_metrics$median_percentile <= obs$median_percentile, na.rm = TRUE)
    pct_delta <- mean(random_metrics$delta <= obs$delta, na.rm = TRUE)
    p_top <- (sum(random_metrics$n_top_rank >= obs$n_top_rank, na.rm = TRUE) + 1) / (generated + 1)
    p_med <- (sum(random_metrics$median_percentile >= obs$median_percentile, na.rm = TRUE) + 1) / (generated + 1)
    p_delta <- (sum(random_metrics$delta >= obs$delta, na.rm = TRUE) + 1) / (generated + 1)
    interpretation <- if (p_top <= 0.05 && p_med <= 0.05 && p_delta <= 0.05) {
      "supports_loop_tal_context_beyond_matched_random"
    } else if (sum(c(p_top, p_med, p_delta) <= 0.10) >= 2 || p_delta <= 0.05) {
      "partial_support"
    } else {
      "not_above_matched_random"
    }
    summary_rows[[nm]] <- data.table(
      module_name = nm,
      observed_loop_tal_mean_rank = obs$mean_rank,
      observed_n_donors_top_rank = obs$n_top_rank,
      observed_n_donors_top_two = obs$n_top_two,
      observed_median_loop_tal_percentile = obs$median_percentile,
      observed_loop_tal_delta = obs$delta,
      random_set_n = generated,
      random_percentile_for_top_rank_count = pct_top,
      random_percentile_for_median_percentile = pct_med,
      random_percentile_for_delta = pct_delta,
      empirical_p_top_rank_count = p_top,
      empirical_p_median_percentile = p_med,
      empirical_p_delta = p_delta,
      benchmark_status = status,
      interpretation = interpretation,
      notes = "Empirical P values are random-set benchmark tail proportions, not donor-level biological replicate tests. Rank metrics can saturate when matched random sets also rank Loop/TAL highly; delta support is interpreted as partial support when rank metrics are not selective."
    )
    obs_ranks <- obs$ranks
    donor_detail_rows[[nm]] <- rbindlist(lapply(unique(obs_ranks$donor_id), function(d) {
      rd <- donor_random[donor_id == d]
      obsd <- obs_ranks[donor_id == d]
      rank_q <- quantile(rd$loop_tal_rank, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      pct_q <- quantile(rd$loop_tal_percentile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      support <- if (obsd$loop_tal_percentile > pct_q[[3]]) {
        "above_random"
      } else if (obsd$loop_tal_percentile < pct_q[[1]]) {
        "below_random"
      } else {
        "within_random"
      }
      data.table(
        module_name = nm,
        donor_id = d,
        observed_loop_tal_rank = obsd$loop_tal_rank,
        observed_loop_tal_percentile = obsd$loop_tal_percentile,
        random_median_loop_tal_rank = rank_q[[2]],
        random_95pct_rank_interval = paste(rank_q[[1]], rank_q[[3]], sep = "-"),
        random_median_percentile = pct_q[[2]],
        random_95pct_percentile_interval = paste(pct_q[[1]], pct_q[[3]], sep = "-"),
        donor_level_support_vs_random = support,
        notes = "Donor-level random comparison uses matched random gene sets summarized at donor x compartment level."
      )
    }))
    random_metric_store[[nm]] <- random_metrics
  } else {
    summary_rows[[nm]] <- data.table(
      module_name = nm,
      observed_loop_tal_mean_rank = obs$mean_rank,
      observed_n_donors_top_rank = obs$n_top_rank,
      observed_n_donors_top_two = obs$n_top_two,
      observed_median_loop_tal_percentile = obs$median_percentile,
      observed_loop_tal_delta = obs$delta,
      random_set_n = 0L,
      random_percentile_for_top_rank_count = NA_real_,
      random_percentile_for_median_percentile = NA_real_,
      random_percentile_for_delta = NA_real_,
      empirical_p_top_rank_count = NA_real_,
      empirical_p_median_percentile = NA_real_,
      empirical_p_delta = NA_real_,
      benchmark_status = status,
      interpretation = "not_evaluable",
      notes = "Random-set generation failed."
    )
  }
}

random_design <- rbindlist(design_rows, fill = TRUE)
random_summary <- rbindlist(summary_rows, fill = TRUE)
random_detail <- rbindlist(donor_detail_rows, fill = TRUE)
small_design <- module_priority[module_role == "small_module_not_for_main_claim", .(
  module_name,
  n_module_genes_detected = n_genes_detected,
  n_random_sets_requested = 0L,
  n_random_sets_generated = 0L,
  matching_variables = "not_attempted",
  binning_strategy = "not_attempted",
  fallback_used = "not_attempted",
  random_seed = random_seed,
  benchmark_status = "not_attempted_small_module",
  notes = "Small module not benchmarked for main claim."
)]
random_design <- rbind(random_design, small_design, fill = TRUE)
fwrite(random_design, file.path(table_dir, "random_set_matching_design.tsv"), sep = "\t")
fwrite(random_summary, file.path(table_dir, "scrna_random_set_benchmark_summary.tsv"), sep = "\t")
fwrite(random_detail, file.path(table_dir, "scrna_random_set_benchmark_donor_detail.tsv"), sep = "\t")

primary_modules <- module_priority[module_role == "primary_main_claim", module_name]
claim_modules <- c(primary_modules, "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy", "Curated_exemplar_panel")
claim_modules <- intersect(unique(claim_modules), names(manifest_genes))

curated_genes <- unique(toupper(exemplar$gene))
tal_marker_panel <- c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16")
calcium_panel <- c("CASR", "CLDN14", "CLDN10", "PKD2")
single_drivers <- c("UMOD", "CLDN14", "CASR", "CLDN10", "SLC12A1", "KCNJ1", "CLDN16")
known_transport <- unique(c(curated_genes, tal_marker_panel, calcium_panel, "SLC12A3", "SLC34A1", "TRPV6", "SLC30A10"))

gene_contrib_rows <- lapply(primary_modules, function(nm) {
  genes <- manifest_genes[[nm]]
  genes <- genes[genes %in% rownames(gene_by_group_mean)]
  lapply(genes, function(g) {
    lt <- universe[gene == g]
    other_mean <- universe[gene == g, mean_expression_global]
    if (nrow(lt) && (other_mean - lt$mean_expression_loop_tal) != 0) {
      other_comp_mean <- group_dt[broad_compartment != loop_label, group_id]
      feature <- feature_lookup[[g]]
      gene_other_mean <- mean(as.numeric(gene_by_group_mean[feature, other_comp_mean]), na.rm = TRUE)
    } else {
      gene_other_mean <- NA_real_
    }
    data.table(
      module_name = nm,
      gene = g,
      mean_expression_loop_tal = lt$mean_expression_loop_tal,
      detection_fraction_loop_tal = lt$detection_fraction_loop_tal,
      mean_expression_other_compartments = gene_other_mean,
      loop_tal_specificity_ratio = (lt$mean_expression_loop_tal + 1e-9) / (gene_other_mean + 1e-9),
      curated_exemplar = ifelse(g %in% curated_genes, "yes", "no"),
      known_tal_marker = ifelse(g %in% tal_marker_panel, "yes", "no"),
      known_transport_or_calcium_gene = ifelse(g %in% known_transport, "yes", "no"),
      notes = "Descriptive contribution proxy from Loop/TAL expression and specificity; not causal evidence."
    )
  }) |> rbindlist(fill = TRUE)
}) |> rbindlist(fill = TRUE)
gene_contrib_rows[, approx_contribution_rank := frank(-mean_expression_loop_tal, ties.method = "min"), by = module_name]
setcolorder(gene_contrib_rows, c("module_name", "gene", "mean_expression_loop_tal", "detection_fraction_loop_tal",
                                 "mean_expression_other_compartments", "loop_tal_specificity_ratio",
                                 "approx_contribution_rank", "curated_exemplar", "known_tal_marker",
                                 "known_transport_or_calcium_gene", "notes"))
fwrite(gene_contrib_rows, file.path(table_dir, "scrna_gene_contribution_to_loop_tal_signal.tsv"), sep = "\t")

top_contrib <- gene_contrib_rows[, .SD[order(approx_contribution_rank)][1:min(.N, 10)], by = module_name]
driver_sensitivity_rows <- list()
removal_specs_for_module <- function(nm) {
  top5 <- top_contrib[module_name == nm][order(approx_contribution_rank)][1:min(.N, 5), gene]
  top10 <- top_contrib[module_name == nm][order(approx_contribution_rank)][1:min(.N, 10), gene]
  c(
    setNames(as.list(single_drivers), paste0("without_", single_drivers)),
    list(
      without_curated_exemplar_panel = curated_genes,
      without_TAL_marker_panel = tal_marker_panel,
      without_calcium_ion_panel = calcium_panel,
      without_top5_contributors = top5,
      without_top10_contributors = top10
    )
  )
}

for (nm in primary_modules) {
  base_genes <- manifest_genes[[nm]]
  base_genes <- base_genes[base_genes %in% rownames(gene_by_group_mean)]
  base_metrics <- summary_metrics(score_module_group(base_genes))
  specs <- removal_specs_for_module(nm)
  for (rtype in names(specs)) {
    requested <- unique(toupper(unlist(specs[[rtype]])))
    present <- intersect(requested, base_genes)
    remain <- setdiff(base_genes, present)
    if (length(remain) < 3) {
      driver_sensitivity_rows[[paste(nm, rtype)]] <- data.table(
        base_module = nm,
        removal_type = rtype,
        removed_genes_requested = collapse0(requested),
        removed_genes_present = collapse0(present),
        n_genes_remaining = length(remain),
        loop_tal_rank_before = base_metrics$mean_rank,
        loop_tal_rank_after = NA_real_,
        n_donors_top_rank_before = base_metrics$n_top_rank,
        n_donors_top_rank_after = NA_integer_,
        median_loop_tal_percentile_before = base_metrics$median_percentile,
        median_loop_tal_percentile_after = NA_real_,
        support_retained = "not_evaluable",
        support_change = "not_evaluable",
        interpretation = "not_evaluable",
        notes = "Too few genes remain after removal."
      )
      next
    }
    after_metrics <- summary_metrics(score_module_group(remain))
    retained <- if (after_metrics$n_top_rank >= 3 && after_metrics$median_percentile >= 80) "yes" else if (after_metrics$n_top_two >= 3 || after_metrics$median_percentile >= 50) "partial" else "no"
    rank_change <- after_metrics$mean_rank - base_metrics$mean_rank
    pct_change <- after_metrics$median_percentile - base_metrics$median_percentile
    support_change <- if (retained == "no" && base_metrics$n_top_rank >= 3) {
      "collapsed"
    } else if (rank_change > 0.5 || pct_change < -20) {
      "weakened"
    } else if (rank_change < -0.5 || pct_change > 20) {
      "strengthened"
    } else {
      "unchanged"
    }
    interpretation <- if (support_change %in% c("unchanged", "strengthened") && retained == "yes") {
      "robust_to_driver_removal"
    } else if (retained %in% c("yes", "partial")) {
      "partly_driver_dependent"
    } else {
      "driver_dependent"
    }
    driver_sensitivity_rows[[paste(nm, rtype)]] <- data.table(
      base_module = nm,
      removal_type = rtype,
      removed_genes_requested = collapse0(requested),
      removed_genes_present = collapse0(present),
      n_genes_remaining = length(remain),
      loop_tal_rank_before = base_metrics$mean_rank,
      loop_tal_rank_after = after_metrics$mean_rank,
      n_donors_top_rank_before = base_metrics$n_top_rank,
      n_donors_top_rank_after = after_metrics$n_top_rank,
      median_loop_tal_percentile_before = base_metrics$median_percentile,
      median_loop_tal_percentile_after = after_metrics$median_percentile,
      support_retained = retained,
      support_change = support_change,
      interpretation = interpretation,
      notes = "Driver-removal sensitivity is descriptive and based on donor x compartment module ranks."
    )
  }
}
driver_sensitivity <- rbindlist(driver_sensitivity_rows, fill = TRUE)
fwrite(driver_sensitivity, file.path(table_dir, "scrna_known_driver_removal_sensitivity.tsv"), sep = "\t")

driver_module_summary <- driver_sensitivity[, .(
  any_collapsed = any(support_change == "collapsed", na.rm = TRUE),
  any_weakened = any(support_change %in% c("weakened", "collapsed"), na.rm = TRUE),
  retained_yes_fraction = mean(support_retained == "yes", na.rm = TRUE),
  panel_curated = interpretation[removal_type == "without_curated_exemplar_panel"][1],
  panel_tal = interpretation[removal_type == "without_TAL_marker_panel"][1],
  panel_calcium = interpretation[removal_type == "without_calcium_ion_panel"][1]
), by = base_module]
driver_module_summary[, driver_removal_support := fifelse(any_collapsed, "driver_dependent",
  fifelse(any_weakened | retained_yes_fraction < 0.8, "partly_driver_dependent", "robust"))]

random_support_for_module <- function(nm) {
  x <- random_summary[module_name == nm]
  if (!nrow(x)) return("not_evaluable")
  switch(x$interpretation[1],
    supports_loop_tal_context_beyond_matched_random = "supports_beyond_matched_random",
    partial_support = "partial_support",
    not_above_matched_random = "not_supported",
    "not_evaluable"
  )
}
driver_support_for_module <- function(nm) {
  x <- driver_module_summary[base_module == nm, driver_removal_support]
  if (!length(x) || is.na(x)) "not_evaluable" else x
}

claim_rows <- rbindlist(lapply(c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100",
                                  "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy", "Curated_exemplar_panel"), function(nm) {
  donor_support <- stage4b1_consistency[module_name == nm, support_summary]
  rand_support <- random_support_for_module(nm)
  drv_support <- driver_support_for_module(nm)
  role <- module_priority[module_name == nm, module_role]
  strength <- if (role == "primary_main_claim" && donor_support == "strong_descriptive_support" &&
                  rand_support == "supports_beyond_matched_random" && drv_support == "robust") {
    "strong_main_claim_allowed"
  } else if (role == "primary_main_claim" && donor_support %in% c("strong_descriptive_support", "moderate_descriptive_support") &&
             rand_support %in% c("supports_beyond_matched_random", "partial_support") &&
             drv_support %in% c("robust", "partly_driver_dependent")) {
    "moderate_main_claim_allowed"
  } else if (role %in% c("proxy_TWAS_caution", "secondary_supporting")) {
    "supplementary_only"
  } else {
    "not_for_claim"
  }
  wording <- switch(strength,
    strong_main_claim_allowed = "MAGMA-prioritized modules showed donor-level and robustness-supported localization to a Loop/TAL-associated papillary single-nucleus expression context.",
    moderate_main_claim_allowed = "MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal.",
    supplementary_only = "TWAS-proxy or curated exemplar modules were retained as supplementary interpretive context.",
    not_for_claim = "This module was not used for the main Loop/TAL context claim."
  )
  data.table(
    module_name = nm,
    stage4B1_donor_support = ifelse(length(donor_support), donor_support, "not_available"),
    random_benchmark_support = rand_support,
    driver_removal_support = drv_support,
    overall_claim_strength = strength,
    recommended_manuscript_claim = wording,
    claim_not_allowed = "causal mediation; causal cell type; papilla-specific genetic regulation; plaque nucleation site; cell-level inferential enrichment",
    notes = "Claim decision separates primary MAGMA modules from TWAS-proxy and curated exemplar modules."
  )
}), fill = TRUE)
fwrite(claim_rows, file.path(table_dir, "loop_tal_claim_decision_table.tsv"), sep = "\t")

stage4b2_evidence <- copy(stage4b1_evidence)
main_claim <- claim_rows[module_name %in% primary_modules]
global_random <- if (main_claim[random_benchmark_support == "supports_beyond_matched_random", .N] >= 2) "supports_beyond_matched_random" else if (main_claim[random_benchmark_support %in% c("supports_beyond_matched_random", "partial_support"), .N] >= 2) "partial_support" else "not_supported"
global_driver <- if (main_claim[driver_removal_support == "robust", .N] >= 2) "robust" else if (main_claim[driver_removal_support %in% c("robust", "partly_driver_dependent"), .N] >= 2) "partly_driver_dependent" else "driver_dependent"
stage4b2_evidence[, random_benchmark_support := global_random]
stage4b2_evidence[, driver_removal_relevance := fifelse(toupper(gene) %in% known_transport, "known_transport_or_calcium_or_TAL_marker", "not_flagged_driver")]
stage4b2_evidence[, driver_removal_support := global_driver]
stage4b2_evidence[, loop_tal_context_final_support := fifelse(
  snrna_context_label == "Loop/TAL-associated" & global_random == "supports_beyond_matched_random" & global_driver == "robust",
  "strong_context_support",
  fifelse(snrna_context_label == "Loop/TAL-associated" & global_random %in% c("supports_beyond_matched_random", "partial_support"),
          "moderate_context_support",
          fifelse(snrna_context_label == "Loop/TAL-associated", "weak_context_support", "not_resolved"))
)]
stage4b2_evidence[, stage4B2_claim_modifier := fifelse(loop_tal_context_final_support %in% c("strong_context_support", "moderate_context_support"),
  "supports_robust_loop_tal_context", "do_not_use_for_main_loop_tal_claim")]
stage4b2_evidence[, stage4B2_notes := "Stage 4B2 robustness is module-level context evidence; it does not upgrade genes to causal status."]
fwrite(stage4b2_evidence, file.path(table_dir, "evidence_model_with_stage4B2_scrna_robustness_v0.1.tsv"), sep = "\t")

fwrite(random_summary, file.path(table_dir, "figure2_source_random_benchmark.tsv"), sep = "\t")
fwrite(driver_sensitivity, file.path(table_dir, "figure2_source_driver_removal.tsv"), sep = "\t")
fwrite(claim_rows, file.path(table_dir, "figure2_source_claim_decision.tsv"), sep = "\t")

main_strength <- if (claim_rows[overall_claim_strength == "strong_main_claim_allowed", .N] >= 2) {
  "strong"
} else if (claim_rows[overall_claim_strength %in% c("strong_main_claim_allowed", "moderate_main_claim_allowed"), .N] >= 2) {
  "moderate"
} else {
  "downgraded"
}
main_wording <- if (main_strength == "strong") {
  "MAGMA-prioritized modules showed donor-level and robustness-supported localization to a Loop/TAL-associated papillary single-nucleus expression context."
} else if (main_strength == "moderate") {
  "MAGMA-prioritized modules showed donor-level descriptive Loop/TAL-associated patterns with partial support beyond matched random expectations and robustness to known-driver removal."
} else {
  "Single-nucleus projection provided descriptive Loop/TAL-associated patterns, but robustness analyses did not support a main Loop/TAL convergence claim."
}

writeLines(c(
  "# Manuscript replacement text for Stage 4B2",
  "",
  "## Methods: matched random-set benchmarking and driver-removal sensitivity",
  "",
  paste0("Stage 4B2 used the audited GSE231569 Seurat object (`", obj_path, "`) and Stage 3R/4B1 module definitions. For each analyzable module, genes were matched to random gene sets using binned global expression and detection fraction across the single-nucleus atlas. Random-set module scores were computed from gene x donor-compartment summaries and evaluated using Loop/TAL rank, number of donors with Loop/TAL top rank, median Loop/TAL percentile rank, and Loop/TAL-versus-other-compartment score delta. Known-driver sensitivity removed individual transport/ion-handling genes and panels including curated exemplars, TAL markers, calcium/ion genes, and top module contributors. All summaries were descriptive and used donor x broad-compartment as the primary unit."),
  "",
  "## Results: robustness of the Loop/TAL context",
  "",
  paste0("Stage 4B1 showed donor-level Loop/TAL support for the primary MAGMA-prioritized modules. Stage 4B2 further evaluated whether this pattern exceeded expression/detection-matched random expectation and whether it persisted after known-driver removal. Overall claim decision: ", main_strength, ". Recommended wording: ", main_wording),
  "",
  "## Limitations",
  "",
  "These analyses remain limited by four donors, imbalanced Loop/TAL nuclei, lack of gene-length/SNP-count matching annotations for the random benchmark, and possible dependence on high-expression renal transport genes. The results do not establish causal mediation, a causal cell type, papilla-specific genetic regulation, or a plaque nucleation site."
), file.path(doc_dir, "manuscript_replacement_text_stage4B2.md"))

reviewer <- c(
  "# Stage 4B2 simulated reviewer check",
  "",
  "1. Does the Loop/TAL signal exceed expression/detection-matched random expectation?",
  paste0("   Primary-module random benchmark support: ", paste(claim_rows[module_name %in% primary_modules, paste0(module_name, "=", random_benchmark_support)], collapse = "; "), "."),
  "",
  "2. Is the signal robust to removal of UMOD, CLDN14, CASR, CLDN10 and the curated exemplar panel?",
  paste0("   Primary-module driver-removal support: ", paste(claim_rows[module_name %in% primary_modules, paste0(module_name, "=", driver_removal_support)], collapse = "; "), "."),
  "",
  "3. Is the signal driven by one or two TAL marker genes?",
  "   The driver-removal table tests individual genes, TAL marker panel, curated exemplar panel, calcium/ion panel, and top contributors. Any weakened/collapsed rows should be treated as driver dependence.",
  "",
  "4. Are donor-level summaries still the primary evidence?",
  "   Yes. Random and removal metrics are derived from donor x broad-compartment summaries.",
  "",
  "5. Are TWAS-proxy modules separated from main MAGMA modules?",
  "   Yes. R3 and R2_R3 are marked proxy_TWAS_caution/supplementary rather than main-claim modules.",
  "",
  "6. Does any result still risk overclaiming causality?",
  "   Yes if described as causal mediation, causal cell type, papilla-specific regulation, or plaque nucleation. These phrases remain disallowed.",
  "",
  "7. What exact wording should be used for the main Results claim?",
  paste0("   ", main_wording),
  "",
  "8. What exact wording should be used for the Limitations?",
  "   The analysis is donor-level and descriptive, limited by four donors, imbalanced Loop/TAL nuclei, expression/detection-only random matching, and inability to infer causal cell type.",
  "",
  "9. What should remain supplementary only?",
  "   TWAS-proxy modules, curated exemplar-only analyses, small modules, and gene contribution rankings.",
  "",
  "10. Should Stage 4C figure planning begin?",
  paste0("   ", ifelse(main_strength %in% c("strong", "moderate"), "Yes, with conservative source-data-backed panels.", "Only after downgrading the claim and deciding whether a supplementary figure is more appropriate."))
)
writeLines(reviewer, file.path(doc_dir, "stage4B2_simulated_reviewer_check.md"))

report <- c(
  "# Stage 4B2 report: snRNA robustness benchmark and driver-removal sensitivity",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Random benchmark design and completion",
  "",
  paste0("- Random seed: ", random_seed),
  paste0("- Default random sets per analyzable module: ", n_random_default),
  "- Matching variables: global mean expression and global detection fraction.",
  "- Gene length and SNP-count annotations were unavailable; they were not invented.",
  "",
  "## Benchmark support summary",
  "",
  paste0("- Primary modules with beyond-random support: ", claim_rows[module_name %in% primary_modules & random_benchmark_support == "supports_beyond_matched_random", .N], " / ", length(primary_modules)),
  paste0("- Primary modules with partial random support: ", claim_rows[module_name %in% primary_modules & random_benchmark_support == "partial_support", .N], " / ", length(primary_modules)),
  "",
  "## Known-driver removal summary",
  "",
  paste0("- Primary modules robust to driver removal: ", claim_rows[module_name %in% primary_modules & driver_removal_support == "robust", .N], " / ", length(primary_modules)),
  paste0("- Primary modules partly driver-dependent: ", claim_rows[module_name %in% primary_modules & driver_removal_support == "partly_driver_dependent", .N], " / ", length(primary_modules)),
  paste0("- Primary modules driver-dependent: ", claim_rows[module_name %in% primary_modules & driver_removal_support == "driver_dependent", .N], " / ", length(primary_modules)),
  "",
  "## Gene contribution audit summary",
  "",
  paste0("- Primary-module gene contribution rows: ", nrow(gene_contrib_rows)),
  "- Contribution ranks are descriptive and not causal evidence.",
  "",
  "## Updated Loop/TAL claim decision",
  "",
  paste0("- Manuscript main claim level: ", main_strength),
  paste0("- Recommended wording: ", main_wording),
  "",
  "## Limitations",
  "",
  "- Four donors and imbalanced Loop/TAL nuclei remain limiting.",
  "- Random matching used expression/detection only because gene length and SNP-count annotations were unavailable.",
  "- Known-driver removal is descriptive and module-level.",
  "- No causal cell-type or plaque-nucleation claim is allowed.",
  "",
  "## Stage 4C decision",
  "",
  ifelse(main_strength %in% c("strong", "moderate"),
         "Stage 4C figure planning can begin after human review, using conservative source-data-backed panels.",
         "Stage 4C should begin only after the Loop/TAL claim is downgraded or moved to supplementary context.")
)
writeLines(report, file.path(doc_dir, "stage4B2_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4B2_completed",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 4A, 4B1, and 4B2 completed; matched random-set benchmark, known-driver removal sensitivity, gene contribution audit, claim decision, and manuscript replacement text generated",
    blocking_issues = "Stage 4C final figure planning not started; full Stage 4 not complete until figure/source-data integration is reviewed",
    next_stage_ready = "stage4C_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Completed Stage 4B2\n")
cat("Primary modules:", paste(primary_modules, collapse = ", "), "\n")
cat("Main claim level:", main_strength, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
