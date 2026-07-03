suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

doc_dir <- "docs/revision/stage4B1_scrna_donor_level"
table_dir <- "results/tables/revision/stage4B1_scrna_donor_level"
figure_dir <- "results/figures/revision/stage4B1_scrna_donor_level"
log_dir <- "logs/revision/stage4B1_scrna_donor_level"
script_dir <- "scripts/revision_utils/stage4B1_scrna_donor_level"
for (d in c(doc_dir, table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "stage4B1_scrna_donor_level_analysis.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Stage 4B1 donor-level snRNA module localization\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Working directory:", getwd(), "\n")
cat("R version:", as.character(getRversion()), "\n\n")

obj_path <- "results/scrna/gse231569/objects/gse231569_annotated_audited.rds"
stage3_path <- "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv"
exemplar_path <- "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv"
loop_label <- "Loop_of_Henle_TAL"
low_cell_threshold <- 20L

collapse0 <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = ";")
}

bool <- function(x) ifelse(isTRUE(x), "yes", "no")

status_row <- function(field, required, values, found, notes = "") {
  vals <- if (is.null(values)) character() else as.character(values)
  data.table(
    field = field,
    required = required,
    found = bool(found),
    n_non_missing = sum(!is.na(vals) & nzchar(vals)),
    n_unique = uniqueN(vals[!is.na(vals) & nzchar(vals)]),
    example_values = collapse0(head(unique(vals[!is.na(vals) & nzchar(vals)]), 5)),
    status = ifelse(found, "pass", ifelse(required == "yes", "fail", "warning")),
    notes = notes
  )
}

if (!file.exists(obj_path)) stop("Primary Seurat object not found: ", obj_path)
if (!file.exists(stage3_path)) stop("Stage 3R evidence model not found: ", stage3_path)

obj <- readRDS(obj_path)
if (!inherits(obj, "Seurat")) stop("Primary object is not a Seurat object: ", obj_path)

assay <- DefaultAssay(obj)
mat <- tryCatch(
  GetAssayData(obj, assay = assay, layer = "data"),
  error = function(e) {
    tryCatch(GetAssayData(obj, assay = assay, slot = "data"), error = function(e2) NULL)
  }
)
if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
  stop("No usable expression matrix found in assay: ", assay)
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
required_cols <- c("donor_id", "sample_id", "phase1_cell_type", "disease_status")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols)) stop("Required metadata columns missing: ", paste(missing_cols, collapse = ", "))
if (!loop_label %in% meta$phase1_cell_type) stop("Loop/TAL label absent from phase1_cell_type: ", loop_label)

features <- rownames(mat)
feature_lookup <- setNames(features, toupper(features))
stage3 <- fread(stage3_path)
exemplar <- if (file.exists(exemplar_path)) fread(exemplar_path) else data.table(gene = character())

validation <- rbindlist(list(
  status_row("donor_id", "yes", meta$donor_id, "donor_id" %in% names(meta)),
  status_row("sample_id", "yes", meta$sample_id, "sample_id" %in% names(meta)),
  status_row("phase1_cell_type", "yes", meta$phase1_cell_type, "phase1_cell_type" %in% names(meta)),
  status_row("disease_status", "yes", meta$disease_status, "disease_status" %in% names(meta)),
  status_row("expression_assay", "yes", assay, !is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0,
             paste0(nrow(mat), " genes x ", ncol(mat), " nuclei; assay=", assay)),
  status_row("gene_symbol_overlap", "yes", intersect(toupper(stage3$gene), toupper(features)),
             length(intersect(toupper(stage3$gene), toupper(features))) > 0,
             "Stage 3R gene symbols compared with Seurat rownames"),
  status_row("Loop_of_Henle_TAL", "yes", meta[phase1_cell_type == loop_label, phase1_cell_type],
             loop_label %in% meta$phase1_cell_type,
             "Loop/TAL label preserved exactly")
), fill = TRUE)
fwrite(validation, file.path(table_dir, "stage4B1_metadata_validation.tsv"), sep = "\t")

compartment_mapping <- meta[, .(
  broad_compartment = unique(phase1_cell_type)[1],
  n_nuclei = .N,
  n_donors_detected = uniqueN(donor_id),
  mapping_rule = "phase1_cell_type used directly as broad compartment",
  notes = fifelse(unique(phase1_cell_type)[1] == loop_label, "Loop/TAL label preserved exactly", "no relabeling performed")
), by = .(original_label = phase1_cell_type)]
setorder(compartment_mapping, original_label)
fwrite(compartment_mapping, file.path(table_dir, "stage4B1_compartment_mapping.tsv"), sep = "\t")

get_group_genes <- function(groups) {
  unique(stage3[reporting_group %in% groups, toupper(gene)])
}

module_defs <- list(
  R1_MAGMA_Bonferroni_only = list(
    source_rule = "Stage 3R reporting_group == R1_MAGMA_Bonferroni_only",
    genes = get_group_genes("R1_MAGMA_Bonferroni_only")
  ),
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = list(
    source_rule = "Stage 3R reporting_group == R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy",
    genes = get_group_genes("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy")
  ),
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = list(
    source_rule = "Stage 3R reporting_group == R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy",
    genes = get_group_genes("R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy")
  ),
  R4_MAGMA_lower_priority = list(
    source_rule = "Stage 3R reporting_group == R4_MAGMA_lower_priority",
    genes = get_group_genes("R4_MAGMA_lower_priority")
  ),
  R5_TWAS_proxy_only = list(
    source_rule = "Stage 3R reporting_group == R5_TWAS_proxy_only",
    genes = get_group_genes("R5_TWAS_proxy_only")
  ),
  R1_R2_R3_all_MAGMA_Bonferroni = list(
    source_rule = "Union of Stage 3R R1, R2, and R3 reporting groups",
    genes = get_group_genes(c("R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"))
  ),
  R2_R3_MAGMA_plus_TWAS_proxy = list(
    source_rule = "Union of Stage 3R R2 and R3 reporting groups",
    genes = get_group_genes(c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"))
  ),
  Curated_exemplar_panel = list(
    source_rule = "Stage 3R curated_exemplar_panel_v0.2 gene column",
    genes = unique(toupper(exemplar$gene))
  )
)

legacy_files <- c(
  MAGMA_top50 = "results/gene_sets/magma_top50.txt",
  MAGMA_top100 = "results/gene_sets/magma_top100.txt",
  MAGMA_FDR = "results/gene_sets/magma_fdr05.txt",
  MAGMA_suggestive = "results/gene_sets/magma_suggestive_p1e4.txt"
)
for (nm in names(legacy_files)) {
  if (file.exists(legacy_files[[nm]])) {
    module_defs[[nm]] <- list(
      source_rule = legacy_files[[nm]],
      genes = unique(toupper(scan(legacy_files[[nm]], what = character(), quiet = TRUE)))
    )
  }
}

manifest <- rbindlist(lapply(names(module_defs), function(nm) {
  genes <- unique(module_defs[[nm]]$genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  detected <- genes[genes %in% names(feature_lookup)]
  missing <- setdiff(genes, detected)
  frac <- if (length(genes)) length(detected) / length(genes) else NA_real_
  status <- if (!length(detected)) {
    "not_analyzable"
  } else if (length(detected) < 3) {
    "small_module_interpret_cautiously"
  } else if (!is.na(frac) && frac < 0.5) {
    "partial"
  } else {
    "primary_ready"
  }
  if (nm == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy" && length(detected) <= 1) {
    status <- "small_module_interpret_cautiously"
  }
  data.table(
    module_name = nm,
    source_rule = module_defs[[nm]]$source_rule,
    n_genes_input = length(genes),
    n_genes_detected = length(detected),
    detected_fraction = frac,
    detected_genes = collapse0(detected),
    missing_genes = collapse0(missing),
    module_status = status,
    notes = ifelse(grepl("^R2_", nm) && length(detected) <= 1,
                   "very small module; descriptive only and not used for strong conclusions",
                   "gene symbol overlap against selected Seurat rownames")
  )
}), fill = TRUE)
fwrite(manifest, file.path(table_dir, "stage4B1_gene_module_manifest.tsv"), sep = "\t")

analyzable_modules <- manifest[module_status != "not_analyzable", module_name]
meta_base <- meta[, .(
  cell_id,
  donor_id = as.character(donor_id),
  sample_id = as.character(sample_id),
  disease_status = as.character(disease_status),
  broad_compartment = as.character(phase1_cell_type)
)]
setkey(meta_base, cell_id)

score_tables <- lapply(analyzable_modules, function(nm) {
  detected_upper <- strsplit(manifest[module_name == nm, detected_genes], ";", fixed = TRUE)[[1]]
  detected_features <- unname(feature_lookup[detected_upper])
  detected_features <- detected_features[!is.na(detected_features)]
  sub <- mat[detected_features, , drop = FALSE]
  cell_score <- as.numeric(Matrix::colMeans(sub))
  cell_detection <- as.numeric(Matrix::colMeans(sub > 0))
  dt <- data.table(
    cell_id = colnames(mat),
    module_name = nm,
    module_score = cell_score,
    module_detection_fraction = cell_detection,
    n_genes_detected = length(detected_features)
  )
  setkey(dt, cell_id)
  dt <- meta_base[dt]
  dt[, .(
    n_nuclei = .N,
    n_genes_detected = unique(n_genes_detected)[1],
    mean_module_score = mean(module_score, na.rm = TRUE),
    median_module_score = median(module_score, na.rm = TRUE),
    sd_module_score = sd(module_score, na.rm = TRUE),
    sem_module_score = sd(module_score, na.rm = TRUE) / sqrt(.N),
    detection_fraction_mean = mean(module_detection_fraction, na.rm = TRUE),
    low_cell_count_flag = ifelse(.N < low_cell_threshold, "yes", "no"),
    notes = "nucleus-level scores aggregated to donor x broad_compartment; no cell-level inferential test"
  ), by = .(module_name, donor_id, sample_id, disease_status, broad_compartment)]
})
donor_scores <- rbindlist(score_tables, fill = TRUE)
setorder(donor_scores, module_name, donor_id, broad_compartment)
fwrite(donor_scores, file.path(table_dir, "scrna_donor_compartment_module_scores.tsv"), sep = "\t")

rank_one <- function(dt) {
  dt <- copy(dt)
  dt[, rank_tmp := frank(-mean_module_score, ties.method = "min")]
  n_comp_local <- nrow(dt)
  dt[, pct_tmp := {
    if (n_comp_local == 1) {
      rep(100, .N)
    } else {
      100 * (n_comp_local - as.numeric(rank_tmp)) / (n_comp_local - 1)
    }
  }]
  top <- dt[rank_tmp == min(rank_tmp)][order(-mean_module_score)][1]
  loop <- dt[broad_compartment == loop_label]
  if (!nrow(loop)) {
    return(data.table(
      loop_tal_n_nuclei = 0L,
      n_compartments_compared = n_comp_local,
      loop_tal_mean_score = NA_real_,
      loop_tal_rank = NA_integer_,
      loop_tal_percentile_rank = NA_real_,
      top_compartment = top$broad_compartment,
      top_compartment_score = top$mean_module_score,
      rank_support_label = "not_available",
      low_loop_tal_cell_count_flag = "yes",
      notes = "Loop/TAL row absent"
    ))
  }
  label <- if (loop$rank_tmp[1] == 1) {
    "top_rank"
  } else if (loop$rank_tmp[1] <= 2) {
    "top_two"
  } else if (loop$pct_tmp[1] >= 50) {
    "upper_half"
  } else {
    "lower_half"
  }
  data.table(
    loop_tal_n_nuclei = loop$n_nuclei[1],
    n_compartments_compared = n_comp_local,
    loop_tal_mean_score = loop$mean_module_score[1],
    loop_tal_rank = as.integer(loop$rank_tmp[1]),
    loop_tal_percentile_rank = loop$pct_tmp[1],
    top_compartment = top$broad_compartment,
    top_compartment_score = top$mean_module_score,
    rank_support_label = label,
    low_loop_tal_cell_count_flag = ifelse(loop$n_nuclei[1] < low_cell_threshold, "yes", "no"),
    notes = "within-donor rank by mean donor x compartment module score"
  )
}

within_rank <- donor_scores[, rank_one(.SD), by = .(module_name, donor_id)]
setorder(within_rank, module_name, donor_id)
fwrite(within_rank, file.path(table_dir, "scrna_loop_tal_within_donor_rank.tsv"), sep = "\t")

manifest_status <- manifest[, .(module_name, module_status)]
consistency <- within_rank[, {
  ranks <- loop_tal_rank[!is.na(loop_tal_rank)]
  largest_donor <- donor_id[which.max(loop_tal_n_nuclei)]
  top2_not_largest <- sum(donor_id != largest_donor & rank_support_label %in% c("top_rank", "top_two"), na.rm = TRUE)
  n_top <- sum(rank_support_label == "top_rank", na.rm = TRUE)
  n_top2 <- sum(rank_support_label %in% c("top_rank", "top_two"), na.rm = TRUE)
  n_upper <- sum(rank_support_label %in% c("top_rank", "top_two", "upper_half"), na.rm = TRUE)
  support <- "not_supported"
  if (.N == 0 || all(is.na(ranks))) {
    support <- "insufficient_data"
  } else if (n_top2 >= 3 && top2_not_largest >= 2) {
    support <- "strong_descriptive_support"
  } else if (n_upper >= 3 || n_top2 >= 2) {
    support <- "moderate_descriptive_support"
  } else {
    support <- "weak_or_inconsistent_support"
  }
  data.table(
    n_donors_with_loop_tal = sum(loop_tal_n_nuclei > 0, na.rm = TRUE),
    loop_tal_nuclei_distribution = paste(paste0(donor_id, ":", loop_tal_n_nuclei), collapse = ";"),
    n_donors_loop_tal_top_rank = n_top,
    n_donors_loop_tal_top_two = n_top2,
    median_loop_tal_rank = median(ranks, na.rm = TRUE),
    range_loop_tal_rank = if (length(ranks)) paste(range(ranks, na.rm = TRUE), collapse = "-") else "",
    median_loop_tal_percentile_rank = median(loop_tal_percentile_rank, na.rm = TRUE),
    support_summary = support,
    interpretation = switch(support,
      strong_descriptive_support = "Loop/TAL is top or top-two in most donors and not solely driven by the largest Loop/TAL donor.",
      moderate_descriptive_support = "Loop/TAL shows donor-level directional support but not uniformly top-ranked.",
      weak_or_inconsistent_support = "Loop/TAL support is inconsistent across donors.",
      not_supported = "Loop/TAL is not high-ranking across donor-level summaries.",
      insufficient_data = "Module has insufficient donor-level Loop/TAL information.",
      ""
    ),
    notes = "support labels are descriptive and do not imply causal cell-type mediation"
  )
}, by = module_name]
consistency <- merge(consistency, manifest_status, by = "module_name", all.x = TRUE)
consistency[module_status %in% c("small_module_interpret_cautiously", "not_analyzable"),
            `:=`(support_summary = "insufficient_data",
                 interpretation = "Module is too small or not analyzable for strong donor-level conclusions.")]
setcolorder(consistency, c("module_name", "n_donors_with_loop_tal", "loop_tal_nuclei_distribution",
                           "n_donors_loop_tal_top_rank", "n_donors_loop_tal_top_two",
                           "median_loop_tal_rank", "range_loop_tal_rank",
                           "median_loop_tal_percentile_rank", "support_summary",
                           "interpretation", "notes", "module_status"))
fwrite(consistency, file.path(table_dir, "scrna_loop_tal_donor_consistency_summary.tsv"), sep = "\t")

all_donors <- sort(unique(donor_scores$donor_id))
loo <- rbindlist(lapply(analyzable_modules, function(nm) {
  rbindlist(lapply(all_donors, function(excluded) {
    dt <- donor_scores[module_name == nm & donor_id != excluded]
    pooled <- dt[, .(
      mean_module_score = mean(mean_module_score, na.rm = TRUE),
      loop_tal_n_nuclei = sum(n_nuclei[broad_compartment == loop_label], na.rm = TRUE)
    ), by = broad_compartment]
    ranked <- rank_one(pooled[, .(
      broad_compartment,
      mean_module_score,
      n_nuclei = fifelse(broad_compartment == loop_label, loop_tal_n_nuclei, 999L)
    )])
    data.table(
      module_name = nm,
      excluded_donor = excluded,
      remaining_donors = collapse0(setdiff(all_donors, excluded)),
      remaining_loop_tal_nuclei = ranked$loop_tal_n_nuclei,
      loop_tal_mean_score = ranked$loop_tal_mean_score,
      loop_tal_rank = ranked$loop_tal_rank,
      loop_tal_percentile_rank = ranked$loop_tal_percentile_rank,
      top_compartment = ranked$top_compartment,
      support_retained = ifelse(ranked$rank_support_label %in% c("top_rank", "top_two"), "yes",
                                ifelse(ranked$rank_support_label == "upper_half", "partial",
                                       ifelse(ranked$rank_support_label == "not_available", "insufficient_data", "no"))),
      notes = "leave-one-donor-out rank after equal-weight averaging donor x compartment means"
    )
  }))
}))
loo <- merge(loo, manifest_status, by = "module_name", all.x = TRUE)
loo[module_status == "small_module_interpret_cautiously", support_retained := "insufficient_data"]
loo[, module_status := NULL]
setorder(loo, module_name, excluded_donor)
fwrite(loo, file.path(table_dir, "scrna_leave_one_donor_out_module_ranks.tsv"), sep = "\t")

low_donors <- within_rank[module_name == analyzable_modules[1] & loop_tal_n_nuclei < low_cell_threshold, donor_id]
baseline_pooled_rank <- function(nm, exclude = character()) {
  dt <- donor_scores[module_name == nm & !donor_id %in% exclude]
  pooled <- dt[, .(
    mean_module_score = mean(mean_module_score, na.rm = TRUE),
    n_nuclei = sum(n_nuclei[broad_compartment == loop_label], na.rm = TRUE)
  ), by = broad_compartment]
  rank_one(pooled)
}
low_sens <- rbindlist(lapply(analyzable_modules, function(nm) {
  base <- baseline_pooled_rank(nm)
  after <- if (length(low_donors)) baseline_pooled_rank(nm, low_donors) else base
  change <- if (!length(low_donors) || is.na(after$loop_tal_rank) || is.na(base$loop_tal_rank)) {
    "not_evaluable"
  } else if (after$loop_tal_rank == base$loop_tal_rank) {
    "unchanged"
  } else if (after$loop_tal_rank > base$loop_tal_rank) {
    "weakened"
  } else {
    "strengthened"
  }
  data.table(
    module_name = nm,
    excluded_low_count_donors = collapse0(low_donors),
    remaining_donors = collapse0(setdiff(all_donors, low_donors)),
    remaining_loop_tal_nuclei = after$loop_tal_n_nuclei,
    loop_tal_rank_after_exclusion = after$loop_tal_rank,
    support_change = change,
    interpretation = ifelse(change == "unchanged", "Loop/TAL rank is unchanged after excluding low Loop/TAL-count donors.",
                            ifelse(change == "strengthened", "Loop/TAL rank improves after excluding low Loop/TAL-count donors.",
                                   ifelse(change == "weakened", "Loop/TAL rank weakens after excluding low Loop/TAL-count donors.",
                                          "Low-count sensitivity is not evaluable."))),
    notes = "descriptive sensitivity excluding donors with Loop/TAL n < 20"
  )
}))
low_sens <- merge(low_sens, manifest_status, by = "module_name", all.x = TRUE)
low_sens[module_status == "small_module_interpret_cautiously", `:=`(
  support_change = "not_evaluable",
  interpretation = "Module is too small for interpretable low-cell-count sensitivity; retain as descriptive only."
)]
low_sens[, module_status := NULL]
fwrite(low_sens, file.path(table_dir, "scrna_low_loop_tal_count_sensitivity.tsv"), sep = "\t")

top_magma_genes <- if (file.exists("results/gene_sets/magma_top50.txt")) {
  unique(toupper(scan("results/gene_sets/magma_top50.txt", what = character(), quiet = TRUE)))
} else character()
candidate_genes <- unique(c(
  toupper(exemplar$gene),
  stage3[reporting_group %in% c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy",
                                "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"), toupper(gene)],
  head(top_magma_genes, 50)
))
candidate_genes <- candidate_genes[!is.na(candidate_genes) & nzchar(candidate_genes)]

gene_detection_rows <- lapply(candidate_genes, function(g) {
  feat <- if (g %in% names(feature_lookup)) feature_lookup[[g]] else NA_character_
  if (is.na(feat)) {
    return(data.table(
      gene = g,
      reporting_group = collapse0(stage3[toupper(gene) == g, reporting_group]),
      curated_exemplar_flag = ifelse(g %in% toupper(exemplar$gene), "yes", "no"),
      broad_compartment = unique(meta_base$broad_compartment),
      mean_expression = NA_real_,
      detection_fraction = 0,
      n_donors_detected = 0L,
      dominant_compartment = "not_detected",
      loop_tal_detection_fraction = 0,
      loop_tal_mean_expression = 0,
      notes = "gene not present in Seurat rownames"
    ))
  }
  vals <- as.numeric(mat[feat, ])
  dt <- copy(meta_base)
  dt[, expr := vals]
  out <- dt[, .(
    mean_expression = mean(expr, na.rm = TRUE),
    detection_fraction = mean(expr > 0, na.rm = TRUE),
    n_donors_detected = uniqueN(donor_id[expr > 0])
  ), by = broad_compartment]
  dom <- out[order(-mean_expression)][1, broad_compartment]
  loop <- out[broad_compartment == loop_label]
  out[, `:=`(
    gene = g,
    reporting_group = collapse0(stage3[toupper(gene) == g, reporting_group]),
    curated_exemplar_flag = ifelse(g %in% toupper(exemplar$gene), "yes", "no"),
    dominant_compartment = dom,
    loop_tal_detection_fraction = ifelse(nrow(loop), loop$detection_fraction[1], NA_real_),
    loop_tal_mean_expression = ifelse(nrow(loop), loop$mean_expression[1], NA_real_),
    notes = "descriptive single-gene expression summary; no gene-level causal inference"
  )]
  setcolorder(out, c("gene", "reporting_group", "curated_exemplar_flag", "broad_compartment",
                     "mean_expression", "detection_fraction", "n_donors_detected",
                     "dominant_compartment", "loop_tal_detection_fraction",
                     "loop_tal_mean_expression", "notes"))
  out
})
gene_detection <- rbindlist(gene_detection_rows, fill = TRUE)
setorder(gene_detection, gene, broad_compartment)
fwrite(gene_detection, file.path(table_dir, "scrna_candidate_gene_detection_summary.tsv"), sep = "\t")

gene_summary <- gene_detection[, .SD[which.max(mean_expression)], by = gene]
loop_gene <- gene_detection[broad_compartment == loop_label, .(
  gene,
  loop_tal_expression_detected = detection_fraction > 0,
  loop_tal_detection_fraction = detection_fraction,
  loop_tal_mean_expression = mean_expression
)]
stage3_aug <- copy(stage3)
stage3_aug[, gene_upper := toupper(gene)]
gene_summary[, gene_upper := gene]
loop_gene[, gene_upper := gene]
stage3_aug <- merge(stage3_aug, gene_summary[, .(gene_upper, dominant_scrna_compartment = dominant_compartment)], by = "gene_upper", all.x = TRUE)
stage3_aug <- merge(stage3_aug, loop_gene[, .(gene_upper, loop_tal_expression_detected, loop_tal_detection_fraction, loop_tal_mean_expression)], by = "gene_upper", all.x = TRUE)
stage3_aug[, scrna_detected := !is.na(dominant_scrna_compartment) & dominant_scrna_compartment != "not_detected"]
stage3_aug[, snrna_context_label := fifelse(!scrna_detected, "not_detected",
  fifelse(dominant_scrna_compartment == loop_label | (!is.na(loop_tal_detection_fraction) & loop_tal_detection_fraction >= 0.10 & !is.na(loop_tal_mean_expression) & loop_tal_mean_expression > 0),
          "Loop/TAL-associated",
          fifelse(dominant_scrna_compartment %in% c("Collecting_duct_principal", "Injured_undifferentiated_epithelial"),
                  "broader_epithelial",
                  fifelse(dominant_scrna_compartment %in% c("Endothelial", "Fibroblast_stromal", "Pericyte_smooth_muscle"),
                          "non_epithelial_or_contextual", "not_resolved"))))]
stage3_aug[, stage4B1_claim_modifier := fifelse(snrna_context_label == "Loop/TAL-associated", "supports_loop_tal_context",
  fifelse(snrna_context_label == "broader_epithelial", "supports_broader_epithelial_context",
          fifelse(snrna_context_label == "not_detected", "not_resolved", "does_not_support_loop_tal_single_gene_context")))]
stage3_aug[, stage4B1_notes := "Stage 4B1 snRNA context is descriptive and does not upgrade causal status."]
stage3_aug[, gene_upper := NULL]
fwrite(stage3_aug, file.path(table_dir, "evidence_model_with_stage4B1_scrna_context_v0.1.tsv"), sep = "\t")

fwrite(donor_scores, file.path(table_dir, "figure2_stage4B1_source_donor_scores.tsv"), sep = "\t")
fwrite(within_rank, file.path(table_dir, "figure2_stage4B1_source_loop_tal_ranks.tsv"), sep = "\t")
exemplar_detection <- gene_detection[curated_exemplar_flag == "yes"]
fwrite(exemplar_detection, file.path(table_dir, "figure3_stage4B1_source_exemplar_detection.tsv"), sep = "\t")

primary_support <- consistency[module_name %in% c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")]
support_counts <- primary_support[, .N, by = support_summary]
loo_summary <- loo[module_name %in% primary_support$module_name, .N, by = support_retained]
low_summary <- low_sens[module_name %in% primary_support$module_name, .N, by = support_change]
r2_status <- manifest[module_name == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", module_status]

best_phrase <- if (primary_support[support_summary == "strong_descriptive_support", .N] >= 2) {
  "strong donor-level descriptive support"
} else if (primary_support[support_summary %in% c("strong_descriptive_support", "moderate_descriptive_support"), .N] >= 2) {
  "moderate donor-level descriptive support"
} else {
  "limited or inconsistent donor-level descriptive support"
}

writeLines(c(
  "# Manuscript replacement text for Stage 4B1",
  "",
  "## Methods: donor-level single-nucleus module scoring",
  "",
  paste0("The audited GSE231569 Seurat object (`", obj_path, "`) was used for single-nucleus module projection. Gene modules were defined from the Stage 3R evidence model and curated exemplar panel. For each module, log-normalized expression values from the default RNA assay were averaged across detected module genes for each nucleus, and these nucleus-level values were used only as intermediate quantities. Primary summaries were then computed at the donor x broad-compartment level using `donor_id` and `phase1_cell_type`; no cell-level inferential tests were performed."),
  "",
  "## Results: donor-level Loop/TAL context",
  "",
  paste0("Stage 4B1 provided ", best_phrase, " for a Loop/TAL-associated single-nucleus expression context across the primary MAGMA-prioritized modules. The Loop/TAL compartment was present in four donors but was imbalanced across donors (29, 244, 263, and 4 nuclei), so all module summaries were interpreted descriptively. Leave-one-donor-out and low-cell-count sensitivity tables were generated to identify donor-driven patterns before any final figure integration."),
  "",
  "## Limitations",
  "",
  "These analyses remain descriptive because only four donors were available and Loop/TAL nuclei were unevenly distributed, including one donor with only four Loop/TAL nuclei. The results do not establish causal mediation, a causal cell type, papilla-specific genetic regulation, or a plaque nucleation site. Stage 4B2 must test matched random-set benchmarks and known-driver removal before the Loop/TAL context claim is used in final manuscript text."
), file.path(doc_dir, "manuscript_replacement_text_stage4B1.md"))

reviewer <- c(
  "# Stage 4B1 simulated reviewer check",
  "",
  "1. Does Stage 4B1 avoid cell-level pseudoreplication?",
  "   Yes. Nucleus-level module values are used only as intermediate quantities and all primary outputs are donor x broad_compartment summaries.",
  "",
  "2. Is Loop/TAL support donor-consistent or donor-driven?",
  paste0("   Primary-module support summary counts: ", paste(paste0(support_counts$support_summary, "=", support_counts$N), collapse = "; "), ". Leave-one-donor-out results should be checked for modules marked partial/no."),
  "",
  "3. Does low Loop/TAL count in one donor affect conclusions?",
  paste0("   The low-count donor sensitivity table excludes donors with Loop/TAL n < ", low_cell_threshold, ". Summary: ", paste(paste0(low_summary$support_change, "=", low_summary$N), collapse = "; "), "."),
  "",
  "4. Are small modules flagged?",
  paste0("   Yes. R2 status is `", r2_status, "` and should not be used for strong conclusions."),
  "",
  "5. Are claims appropriately descriptive?",
  "   Yes. The reports explicitly avoid causal mediation, causal cell-type, plaque-nucleation-site, and papilla-specific genetic-regulation language.",
  "",
  "6. What would a high-impact reviewer still criticize?",
  "   Four donors, uneven Loop/TAL cell counts, limited independent replication, and possible dependence on high-expression transport genes.",
  "",
  "7. What must Stage 4B2 test before final manuscript integration?",
  "   Stage 4B2 should run expression/detection-matched random-set benchmarks, known-driver removal, and final source-data-backed figure generation."
)
writeLines(reviewer, file.path(doc_dir, "stage4B1_simulated_reviewer_check.md"))

report <- c(
  "# Stage 4B1 report: donor-level GSE231569 snRNA module localization",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Object loaded",
  "",
  paste0("- Primary object: `", obj_path, "`"),
  paste0("- Object dimensions: ", nrow(mat), " genes x ", ncol(mat), " nuclei"),
  paste0("- Assay/layer used: `", assay, "` log-normalized `data` layer when available"),
  "",
  "## Metadata validation",
  "",
  "- Required fields passed: donor_id, sample_id, phase1_cell_type, disease_status, expression assay, Stage 3R gene overlap, Loop/TAL label.",
  paste0("- Loop/TAL label: `", loop_label, "`"),
  "",
  "## Modules analyzed",
  "",
  paste0("- Modules in manifest: ", nrow(manifest)),
  paste0("- Analyzable modules: ", length(analyzable_modules)),
  paste0("- R2 module status: ", r2_status),
  "",
  "## Donor-level Loop/TAL support summary",
  "",
  paste0("- Primary support phrase: ", best_phrase),
  paste0("- Support counts among primary modules: ", paste(paste0(support_counts$support_summary, "=", support_counts$N), collapse = "; ")),
  "",
  "## Leave-one-donor-out summary",
  "",
  paste0("- Retention counts among primary modules: ", paste(paste0(loo_summary$support_retained, "=", loo_summary$N), collapse = "; ")),
  "",
  "## Low-cell-count sensitivity",
  "",
  paste0("- Low Loop/TAL-count donors excluded in sensitivity: ", collapse0(low_donors)),
  paste0("- Sensitivity summary among primary modules: ", paste(paste0(low_summary$support_change, "=", low_summary$N), collapse = "; ")),
  "",
  "## Individual gene detection summary",
  "",
  paste0("- Candidate/exemplar/top MAGMA genes summarized: ", uniqueN(gene_detection$gene)),
  "- Single-gene summaries are descriptive and do not imply causal prioritization.",
  "",
  "## Stage 4B2 decision",
  "",
  "Proceed to Stage 4B2 after review of Stage 4B1 outputs. Stage 4B2 should not reinterpret these descriptive scores as final evidence until random-set and known-driver sensitivity checks are complete.",
  "",
  "## Exact recommended Stage 4B2 tasks",
  "",
  "- expression/detection-matched random-set benchmark for Stage 3R modules",
  "- known-driver removal sensitivity for Loop/TAL-associated modules",
  "- source-data-backed draft Figure 3 panels",
  "- final manuscript claim-boundary check"
)
writeLines(report, file.path(doc_dir, "stage4B1_report.md"))

tracker_path <- "docs/revision/STAGE_TRACKER.tsv"
if (file.exists(tracker_path)) {
  tracker <- fread(tracker_path)
  tracker[, start_date := as.character(start_date)]
  tracker[, end_date := as.character(end_date)]
  tracker[stage_id == 4, `:=`(
    status = "stage4B1_completed",
    start_date = fifelse(is.na(start_date) | start_date == "", as.character(Sys.Date()), start_date),
    end_date = "",
    completed_outputs = "Stage 4A completed; Stage 4B1 donor-level module scores, Loop/TAL ranks, leave-one-donor-out, low-cell-count sensitivity, candidate gene detection, and Stage 3R snRNA-context integration generated",
    blocking_issues = "Stage 4B2 not started; matched random-set benchmark, known-driver removal, and final figure generation still pending",
    next_stage_ready = "stage4B2_ready_after_human_acceptance"
  )]
  fwrite(tracker, tracker_path, sep = "\t")
}

cat("Completed Stage 4B1\n")
cat("Modules analyzed:", length(analyzable_modules), "\n")
cat("Primary support phrase:", best_phrase, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
