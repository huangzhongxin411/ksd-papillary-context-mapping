#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
})

stage_id <- "stage6B1_spatial_context_projection"
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
fig_dir <- file.path("results/figures/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
for (d in c(table_dir, doc_dir, fig_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage6B1_spatial_context_projection.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 6B1 conservative spatial tissue-context projection\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

object_dir <- "data/processed/spatial/gse206306/seurat_objects"
sample_ids <- c("GSM6250307", "GSM6250308", "GSM6250309", "GSM6250310", "GSM7166170")
object_paths <- file.path(object_dir, paste0(sample_ids, "_spatial_qc.rds"))
paths <- list(
  stage3 = "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  exemplar = "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
  stage4_claim = "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
  stage5_claim = "results/tables/revision/stage5C1_gse73680_figure4_draft/gse73680_final_claim_lock_stage5C1.tsv",
  stage6A_feasibility = "results/tables/revision/stage6A_spatial_twas_audit/spatial_module_feasibility_audit.tsv",
  stage6A_primary = "results/tables/revision/stage6A_spatial_twas_audit/spatial_primary_input_audit.tsv",
  stage6A_roi = "docs/revision/stage6A_spatial_twas_audit/spatial_roi_plaque_annotation_audit.md",
  top50 = "results/gene_sets/magma_top50.txt",
  top100 = "results/gene_sets/magma_top100.txt"
)

required_files <- c(object_paths, unlist(paths))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Required Stage 6B1 input missing: ", paste(missing_files, collapse = ", "))

check_cols <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop("Missing columns in ", label, ": ", paste(missing, collapse = ", "))
}

stage3 <- fread(paths$stage3)
exemplar <- fread(paths$exemplar)
stage6A_feasibility <- fread(paths$stage6A_feasibility)
check_cols(stage3, c("gene", "reporting_group"), "Stage 3R evidence model")
check_cols(exemplar, c("gene"), "Stage 3R curated exemplar panel")
check_cols(stage6A_feasibility, c("module_name", "feasible_for_stage6B"), "Stage 6A feasibility table")

read_set <- function(path) unique(toupper(scan(path, what = character(), quiet = TRUE)))
get_group <- function(groups) unique(toupper(stage3[reporting_group %in% groups, gene]))

module_defs <- list(
  R1_MAGMA_Bonferroni_only = get_group("R1_MAGMA_Bonferroni_only"),
  R1_R2_R3_all_MAGMA_Bonferroni = get_group(c("R1_MAGMA_Bonferroni_only", "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy")),
  MAGMA_top50 = read_set(paths$top50),
  MAGMA_top100 = read_set(paths$top100),
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = get_group("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy"),
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = get_group("R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy"),
  R2_R3_MAGMA_plus_TWAS_proxy = get_group(c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy")),
  R5_TWAS_proxy_only = get_group("R5_TWAS_proxy_only"),
  Curated_exemplar_panel = unique(toupper(exemplar$gene)),
  LoopTAL_signature = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2"),
  epithelial_general = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  proximal_tubule = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "AQP1", "ALDOB"),
  collecting_duct = c("AQP2", "AQP3", "SCNN1A", "SCNN1B", "SCNN1G", "KRT19"),
  injury_epithelial = c("HAVCR1", "LCN2", "SPP1", "VCAM1", "KRT8", "KRT18", "VIM"),
  ECM_fibrosis = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "TGFBI", "MMP2"),
  mineralization_remodeling = c("SPP1", "ALPL", "MMP7", "MMP9", "RUNX2", "BGLAP", "POSTN"),
  immune_myeloid = c("PTPRC", "AIF1", "LST1", "CD68", "CD14", "C1QA", "C1QB"),
  endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CLDN5"),
  fibroblast_stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA")
)

primary_modules <- c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")
twas_modules <- c("R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy", "R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy", "R2_R3_MAGMA_plus_TWAS_proxy", "R5_TWAS_proxy_only")
context_signatures <- c("LoopTAL_signature", "epithelial_general", "proximal_tubule", "collecting_duct", "injury_epithelial", "ECM_fibrosis", "mineralization_remodeling", "immune_myeloid", "endothelial", "fibroblast_stromal")
tested_contexts <- c("LoopTAL_signature", "injury_epithelial", "ECM_fibrosis", "mineralization_remodeling", "epithelial_general", "immune_myeloid", "endothelial", "fibroblast_stromal")

module_roles <- c(
  R1_MAGMA_Bonferroni_only = "primary_MAGMA_module",
  R1_R2_R3_all_MAGMA_Bonferroni = "primary_MAGMA_module",
  MAGMA_top50 = "primary_MAGMA_module", MAGMA_top100 = "primary_MAGMA_module",
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy = "secondary_TWAS_proxy_module",
  R2_R3_MAGMA_plus_TWAS_proxy = "secondary_TWAS_proxy_module",
  R5_TWAS_proxy_only = "secondary_TWAS_proxy_module",
  Curated_exemplar_panel = "curated_exemplar_context",
  LoopTAL_signature = "spatial_context_signature", epithelial_general = "spatial_context_signature",
  proximal_tubule = "spatial_context_signature", collecting_duct = "spatial_context_signature",
  injury_epithelial = "injury_or_remodeling_signature", ECM_fibrosis = "injury_or_remodeling_signature",
  mineralization_remodeling = "injury_or_remodeling_signature", immune_myeloid = "spatial_context_signature",
  endothelial = "spatial_context_signature", fibroblast_stromal = "spatial_context_signature"
)

objects <- list()
input_rows <- list()
qc_rows <- list()
section_detected <- list()
spot_rows <- list()

add_input <- function(item, expected, found, status, path, notes) {
  input_rows[[length(input_rows) + 1]] <<- data.table(
    input_item = item, expected = expected, found = as.character(found), status = status,
    file_path = path, notes = notes
  )
}

add_input("five_spatial_RDS_objects", "5 readable QC-passed RDS objects", sum(file.exists(object_paths)), ifelse(all(file.exists(object_paths)), "pass", "fail"), paste(object_paths, collapse = ";"), "Stage 6A recommended primary inputs.")
add_input("Stage3R_module_table", "exists/readable", file.exists(paths$stage3), ifelse(file.exists(paths$stage3), "pass", "fail"), paths$stage3, "Current module definitions, not legacy P1 definitions.")
add_input("Stage6A_module_feasibility", "exists/readable", file.exists(paths$stage6A_feasibility), ifelse(file.exists(paths$stage6A_feasibility), "pass", "fail"), paths$stage6A_feasibility, "R5 TWAS-only is expected to remain not score-feasible.")
add_input("ROI_plaque_annotation", "absent", FALSE, "pass_boundary_locked", paths$stage6A_roi, "No ROI-aware analysis permitted; tissue-context projection only.")

mean_z_score <- function(mat) {
  dense <- as.matrix(mat)
  gene_sd <- apply(dense, 1, sd)
  keep <- is.finite(gene_sd) & gene_sd > 0
  if (!any(keep)) return(rep(NA_real_, ncol(dense)))
  gene_mean <- rowMeans(dense[keep, , drop = FALSE])
  z <- sweep(dense[keep, , drop = FALSE], 1, gene_mean, "-")
  z <- sweep(z, 1, gene_sd[keep], "/")
  colMeans(z, na.rm = TRUE)
}

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  path <- object_paths[i]
  cat("Loading", sid, "\n")
  obj <- tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
  if (!inherits(obj, "Seurat")) stop(path, " is not a Seurat object")
  if (!"Spatial" %in% names(obj@assays)) stop("Spatial assay missing in ", sid)
  if (!all(c("nCount_Spatial", "nFeature_Spatial") %in% colnames(obj@meta.data))) stop("Required QC metadata missing in ", sid)
  if (!length(obj@images)) stop("Tissue image missing in ", sid)
  coords <- tryCatch(as.data.table(Seurat::GetTissueCoordinates(obj)), error = function(e) stop("Coordinates unavailable in ", sid, ": ", conditionMessage(e)))
  if (nrow(coords) != ncol(obj)) stop("Coordinate/spot count mismatch in ", sid)

  add_input(paste0(sid, "_Spatial_assay"), "present", TRUE, "pass", path, "Spatial assay verified.")
  add_input(paste0(sid, "_coordinates"), "one coordinate row per spot", nrow(coords), ifelse(nrow(coords) == ncol(obj), "pass", "fail"), path, "Coordinates verified.")
  add_input(paste0(sid, "_tissue_image"), "present", length(obj@images) > 0, "pass", path, "Tissue image verified.")
  add_input(paste0(sid, "_nCount_Spatial"), "present", "nCount_Spatial" %in% colnames(obj@meta.data), "pass", path, "Complexity covariate verified.")
  add_input(paste0(sid, "_nFeature_Spatial"), "present", "nFeature_Spatial" %in% colnames(obj@meta.data), "pass", path, "Complexity covariate verified.")

  DefaultAssay(obj) <- "Spatial"
  if (!"percent_mt" %in% colnames(obj@meta.data)) {
    obj[["percent_mt"]] <- PercentageFeatureSet(obj, assay = "Spatial", pattern = "^MT-")
    pct_status <- "computed"
  } else {
    pct_status <- "present"
  }
  if (!all(is.finite(obj$percent_mt))) stop("percent_mt contains non-finite values in ", sid)
  add_input(paste0(sid, "_percent_mt"), "present or computable", pct_status, "pass", path, "Computed from Spatial counts when absent.")

  obj <- NormalizeData(obj, assay = "Spatial", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  data_mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"), error = function(e) GetAssayData(obj, assay = "Spatial", slot = "data"))
  counts_mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "counts"), error = function(e) GetAssayData(obj, assay = "Spatial", slot = "counts"))

  feature_names <- rownames(data_mat)
  feature_upper <- toupper(feature_names)
  feature_lookup <- setNames(feature_names[!duplicated(feature_upper)], feature_upper[!duplicated(feature_upper)])
  nonzero <- Matrix::rowSums(counts_mat) > 0
  nonzero_upper <- unique(feature_upper[nonzero])

  meta <- as.data.table(obj@meta.data, keep.rownames = "spot_barcode")
  if ("cell" %in% names(coords)) setnames(coords, "cell", "spot_barcode")
  if (!"spot_barcode" %in% names(coords)) coords[, spot_barcode := rownames(Seurat::GetTissueCoordinates(obj))]
  base <- merge(meta, coords[, .(spot_barcode, x_coordinate = x, y_coordinate = y)], by = "spot_barcode", all.x = TRUE, sort = FALSE)
  if (anyNA(base$x_coordinate) || anyNA(base$y_coordinate)) stop("Coordinate merge produced missing values in ", sid)
  base[, `:=`(
    global_spot_id = paste(sid, spot_barcode, sep = ":"),
    section_id = sid,
    sample_id = sid
  )]
  if (uniqueN(base$global_spot_id) != nrow(base)) stop("Global spot IDs are not unique in ", sid)
  add_input(paste0(sid, "_global_spot_id"), "unique section_id:barcode", uniqueN(base$global_spot_id), "pass", path, "Global spot IDs verified.")

  qc_rows[[sid]] <- data.table(
    dataset_id = "GSE206306", section_id = sid, sample_id = sid,
    n_spots_total = ncol(obj), n_tissue_spots = ncol(obj),
    median_nCount_Spatial = median(obj$nCount_Spatial),
    median_nFeature_Spatial = median(obj$nFeature_Spatial),
    median_percent_mt = median(obj$percent_mt),
    image_available = TRUE, roi_available = FALSE,
    notes = "QC-passed tissue spots; no plaque/lesion ROI; tissue-context projection only."
  )

  for (nm in names(module_defs)) {
    genes <- unique(module_defs[[nm]])
    detected_upper <- intersect(genes, nonzero_upper)
    detected_actual <- unname(feature_lookup[detected_upper])
    detected_actual <- detected_actual[!is.na(detected_actual)]
    section_detected[[length(section_detected) + 1]] <- data.table(
      section_id = sid, module_name = nm, n_genes_detected = length(detected_actual),
      detected_genes = paste(detected_upper, collapse = ";")
    )

    stage6A_status <- stage6A_feasibility[module_name == nm, feasible_for_stage6B][1]
    score_allowed <- !is.na(stage6A_status) && stage6A_status != "no" && length(detected_actual) >= 1
    if (!score_allowed) next

    score <- mean_z_score(data_mat[detected_actual, , drop = FALSE])
    score_map <- data.table(spot_barcode = colnames(data_mat), module_score = score)
    out <- merge(base, score_map, by = "spot_barcode", all.x = TRUE, sort = FALSE)
    if (anyNA(out$module_score)) stop("Module score merge failed for ", sid, " / ", nm)
    spot_rows[[length(spot_rows) + 1]] <- out[, .(
      global_spot_id, section_id, sample_id, spot_barcode,
      x_coordinate, y_coordinate,
      nCount_Spatial, nFeature_Spatial, percent_mt,
      module_name = nm, module_role = unname(module_roles[nm]),
      n_genes_detected = length(detected_actual), module_score,
      score_method = "within-section mean z-scored log-normalized expression across non-zero detected genes",
      notes = "Spot-level score for visualization/descriptive co-distribution only; spot is not a biological replicate."
    )]
  }
  objects[[sid]] <- obj
}

input_check <- rbindlist(input_rows, fill = TRUE)
fwrite(input_check, file.path(table_dir, "stage6B1_input_consistency_check.tsv"), sep = "\t")
section_qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(section_qc, file.path(table_dir, "spatial_section_qc_stage6B1.tsv"), sep = "\t")

section_detected_dt <- rbindlist(section_detected)
manifest <- rbindlist(lapply(names(module_defs), function(nm) {
  genes <- unique(module_defs[[nm]])
  rows <- section_detected_dt[module_name == nm]
  union_detected <- unique(unlist(strsplit(rows$detected_genes[nzchar(rows$detected_genes)], ";", fixed = TRUE)))
  union_detected <- union_detected[nzchar(union_detected)]
  frac <- length(union_detected) / length(genes)
  n_sections <- sum(rows$n_genes_detected > 0)
  status <- if (nm == "R5_TWAS_proxy_only" || stage6A_feasibility[module_name == nm, feasible_for_stage6B][1] == "no") {
    "not_score_feasible"
  } else if (nm == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy" || min(rows$n_genes_detected) < 3) {
    "small_module_interpret_cautiously"
  } else if (nm %in% primary_modules || nm %in% context_signatures) {
    "primary_ready"
  } else {
    "supplementary_ready"
  }
  data.table(
    module_name = nm, module_role = unname(module_roles[nm]),
    n_genes_input = length(genes), n_genes_detected_union = length(union_detected),
    detected_fraction_union = frac, n_sections_detected = n_sections,
    module_status = status,
    notes = ifelse(nm %in% twas_modules, "TWAS-proxy module; supplementary only and never used for main spatial claims.", ifelse(nm == "Curated_exemplar_panel", "Curated context only; no single-gene validation.", "Current Stage 3R/Stage 5 definition; section consistency required for interpretation."))
  )
}))
fwrite(manifest, file.path(table_dir, "spatial_module_manifest_stage6B1.tsv"), sep = "\t")

spot_scores <- rbindlist(spot_rows, fill = TRUE)
if (anyDuplicated(spot_scores[, .(global_spot_id, module_name)])) stop("Duplicate global_spot_id/module_name rows in spot score table")
fwrite(spot_scores, file.path(table_dir, "spatial_spot_module_scores.tsv"), sep = "\t")
cat("Spot score rows:", nrow(spot_scores), "\n")

# Complexity adjustment is performed within each section and module.
residualized <- spot_scores[, {
  d <- copy(.SD)
  d[, `:=`(log_nCount = log1p(nCount_Spatial), log_nFeature = log1p(nFeature_Spatial))]
  use_pct <- all(is.finite(d$percent_mt)) && sd(d$percent_mt) > 0
  covars <- c("log_nCount", "log_nFeature", if (use_pct) "percent_mt")
  keep <- complete.cases(d[, covars, with = FALSE]) & is.finite(d$module_score)
  status <- if (sum(keep) < 20 || sd(d$module_score[keep]) == 0) "failed_insufficient_data" else if (use_pct) "success" else "success_without_percent_mt"
  resid <- rep(NA_real_, nrow(d))
  if (status != "failed_insufficient_data") {
    fit <- lm(reformulate(covars, response = "module_score"), data = d[keep])
    resid[keep] <- residuals(fit)
  }
  data.table(
    global_spot_id = d$global_spot_id,
    section_id = d$section_id,
    sample_id = d$sample_id,
    module_name = d$module_name,
    raw_module_score = d$module_score,
    residualized_module_score = resid,
    nCount_Spatial = d$nCount_Spatial,
    nFeature_Spatial = d$nFeature_Spatial,
    percent_mt = d$percent_mt,
    residual_model_status = status,
    notes = paste0("Within-section residualization against ", paste(covars, collapse = "+"), "; spot-level residual for descriptive projection only.")
  )
}, by = .(section_id, module_name)]
# Remove grouping columns duplicated by the returned table.
for (nm in c("section_id", "module_name")) {
  dup <- which(names(residualized) == nm)
  if (length(dup) > 1) residualized[, (dup[-1]) := NULL]
}
fwrite(residualized, file.path(table_dir, "spatial_spot_module_scores_residualized.tsv"), sep = "\t")

score_summary <- merge(
  spot_scores[, .(
    module_role = module_role[1], n_spots = .N, n_genes_detected = n_genes_detected[1],
    mean_raw_score = mean(module_score), median_raw_score = median(module_score), sd_raw_score = sd(module_score)
  ), by = .(section_id, sample_id, module_name)],
  residualized[, .(
    mean_residualized_score = mean(residualized_module_score, na.rm = TRUE),
    median_residualized_score = median(residualized_module_score, na.rm = TRUE),
    score_status = ifelse(all(residual_model_status %in% c("success", "success_without_percent_mt")), "score_and_adjustment_ready", "adjustment_failed"),
    residual_status = residual_model_status[1]
  ), by = .(section_id, sample_id, module_name)],
  by = c("section_id", "sample_id", "module_name"), all = TRUE
)
score_summary[, notes := "Section-level descriptive summary; section, not spot, is the biological unit."]
score_summary[, residual_status := NULL]
setcolorder(score_summary, c("section_id", "sample_id", "module_name", "module_role", "n_spots", "n_genes_detected", "mean_raw_score", "median_raw_score", "sd_raw_score", "mean_residualized_score", "median_residualized_score", "score_status", "notes"))
fwrite(score_summary, file.path(table_dir, "spatial_section_module_score_summary.tsv"), sep = "\t")

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(list(n = sum(keep), r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  list(n = sum(keep), r = unname(ct$estimate), p = ct$p.value)
}

codist_rows <- list()
for (sid in sample_ids) {
  raw_sid <- spot_scores[section_id == sid]
  res_sid <- residualized[section_id == sid]
  for (mod in primary_modules) {
    mod_raw <- raw_sid[module_name == mod, .(global_spot_id, mod_raw = module_score)]
    mod_res <- res_sid[module_name == mod, .(global_spot_id, mod_res = residualized_module_score)]
    for (ctx in tested_contexts) {
      ctx_raw <- raw_sid[module_name == ctx, .(global_spot_id, ctx_raw = module_score)]
      ctx_res <- res_sid[module_name == ctx, .(global_spot_id, ctx_res = residualized_module_score)]
      d <- Reduce(function(x, y) merge(x, y, by = "global_spot_id", all = FALSE), list(mod_raw, mod_res, ctx_raw, ctx_res))
      cr <- safe_cor(d$mod_raw, d$ctx_raw)
      ca <- safe_cor(d$mod_res, d$ctx_res)
      direction_consistent <- is.finite(cr$r) && is.finite(ca$r) && sign(cr$r) == sign(ca$r)
      adj <- if (!is.finite(cr$r) || !is.finite(ca$r)) {
        "not_evaluable"
      } else if (sign(cr$r) != sign(ca$r) || (abs(cr$r) >= 0.10 && abs(ca$r) < 0.03)) {
        "lost_after_complexity_adjustment"
      } else if (abs(ca$r) < 0.70 * abs(cr$r)) {
        "attenuated_after_complexity_adjustment"
      } else {
        "retained_after_complexity_adjustment"
      }
      interp <- if (!is.finite(ca$r)) "not_evaluable" else if (ca$r >= 0.10) "positive_codistribution" else if (ca$r <= -0.10) "negative_codistribution" else "weak_or_inconsistent"
      codist_rows[[length(codist_rows) + 1]] <- data.table(
        section_id = sid, module_name = mod, context_signature = ctx,
        n_spots = min(cr$n, ca$n), spearman_r_raw = cr$r, spearman_p_raw = cr$p,
        spearman_r_residualized = ca$r, spearman_p_residualized = ca$p,
        direction_consistent = direction_consistent, adjustment_effect = adj,
        interpretation = interp,
        notes = "Within-section spot-level Spearman P values are descriptive only; section-level consistency determines claim strength."
      )
    }
  }
}
codist <- rbindlist(codist_rows)
fwrite(codist, file.path(table_dir, "spatial_within_section_codistribution.tsv"), sep = "\t")

consistency <- codist[, {
  eval <- is.finite(spearman_r_raw) & is.finite(spearman_r_residualized)
  n_eval <- sum(eval)
  n_pos_raw <- sum(spearman_r_raw[eval] >= 0.10)
  n_pos_res <- sum(spearman_r_residualized[eval] >= 0.10)
  med_raw <- if (n_eval) median(spearman_r_raw[eval]) else NA_real_
  med_res <- if (n_eval) median(spearman_r_residualized[eval]) else NA_real_
  range_res <- if (n_eval) paste0(formatC(min(spearman_r_residualized[eval]), digits = 3, format = "f"), " to ", formatC(max(spearman_r_residualized[eval]), digits = 3, format = "f")) else ""
  label <- if (!n_eval) {
    "not_evaluable"
  } else if (n_pos_res >= ceiling(0.8 * n_eval) && med_res >= 0.10 && abs(med_res) >= 0.70 * max(abs(med_raw), 1e-8)) {
    "consistent_positive"
  } else if (n_pos_res >= ceiling(0.6 * n_eval) && med_res >= 0.05) {
    "directionally_positive_but_attenuated"
  } else if (med_res <= 0 && n_pos_res <= floor(0.4 * n_eval)) {
    "not_supported"
  } else {
    "mixed_or_weak"
  }
  allowed <- switch(label,
    consistent_positive = "papillary_tissue_context_projection",
    directionally_positive_but_attenuated = "supplementary_spatial_codistribution",
    mixed_or_weak = "supplementary_only",
    not_supported = "not_for_claim",
    not_evaluable = "not_for_claim"
  )
  .(
    n_sections_evaluable = n_eval,
    n_sections_positive_raw = n_pos_raw,
    n_sections_positive_residualized = n_pos_res,
    median_spearman_raw = med_raw,
    median_spearman_residualized = med_res,
    range_spearman_residualized = range_res,
    consistency_label = label,
    spatial_claim_allowed = allowed,
    notes = "Positive threshold rho >= 0.10; section consistency, not spot-level P values, controls interpretation."
  )
}, by = .(module_name, context_signature)]
fwrite(consistency, file.path(table_dir, "spatial_section_level_consistency_summary.tsv"), sep = "\t")

strength_rank <- c(not_supported = 0L, descriptive_only = 1L, weak_supplementary_context = 2L, moderate_supplementary_context = 3L)
context_decision <- function(context, component) {
  d <- consistency[context_signature == context]
  counts <- d[, .N, by = consistency_label]
  n_consistent <- d[consistency_label == "consistent_positive", .N]
  n_attenuated <- d[consistency_label == "directionally_positive_but_attenuated", .N]
  n_mixed <- d[consistency_label == "mixed_or_weak", .N]
  strength <- if (n_consistent >= 2) "moderate_supplementary_context" else if (n_consistent + n_attenuated >= 2) "weak_supplementary_context" else if (n_mixed > 0 || n_consistent + n_attenuated > 0) "descriptive_only" else "not_supported"
  summary <- paste0(n_consistent, "/4 primary modules consistent_positive; ", n_attenuated, "/4 directionally_positive_but_attenuated; ", n_mixed, "/4 mixed_or_weak.")
  data.table(
    evidence_component = component,
    stage6B1_result_summary = summary,
    allowed_claim = ifelse(strength == "moderate_supplementary_context", "section-level spatial co-distribution within papillary tissue context", ifelse(strength == "weak_supplementary_context", "supplementary spatial co-distribution with attenuation caveat", ifelse(strength == "descriptive_only", "descriptive papillary tissue-context projection only", "no spatial support claim"))),
    disallowed_claim = "plaque-specific validation; spatial validation; plaque nucleation site; causal spatial niche; papilla-specific genetic regulation; cell-type-specific genetic mechanism",
    claim_strength = strength,
    manuscript_use = ifelse(strength == "moderate_supplementary_context", "supplementary_or_main_caveat", ifelse(strength == "weak_supplementary_context", "supplementary", ifelse(strength == "descriptive_only", "supplementary_descriptive", "do_not_claim"))),
    notes = "No plaque/lesion ROI; spot-level P values are descriptive only."
  )
}

claim_rows <- rbindlist(list(
  context_decision("LoopTAL_signature", "primary_MAGMA_LoopTAL_codistribution"),
  context_decision("injury_epithelial", "primary_MAGMA_injury_codistribution"),
  context_decision("ECM_fibrosis", "primary_MAGMA_ECM_codistribution"),
  context_decision("mineralization_remodeling", "primary_MAGMA_mineralization_codistribution")
))
claim_rows <- rbind(claim_rows, data.table(
  evidence_component = c("TWAS_proxy_spatial_projection", "curated_exemplar_spatial_projection"),
  stage6B1_result_summary = c(
    "TWAS-proxy modules were scored only when feasible; R2 is a one-gene descriptive module and R5 is not score-feasible.",
    "The six-gene curated exemplar panel is spatially detectable but remains interpretive context only."
  ),
  allowed_claim = c("supplementary Kidney_Cortex-proxy spatial context only", "descriptive curated-exemplar spatial projection"),
  disallowed_claim = c(
    "papilla-specific genetic regulation; causal expression; spatial validation; SMR/coloc support",
    "validated disease genes; single-gene validation; therapeutic target validation"
  ),
  claim_strength = c("descriptive_only", "descriptive_only"),
  manuscript_use = c("supplementary_only", "supplementary_descriptive"),
  notes = c("TWAS remains GTEx Kidney_Cortex proxy evidence.", "Curation does not upgrade evidence strength.")
), fill = TRUE)

max_rank <- max(strength_rank[claim_rows[1:4]$claim_strength])
overall_strength <- names(strength_rank)[match(max_rank, strength_rank)]
if (overall_strength == "moderate_supplementary_context") {
  overall_allowed <- "papillary tissue-context projection with moderate supplementary section-level spatial co-distribution"
  overall_use <- "supplementary_or_main_caveat"
} else if (overall_strength == "weak_supplementary_context") {
  overall_allowed <- "weak supplementary papillary spatial co-distribution"
  overall_use <- "supplementary"
} else if (overall_strength == "descriptive_only") {
  overall_allowed <- "descriptive papillary tissue-context projection only"
  overall_use <- "supplementary_descriptive"
} else {
  overall_allowed <- "no spatial support claim"
  overall_use <- "do_not_claim"
}
claim_rows <- rbind(claim_rows, data.table(
  evidence_component = "overall_spatial_support",
  stage6B1_result_summary = paste(claim_rows[1:4, paste(evidence_component, claim_strength, sep = "=")], collapse = "; "),
  allowed_claim = overall_allowed,
  disallowed_claim = "plaque-specific validation; spatial validation; genetic risk localizes to plaque; causal spatial niche; plaque nucleation site; papilla-specific regulatory mechanism",
  claim_strength = overall_strength,
  manuscript_use = overall_use,
  notes = "Overall claim is capped at supplementary context because no plaque/lesion ROI exists and only five sections are available."
), fill = TRUE)
fwrite(claim_rows, file.path(table_dir, "spatial_stage6B1_claim_decision_table.tsv"), sep = "\t")

cat("\nSection-level consistency counts:\n")
print(consistency[, .N, by = consistency_label][order(consistency_label)])
cat("\nClaim decisions:\n")
print(claim_rows[, .(evidence_component, claim_strength)])
cat("\nSpatial autocorrelation was not run: no preregistered/validated neighbor graph; see audit memo.\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
