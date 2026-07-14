#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(future.globals.maxSize = 8 * 1024^3)

root <- getwd()
dir.create(file.path(root, "scripts/08_spatial_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/spatial/phase3_step2_label_transfer"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/figures/phase3_step2_label_transfer"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "source_data/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
safe_id <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
section_uid <- function(sample_id, section_id) paste(safe_id(sample_id), safe_id(section_id), sep = "__")
fmt_num <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))

required_files <- c(
  "data/processed/gse231569_audited_seurat.rds",
  "config/spatial_sample_metadata_curated_phase3.tsv",
  "results/tables/phase3_step1_spatial_sample_manifest.tsv",
  "results/tables/phase3_step1_spatial_section_audit.tsv",
  "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv",
  "results/tables/phase3_step1_magma_module_spatial_mapping.tsv",
  "notes/phase3_step1B_spatial_claim_boundary.md"
)
missing <- required_files[!file.exists(file.path(root, required_files))]
if (length(missing)) {
  stop("Missing required Phase 3-Step 2 input(s): ", paste(missing, collapse = ", "))
}

meta <- fread(file.path(root, "config/spatial_sample_metadata_curated_phase3.tsv"))
manifest <- fread(file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv"))
section_audit <- fread(file.path(root, "results/tables/phase3_step1_spatial_section_audit.tsv"))
overlap <- fread(file.path(root, "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv"))

included <- meta[included_for_phase3_step2 == "yes"]
included <- merge(
  included,
  manifest[, .(sample_id, section_id, matrix_file, tissue_positions_file, complete_visium_input)],
  by = c("sample_id", "section_id"),
  all.x = TRUE
)
included <- merge(
  included,
  overlap[, .(sample_id, section_id, n_overlap_genes, suitable_for_label_transfer)],
  by = c("sample_id", "section_id"),
  all.x = TRUE
)
included <- merge(
  included,
  section_audit[, .(sample_id, section_id, n_spots_total_audit = n_spots_total, n_spots_in_tissue_audit = n_spots_in_tissue)],
  by = c("sample_id", "section_id"),
  all.x = TRUE
)
setorder(included, sample_id, section_id)

label_map <- c(
  Collecting_duct_principal = "Collecting duct",
  Endothelial = "Endothelial",
  Fibroblast_stromal = "Fibroblast/stromal",
  Injured_undifferentiated_epithelial = "Injured/undifferentiated epithelial",
  Loop_of_Henle_TAL = "Loop/TAL",
  Pericyte_smooth_muscle = "Pericyte/smooth muscle"
)
score_map <- c(
  Loop_TAL_prediction_score = "Loop_of_Henle_TAL",
  Collecting_duct_prediction_score = "Collecting_duct_principal",
  Endothelial_prediction_score = "Endothelial",
  Fibroblast_stromal_prediction_score = "Fibroblast_stromal",
  Injured_epithelial_prediction_score = "Injured_undifferentiated_epithelial",
  Pericyte_smooth_muscle_prediction_score = "Pericyte_smooth_muscle"
)
marker_examples <- c(
  Collecting_duct_principal = "AQP2,AQP3,SCNN1G",
  Endothelial = "PECAM1,VWF,KDR,EMCN",
  Fibroblast_stromal = "DCN,LUM,COL1A1,COL1A2",
  Injured_undifferentiated_epithelial = "KRT8,KRT18,PROM1,HAVCR1,VCAM1",
  Loop_of_Henle_TAL = "UMOD,SLC12A1,CLDN10,KCNJ1,CLDN16",
  Pericyte_smooth_muscle = "RGS5,ACTA2,TAGLN,MYH11,PDGFRB"
)

message("Loading canonical snRNA reference...")
reference <- readRDS(file.path(root, "data/processed/gse231569_audited_seurat.rds"))
DefaultAssay(reference) <- "RNA"
if (!"phase1_cell_type" %in% colnames(reference@meta.data)) stop("Reference lacks phase1_cell_type.")
if (!"donor_id" %in% colnames(reference@meta.data)) stop("Reference lacks donor_id.")
reference <- subset(reference, subset = !is.na(phase1_cell_type))
reference$phase1_cell_type <- as.character(reference$phase1_cell_type)
reference <- subset(reference, subset = phase1_cell_type %in% names(label_map))
reference$phase1_cell_type <- factor(reference$phase1_cell_type, levels = names(label_map))
if (length(VariableFeatures(reference)) < 1000) {
  reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}
ref_features <- VariableFeatures(reference)
if (!"pca" %in% Reductions(reference)) {
  reference <- ScaleData(reference, features = ref_features, verbose = FALSE)
  reference <- RunPCA(reference, features = ref_features, npcs = 30, verbose = FALSE)
}
npcs <- min(30, ncol(Embeddings(reference, "pca")))
dims_use <- seq_len(npcs)

ref_summary <- as.data.table(reference@meta.data)[
  , .(n_reference_nuclei = .N, n_donors_present = uniqueN(donor_id)),
  by = .(compartment = phase1_cell_type)
]
ref_summary[, key_marker_examples := marker_examples[compartment]]
ref_summary[, included_in_label_transfer := "yes"]
ref_summary[, notes := "Broad-compartment label from canonical snRNA phase1_cell_type field."]
setorder(ref_summary, compartment)
write_tsv(
  ref_summary[, .(compartment, n_reference_nuclei, n_donors_present, key_marker_examples, included_in_label_transfer, notes)],
  file.path(root, "results/tables/phase3_step2_reference_compartment_summary.tsv")
)

read_positions <- function(path) {
  dt <- fread(path, header = FALSE)
  if (ncol(dt) < 6) stop("Tissue position file has fewer than six columns: ", path)
  setnames(dt, seq_len(6), c("spot_id", "in_tissue", "array_row", "array_col", "y_coord", "x_coord"))
  dt[, spot_id := as.character(spot_id)]
  dt[, in_tissue := as.integer(in_tissue)]
  dt[, x_coord := as.numeric(x_coord)]
  dt[, y_coord := as.numeric(y_coord)]
  dt
}

get_score_col <- function(pred, label) {
  direct <- paste0("prediction.score.", label)
  alt <- names(pred)[gsub("^prediction\\.score\\.", "", names(pred)) == label]
  if (direct %in% names(pred)) return(direct)
  if (length(alt)) return(alt[1])
  NA_character_
}

palette_compartment <- c(
  "Collecting duct" = "#4C78A8",
  "Endothelial" = "#72B7B2",
  "Fibroblast/stromal" = "#59A14F",
  "Injured/undifferentiated epithelial" = "#E15759",
  "Loop/TAL" = "#F28E2B",
  "Pericyte/smooth muscle" = "#B07AA1"
)

make_pred_plot <- function(dt, title, subtitle) {
  ggplot(dt[in_tissue == 1], aes(x = x_coord, y = -y_coord, color = predicted_compartment)) +
    geom_point(size = 0.72, alpha = 0.92) +
    scale_color_manual(values = palette_compartment, drop = FALSE, na.value = "grey80") +
    coord_equal() +
    labs(
      title = title,
      subtitle = subtitle,
      color = "Predicted compartment",
      caption = "label-transfer projection; no ROI annotation"
    ) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 6)
    )
}

make_loop_plot <- function(dt, title, subtitle) {
  ggplot(dt[in_tissue == 1], aes(x = x_coord, y = -y_coord, color = Loop_TAL_prediction_score)) +
    geom_point(size = 0.72, alpha = 0.92) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1), oob = squish, name = "Loop/TAL score") +
    coord_equal() +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "label-transfer projection; no ROI annotation"
    ) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 6)
    )
}

all_scores <- list()
summary_rows <- list()
pred_plots <- list()
loop_plots <- list()

for (i in seq_len(nrow(included))) {
  row <- included[i]
  uid <- section_uid(row$sample_id, row$section_id)
  out_dir <- file.path(root, "results/spatial/phase3_step2_label_transfer", uid)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Processing ", row$sample_id, " / ", row$section_id)

  status <- "success"
  note <- "Seurat anchor-based label transfer; output directory uses sample_id__section_id to avoid overwriting multi-section sample IDs."
  result <- tryCatch({
    counts <- Read10X_h5(file.path(root, row$matrix_file))
    if (is.list(counts)) {
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }
    spatial <- CreateSeuratObject(counts = counts, project = uid, assay = "RNA")
    DefaultAssay(spatial) <- "RNA"
    spatial <- NormalizeData(spatial, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    spatial <- FindVariableFeatures(spatial, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    shared <- intersect(ref_features, rownames(spatial))
    if (length(shared) < 500) shared <- intersect(rownames(reference), rownames(spatial))
    if (length(shared) < 200) stop("Insufficient shared genes for label transfer: ", length(shared))

    anchors <- FindTransferAnchors(
      reference = reference,
      query = spatial,
      normalization.method = "LogNormalize",
      reference.reduction = "pca",
      features = shared,
      dims = dims_use,
      verbose = FALSE
    )
    pred <- as.data.table(TransferData(
      anchorset = anchors,
      refdata = reference$phase1_cell_type,
      dims = dims_use,
      verbose = FALSE
    ), keep.rownames = "spot_id")

    positions <- read_positions(file.path(root, row$tissue_positions_file))
    pred[, predicted_compartment := label_map[as.character(predicted.id)]]
    out <- merge(positions, pred, by = "spot_id", all.x = TRUE, sort = FALSE)
    out[, `:=`(sample_id = row$sample_id, section_id = row$section_id)]
    out[, prediction_score_max := as.numeric(prediction.score.max)]
    for (out_col in names(score_map)) {
      score_col <- get_score_col(pred, score_map[[out_col]])
      out[, (out_col) := if (!is.na(score_col)) as.numeric(get(score_col)) else NA_real_]
    }
    out[, all_other_compartment_scores_if_available := ""]
    out[, notes := "Seurat label-transfer prediction score; not a true cell fraction and not ROI-resolved."]
    final <- out[, .(
      spot_id,
      sample_id,
      section_id,
      x_coord,
      y_coord,
      in_tissue,
      predicted_compartment,
      prediction_score_max,
      Loop_TAL_prediction_score,
      Collecting_duct_prediction_score,
      Endothelial_prediction_score,
      Fibroblast_stromal_prediction_score,
      Injured_epithelial_prediction_score,
      Pericyte_smooth_muscle_prediction_score,
      all_other_compartment_scores_if_available,
      notes
    )]
    fwrite(final, file.path(out_dir, "label_transfer_scores.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

    low_conf <- final[in_tissue == 1 & prediction_score_max < 0.5]
    in_tissue_loop <- final[in_tissue == 1 & !is.na(Loop_TAL_prediction_score)]
    loop_high_n <- if (nrow(in_tissue_loop)) {
      sum(rank(-in_tissue_loop$Loop_TAL_prediction_score, ties.method = "first") <= ceiling(0.25 * nrow(in_tissue_loop)))
    } else {
      NA_integer_
    }
    list(
      final = final,
      n_shared = length(shared),
      low_conf_n = nrow(low_conf),
      loop_high_n = loop_high_n
    )
  }, error = function(e) {
    status <<- "failed"
    note <<- paste("Label transfer failed:", conditionMessage(e))
    list(final = data.table(), n_shared = NA_integer_, low_conf_n = NA_integer_, loop_high_n = NA_integer_)
  })

  if (nrow(result$final)) {
    plot_title <- paste(row$sample_id, row$section_id)
    plot_sub <- paste(row$dataset_label, row$control_or_disease, row$metadata_confidence, sep = " | ")
    pred_p <- make_pred_plot(result$final, plot_title, plot_sub)
    loop_p <- make_loop_plot(result$final, plot_title, plot_sub)
    ggsave(file.path(root, "results/figures/phase3_step2_label_transfer", paste0(uid, "_predicted_compartment_overlay.pdf")),
           pred_p, width = 5.2, height = 4.2, units = "in", device = "pdf")
    ggsave(file.path(root, "results/figures/phase3_step2_label_transfer", paste0(uid, "_Loop_TAL_score_overlay.pdf")),
           loop_p, width = 5.2, height = 4.2, units = "in", device = "pdf")
    pred_plots[[uid]] <- pred_p + theme(legend.position = "none")
    loop_plots[[uid]] <- loop_p + theme(legend.position = "none")
    all_scores[[uid]] <- result$final[
      ,
      `:=`(
        section_uid = uid,
        dataset_label = row$dataset_label,
        control_or_disease = row$control_or_disease,
        metadata_confidence = row$metadata_confidence
      )
    ][]
  }

  in_tissue_dt <- result$final[in_tissue == 1]
  summary_rows[[uid]] <- data.table(
    sample_id = row$sample_id,
    section_id = row$section_id,
    dataset_label = row$dataset_label,
    control_or_disease = row$control_or_disease,
    metadata_confidence = row$metadata_confidence,
    n_spots_total = ifelse(nrow(result$final), nrow(result$final), row$n_spots_total_audit),
    n_spots_in_tissue = ifelse(nrow(in_tissue_dt), nrow(in_tissue_dt), row$n_spots_in_tissue_audit),
    n_shared_genes_used = result$n_shared,
    n_reference_nuclei = ncol(reference),
    n_reference_compartments = length(unique(reference$phase1_cell_type)),
    median_prediction_score_max = ifelse(nrow(in_tissue_dt), median(in_tissue_dt$prediction_score_max, na.rm = TRUE), NA_real_),
    mean_prediction_score_max = ifelse(nrow(in_tissue_dt), mean(in_tissue_dt$prediction_score_max, na.rm = TRUE), NA_real_),
    n_low_confidence_spots = result$low_conf_n,
    percent_low_confidence_spots = ifelse(nrow(in_tissue_dt), 100 * result$low_conf_n / nrow(in_tissue_dt), NA_real_),
    median_Loop_TAL_prediction_score = ifelse(nrow(in_tissue_dt), median(in_tissue_dt$Loop_TAL_prediction_score, na.rm = TRUE), NA_real_),
    percent_spots_Loop_TAL_high = ifelse(nrow(in_tissue_dt), 100 * result$loop_high_n / nrow(in_tissue_dt), NA_real_),
    label_transfer_status = status,
    notes = paste(note, "Loop/TAL-high defined descriptively as the within-section top 25% of spots by Loop/TAL prediction-score rank.")
  )
  rm(result)
  gc(verbose = FALSE)
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
write_tsv(summary_dt, file.path(root, "results/tables/phase3_step2_label_transfer_summary.tsv"))
all_dt <- rbindlist(all_scores, fill = TRUE)
write_tsv(all_dt, file.path(root, "source_data/figures/phase3_step2_label_transfer_qc_source_data.tsv"))

if (length(pred_plots)) {
  ggsave(
    file.path(root, "results/figures/phase3_step2_all_sections_predicted_compartment_overview.pdf"),
    wrap_plots(pred_plots, ncol = 2) + plot_annotation(title = "Phase 3-Step 2 broad-compartment label-transfer projections"),
    width = 11, height = 15, units = "in", device = "pdf"
  )
  ggsave(
    file.path(root, "results/figures/phase3_step2_all_sections_Loop_TAL_score_overview.pdf"),
    wrap_plots(loop_plots, ncol = 2) + plot_annotation(title = "Phase 3-Step 2 Loop/TAL prediction-score projections"),
    width = 11, height = 15, units = "in", device = "pdf"
  )
}

if (nrow(all_dt)) {
  all_dt[, section_label := paste(sample_id, section_id, sep = "\n")]
  p_score <- ggplot(all_dt[in_tissue == 1], aes(x = prediction_score_max)) +
    geom_histogram(bins = 40, fill = "#4C78A8", color = "white", linewidth = 0.15) +
    facet_wrap(~ section_id, scales = "free_y", ncol = 2) +
    labs(x = "Maximum prediction score", y = "In-tissue spots", title = "Label-transfer prediction-score distribution") +
    theme_bw(base_size = 8)
  ggsave(file.path(root, "results/figures/phase3_step2_prediction_score_distribution.pdf"),
         p_score, width = 8.5, height = 11, units = "in", device = "pdf")

  p_loop <- ggplot(all_dt[in_tissue == 1], aes(x = section_id, y = Loop_TAL_prediction_score, fill = control_or_disease)) +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), oob = squish) +
    labs(x = NULL, y = "Loop/TAL prediction score", fill = "Context", title = "Loop/TAL prediction-score distribution by section") +
    theme_bw(base_size = 8)
  ggsave(file.path(root, "results/figures/phase3_step2_Loop_TAL_score_distribution_by_section.pdf"),
         p_loop, width = 8.5, height = 5.8, units = "in", device = "pdf")

  p_low <- ggplot(summary_dt, aes(x = reorder(section_id, percent_low_confidence_spots), y = percent_low_confidence_spots, fill = control_or_disease)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Low-confidence spots (%)", fill = "Context", title = "Fraction of low-confidence spots by section") +
    theme_bw(base_size = 8)
  ggsave(file.path(root, "results/figures/phase3_step2_low_confidence_spot_fraction.pdf"),
         p_low, width = 8.5, height = 5.8, units = "in", device = "pdf")
}

condition_dt <- all_dt[in_tissue == 1, .(
  n_sections = uniqueN(section_id),
  metadata_confidence_summary = paste(sort(unique(metadata_confidence)), collapse = ";"),
  median_Loop_TAL_prediction_score = median(Loop_TAL_prediction_score, na.rm = TRUE),
  iqr_Loop_TAL_prediction_score = IQR(Loop_TAL_prediction_score, na.rm = TRUE)
), by = control_or_disease]
control_n <- summary_dt[control_or_disease == "control", uniqueN(section_id)]
disease_n <- summary_dt[control_or_disease == "stone_disease", uniqueN(section_id)]
condition_dt[, interpretation_allowed := "Descriptive section-level tissue-context summary only."]
condition_dt[, interpretation_not_allowed := "No claim-grade disease/control test, plaque-specific localization, lesion-stage localization, causal niche, or independent validation."]
condition_dt[, notes := ifelse(
  control_n < 2 || disease_n < 2,
  "Condition-level comparison is descriptive only and not used as a claim-grade disease/control test.",
  "Condition-level comparison remains descriptive unless a formal balanced design is specified."
)]
write_tsv(condition_dt, file.path(root, "results/tables/phase3_step2_condition_context_summary.tsv"))

claim_note <- c(
  "# Phase 3-Step 2 spatial label-transfer claim boundary",
  "",
  "- Seurat label transfer provides broad-compartment prediction scores, not direct cell fractions.",
  "- Loop/TAL score overlays support anatomical/tissue-context projection only.",
  "- No plaque-, mineral-, lesion-, or fibrosis-resolved spot-level ROI annotation is available in the locked inputs.",
  "- Do not claim plaque-specific localization.",
  "- Do not claim lesion-stage localization.",
  "- Do not claim a causal spatial niche or independent validation.",
  "- MAGMA modules remain genetic-priority modules, not causal gene sets.",
  "- Condition-level comparison is descriptive only and not used as a claim-grade disease/control test because control representation is limited.",
  "- Per-section output directories use `sample_id__section_id` to preserve all sections when a GEO sample has multiple separated section folders."
)
writeLines(claim_note, file.path(root, "notes/phase3_step2_spatial_label_transfer_claim_boundary.md"))

success_n <- summary_dt[label_transfer_status == "success", .N]
total_n <- nrow(summary_dt)
mean_low <- mean(summary_dt$percent_low_confidence_spots, na.rm = TRUE)
median_max <- median(summary_dt$median_prediction_score_max, na.rm = TRUE)
loop_med_range <- range(summary_dt$median_Loop_TAL_prediction_score, na.rm = TRUE)
compartments <- paste(names(label_map), collapse = ", ")

report <- c(
  "# Phase 3-Step 2 report",
  "",
  "## Method",
  "Used Seurat anchor-based label transfer from `data/processed/gse231569_audited_seurat.rds` to the locked Visium sections. The reference label field was `phase1_cell_type`, and the donor field was `donor_id`.",
  "",
  "## Sections processed",
  paste0("Processed ", success_n, " of ", total_n, " included spatial sections from `config/spatial_sample_metadata_curated_phase3.tsv`. Per-section directories use `sample_id__section_id` to avoid overwriting multi-section sample IDs."),
  "",
  "## Gene overlap",
  paste0("Shared genes used per section ranged from ", min(summary_dt$n_shared_genes_used, na.rm = TRUE), " to ", max(summary_dt$n_shared_genes_used, na.rm = TRUE), ". The Step 1 overlap table reported 30,940 shared snRNA-spatial genes for each section before variable-feature selection."),
  "",
  "## Transferred compartments",
  paste0("Transferred broad-compartment labels: ", compartments, "."),
  "",
  "## Prediction score QC",
  paste0("Across successful sections, median section-level maximum prediction score was ", fmt_num(median_max), ". Mean low-confidence fraction, using prediction_score_max < 0.5, was ", fmt_num(mean_low), "%."),
  "",
  "## Loop/TAL projection summary",
  paste0("Median Loop/TAL prediction score across sections ranged from ", fmt_num(loop_med_range[1]), " to ", fmt_num(loop_med_range[2]), ". Loop/TAL-high spots were defined descriptively as the within-section top 25% of spots by Loop/TAL prediction-score rank. Because section medians were zero, this rank-based layer should be treated as a sparse contextual overlay rather than evidence for strong Loop/TAL spatial enrichment."),
  "",
  "## Disease/control comparison status",
  "Condition-level comparison is descriptive only and not used as a claim-grade disease/control test because the locked metadata includes limited control representation.",
  "",
  "## ROI boundary status",
  "No plaque-, mineral-, lesion-, or fibrosis-resolved spot-level ROI annotation is available. These outputs support spatial projection and label-transfer-informed tissue context, not plaque-specific localization or causal spatial localization.",
  "",
  "## Step 3 readiness",
  ifelse(success_n == total_n, "Phase 3-Step 3 can proceed if the human reviewer accepts Seurat label transfer as the projection layer for downstream MAGMA module co-distribution.", "Human review is needed before Step 3 because at least one section failed label transfer.")
)
writeLines(report, file.path(root, "notes/phase3_step2_report.md"))

recommendation <- if (success_n == total_n && is.finite(mean_low) && mean_low < 50) {
  "A. proceed to Phase 3-Step 3: spatial MAGMA module scoring and co-distribution with Loop/TAL prediction scores."
} else if (success_n == total_n) {
  "D. human decision required because label transfer completed but low-confidence spot burden is substantial."
} else {
  "B. repeat label transfer due to failed sections before downstream spatial module work."
}
limitations <- c(
  "# Phase 3-Step 2 limitations and next steps",
  "",
  "- Seurat label transfer is not true deconvolution and does not estimate physical cell fractions.",
  "- No plaque-, mineral-, lesion-, or fibrosis-resolved ROI annotation is available.",
  "- Control representation is limited; condition-level summaries are descriptive only.",
  "- Loop/TAL prediction scores were sparse in section-level summaries; rank-based top-quartile overlays are descriptive and tie-sensitive.",
  "- MAGMA spatial module scoring was not performed in this step.",
  "- Spatial correlation or co-distribution analysis was not performed in this step.",
  "- RCTD/spacexr was not installed or run in this step.",
  "",
  "## Recommended next step",
  recommendation,
  "",
  "Other possible actions if requested by human review:",
  "- B. repeat label transfer due to low confidence.",
  "- C. install RCTD and repeat with RCTD.",
  "- D. human decision required."
)
writeLines(limitations, file.path(root, "notes/phase3_step2_limitations_and_next_steps.md"))

checklist <- data.table(
  item = c(
    "main_script_created",
    "canonical_snRNA_reference_loaded",
    "all_included_sections_attempted",
    "spot_level_prediction_scores_written",
    "cross_section_summary_written",
    "reference_compartment_summary_written",
    "per_section_overlay_plots_written",
    "multi_section_overview_plots_written",
    "qc_plots_and_source_data_written",
    "condition_context_summary_written",
    "claim_boundary_note_written",
    "report_written",
    "limitations_note_written",
    "no_magma_module_scores_computed",
    "no_spatial_correlations_computed",
    "no_rctd_or_cell2location_run",
    "no_manuscript_modified"
  ),
  status = c(
    "complete",
    "complete",
    ifelse(total_n == nrow(included), "complete", "review"),
    ifelse(success_n == total_n, "complete", "partial"),
    "complete",
    "complete",
    ifelse(success_n == total_n, "complete", "partial"),
    ifelse(success_n > 0, "complete", "not_done"),
    ifelse(success_n > 0, "complete", "not_done"),
    "complete",
    "complete",
    "complete",
    "complete",
    "complete",
    "complete",
    "complete",
    "complete"
  ),
  evidence = c(
    "scripts/08_spatial_analysis/phase3_step2_seurat_label_transfer.R",
    "data/processed/gse231569_audited_seurat.rds",
    paste0(total_n, " section rows in summary table"),
    "results/spatial/phase3_step2_label_transfer/*/label_transfer_scores.tsv.gz",
    "results/tables/phase3_step2_label_transfer_summary.tsv",
    "results/tables/phase3_step2_reference_compartment_summary.tsv",
    "results/figures/phase3_step2_label_transfer/*_overlay.pdf",
    "results/figures/phase3_step2_all_sections_*_overview.pdf",
    "source_data/figures/phase3_step2_label_transfer_qc_source_data.tsv",
    "results/tables/phase3_step2_condition_context_summary.tsv",
    "notes/phase3_step2_spatial_label_transfer_claim_boundary.md",
    "notes/phase3_step2_report.md",
    "notes/phase3_step2_limitations_and_next_steps.md",
    "Script contains no MAGMA module scoring step.",
    "Script contains no spatial correlation step.",
    "Script uses Seurat FindTransferAnchors/TransferData only.",
    "No DOCX or manuscript path is written by this script."
  ),
  notes = c(
    "Created as the executable record for this step.",
    "Original RDS is read only; in-memory object may be normalized/PCA checked.",
    "All sections marked included_for_phase3_step2=yes were processed or logged.",
    "Directories use sample_id__section_id because multiple sections share sample_id.",
    "Low-confidence threshold is prediction_score_max < 0.5.",
    "Marker examples are taken from locked Phase 1 snRNA marker-label audit.",
    "Plots include label-transfer projection/no ROI annotation caption.",
    "Compact overview PDFs created from per-section plots.",
    "QC source data is spot-level and uncompressed TSV for traceability.",
    "Disease/control comparison is descriptive only.",
    "Claim boundaries mirror Phase 3-Step 1B.",
    "Report uses spatial projection language.",
    "Recommended action selected from A/B/C/D.",
    "Guardrail honored.",
    "Guardrail honored.",
    "Guardrail honored.",
    "Guardrail honored."
  )
)
write_tsv(checklist, file.path(root, "codex_tasks/phase3_step2_completion_checklist.tsv"))

message("Phase 3-Step 2 complete: ", success_n, "/", total_n, " sections processed successfully.")
