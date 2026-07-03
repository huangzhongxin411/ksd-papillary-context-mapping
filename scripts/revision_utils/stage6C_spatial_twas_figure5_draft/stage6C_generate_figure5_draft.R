#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
})

stage_id <- "stage6C_spatial_twas_figure5_draft"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
for (d in c(fig_dir, table_dir, doc_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "stage6C_generate_figure5_draft.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Stage 6C conservative Figure 5 draft generation\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

paths <- list(
  section_qc = "results/tables/revision/stage6B1_spatial_context_projection/spatial_section_qc_stage6B1.tsv",
  spot_scores = "results/tables/revision/stage6B1_spatial_context_projection/spatial_spot_module_scores.tsv",
  consistency = "results/tables/revision/stage6B1_spatial_context_projection/spatial_section_level_consistency_summary.tsv",
  within_section = "results/tables/revision/stage6B1_spatial_context_projection/spatial_within_section_codistribution.tsv",
  claim_decision = "results/tables/revision/stage6B1_spatial_context_projection/spatial_stage6B1_claim_decision_table.tsv",
  module_manifest = "results/tables/revision/stage6B1_spatial_context_projection/spatial_module_manifest_stage6B1.tsv",
  twas_boundary = "results/tables/revision/stage6A_spatial_twas_audit/twas_boundary_integration_audit.tsv",
  integrated_boundary = "results/tables/revision/stage6A_spatial_twas_audit/stage6_integrated_claim_boundary_table.tsv",
  roi_audit = "docs/revision/stage6A_spatial_twas_audit/spatial_roi_plaque_annotation_audit.md"
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Required Stage 6A/6B1 source missing: ", paste(missing, collapse = ", "))

check_cols <- function(d, cols, label) {
  absent <- setdiff(cols, names(d))
  if (length(absent)) stop("Missing required columns in ", label, ": ", paste(absent, collapse = ", "))
}

section_qc <- fread(paths$section_qc)
spot_scores <- fread(paths$spot_scores)
consistency <- fread(paths$consistency)
within_section <- fread(paths$within_section)
claim_decision <- fread(paths$claim_decision)
module_manifest <- fread(paths$module_manifest)
twas_boundary <- fread(paths$twas_boundary)

check_cols(section_qc, c("section_id", "n_tissue_spots", "median_nFeature_Spatial", "image_available", "roi_available"), "section QC")
check_cols(spot_scores, c("global_spot_id", "section_id", "x_coordinate", "y_coordinate", "module_name", "module_score"), "spot scores")
check_cols(consistency, c("module_name", "context_signature", "median_spearman_raw", "median_spearman_residualized", "consistency_label", "spatial_claim_allowed"), "consistency summary")
check_cols(within_section, c("section_id", "module_name", "context_signature", "spearman_r_raw", "spearman_r_residualized", "adjustment_effect"), "within-section co-distribution")
check_cols(claim_decision, c("evidence_component", "claim_strength", "allowed_claim", "disallowed_claim"), "claim decision")
check_cols(module_manifest, c("module_name", "n_genes_input", "n_genes_detected_union", "module_status"), "module manifest")
check_cols(twas_boundary, c("twas_component", "n_genes_or_models", "current_status", "key_limitation"), "TWAS boundary")

primary_modules <- c("R1_MAGMA_Bonferroni_only", "R1_R2_R3_all_MAGMA_Bonferroni", "MAGMA_top50", "MAGMA_top100")
key_contexts <- c("LoopTAL_signature", "injury_epithelial", "ECM_fibrosis", "mineralization_remodeling")
module_labels <- c(
  R1_MAGMA_Bonferroni_only = "R1 MAGMA Bonf.",
  R1_R2_R3_all_MAGMA_Bonferroni = "R1-R3 MAGMA Bonf.",
  MAGMA_top50 = "MAGMA top 50",
  MAGMA_top100 = "MAGMA top 100"
)
context_labels <- c(
  LoopTAL_signature = "Loop / TAL",
  injury_epithelial = "Injury epithelial",
  ECM_fibrosis = "ECM / fibrosis",
  mineralization_remodeling = "Mineralization / remodeling"
)

write_source <- function(d, filename) {
  path <- file.path(table_dir, filename)
  out <- copy(d)
  char_cols <- names(out)[vapply(out, is.character, logical(1))]
  for (nm in char_cols) set(out, j = nm, value = gsub("[\r\n]+", " | ", out[[nm]]))
  fwrite(out, path, sep = "\t", na = "")
  cat("Wrote source:", path, "\n")
  path
}

# Panel A: resource overview and representative-section selection independent of correlation.
eligible <- section_qc[image_available == TRUE]
if (!nrow(eligible)) stop("No image-available section for representative visualization")
target_nfeature <- median(eligible$median_nFeature_Spatial)
eligible[, display_distance := abs(median_nFeature_Spatial - target_nfeature)]
setorder(eligible, display_distance, -n_tissue_spots, section_id)
representative_section <- eligible$section_id[1]
panelA <- copy(section_qc)
panelA[, `:=`(
  representative_for_display = section_id == representative_section,
  representative_selection_rule = "image available; median nFeature closest to the five-section median; ties resolved by larger tissue-spot count; correlation not used",
  no_roi_boundary = "No plaque/mineral/lesion ROI annotation"
)]
panelA_path <- write_source(panelA, "figure5_panelA_source.tsv")
cat("Representative display section:", representative_section, "(QC rule; not correlation)\n")

# Panel B: one representative section, one primary module and Loop/TAL signature.
panelB <- spot_scores[
  section_id == representative_section & module_name %in% c("MAGMA_top50", "LoopTAL_signature")
]
panelB[, module_label := factor(
  fifelse(module_name == "MAGMA_top50", "MAGMA top 50", "Loop / TAL signature"),
  levels = c("MAGMA top 50", "Loop / TAL signature")
)]
panelB[, `:=`(
  representative_selection_rule = "closest median nFeature among image-available sections; correlation not used",
  interpretation_boundary = "visualization only; not lesion localization or plaque enrichment"
)]
if (nrow(panelB) != 2L * section_qc[section_id == representative_section, n_tissue_spots]) stop("Panel B expected two score rows per representative spot")
panelB_path <- write_source(panelB, "figure5_panelB_source.tsv")

# Resolve representative histology image and scaling information.
spatial_dir <- file.path("data/processed/spatial/gse206306", representative_section, "spatial")
image_path <- file.path(spatial_dir, "tissue_hires_image.png")
scale_path <- file.path(spatial_dir, "scalefactors_json.json")
if (!file.exists(image_path) || !file.exists(scale_path)) stop("Representative histology image or scale factor missing: ", representative_section)
img <- png::readPNG(image_path)
scale_json <- jsonlite::fromJSON(scale_path)
hires_scale <- as.numeric(scale_json$tissue_hires_scalef)
if (!is.finite(hires_scale)) stop("Invalid tissue_hires_scalef for ", representative_section)
if (length(dim(img)) == 3L) {
  img[, , 1:3] <- 0.70 * img[, , 1:3] + 0.30
} else {
  img <- 0.70 * img + 0.30
}
image_width <- dim(img)[2]
image_height <- dim(img)[1]
panelB[, `:=`(
  x_hires = x_coordinate * hires_scale,
  y_hires_plot = image_height - y_coordinate * hires_scale,
  histology_image = image_path,
  tissue_hires_scalef = hires_scale
)]
# Re-write after adding display coordinates and provenance.
panelB_path <- write_source(panelB, "figure5_panelB_source.tsv")

# Panel C: adjusted section-level summary for four primary modules x four contexts.
panelC <- consistency[module_name %in% primary_modules & context_signature %in% key_contexts]
panelC[, `:=`(
  module_label = factor(unname(module_labels[module_name]), levels = rev(unname(module_labels))),
  context_label = factor(unname(context_labels[context_signature]), levels = unname(context_labels)),
  consistency_short = fifelse(consistency_label == "consistent_positive", "consistent",
    fifelse(consistency_label == "directionally_positive_but_attenuated", "attenuated",
      fifelse(consistency_label == "mixed_or_weak", "mixed", "not supported"))),
  cell_label = paste0(formatC(median_spearman_residualized, digits = 2, format = "f"), "\n", fifelse(consistency_label == "consistent_positive", "CP", fifelse(consistency_label == "directionally_positive_but_attenuated", "ATT", fifelse(consistency_label == "mixed_or_weak", "MIX", "NS"))))
)]
if (nrow(panelC) != 16L) stop("Panel C expected 16 rows; found ", nrow(panelC))
panelC_path <- write_source(panelC, "figure5_panelC_source.tsv")

# Panel D: raw versus complexity-adjusted values for every module-section pair.
panelD <- within_section[module_name %in% primary_modules & context_signature %in% key_contexts]
panelD[, `:=`(
  module_label = unname(module_labels[module_name]),
  context_label = factor(unname(context_labels[context_signature]), levels = unname(context_labels)),
  context_label_short = factor(
    fifelse(context_signature == "LoopTAL_signature", "Loop / TAL",
      fifelse(context_signature == "injury_epithelial", "Injury",
        fifelse(context_signature == "ECM_fibrosis", "ECM / fibrosis", "Mineralization"))),
    levels = c("Loop / TAL", "Injury", "ECM / fibrosis", "Mineralization")
  ),
  pair_id = paste(section_id, module_name, sep = ":")
)]
panelD_path <- write_source(panelD, "figure5_panelD_source.tsv")

# Panel E: TWAS count and module-boundary source.
one_row <- twas_boundary[twas_component == "one-SNP TWAS models"][1]
multi_row <- twas_boundary[twas_component == "multi-SNP TWAS models"][1]
r2_row <- module_manifest[module_name == "R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy"][1]
r5_row <- module_manifest[module_name == "R5_TWAS_proxy_only"][1]
panelE <- rbindlist(list(
  data.table(component = "FDR-supported TWAS", category = "one-SNP", value = one_row$n_genes_or_models, status = one_row$current_status, display_note = "weaker proxy"),
  data.table(component = "FDR-supported TWAS", category = "multi-SNP", value = multi_row$n_genes_or_models, status = multi_row$current_status, display_note = "supplementary proxy"),
  data.table(component = "R2 spatial module", category = "one-gene descriptive", value = r2_row$n_genes_detected_union, status = r2_row$module_status, display_note = "not for main claim"),
  data.table(component = "R5 spatial module", category = "not score-feasible", value = r5_row$n_genes_detected_union, status = r5_row$module_status, display_note = "1/4 genes detected"),
  data.table(component = "TWAS boundary", category = "GTEx v8 Kidney_Cortex proxy", value = NA_real_, status = "boundary_only", display_note = "not papilla-specific regulation; not SMR/coloc support")
), fill = TRUE)
panelE_path <- write_source(panelE, "figure5_panelE_source.tsv")

# Panel F: explicit final boundary strip.
panelF <- data.table(
  display_order = 1:6,
  statement = c(
    "Moderate supplementary spatial context",
    "Loop/TAL-associated tissue-context projection",
    "No plaque ROI",
    "Not plaque-specific validation",
    "Not causal spatial niche",
    "TWAS = Kidney_Cortex proxy"
  ),
  statement_display = c(
    "Moderate supplementary\nspatial context",
    "Loop/TAL-associated\ntissue-context projection",
    "No plaque ROI",
    "Not plaque-specific\nvalidation",
    "Not causal\nspatial niche",
    "TWAS = Kidney_Cortex\nproxy"
  ),
  statement_type = c("allowed", "allowed", "boundary", "forbidden", "forbidden", "boundary"),
  notes = c(
    "Overall Stage 6C support ceiling.", "Only spatial context retained after adjustment.",
    "No plaque/mineral/lesion/fibrosis ROI annotation.", "Spatial projection is not validation.",
    "No causal or plaque-nucleation inference.", "GTEx tissue proxy, not papilla-specific regulation."
  )
)
panelF_path <- write_source(panelF, "figure5_panelF_source.tsv")

# Final claim lock.
claim_lock <- data.table(
  evidence_component = c(
    "spatial_resource_layer", "spatial_no_roi_boundary",
    "primary_MAGMA_LoopTAL_codistribution", "primary_MAGMA_injury_codistribution",
    "primary_MAGMA_ECM_codistribution", "primary_MAGMA_mineralization_codistribution",
    "TWAS_proxy_layer", "curated_exemplar_spatial_projection", "overall_spatial_TWAS_claim"
  ),
  support_level = c(
    "moderate_supplementary_context", "boundary_only",
    "moderate_supplementary_context", "not_supported", "not_supported",
    "weak_or_mixed_context", "descriptive_only", "descriptive_only",
    "moderate_supplementary_context"
  ),
  allowed_wording = c(
    "Five QC-passed GSE206306 sections support papillary tissue-context projection.",
    "No plaque/mineral/lesion/fibrosis ROI annotation was available; spatial interpretation is not plaque-specific.",
    "Moderate supplementary Loop/TAL-associated section-level spatial co-distribution for MAGMA-prioritized modules.",
    "Injury co-distribution was not retained after complexity adjustment.",
    "ECM/fibrosis co-distribution was not retained after complexity adjustment.",
    "Mineralization/remodeling spatial co-distribution was mixed/weak and descriptive only.",
    "GTEx v8 Kidney_Cortex TWAS proxy evidence is supplementary only.",
    "Curated exemplars provide descriptive spatial projection without evidence upgrading.",
    "moderate supplementary Loop/TAL-associated papillary tissue-context projection for MAGMA-prioritized modules"
  ),
  forbidden_wording = c(
    "validation cohort; plaque-specific sampling",
    "plaque-specific validation; genetic risk localizes to plaque; plaque nucleation site",
    "spatial validation; causal cell type; causal papillary niche",
    "spatial injury mechanism; plaque injury niche",
    "ECM mechanism; fibrotic plaque localization",
    "mineralization mechanism; calcification or plaque nucleation site",
    "papilla-specific genetic regulation; causal expression; SMR/coloc-supported gene",
    "validated disease gene; therapeutic target; single-gene validation",
    "plaque-specific validation; spatial validation; genetic risk localizes to plaque; plaque nucleation site; causal papillary niche; papilla-specific genetic regulation; therapeutic target; SMR/coloc-supported gene"
  ),
  main_or_supplement = c("supplementary", "boundary", "supplementary", "supplementary_null", "supplementary_null", "supplementary_descriptive", "supplementary", "supplementary_descriptive", "supplementary"),
  notes = c(
    "Five sections; 7,747 spots; spots are visualization units.",
    "Hard interpretation boundary.",
    "Two modules consistent positive and two directionally positive but attenuated.",
    "Adjusted median rho near zero/slightly negative.",
    "Adjusted median rho near zero/slightly negative.",
    "All four modules mixed/weak after adjustment.",
    "51 FDR-supported genes: 42 one-SNP and 9 multi-SNP models.",
    "Six-gene panel detectable but interpretive only.",
    "Locked Stage 6C wording."
  )
)
claim_lock_path <- file.path(table_dir, "spatial_twas_final_claim_lock_stage6C.tsv")
fwrite(claim_lock, claim_lock_path, sep = "\t")

theme_stage <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#333333"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = max(base_size - 0.2, 8), color = "#555555", lineheight = 0.95),
      axis.title = element_text(size = max(base_size, 9.5)),
      axis.text = element_text(size = max(base_size - 0.2, 8.5), color = "#333333"),
      legend.title = element_text(face = "bold", size = max(base_size - 0.1, 8.5)),
      legend.text = element_text(size = max(base_size - 0.3, 8)),
      plot.margin = margin(6, 6, 6, 6)
    )
}

panel_label <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = 14, color = "#222222"),
          plot.tag.position = c(0, 1))
}

# A: section resources and no-ROI boundary.
panelA[, section_label := factor(section_id, levels = rev(section_id))]
pA <- ggplot(panelA, aes(y = section_label, x = n_tissue_spots, fill = representative_for_display)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = format(n_tissue_spots, big.mark = ",")), hjust = -0.12, size = 3.0, color = "#333333") +
  scale_fill_manual(values = c("TRUE" = "#B99B5A", "FALSE" = "#7F9DA6"), guide = "none") +
  coord_cartesian(xlim = c(0, max(panelA$n_tissue_spots) * 1.22), clip = "off") +
  labs(
    title = "Five spatial sections; no ROI",
    subtitle = paste0("5 sections | 7,747 spots | No plaque/mineral/lesion ROI\n", representative_section, " marked in amber: QC-selected for display, not correlation"),
    caption = "Papillary tissue-context projection only; not plaque-specific validation",
    x = "Tissue spots", y = NULL
  ) +
  theme_stage(9)

# B: representative overlay. Evidence remains in C/D.
panelB_color_limits <- as.numeric(quantile(panelB$module_score, c(0.02, 0.98), na.rm = TRUE))
pB <- ggplot(panelB, aes(x = x_hires, y = y_hires_plot)) +
  annotation_raster(img, xmin = 0, xmax = image_width, ymin = 0, ymax = image_height) +
  geom_point(aes(color = module_score), size = 0.75, alpha = 0.82) +
  facet_wrap(~module_label, nrow = 1) +
  scale_color_gradient2(low = "#7F9DA6", mid = "#E6E9EA", high = "#9B5C4D", midpoint = 0, limits = panelB_color_limits, oob = scales::squish, name = "Mean-z\nscore") +
  coord_fixed(xlim = c(0, image_width), ylim = c(0, image_height), expand = FALSE) +
  labs(title = "Representative papillary spatial projection", subtitle = paste0(representative_section, "; visualization only, not lesion localization"), x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8.5, color = "#555555"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right", legend.title = element_text(face = "bold", size = 8.5), legend.text = element_text(size = 8),
    plot.margin = margin(6, 6, 6, 6)
  )

# C: main adjusted evidence matrix.
pC <- ggplot(panelC, aes(x = context_label, y = module_label, fill = median_spearman_residualized)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.9, color = "#222222") +
  scale_fill_gradient2(low = "#9B5C4D", mid = "#F2F2F2", high = "#0F4C5C", midpoint = 0, limits = c(-0.15, 0.15), oob = scales::squish, name = "Adjusted\nmedian rho") +
  labs(title = "Section-level adjusted co-distribution", subtitle = "CP = consistent; ATT = positive but attenuated; MIX = mixed; NS = not supported", x = NULL, y = NULL) +
  theme_stage(8.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

# D: raw versus adjusted values, with all module-section pairs visible.
panelD_long <- melt(
  panelD,
  id.vars = c("section_id", "module_name", "module_label", "context_signature", "context_label", "context_label_short", "pair_id"),
  measure.vars = c("spearman_r_raw", "spearman_r_residualized"),
  variable.name = "score_type", value.name = "spearman_r"
)
panelD_long[, score_type := factor(score_type, levels = c("spearman_r_raw", "spearman_r_residualized"), labels = c("Raw", "Complexity-adjusted"))]
panelD_median <- panelD_long[, .(median_rho = median(spearman_r)), by = .(context_label_short, score_type)]
pD <- ggplot(panelD_long, aes(x = score_type, y = spearman_r, group = pair_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_line(color = "#AEB7BA", alpha = 0.45, linewidth = 0.35) +
  geom_point(aes(color = score_type), size = 1.45, alpha = 0.68) +
  geom_line(data = panelD_median, aes(x = score_type, y = median_rho, group = 1), inherit.aes = FALSE, color = "#222222", linewidth = 1.0) +
  geom_point(data = panelD_median, aes(x = score_type, y = median_rho), inherit.aes = FALSE, color = "#222222", fill = "white", shape = 21, size = 3.0, stroke = 0.8) +
  facet_wrap(~context_label_short, nrow = 1) +
  scale_color_manual(values = c("Raw" = "#7F9DA6", "Complexity-adjusted" = "#245A64"), guide = "none") +
  labs(title = "Complexity adjustment changes the context signal", subtitle = "Thin lines = section x primary-module pairs; black line = median; no spot-level P values shown", x = NULL, y = "Within-section Spearman rho") +
  theme_stage(8.3) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), strip.background = element_rect(fill = "#F4F5F5", color = NA), strip.text = element_text(face = "bold", size = 8))

# E: TWAS proxy boundary.
twas_plot <- panelE[component == "FDR-supported TWAS"]
twas_plot[, category := factor(category, levels = c("one-SNP", "multi-SNP"))]
pE <- ggplot(twas_plot, aes(x = "GTEx v8 Kidney_Cortex", y = value, fill = category)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(category, "\n", value)), position = position_stack(vjust = 0.5), size = 3.1, fontface = "bold", lineheight = 0.9, color = "#222222") +
  scale_fill_manual(values = c("one-SNP" = "#C88A78", "multi-SNP" = "#7F9DA6"), guide = "none") +
  annotate("text", x = 0.18, y = 0, hjust = 0, label = "51 FDR-supported genes: 42 one-SNP (weaker proxy) | 9 multi-SNP\nR2: one-gene descriptive | R5: not score-feasible\nNot papilla-specific regulation; not SMR/coloc support", size = 3.0, lineheight = 0.98, color = "#444444") +
  coord_flip(xlim = c(0.02, 1.45), ylim = c(0, 55), clip = "off") +
  labs(title = "TWAS boundary", subtitle = "Supplementary GTEx Kidney_Cortex proxy only", x = NULL, y = "FDR-supported genes") +
  theme_stage(8.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(6, 12, 6, 6))

# F: final interpretation strip.
pF <- ggplot(panelF, aes(x = display_order, y = 1, fill = statement_type)) +
  geom_tile(width = 0.94, height = 0.72, color = "white", linewidth = 0.8) +
  geom_text(aes(label = statement_display), size = 3.0, lineheight = 0.92, fontface = "bold", color = "#2F2F2F") +
  scale_fill_manual(values = c(allowed = "#DDE8EA", boundary = "#F3E9D2", forbidden = "#F1DDD7"), guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = NULL, expand = expansion(add = 0.05)) +
  coord_cartesian(ylim = c(0.55, 1.45), clip = "off") +
  labs(title = "Locked spatial/TWAS interpretation") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(6, 6, 6, 6))

figure5 <- (
  panel_label(pA, "A") | panel_label(pB, "B")
) / (
  panel_label(pC, "C") | panel_label(pD, "D")
) / (
  panel_label(pE, "E")
) / panel_label(pF, "F") +
  plot_layout(heights = c(0.92, 1.08, 0.58, 0.42), widths = c(0.92, 1.28), guides = "keep") +
  plot_annotation(
    title = "Figure 5 draft v0.1. Spatial and TWAS layers provide supplementary papillary context",
    subtitle = "Loop/TAL co-distribution is retained or attenuated after complexity adjustment; injury/ECM are not retained; no plaque ROI is available",
    theme = theme(
      plot.title = element_text(face = "bold", family = "sans", size = 17, color = "#222222"),
      plot.subtitle = element_text(family = "sans", size = 10.5, color = "#555555"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

pdf_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.1.pdf")
png_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.1.png")
svg_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.1.svg")

grDevices::pdf(pdf_path, width = 15, height = 12.5, useDingbats = FALSE, bg = "white")
print(figure5)
grDevices::dev.off()
ggsave(png_path, figure5, width = 15, height = 12.5, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(svg_path, figure5, width = 15, height = 12.5, units = "in", device = svglite, bg = "white")
cat("Figure files written:\n", pdf_path, "\n", png_path, "\n", svg_path, "\n")

source_paths <- c(panelA_path, panelB_path, panelC_path, panelD_path, panelE_path, panelF_path)
source_tables <- c(
  paths$section_qc,
  paste(paths$spot_scores, image_path, scale_path, sep = ";"),
  paths$consistency,
  paths$within_section,
  paste(paths$twas_boundary, paths$module_manifest, sep = ";"),
  paste(paths$claim_decision, claim_lock_path, sep = ";")
)
manifest <- rbindlist(lapply(seq_along(source_paths), function(i) {
  d <- fread(source_paths[i])
  data.table(
    figure = "Figure 5 draft v0.1", panel = LETTERS[i], source_table = source_tables[i],
    n_rows = nrow(d), n_columns = ncol(d), ready_for_publication_source_data = "yes_draft",
    notes = c(
      "Five-section QC and no-ROI boundary; representative selection rule recorded.",
      "QC-selected representative section with two spot-level overlays; visualization only; image and scale factor listed.",
      "Four primary modules x four key contexts; adjusted section-level evidence matrix.",
      "All raw/adjusted module-section pairs for four key contexts; no spot-level P values plotted.",
      "TWAS proxy counts plus R2/R5 spatial-module status.",
      "Required claim strip and locked spatial/TWAS wording."
    )[i]
  )
}))
manifest_path <- file.path(table_dir, "figure5_source_data_manifest_v0.1.tsv")
fwrite(manifest, manifest_path, sep = "\t")

visual_qc <- data.table(
  figure_id = "Figure 5", version = "v0.1",
  pdf_exists = file.exists(pdf_path), png_exists = file.exists(png_path), svg_exists = file.exists(svg_path),
  png_intended_dpi = 600, minimum_configured_font_pt = 8.0,
  panel_labels = "A-F present", legend_placement = "outside dense panels",
  palette_check = "pass_semantic_teal_bluegrey_amber_terracotta",
  claim_boundary_check = "pass_no_ROI_no_spatial_validation_TWAS_proxy",
  source_data_check = ifelse(all(file.exists(source_paths)), "pass", "fail"),
  legend_file = file.path(doc_dir, "draft_figure5_legend_v0.1.md"),
  legend_exists = file.exists(file.path(doc_dir, "draft_figure5_legend_v0.1.md")),
  rasterization_status = "controlled_600dpi_export; representative spot overlay rasterized only in PNG, vector in PDF/SVG",
  visual_status = "requires_rendered_agent_and_human_review",
  action_required = "Inspect PNG and independently rendered PDF at full size and 50% scale."
)
fwrite(visual_qc, file.path(table_dir, "figure5_visual_qc_v0.1.tsv"), sep = "\t")

cat("Claim lock:", claim_lock_path, "\n")
cat("Source manifest:", manifest_path, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
