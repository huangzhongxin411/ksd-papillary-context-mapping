#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
  library(png)
  library(ragg)
  library(svglite)
})

stage_id <- "stage7E_figure5_evidence_weighted_v0.4"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
in_dir <- "results/tables/revision/stage6C_spatial_twas_figure5_draft"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "stage7E_figure5_evidence_weighted_v0.4.log")
sink(log_path, split = TRUE)
on.exit({
  cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  sink()
}, add = TRUE)

cat("Stage:", stage_id, "\n")
cat("Rule: frozen Stage 6C panel sources only; no analysis rerun.\n")

read_panel <- function(id) {
  path <- file.path(in_dir, paste0("figure5_panel", id, "_source.tsv"))
  if (!file.exists(path)) stop("Missing frozen source: ", path)
  fread(path)
}

panelA <- read_panel("A")
panelB <- read_panel("B")
panelC <- read_panel("C")
panelD <- read_panel("D")
panelE <- read_panel("E")
panelF_source <- read_panel("F")
claim_lock_path <- file.path(in_dir, "spatial_twas_final_claim_lock_stage6C.tsv")
claim_lock <- fread(claim_lock_path)

stopifnot(
  nrow(panelA) == 5,
  sum(panelA$n_tissue_spots) == 7747,
  sum(panelA$representative_for_display) == 1,
  unique(panelA$section_id[panelA$representative_for_display]) == "GSM6250309",
  uniqueN(panelB$section_id) == 1,
  uniqueN(panelB$module_name) == 2,
  nrow(panelC) == 16,
  uniqueN(panelC$module_name) == 4,
  uniqueN(panelC$context_signature) == 4,
  nrow(panelE[component == "FDR-supported TWAS"]) == 2,
  sum(panelE[component == "FDR-supported TWAS"]$value) == 51
)

font_family <- "Helvetica"
ink <- "#252A2D"
muted <- "#667176"
teal <- "#2C6872"
teal_mid <- "#7E9FA6"
teal_light <- "#E5EEF0"
gold <- "#B7954D"
gold_light <- "#F4EBD6"
brick <- "#A55F50"
brick_mid <- "#C88A78"
brick_light <- "#F3E2DD"
line_col <- "#D5DADC"

theme_pub <- function(base_size = 9.1) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = ink),
      plot.title = element_text(face = "bold", size = base_size + 1.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.1, color = muted, lineheight = 0.95,
                                   margin = margin(b = 5)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, color = ink),
      axis.line = element_line(linewidth = 0.45, color = ink),
      axis.ticks = element_line(linewidth = 0.4, color = ink),
      legend.title = element_text(face = "bold", size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      plot.margin = margin(6, 7, 6, 7)
    )
}

tag_panel <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(family = font_family, face = "bold", size = 14, color = ink),
          plot.tag.position = c(0, 1),
          plot.title = element_text(margin = margin(l = 20, b = 2)))
}

module_short <- c(
  "R1_MAGMA_Bonferroni_only" = "R1",
  "R1_R2_R3_all_MAGMA_Bonferroni" = "R1-R3",
  "MAGMA_top50" = "Top 50",
  "MAGMA_top100" = "Top 100"
)
context_short <- c(
  "LoopTAL_signature" = "Loop / TAL",
  "injury_epithelial" = "Injury",
  "ECM_fibrosis" = "ECM / fibrosis",
  "mineralization_remodeling" = "Mineralization / remodeling"
)

# A. Resource plus hard no-ROI boundary.
panelA[, section_label := factor(section_id, levels = rev(section_id))]
pA <- ggplot(panelA, aes(y = section_label, x = n_tissue_spots, fill = representative_for_display)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.4) +
  geom_text(aes(label = comma(n_tissue_spots)), hjust = -0.12, size = 3.0,
            family = font_family, color = ink) +
  annotate("rect", xmin = 1940, xmax = 2980, ymin = 5.52, ymax = 5.92,
           fill = brick_light, color = NA) +
  annotate("text", x = 2460, y = 5.72, label = "NO PLAQUE / MINERAL / LESION ROI",
           family = font_family, fontface = "bold", size = 2.70, color = brick) +
  scale_fill_manual(values = c("TRUE" = gold, "FALSE" = teal_mid), guide = "none") +
  coord_cartesian(xlim = c(0, 3050), ylim = c(0.5, 6.05), clip = "off") +
  labs(title = "Spatial resource and no-ROI boundary",
       subtitle = "5 QC-passed sections · 7,747 tissue spots",
       caption = "GSM6250309: QC-selected display section, not correlation-selected\nPapillary tissue-context projection only",
       x = "Tissue spots", y = NULL) +
  theme_pub(8.7) +
  theme(plot.caption = element_text(size = 7.8, color = muted, hjust = 0, lineheight = 0.95),
        plot.margin = margin(6, 14, 5, 7))

# B. Representative histology-aligned overlays; visualization only.
image_path <- unique(panelB$histology_image)
if (length(image_path) != 1 || !file.exists(image_path)) stop("Panel B histology image unavailable")
img <- png::readPNG(image_path)
image_width <- dim(img)[2]
image_height <- dim(img)[1]
if (length(dim(img)) == 3L) img[, , 1:3] <- 0.76 * img[, , 1:3] + 0.24
limits_B <- as.numeric(quantile(panelB$module_score, c(0.02, 0.98), na.rm = TRUE))
panelB[, module_display := factor(module_label, levels = c("MAGMA top 50", "Loop / TAL signature"))]
pB <- ggplot(panelB, aes(x = x_hires, y = y_hires_plot)) +
  annotation_raster(img, xmin = 0, xmax = image_width, ymin = 0, ymax = image_height) +
  geom_point(aes(color = module_score), size = 0.65, alpha = 0.75) +
  facet_wrap(~module_display, nrow = 1) +
  scale_color_gradient2(low = "#9AAFB4", mid = "#F3F4F4", high = "#C39A90", midpoint = 0,
                        limits = limits_B, oob = squish, name = "Mean-z\nscore") +
  coord_fixed(xlim = c(0, image_width), ylim = c(0, image_height), expand = FALSE) +
  labs(title = "Representative papillary projection",
       subtitle = "GSM6250309 · Visualization only; not lesion localization", x = NULL, y = NULL) +
  theme_void(base_family = font_family, base_size = 8.6) +
  theme(
    plot.title = element_text(face = "bold", size = 10.3, color = ink, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 8.2, color = muted, margin = margin(b = 4)),
    strip.text = element_text(face = "bold", size = 8.5, color = ink),
    legend.position = "right", legend.title = element_text(face = "bold", size = 8.0),
    legend.text = element_text(size = 7.7), legend.key.height = unit(0.48, "cm"),
    plot.margin = margin(6, 6, 6, 7)
  )

# C. Main adjusted section-level evidence matrix.
panelC[, module_display := factor(module_short[module_name], levels = rev(unname(module_short)))]
panelC[, context_display := factor(context_short[context_signature], levels = unname(context_short))]
pC <- ggplot(panelC, aes(x = context_display, y = module_display,
                         fill = median_spearman_residualized)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label), family = font_family, size = 3.25,
            lineheight = 0.9, color = ink) +
  scale_fill_gradient2(low = "#B9877A", mid = "#F4F4F3", high = teal,
                       midpoint = 0, limits = c(-0.15, 0.15), oob = squish,
                       name = "Adjusted\nmedian rho") +
  labs(title = "Adjusted section-level co-distribution",
       subtitle = "Five-section median after complexity adjustment",
       caption = "CP = context-positive  ·  ATT = attenuated  ·  NS = not supported  ·  MIX = mixed/descriptive",
       x = NULL, y = NULL) +
  theme_pub(8.8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "right", legend.key.height = unit(0.55, "cm"),
        plot.caption = element_text(size = 7.6, color = muted, hjust = 0, margin = margin(t = 3)),
        plot.margin = margin(6, 8, 5, 7))

# D. Raw-to-adjusted changes for section x module pairs.
panelD[, module_display := unname(module_short[module_name])]
panelD[, context_display := factor(context_short[context_signature], levels = unname(context_short))]
panelD[, pair_id := paste(section_id, module_name, sep = ":")]
panelD_long <- melt(
  panelD,
  id.vars = c("section_id", "module_name", "module_display", "context_signature", "context_display", "pair_id"),
  measure.vars = c("spearman_r_raw", "spearman_r_residualized"),
  variable.name = "score_type", value.name = "spearman_r"
)
panelD_long[, score_type := factor(score_type,
  levels = c("spearman_r_raw", "spearman_r_residualized"),
  labels = c("Raw", "Adjusted"))]
panelD_median <- panelD_long[, .(median_rho = median(spearman_r)), by = .(context_display, score_type)]
pD <- ggplot(panelD_long, aes(x = score_type, y = spearman_r, group = pair_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#A3A8AA", linewidth = 0.38) +
  geom_line(color = "#AEB7BA", alpha = 0.20, linewidth = 0.30) +
  geom_point(color = teal_mid, size = 1.15, alpha = 0.43) +
  geom_line(data = panelD_median, aes(x = score_type, y = median_rho, group = 1),
            inherit.aes = FALSE, color = ink, linewidth = 1.20) +
  geom_point(data = panelD_median, aes(x = score_type, y = median_rho),
             inherit.aes = FALSE, color = ink, fill = "white", shape = 21,
             size = 2.8, stroke = 0.75) +
  facet_wrap(~context_display, nrow = 1) +
  labs(title = "Complexity-adjustment changes",
       subtitle = "Section-module pairs; black line = median",
       x = NULL, y = "Within-section rho") +
  theme_pub(8.2) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.background = element_rect(fill = "#F3F5F5", color = NA),
        strip.text = element_text(face = "bold", size = 7.8),
        panel.spacing.x = unit(0.15, "cm"), plot.margin = margin(6, 7, 5, 7))

# E. Compact TWAS proxy burden and boundary card.
twas_counts <- panelE[component == "FDR-supported TWAS"]
one_n <- twas_counts[category == "one-SNP", value]
multi_n <- twas_counts[category == "multi-SNP", value]
pE <- ggplot() +
  annotate("rect", xmin = 0.03, xmax = 0.44, ymin = 0.10, ymax = 0.90,
           fill = "#F5F7F7", color = line_col, linewidth = 0.5) +
  annotate("rect", xmin = 0.47, xmax = 0.97, ymin = 0.10, ymax = 0.90,
           fill = "#FAF8F6", color = line_col, linewidth = 0.5) +
  annotate("text", x = 0.06, y = 0.81, label = "TWAS PROXY BURDEN", hjust = 0,
           family = font_family, fontface = "bold", size = 3.2, color = teal) +
  annotate("text", x = 0.06, y = 0.63, label = "51", hjust = 0,
           family = font_family, fontface = "bold", size = 7.0, color = ink) +
  annotate("text", x = 0.145, y = 0.63, label = "FDR-supported genes", hjust = 0,
           family = font_family, size = 3.15, color = ink) +
  annotate("rect", xmin = 0.06, xmax = 0.06 + 0.32 * one_n / 51,
           ymin = 0.36, ymax = 0.47, fill = brick_mid, color = "white", linewidth = 0.4) +
  annotate("rect", xmin = 0.06 + 0.32 * one_n / 51, xmax = 0.38,
           ymin = 0.36, ymax = 0.47, fill = teal_mid, color = "white", linewidth = 0.4) +
  annotate("text", x = 0.06 + 0.16 * one_n / 51, y = 0.415, label = "42",
           family = font_family, fontface = "bold", size = 3.0, color = ink) +
  annotate("text", x = 0.06 + 0.32 * one_n / 51 + 0.16 * multi_n / 51,
           y = 0.415, label = "9", family = font_family, fontface = "bold",
           size = 3.0, color = ink) +
  annotate("text", x = 0.06, y = 0.27,
           label = "42 one-SNP (weaker proxy)  ·  9 multi-SNP",
           hjust = 0, family = font_family, size = 2.9, color = muted) +
  annotate("text", x = 0.50, y = 0.81, label = "BOUNDARY NOTES", hjust = 0,
           family = font_family, fontface = "bold", size = 3.2, color = brick) +
  annotate("rect", xmin = 0.50, xmax = 0.94, ymin = 0.66, ymax = 0.74,
           fill = "#F4E9E5", color = NA) +
  annotate("rect", xmin = 0.50, xmax = 0.94, ymin = 0.54, ymax = 0.62,
           fill = "#F6EFEC", color = NA) +
  annotate("rect", xmin = 0.50, xmax = 0.94, ymin = 0.42, ymax = 0.50,
           fill = "#F6EFEC", color = NA) +
  annotate("rect", xmin = 0.50, xmax = 0.94, ymin = 0.30, ymax = 0.38,
           fill = "#F6EFEC", color = NA) +
  annotate("rect", xmin = 0.50, xmax = 0.94, ymin = 0.18, ymax = 0.26,
           fill = "#F4E9E5", color = NA) +
  annotate("text", x = 0.52, y = 0.70, label = "GTEx kidney cortex proxy only",
           hjust = 0, family = font_family, size = 2.75, color = ink) +
  annotate("text", x = 0.52, y = 0.58, label = "R2 descriptive",
           hjust = 0, family = font_family, size = 2.75, color = ink) +
  annotate("text", x = 0.52, y = 0.46, label = "R5 not score-feasible",
           hjust = 0, family = font_family, size = 2.75, color = ink) +
  annotate("text", x = 0.52, y = 0.34, label = "No papilla-specific eQTL",
           hjust = 0, family = font_family, size = 2.75, color = ink) +
  annotate("text", x = 0.52, y = 0.22, label = "No SMR/coloc support",
           hjust = 0, family = font_family, size = 2.75, color = ink) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "TWAS proxy boundary",
       subtitle = "Supplementary GTEx kidney cortex context only") +
  theme_void(base_family = font_family, base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10.5, color = ink, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8.5, color = muted, margin = margin(b = 3)),
        plot.margin = margin(6, 8, 5, 7))

# F. Two-zone interpretation boundary banner.
pF <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.48, ymin = 0.12, ymax = 0.88,
           fill = teal_light, color = "white", linewidth = 0.7) +
  annotate("rect", xmin = 0.52, xmax = 0.98, ymin = 0.12, ymax = 0.88,
           fill = brick_light, color = "white", linewidth = 0.7) +
  annotate("text", x = 0.25, y = 0.76, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 3.0, color = teal) +
  annotate("text", x = 0.75, y = 0.76, label = "NOT CLAIMED",
           family = font_family, fontface = "bold", size = 3.0, color = brick) +
  annotate("text", x = 0.25, y = 0.38,
           label = "Moderate supplementary spatial context\nLoop/TAL-associated tissue-context projection\nGTEx kidney cortex proxy context",
           family = font_family, fontface = "bold", size = 2.75,
           lineheight = 1.22, color = ink) +
  annotate("text", x = 0.75, y = 0.38,
           label = "No plaque ROI\nNot plaque-specific validation\nNot causal spatial niche\nNot plaque nucleation site",
           family = font_family, fontface = "bold", size = 2.75,
           lineheight = 1.15, color = ink) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Interpretation boundary") +
  theme_void(base_family = font_family, base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10.5, color = ink, margin = margin(b = 0)),
        plot.margin = margin(5, 8, 4, 8))

figure5 <- (
  tag_panel(pA, "A") | tag_panel(pB, "B")
) / (
  tag_panel(pC, "C") | tag_panel(pD, "D")
) / tag_panel(pE, "E") / tag_panel(pF, "F") +
  plot_layout(heights = c(0.70, 1.40, 0.68, 0.43), widths = c(0.96, 1.04)) +
  plot_annotation(
    title = "Spatial and TWAS layers provide supplementary papillary context",
    subtitle = "Loop/TAL co-distribution is retained or attenuated after complexity adjustment; no plaque ROI is available",
    theme = theme(
      text = element_text(family = font_family, color = ink),
      plot.title = element_text(face = "bold", size = 16.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10.1, color = muted, margin = margin(b = 4)),
      plot.margin = margin(8, 8, 7, 8)
    )
  )

pdf_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.4.pdf")
png_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.4.png")
svg_path <- file.path(fig_dir, "figure5_spatial_twas_context_draft_v0.4.svg")

grDevices::pdf(pdf_path, width = 14.5, height = 11.0, family = font_family,
               useDingbats = FALSE, bg = "white")
print(figure5)
grDevices::dev.off()
ggsave(png_path, figure5, width = 14.5, height = 11.0, units = "in", dpi = 600,
       device = ragg::agg_png, bg = "white")
ggsave(svg_path, figure5, width = 14.5, height = 11.0, units = "in",
       device = svglite::svglite, bg = "white")

source_files <- list(
  file.path(in_dir, "figure5_panelA_source.tsv"),
  file.path(in_dir, "figure5_panelB_source.tsv"),
  file.path(in_dir, "figure5_panelC_source.tsv"),
  file.path(in_dir, "figure5_panelD_source.tsv"),
  file.path(in_dir, "figure5_panelE_source.tsv"),
  file.path(in_dir, c("figure5_panelF_source.tsv", "spatial_twas_final_claim_lock_stage6C.tsv"))
)
manifest <- data.table(
  panel = LETTERS[1:6],
  source_data_file = vapply(source_files, paste, collapse = ";", FUN.VALUE = character(1)),
  source_exists = vapply(source_files, function(x) ifelse(all(file.exists(x)), "yes", "no"), FUN.VALUE = character(1)),
  row_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) nrow(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  column_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) ncol(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  data_changed_from_v0.3 = "no",
  claim_changed_from_v0.3 = "no",
  ready_for_publication_source_data = "yes",
  notes = c(
    "Five-section QC, section spot counts, representative-selection rule, and no-ROI boundary.",
    "Frozen representative-section spot scores and histology coordinates; visualization only.",
    "Four modules x four contexts with adjusted median rho and locked status labels.",
    "All raw/adjusted section-module pairs; no spot-level P values.",
    "TWAS proxy counts and R2/R5 feasibility boundaries.",
    "Stage 6C boundary source plus locked claim table; wording condensed without claim upgrade."
  )
)
fwrite(manifest, file.path(table_dir, "figure5_source_data_manifest_v0.4.tsv"), sep = "\t")

changes <- data.table(
  panel = c("overall", LETTERS[1:6]),
  issue_in_v0.3 = c(
    "Evidence hierarchy remained too even",
    "Resource band remained relatively tall",
    "Representative image retained too much visual weight",
    "Main adjusted matrix needed greater central weight",
    "Adjustment trajectories needed greater central weight",
    "TWAS panel still read as two text boxes",
    "Interpretation boundary still resembled a checklist"
  ),
  redesign_made_in_v0.4 = c(
    "Reweighted four bands around dominant C/D evidence row",
    "Compressed first band and retained compact no-ROI badge",
    "Reduced band height and further subdued overlay palette",
    "Enlarged central heatmap and cell labels",
    "Enlarged trajectory panel and strengthened medians",
    "Rebuilt as visual burden card plus five compact boundary labels",
    "Rebuilt as two large supported/not-claimed zones"
  ),
  data_changed = "no",
  claim_changed = "no",
  reason = c("Improve scientific hierarchy", "Make resource boundary immediate", "Prevent illustrative overlay overreading", "Make adjusted spatial evidence primary", "Clarify attenuation without spot-level inference", "Downgrade one-SNP-heavy TWAS appropriately", "Make claim ceiling explicit"),
  human_review_needed = "yes",
  notes = c("Six panels preserved", "All five sections and counts retained", "Same section, scores, coordinates, and display-scale rule", "All 16 values/status labels retained", "All frozen raw/adjusted pairs retained", "Counts 51/42/9 and R2/R5 boundaries retained", "No plaque ROI, validation, causal niche, or nucleation claim")
)
fwrite(changes, file.path(table_dir, "figure5_visual_redesign_log_v0.4.tsv"), sep = "\t")

cat("Wrote:", pdf_path, "\n")
cat("Wrote:", png_path, "\n")
cat("Wrote:", svg_path, "\n")
cat("Rows: A=", nrow(panelA), ", B=", nrow(panelB), ", C=", nrow(panelC),
    ", D=", nrow(panelD), ", E=", nrow(panelE), "\n", sep = "")
cat("data_changed=no; claim_changed=no\n")
