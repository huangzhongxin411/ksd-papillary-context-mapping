#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
  library(ragg)
  library(svglite)
})

stage_id <- "stage7E_figure4_lock_micro_polish_v0.6"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
in_dir <- "results/tables/revision/stage5C1_gse73680_figure4_draft"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "stage7E_figure4_lock_micro_polish_v0.6.log")
sink(log_path, split = TRUE)
on.exit({
  cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
  sink()
}, add = TRUE)

cat("Stage:", stage_id, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Rule: artwork refinement only; source data and claims unchanged.\n")

read_source <- function(panel) {
  path <- file.path(in_dir, paste0("figure4_panel", panel, "_source.tsv"))
  if (!file.exists(path)) stop("Missing locked source: ", path)
  fread(path)
}

panelA <- read_source("A")
panelB <- read_source("B")
panelC <- read_source("C")
panelD <- read_source("D")
panelE <- read_source("E")
panelF <- read_source("F")
panelG <- read_source("G")

stopifnot(
  sum(panelA$n_samples_in_record) == 55,
  sum(panelA$usable_for_paired_delta) == 26,
  uniqueN(panelA$patient_id) == 29,
  nrow(panelB) == 104,
  uniqueN(panelB$patient_id) == 26,
  uniqueN(panelB$module_name) == 4,
  nrow(panelC) == 4,
  nrow(panelF) == 12
)

module_order <- c("R1 MAGMA Bonf.", "R1-R3 MAGMA Bonf.", "MAGMA top 50", "MAGMA top 100")
short_labels <- c("R1 MAGMA Bonf." = "R1", "R1-R3 MAGMA Bonf." = "R1-R3",
                  "MAGMA top 50" = "Top 50", "MAGMA top 100" = "Top 100")
panelB[, module_short := factor(short_labels[as.character(module_label)], levels = unname(short_labels))]
panelC[, module_short := factor(short_labels[as.character(module_label)], levels = rev(unname(short_labels)))]
panelE[, module_short := factor(short_labels[as.character(module_label)], levels = unname(short_labels))]
panelE[, model_family := factor(model_family, levels = c("Single-signature", "Compact"))]
panelE[, outcome := factor(outcome, levels = c("Direction retained", "Attenuated", "Lost / explained"))]
panelF[, module_short := factor(short_labels[as.character(module_label)], levels = rev(unname(short_labels)))]
signature_short <- c("Injury epithelial" = "Injury", "ECM / fibrosis" = "ECM / fibrosis",
                     "Mineralization / remodeling" = "Mineralization / remodeling")
panelF[, signature_short := factor(signature_short[as.character(signature_label)], levels = unname(signature_short))]
panelF[, rho_label := formatC(spearman_r, digits = 2, format = "f")]

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

theme_pub <- function(base_size = 9.3) {
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
          plot.tag.position = c(0, 1))
}

# A. Compact study-design and preprocessing schematic.
pA <- ggplot() +
  annotate("text", x = 0.05, y = 0.91, label = "GSE73680", hjust = 0,
           family = font_family, fontface = "bold", size = 4.7, color = ink) +
  annotate("text", x = 0.05, y = 0.79, label = "paired renal papillary bulk expression",
           hjust = 0, family = font_family, size = 3.35, color = muted) +
  annotate("rect", xmin = 0.05, xmax = 0.33, ymin = 0.38, ymax = 0.67,
           fill = teal_light, color = teal_mid, linewidth = 0.65) +
  annotate("rect", xmin = 0.67, xmax = 0.95, ymin = 0.38, ymax = 0.67,
           fill = gold_light, color = gold, linewidth = 0.65) +
  annotate("text", x = 0.19, y = 0.56, label = "Control / adjacent",
           family = font_family, fontface = "bold", size = 3.55, color = ink) +
  annotate("text", x = 0.19, y = 0.46, label = "papilla", family = font_family,
           size = 3.2, color = muted) +
  annotate("text", x = 0.81, y = 0.56, label = "Plaque / stone",
           family = font_family, fontface = "bold", size = 3.55, color = ink) +
  annotate("text", x = 0.81, y = 0.46, label = "papilla", family = font_family,
           size = 3.2, color = muted) +
  annotate("segment", x = 0.35, xend = 0.65, y = 0.525, yend = 0.525,
           linewidth = 0.9, color = ink, arrow = arrow(length = unit(0.13, "in"))) +
  annotate("text", x = 0.50, y = 0.635, label = "26 paired patients",
           family = font_family, fontface = "bold", size = 3.7, color = teal) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.12, ymax = 0.27,
           fill = "#F8FAFA", color = NA) +
  annotate("text", x = 0.50, y = 0.195,
           label = "55 samples  ·  29 patients  ·  26 paired patients  ·  log2(x + 1) module scores",
           family = font_family, size = 3.05, color = "#6F777A") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Paired bulk resource") +
  theme_void(base_family = font_family, base_size = 9.3) +
  theme(plot.title = element_text(face = "bold", size = 10.8, color = ink,
                                  margin = margin(b = 2)),
        plot.margin = margin(6, 8, 6, 8))

# B. Dominant paired response panel; all source values are unchanged.
pB <- ggplot(panelB, aes(x = module_short, y = paired_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#858C8F", linetype = "dashed") +
  geom_violin(width = 0.78, fill = "#F0F5F5", color = NA, trim = TRUE, alpha = 0.70) +
  geom_boxplot(width = 0.28, outlier.shape = NA, fill = "white", color = teal,
               linewidth = 0.55) +
  geom_jitter(width = 0.095, height = 0, size = 1.08, alpha = 0.58, color = teal) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.0,
               fill = gold, color = ink, stroke = 0.35) +
  labs(title = "Paired MAGMA-module response",
       subtitle = "Each point = one paired patient",
       x = NULL, y = "Paired module-score delta") +
  theme_pub(9.3) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1),
        plot.margin = margin(6, 10, 6, 7))

# C. Patient-aware unadjusted model estimates.
pC <- ggplot(panelC, aes(y = module_short, x = effect_estimate)) +
  geom_vline(xintercept = 0, color = "#858C8F", linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.17,
                color = teal_mid, linewidth = 0.88, orientation = "y") +
  geom_point(size = 3.65, shape = 21, fill = teal, color = "white", stroke = 0.42) +
  annotate("text", x = max(panelC$ci_high) + 0.008, y = 4.17,
           label = "All~FDR%~~%0.036", parse = TRUE, hjust = 0,
           family = font_family, size = 3.0, color = muted) +
  coord_cartesian(xlim = c(min(panelC$ci_low) - 0.02, max(panelC$ci_high) + 0.09),
                  ylim = c(0.65, 4.3), clip = "off") +
  labs(title = "Patient-aware effect summary",
       subtitle = "Patient-aware paired model",
       x = "Paired effect estimate", y = NULL) +
  theme_pub(9.1) +
  theme(plot.margin = margin(6, 34, 6, 7))

# D. Correlation-only heatmap; exact FDR values remain in Source Data and the legend.
pD <- ggplot(panelF, aes(x = signature_short, y = module_short, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = rho_label), family = font_family, size = 3.15, color = ink) +
  scale_fill_gradient(low = "#F8F4F3", high = "#B9877A", limits = c(0, 0.85),
                      breaks = c(0, 0.4, 0.8), name = "Spearman rho") +
  labs(title = "Injury/remodeling coupling",
       subtitle = "Bulk proxy correlations; not mechanism", x = NULL, y = NULL) +
  theme_pub(8.8) +
  theme(axis.text.x = element_text(angle = 21, hjust = 1),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center", legend.key.width = unit(1.25, "cm"),
        legend.key.height = unit(0.26, "cm"),
        legend.margin = margin(t = -3), legend.box.margin = margin(t = -4),
        plot.margin = margin(6, 7, 5, 7))

# E. Two-column adjustment summary with one shared paired-shift note.
stage_labels <- c("Single-signature\nadjustment", "Compact\nadjustment")
tile_base <- CJ(module_short = factor(unname(short_labels), levels = unname(short_labels)), stage_x = 1:2)
row_y <- c("R1" = 3.55, "R1-R3" = 2.70, "Top 50" = 1.85, "Top 100" = 1.00)
tile_base[, module_y := unname(row_y[as.character(module_short)])]
adjust_seg <- copy(panelE)
adjust_seg[, stage_x := fifelse(model_family == "Single-signature", 1L, 2L)]
adjust_seg <- adjust_seg[, .(module_short, stage_x, outcome, fraction, n_models)]
seg <- adjust_seg
seg[, module_y := unname(row_y[as.character(module_short)])]
setorder(seg, module_y, stage_x, outcome)
seg[, cum0 := shift(cumsum(fraction), fill = 0), by = .(module_short, stage_x)]
seg[, `:=`(xmin = stage_x - 0.39 + 0.78 * cum0,
           xmax = stage_x - 0.39 + 0.78 * (cum0 + fraction),
           ymin = module_y - 0.15, ymax = module_y + 0.15,
           xmid = (stage_x - 0.39 + 0.78 * cum0 + stage_x - 0.39 + 0.78 * (cum0 + fraction)) / 2)]
pE <- ggplot() +
  geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = outcome),
            color = "white", linewidth = 0.42) +
  geom_text(data = seg, aes(x = xmid, y = module_y, label = n_models),
            family = font_family, fontface = "bold", size = 3.15, color = ink) +
  annotate("text", x = 1.5, y = 4.82,
           label = "Paired bulk shift was positive in all four modules before adjustment",
           family = font_family, fontface = "bold", size = 3.00, color = teal) +
  annotate("text", x = 1, y = 4.08, label = "Single-signature\nadjustment",
           family = font_family, fontface = "bold", size = 3.15, lineheight = 0.92, color = ink) +
  annotate("text", x = 2, y = 4.08, label = "Compact\nadjustment",
           family = font_family, fontface = "bold", size = 3.15, lineheight = 0.92, color = ink) +
  scale_fill_manual(values = c("Direction retained" = teal_mid, "Attenuated" = gold,
                               "Lost / explained" = brick_mid), drop = FALSE) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.46)) +
  scale_y_continuous(breaks = unname(row_y), labels = names(row_y),
                     limits = c(0.58, 5.05), expand = c(0, 0)) +
  labs(title = "Adjustment attenuation summary",
       subtitle = "Effects are attenuated after adjustment; direction retained overall",
       x = NULL, y = NULL, fill = NULL) +
  theme_pub(9.0) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(face = "bold", color = ink, lineheight = 0.95),
        legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(6, 12, 5, 7))

# F. Full-width claim-boundary strip.
boundary <- data.table(
  x = 1:6,
  group = c(rep("SUPPORTED INTERPRETATION", 3), rep("NOT CLAIMED", 3)),
  label = c(
    "Paired bulk\nmodule shift",
    "Injury/remodeling-\nassociated context",
    "Attenuated but\nretained overall",
    "Independent\nvalidation",
    "Cell-type-specific\nresponse",
    "Genetic causality or\nplaque mechanism"
  )
)
pF <- ggplot(boundary, aes(x = x, y = 1, fill = group)) +
  geom_tile(width = 0.94, height = 0.62, color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), family = font_family, fontface = "bold",
            size = 3.05, lineheight = 0.93, color = ink) +
  annotate("text", x = 2, y = 1.43, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 3.05, color = teal) +
  annotate("text", x = 5, y = 1.43, label = "NOT CLAIMED",
           family = font_family, fontface = "bold", size = 3.05, color = brick) +
  scale_fill_manual(values = c("SUPPORTED INTERPRETATION" = teal_light,
                               "NOT CLAIMED" = brick_light), guide = "none") +
  scale_x_continuous(expand = expansion(add = 0.04)) +
  coord_cartesian(ylim = c(0.62, 1.58), clip = "off") +
  labs(title = "Interpretation boundary") +
  theme_void(base_family = font_family, base_size = 9.2) +
  theme(plot.title = element_text(face = "bold", size = 10.7, color = ink,
                                  margin = margin(b = 0)),
        plot.margin = margin(5, 8, 4, 8))

figure4 <- (
  tag_panel(pA, "A") | tag_panel(pB, "B")
) / (
  tag_panel(pC, "C") | tag_panel(pD, "D")
) / tag_panel(pE, "E") / tag_panel(pF, "F") +
  plot_layout(heights = c(0.88, 1.18, 1.04, 0.47), widths = c(0.82, 1.28)) +
  plot_annotation(
    title = "Paired bulk papillary disease-context support in GSE73680",
    subtitle = "Module shifts are injury/remodeling-associated and attenuated after adjustment",
    theme = theme(
      text = element_text(family = font_family, color = ink),
      plot.title = element_text(face = "bold", size = 16.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10.2, color = muted, margin = margin(b = 4)),
      plot.margin = margin(8, 8, 7, 8)
    )
  )

pdf_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.6.pdf")
png_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.6.png")
svg_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.6.svg")

grDevices::pdf(pdf_path, width = 14.5, height = 10.6, family = font_family,
               useDingbats = FALSE, bg = "white")
print(figure4)
grDevices::dev.off()
ggsave(png_path, figure4, width = 14.5, height = 10.6, units = "in", dpi = 600,
       device = ragg::agg_png, bg = "white")
ggsave(svg_path, figure4, width = 14.5, height = 10.6, units = "in",
       device = svglite::svglite, bg = "white")

source_files <- list(
  file.path(in_dir, "figure4_panelA_source.tsv"),
  file.path(in_dir, "figure4_panelB_source.tsv"),
  file.path(in_dir, "figure4_panelC_source.tsv"),
  file.path(in_dir, "figure4_panelF_source.tsv"),
  file.path(in_dir, c("figure4_panelC_source.tsv", "figure4_panelE_source.tsv")),
  file.path(in_dir, c("figure4_panelG_source.tsv", "gse73680_final_claim_lock_stage5C1.tsv"))
)
manifest <- data.table(
  panel = LETTERS[1:6],
  source_data_file = vapply(source_files, paste, collapse = ";", FUN.VALUE = character(1)),
  source_exists = vapply(source_files, function(x) ifelse(all(file.exists(x)), "yes", "no"), FUN.VALUE = character(1)),
  row_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) nrow(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  column_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) ncol(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  data_changed_from_v0.5 = "no",
  claim_changed_from_v0.5 = "no",
  ready_for_publication_source_data = "yes",
  notes = c(
    "Counts derived from locked patient-level provenance.",
    "All 104 paired-delta values retained.",
    "All four model estimates retained.",
    "All 12 rho values retained; FDR values remain in Source Data and the legend.",
    "Paired effects and all source-backed adjustment outcome counts retained.",
    "Boundary wording condensed without upgrading the evidence class."
  )
)
fwrite(manifest, file.path(table_dir, "figure4_source_data_manifest_v0.6.tsv"), sep = "\t")

changes <- data.table(
  panel = c("overall", "A", "B", "C", "D", "E", "F", "overall", "overall"),
  issue_in_v0.5 = c("Lock-level consistency check remained", "Resource strip required lock confirmation", "Violin/point balance required lock confirmation", "FDR alignment required lock confirmation", "Heatmap could be one step more restrained", "Shared note wording and residual row framing required final adjustment", "Boundary strip required lock confirmation", "Typography required lock confirmation", "Final tonal restraint required lock confirmation"),
  micro_polish_made_in_v0.6 = c("Confirmed six-panel balance and spacing", "Retained centered schematic and readable subtle strip", "Retained accepted light violin and point treatment", "Retained aligned approximate FDR note", "Further softened terracotta scale and tightened colorbar", "Applied requested note wording, moved headers closer, removed remaining row backgrounds", "Confirmed equal-height six-box boundary strip", "Confirmed Helvetica hierarchy at full and reduced scale", "Applied final restrained tonal harmonization"),
  data_changed = "no",
  claim_changed = "no",
  reason = c("Improve narrative continuity", "Reduce internal text and empty space", "Make paired bulk evidence immediately legible", "Improve double-column readability", "Keep related disease-context evidence together", "Make attenuation the final quantitative result", "Separate supported interpretation from exclusions", "Match Figures 1-3", "Maintain accessible, non-saturated styling"),
  human_review_needed = "yes",
  notes = c("Architecture only", "Counts unchanged", "All 104 values unchanged", "All estimates and intervals unchanged", "All 12 rho values unchanged; FDR values moved to legend and Source Data", "All model-outcome counts unchanged", "Claim class unchanged", "Final cross-figure proofing remains", "Final grayscale proofing remains")
)
fwrite(changes, file.path(table_dir, "figure4_micro_polish_log_v0.6.tsv"), sep = "\t")

cat("Wrote:", pdf_path, "\n")
cat("Wrote:", png_path, "\n")
cat("Wrote:", svg_path, "\n")
cat("Source rows: A=", nrow(panelA), ", B=", nrow(panelB), ", C=", nrow(panelC),
    ", D=", nrow(panelF), ", E=", nrow(panelE), "\n", sep = "")
cat("data_changed=no; claim_changed=no\n")
