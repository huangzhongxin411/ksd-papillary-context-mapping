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

stage_id <- "stage7E_figure4_refinement_v0.2"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
in_dir <- "results/tables/revision/stage5C1_gse73680_figure4_draft"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "stage7E_figure4_refinement_v0.2.log")
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
panelB[, module_label := factor(module_label, levels = module_order)]
panelC[, module_label := factor(module_label, levels = rev(module_order))]
panelE[, module_label := factor(module_label, levels = module_order)]
panelE[, model_family := factor(model_family, levels = c("Single-signature", "Compact"))]
panelE[, outcome := factor(outcome, levels = c("Direction retained", "Attenuated", "Lost / explained"))]
panelF[, module_label := factor(module_label, levels = rev(module_order))]
panelF[, signature_label := factor(signature_label,
  levels = c("Injury epithelial", "ECM / fibrosis", "Mineralization / remodeling"))]
panelF[, cell_label := paste0("rho ", formatC(spearman_r, digits = 2, format = "f"),
                              "\nFDR ", formatC(fdr, digits = 3, format = "g"))]

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
  annotate("segment", x = 0.35, xend = 0.65, y = 0.535, yend = 0.535,
           linewidth = 0.9, color = ink, arrow = arrow(length = unit(0.13, "in"))) +
  annotate("text", x = 0.50, y = 0.66, label = "26 paired patients",
           family = font_family, fontface = "bold", size = 3.7, color = teal) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.12, ymax = 0.27,
           fill = "#F5F7F7", color = NA) +
  annotate("text", x = 0.50, y = 0.195,
           label = "55 samples  |  29 patients  |  log2(x + 1) module scores",
           family = font_family, fontface = "bold", size = 3.2, color = ink) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Paired bulk disease-context resource") +
  theme_void(base_family = font_family, base_size = 9.3) +
  theme(plot.title = element_text(face = "bold", size = 10.8, color = ink,
                                  margin = margin(b = 2)),
        plot.margin = margin(6, 8, 6, 8))

# B. Dominant paired response panel; all source values are unchanged.
pB <- ggplot(panelB, aes(x = module_label, y = paired_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#858C8F", linetype = "dashed") +
  geom_violin(width = 0.78, fill = teal_light, color = NA, trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.28, outlier.shape = NA, fill = "white", color = teal,
               linewidth = 0.55) +
  geom_jitter(width = 0.095, height = 0, size = 1.35, alpha = 0.62, color = teal) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.0,
               fill = gold, color = ink, stroke = 0.35) +
  labs(title = "Paired MAGMA-module response",
       subtitle = "Paired plaque/stone bulk shift; each point is one patient (n = 26)",
       x = NULL, y = "Paired module-score delta") +
  theme_pub(9.3) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1),
        plot.margin = margin(6, 10, 6, 7))

# C. Patient-aware unadjusted model estimates.
pC <- ggplot(panelC, aes(y = module_label, x = effect_estimate)) +
  geom_vline(xintercept = 0, color = "#858C8F", linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.16,
                color = teal_mid, linewidth = 0.75, orientation = "y") +
  geom_point(size = 3.35, shape = 21, fill = teal, color = "white", stroke = 0.4) +
  geom_text(aes(x = ci_high + 0.014, label = fdr_label), hjust = 0,
            family = font_family, size = 3.0, color = muted) +
  coord_cartesian(xlim = c(min(panelC$ci_low) - 0.02, max(panelC$ci_high) + 0.09), clip = "off") +
  labs(title = "Patient-aware effect summary",
       subtitle = "Patient-aware paired model of the bulk module shift",
       x = "Group effect estimate (95% CI)", y = NULL) +
  theme_pub(9.1) +
  theme(plot.margin = margin(6, 34, 6, 7))

# D. Injury/remodeling coupling, with the source-backed paired signature shifts as a compact header.
injury_summary <- panelD[signature_name %in% c("injury_epithelial", "ECM_fibrosis", "mineralization_remodeling")]
injury_summary[, signature_short := factor(signature_label,
  levels = c("Injury epithelial", "ECM / fibrosis", "Mineralization / remodeling"))]
pD_top <- ggplot(injury_summary, aes(x = signature_short, y = mean_delta)) +
  geom_hline(yintercept = 0, color = "#858C8F", linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.58, fill = brick_mid) +
  geom_text(aes(label = sprintf("%+.2f", mean_delta)), vjust = -0.35,
            family = font_family, size = 2.75, color = ink) +
  scale_y_continuous(limits = c(0, 0.42), breaks = c(0, 0.2, 0.4), expand = expansion(mult = c(0, 0))) +
  labs(title = "Injury/remodeling signature coupling",
       subtitle = "Injury/remodeling-associated bulk proxies; correlation is not mechanism",
       tag = "D", y = NULL, x = NULL) +
  theme_pub(7.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10.6),
        plot.subtitle = element_text(size = 8.8),
        plot.tag = element_text(family = font_family, face = "bold", size = 14, color = ink),
        plot.tag.position = c(-0.12, 1.02),
        plot.margin = margin(5, 34, 0, 7))

pD_heat <- ggplot(panelF, aes(x = signature_label, y = module_label, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = cell_label), family = font_family, size = 2.65,
            lineheight = 0.9, color = ink) +
  scale_fill_gradient(low = "#F5F1F0", high = brick, limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "Spearman\nrho") +
  labs(x = NULL, y = NULL) +
  theme_pub(8.2) +
  theme(axis.text.x = element_text(angle = 21, hjust = 1),
        legend.position = "right", legend.key.height = unit(0.56, "cm"),
        plot.margin = margin(0, 7, 5, 7))

pD <- (pD_top / pD_heat) +
  plot_layout(heights = c(0.48, 1))

# E. Adjustment outcomes remain fully visible and explicitly attenuated.
pE <- ggplot(panelE, aes(x = module_label, y = fraction, fill = outcome)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.35) +
  geom_text(aes(label = n_models), position = position_stack(vjust = 0.5),
            family = font_family, size = 2.9, color = ink) +
  facet_wrap(~model_family, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c("Direction retained" = teal_mid,
                               "Attenuated" = gold,
                               "Lost / explained" = brick_mid), drop = FALSE) +
  labs(title = "Composition/injury/remodeling-aware sensitivity",
       subtitle = "Effects attenuated after adjustment; direction retained overall",
       x = NULL, y = "Sensitivity-model outcomes", fill = NULL) +
  theme_pub(8.9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom", legend.direction = "horizontal",
        strip.background = element_rect(fill = "#F2F4F4", color = NA),
        strip.text = element_text(face = "bold", size = 8.6),
        panel.spacing.x = unit(0.55, "cm"))

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
  tag_panel(pC, "C") | pD
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

pdf_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.2.pdf")
png_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.2.png")
svg_path <- file.path(fig_dir, "figure4_gse73680_bulk_context_draft_v0.2.svg")

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
  file.path(in_dir, c("figure4_panelD_source.tsv", "figure4_panelF_source.tsv")),
  file.path(in_dir, "figure4_panelE_source.tsv"),
  file.path(in_dir, c("figure4_panelG_source.tsv", "gse73680_final_claim_lock_stage5C1.tsv"))
)
manifest <- data.table(
  panel = LETTERS[1:6],
  source_data_file = vapply(source_files, paste, collapse = ";", FUN.VALUE = character(1)),
  source_exists = vapply(source_files, function(x) ifelse(all(file.exists(x)), "yes", "no"), FUN.VALUE = character(1)),
  row_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) nrow(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  column_count_if_readable = vapply(source_files, function(x) paste(vapply(x, function(y) ncol(fread(y)), integer(1)), collapse = ";"), FUN.VALUE = character(1)),
  data_changed_from_v0.1 = "no",
  claim_changed_from_v0.1 = "no",
  ready_for_publication_source_data = "yes",
  notes = c(
    "Counts derived from locked patient-level provenance.",
    "All 104 paired-delta values retained.",
    "All four model estimates retained.",
    "Marker signatures remain bulk proxies; correlations do not imply mechanism.",
    "All source-backed outcome counts retained.",
    "Boundary wording condensed without upgrading the evidence class."
  )
)
fwrite(manifest, file.path(table_dir, "figure4_source_data_manifest_v0.2.tsv"), sep = "\t")

changes <- data.table(
  panel = c("overall", "A", "B", "C", "D", "E", "F", "overall", "overall"),
  issue_in_v0.1 = c("Seven panels dispersed the evidence flow", "Study schematic used excess space", "Paired response had limited visual priority", "Forest plot had excess right margin", "Signature response and coupling were split", "Adjustment summary competed with other panels", "Coupling and boundary were separate panels", "Default sans typography", "Palette needed production consistency"),
  change_made_in_v0.2 = c("Consolidated to six panels A-F", "Compact design with preprocessing metric", "Dominant violin/box/point response", "Tightened coefficient display", "Combined injury/remodeling response and coupling", "Full-width two-facet adjustment summary", "Full-width claim-boundary strip", "Applied Helvetica hierarchy", "Retained restrained semantic teal/gold/terracotta"),
  data_changed = "no",
  claim_changed = "no",
  reason = c("Improve narrative continuity", "Reduce internal text and empty space", "Make paired bulk evidence immediately legible", "Improve double-column readability", "Keep related disease-context evidence together", "Make attenuation the final quantitative result", "Separate supported interpretation from exclusions", "Match Figures 1-3", "Maintain accessible, non-saturated styling"),
  human_review_needed = "yes",
  notes = c("Architecture only", "Counts unchanged", "All 104 values unchanged", "All estimates and intervals unchanged", "Only three injury/remodeling proxies shown above the unchanged 12-cell heatmap", "All model-outcome counts unchanged", "Claim class unchanged", "Final cross-figure proofing remains", "Final grayscale proofing remains")
)
fwrite(changes, file.path(table_dir, "figure4_visual_change_log_v0.2.tsv"), sep = "\t")

cat("Wrote:", pdf_path, "\n")
cat("Wrote:", png_path, "\n")
cat("Wrote:", svg_path, "\n")
cat("Source rows: A=", nrow(panelA), ", B=", nrow(panelB), ", C=", nrow(panelC),
    ", D=", nrow(injury_summary), "+", nrow(panelF), ", E=", nrow(panelE), "\n", sep = "")
cat("data_changed=no; claim_changed=no\n")
