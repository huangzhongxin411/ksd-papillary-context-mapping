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

stage_id <- "stage7E_figure4_design_handoff"
fig_dir <- file.path("results/figures/revision", stage_id)
table_dir <- file.path("results/tables/revision", stage_id)
doc_dir <- file.path("docs/revision", stage_id)
log_dir <- file.path("logs/revision", stage_id)
in_dir <- "results/tables/revision/stage5C1_gse73680_figure4_draft"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "stage7E_figure4_design_handoff.log")
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
           linewidth = 1.05, color = ink, arrow = arrow(length = unit(0.13, "in"))) +
  annotate("text", x = 0.50, y = 0.635, label = "26 paired patients",
           family = font_family, fontface = "bold", size = 3.7, color = teal) +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.12, ymax = 0.27,
           fill = "#F8FAFA", color = NA) +
  annotate("text", x = 0.50, y = 0.195,
           label = "55 samples  ·  29 patients  ·  26 paired patients  ·  log2(x + 1) module scores",
           family = font_family, size = 3.15, color = "#626C70") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Paired bulk resource") +
  theme_void(base_family = font_family, base_size = 9.3) +
  theme(plot.title = element_text(face = "bold", size = 10.8, color = ink,
                                  margin = margin(b = 2)),
        plot.margin = margin(6, 8, 6, 8))

# B. Dominant paired response panel; all source values are unchanged.
pB <- ggplot(panelB, aes(x = module_short, y = paired_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#858C8F", linetype = "dashed") +
  geom_violin(width = 0.78, fill = "#F3F7F7", color = NA, trim = TRUE, alpha = 0.68) +
  geom_boxplot(width = 0.28, outlier.shape = NA, fill = "white", color = teal,
               linewidth = 0.55) +
  geom_jitter(width = 0.095, height = 0, size = 1.08, alpha = 0.58, color = teal) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.0,
               fill = gold, color = ink, stroke = 0.35) +
  labs(title = "Paired MAGMA-module response",
       subtitle = "Each point = one paired patient",
       x = NULL, y = "Paired module-score delta") +
  theme_pub(9.3) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1),
        plot.margin = margin(6, 10, 6, 7))

# C. Patient-aware unadjusted model estimates.
pC <- ggplot(panelC, aes(y = module_short, x = effect_estimate)) +
  geom_vline(xintercept = 0, color = "#A4AAAC", linetype = "dashed", linewidth = 0.40) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.17,
                color = teal_mid, linewidth = 0.98, orientation = "y") +
  geom_point(size = 3.85, shape = 21, fill = teal, color = "white", stroke = 0.44) +
  annotate("text", x = max(panelC$ci_high) + 0.006, y = 4.10,
           label = "All FDR = 0.036", hjust = 0,
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
  geom_text(aes(label = rho_label), family = font_family, size = 3.35, color = ink) +
  scale_fill_gradient(low = "#F9F5F4", high = "#BE9185", limits = c(0, 0.85),
                      breaks = c(0, 0.4, 0.8), name = "Spearman rho") +
  labs(title = "Injury/remodeling coupling",
       subtitle = "Bulk proxy correlations; not mechanism", x = NULL, y = NULL) +
  theme_pub(8.8) +
  theme(axis.text.x = element_text(angle = 16, hjust = 1),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center", legend.key.width = unit(1.25, "cm"),
        legend.key.height = unit(0.26, "cm"),
        legend.margin = margin(t = -5), legend.box.margin = margin(t = -7),
        plot.margin = margin(6, 7, 2, 7))

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
           ymin = module_y - 0.175, ymax = module_y + 0.175,
           xmid = (stage_x - 0.39 + 0.78 * cum0 + stage_x - 0.39 + 0.78 * (cum0 + fraction)) / 2)]
pE <- ggplot() +
  geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = outcome),
            color = "white", linewidth = 0.42) +
  geom_text(data = seg, aes(x = xmid, y = module_y, label = n_models),
            family = font_family, fontface = "bold", size = 3.15, color = ink) +
  annotate("text", x = 1.5, y = 4.66,
           label = "Paired bulk shift was positive in all four modules before adjustment",
           family = font_family, fontface = "bold", size = 2.85, color = teal) +
  annotate("text", x = 1, y = 4.00, label = "Single-signature\nadjustment",
           family = font_family, fontface = "bold", size = 3.15, lineheight = 0.92, color = ink) +
  annotate("text", x = 2, y = 4.00, label = "Compact\nadjustment",
           family = font_family, fontface = "bold", size = 3.15, lineheight = 0.92, color = ink) +
  scale_fill_manual(values = c("Direction retained" = teal_mid, "Attenuated" = gold,
                               "Lost / explained" = brick_mid), drop = FALSE) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.46)) +
  scale_y_continuous(breaks = unname(row_y), labels = names(row_y),
                     limits = c(0.58, 4.84), expand = c(0, 0)) +
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
  geom_tile(width = 0.94, height = 0.66, color = "white", linewidth = 0.65) +
  geom_text(aes(label = label), family = font_family, fontface = "bold",
            size = 3.05, lineheight = 0.93, color = ink) +
  annotate("text", x = 2, y = 1.38, label = "SUPPORTED INTERPRETATION",
           family = font_family, fontface = "bold", size = 3.05, color = teal) +
  annotate("text", x = 5, y = 1.38, label = "NOT CLAIMED",
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

# Export independent, editable components only. No assembled Figure 4 is generated.
panel_specs <- list(
  A = list(plot = pA, stem = "panelA_paired_bulk_resource", width = 6.4, height = 3.0),
  B = list(plot = pB, stem = "panelB_paired_module_response", width = 6.8, height = 3.4),
  C = list(plot = pC, stem = "panelC_patient_aware_effect_summary", width = 6.5, height = 3.5),
  D = list(plot = pD, stem = "panelD_injury_remodeling_coupling_heatmap", width = 6.8, height = 3.8),
  E = list(plot = pE, stem = "panelE_adjustment_attenuation_summary", width = 12.5, height = 3.8),
  F = list(plot = pF, stem = "panelF_interpretation_boundary", width = 12.5, height = 2.0)
)

save_pdf <- function(path, plot, width, height) {
  grDevices::pdf(path, width = width, height = height, family = font_family,
                 useDingbats = FALSE, bg = "white")
  print(plot)
  grDevices::dev.off()
}

tag_panel_handoff <- function(p, tag) {
  tag_panel(p, tag) +
    theme(plot.title = element_text(margin = margin(l = 22, b = 2)))
}

for (panel_id in names(panel_specs)) {
  spec <- panel_specs[[panel_id]]
  labeled <- tag_panel_handoff(spec$plot, panel_id)
  save_pdf(file.path(fig_dir, paste0(spec$stem, ".pdf")), labeled, spec$width, spec$height)
  ggsave(file.path(fig_dir, paste0(spec$stem, ".svg")), labeled,
         width = spec$width, height = spec$height, units = "in",
         device = svglite::svglite, bg = "white")
  ggsave(file.path(fig_dir, paste0("panel", panel_id, "_clean_no_letter.svg")), spec$plot,
         width = spec$width, height = spec$height, units = "in",
         device = svglite::svglite, bg = "white")
  cat("Exported panel", panel_id, "as labeled PDF/SVG and clean SVG\n")
}

source_A <- file.path(in_dir, "figure4_panelA_source.tsv")
source_B <- file.path(in_dir, "figure4_panelB_source.tsv")
source_C <- file.path(in_dir, "figure4_panelC_source.tsv")
source_D <- file.path(in_dir, "figure4_panelF_source.tsv")
source_E <- file.path(in_dir, "figure4_panelE_source.tsv")
source_F <- file.path(in_dir, "figure4_panelG_source.tsv")
claim_lock <- file.path(in_dir, "gse73680_final_claim_lock_stage5C1.tsv")

lock_rows <- list(
  data.table(panel = "A", display_element = c("Dataset", "Sample count", "Patient count", "Paired-patient count", "Expression scale", "Reference condition", "Plaque/stone condition"),
             locked_value = c("GSE73680", "55 samples", "29 patients", "26 paired patients", "log2(x + 1) module scores", "Control / adjacent papilla", "Plaque / stone papilla"),
             source_file = source_A, claim_role = "paired bulk resource", must_not_change = "yes",
             notes = "Locked cohort/design text from v0.7"),
  data.table(panel = "B", display_element = c("Module label", "Module label", "Module label", "Module label", "Paired-delta observations"),
             locked_value = c("R1", "R1-R3", "Top 50", "Top 100", "104 observations: 26 paired patients x 4 modules"),
             source_file = source_B, claim_role = "paired bulk module shift", must_not_change = "yes",
             notes = c(rep("Short display label; full definition remains in legend", 4), "All source-backed paired deltas remain locked"))
)

lock_rows[[length(lock_rows) + 1L]] <- panelC[, .(
  panel = "C",
  display_element = paste0(as.character(module_short), " patient-aware estimate"),
  locked_value = sprintf("effect=%.15g; 95%% CI=[%.15g, %.15g]; FDR=%.15g", effect_estimate, ci_low, ci_high, fdr),
  source_file = source_C,
  claim_role = "patient-aware paired effect summary",
  must_not_change = "yes",
  notes = "Displayed common rounded label: All FDR = 0.036"
)]

lock_rows[[length(lock_rows) + 1L]] <- panelF[, .(
  panel = "D",
  display_element = paste(as.character(module_short), "x", as.character(signature_short)),
  locked_value = sprintf("Spearman rho=%.15g", spearman_r),
  source_file = source_D,
  claim_role = "bulk proxy correlation; not mechanism",
  must_not_change = "yes",
  notes = "FDR remains in Source Data and must not be added to the heatmap body"
)]

lock_rows[[length(lock_rows) + 1L]] <- rbindlist(lapply(unname(short_labels), function(m) {
  data.table(
    panel = "E",
    display_element = c(paste(m, "single-signature adjustment"), paste(m, "compact adjustment")),
    locked_value = c("3 retained / 3 attenuated / 2 lost", "3 retained / 2 lost"),
    source_file = source_E,
    claim_role = "adjustment attenuation summary",
    must_not_change = "yes",
    notes = "State colors are locked: teal retained; gold attenuated; terracotta lost/explained"
  )
}))

lock_rows[[length(lock_rows) + 1L]] <- boundary[, .(
  panel = "F",
  display_element = paste0(group, " box"),
  locked_value = gsub("- ", "-", gsub("\\n", " ", label)),
  source_file = paste(source_F, claim_lock, sep = ";"),
  claim_role = ifelse(group == "SUPPORTED INTERPRETATION", "supported claim", "claim boundary"),
  must_not_change = "yes",
  notes = "Locked v0.7 boundary wording"
)]

lock_table <- rbindlist(lock_rows, use.names = TRUE, fill = TRUE)
fwrite(lock_table, file.path(table_dir, "figure4_source_value_lock.tsv"), sep = "\t")

cat("Wrote source-value lock with", nrow(lock_table), "rows\n")
cat("No full Figure 4 was generated. data_changed=no; claim_changed=no\n")
