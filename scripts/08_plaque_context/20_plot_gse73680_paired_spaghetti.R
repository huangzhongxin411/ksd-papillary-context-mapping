suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  control = "#8AA0A8",
  plaque = "#B08A45",
  magma = "#3E6672",
  p1 = "#B59A5B",
  muted = "#C9C9C9",
  text = "#303030"
)

label_module <- function(x) {
  map <- c(
    MAGMA_top50 = "MAGMA top 50",
    MAGMA_top100 = "MAGMA top 100",
    MAGMA_FDR = "MAGMA FDR",
    MAGMA_suggestive = "MAGMA suggestive",
    P1_core_TAL_candidates = "P1 core"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

scores <- fread("results/gse73680/tables/gse73680_patient_level_module_score_matrix.tsv")
cons <- fread("results/gse73680/tables/gse73680_paired_direction_consistency.tsv")

modules_keep <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")
scores <- scores[module_name %in% modules_keep]

wide <- dcast(scores, module_name + patient_id ~ group_curated, value.var = "patient_level_module_score")
paired <- wide[!is.na(control_or_adjacent) & !is.na(plaque_or_stone_papilla)]
paired[, paired_delta := plaque_or_stone_papilla - control_or_adjacent]
paired[, direction := fifelse(paired_delta > 0, "positive", fifelse(paired_delta < 0, "negative", "zero"))]
paired[, module_label := factor(label_module(module_name),
                                levels = label_module(modules_keep))]

long <- melt(
  paired,
  id.vars = c("module_name", "module_label", "patient_id", "paired_delta", "direction"),
  measure.vars = c("control_or_adjacent", "plaque_or_stone_papilla"),
  variable.name = "group_curated",
  value.name = "patient_level_module_score"
)
long[, group_label := factor(
  fifelse(group_curated == "control_or_adjacent", "Control/adjacent", "Plaque/stone papilla"),
  levels = c("Control/adjacent", "Plaque/stone papilla")
)]

summary <- paired[, .(
  n_paired_patients = .N,
  n_positive_delta = sum(paired_delta > 0, na.rm = TRUE),
  positive_fraction = mean(paired_delta > 0, na.rm = TRUE),
  median_delta = median(paired_delta, na.rm = TRUE)
), by = .(module_name, module_label)]
summary <- merge(summary, cons[feature_type == "module", .(module_name = module_or_gene, sign_test_p, interpretation)],
                 by = "module_name", all.x = TRUE)

fwrite(long, "results/tables/gse73680_paired_module_delta_long.tsv", sep = "\t")
fwrite(summary, "results/tables/gse73680_paired_module_delta_summary.tsv", sep = "\t")

fig <- ggplot(long, aes(group_label, patient_level_module_score, group = patient_id)) +
  geom_line(aes(color = direction), alpha = 0.52, linewidth = 0.38) +
  geom_point(aes(fill = group_label), shape = 21, color = "#555555", stroke = 0.2, size = 1.8, alpha = 0.9) +
  facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c(positive = pal$magma, negative = "#9A5F52", zero = pal$muted)) +
  scale_fill_manual(values = c("Control/adjacent" = pal$control, "Plaque/stone papilla" = pal$plaque)) +
  labs(
    title = "GSE73680 paired patient-level MAGMA module shifts",
    subtitle = "Each line represents one paired patient; positive direction supports module-level disease-context association",
    x = NULL,
    y = "Patient-level module score",
    color = "Delta",
    fill = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "#555555"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/figure4_paired_module_spaghetti_v0.1.pdf", fig,
       width = 10.8, height = 6.7, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure4_paired_module_spaghetti_v0.1.png", fig,
       width = 10.8, height = 6.7, units = "in", dpi = 260, bg = "white")

dot <- summary[, module_label := factor(as.character(module_label), levels = rev(label_module(modules_keep)))]
dot[, support_class := fifelse(module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive") &
                                 positive_fraction >= 0.65,
                               "MAGMA directional support",
                               "Exploratory/background")]
fig_dot <- ggplot(dot, aes(median_delta, module_label, fill = support_class)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.3) +
  geom_point(shape = 21, color = "#555555", stroke = 0.25, size = 4.2) +
  geom_text(aes(label = sprintf("%d/%d positive", n_positive_delta, n_paired_patients)),
            nudge_x = 0.045, size = 2.8, hjust = 0, color = pal$text) +
  scale_fill_manual(values = c("MAGMA directional support" = pal$magma,
                               "Exploratory/background" = pal$muted)) +
  coord_cartesian(xlim = c(min(dot$median_delta, na.rm = TRUE) - 0.05,
                           max(dot$median_delta, na.rm = TRUE) + 0.28)) +
  labs(
    title = "Paired module delta summary",
    x = "Median paired delta",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("results/figures/gse73680_paired_module_delta_dotplot.pdf", fig_dot,
       width = 7.6, height = 4.6, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/gse73680_paired_module_delta_dotplot.png", fig_dot,
       width = 7.6, height = 4.6, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# GSE73680 Paired Module Spaghetti v0.1",
  "",
  "The paired spaghetti plot visualizes patient-level module scores for control/adjacent and plaque/stone papilla samples.",
  "",
  "Interpretation boundary: this supports paired module-level disease-context association for MAGMA-prioritized modules when most paired patients shift in the same direction. It does not establish genetic causality or cell-type-specific disease response."
), "docs/gse73680_paired_module_spaghetti_v0.1.md")

message("wrote GSE73680 paired module spaghetti outputs")
