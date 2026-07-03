# Template: paired patient spaghetti/slopegraph.
# Required columns: group_label, score, patient_id, module_label.

make_spaghetti_paired_template <- function(dt, title = "Paired disease-context shifts") {
  ggplot(dt, aes(group_label, score, group = patient_id)) +
    geom_line(color = project_palette[["medium_grey"]], alpha = 0.22, linewidth = 0.28) +
    geom_point(aes(fill = group_label), shape = 21, color = "white", stroke = 0.10, size = 1.20, alpha = 0.82) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.25, color = project_palette[["ink"]]) +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2.25, color = project_palette[["ink"]]) +
    facet_wrap(~ module_label, ncol = 3, scales = "free_y") +
    scale_fill_project(values = c("Control/adjacent" = project_palette[["bluegrey"]],
      "Plaque/stone papilla" = project_palette[["sand_gold"]])) +
    labs(title = title, x = NULL, y = "Module score", fill = NULL) +
    theme_highimpact(8.4) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
}
