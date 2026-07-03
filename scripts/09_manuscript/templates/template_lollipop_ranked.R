# Template: ranked lollipop contributor plot.
# Required columns: label, score, group. Optional: size_value.

make_lollipop_ranked_template <- function(dt, title, xlab = "Score") {
  if (!"size_value" %in% names(dt)) dt$size_value <- 1
  dt$label <- factor(dt$label, levels = rev(dt$label))
  ggplot(dt, aes(score, label)) +
    geom_segment(aes(x = 0, xend = score, yend = label), color = default_project_palette[["medium_grey"]], linewidth = 0.7) +
    geom_point(aes(fill = group, size = size_value), shape = 21, color = "#5A6062", stroke = 0.35) +
    scale_fill_project() +
    labs(title = title, x = xlab, y = NULL, fill = NULL, size = NULL) +
    theme_highimpact(9)
}
