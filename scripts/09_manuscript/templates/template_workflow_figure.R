# Template: workflow schematic with semantic steps.
# Required columns: x, y, label, fill.

make_workflow_figure_template <- function(dt, title = "Workflow") {
  seg <- NULL
  if (nrow(dt) > 1) {
    seg <- data.frame(
      x = dt$x[-nrow(dt)] + 0.08,
      xend = dt$x[-1] - 0.08,
      y = dt$y[-nrow(dt)],
      yend = dt$y[-1]
    )
  }
  ggplot(dt, aes(x, y)) +
    geom_label(aes(label = label, fill = fill), color = "white", fontface = "bold", size = 3.0,
      label.padding = unit(0.25, "lines"), linewidth = 0) +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
      arrow = arrow(length = unit(0.10, "in")), color = default_project_palette[["dark_grey"]], linewidth = 0.55) +
    scale_fill_identity() +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 16, color = default_project_palette[["ink"]]))
}
