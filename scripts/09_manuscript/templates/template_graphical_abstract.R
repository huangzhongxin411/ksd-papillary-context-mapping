# Template: graphical abstract / Figure 1 framework.
# Expects `figure_theme_highimpact_v0.2.R` to be sourced.

make_graphical_abstract_template <- function(cards, claim_box, title, stem) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.94, label = title, fontface = "bold", size = 4.5, color = project_palette[["ink"]]) +
    geom_label(data = cards, aes(x = x, y = 0.68, label = layer, fill = fill), color = "white",
      fontface = "bold", size = 3.0, label.padding = unit(0.24, "lines"), linewidth = 0) +
    geom_label(data = cards, aes(x = x, y = 0.52, label = badge), fill = "white", color = project_palette[["ink"]],
      size = 2.6, label.padding = unit(0.24, "lines"), linewidth = 0.35, lineheight = 0.92) +
    annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.14, ymax = 0.35, fill = "#F8FAFA",
      color = project_palette[["medium_grey"]], linewidth = 0.55) +
    annotate("text", x = 0.12, y = 0.30, label = claim_box, hjust = 0, vjust = 1,
      size = 2.7, color = project_palette[["ink"]], lineheight = 0.95) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.03, 0.97), ylim = c(0.08, 0.98), clip = "off") +
    theme_void(base_size = 10)
  save_figure_pdf_png(p, stem, width = 13.2, height = 5.9)
  p
}
