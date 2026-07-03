# Template: compact heatmap with numeric labels.
# Required columns: x, y, value.

make_heatmap_compact_template <- function(dt, title, legend_title = "Value") {
  dt$text_color <- ifelse(dt$value >= 0.70, "white", project_palette[["ink"]])
  ggplot(dt, aes(x, y, fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", value), color = text_color), fontface = "bold", size = 2.9) +
    scale_color_identity() +
    scale_fill_gradientn(colors = c("#F1F3F3", project_palette[["bluegrey"]], project_palette[["deep_teal"]]),
      name = legend_title) +
    labs(title = title, x = NULL, y = NULL) +
    theme_highimpact(8.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
}
