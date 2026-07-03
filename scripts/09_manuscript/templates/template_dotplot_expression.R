# Template: expression dotplot.
# Required columns: cell_type, gene, pct_expressed, avg_expression.

make_dotplot_expression_template <- function(dt, title = "Expression across audited cell types") {
  ggplot(dt, aes(cell_type, gene)) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21,
      color = project_palette[["medium_grey"]], stroke = 0.35) +
    scale_fill_gradientn(colors = c("#EEF3F4", project_palette[["bluegrey"]], project_palette[["deep_teal"]]),
      name = "Average\nexpression") +
    scale_size_continuous(range = c(0.9, 5.0), breaks = c(0.2, 0.5, 0.8), labels = c("20", "50", "80"),
      name = "Detected (%)") +
    labs(title = title, x = NULL, y = NULL) +
    theme_highimpact(8.8) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
}
