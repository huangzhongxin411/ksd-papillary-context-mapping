# Template: evidence matrix with glyphs/text.
# Required columns: row_label, evidence, call.

make_evidence_matrix_template <- function(dt, title = "Evidence matrix") {
  dt$evidence <- factor(dt$evidence, levels = unique(dt$evidence))
  dt$row_label <- factor(dt$row_label, levels = rev(unique(dt$row_label)))
  glyph <- dt[dt$call %in% c("strong", "moderate", "contextual"), ]
  text_dt <- dt[!dt$call %in% c("strong", "moderate", "contextual"), ]
  ggplot() +
    geom_tile(data = dt, aes(evidence, row_label), fill = "#F7F8F8", color = "white", linewidth = 0.6) +
    geom_point(data = glyph[glyph$call == "strong", ], aes(evidence, row_label), shape = 21,
      fill = default_project_palette[["deep_teal"]], color = default_project_palette[["deep_teal"]], size = 4.4) +
    geom_point(data = glyph[glyph$call == "moderate", ], aes(evidence, row_label), shape = 21,
      fill = default_project_palette[["sand_gold"]], color = default_project_palette[["dark_grey"]], size = 4.4, stroke = 0.35) +
    geom_point(data = glyph[glyph$call == "contextual", ], aes(evidence, row_label), shape = 21,
      fill = "white", color = default_project_palette[["dark_grey"]], size = 4.4, stroke = 0.7) +
    geom_text(data = text_dt, aes(evidence, row_label, label = call), size = 2.6, color = default_project_palette[["ink"]]) +
    labs(title = title, x = NULL, y = NULL) +
    theme_compact_heatmap(8.8) +
    theme(axis.text.y = element_text(face = "italic"))
}
