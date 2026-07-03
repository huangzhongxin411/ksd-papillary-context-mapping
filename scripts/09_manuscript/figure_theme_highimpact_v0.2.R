suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

project_palette <- c(
  deep_teal = "#245A64",
  loop_teal = "#0F4C5C",
  sand_gold = "#B99B5A",
  terracotta = "#9B5C4D",
  light_grey = "#E6E9EA",
  bluegrey = "#7F9DA6",
  ink = "#2F3335",
  medium_grey = "#AEB8BC",
  pale_bluegrey = "#E9EFF1",
  pale_sand = "#EFE5CF",
  white = "#FFFFFF"
)

project_color_values <- c(
  "MAGMA" = project_palette[["deep_teal"]],
  "MAGMA top 50" = project_palette[["deep_teal"]],
  "MAGMA top 100" = "#47757E",
  "MAGMA FDR" = "#6B8F98",
  "MAGMA suggestive" = "#9FB3BA",
  "Loop/TAL" = project_palette[["loop_teal"]],
  "P1" = project_palette[["sand_gold"]],
  "P1 core" = project_palette[["sand_gold"]],
  "P1 candidate" = project_palette[["sand_gold"]],
  "Disease context" = project_palette[["terracotta"]],
  "Plaque/stone papilla" = project_palette[["terracotta"]],
  "Background" = project_palette[["light_grey"]],
  "No support" = "#CFCFCF",
  "Audited context" = project_palette[["bluegrey"]]
)

theme_highimpact <- function(base_size = 9, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, color = project_palette[["ink"]], margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9.2, color = "#5A6062", margin = margin(b = 4)),
      axis.title = element_text(size = 9.2, color = project_palette[["ink"]]),
      axis.text = element_text(size = 8.5, color = project_palette[["ink"]]),
      legend.title = element_text(size = 8.8, face = "bold", color = project_palette[["ink"]]),
      legend.text = element_text(size = 8.5, color = project_palette[["ink"]]),
      strip.background = element_rect(fill = project_palette[["light_grey"]], color = NA),
      strip.text = element_text(size = 8.8, face = "bold", color = project_palette[["ink"]]),
      axis.line = element_line(linewidth = 0.5, color = "#5A6062"),
      axis.ticks = element_line(linewidth = 0.45, color = "#5A6062"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

scale_fill_project <- function(..., values = project_color_values, na.value = "#CFCFCF") {
  scale_fill_manual(..., values = values, na.value = na.value)
}

scale_color_project <- function(..., values = project_color_values, na.value = "#CFCFCF") {
  scale_color_manual(..., values = values, na.value = na.value)
}

save_figure_pdf_png <- function(plot, stem, width, height, dpi = 600, bg = "white") {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = "pdf", bg = bg)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = bg)
  invisible(c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png")))
}

italicize_gene_labels <- function(labels) {
  parse(text = paste0("italic('", labels, "')"))
}

panel_label_style <- function(label, x = -Inf, y = Inf) {
  annotate("text", x = x, y = y, label = label, hjust = -0.2, vjust = 1.2,
    fontface = "bold", size = 13 / ggplot2::.pt, color = project_palette[["ink"]])
}

figure_qc_check <- function(figure_id, panel_id = "all", pdf_path, png_path,
                            panel_label_present = TRUE, legend_position_ok = TRUE,
                            uses_project_palette = TRUE, claim_boundary_ok = TRUE,
                            min_font_size = 8.5, notes = "") {
  data.frame(
    figure_id = figure_id,
    panel_id = panel_id,
    min_font_size = min_font_size,
    panel_label_present = panel_label_present,
    legend_position_ok = legend_position_ok,
    uses_project_palette = uses_project_palette,
    claim_boundary_ok = claim_boundary_ok,
    png_600dpi_exists = file.exists(png_path),
    pdf_exists = file.exists(pdf_path),
    visual_review_status = ifelse(file.exists(png_path) && file.exists(pdf_path) && claim_boundary_ok, "pass", "review_needed"),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

# Backward-compatible aliases used by Phase 15 scripts.
hi_pal <- as.list(project_palette)
hi_pal$loop_tal <- project_palette[["loop_teal"]]
hi_pal$deep_teal <- project_palette[["deep_teal"]]
hi_pal$sand <- project_palette[["sand_gold"]]
hi_pal$bg_grey <- project_palette[["light_grey"]]
hi_pal$dark_grey <- "#5A6062"
theme_hi <- theme_highimpact
hi_save_fig <- save_figure_pdf_png
hi_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, path, useBytes = TRUE)
}
