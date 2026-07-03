suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

hi_pal <- list(
  ink = "#243038",
  deep_teal = "#0F5A64",
  loop_tal = "#0F5A64",
  bluegrey = "#6F929B",
  pale_bluegrey = "#EEF2F3",
  bg_grey = "#EEF2F3",
  light_grey = "#EEF2F3",
  medium_grey = "#6F929B",
  dark_grey = "#243038",
  sand = "#C49A3A",
  pale_sand = "#EEF2F3",
  terracotta = "#A75A49",
  pale_terracotta = "#EEF2F3",
  sage = "#6F929B",
  white = "#FFFFFF"
)

hi_module_cols <- c(
  "MAGMA top 50" = hi_pal$deep_teal,
  "MAGMA top 100" = hi_pal$deep_teal,
  "MAGMA FDR" = hi_pal$bluegrey,
  "MAGMA suggestive" = hi_pal$bluegrey,
  "P1 core" = hi_pal$sand,
  "P1 candidate" = hi_pal$sand,
  "Other MAGMA gene" = hi_pal$deep_teal
)

hi_cell_cols <- c(
  "Loop/TAL" = hi_pal$loop_tal,
  "Collecting duct" = hi_pal$bluegrey,
  "Endothelial" = hi_pal$bluegrey,
  "Fibroblast/stromal" = hi_pal$bluegrey,
  "Fib/stromal" = hi_pal$bluegrey,
  "Injured epithelial" = hi_pal$terracotta,
  "Injured epi." = hi_pal$terracotta,
  "Perivascular/mural-like" = hi_pal$bluegrey,
  "Perivasc." = hi_pal$bluegrey
)

theme_hi <- function(base_size = 9, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.5, color = hi_pal$ink, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size, color = hi_pal$dark_grey, margin = margin(b = 4)),
      axis.title = element_text(size = base_size + 0.8, color = hi_pal$ink),
      axis.text = element_text(size = base_size, color = hi_pal$ink),
      legend.title = element_text(size = base_size, face = "bold", color = hi_pal$ink),
      legend.text = element_text(size = base_size - 0.2, color = hi_pal$ink),
      strip.background = element_rect(fill = hi_pal$bg_grey, color = NA),
      strip.text = element_text(size = base_size, face = "bold", color = hi_pal$ink),
      axis.line = element_line(linewidth = 0.5, color = hi_pal$dark_grey),
      axis.ticks = element_line(linewidth = 0.45, color = hi_pal$dark_grey),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = hi_pal$white, color = NA),
      legend.key = element_rect(fill = hi_pal$white, color = NA)
    )
}

theme_hi_void <- function(base_size = 10, base_family = "sans") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.2, color = hi_pal$ink),
      plot.subtitle = element_text(size = base_size, color = hi_pal$dark_grey),
      plot.background = element_rect(fill = hi_pal$white, color = NA)
    )
}

hi_save_fig <- function(plot, stem, width, height, dpi = 600) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = "pdf", bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  if (!requireNamespace("svglite", quietly = TRUE)) stop("svglite is required for editable SVG output")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
}

hi_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, path, useBytes = TRUE)
}

hi_label_cell <- function(x, short = FALSE) {
  full <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like",
    Pericyte_smooth_muscle = "Perivascular/mural-like"
  )
  short_map <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fib/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epi.",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivasc.",
    Pericyte_smooth_muscle = "Perivasc."
  )
  map <- if (short) short_map else full
  unname(ifelse(x %in% names(map), map[x], x))
}

hi_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "P<0.001", sprintf("P=%.3f", x)))
}
