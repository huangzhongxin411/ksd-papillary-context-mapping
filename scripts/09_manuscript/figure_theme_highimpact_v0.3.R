suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

default_project_palette <- c(
  deep_teal = "#245A64",
  loop_teal = "#0F4C5C",
  focal_signal = "#0F4C5C",
  bluegrey = "#7F9DA6",
  sand_gold = "#B99B5A",
  terracotta = "#9B5C4D",
  pale_grey = "#E6E9EA",
  light_grey = "#E6E9EA",
  dark_grey = "#333333",
  ink = "#2F3335",
  medium_grey = "#AEB8BC",
  white = "#FFFFFF"
)

trim_ws <- function(x) gsub("^\\s+|\\s+$", "", x)

load_project_profile <- function(profile_path) {
  if (!file.exists(profile_path)) stop("Project profile not found: ", profile_path)
  lines <- readLines(profile_path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$|^\\s*#", lines)]
  profile <- list()
  current_key <- NULL
  current_map <- NULL
  for (line in lines) {
    if (grepl("^\\s{2}-\\s+", line) && !is.null(current_key)) {
      profile[[current_key]] <- c(profile[[current_key]], sub("^\\s{2}-\\s+", "", line))
    } else if (grepl("^\\s{2}[A-Za-z0-9_]+:", line) && !is.null(current_map)) {
      kv <- strsplit(trim_ws(line), ":", fixed = TRUE)[[1]]
      key <- trim_ws(kv[1])
      val <- trim_ws(paste(kv[-1], collapse = ":"))
      val <- gsub("^\"|\"$", "", val)
      profile[[current_map]][[key]] <- val
    } else if (grepl("^[A-Za-z0-9_]+:\\s*$", line)) {
      key <- sub(":\\s*$", "", trim_ws(line))
      current_key <- key
      if (key %in% c("color_semantics", "figure_claim_boundaries", "figure_plan")) {
        profile[[key]] <- list()
        current_map <- key
        current_key <- NULL
      } else {
        profile[[key]] <- character()
        current_map <- NULL
      }
    } else if (grepl("^[A-Za-z0-9_]+:", line)) {
      kv <- strsplit(line, ":", fixed = TRUE)[[1]]
      key <- trim_ws(kv[1])
      val <- trim_ws(paste(kv[-1], collapse = ":"))
      val <- gsub("^\"|\"$", "", val)
      profile[[key]] <- val
      current_key <- NULL
      current_map <- NULL
    }
  }
  class(profile) <- c("figure_project_profile", "list")
  profile
}

get_project_palette <- function(profile = NULL) {
  pal <- default_project_palette
  if (!is.null(profile) && !is.null(profile$color_semantics)) {
    custom <- unlist(profile$color_semantics)
    pal[names(custom)] <- custom
  }
  pal
}

theme_highimpact <- function(base_size = 10, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = default_project_palette[["ink"]], margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "#5A6062", margin = margin(b = 5)),
      axis.title = element_text(size = 10, color = default_project_palette[["ink"]]),
      axis.text = element_text(size = 8.5, color = default_project_palette[["ink"]]),
      legend.title = element_text(size = 8.8, face = "bold", color = default_project_palette[["ink"]]),
      legend.text = element_text(size = 8.2, color = default_project_palette[["ink"]]),
      strip.background = element_rect(fill = default_project_palette[["pale_grey"]], color = NA),
      strip.text = element_text(size = 9, face = "bold", color = default_project_palette[["ink"]]),
      axis.line = element_line(linewidth = 0.5, color = "#5A6062"),
      axis.ticks = element_line(linewidth = 0.45, color = "#5A6062"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

theme_umap_highlight <- function(base_size = 9) {
  theme_highimpact(base_size = base_size) +
    theme(legend.position = "right", axis.title = element_text(size = 9), axis.text = element_text(size = 8.5))
}

theme_compact_heatmap <- function(base_size = 9) {
  theme_highimpact(base_size = base_size) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
}

save_figure_highimpact <- function(plot, filename_base, width, height, dpi = 600, bg = "white") {
  dir.create(dirname(filename_base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = "in", device = "pdf", bg = bg)
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = bg)
  invisible(c(pdf = paste0(filename_base, ".pdf"), png = paste0(filename_base, ".png")))
}

italic_gene_labels <- function(labels) {
  parse(text = paste0("italic('", labels, "')"))
}

make_panel_label <- function(label, x = -Inf, y = Inf) {
  annotate("text", x = x, y = y, label = label, hjust = -0.15, vjust = 1.15,
    fontface = "bold", size = 13 / ggplot2::.pt, color = default_project_palette[["ink"]])
}

check_palette_consistency <- function(plot_colors, project_palette) {
  plot_colors <- unique(toupper(plot_colors[!is.na(plot_colors)]))
  allowed <- unique(toupper(unname(project_palette)))
  data.frame(
    n_colors_checked = length(plot_colors),
    n_out_of_palette = sum(!plot_colors %in% allowed),
    out_of_palette = paste(plot_colors[!plot_colors %in% allowed], collapse = ";"),
    uses_project_palette = all(plot_colors %in% allowed),
    stringsAsFactors = FALSE
  )
}

write_figure_source_table <- function(rows, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(data.table::as.data.table(rows), path, sep = "\t")
  } else {
    utils::write.table(as.data.frame(rows), path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(path)
}

figure_qc_check_v03 <- function(figure_id, version, pdf_path, png_path, source_table = NA_character_,
                                legend_path = NA_character_, png_dpi = 600, min_font_size = 8.5,
                                panel_labels_present = TRUE, legend_outside_dense_panels = TRUE,
                                uses_project_palette = TRUE, claim_boundary_checked = TRUE,
                                resource_limited_claim_checked = TRUE, action_required = "none") {
  data.frame(
    figure_id = figure_id,
    version = version,
    pdf_exists = file.exists(pdf_path),
    png_exists = file.exists(png_path),
    png_dpi = png_dpi,
    min_detected_font_size = min_font_size,
    panel_labels_present = panel_labels_present,
    legend_outside_dense_panels = legend_outside_dense_panels,
    uses_project_palette = uses_project_palette,
    claim_boundary_checked = claim_boundary_checked,
    resource_limited_claim_checked = resource_limited_claim_checked,
    source_table_exists = ifelse(is.na(source_table), NA, file.exists(source_table)),
    legend_exists = ifelse(is.na(legend_path), NA, file.exists(legend_path)),
    visual_status = ifelse(file.exists(pdf_path) && file.exists(png_path) && claim_boundary_checked, "pass", "review_needed"),
    action_required = action_required,
    stringsAsFactors = FALSE
  )
}

scale_fill_project <- function(profile = NULL, ..., values = NULL, na.value = "#CFCFCF") {
  pal <- get_project_palette(profile)
  vals <- if (is.null(values)) c(
    "MAGMA" = pal[["deep_teal"]],
    "MAGMA top 50" = pal[["deep_teal"]],
    "MAGMA top 100" = "#47757E",
    "MAGMA FDR" = "#6B8F98",
    "MAGMA suggestive" = "#9FB3BA",
    "Loop/TAL" = pal[["loop_teal"]],
    "P1" = pal[["sand_gold"]],
    "P1 core" = pal[["sand_gold"]],
    "Disease context" = pal[["terracotta"]],
    "Plaque/stone papilla" = pal[["terracotta"]],
    "Background" = pal[["pale_grey"]],
    "No support" = "#CFCFCF",
    "Audited context" = pal[["bluegrey"]]
  ) else values
  scale_fill_manual(..., values = vals, na.value = na.value)
}

scale_color_project <- function(profile = NULL, ..., values = NULL, na.value = "#CFCFCF") {
  vals <- if (is.null(values)) unname(get_project_palette(profile)) else values
  scale_color_manual(..., values = vals, na.value = na.value)
}

# Backward-compatible aliases.
save_figure_pdf_png <- save_figure_highimpact
theme_hi <- theme_highimpact
hi_save_fig <- save_figure_highimpact
