suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

r_packages <- data.table(
  name = c("ggplot2", "scales", "data.table", "cowplot", "patchwork", "ggrepel", "ggrastr", "ragg",
    "svglite", "colorspace", "ComplexHeatmap", "circlize", "grid", "ggforce", "ggnewscale", "ggtext"),
  importance = c("required", "required", "required", "recommended", "recommended", "recommended", "recommended",
    "recommended", "recommended", "recommended", "recommended", "recommended", "required", "optional", "optional", "optional")
)

system_tools <- data.table(
  name = c("inkscape", "magick", "gs", "pdftoppm", "qpdf"),
  importance = c("optional", "recommended", "recommended", "recommended", "recommended")
)

pkg_dt <- r_packages[, .(
  dependency_type = "R package",
  name,
  importance,
  available = requireNamespace(name, quietly = TRUE),
  version = if (requireNamespace(name, quietly = TRUE)) as.character(utils::packageVersion(name)) else NA_character_,
  path_or_note = NA_character_
), by = seq_len(nrow(r_packages))][, seq_len := NULL]

tool_dt <- system_tools[, .(
  dependency_type = "system tool",
  name,
  importance,
  available = nzchar(Sys.which(name)),
  version = NA_character_,
  path_or_note = unname(Sys.which(name))
), by = seq_len(nrow(system_tools))][, seq_len := NULL]
tool_dt[available == FALSE, path_or_note := NA_character_]

out <- rbind(pkg_dt, tool_dt, fill = TRUE)
out[, status := fifelse(available, "available",
  fifelse(importance == "required", "missing_required",
    fifelse(importance == "recommended", "missing_recommended", "missing_optional")))]
out[, fallback := fifelse(available, "none",
  fifelse(name == "ggtext", "use plain text labels instead of markdown labels",
    fifelse(dependency_type == "system tool", "skip system-level PDF/SVG conversion or compression QC and record warning",
      "use base plotting/device fallback")))]

fwrite(out, "results/tables/figure_dependency_check_v0.2.tsv", sep = "\t")
message("Wrote results/tables/figure_dependency_check_v0.2.tsv")
