suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

r_packages <- c(
  "ggplot2", "patchwork", "cowplot", "ggtext", "ggrepel", "ggrastr", "ragg",
  "svglite", "scales", "colorspace", "ComplexHeatmap", "circlize", "grid",
  "ggforce", "ggnewscale"
)

system_tools <- c("inkscape", "magick", "gs", "pdftoppm", "qpdf")

pkg_dt <- data.table(
  dependency_type = "R package",
  name = r_packages,
  available = vapply(r_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
  version = vapply(r_packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, FUN.VALUE = character(1)),
  path_or_note = NA_character_
)

tool_dt <- data.table(
  dependency_type = "system tool",
  name = system_tools,
  available = nzchar(Sys.which(system_tools)),
  version = NA_character_,
  path_or_note = unname(Sys.which(system_tools))
)
tool_dt[available == FALSE, path_or_note := NA_character_]

out <- rbind(pkg_dt, tool_dt, fill = TRUE)
out[, status := ifelse(available, "available", "missing")]

fwrite(out, "results/tables/figure_dependency_check_v0.1.tsv", sep = "\t")
message("Wrote results/tables/figure_dependency_check_v0.1.tsv")
