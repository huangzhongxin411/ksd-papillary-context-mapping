# Template: audited cell atlas UMAP highlight.
# Required columns: UMAP_1, UMAP_2, cell_label, is_target.

make_umap_highlight_template <- function(dt, title = "Audited snRNA-seq atlas") {
  ggplot() +
    geom_point(data = dt, aes(UMAP_1, UMAP_2), color = project_palette[["light_grey"]], alpha = 0.20, size = 0.13) +
    geom_point(data = dt[!dt$is_target, ], aes(UMAP_1, UMAP_2, color = cell_label), alpha = 0.40, size = 0.15) +
    geom_point(data = dt[dt$is_target, ], aes(UMAP_1, UMAP_2), color = project_palette[["loop_teal"]], alpha = 0.95, size = 0.46) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    scale_color_project(name = "Audited cell type") +
    theme_highimpact(9)
}
