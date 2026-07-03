suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
fig_dir <- "results/figures/supplementary/spatial_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

p1_genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
pal <- c(light = "#EAF1F2", sand = "#C99A2E", slate = "#6F929B", ink = "#22313B", terracotta = "#A85B4B")
objects <- list.files(obj_dir, pattern = "_spatial_scores_v0\\.1\\.rds$", full.names = TRUE)

get_data_matrix <- function(obj) {
  tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
           error = function(e) GetAssayData(obj, assay = "Spatial", slot = "data"))
}
get_count_matrix <- function(obj) {
  tryCatch(GetAssayData(obj, assay = "Spatial", layer = "counts"),
           error = function(e) GetAssayData(obj, assay = "Spatial", slot = "counts"))
}

rows <- list()
plot_gene_sample <- function(path, gene) {
  sample_id <- sub("_spatial_scores_v0\\.1\\.rds$", "", basename(path))
  obj <- readRDS(path)
  if (!gene %in% rownames(obj)) {
    return(ggplot() + annotate("text", 0.5, 0.5, label = paste(sample_id, "\nnot detected"), size = 3) +
             xlim(0, 1) + ylim(0, 1) + theme_void())
  }
  SpatialFeaturePlot(obj, features = gene, images = names(obj@images)[1], pt.size.factor = 1.6,
                     alpha = c(0.35, 1), combine = FALSE)[[1]] +
    scale_fill_gradient(low = pal["light"], high = pal["sand"]) +
    labs(title = sample_id, fill = "log expr") +
    theme(plot.title = element_text(family = "Arial", face = "bold", size = 10, color = pal["ink"]),
          legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}

for (path in objects) {
  sample_id <- sub("_spatial_scores_v0\\.1\\.rds$", "", basename(path))
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Spatial"
  if (!"data" %in% Layers(obj[["Spatial"]])) {
    obj <- NormalizeData(obj, assay = "Spatial", normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
  }
  counts <- get_count_matrix(obj)
  data <- get_data_matrix(obj)
  for (gene in p1_genes) {
    present <- gene %in% rownames(obj)
    if (present) {
      detected_fraction <- mean(counts[gene, ] > 0)
      mean_expr <- mean(data[gene, ])
      max_expr <- max(data[gene, ])
    } else {
      detected_fraction <- 0
      mean_expr <- NA_real_
      max_expr <- NA_real_
    }
    rows[[length(rows) + 1]] <- tibble(
      sample_id = sample_id,
      gene = gene,
      gene_present = present,
      detected_spot_fraction = detected_fraction,
      mean_log_expression = mean_expr,
      max_log_expression = max_expr,
      spatial_plot_allowed = present && detected_fraction >= 0.05
    )
  }
}

audit <- bind_rows(rows) %>% arrange(gene, sample_id)
write_tsv(audit, file.path(out_dir, "p1_gene_spatial_expression_audit_v0.1.tsv"))

for (gene in p1_genes) {
  if (!any(audit$gene == gene & audit$spatial_plot_allowed)) next
  plots <- lapply(objects, plot_gene_sample, gene = gene)
  fig <- wrap_plots(plots, ncol = length(objects)) +
    plot_annotation(
      title = paste0("P1 gene spatial expression audit: ", gene),
      subtitle = "Individual gene overlay; not P1 single-gene disease validation",
      theme = theme(
        plot.title = element_text(family = "Arial", face = "bold", size = 14, color = pal["ink"]),
        plot.subtitle = element_text(family = "Arial", size = 9.5, color = pal["terracotta"])
      )
    )
  ggsave(file.path(fig_dir, paste0("p1_gene_spatial_overlay_", gene, "_v0.1.pdf")),
         fig, width = 15.5, height = 3.8, device = "pdf")
  ggsave(file.path(fig_dir, paste0("p1_gene_spatial_overlay_", gene, "_v0.1.png")),
         fig, width = 15.5, height = 3.8, dpi = 600)
}
