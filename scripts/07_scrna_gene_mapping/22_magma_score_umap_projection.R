suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Seurat)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

audit_path <- "results/tables/gse231569_celllevel_resource_audit_v0.1.tsv"
if (!file.exists(audit_path)) {
  stop("Run scripts/06_scrna_processing/30_rebuild_gse231569_seurat.R first.")
}
audit <- fread(audit_path)
object_path <- "data/processed/gse231569_audited_seurat.rds"
resource_ok <- file.exists(object_path)

if (!resource_ok) {
  donor <- fread("results/tables/gse231569_magma_score_by_celltype_donor.tsv")
  fwrite(donor, "results/tables/gse231569_donor_celltype_magma_scores_v0.2.tsv", sep = "\t")
  blocker <- data.table(
    requested_output = c("gse231569_celllevel_magma_scores.tsv", "figure2_umap_audited_celltypes_v0.1", "figure2_magma_score_umap_v0.1"),
    status = "not_generated_resource_missing",
    reason = "No usable cell-level GSE231569 expression object and embedding table in current workspace"
  )
  fwrite(blocker, "results/tables/gse231569_celllevel_magma_score_blocker_v0.1.tsv", sep = "\t")

  p <- ggplot(audit, aes(status, factor(resource, levels = rev(resource)), fill = status)) +
    geom_col(width = 0.58, color = "#666666", linewidth = 0.25) +
    geom_text(aes(label = status), hjust = 1.05, size = 3, color = "white") +
    scale_fill_manual(values = c(available = "#3E6672", available_pdf_only = "#6F8F98",
                                 candidate_found = "#B59A5B", missing = "#9A5F52")) +
    labs(title = "GSE231569 cell-level UMAP projection readiness",
         subtitle = "Per-cell MAGMA score UMAP was not generated because cell-level resources are missing",
         x = NULL, y = NULL, fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, color = "#555555"),
          legend.position = "none", panel.grid.minor = element_blank())
  ggsave("results/figures/figure2_umap_audited_celltypes_v0.1.pdf", p, width = 8.0, height = 3.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_umap_audited_celltypes_v0.1.png", p, width = 8.0, height = 3.8, units = "in", dpi = 260, bg = "white")
  ggsave("results/figures/figure2_magma_score_umap_v0.1.pdf", p, width = 8.0, height = 3.8, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_score_umap_v0.1.png", p, width = 8.0, height = 3.8, units = "in", dpi = 260, bg = "white")
  message("cell-level MAGMA score UMAP not generated; wrote blocker outputs")
} else {
  obj <- readRDS(object_path)
  if (!"umap" %in% names(obj@reductions)) stop("Seurat object lacks a UMAP reduction.")
  meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  celltype_col <- if ("audited_broad_cell_type" %in% names(meta)) "audited_broad_cell_type" else "phase1_cell_type"
  donor_col <- if ("donor_id" %in% names(meta)) "donor_id" else "sample_id"
  sample_col <- if ("sample_id" %in% names(meta)) "sample_id" else donor_col
  umap <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "cell_id")
  setnames(umap, old = names(umap)[2:3], new = c("UMAP_1", "UMAP_2"))
  meta <- merge(meta, umap, by = "cell_id")
  meta[, audited_broad_cell_type := get(celltype_col)]
  meta[, donor_id := get(donor_col)]
  meta[, sample_id := get(sample_col)]

  read_set <- function(path) unique(fread(path, header = FALSE)[[1]])
  sets <- list(
    MAGMA_top50 = read_set("results/gene_sets/magma_top50.txt"),
    MAGMA_top100 = read_set("results/gene_sets/magma_top100.txt"),
    MAGMA_FDR = read_set("results/gene_sets/magma_fdr05.txt"),
    MAGMA_suggestive = read_set("results/gene_sets/magma_suggestive_p1e4.txt")
  )

  expr <- tryCatch(
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
    error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
  )
  score_wide <- meta[, .(cell_id, sample_id, donor_id, disease_status, audited_broad_cell_type, UMAP_1, UMAP_2)]
  score_meta <- rbindlist(lapply(names(sets), function(nm) {
    genes <- intersect(sets[[nm]], rownames(expr))
    if (length(genes) == 0) {
      score_wide[, (paste0(nm, "_score")) := NA_real_]
    } else {
      sub <- as.matrix(expr[genes, , drop = FALSE])
      z <- t(scale(t(sub)))
      z[!is.finite(z)] <- NA_real_
      score_wide[, (paste0(nm, "_score")) := as.numeric(colMeans(z, na.rm = TRUE))]
    }
    data.table(module_name = nm, n_genes_input = length(sets[[nm]]),
               n_genes_detected = length(genes),
               detected_fraction = length(genes) / max(1, length(sets[[nm]])))
  }))

  score_long <- melt(
    score_wide,
    id.vars = c("cell_id", "sample_id", "donor_id", "disease_status", "audited_broad_cell_type", "UMAP_1", "UMAP_2"),
    measure.vars = paste0(names(sets), "_score"),
    variable.name = "module_name",
    value.name = "celllevel_module_score"
  )
  score_long[, module_name := sub("_score$", "", module_name)]
  score_long <- merge(score_long, score_meta, by = "module_name", all.x = TRUE)
  fwrite(score_long, "results/tables/gse231569_celllevel_magma_scores.tsv", sep = "\t")
  fwrite(data.table(
    requested_output = c("gse231569_celllevel_magma_scores.tsv", "figure2_umap_audited_celltypes_v0.1", "figure2_magma_score_umap_v0.1"),
    status = "generated",
    reason = "Audited GSE231569 Seurat object with UMAP reduction was available locally"
  ), "results/tables/gse231569_celllevel_magma_score_blocker_v0.1.tsv", sep = "\t")

  donor_summary <- score_long[, .(
    n_cells = .N,
    mean_score = mean(celllevel_module_score, na.rm = TRUE),
    median_score = median(celllevel_module_score, na.rm = TRUE),
    n_genes_detected = n_genes_detected[1],
    detected_fraction = detected_fraction[1]
  ), by = .(module_name, donor_id, sample_id, disease_status, audited_broad_cell_type)]
  fwrite(donor_summary, "results/tables/gse231569_donor_celltype_magma_scores_v0.2.tsv", sep = "\t")

  label_cell <- function(x) {
    map <- c(
      Collecting_duct_principal = "Collecting duct",
      Fibroblast_stromal = "Fibroblast/stromal",
      Endothelial = "Endothelial",
      Injured_undifferentiated_epithelial = "Injured epithelial",
      Loop_of_Henle_TAL = "Loop/TAL",
      Perivascular_mural_like = "Perivascular/mural-like",
      Pericyte_smooth_muscle = "Perivascular/mural-like"
    )
    unname(ifelse(x %in% names(map), map[x], x))
  }
  score_wide[, cell_label := label_cell(audited_broad_cell_type)]
  score_wide[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]

  set.seed(1)
  plot_cells <- score_wide
  if (nrow(plot_cells) > 22000) plot_cells <- plot_cells[sample(.N, 22000)]
  p_atlas <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_label, alpha = is_tal), size = 0.18) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.38), guide = "none") +
    labs(title = "Audited GSE231569 single-nucleus atlas", x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right",
          panel.grid = element_blank())
  ggsave("results/figures/figure2_umap_audited_celltypes_v0.1.pdf", p_atlas,
         width = 7.4, height = 5.4, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_umap_audited_celltypes_v0.1.png", p_atlas,
         width = 7.4, height = 5.4, units = "in", dpi = 260, bg = "white")

  score_plot <- score_long[module_name == "MAGMA_top50"]
  if (nrow(score_plot) > 22000) score_plot <- score_plot[sample(.N, 22000)]
  qlim <- quantile(score_plot$celllevel_module_score, c(0.02, 0.98), na.rm = TRUE)
  p_score <- ggplot(score_plot, aes(UMAP_1, UMAP_2, color = celllevel_module_score)) +
    geom_point(size = 0.18, alpha = 0.75) +
    scale_color_gradient(low = "#E8ECEE", high = "#3E6672", limits = qlim, oob = scales::squish) +
    labs(title = "MAGMA top 50 module score projected on GSE231569 UMAP",
         x = "UMAP 1", y = "UMAP 2", color = "Score") +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
  ggsave("results/figures/figure2_magma_score_umap_v0.1.pdf", p_score,
         width = 7.4, height = 5.4, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_score_umap_v0.1.png", p_score,
         width = 7.4, height = 5.4, units = "in", dpi = 260, bg = "white")

  donor_summary[, cell_label := label_cell(audited_broad_cell_type)]
  donor_summary[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  donor_summary[, module_label := factor(module_name, levels = names(sets),
                                         labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  p_violin <- ggplot(donor_summary, aes(cell_label, mean_score, fill = is_tal)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "#555555", linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.09, height = 0), size = 1.7,
               shape = 21, color = "#555555", stroke = 0.2) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("TRUE" = "#3E6672", "FALSE" = "#D8D8D8"), labels = c("Other", "Loop/TAL")) +
    labs(title = "Donor-level cell-type MAGMA module scores",
         subtitle = "Each point = donor-level cell-type score",
         x = NULL, y = "Mean per-cell module score", fill = NULL) +
    theme_bw(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, color = "#555555"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave("results/figures/figure2_magma_score_violin_v0.1.pdf", p_violin,
         width = 9.4, height = 6.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_score_violin_v0.1.png", p_violin,
         width = 9.4, height = 6.6, units = "in", dpi = 260, bg = "white")

  message("wrote GSE231569 cell-level MAGMA score UMAP outputs")
}
