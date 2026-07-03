suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
table_dir <- "results/tables"
figure_dir <- "results/figures"
scrna_table_dir <- "results/scrna/gse231569/tables"
scrna_figure_dir <- "results/scrna/gse231569/figures"
config_dir <- "config"
for (d in c(table_dir, figure_dir, scrna_table_dir, scrna_figure_dir, config_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

obj <- readRDS(obj_path)
if (!"sample_id" %in% colnames(obj@meta.data)) obj$sample_id <- obj$orig.ident
if (!"donor_id" %in% colnames(obj@meta.data)) obj$donor_id <- obj$sample_id
if (!"phase1_cell_type" %in% colnames(obj@meta.data)) {
  stop("phase1_cell_type is missing; run phase1_gse231569_quick_projection.R first")
}

mat <- GetAssayData(obj, assay = "RNA", layer = "data")
clusters <- as.character(obj$seurat_clusters)
cell_types <- as.character(obj$phase1_cell_type)

marker_sets <- list(
  Proximal_tubule = c("LRP2", "SLC34A1", "CUBN", "ALDOB"),
  Loop_of_Henle_TAL = c("UMOD", "SLC12A1", "CLDN10", "KCNJ1", "CLDN16"),
  Collecting_duct_principal = c("AQP2", "AQP3", "SCNN1G"),
  Intercalated = c("ATP6V1B1", "SLC4A1", "FOXI1"),
  Injured_undifferentiated_epithelial = c("KRT8", "KRT18", "PROM1", "HAVCR1", "VCAM1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
  Fibroblast_stromal = c("DCN", "LUM", "COL1A1", "COL1A2"),
  Pericyte_smooth_muscle = c("RGS5", "ACTA2", "TAGLN", "MYH11", "PDGFRB", "CSPG4", "MCAM"),
  Macrophage_myeloid = c("LST1", "C1QA", "C1QB", "TYROBP", "CD68"),
  T_cell = c("CD3D", "TRAC", "IL7R"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  B_cell_plasma = c("MS4A1", "CD79A", "MZB1", "JCHAIN")
)

audit_one_grouping <- function(group_values, group_name) {
  groups <- sort(unique(group_values))
  rbindlist(lapply(names(marker_sets), function(marker_set) {
    genes <- intersect(marker_sets[[marker_set]], rownames(mat))
    rbindlist(lapply(groups, function(grp) {
      cells <- colnames(obj)[group_values == grp]
      if (length(genes) == 0 || length(cells) == 0) {
        mean_expr <- NA_real_
        pct_any <- NA_real_
      } else {
        sub <- mat[genes, cells, drop = FALSE]
        mean_expr <- mean(Matrix::colMeans(sub))
        pct_any <- mean(Matrix::colSums(sub > 0) > 0)
      }
      data.table(
        grouping = group_name,
        group = grp,
        candidate_cell_type = marker_set,
        marker_genes_requested = paste(marker_sets[[marker_set]], collapse = ","),
        marker_genes_present = paste(genes, collapse = ","),
        n_markers_present = length(genes),
        mean_marker_expression = mean_expr,
        pct_cells_with_any_marker = pct_any,
        n_cells = length(cells)
      )
    }))
  }))
}

audit <- rbind(
  audit_one_grouping(clusters, "seurat_cluster"),
  audit_one_grouping(cell_types, "phase1_cell_type")
)
fwrite(audit, file.path(table_dir, "gse231569_marker_audit.tsv"), sep = "\t")
fwrite(audit, file.path(scrna_table_dir, "gse231569_marker_audit.tsv"), sep = "\t")

cluster_audit <- audit[grouping == "seurat_cluster"]
cluster_best <- cluster_audit[
  order(group, -mean_marker_expression, -pct_cells_with_any_marker)
][
  ,
  .SD[1],
  by = group
]
setnames(cluster_best, "group", "cluster")
cluster_current <- unique(data.table(cluster = clusters, current_phase1_cell_type = cell_types))
cluster_current <- cluster_current[, .(
  current_phase1_cell_type = paste(sort(unique(current_phase1_cell_type)), collapse = ";")
), by = cluster]
cluster_assignment <- merge(cluster_current, cluster_best, by = "cluster", all.x = TRUE)
immune_best <- cluster_audit[
  candidate_cell_type %in% c("Macrophage_myeloid", "T_cell", "Mast_cell", "B_cell_plasma")
][
  order(group, -mean_marker_expression, -pct_cells_with_any_marker)
][
  ,
  .SD[1],
  by = group
][
  ,
  .(
    cluster = group,
    top_immune_marker_set = candidate_cell_type,
    top_immune_mean_marker_expression = mean_marker_expression,
    top_immune_pct_cells_with_any_marker = pct_cells_with_any_marker
  )
]
cluster_assignment <- merge(cluster_assignment, immune_best, by = "cluster", all.x = TRUE)
cluster_assignment[, audit_flag := fifelse(
  grepl("Macrophage|T_cell|Mast_cell|B_cell", candidate_cell_type) &
    !grepl("Macrophage|T_cell|Mast|B_cell|Immune", current_phase1_cell_type),
  "possible_immune_underannotation",
  fifelse(
    candidate_cell_type != current_phase1_cell_type,
    "marker_top_label_differs_from_phase1",
    "phase1_label_consistent_with_top_marker_set"
  )
)]
cluster_assignment[, immune_review_flag := fifelse(
  top_immune_mean_marker_expression >= 0.2 & top_immune_pct_cells_with_any_marker >= 0.4,
  "immune_marker_signal_requires_manual_review",
  "no_strong_cluster_level_immune_signal"
)]
cluster_assignment[, confidence := fifelse(
  n_cells < 50,
  "exploratory_low_cell_count",
  fifelse(mean_marker_expression > 0.25 & pct_cells_with_any_marker > 0.25, "medium_high", "low_medium")
)]
fwrite(
  cluster_assignment,
  file.path(table_dir, "gse231569_cluster_marker_assignment_audit.tsv"),
  sep = "\t"
)
fwrite(
  cluster_assignment,
  file.path(scrna_table_dir, "gse231569_cluster_marker_assignment_audit.tsv"),
  sep = "\t"
)

cell_counts <- as.data.table(obj@meta.data)[
  ,
  .(
    n_cells = .N,
    n_donors = uniqueN(donor_id),
    n_samples = uniqueN(sample_id),
    median_nCount_RNA = as.numeric(median(nCount_RNA, na.rm = TRUE)),
    median_nFeature_RNA = as.numeric(median(nFeature_RNA, na.rm = TRUE)),
    median_percent_mt = as.numeric(median(percent.mt, na.rm = TRUE))
  ),
  by = .(broad_cell_type = phase1_cell_type)
][order(-n_cells)]
fwrite(cell_counts, file.path(table_dir, "gse231569_cell_counts.tsv"), sep = "\t")
fwrite(cell_counts, file.path(scrna_table_dir, "gse231569_cell_counts.tsv"), sep = "\t")

donor_counts <- as.data.table(obj@meta.data)[
  ,
  .(n_cells = .N),
  by = .(donor_id, sample_id, disease_status, broad_cell_type = phase1_cell_type)
][order(donor_id, broad_cell_type)]
fwrite(donor_counts, file.path(table_dir, "gse231569_donor_cell_counts.tsv"), sep = "\t")
fwrite(donor_counts, file.path(scrna_table_dir, "gse231569_donor_cell_counts.tsv"), sep = "\t")

all_markers <- intersect(unique(unlist(marker_sets)), rownames(obj))
immune_markers <- intersect(unlist(marker_sets[c("Macrophage_myeloid", "T_cell", "Mast_cell", "B_cell_plasma")]), rownames(obj))
mural_markers <- intersect(c(marker_sets$Pericyte_smooth_muscle, marker_sets$Fibroblast_stromal, marker_sets$Endothelial), rownames(obj))
transport_markers <- intersect(c(marker_sets$Loop_of_Henle_TAL, marker_sets$Collecting_duct_principal, marker_sets$Intercalated), rownames(obj))

pdf(file.path(figure_dir, "gse231569_marker_dotplot_by_cluster.pdf"), width = 16, height = 8)
print(DotPlot(obj, features = all_markers, group.by = "seurat_clusters") + RotatedAxis())
dev.off()
file.copy(file.path(figure_dir, "gse231569_marker_dotplot_by_cluster.pdf"), file.path(scrna_figure_dir, "gse231569_marker_dotplot_by_cluster.pdf"), overwrite = TRUE)

if (length(immune_markers) > 0) {
  pdf(file.path(figure_dir, "gse231569_immune_marker_featureplots.pdf"), width = 14, height = 10)
  print(FeaturePlot(obj, features = immune_markers, ncol = 4))
  dev.off()
  file.copy(file.path(figure_dir, "gse231569_immune_marker_featureplots.pdf"), file.path(scrna_figure_dir, "gse231569_immune_marker_featureplots.pdf"), overwrite = TRUE)
}

if (length(mural_markers) > 0) {
  pdf(file.path(figure_dir, "gse231569_mural_audit_dotplot.pdf"), width = 12, height = 6)
  print(DotPlot(obj, features = mural_markers, group.by = "seurat_clusters") + RotatedAxis())
  dev.off()
  file.copy(file.path(figure_dir, "gse231569_mural_audit_dotplot.pdf"), file.path(scrna_figure_dir, "gse231569_mural_audit_dotplot.pdf"), overwrite = TRUE)
}

if (length(transport_markers) > 0) {
  pdf(file.path(figure_dir, "gse231569_transport_marker_audit.pdf"), width = 12, height = 6)
  print(DotPlot(obj, features = transport_markers, group.by = "seurat_clusters") + RotatedAxis())
  dev.off()
  file.copy(file.path(figure_dir, "gse231569_transport_marker_audit.pdf"), file.path(scrna_figure_dir, "gse231569_transport_marker_audit.pdf"), overwrite = TRUE)
}

cluster_to_celltype <- cluster_assignment[, .(
  cluster,
  broad_cell_type = current_phase1_cell_type,
  fine_cell_type = current_phase1_cell_type,
  confidence,
  main_positive_markers = marker_genes_present,
  negative_markers_checked = "PECAM1,VWF,DCN,COL1A1,LST1,C1QA,CD3D,TRAC,TPSAB1,KIT",
  audit_top_marker_set = candidate_cell_type,
  top_immune_marker_set,
  top_immune_mean_marker_expression,
  top_immune_pct_cells_with_any_marker,
  immune_review_flag,
  notes = audit_flag
)]
fwrite(cluster_to_celltype, file.path(config_dir, "gse231569_cluster_to_celltype.tsv"), sep = "\t")

pdf(file.path(figure_dir, "gse231569_umap_audited_broad_celltype.pdf"), width = 8, height = 7)
print(DimPlot(obj, group.by = "phase1_cell_type", label = TRUE, repel = TRUE) + NoLegend())
dev.off()
file.copy(file.path(figure_dir, "gse231569_umap_audited_broad_celltype.pdf"), file.path(scrna_figure_dir, "gse231569_umap_audited_broad_celltype.pdf"), overwrite = TRUE)

saveRDS(obj, "results/scrna/gse231569/objects/gse231569_annotated_audited.rds")

cat("wrote\t", file.path(table_dir, "gse231569_marker_audit.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "gse231569_cluster_marker_assignment_audit.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "gse231569_cell_counts.tsv"), "\n", sep = "")
cat("wrote\t", file.path(table_dir, "gse231569_donor_cell_counts.tsv"), "\n", sep = "")
cat("wrote\t", file.path(config_dir, "gse231569_cluster_to_celltype.tsv"), "\n", sep = "")
