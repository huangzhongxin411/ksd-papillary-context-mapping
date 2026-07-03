suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

raw_dir <- "data/raw/GSE231569/raw_mtx"
out_dir <- "results"
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
processed_dir <- "data/processed/GSE231569"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.frame(
  sample_id = c("GSM7290910", "GSM7290912", "GSM7290913", "GSM7290914"),
  prefix = c(
    "GSM7290910_3778_K2200020_1_D1_N1",
    "GSM7290912_3693_K2100055_1_D1_N1",
    "GSM7290913_IU-KRP-429_D1_N1",
    "GSM7290914_IU-KRP-473_D1_N1"
  ),
  disease_status = c("healthy_control", "healthy_control", "stone_disease", "stone_disease"),
  tissue = "Papilla",
  stringsAsFactors = FALSE
)

read_one <- function(prefix, sample_id, disease_status) {
  mat <- ReadMtx(
    mtx = file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz")),
    cells = file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz")),
    features = file.path(raw_dir, paste0(prefix, "_features.tsv.gz")),
    feature.column = 2,
    unique.features = TRUE
  )
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  obj$sample_id <- sample_id
  obj$disease_status <- disease_status
  obj$tissue <- "Papilla"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

objs <- lapply(seq_len(nrow(samples)), function(i) {
  read_one(samples$prefix[i], samples$sample_id[i], samples$disease_status[i])
})

obj <- Reduce(function(x, y) merge(x, y), objs)
qc_before <- data.frame(
  metric = c("cells_before_qc", "features"),
  value = c(ncol(obj), nrow(obj))
)

obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 20)
qc_after <- data.frame(
  metric = c("cells_after_qc", "median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt"),
  value = c(
    ncol(obj),
    median(obj$nFeature_RNA),
    median(obj$nCount_RNA),
    median(obj$percent.mt)
  )
)
write.table(
  rbind(qc_before, qc_after),
  file.path(table_dir, "phase1_gse231569_qc_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- JoinLayers(obj, assay = "RNA")
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:25, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:25, verbose = FALSE)

marker_sets <- list(
  Proximal_tubule = c("LRP2", "SLC34A1", "CUBN", "ALDOB", "SLC5A12"),
  Loop_of_Henle_TAL = c("UMOD", "SLC12A1", "CLDN10", "KCNJ1"),
  Collecting_duct_principal = c("AQP2", "AQP3", "SCNN1G", "PDE1A"),
  Intercalated_cell = c("ATP6V1B1", "SLC4A1", "FOXI1", "ATP6V0D2"),
  Injured_undifferentiated_epithelial = c("KRT8", "KRT18", "PROM1", "VCAM1", "MMP7"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
  Fibroblast_stromal = c("DCN", "LUM", "COL1A1", "COL1A2"),
  Macrophage_myeloid = c("LST1", "C1QA", "C1QB", "CD68", "TYROBP"),
  T_cell = c("CD3D", "CD3E", "TRAC", "IL7R"),
  Pericyte_smooth_muscle = c("RGS5", "ACTA2", "PDGFRB", "MYH11")
)

data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
score_one <- function(genes) {
  present <- intersect(genes, rownames(data_mat))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(data_mat)))
  }
  Matrix::colMeans(data_mat[present, , drop = FALSE])
}
for (nm in names(marker_sets)) {
  obj[[paste0("marker_", nm)]] <- score_one(marker_sets[[nm]])
}

marker_cols <- paste0("marker_", names(marker_sets))
cluster_scores <- aggregate(
  obj@meta.data[, marker_cols, drop = FALSE],
  by = list(seurat_clusters = obj$seurat_clusters),
  FUN = mean,
  na.rm = TRUE
)
best <- apply(cluster_scores[, marker_cols, drop = FALSE], 1, function(x) {
  names(marker_sets)[which.max(x)]
})
cluster_anno <- data.frame(
  seurat_clusters = cluster_scores$seurat_clusters,
  phase1_cell_type = best,
  cluster_scores[, marker_cols, drop = FALSE],
  check.names = FALSE
)
write.table(
  cluster_anno,
  file.path(table_dir, "phase1_gse231569_cluster_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
anno_map <- setNames(cluster_anno$phase1_cell_type, cluster_anno$seurat_clusters)
obj$phase1_cell_type <- unname(anno_map[as.character(obj$seurat_clusters)])

ranked <- read.delim("results/tables/phase1_candidate_genes_ranked_protein_coding.tsv")
top50 <- scan("results/tables/phase1_candidate_genes_top50.txt", what = character(), quiet = TRUE)
top100 <- scan("results/tables/phase1_candidate_genes_top100.txt", what = character(), quiet = TRUE)

risk_score <- function(genes) {
  present <- intersect(genes, rownames(data_mat))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(data_mat)))
  }
  Matrix::colMeans(data_mat[present, , drop = FALSE])
}
obj$risk_top50_score <- risk_score(top50)
obj$risk_top100_score <- risk_score(top100)

rank_lookup <- ranked[, c("gene", "rank_in_candidate_list")]
projection_genes <- intersect(ranked$gene[ranked$rank_in_candidate_list <= 100], rownames(data_mat))
loc_rows <- lapply(projection_genes, function(g) {
  vals <- data_mat[g, ]
  pct <- tapply(vals > 0, obj$phase1_cell_type, mean)
  avg <- tapply(vals, obj$phase1_cell_type, mean)
  data.frame(
    gene = g,
    cell_type = names(avg),
    avg_expression = as.numeric(avg),
    pct_expressed = as.numeric(pct),
    stringsAsFactors = FALSE
  )
})
loc <- do.call(rbind, loc_rows)
loc <- merge(loc, rank_lookup, by = "gene", all.x = TRUE)
loc$module_group <- ifelse(loc$rank_in_candidate_list <= 50, "top50", "top100")
loc <- loc[order(loc$rank_in_candidate_list, loc$cell_type), ]
write.table(
  loc,
  file.path(table_dir, "phase1_celltype_localization.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

module_summary <- aggregate(
  obj@meta.data[, c("risk_top50_score", "risk_top100_score")],
  by = list(cell_type = obj$phase1_cell_type, disease_status = obj$disease_status),
  FUN = mean,
  na.rm = TRUE
)
write.table(
  module_summary,
  file.path(table_dir, "phase1_module_score_by_celltype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file.path(figure_dir, "phase1_gse231569_umap_celltype.pdf"), width = 8, height = 6)
print(DimPlot(obj, group.by = "phase1_cell_type", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

top_dot <- intersect(ranked$gene[ranked$rank_in_candidate_list <= 30], rownames(obj))
pdf(file.path(figure_dir, "phase1_scrna_dotplot_top30.pdf"), width = 12, height = 5)
print(DotPlot(obj, features = top_dot, group.by = "phase1_cell_type") + RotatedAxis())
dev.off()

pdf(file.path(figure_dir, "phase1_scrna_module_score.pdf"), width = 8, height = 5)
print(VlnPlot(obj, features = c("risk_top50_score", "risk_top100_score"), group.by = "phase1_cell_type", pt.size = 0, stack = FALSE) + NoLegend())
dev.off()

saveRDS(obj, file.path(processed_dir, "phase1_gse231569_quick_seurat.rds"))

cat("cells_after_qc\t", ncol(obj), "\n", sep = "")
cat("projection_genes_present_top100\t", length(projection_genes), "\n", sep = "")
cat("outputs\t", file.path(table_dir, "phase1_celltype_localization.tsv"), "\n", sep = "")
