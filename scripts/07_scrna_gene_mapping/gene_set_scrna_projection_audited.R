suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
gene_set_dir <- if (length(args) >= 1) args[[1]] else "results/gene_sets"
out_prefix <- if (length(args) >= 2) args[[2]] else "custom_gene_sets"

obj_path <- "data/processed/GSE231569/phase1_gse231569_quick_seurat.rds"
cfg_path <- "config/gse231569_cluster_to_celltype_audited.tsv"
table_dir <- "results/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

gene_set_files <- list.files(gene_set_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(gene_set_files) == 0) {
  stop("No .txt gene set files found in ", gene_set_dir)
}

obj <- readRDS(obj_path)
cfg <- fread(cfg_path)
cfg[, cluster := as.character(cluster)]
cluster_to_broad <- setNames(cfg$broad_cell_type, cfg$cluster)
cluster_to_conf <- setNames(cfg$annotation_confidence, cfg$cluster)
cluster_to_immune <- setNames(cfg$immune_review_flag, cfg$cluster)
cluster_to_mural <- setNames(cfg$mural_review_flag, cfg$cluster)
obj$audited_broad_cell_type <- unname(cluster_to_broad[as.character(obj$seurat_clusters)])
obj$audited_annotation_confidence <- unname(cluster_to_conf[as.character(obj$seurat_clusters)])
obj$audited_immune_review_flag <- unname(cluster_to_immune[as.character(obj$seurat_clusters)])
obj$audited_mural_review_flag <- unname(cluster_to_mural[as.character(obj$seurat_clusters)])
if (!"donor_id" %in% colnames(obj@meta.data)) obj$donor_id <- obj$sample_id

mat <- GetAssayData(obj, assay = "RNA", layer = "data")
meta <- as.data.table(obj@meta.data)

score_genes <- function(genes) {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[present, , drop = FALSE])
}

score_cols <- character()
set_info <- lapply(gene_set_files, function(path) {
  label <- tools::file_path_sans_ext(basename(path))
  genes <- scan(path, what = character(), quiet = TRUE)
  col <- paste0(label, "_score")
  meta[[col]] <<- score_genes(genes)
  score_cols <<- c(score_cols, col)
  data.table(
    gene_set = label,
    file = path,
    n_genes_requested = length(unique(genes)),
    n_genes_present = length(intersect(unique(genes), rownames(mat)))
  )
})
fwrite(rbindlist(set_info), file.path(table_dir, paste0(out_prefix, "_gene_set_presence.tsv")), sep = "\t")

donor_score <- meta[
  ,
  c(
    .(
      n_cells = .N,
      annotation_confidence = paste(sort(unique(audited_annotation_confidence)), collapse = ";"),
      immune_review_cells = sum(audited_immune_review_flag %in% TRUE),
      mural_review_cells = sum(audited_mural_review_flag %in% TRUE)
    ),
    as.list(vapply(score_cols, function(x) mean(get(x), na.rm = TRUE), numeric(1)))
  ),
  by = .(donor_id, sample_id, disease_status, audited_broad_cell_type)
]
fwrite(donor_score, file.path(table_dir, paste0(out_prefix, "_module_score_by_donor_celltype.tsv")), sep = "\t")

cell_score <- meta[
  ,
  c(
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_id),
      annotation_confidence = paste(sort(unique(audited_annotation_confidence)), collapse = ";"),
      immune_review_cells = sum(audited_immune_review_flag %in% TRUE),
      mural_review_cells = sum(audited_mural_review_flag %in% TRUE)
    ),
    as.list(vapply(score_cols, function(x) mean(get(x), na.rm = TRUE), numeric(1)))
  ),
  by = audited_broad_cell_type
]
fwrite(cell_score, file.path(table_dir, paste0(out_prefix, "_module_score_by_celltype.tsv")), sep = "\t")

cat("wrote\t", file.path(table_dir, paste0(out_prefix, "_module_score_by_donor_celltype.tsv")), "\n", sep = "")
cat("wrote\t", file.path(table_dir, paste0(out_prefix, "_module_score_by_celltype.tsv")), "\n", sep = "")
