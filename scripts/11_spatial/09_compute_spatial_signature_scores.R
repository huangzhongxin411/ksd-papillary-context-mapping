suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(dplyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
defs <- read_tsv(file.path(out_dir, "spatial_signature_definitions_v0.1.tsv"), show_col_types = FALSE)
overlap <- read_tsv(file.path(out_dir, "spatial_signature_gene_overlap_v0.1.tsv"), show_col_types = FALSE)
objects <- list.files(obj_dir, pattern = "_spatial_qc\\.rds$", full.names = TRUE)

get_data_matrix <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "Spatial", layer = "data"),
    error = function(e) GetAssayData(obj, assay = "Spatial", slot = "data")
  )
}

mean_z_score <- function(mat) {
  dense <- as.matrix(mat)
  gene_means <- rowMeans(dense)
  gene_sds <- apply(dense, 1, sd)
  keep <- is.finite(gene_sds) & gene_sds > 0
  if (!any(keep)) return(rep(NA_real_, ncol(dense)))
  z <- sweep(sweep(dense[keep, , drop = FALSE], 1, gene_means[keep], "-"), 1, gene_sds[keep], "/")
  colMeans(z, na.rm = TRUE)
}

qc_rows <- list()
for (path in objects) {
  sample_id <- sub("_spatial_qc\\.rds$", "", basename(path))
  message("Scoring ", sample_id)
  obj <- readRDS(path)
  DefaultAssay(obj) <- "Spatial"
  obj <- NormalizeData(obj, assay = "Spatial", normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  }
  data_mat <- get_data_matrix(obj)
  sample_overlap <- overlap %>% filter(sample_id == !!sample_id)

  for (i in seq_len(nrow(sample_overlap))) {
    sig <- sample_overlap$signature[i]
    allowed <- isTRUE(sample_overlap$score_allowed[i])
    genes <- strsplit(sample_overlap$detected_genes[i], ";", fixed = TRUE)[[1]]
    genes <- genes[nzchar(genes)]
    score_col <- paste0("score_", sig)

    if (allowed && length(genes) > 0) {
      score <- mean_z_score(data_mat[genes, , drop = FALSE])
      obj[[score_col]] <- score
      qc_rows[[length(qc_rows) + 1]] <- tibble(
        sample_id = sample_id,
        signature = sig,
        score_allowed = TRUE,
        n_genes_used = length(genes),
        score_mean = mean(score, na.rm = TRUE),
        score_sd = sd(score, na.rm = TRUE),
        score_min = min(score, na.rm = TRUE),
        score_q25 = unname(quantile(score, 0.25, na.rm = TRUE)),
        score_median = median(score, na.rm = TRUE),
        score_q75 = unname(quantile(score, 0.75, na.rm = TRUE)),
        score_max = max(score, na.rm = TRUE),
        notes = "deterministic mean-z score"
      )
    } else {
      qc_rows[[length(qc_rows) + 1]] <- tibble(
        sample_id = sample_id,
        signature = sig,
        score_allowed = FALSE,
        n_genes_used = length(genes),
        score_mean = NA_real_,
        score_sd = NA_real_,
        score_min = NA_real_,
        score_q25 = NA_real_,
        score_median = NA_real_,
        score_q75 = NA_real_,
        score_max = NA_real_,
        notes = "not scored: gene-overlap threshold not met"
      )
    }
  }
  saveRDS(obj, file.path(obj_dir, paste0(sample_id, "_spatial_scores_v0.1.rds")))
}

qc <- bind_rows(qc_rows) %>% arrange(sample_id, signature)
write_tsv(qc, file.path(out_dir, "spatial_module_score_qc_v0.1.tsv"))
message("Wrote score QC: ", nrow(qc), " rows.")
