suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
risk_signatures <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core")
context_signatures <- c("Loop_TAL", "Collecting_duct", "Injured_epithelial", "Injury_remodeling",
                        "Fibrosis_ECM", "Immune_myeloid", "Endothelial", "Ion_mineral_handling")

objects <- list.files(obj_dir, pattern = "_spatial_scores_v0\\.1\\.rds$", full.names = TRUE)
if (!length(objects)) stop("No scored objects found.")

rows <- list()
for (path in objects) {
  sample_id <- sub("_spatial_scores_v0\\.1\\.rds$", "", basename(path))
  obj <- readRDS(path)
  meta <- obj@meta.data
  for (risk in risk_signatures) {
    risk_col <- paste0("score_", risk)
    for (context in context_signatures) {
      context_col <- paste0("score_", context)
      if (!all(c(risk_col, context_col) %in% colnames(meta))) {
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = NA_integer_, rho_raw = NA_real_, p_raw = NA_real_,
          direction = NA_character_, notes = "score_missing"
        )
        next
      }
      x <- meta[[risk_col]]
      y <- meta[[context_col]]
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 10) {
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = sum(keep), rho_raw = NA_real_, p_raw = NA_real_,
          direction = NA_character_, notes = "too_few_spots"
        )
      } else {
        ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
        rho <- unname(ct$estimate)
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = sum(keep), rho_raw = rho, p_raw = ct$p.value,
          direction = ifelse(rho > 0, "positive", ifelse(rho < 0, "negative", "zero")),
          notes = "spot-level p descriptive only; spots are not independent replicates"
        )
      }
    }
  }
}

write_tsv(bind_rows(rows), file.path(out_dir, "spatial_score_correlations_by_sample_v0.1.tsv"))
