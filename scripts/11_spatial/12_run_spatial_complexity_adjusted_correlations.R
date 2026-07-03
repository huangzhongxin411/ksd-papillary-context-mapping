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

residualize <- function(meta, score_col) {
  covars <- c("log_nCount", "nFeature_Spatial")
  if ("percent.mt" %in% colnames(meta) && any(is.finite(meta$percent.mt)) && sd(meta$percent.mt, na.rm = TRUE) > 0) {
    covars <- c(covars, "percent.mt")
  }
  dat <- meta %>%
    mutate(log_nCount = log1p(nCount_Spatial)) %>%
    select(all_of(c(score_col, covars))) %>%
    filter(if_all(everything(), is.finite))
  if (nrow(dat) < 20) return(list(resid = NULL, index = integer(), covars = paste(covars, collapse = "+")))
  names(dat)[1] <- "score"
  fit <- lm(reformulate(covars, response = "score"), data = dat)
  list(resid = residuals(fit), index = as.integer(rownames(dat)), covars = paste(covars, collapse = "+"))
}

rows <- list()
for (path in objects) {
  sample_id <- sub("_spatial_scores_v0\\.1\\.rds$", "", basename(path))
  obj <- readRDS(path)
  meta <- obj@meta.data
  meta$.spot_index <- seq_len(nrow(meta))
  residual_cache <- list()

  for (sig in unique(c(risk_signatures, context_signatures))) {
    col <- paste0("score_", sig)
    if (col %in% colnames(meta)) {
      tmp <- meta
      rownames(tmp) <- as.character(tmp$.spot_index)
      residual_cache[[sig]] <- residualize(tmp, col)
    }
  }

  for (risk in risk_signatures) {
    for (context in context_signatures) {
      r1 <- residual_cache[[risk]]
      r2 <- residual_cache[[context]]
      if (is.null(r1) || is.null(r2) || is.null(r1$resid) || is.null(r2$resid)) {
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = NA_integer_, rho_adjusted = NA_real_, p_adjusted_descriptive = NA_real_,
          covariates_used = NA_character_, notes = "score_missing_or_residualization_failed"
        )
        next
      }
      common <- intersect(r1$index, r2$index)
      x <- r1$resid[match(common, r1$index)]
      y <- r2$resid[match(common, r2$index)]
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 10) {
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = sum(keep), rho_adjusted = NA_real_, p_adjusted_descriptive = NA_real_,
          covariates_used = r1$covars, notes = "too_few_spots"
        )
      } else {
        ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
        rows[[length(rows) + 1]] <- tibble(
          sample_id = sample_id, risk_signature = risk, context_signature = context,
          n_spots = sum(keep), rho_adjusted = unname(ct$estimate),
          p_adjusted_descriptive = ct$p.value, covariates_used = r1$covars,
          notes = "residual Spearman; spot-level p descriptive only"
        )
      }
    }
  }
}

write_tsv(bind_rows(rows), file.path(out_dir, "spatial_score_correlations_adjusted_v0.1.tsv"))
