suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
obj_dir <- "data/processed/spatial/gse206306/seurat_objects"
defs_path <- file.path(out_dir, "spatial_signature_definitions_v0.1.tsv")
if (!file.exists(defs_path)) stop("Missing signature definitions: ", defs_path)

defs <- read_tsv(defs_path, show_col_types = FALSE)
objects <- list.files(obj_dir, pattern = "_spatial_qc\\.rds$", full.names = TRUE)
if (!length(objects)) stop("No Phase 28B QC Seurat objects found.")

threshold_for <- function(signature) {
  if (grepl("^MAGMA_", signature)) return(10L)
  if (signature == "P1_core") return(3L)
  if (signature == "Loop_TAL") return(3L)
  4L
}

audit_one <- function(path) {
  sample_id <- sub("_spatial_qc\\.rds$", "", basename(path))
  obj <- readRDS(path)
  features <- rownames(obj)
  defs %>%
    group_by(signature) %>%
    summarise(input_genes = list(unique(gene)), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      sample_id = sample_id,
      input_gene_n = length(input_genes),
      detected = list(intersect(input_genes, features)),
      missing = list(setdiff(input_genes, features)),
      detected_gene_n = length(detected),
      detected_gene_fraction = ifelse(input_gene_n > 0, detected_gene_n / input_gene_n, NA_real_),
      detected_genes = paste(detected, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      score_allowed = detected_gene_n >= threshold_for(signature)
    ) %>%
    ungroup() %>%
    select(sample_id, signature, input_gene_n, detected_gene_n, detected_gene_fraction,
           detected_genes, missing_genes, score_allowed)
}

overlap <- bind_rows(lapply(objects, audit_one)) %>%
  arrange(sample_id, signature)

write_tsv(overlap, file.path(out_dir, "spatial_signature_gene_overlap_v0.1.tsv"))
message("Wrote overlap audit for ", length(unique(overlap$sample_id)), " samples.")
