suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
processed_dir <- "data/processed/gse73680"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

gene_matrix_path <- file.path(processed_dir, "gse73680_gene_expression_matrix.tsv.gz")
p1 <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
tal <- c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2")
injury <- c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "HAVCR1", "CCL2", "CCL7", "VCAM1", "KRT8", "KRT18")
read_gene_set <- function(path, fallback = character()) if (file.exists(path)) fread(path, header = FALSE)[[1]] else fallback
sets <- list(
  P1_core_TAL_candidates = p1,
  MAGMA_top50 = read_gene_set("results/gene_sets/magma_top50.txt"),
  MAGMA_top100 = read_gene_set("results/gene_sets/magma_top100.txt"),
  MAGMA_top200 = read_gene_set("results/gene_sets/magma_top200.txt"),
  MAGMA_suggestive = read_gene_set("results/gene_sets/magma_suggestive_p1e4.txt"),
  MAGMA_FDR = read_gene_set("results/gene_sets/magma_fdr05.txt"),
  TAL_marker_set = tal,
  injury_remodeling_marker_set = injury
)

if (!file.exists(gene_matrix_path)) {
  out <- rbindlist(lapply(names(sets), function(nm) {
    genes <- unique(sets[[nm]])
    data.table(gene_set = nm, n_input_genes = length(genes), n_detected = 0L, detected_fraction = 0,
               detected_genes = "", missing_genes = paste(genes, collapse = ";"),
               usable_for_single_gene_analysis = FALSE, usable_for_module_score = "do_not_use_for_module_score",
               notes = "Gene expression matrix missing; availability not assessable.")
  }))
  fwrite(out, file.path(table_dir, "gse73680_gene_availability.tsv"), sep = "\t")
  fwrite(out[gene_set == "P1_core_TAL_candidates"], file.path(table_dir, "gse73680_p1_gene_availability.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
mat <- fread(gene_matrix_path, select = "gene")
available <- unique(mat$gene)
out <- rbindlist(lapply(names(sets), function(nm) {
  genes <- unique(sets[[nm]])
  detected <- intersect(genes, available)
  frac <- length(detected) / max(1, length(genes))
  data.table(gene_set = nm, n_input_genes = length(genes), n_detected = length(detected), detected_fraction = frac,
             detected_genes = paste(detected, collapse = ";"),
             missing_genes = paste(setdiff(genes, detected), collapse = ";"),
             usable_for_single_gene_analysis = if (nm == "P1_core_TAL_candidates") length(detected) > 0 else NA,
             usable_for_module_score = fifelse(frac >= 0.70, "TRUE", fifelse(frac >= 0.40, "exploratory_module_only", "do_not_use_for_module_score")),
             notes = "Availability assessed against gene-level expression matrix.")
}))
fwrite(out, file.path(table_dir, "gse73680_gene_availability.tsv"), sep = "\t")
fwrite(out[gene_set == "P1_core_TAL_candidates"], file.path(table_dir, "gse73680_p1_gene_availability.tsv"), sep = "\t")
message("wrote GSE73680 gene availability outputs")
