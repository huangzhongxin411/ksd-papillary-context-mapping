suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_gene_list <- function(path, fallback = character()) {
  if (!file.exists(path)) return(unique(fallback))
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x) & !grepl("^#", x)]
  unique(x)
}

signature_sources <- list(
  MAGMA_top50 = list(
    genes = read_gene_list("results/gene_sets/magma_top50.txt"),
    source = "project_MAGMA_gene_set",
    notes = "MAGMA top 50 genes from project gene-set file"
  ),
  MAGMA_top100 = list(
    genes = read_gene_list("results/gene_sets/magma_top100.txt"),
    source = "project_MAGMA_gene_set",
    notes = "MAGMA top 100 genes from project gene-set file"
  ),
  MAGMA_FDR = list(
    genes = read_gene_list("results/gene_sets/magma_fdr05.txt"),
    source = "project_MAGMA_gene_set",
    notes = "MAGMA FDR-significant genes from project gene-set file"
  ),
  MAGMA_suggestive = list(
    genes = read_gene_list("results/gene_sets/magma_suggestive_p1e4.txt"),
    source = "project_MAGMA_gene_set",
    notes = "MAGMA suggestive genes from project gene-set file"
  ),
  P1_core = list(
    genes = c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2"),
    source = "phase28c_protocol",
    notes = "Six P1 candidate genes; interpreted gene-by-gene only"
  ),
  Loop_TAL = list(
    genes = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CASR"),
    source = "phase28c_protocol",
    notes = "Loop/TAL cellular-context markers"
  ),
  Collecting_duct = list(
    genes = c("AQP2", "AQP3", "SCNN1G", "ATP6V1B1", "SLC4A1", "FOXI1"),
    source = "phase28c_protocol",
    notes = "Collecting duct markers"
  ),
  Injured_epithelial = list(
    genes = c("KRT8", "KRT18", "PROM1", "HAVCR1", "VCAM1", "LCN2"),
    source = "phase28c_protocol",
    notes = "Injured epithelial markers"
  ),
  Injury_remodeling = list(
    genes = c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "CCL2", "CCL7", "HAVCR1"),
    source = "phase28c_protocol",
    notes = "Injury/remodeling markers"
  ),
  Fibrosis_ECM = list(
    genes = c("COL1A1", "COL1A2", "COL3A1", "FN1", "DCN", "LUM", "ACTA2"),
    source = "phase28c_protocol",
    notes = "Fibrosis and extracellular matrix markers"
  ),
  Immune_myeloid = list(
    genes = c("PTPRC", "LST1", "C1QA", "C1QB", "CD68", "TYROBP"),
    source = "phase28c_protocol",
    notes = "Immune/myeloid markers"
  ),
  Endothelial = list(
    genes = c("PECAM1", "VWF", "KDR", "EMCN"),
    source = "phase28c_protocol",
    notes = "Endothelial markers"
  ),
  Ion_mineral_handling = list(
    genes = c("UMOD", "SLC12A1", "CLDN10", "CLDN14", "CASR", "SLC34A1", "SLC34A3", "TRPV5"),
    source = "phase28c_protocol",
    notes = "Ion/mineral handling marker set"
  )
)

defs <- bind_rows(lapply(names(signature_sources), function(sig) {
  item <- signature_sources[[sig]]
  tibble(
    signature = sig,
    gene = unique(item$genes),
    source = item$source,
    notes = item$notes
  )
})) %>%
  filter(!is.na(gene), nzchar(gene)) %>%
  distinct(signature, gene, .keep_all = TRUE)

write_tsv(defs, file.path(out_dir, "spatial_signature_definitions_v0.1.tsv"))
message("Wrote ", nrow(defs), " signature-gene rows.")
