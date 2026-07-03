suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

candidate_patterns <- "\\.(rds|RDS|rda|RData|h5ad|h5seurat|loom|mtx|mtx.gz|h5)$|barcodes|features|genes|matrix"
paths <- list.files(c("data", "results", "config"), recursive = TRUE, full.names = TRUE, all.files = FALSE)
hits <- paths[grepl("GSE231569|gse231569", paths) & grepl(candidate_patterns, paths, ignore.case = TRUE)]
embedding_hits <- paths[grepl("GSE231569|gse231569", paths) & grepl("umap|embedding|coord|cell_metadata", paths, ignore.case = TRUE)]

audit <- data.table(
  resource = c("cell_level_expression_object", "cell_level_embedding_table", "audited_annotation_config", "existing_umap_pdf", "aggregated_module_tables"),
  required_for = c("Seurat rebuild and per-cell module scores", "UMAP atlas and score projection",
                   "audited labels", "visual reference only", "fallback donor/cell-type summaries"),
  status = c(
    ifelse(any(grepl("\\.(rds|RDS|rda|RData|h5ad|h5seurat|loom|mtx|mtx.gz|h5)$", hits, ignore.case = TRUE)), "candidate_found", "missing"),
    ifelse(length(embedding_hits) > 0, "candidate_found", "missing"),
    ifelse(file.exists("config/gse231569_cluster_to_celltype_audited.tsv"), "available", "missing"),
    ifelse(file.exists("results/scrna/gse231569/figures/gse231569_umap_audited_broad_celltype.pdf"), "available_pdf_only", "missing"),
    ifelse(file.exists("results/tables/magma_scrna_module_score_by_donor_celltype.tsv"), "available", "missing")
  ),
  path_or_note = c(
    paste(hits, collapse = "; "),
    paste(embedding_hits, collapse = "; "),
    "config/gse231569_cluster_to_celltype_audited.tsv",
    "results/scrna/gse231569/figures/gse231569_umap_audited_broad_celltype.pdf",
    "results/tables/magma_scrna_module_score_by_donor_celltype.tsv"
  )
)
audit[path_or_note == "", path_or_note := "not found in current workspace"]
fwrite(audit, "results/tables/gse231569_celllevel_resource_audit_v0.1.tsv", sep = "\t")

preferred_object <- "results/scrna/gse231569/objects/gse231569_annotated_audited.rds"
can_rebuild <- file.exists(preferred_object)

if (!can_rebuild) {
  writeLines(c(
    "# GSE231569 Cell-level Resource Blocker Memo v0.1",
    "",
    "A rebuild of a GSE231569 audited Seurat object was not attempted because the current workspace does not contain both a usable cell-level expression object/matrix and a cell-level UMAP or embedding table.",
    "",
    "Available resources support aggregated audited cell-type and donor-celltype module summaries, but not a compliant per-cell MAGMA score UMAP.",
    "",
    "Required files for a future rebuild:",
    "",
    "- Cell-by-gene count or normalized expression matrix, with cell barcodes.",
    "- Gene feature table.",
    "- Cell-level metadata with donor/sample labels.",
    "- Cell-level UMAP coordinates or enough raw data to recompute UMAP.",
    "- Audited cell-type labels or cluster-to-celltype mapping.",
    "",
    "Manuscript boundary: main figures should not display text implying that a per-cell UMAP score projection was performed."
  ), "docs/gse231569_celllevel_resource_blocker_memo.md", useBytes = TRUE)
  message("GSE231569 cell-level rebuild blocked by missing resources")
} else {
  obj <- readRDS(preferred_object)
  saveRDS(obj, "data/processed/gse231569_audited_seurat.rds")
  writeLines(c(
    "# GSE231569 Cell-level Resource Memo v0.1",
    "",
    "A usable audited GSE231569 Seurat object was found in the current workspace and copied to `data/processed/gse231569_audited_seurat.rds` for Phase 9 cell-level score projection.",
    "",
    "This object contains cell-level metadata and a UMAP reduction. It can support audited cell-type UMAP and MAGMA module-score projection."
  ), "docs/gse231569_celllevel_resource_blocker_memo.md", useBytes = TRUE)
  message("wrote data/processed/gse231569_audited_seurat.rds")
}
