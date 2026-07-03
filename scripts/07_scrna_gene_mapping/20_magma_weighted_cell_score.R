suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  tal = "#3E6672",
  scrna = "#6F8F98",
  p1 = "#B59A5B",
  muted = "#D8D8D8",
  text = "#303030",
  grid = "#A8A8A8"
)

label_cell <- function(x) {
  map <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

label_module <- function(x) {
  map <- c(
    magma_top50_score = "MAGMA top 50",
    magma_top100_score = "MAGMA top 100",
    magma_fdr05_score = "MAGMA FDR",
    magma_suggestive_p1e4_score = "MAGMA suggestive"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

object_candidates <- c(
  "data/processed/GSE231569/GSE231569_seurat.rds",
  "data/processed/gse231569/gse231569_seurat.rds",
  "data/raw/GSE231569/GSE231569_seurat.rds",
  "results/scrna/gse231569/gse231569_seurat.rds"
)

embedding_candidates <- list.files(
  path = c("data", "results", "config"),
  pattern = "umap|embedding|coord|cell_metadata",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
embedding_candidates <- embedding_candidates[grepl("\\.(tsv|csv|txt|rds|RDS)$", embedding_candidates)]

readiness <- data.table(
  item = c("GSE231569 Seurat object", "cell-level UMAP or embedding table", "audited donor-celltype module table"),
  required_for = c("per-cell MAGMA-weighted module score", "UMAP score projection", "donor-level summary"),
  status = c(
    ifelse(any(file.exists(object_candidates)), "available", "missing"),
    ifelse(length(embedding_candidates) > 0, "candidate_found", "missing"),
    ifelse(file.exists("results/tables/magma_scrna_module_score_by_donor_celltype.tsv"), "available", "missing")
  ),
  path_or_note = c(
    paste(object_candidates[file.exists(object_candidates)], collapse = "; "),
    paste(embedding_candidates, collapse = "; "),
    "results/tables/magma_scrna_module_score_by_donor_celltype.tsv"
  )
)
readiness[path_or_note == "", path_or_note := "not found in current workspace"]
fwrite(readiness, "results/tables/gse231569_magma_score_projection_readiness_v0.1.tsv", sep = "\t")

donor_scores <- fread("results/tables/magma_scrna_module_score_by_donor_celltype.tsv")
score_cols <- c("magma_top50_score", "magma_top100_score", "magma_fdr05_score", "magma_suggestive_p1e4_score")
score_long <- melt(
  donor_scores,
  id.vars = c("donor_id", "sample_id", "disease_status", "audited_broad_cell_type", "n_cells", "annotation_confidence"),
  measure.vars = score_cols,
  variable.name = "module",
  value.name = "module_score"
)
score_long[, module_label := factor(label_module(module),
                                    levels = c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive"))]
score_long[, cell_label := factor(label_cell(audited_broad_cell_type))]
score_long[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
score_long[, interpretation_scope := "donor-celltype module score summary; not a per-cell scDRS analysis"]
fwrite(score_long, "results/tables/gse231569_magma_score_by_celltype_donor.tsv", sep = "\t")

celltype_summary <- score_long[, .(
  n_donors = uniqueN(donor_id),
  median_score = as.numeric(median(module_score, na.rm = TRUE)),
  mean_score = as.numeric(mean(module_score, na.rm = TRUE)),
  min_score = as.numeric(min(module_score, na.rm = TRUE)),
  max_score = as.numeric(max(module_score, na.rm = TRUE)),
  median_n_cells = as.numeric(median(n_cells, na.rm = TRUE))
), by = .(audited_broad_cell_type, cell_label, module_label)]
setorder(celltype_summary, module_label, -median_score)
fwrite(celltype_summary, "results/tables/gse231569_magma_score_by_celltype_summary.tsv", sep = "\t")

fig_violin <- ggplot(score_long, aes(cell_label, module_score, fill = is_tal)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "#555555", linewidth = 0.25) +
  geom_point(aes(size = n_cells), position = position_jitter(width = 0.10, height = 0),
             shape = 21, color = "#555555", stroke = 0.2, alpha = 0.85) +
  facet_wrap(~ module_label, ncol = 2) +
  scale_fill_manual(values = c("TRUE" = pal$tal, "FALSE" = pal$muted), labels = c("Other", "Loop/TAL")) +
  scale_size_continuous(range = c(1.4, 4.2)) +
  labs(
    title = "Donor-level MAGMA module scores across audited GSE231569 cell types",
    subtitle = "Cell-level UMAP score projection requires a Seurat object or embedding table not present in the current workspace",
    x = NULL,
    y = "Donor-celltype module score",
    fill = NULL,
    size = "Cells"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "#555555"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/figure2_magma_score_violin_v0.1.pdf", fig_violin,
       width = 10.8, height = 7.2, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure2_magma_score_violin_v0.1.png", fig_violin,
       width = 10.8, height = 7.2, units = "in", dpi = 260, bg = "white")

fig_readiness <- ggplot(readiness, aes(status, factor(item, levels = rev(item)), fill = status)) +
  geom_col(width = 0.58, color = "#666666", linewidth = 0.25) +
  geom_text(aes(label = status), hjust = 1.08, color = "white", size = 3) +
  scale_fill_manual(values = c(available = pal$tal, candidate_found = pal$p1, missing = "#9A5F52")) +
  labs(
    title = "GSE231569 MAGMA score UMAP readiness",
    subtitle = "Per-cell score UMAP was not generated because required cell-level resources are missing",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank())

ggsave("results/figures/figure2_magma_score_umap_v0.1.pdf", fig_readiness,
       width = 8.2, height = 3.6, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure2_magma_score_umap_v0.1.png", fig_readiness,
       width = 8.2, height = 3.6, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# GSE231569 MAGMA Score Projection Resource Note v0.1",
  "",
  "A formal per-cell MAGMA-weighted KSD module score and UMAP projection was not generated in this run because no usable GSE231569 Seurat object, h5ad object, or cell-level UMAP coordinate table was found in the current workspace.",
  "",
  "Generated outputs therefore use the already audited donor-celltype MAGMA module score table and should be described as donor-celltype module score summaries, not as scDRS or per-cell disease relevance scores.",
  "",
  "Required resources for a future per-cell projection:",
  "",
  "- Cell-by-gene normalized expression matrix or Seurat object.",
  "- Cell-level UMAP coordinates.",
  "- Audited cell type labels.",
  "- Donor/sample metadata.",
  "- MAGMA gene-level weights or explicitly unweighted module-score definition."
), "docs/gse231569_magma_score_umap_resource_note_v0.1.md")

message("wrote GSE231569 MAGMA score readiness and donor-celltype summaries")
