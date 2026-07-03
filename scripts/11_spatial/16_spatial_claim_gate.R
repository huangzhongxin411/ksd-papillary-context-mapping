suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

out_dir <- "results/spatial/phase28c_projection"
doc_dir <- "docs/spatial"
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

load_qc <- read_tsv("results/spatial/phase28b_loading_rescue/spatial_load_qc_v0.2.tsv", show_col_types = FALSE)
gate_prev <- read_tsv("results/spatial/phase28b_loading_rescue/spatial_projection_gate_v0.1.tsv", show_col_types = FALSE)
overlap <- read_tsv(file.path(out_dir, "spatial_signature_gene_overlap_v0.1.tsv"), show_col_types = FALSE)
score_qc <- read_tsv(file.path(out_dir, "spatial_module_score_qc_v0.1.tsv"), show_col_types = FALSE)
summary <- read_tsv(file.path(out_dir, "spatial_projection_summary_v0.1.tsv"), show_col_types = FALSE)
p1 <- read_tsv(file.path(out_dir, "p1_gene_spatial_expression_audit_v0.1.tsv"), show_col_types = FALSE)

all_loaded <- nrow(load_qc) >= 5 && all(load_qc$qc_status == "pass")
risk_modules <- c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")
risk_projected <- score_qc %>%
  filter(signature %in% risk_modules, score_allowed) %>%
  count(signature) %>%
  filter(n >= 4) %>%
  pull(signature)

loop_support <- summary %>%
  filter(risk_signature %in% risk_modules, context_signature == "Loop_TAL",
         support_class %in% c("strong_spatial_context_support", "moderate_spatial_context_support"))
injury_support <- summary %>%
  filter(risk_signature %in% risk_modules, context_signature == "Injury_remodeling",
         support_class %in% c("strong_spatial_context_support", "moderate_spatial_context_support"))
p1_resolved <- p1 %>%
  group_by(gene) %>%
  summarise(n_samples_resolved = sum(spatial_plot_allowed), .groups = "drop") %>%
  filter(n_samples_resolved > 0)

claim_rows <- tribble(
  ~claim, ~allowed, ~evidence_required, ~evidence_observed, ~recommended_wording, ~forbidden_wording,
  "GSE206306 spatial resources were loaded and QC-passed.",
  all_loaded,
  "Five samples with matrix, image, coordinates and QC pass.",
  paste0(sum(load_qc$qc_status == "pass"), "/", nrow(load_qc), " samples QC pass."),
  "GSE206306 spatial objects were loaded and passed technical QC.",
  "Spatial validation of genetic risk.",

  "MAGMA risk modules can be projected onto papillary spatial sections.",
  length(risk_projected) > 0,
  "Risk module scores allowed in at least four QC-passed samples.",
  paste(risk_projected, collapse = "; "),
  "MAGMA-prioritized modules were projected as spatial context scores in QC-passed papillary sections.",
  "Risk genes spatially validated in papilla.",

  "MAGMA risk modules show spatial co-distribution with Loop/TAL context.",
  nrow(loop_support) > 0,
  "At least one MAGMA module has moderate/strong adjusted cross-sample support with Loop/TAL.",
  paste(loop_support$risk_signature, loop_support$support_class, sep = ":", collapse = "; "),
  "MAGMA module projections showed cross-sample co-distribution with Loop/TAL context signatures.",
  "Loop/TAL spatial validation of genetic risk.",

  "MAGMA risk modules show spatial co-distribution with injury/remodeling context.",
  nrow(injury_support) > 0,
  "At least one MAGMA module has moderate/strong adjusted cross-sample support with injury/remodeling.",
  paste(injury_support$risk_signature, injury_support$support_class, sep = ":", collapse = "; "),
  "MAGMA module projections showed cross-sample co-distribution with injury/remodeling signatures.",
  "Causal injury niche or plaque-specific localization.",

  "P1 genes show spatially resolvable expression.",
  nrow(p1_resolved) > 0,
  "Gene-by-gene detection in at least one sample with >=5% detected spots.",
  paste(p1_resolved$gene, p1_resolved$n_samples_resolved, sep = ":", collapse = "; "),
  "P1 genes were assessed gene-by-gene; only genes passing detection thresholds were visualized.",
  "Spatially validated P1 genes.",

  "Spatial validation of genetic risk.",
  FALSE,
  "Would require independent spatial validation design and stronger causal/localization evidence.",
  "Not observed; Phase 28C is spatial projection only.",
  "Do not use.",
  "Spatial validation of genetic risk.",

  "Plaque-specific localization.",
  FALSE,
  "Would require plaque ROI annotations and lesion-specific tests.",
  "No plaque ROI annotations used.",
  "Do not use.",
  "Plaque-specific localization.",

  "Causal spatial niche.",
  FALSE,
  "Would require causal or perturbational evidence.",
  "No causal evidence generated.",
  "Do not use.",
  "Causal spatial niche."
)

write_tsv(claim_rows, file.path(out_dir, "spatial_claim_gate_v0.1.tsv"))

passed <- overlap %>%
  filter(score_allowed) %>%
  group_by(signature) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  arrange(signature)

scores <- score_qc %>%
  filter(score_allowed) %>%
  group_by(signature) %>%
  summarise(n_samples_scored = n_distinct(sample_id), median_genes_used = median(n_genes_used), .groups = "drop") %>%
  arrange(signature)

strongest <- summary %>%
  filter(is.finite(median_rho_adjusted)) %>%
  arrange(desc(median_rho_adjusted)) %>%
  slice_head(n = 10)

fig_files <- list.files("results/figures/supplementary/spatial_projection", pattern = "\\.(pdf|png)$", full.names = FALSE)
safe_wording <- if (nrow(loop_support) > 0 || nrow(injury_support) > 0) {
  "GSE206306 spatial projection provided additional papillary context support: MAGMA-prioritized modules were spatially projected in QC-passed human renal papillary sections and showed consistent co-distribution with Loop/TAL and/or injury-remodeling signatures. These analyses were interpreted as spatial context mapping rather than plaque-specific validation or causal mediation."
} else {
  "Although GSE206306 spatial objects were loaded successfully, signature projection did not show consistent cross-sample spatial co-distribution; spatial data were therefore retained as a resource and QC extension rather than a support layer."
}

memo <- c(
  "# Phase 28C spatial projection completion memo v0.1",
  "",
  "## 1. Signatures passing overlap thresholds",
  paste0("- ", passed$signature, ": ", passed$n_samples, " samples", collapse = "\n"),
  "",
  "## 2. Scores calculated",
  paste0("- ", scores$signature, ": ", scores$n_samples_scored, " samples; median genes used=", scores$median_genes_used, collapse = "\n"),
  "",
  "## 3. Samples used",
  paste0("- ", load_qc$sample_id, ": ", load_qc$n_spots_matrix, " spots; QC=", load_qc$qc_status, collapse = "\n"),
  "",
  "## 4. Spatial overlays generated",
  paste0("- ", fig_files, collapse = "\n"),
  "",
  "## 5. Strongest adjusted correlations",
  paste0("- ", strongest$risk_signature, " vs ", strongest$context_signature,
         ": median adjusted rho=", round(strongest$median_rho_adjusted, 3),
         "; support=", strongest$support_class, collapse = "\n"),
  "",
  "## 6. Consistency across samples",
  paste0("- Loop/TAL supported pairs: ", ifelse(nrow(loop_support) > 0, paste(loop_support$risk_signature, collapse = ", "), "none")),
  paste0("- Injury/remodeling supported pairs: ", ifelse(nrow(injury_support) > 0, paste(injury_support$risk_signature, collapse = ", "), "none")),
  "",
  "## 7. P1 resolvability",
  if (nrow(p1_resolved) > 0) paste0("- ", p1_resolved$gene, ": resolved in ", p1_resolved$n_samples_resolved, " samples", collapse = "\n") else "- No P1 gene passed the >=5% detected-spot plotting threshold.",
  "",
  "## 8. Supplementary Figure S7 readiness",
  "- Supplementary Figure S7 candidate was generated if all required plotting inputs were available.",
  "",
  "## 9. Should manuscript change?",
  "- Do not change the main conclusion or Figures 1-5 before review. This phase can be considered as a supplementary spatial context mapping layer only after claim-gate review.",
  "",
  "## 10. Exact claim-safe wording",
  safe_wording,
  "",
  "## Forbidden wording",
  "- spatial validation",
  "- plaque-specific validation",
  "- genetic risk localizes to plaque",
  "- causal papillary niche",
  "- spatially validated P1 genes"
)

writeLines(memo, file.path(doc_dir, "phase28c_spatial_projection_completion_memo_v0.1.md"))
