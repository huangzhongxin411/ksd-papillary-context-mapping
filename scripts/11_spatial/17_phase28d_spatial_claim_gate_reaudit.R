suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

in_dir <- "results/spatial/phase28c_projection"
out_dir <- "results/spatial/phase28d_spatial_hardening"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary <- read_tsv(file.path(in_dir, "spatial_projection_summary_v0.1.tsv"), show_col_types = FALSE)
gate <- read_tsv(file.path(in_dir, "spatial_claim_gate_v0.1.tsv"), show_col_types = FALSE)
p1 <- read_tsv(file.path(in_dir, "p1_gene_spatial_expression_audit_v0.1.tsv"), show_col_types = FALSE)

supported_pair <- function(risk, context) {
  x <- summary %>%
    filter(risk_signature %in% risk, context_signature == context,
           support_class %in% c("moderate_spatial_context_support", "strong_spatial_context_support"))
  if (!nrow(x)) return("none")
  paste0(x$risk_signature, " median adjusted rho=", round(x$median_rho_adjusted, 3),
         " (", x$support_class, ")", collapse = "; ")
}

p1_resolved <- p1 %>%
  group_by(gene) %>%
  summarise(n_samples_resolved = sum(spatial_plot_allowed), .groups = "drop") %>%
  arrange(match(gene, c("CLDN10", "HIBADH", "PKD2", "UMOD", "CASR", "CLDN14")))

out <- tribble(
  ~claim, ~status, ~supporting_result, ~allowed_wording, ~forbidden_wording, ~manuscript_use,
  "GSE206306 spatial resources were recovered, loaded and QC-passed.",
  "allowed",
  "Phase 28B/28C: 5/5 samples loaded with image and coordinates; QC pass.",
  "GSE206306 spatial resources were recovered, loaded and passed technical QC.",
  "spatial validation of genetic risk",
  "Methods and supplementary Results",

  "MAGMA risk modules can be projected onto human renal papillary spatial sections.",
  "allowed",
  "MAGMA_top50, MAGMA_top100, MAGMA_FDR and MAGMA_suggestive were scoreable in 5/5 samples.",
  "MAGMA-prioritized modules were projected onto QC-passed human renal papillary spatial sections.",
  "risk genes were spatially validated",
  "Supplementary Results only",

  "MAGMA_top50/top100 showed moderate spatial co-distribution with Loop/TAL context.",
  "allowed",
  supported_pair(c("MAGMA_top50", "MAGMA_top100"), "Loop_TAL"),
  "MAGMA_top50 and MAGMA_top100 showed moderate cross-sample co-distribution with Loop/TAL context signatures.",
  "Loop/TAL spatial validation of genetic risk",
  "Supplementary Results with boundary sentence",

  "MAGMA_top50/top100 showed moderate spatial co-distribution with ion/mineral-handling context.",
  "allowed",
  supported_pair(c("MAGMA_top50", "MAGMA_top100"), "Ion_mineral_handling"),
  "MAGMA_top50 and MAGMA_top100 showed moderate cross-sample co-distribution with ion/mineral-handling signatures.",
  "ion-handling mechanism or pathway activity",
  "Supplementary Results with boundary sentence",

  "Injury/remodeling spatial support.",
  "not_allowed",
  supported_pair(c("MAGMA_top50", "MAGMA_top100"), "Injury_remodeling"),
  "Injury/remodeling signatures did not meet the predefined support threshold.",
  "injury/remodeling spatial support",
  "Limitations / boundary only",

  "P1 spatial validation.",
  "not_allowed",
  paste0(p1_resolved$gene, " resolved in ", p1_resolved$n_samples_resolved, "/5 samples", collapse = "; "),
  "P1 genes were assessed only for spatial assay resolvability.",
  "spatially validated P1 genes; P1 disease validation",
  "Supplementary Results boundary only",

  "Plaque-specific localization.",
  "not_allowed",
  "No plaque ROI annotations or lesion-specific spatial tests were used.",
  "No plaque-specific localization claim is made.",
  "plaque-specific localization; genetic risk localizes to plaque",
  "Limitations only",

  "Causal mediation.",
  "not_allowed",
  "No causal or perturbational spatial evidence was generated.",
  "No causal spatial inference is made.",
  "causal papillary niche; causal mediation",
  "Limitations only",

  "Spatial validation of genetic risk.",
  "not_allowed",
  "Phase 28D is claim-gated spatial context mapping, not validation.",
  "Spatial projection was interpreted as context mapping only.",
  "spatial validation",
  "Limitations and claim-boundary audit"
)

write_tsv(out, file.path(out_dir, "spatial_claim_gate_final_v0.1.tsv"))
