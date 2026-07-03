suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results/plaque", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

p1_resp <- fread("results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv")
module_resp <- fread("results/gse73680/tables/gse73680_patient_level_module_response.tsv")
coupling <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
fwrite(p1_resp, "results/plaque/gse73680_candidate_gene_response.tsv", sep = "\t")

twas_smr_placeholder <- data.table(
  gene = character(),
  source = character(),
  plaque_response_status = character(),
  note = "TWAS/SMR-coloc-supported genes are pending resource-confirmed analyses"
)
fwrite(twas_smr_placeholder, "results/plaque/gse73680_twas_smr_gene_response.tsv", sep = "\t")
fwrite(coupling, "results/plaque/gse73680_risk_injury_coupling.tsv", sep = "\t")

integrated <- module_resp[, .(module_name, n_paired_patients, paired_delta, p_value, fdr, interpretation)]
integrated <- merge(integrated, coupling[analysis == "Paired delta" & injury_module == "injury_remodeling",
                                        .(module_name, injury_delta_rho = rho, injury_delta_fdr = fdr, robustness_summary)],
                    by = "module_name", all.x = TRUE)
integrated[, plaque_context_interpretation := fcase(
  fdr < 0.05 & injury_delta_rho > 0 & injury_delta_fdr < 0.10, "module plaque response with paired injury-coupling support",
  fdr < 0.05, "module plaque response support",
  !is.na(injury_delta_rho) & injury_delta_rho > 0 & injury_delta_fdr < 0.10, "injury-coupling support only",
  default = "limited plaque context support"
)]
fwrite(integrated, "results/plaque/gse73680_plaque_context_integrated.tsv", sep = "\t")
writeLines(c(
  "# GSE73680 Integrated Plaque Context v0.1",
  "",
  "This table integrates existing module plaque response and risk-injury coupling evidence.",
  "It supports plaque/stone papilla disease-context association only and should not be described as disease-gene validation."
), "docs/gse73680_integrated_plaque_context_v0.1.md", useBytes = TRUE)
message("wrote integrated plaque context outputs")
