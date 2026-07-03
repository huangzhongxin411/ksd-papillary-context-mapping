suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

design <- fread(file.path(table_dir, "gse73680_analysis_design.tsv"))
p1 <- fread(file.path(table_dir, "gse73680_p1_gene_response.tsv"))
modules <- fread(file.path(table_dir, "gse73680_module_score_response.tsv"))
p1_patient <- fread(file.path(table_dir, "gse73680_patient_level_p1_gene_response.tsv"))
module_patient <- fread(file.path(table_dir, "gse73680_patient_level_module_response.tsv"))
inj <- fread(file.path(table_dir, "gse73680_module_injury_correlation.tsv"))
bench <- fread(file.path(table_dir, "gse73680_random_module_benchmark.tsv"))

get_design <- function(item) design[design_item == item, value][1]
p1_fdr <- p1[fdr < 0.05]
p1_nom <- p1[p_value < 0.05 & fdr >= 0.05]
module_fdr <- modules[fdr < 0.05]
module_nom <- modules[p_value < 0.05 & fdr >= 0.05]
p1_patient_fdr <- p1_patient[fdr < 0.05]
p1_patient_nom <- p1_patient[p_value < 0.05 & fdr >= 0.05]
module_patient_fdr <- module_patient[fdr < 0.05]
module_patient_nom <- module_patient[p_value < 0.05 & fdr >= 0.05]
inj_strong <- inj[interpretation == "strong disease-context association"]
bench_hit <- bench[benchmark_interpretation == "module response exceeds random expectation"]

summary <- data.table(
  evidence_item = c("analysis_design", "P1_gene_response", "module_score_response",
                    "patient_level_P1_sensitivity", "patient_level_module_sensitivity",
                    "module_injury_correlation", "random_gene_set_benchmark", "overall_result4_readiness"),
  result = c(
    paste0(get_design("recommended_model"), "; ", get_design("n_samples_included"), " included samples; ",
           get_design("n_patients_with_both_groups"), " paired patients"),
    if (nrow(p1_fdr)) paste(p1_fdr$gene, collapse = ";") else if (nrow(p1_nom)) paste(p1_nom$gene, collapse = ";") else "no P1 gene reached nominal P < 0.05",
    if (nrow(module_fdr)) paste(module_fdr$module_name, collapse = ";") else if (nrow(module_nom)) paste(module_nom$module_name, collapse = ";") else "no module reached nominal P < 0.05",
    if (nrow(p1_patient_fdr)) paste(p1_patient_fdr$gene, collapse = ";") else if (nrow(p1_patient_nom)) paste(p1_patient_nom$gene, collapse = ";") else "no patient-level P1 gene reached nominal P < 0.05",
    if (nrow(module_patient_fdr)) paste(module_patient_fdr$module_name, collapse = ";") else if (nrow(module_patient_nom)) paste(module_patient_nom$module_name, collapse = ";") else "no patient-level module reached nominal P < 0.05",
    if (nrow(inj_strong)) paste(unique(inj_strong$module_name), collapse = ";") else "no strong module-injury association by FDR",
    if (nrow(bench_hit)) paste(bench_hit$module_name, collapse = ";") else "no tested module exceeded random expectation",
    "GSE73680 disease-context analysis completed with patient-aware design"
  ),
  strength = c(
    "moderate",
    if (nrow(p1_fdr)) "strong" else if (nrow(p1_nom)) "exploratory" else "negative_or_no_signal",
    if (nrow(module_fdr)) "strong" else if (nrow(module_nom)) "exploratory" else "negative_or_no_signal",
    if (nrow(p1_patient_fdr)) "strong" else if (nrow(p1_patient_nom)) "exploratory" else "negative_or_no_signal",
    if (nrow(module_patient_fdr)) "strong" else if (nrow(module_patient_nom)) "exploratory" else "negative_or_no_signal",
    if (nrow(inj_strong)) "strong" else "negative_or_no_signal",
    if (nrow(bench_hit)) "strong" else "negative_or_no_signal",
    if (nrow(p1_fdr) || nrow(module_fdr) || nrow(bench_hit)) "moderate" else "exploratory"
  ),
  interpretation = c(
    "Repeated patient structure requires patient-aware modeling for disease-context comparisons.",
    if (nrow(p1_fdr)) "At least one P1 gene shows FDR-supported disease-context expression association." else if (nrow(p1_nom)) "Selected P1 genes show nominal exploratory disease-context associations." else "P1 genes show no detectable disease-context differential expression under this design.",
    if (nrow(module_fdr)) "At least one module shows FDR-supported disease-context score association." else if (nrow(module_nom)) "Selected modules show nominal exploratory disease-context associations." else "Module scores show no detectable disease-context difference under this design.",
    if (nrow(p1_patient_fdr)) "Patient-level paired analysis supports at least one P1 gene." else if (nrow(p1_patient_nom)) "Patient-level paired analysis shows nominal P1 signal." else "Patient-level paired P1 analysis does not support a robust single-gene response.",
    if (nrow(module_patient_fdr)) "Patient-level paired analysis supports MAGMA/injury module disease-context shifts." else if (nrow(module_patient_nom)) "Patient-level paired analysis shows nominal module signal." else "Patient-level paired module analysis does not support robust module shifts.",
    if (nrow(inj_strong)) "Some modules correlate strongly with injury-remodeling context." else "Modules do not show strong FDR-supported correlation with injury remodeling.",
    if (nrow(bench_hit)) "At least one module response exceeds random gene-set expectation." else "Observed module responses are not stronger than matched random gene sets.",
    "Use as Result 4 only if the combined evidence is strong or moderate; otherwise place in Supplementary or report as a bounded negative/exploratory analysis."
  ),
  claim_boundary = c(
    "Patient structure is filename-derived and must not be treated as randomized independent sampling.",
    "Expression association only; no causal plaque formation claim.",
    "Module score association only; no genetic mediation or mechanism proof.",
    "Patient-level sensitivity is supportive but remains expression association only.",
    "Patient-level sensitivity is supportive but remains expression association only.",
    "Correlation only; no causal injury-remodeling directionality.",
    "Benchmark is empirical and expression-background limited.",
    "No causal validation, spatial validation, TWAS/coloc support or TAL-mediated mechanism proof."
  ),
  use_in_main_text = c("yes_methods", ifelse(nrow(p1_fdr), "yes", ifelse(nrow(p1_nom), "maybe_exploratory", "no_or_negative")),
                       ifelse(nrow(module_fdr), "yes", ifelse(nrow(module_nom), "maybe_exploratory", "no_or_negative")),
                       ifelse(nrow(p1_patient_fdr), "yes_sensitivity", ifelse(nrow(p1_patient_nom), "maybe_exploratory", "no_or_negative")),
                       ifelse(nrow(module_patient_fdr), "yes_sensitivity", ifelse(nrow(module_patient_nom), "maybe_exploratory", "no_or_negative")),
                       ifelse(nrow(inj_strong), "yes", "supplement_or_negative"),
                       ifelse(nrow(bench_hit), "yes", "supplement_or_negative"), "conditional"),
  use_in_figure4 = c("panel_A_design", ifelse(nrow(p1_fdr) || nrow(p1_nom), "candidate_panel_B", "supplement"),
                     ifelse(nrow(module_fdr) || nrow(module_nom), "candidate_panel_C", "supplement"),
                     ifelse(nrow(p1_patient_fdr) || nrow(p1_patient_nom), "supplement_or_panel_B_annotation", "supplement"),
                     ifelse(nrow(module_patient_fdr) || nrow(module_patient_nom), "candidate_panel_C_annotation", "supplement"),
                     ifelse(nrow(inj_strong), "candidate_panel_D", "supplement"),
                     ifelse(nrow(bench_hit), "candidate_panel_D", "supplement"),
                     "decide_after_review")
)
fwrite(summary, file.path(table_dir, "gse73680_disease_context_summary.tsv"), sep = "\t")

writeLines(c(
  "# GSE73680 Disease-Context Results Memo",
  "",
  "## Status",
  "",
  "Phase 5B disease-context analysis completed. The analysis uses patient-aware modeling where repeated patient samples are present and keeps claims limited to expression associations.",
  "",
  "## Design",
  "",
  paste0("- Recommended model: ", get_design("recommended_model")),
  paste0("- Included samples: ", get_design("n_samples_included")),
  paste0("- Included patients: ", get_design("n_patients_included")),
  paste0("- Patients with both groups: ", get_design("n_patients_with_both_groups")),
  "",
  "## Evidence Summary",
  "",
  paste0("- P1 gene response: ", summary[evidence_item == "P1_gene_response", result], " (", summary[evidence_item == "P1_gene_response", strength], ")"),
  paste0("- Module response: ", summary[evidence_item == "module_score_response", result], " (", summary[evidence_item == "module_score_response", strength], ")"),
  paste0("- Patient-level P1 sensitivity: ", summary[evidence_item == "patient_level_P1_sensitivity", result], " (", summary[evidence_item == "patient_level_P1_sensitivity", strength], ")"),
  paste0("- Patient-level module sensitivity: ", summary[evidence_item == "patient_level_module_sensitivity", result], " (", summary[evidence_item == "patient_level_module_sensitivity", strength], ")"),
  paste0("- Injury correlation: ", summary[evidence_item == "module_injury_correlation", result], " (", summary[evidence_item == "module_injury_correlation", strength], ")"),
  paste0("- Random benchmark: ", summary[evidence_item == "random_gene_set_benchmark", result], " (", summary[evidence_item == "random_gene_set_benchmark", strength], ")"),
  "",
  "## Claim Boundary",
  "",
  "Allowed wording: GSE73680 disease-context analysis reports expression associations for P1 genes or MAGMA/scRNA-linked modules.",
  "",
  "Not allowed wording: causal validation, spatial validation, TWAS/coloc support, or proof of a TAL-mediated mechanism.",
  "",
  "## Figure 4 Decision",
  "",
  "Use Figure 4 in the main text only if the combined P1/module/random-benchmark evidence is strong or coherent enough after review. Otherwise, report GSE73680 as supplementary disease-context evidence or as a bounded negative/exploratory analysis."
), "docs/gse73680_disease_context_results_memo.md")
message("wrote GSE73680 disease-context summary")
