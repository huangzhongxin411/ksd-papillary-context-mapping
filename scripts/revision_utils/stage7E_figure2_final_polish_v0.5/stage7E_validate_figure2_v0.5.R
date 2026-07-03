#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c(
  "results/figures/revision/stage7E_figure2_final_polish_v0.5/figure2_snRNA_context_draft_v0.5.pdf",
  "results/figures/revision/stage7E_figure2_final_polish_v0.5/figure2_snRNA_context_draft_v0.5.png",
  "results/figures/revision/stage7E_figure2_final_polish_v0.5/figure2_snRNA_context_draft_v0.5.svg",
  "docs/revision/stage7E_figure2_final_polish_v0.5/figure2_qc_audit_v0.5.md",
  "results/tables/revision/stage7E_figure2_final_polish_v0.5/figure2_visual_polish_log_v0.5.tsv",
  "results/tables/revision/stage7E_figure2_final_polish_v0.5/figure2_source_data_manifest_v0.5.tsv",
  "docs/revision/stage7E_figure2_final_polish_v0.5/figure2_legend_v0.5.md",
  "docs/revision/stage7E_figure2_final_polish_v0.5/stage7E_figure2_final_polish_report.md"
)

checks <- character()
check_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  checks <<- c(checks, paste0("PASS\t", label))
}

check_true(all(file.exists(file.path(root, required))), "all_8_required_outputs_exist")
check_true(all(file.info(file.path(root, required))$size > 0), "all_required_outputs_nonempty")

manifest <- read.delim(file.path(root, required[6]), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
polish <- read.delim(file.path(root, required[5]), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
manifest_headers <- c(
  "panel", "source_data_file", "source_exists", "row_count_if_readable",
  "column_count_if_readable", "data_changed_from_v0.4", "claim_changed_from_v0.4",
  "ready_for_publication_source_data", "notes"
)
polish_headers <- c(
  "panel", "issue_in_v0.4", "polish_made_in_v0.5", "data_changed",
  "claim_changed", "reason", "human_review_needed", "notes"
)
check_true(identical(names(manifest), manifest_headers), "source_manifest_schema_exact")
check_true(identical(names(polish), polish_headers), "visual_polish_log_schema_exact")
check_true(identical(manifest$panel, LETTERS[1:7]), "source_manifest_panels_A_to_G")
check_true(identical(polish$panel, LETTERS[1:7]), "visual_polish_log_panels_A_to_G")
check_true(all(manifest$data_changed_from_v0.4 == "no"), "manifest_data_changed_all_no")
check_true(all(manifest$claim_changed_from_v0.4 == "no"), "manifest_claim_changed_all_no")
check_true(all(polish$data_changed == "no"), "polish_data_changed_all_no")
check_true(all(polish$claim_changed == "no"), "polish_claim_changed_all_no")
check_true(all(polish$human_review_needed == "yes"), "human_review_flag_all_yes")

expected_rows <- c(6L, 92L, 16L, 16L, 4L, 20L, 7L)
expected_cols <- c(6L, 17L, 14L, 12L, 19L, 17L, 9L)
check_true(identical(as.integer(manifest$row_count_if_readable), expected_rows), "manifest_row_counts_match_frozen_sources")
check_true(identical(as.integer(manifest$column_count_if_readable), expected_cols), "manifest_column_counts_match_frozen_sources")
check_true(all(file.exists(file.path(root, manifest$source_data_file))), "all_manifest_source_paths_exist")
for (i in seq_len(nrow(manifest))) {
  source <- read.delim(file.path(root, manifest$source_data_file[i]), sep = "\t", check.names = FALSE)
  check_true(nrow(source) == expected_rows[i] && ncol(source) == expected_cols[i], paste0("panel_", manifest$panel[i], "_source_dimensions_verified"))
}

svg <- paste(readLines(file.path(root, required[3]), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
required_text <- c(
  "Donor-level snRNA context mapping of MAGMA-prioritized modules",
  "Moderate Loop/TAL-associated context with partial matched-random support",
  "Resource overview", "Donor x compartment module scores", "Loop/TAL rank",
  "Single-donor exclusion", "Matched-random benchmark", "Driver-panel removal",
  "Interpretation boundary", "Partial matched-random support", "Top-ranked in 4/4 donors",
  "Matched-random percentile", "Overall: partial support", "Not causal", "Not plaque"
)
for (needle in required_text) {
  check_true(grepl(needle, svg, fixed = TRUE), paste0("svg_contains_", gsub("[^A-Za-z0-9]+", "_", needle)))
}
forbidden_positive <- c(
  "Strong enrichment", "Validated", "Independent validation", "Therapeutic target",
  "high-confidence causal", "plaque-specific localization", "R1_MAGMA", "MAGMA_top"
)
for (needle in forbidden_positive) {
  check_true(!grepl(needle, svg, fixed = TRUE), paste0("svg_excludes_", gsub("[^A-Za-z0-9]+", "_", needle)))
}

qc_text <- readLines(file.path(root, required[4]), warn = FALSE)
check_true(sum(grepl("\\| PASS \\|", qc_text)) == 15L, "all_15_requested_qc_items_pass")
check_true(any(grepl("source-backed percentile summary", qc_text, fixed = TRUE)), "panel_E_summary_only_limitation_documented")

figure_dir <- file.path(root, "results/figures/revision/stage7E_figure2_final_polish_v0.5")
check_true(length(list.files(figure_dir, pattern = "figure3", ignore.case = TRUE)) == 0L, "no_figure3_output_created")

qc_preview <- file.path(root, "logs/revision/stage7E_figure2_final_polish_v0.5/figure2_snRNA_context_v0.5_qc_50pct.png")
check_true(file.exists(qc_preview) && file.info(qc_preview)$size > 0, "true_50_percent_QC_preview_exists")

format_info <- system2("file", c(file.path(root, required[c(1, 2, 3)]), qc_preview), stdout = TRUE, stderr = TRUE)
check_true(any(grepl("PDF document", format_info, fixed = TRUE)), "pdf_format_recognized")
check_true(any(grepl("PNG image data, 4323 x 4230", format_info, fixed = TRUE)), "png_dimensions_4323_x_4230")
check_true(any(grepl("PNG image data, 2161 x 2115", format_info, fixed = TRUE)), "qc_preview_dimensions_2161_x_2115")
check_true(any(grepl("SVG Scalable Vector Graphics", format_info, fixed = TRUE)) || any(grepl("XML", format_info, fixed = TRUE)), "svg_format_recognized")
check_true(grepl('font-family: "Helvetica"', svg, fixed = TRUE), "svg_uses_Helvetica")
pdf_font_info <- system2("strings", file.path(root, required[1]), stdout = TRUE, stderr = TRUE)
check_true(any(grepl("/BaseFont /Helvetica", pdf_font_info, fixed = TRUE)), "pdf_uses_Helvetica")
check_true(any(grepl("/BaseFont /Helvetica-Bold", pdf_font_info, fixed = TRUE)), "pdf_uses_Helvetica_Bold")

log_dir <- file.path(root, "logs/revision/stage7E_figure2_final_polish_v0.5")
log_path <- file.path(log_dir, "stage7E_figure2_v0.5_validation.log")
writeLines(c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  checks,
  "RESULT\tPASS",
  "source_data_changed_from_v0.4=no",
  "claim_changed_from_v0.4=no",
  "figure3_started=no"
), log_path, useBytes = TRUE)
cat(paste(c(checks, "RESULT\tPASS"), collapse = "\n"), "\n")
