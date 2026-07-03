#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c(
  "results/figures/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_snRNA_context_draft_v0.6.pdf",
  "results/figures/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_snRNA_context_draft_v0.6.png",
  "results/figures/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_snRNA_context_draft_v0.6.svg",
  "docs/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_qc_audit_v0.6.md",
  "results/tables/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_visual_redesign_log_v0.6.tsv",
  "results/tables/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_source_data_manifest_v0.6.tsv",
  "docs/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_legend_v0.6.md",
  "docs/revision/stage7E_figure2_editorial_redesign_v0.6/stage7E_figure2_editorial_redesign_report.md"
)

checks <- character()
check_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  checks <<- c(checks, paste0("PASS\t", label))
}

check_true(all(file.exists(file.path(root, required))), "all_8_required_outputs_exist")
check_true(all(file.info(file.path(root, required))$size > 0), "all_required_outputs_nonempty")

manifest <- read.delim(file.path(root, required[6]), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
redesign <- read.delim(file.path(root, required[5]), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
manifest_headers <- c(
  "panel", "source_data_file", "source_exists", "row_count_if_readable",
  "column_count_if_readable", "data_changed_from_v0.5", "claim_changed_from_v0.5",
  "main_figure_or_source_detail", "notes"
)
redesign_headers <- c(
  "old_panel", "new_panel", "issue_in_v0.5", "redesign_action_v0.6",
  "data_changed", "claim_changed", "reason", "notes"
)
check_true(identical(names(manifest), manifest_headers), "source_manifest_schema_exact")
check_true(identical(names(redesign), redesign_headers), "visual_redesign_log_schema_exact")
check_true(identical(manifest$panel, LETTERS[1:5]), "source_manifest_panels_A_to_E")
check_true(identical(redesign$new_panel, LETTERS[1:5]), "redesign_log_new_panels_A_to_E")
check_true(all(manifest$data_changed_from_v0.5 == "no"), "manifest_data_changed_all_no")
check_true(all(manifest$claim_changed_from_v0.5 == "no"), "manifest_claim_changed_all_no")
check_true(all(redesign$data_changed == "no"), "redesign_data_changed_all_no")
check_true(all(redesign$claim_changed == "no"), "redesign_claim_changed_all_no")

expected_rows <- list(c(6L), c(92L), c(16L, 16L, 20L), c(4L), c(7L))
expected_cols <- list(c(6L), c(17L), c(14L, 12L, 17L), c(19L), c(9L))
for (i in seq_len(nrow(manifest))) {
  paths <- trimws(strsplit(manifest$source_data_file[i], ";", fixed = TRUE)[[1]])
  row_counts <- as.integer(trimws(strsplit(manifest$row_count_if_readable[i], ";", fixed = TRUE)[[1]]))
  col_counts <- as.integer(trimws(strsplit(manifest$column_count_if_readable[i], ";", fixed = TRUE)[[1]]))
  exists_flags <- trimws(strsplit(manifest$source_exists[i], ";", fixed = TRUE)[[1]])
  check_true(length(paths) == length(expected_rows[[i]]), paste0("panel_", manifest$panel[i], "_source_count_expected"))
  check_true(all(exists_flags == "yes") && all(file.exists(file.path(root, paths))), paste0("panel_", manifest$panel[i], "_all_source_paths_exist"))
  check_true(identical(row_counts, expected_rows[[i]]) && identical(col_counts, expected_cols[[i]]), paste0("panel_", manifest$panel[i], "_manifest_dimensions_expected"))
  for (j in seq_along(paths)) {
    source <- read.delim(file.path(root, paths[j]), sep = "\t", check.names = FALSE)
    check_true(nrow(source) == row_counts[j] && ncol(source) == col_counts[j], paste0("panel_", manifest$panel[i], "_source_", j, "_dimensions_verified"))
  }
}
check_true(grepl("detailed_matrices_in_source_data", manifest$main_figure_or_source_detail[manifest$panel == "C"], fixed = TRUE), "panel_C_detail_deferred_in_manifest")

svg <- paste(readLines(file.path(root, required[3]), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
required_text <- c(
  "Donor-level snRNA context mapping of MAGMA-prioritized modules",
  "Moderate Loop/TAL-associated context with partial matched-random support",
  "Resource and", "analysis unit", "MAGMA-prioritized", "modules",
  "GSE231569 papilla snRNA", "43,878", "540", "4", "6",
  "Donor × compartment module scores", "Consistency and robustness", "summary",
  "Top-ranked in 4/4 donors", "Retained after 4/4 donor exclusions",
  "Retained after 5/5 removal panels", "Detailed matrices in Source Data",
  "Matched-random benchmark", "Overall: partial support", "Interpretation boundary",
  "Not causal", "Not plaque"
)
for (needle in required_text) {
  check_true(grepl(needle, svg, fixed = TRUE), paste0("svg_contains_", gsub("[^A-Za-z0-9]+", "_", needle)))
}
for (panel in LETTERS[1:5]) {
  check_true(grepl(paste0(">", panel, "</text>"), svg, fixed = TRUE), paste0("svg_contains_panel_tag_", panel))
}
check_true(!grepl(">F</text>", svg, fixed = TRUE) && !grepl(">G</text>", svg, fixed = TRUE), "svg_has_only_five_panel_tags")

forbidden_positive <- c(
  "Strong enrichment", "Validated", "Independent validation", "Therapeutic target",
  "high-confidence causal", "plaque-specific localization", "R1_MAGMA", "MAGMA_top",
  "GWAS discovery", "fine mapping", "causal loci", "genetic validation", "Manhattan"
)
for (needle in forbidden_positive) {
  check_true(!grepl(needle, svg, fixed = TRUE), paste0("svg_excludes_", gsub("[^A-Za-z0-9]+", "_", needle)))
}

qc_text <- readLines(file.path(root, required[4]), warn = FALSE)
check_true(sum(grepl("\\| PASS \\|", qc_text)) == 12L, "all_12_requested_qc_items_pass")
check_true(any(grepl("Detailed matrices in Source Data", qc_text, fixed = TRUE)), "detail_deferral_documented")

figure_dir <- file.path(root, "results/figures/revision/stage7E_figure2_editorial_redesign_v0.6")
check_true(length(list.files(figure_dir, pattern = "figure3", ignore.case = TRUE)) == 0L, "no_figure3_output_created")

qc_preview <- file.path(root, "logs/revision/stage7E_figure2_editorial_redesign_v0.6/figure2_snRNA_context_v0.6_qc_50pct.png")
check_true(file.exists(qc_preview) && file.info(qc_preview)$size > 0, "true_50_percent_QC_preview_exists")
format_info <- system2("file", c(file.path(root, required[c(1, 2, 3)]), qc_preview), stdout = TRUE, stderr = TRUE)
check_true(any(grepl("PDF document", format_info, fixed = TRUE)), "pdf_format_recognized")
check_true(any(grepl("PNG image data, 4323 x 3870", format_info, fixed = TRUE)), "png_dimensions_4323_x_3870")
check_true(any(grepl("PNG image data, 2161 x 1935", format_info, fixed = TRUE)), "qc_preview_dimensions_2161_x_1935")
check_true(any(grepl("SVG Scalable Vector Graphics", format_info, fixed = TRUE)) || any(grepl("XML", format_info, fixed = TRUE)), "svg_format_recognized")
check_true(grepl('font-family: "Helvetica"', svg, fixed = TRUE), "svg_uses_Helvetica")
pdf_font_info <- system2("strings", file.path(root, required[1]), stdout = TRUE, stderr = TRUE)
check_true(any(grepl("/BaseFont /Helvetica", pdf_font_info, fixed = TRUE)), "pdf_uses_Helvetica")
check_true(any(grepl("/BaseFont /Helvetica-Bold", pdf_font_info, fixed = TRUE)), "pdf_uses_Helvetica_Bold")

log_dir <- file.path(root, "logs/revision/stage7E_figure2_editorial_redesign_v0.6")
log_path <- file.path(log_dir, "stage7E_figure2_v0.6_validation.log")
writeLines(c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  checks,
  "RESULT\tPASS",
  "source_data_changed_from_v0.5=no",
  "claim_changed_from_v0.5=no",
  "full_manhattan_plot_added=no",
  "figure3_started=no"
), log_path, useBytes = TRUE)
cat(paste(c(checks, "RESULT\tPASS"), collapse = "\n"), "\n")
