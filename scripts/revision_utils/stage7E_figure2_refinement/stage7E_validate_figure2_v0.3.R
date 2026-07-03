#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

required <- c(
  "results/figures/revision/stage7E_figure2_refinement/figure2_snRNA_context_draft_v0.3.pdf",
  "results/figures/revision/stage7E_figure2_refinement/figure2_snRNA_context_draft_v0.3.png",
  "results/figures/revision/stage7E_figure2_refinement/figure2_snRNA_context_draft_v0.3.svg",
  "docs/revision/stage7E_figure2_refinement/figure2_qc_audit_v0.3.md",
  "results/tables/revision/stage7E_figure2_refinement/figure2_visual_change_log_v0.3.tsv",
  "results/tables/revision/stage7E_figure2_refinement/figure2_source_data_manifest_v0.3.tsv",
  "docs/revision/stage7E_figure2_refinement/figure2_legend_v0.3.md",
  "docs/revision/stage7E_figure2_refinement/stage7E_figure2_report.md"
)

checks <- character()
check_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  checks <<- c(checks, paste0("PASS\t", label))
}

check_true(all(file.exists(file.path(root, required))), "all_8_required_outputs_exist")
check_true(all(file.info(file.path(root, required))$size > 0), "all_required_outputs_nonempty")

manifest_path <- file.path(root, required[6])
change_path <- file.path(root, required[5])
manifest <- read.delim(manifest_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
changes <- read.delim(change_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

manifest_headers <- c(
  "panel", "source_data_file", "source_exists", "row_count_if_readable",
  "column_count_if_readable", "data_changed_from_v0.2", "claim_changed_from_v0.2",
  "ready_for_publication_source_data", "notes"
)
change_headers <- c(
  "panel", "visual_issue_in_v0.2", "change_made_in_v0.3", "data_changed",
  "claim_changed", "reason", "human_review_needed", "notes"
)

check_true(identical(names(manifest), manifest_headers), "source_manifest_schema_exact")
check_true(identical(names(changes), change_headers), "visual_change_log_schema_exact")
check_true(identical(manifest$panel, LETTERS[1:7]), "source_manifest_panels_A_to_G")
check_true(identical(changes$panel, LETTERS[1:7]), "visual_change_log_panels_A_to_G")
check_true(all(manifest$data_changed_from_v0.2 == "no"), "manifest_data_changed_all_no")
check_true(all(manifest$claim_changed_from_v0.2 == "no"), "manifest_claim_changed_all_no")
check_true(all(changes$data_changed == "no"), "change_log_data_changed_all_no")
check_true(all(changes$claim_changed == "no"), "change_log_claim_changed_all_no")
check_true(all(changes$human_review_needed == "yes"), "human_review_flag_all_yes")

expected_rows <- c(6L, 92L, 16L, 16L, 4L, 20L, 7L)
expected_cols <- c(6L, 17L, 14L, 12L, 19L, 17L, 9L)
check_true(identical(as.integer(manifest$row_count_if_readable), expected_rows), "manifest_row_counts_match_frozen_sources")
check_true(identical(as.integer(manifest$column_count_if_readable), expected_cols), "manifest_column_counts_match_frozen_sources")
check_true(all(file.exists(file.path(root, manifest$source_data_file))), "all_manifest_source_paths_exist")

for (i in seq_len(nrow(manifest))) {
  source <- read.delim(file.path(root, manifest$source_data_file[i]), sep = "\t", check.names = FALSE)
  check_true(nrow(source) == expected_rows[i] && ncol(source) == expected_cols[i], paste0("panel_", manifest$panel[i], "_source_dimensions_verified"))
}

svg_path <- file.path(root, required[3])
svg <- paste(readLines(svg_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
required_svg_text <- c(
  "Donor-level snRNA context mapping of MAGMA-prioritized modules",
  "Moderate Loop/TAL-associated context with partial matched-random support",
  "Bonferroni core",
  "Partial matched-random support",
  "Not causal",
  "Not plaque"
)
for (needle in required_svg_text) {
  check_true(grepl(needle, svg, fixed = TRUE), paste0("svg_contains_", gsub("[^A-Za-z0-9]+", "_", needle)))
}

forbidden_positive <- c(
  "Strong enrichment", "Validated", "Independent validation", "Therapeutic target",
  "high-confidence causal", "plaque-specific localization", "R1_MAGMA", "MAGMA_top"
)
for (needle in forbidden_positive) {
  check_true(!grepl(needle, svg, fixed = TRUE), paste0("svg_excludes_", gsub("[^A-Za-z0-9]+", "_", needle)))
}

figure_dir <- file.path(root, "results/figures/revision/stage7E_figure2_refinement")
check_true(length(list.files(figure_dir, pattern = "figure3", ignore.case = TRUE)) == 0L, "no_figure3_output_created")

format_info <- system2("file", file.path(root, required[c(1, 2, 3)]), stdout = TRUE, stderr = TRUE)
check_true(any(grepl("PDF document", format_info, fixed = TRUE)), "pdf_format_recognized")
check_true(any(grepl("PNG image data, 4323 x 4470", format_info, fixed = TRUE)), "png_dimensions_4323_x_4470")
check_true(any(grepl("SVG Scalable Vector Graphics", format_info, fixed = TRUE)) || any(grepl("XML", format_info, fixed = TRUE)), "svg_format_recognized")

log_dir <- file.path(root, "logs/revision/stage7E_figure2_refinement")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(log_dir, "stage7E_figure2_validation.log")
writeLines(c(
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  checks,
  "RESULT\tPASS",
  "source_data_changed=no",
  "claim_changed=no",
  "figure3_started=no"
), log_path, useBytes = TRUE)

cat(paste(c(checks, "RESULT\tPASS"), collapse = "\n"), "\n")
