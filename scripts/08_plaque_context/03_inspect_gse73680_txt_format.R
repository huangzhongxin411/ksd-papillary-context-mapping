suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
raw_dir <- "data/raw/gse73680/supplementary"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

validation_path <- file.path(table_dir, "gse73680_extracted_txtgz_validation.tsv")
download_log_path <- file.path(table_dir, "gse73680_supplementary_download_log.tsv")
if (!file.exists(validation_path) && !file.exists(download_log_path)) {
  fwrite(data.table(file_name = character(), local_path = character(), gzip_valid = logical(),
                    n_rows_preview = integer(), n_cols_preview = integer(), first_column_name = character(),
                    first_column_example = character(), header_detected = logical(), numeric_columns_detected = logical(),
                    n_numeric_columns = integer(), sample_id_in_header = character(), patient_id_detected = character(),
                    feature_id_type = character(), matrix_layout = character(), expression_scale_guess = character(),
                    file_type_class = character(), usable_for_matrix = logical(), exclude_reason = character(), notes = character()),
         file.path(table_dir, "gse73680_txt_structure_audit.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
if (file.exists(validation_path)) {
  validation <- fread(validation_path)
  ok <- validation[gzip_valid == TRUE, .(
    remote_name = basename(file),
    local_path = file,
    readable_header = TRUE
  )]
} else {
  dl <- fread(download_log_path)
  ok <- dl[gzip_valid == TRUE & readable_header == TRUE]
}
empty_audit <- function() data.table(
  file_name = character(), local_path = character(), gzip_valid = logical(),
  n_rows_preview = integer(), n_cols_preview = integer(), first_column_name = character(),
  first_column_example = character(), header_detected = logical(), numeric_columns_detected = logical(),
  n_numeric_columns = integer(), sample_id_in_header = character(), patient_id_detected = character(),
  feature_id_type = character(), matrix_layout = character(), expression_scale_guess = character(),
  file_type_class = character(), usable_for_matrix = logical(), exclude_reason = character(), notes = character()
)

read_agilent_features <- function(path, nrows = 100L) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  start <- grep("^FEATURES\t", lines)[1]
  if (is.na(start)) return(NULL)
  end_candidates <- grep("^\\*$", lines)
  end_candidates <- end_candidates[end_candidates > start]
  end <- if (length(end_candidates)) end_candidates[1] - 1L else length(lines)
  section <- lines[start:end]
  if (length(section) < 2L) return(NULL)
  section[1] <- sub("^FEATURES\t", "", section[1])
  data_lines <- grep("^DATA\t", section)
  if (!length(data_lines)) return(NULL)
  section[data_lines] <- sub("^DATA\t", "", section[data_lines])
  section <- section[c(1L, data_lines)]
  if (is.finite(nrows)) section <- section[seq_len(min(length(section), nrows + 1L))]
  fread(text = paste(section, collapse = "\n"), fill = TRUE)
}

if (!nrow(ok)) {
  audit <- empty_audit()
  fwrite(audit, file.path(table_dir, "gse73680_txt_structure_audit.tsv"), sep = "\t")
  fwrite(data.table(file_name = character(), header = character()), file.path(table_dir, "gse73680_example_headers.tsv"), sep = "\t")
  fwrite(data.table(file_type_class = character(), matrix_layout = character(), feature_id_type = character(),
                    expression_scale_guess = character(), usable_for_matrix = logical(), N = integer()),
         file.path(table_dir, "gse73680_file_type_classification.tsv"), sep = "\t")
  message("No valid local TXT.gz files for structure audit.")
  quit(save = "no", status = 0)
}

inspect_one <- function(path, name) {
  preview <- tryCatch(read_agilent_features(path, nrows = 100L), error = function(e) NULL)
  if (is.null(preview) || !ncol(preview)) {
    return(data.table(file_name = name, local_path = path, gzip_valid = TRUE, n_rows_preview = 0L,
                      n_cols_preview = 0L, first_column_name = NA_character_, first_column_example = NA_character_,
                      header_detected = FALSE, numeric_columns_detected = FALSE, n_numeric_columns = 0L,
                      sample_id_in_header = NA_character_, patient_id_detected = NA_character_,
                      feature_id_type = "unknown", matrix_layout = "unknown", expression_scale_guess = "unknown",
                      file_type_class = "unknown", usable_for_matrix = FALSE,
                      exclude_reason = "preview_read_failed", notes = "Could not parse preview."))
  }
  expression_candidates <- intersect(c("gProcessedSignal", "gBGSubSignal", "gMeanSignal", "gMedianSignal"), names(preview))
  numeric_cols <- expression_candidates[vapply(preview[, ..expression_candidates], is.numeric, logical(1))]
  feature_col <- if ("GeneName" %in% names(preview)) "GeneName" else if ("ProbeName" %in% names(preview)) "ProbeName" else names(preview)[1]
  first_vals <- preview[[feature_col]]
  first_ex <- as.character(first_vals[which(!is.na(first_vals) & first_vals != "")[1]])
  feature_id_type <- if (feature_col == "GeneName") "gene_symbol" else if (feature_col == "ProbeName") "probe_id" else "unknown"
  matrix_layout <- if (length(numeric_cols) >= 1) "agilent_single_sample_feature_vector" else "unknown"
  vals <- unlist(preview[, ..numeric_cols], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  scale_guess <- if (!length(vals)) "unknown" else if (max(vals, na.rm = TRUE) <= 25 && min(vals, na.rm = TRUE) >= -5) "log2_like" else if (all(abs(vals - round(vals)) < 1e-8) && max(vals, na.rm = TRUE) > 100) "raw_count_like" else "normalized_continuous"
  usable <- matrix_layout %in% c("single_sample_feature_vector", "multi_sample_expression_matrix", "agilent_single_sample_feature_vector") && length(numeric_cols) >= 1 && feature_id_type != "unknown"
  data.table(file_name = name, local_path = path, gzip_valid = TRUE,
             n_rows_preview = nrow(preview), n_cols_preview = ncol(preview),
             first_column_name = feature_col, first_column_example = first_ex,
             header_detected = TRUE, numeric_columns_detected = length(numeric_cols) > 0,
             n_numeric_columns = length(numeric_cols), sample_id_in_header = sub("_.*$", "", name),
             patient_id_detected = sub("^GSM[0-9]+_(ST[0-9]+)_.*$", "\\1", name),
             feature_id_type, matrix_layout, expression_scale_guess = scale_guess,
             file_type_class = if (usable) "expression_vector" else "unknown",
             usable_for_matrix = usable,
             exclude_reason = if (usable) NA_character_ else "no_numeric_expression_column_or_unknown_feature_id",
             notes = "Preview-based structure audit; full parsing occurs during matrix build.")
}
audit <- rbindlist(Map(inspect_one, ok$local_path, ok$remote_name), fill = TRUE)
fwrite(audit, file.path(table_dir, "gse73680_txt_structure_audit.tsv"), sep = "\t")
headers <- rbindlist(lapply(ok$local_path[seq_len(min(5, nrow(ok)))], function(p) {
  x <- tryCatch(fread(p, nrows = 3, fill = TRUE), error = function(e) data.table())
  data.table(file_name = basename(p), header = paste(names(x), collapse = "|"))
}), fill = TRUE)
fwrite(headers, file.path(table_dir, "gse73680_example_headers.tsv"), sep = "\t")
fwrite(audit[, .N, by = .(file_type_class, matrix_layout, feature_id_type, expression_scale_guess, usable_for_matrix)],
       file.path(table_dir, "gse73680_file_type_classification.tsv"), sep = "\t")
message("wrote GSE73680 TXT structure audit")
