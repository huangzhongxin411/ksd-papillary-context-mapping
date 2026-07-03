suppressPackageStartupMessages(library(data.table))

table_dir <- "results/gse73680/tables"
processed_dir <- "data/processed/gse73680"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

audit_path <- file.path(table_dir, "gse73680_txt_structure_audit.tsv")
if (!file.exists(audit_path)) stop("Missing TXT structure audit.")
audit <- fread(audit_path)
usable <- audit[usable_for_matrix == TRUE]

if (!nrow(usable)) {
  fwrite(data.table(n_remote_files = nrow(audit), n_downloaded_ok = 0L, n_structure_usable = 0L,
                    n_parsed_ok = 0L, n_failed = 0L, n_unique_samples = 0L, n_features = 0L,
                    feature_id_type = "unknown", expression_scale_guess = "unknown",
                    matrix_dimension = "0x0", duplicate_sample_ids = 0L, duplicate_feature_ids = NA_integer_,
                    status = "fail", notes = "No usable local TXT.gz files; matrix not built."),
         file.path(table_dir, "gse73680_matrix_build_log.tsv"), sep = "\t")
  fwrite(data.table(), file.path(table_dir, "gse73680_sample_sheet.tsv"), sep = "\t")
  fwrite(data.table(), file.path(table_dir, "gse73680_failed_files.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

read_agilent_features <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  start <- grep("^FEATURES\t", lines)[1]
  if (is.na(start)) stop("FEATURES section not found")
  end_candidates <- grep("^\\*$", lines)
  end_candidates <- end_candidates[end_candidates > start]
  end <- if (length(end_candidates)) end_candidates[1] - 1L else length(lines)
  section <- lines[start:end]
  if (length(section) < 2L) stop("FEATURES section is empty")
  section[1] <- sub("^FEATURES\t", "", section[1])
  data_lines <- grep("^DATA\t", section)
  if (!length(data_lines)) stop("FEATURES section has no DATA rows")
  section[data_lines] <- sub("^DATA\t", "", section[data_lines])
  fread(text = paste(section[c(1L, data_lines)], collapse = "\n"), fill = TRUE)
}

parse_one <- function(path, sample_id) {
  x <- read_agilent_features(path)
  expression_col <- intersect(c("gProcessedSignal", "gBGSubSignal", "gMeanSignal", "gMedianSignal"), names(x))[1]
  if (is.na(expression_col)) stop("no supported Agilent expression column")
  feature_col <- if ("GeneName" %in% names(x)) "GeneName" else if ("ProbeName" %in% names(x)) "ProbeName" else NA_character_
  if (is.na(feature_col)) stop("no supported feature identifier column")
  y <- x[, .(feature_id = as.character(get(feature_col)), value = as.numeric(get(expression_col)))]
  if ("ControlType" %in% names(x)) y <- y[as.integer(x$ControlType) == 0L]
  y <- y[!is.na(feature_id) & feature_id != "" & !feature_id %in% c("DarkCorner", "GE_BrightCorner")]
  y <- y[, .(value = max(value, na.rm = TRUE)), by = feature_id]
  setnames(y, "value", sample_id)
  y
}
parsed <- list(); failed <- list()
for (i in seq_len(nrow(usable))) {
  sid <- usable$sample_id_in_header[i]
  res <- tryCatch(parse_one(usable$local_path[i], sid), error = function(e) e)
  if (inherits(res, "error")) {
    failed[[length(failed) + 1]] <- data.table(file_name = usable$file_name[i], local_path = usable$local_path[i],
                                               error = conditionMessage(res))
  } else {
    parsed[[length(parsed) + 1]] <- res
  }
}
if (!length(parsed)) {
  fwrite(data.table(n_remote_files = nrow(audit), n_downloaded_ok = nrow(audit), n_structure_usable = nrow(usable),
                    n_parsed_ok = 0L, n_failed = length(failed), n_unique_samples = 0L, n_features = 0L,
                    feature_id_type = "unknown", expression_scale_guess = "unknown", matrix_dimension = "0x0",
                    duplicate_sample_ids = 0L, duplicate_feature_ids = NA_integer_, status = "fail",
                    notes = "No files parsed successfully."),
         file.path(table_dir, "gse73680_matrix_build_log.tsv"), sep = "\t")
  fwrite(rbindlist(failed, fill = TRUE), file.path(table_dir, "gse73680_failed_files.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}
mat <- Reduce(function(a, b) merge(a, b, by = "feature_id", all = TRUE), parsed)
fwrite(mat, file.path(processed_dir, "gse73680_expression_matrix.tsv.gz"), sep = "\t")
sample_cols <- setdiff(names(mat), "feature_id")
sample_sheet <- rbindlist(lapply(sample_cols, function(s) {
  vals <- mat[[s]]
  data.table(sample_id = s, patient_id = sub("^GSM[0-9]+_(ST[0-9]+)_.*$", "\\1", s), source_file = NA_character_,
             n_features = sum(!is.na(vals)), n_missing_values = sum(is.na(vals)),
             expression_min = suppressWarnings(min(vals, na.rm = TRUE)),
             expression_median = suppressWarnings(median(vals, na.rm = TRUE)),
             expression_max = suppressWarnings(max(vals, na.rm = TRUE)),
             expression_scale_guess = "pending_full_qc", include_in_matrix = TRUE,
             exclude_reason = NA_character_, notes = "Parsed into matrix.")
}))
fwrite(sample_sheet, file.path(table_dir, "gse73680_sample_sheet.tsv"), sep = "\t")
dup_features <- sum(duplicated(mat$feature_id))
status <- if (length(sample_cols) >= 6 && nrow(mat) >= 5000 && length(parsed) / max(1, nrow(usable)) >= 0.8) "pass" else if (length(sample_cols) >= 3 && nrow(mat) >= 1000) "warning" else "fail"
fwrite(data.table(n_remote_files = nrow(audit), n_downloaded_ok = nrow(audit), n_structure_usable = nrow(usable),
                  n_parsed_ok = length(parsed), n_failed = length(failed), n_unique_samples = length(sample_cols),
                  n_features = nrow(mat), feature_id_type = paste(unique(usable$feature_id_type), collapse = ";"),
                  expression_scale_guess = paste(unique(usable$expression_scale_guess), collapse = ";"),
                  matrix_dimension = paste(nrow(mat), length(sample_cols), sep = "x"),
                  duplicate_sample_ids = sum(duplicated(sample_cols)), duplicate_feature_ids = dup_features,
                  status, notes = "Matrix build completed from local supplementary TXT.gz files."),
       file.path(table_dir, "gse73680_matrix_build_log.tsv"), sep = "\t")
fwrite(rbindlist(failed, fill = TRUE), file.path(table_dir, "gse73680_failed_files.tsv"), sep = "\t")
message("wrote GSE73680 expression matrix build outputs")
