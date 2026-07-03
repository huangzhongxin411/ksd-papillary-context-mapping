suppressPackageStartupMessages(library(data.table))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

filelist_path <- "data/raw/GSE73680/GSE73680_suppl_filelist.txt"
if (!file.exists(filelist_path)) {
  stop("Missing supplementary filelist. Run 00_download_gse73680.R first.")
}

fl <- fread(filelist_path, fill = TRUE)
setnames(fl, gsub("^#", "", names(fl)))
setnames(fl, names(fl)[1], "record_type")
files <- fl[record_type == "File"]
files[, local_path := file.path("data/raw/GSE73680/suppl", Name)]
files[, local_detected := file.exists(local_path)]
files[, geo_accession := sub("_.*$", "", Name)]
files[, stem := sub("\\.txt\\.gz$", "", Name)]
files[, patient_id := sub("^GSM[0-9]+_(ST[0-9]+)_.*$", "\\1", stem)]
files[, sex := fifelse(grepl("_M_", stem), "M", fifelse(grepl("_F_", stem), "F", NA_character_))]
files[, filename_group_token := fifelse(grepl("_N($|_)", stem), "N",
                                  fifelse(grepl("_P_", stem), "P",
                                          fifelse(grepl("_C_", stem), "C", "other")))]
files[, plaque_status := fifelse(filename_group_token == "P", "plaque_or_papillary_plaque", NA_character_)]
files[, disease_status := fifelse(filename_group_token == "N", "normal_context",
                           fifelse(filename_group_token == "P", "stone_plaque_context",
                                   fifelse(filename_group_token == "C", "control_or_other_context", "uncurated")))]
files[, stone_type_from_filename := fifelse(grepl("CaOx", stem) & grepl("CaP", stem), "CaOx_CaP",
                                     fifelse(grepl("CaOx", stem), "CaOx",
                                             fifelse(grepl("CaP", stem), "CaP",
                                                     fifelse(grepl("cystine", stem, ignore.case = TRUE), "cystine", NA_character_))))]
files[, download_required := !local_detected]
files[, notes := "Remote GEO supplementary TXT.gz listed; filename-derived metadata is preliminary and must be manually curated before analysis."]
fwrite(files[, .(
  Name, Time, Size, Type, geo_accession, patient_id, sex, filename_group_token,
  disease_status, plaque_status, stone_type_from_filename, local_path,
  local_detected, download_required, notes
)], "results/tables/gse73680_supplementary_file_manifest.tsv", sep = "\t")

metadata <- files[, .(
  sample_id = geo_accession,
  geo_accession,
  group_raw = stem,
  group_curated = disease_status,
  tissue_context = fifelse(filename_group_token == "P", "papillary_plaque_or_plaque_context",
                    fifelse(filename_group_token == "N", "normal_papilla_context", "control_or_other_context")),
  disease_status,
  plaque_status,
  stone_status = stone_type_from_filename,
  sample_type = "supplementary_TXT_from_GEO_filename",
  include_in_analysis = FALSE,
  exclude_reason = "requires_manual_metadata_and_expression_format_audit",
  notes = "Curated from supplementary filenames only; not yet approved for differential or disease-context analysis."
)]
fwrite(metadata, "results/tables/gse73680_metadata_audit.tsv", sep = "\t")
fwrite(metadata, "config/gse73680_sample_metadata_curated.tsv", sep = "\t")

series_matrix_ok <- file.exists("data/raw/GSE73680/GSE73680_series_matrix.txt.gz") &&
  file.info("data/raw/GSE73680/GSE73680_series_matrix.txt.gz")$size > 0 &&
  length(system2("gzip", c("-t", "data/raw/GSE73680/GSE73680_series_matrix.txt.gz"), stdout = TRUE, stderr = TRUE)) == 0
local_suppl_txt <- any(files$local_detected)
soft_ok <- file.exists("data/raw/GSE73680/GSE73680_series_full.soft.txt") && file.info("data/raw/GSE73680/GSE73680_series_full.soft.txt")$size > 0

file_decision <- data.table(
  resource_type = c("Series matrix", "remote supplementary filelist", "local supplementary TXT.gz files", "platform annotation", "sample metadata", "raw CEL or FASTQ"),
  expected_file = c(
    "GSE73680_series_matrix.txt.gz",
    "GSE73680_suppl_filelist.txt",
    "GSM*.txt.gz per-sample supplementary files",
    "GPL annotation table or platform metadata",
    "SOFT metadata plus filename-derived sample labels",
    "raw CEL or FASTQ"
  ),
  local_detected = c(series_matrix_ok, file.exists(filelist_path), local_suppl_txt, FALSE, soft_ok || file.exists(filelist_path), FALSE),
  download_required = c(!series_matrix_ok, FALSE, !local_suppl_txt, TRUE, !(soft_ok || file.exists(filelist_path)), TRUE),
  usable_for_analysis = c(series_matrix_ok, FALSE, FALSE, FALSE, FALSE, FALSE),
  notes = c(
    if (series_matrix_ok) "Complete series matrix available locally." else "Series matrix download attempted but not available as a valid gzip.",
    "Remote supplementary filelist is available and has been parsed; it is not itself an expression matrix.",
    "Per-sample supplementary TXT.gz files are listed remotely but are not downloaded locally.",
    "Platform/probe annotation is still required before gene-level availability can be finalized.",
    "Metadata is preliminary and filename-derived; manual curation required before analysis.",
    "Raw data are not prioritized for current scope."
  )
)
fwrite(file_decision, "results/tables/gse73680_file_decision.tsv", sep = "\t")

fwrite(data.table(
  dataset = "GSE73680",
  series_matrix_available = series_matrix_ok,
  supplementary_filelist_available = TRUE,
  n_remote_supplementary_txt = nrow(files),
  n_local_supplementary_txt = sum(files$local_detected),
  n_samples_from_filename = uniqueN(files$geo_accession),
  n_patients_from_filename = uniqueN(files$patient_id),
  expression_matrix_ready = series_matrix_ok || local_suppl_txt,
  metadata_ready_for_analysis = FALSE,
  notes = "Resource landing has remote supplementary index and preliminary filename metadata, but no local expression matrix is ready for analysis."
), "results/tables/gse73680_expression_qc.tsv", sep = "\t")

message("wrote GSE73680 filelist audit and corrected resource decisions")
