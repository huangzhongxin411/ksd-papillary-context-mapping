suppressPackageStartupMessages(library(data.table))
options(timeout = 240)

raw_dir <- "data/raw/gse73680/supplementary"
table_dir <- "results/gse73680/tables"
log_dir <- "results/gse73680/logs"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- if (file.exists(file.path(table_dir, "gse73680_supplementary_file_manifest.tsv"))) {
  file.path(table_dir, "gse73680_supplementary_file_manifest.tsv")
} else {
  "results/tables/gse73680_supplementary_file_manifest.tsv"
}
if (!file.exists(manifest_path)) stop("Missing GSE73680 supplementary manifest.")
manifest <- fread(manifest_path)

max_files <- as.integer(Sys.getenv("GSE73680_MAX_FILES", "0"))
if (!is.na(max_files) && max_files > 0) manifest <- manifest[seq_len(min(.N, max_files))]

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73680/suppl"
raw_tar_url <- paste0(base_url, "/GSE73680_RAW.tar")
raw_tar_path <- "data/raw/gse73680/GSE73680_RAW.tar"
valid_gzip <- function(path) {
  file.exists(path) && file.info(path)$size > 0 &&
    length(system2("gzip", c("-t", path), stdout = TRUE, stderr = TRUE)) == 0
}
readable_header <- function(path) {
  if (!valid_gzip(path)) return(FALSE)
  out <- tryCatch({ readLines(gzfile(path), n = 2, warn = FALSE); TRUE }, error = function(e) FALSE)
  out
}
download_one <- function(name) {
  remote_url <- paste0(base_url, "/", name)
  local_path <- file.path(raw_dir, name)
  tmp <- paste0(local_path, ".tmp")
  if (file.exists(local_path) && valid_gzip(local_path)) {
    return(data.table(remote_name = name, remote_url, local_path, download_attempted = FALSE,
                      download_status = "ok_existing", file_size_bytes = file.info(local_path)$size,
                      gzip_valid = TRUE, readable_header = readable_header(local_path),
                      parsed_patient_id = sub("^GSM[0-9]+_(ST[0-9]+)_.*$", "\\1", name),
                      parsed_sample_id = sub("_.*$", "", name),
                      notes = "Existing valid gzip reused."))
  }
  if (file.exists(tmp)) unlink(tmp)
  status <- "failed"
  notes <- NA_character_
  tryCatch({
    utils::download.file(remote_url, tmp, mode = "wb", quiet = TRUE)
    if (valid_gzip(tmp)) {
      file.rename(tmp, local_path)
      status <- "ok"
      notes <- "Downloaded and gzip-validated."
    } else {
      status <- "partial"
      notes <- "Downloaded file failed gzip validation; excluded from matrix construction."
      if (file.exists(tmp)) unlink(tmp)
    }
  }, error = function(e) {
    status <<- "failed"
    notes <<- conditionMessage(e)
    if (file.exists(tmp)) unlink(tmp)
  })
  data.table(remote_name = name, remote_url, local_path, download_attempted = TRUE,
             download_status = status,
             file_size_bytes = if (file.exists(local_path)) file.info(local_path)$size else 0,
             gzip_valid = valid_gzip(local_path),
             readable_header = readable_header(local_path),
             parsed_patient_id = sub("^GSM[0-9]+_(ST[0-9]+)_.*$", "\\1", name),
             parsed_sample_id = sub("_.*$", "", name),
             notes = notes)
}

log <- rbindlist(lapply(manifest$Name, download_one), fill = TRUE)
if (!any(log$gzip_valid & log$readable_header)) {
  raw_tar_exists <- file.exists(raw_tar_path) && file.info(raw_tar_path)$size > 0
  raw_tar_size <- if (raw_tar_exists) file.info(raw_tar_path)$size else 0
  log[, notes := paste(notes, "Individual TXT.gz files appear to be archive members; use GSE73680_RAW.tar for full extraction.")]
  log <- rbind(
    log,
    data.table(remote_name = "GSE73680_RAW.tar", remote_url = raw_tar_url,
               local_path = raw_tar_path, download_attempted = FALSE,
               download_status = if (raw_tar_exists) "raw_tar_present_not_extracted" else "requires_raw_tar_download",
               file_size_bytes = raw_tar_size, gzip_valid = NA, readable_header = NA,
               parsed_patient_id = NA_character_, parsed_sample_id = NA_character_,
               notes = "NCBI filelist entries are inside this 803399680-byte TAR archive; download/extract TAR before matrix construction."),
    fill = TRUE
  )
}
fwrite(log, file.path(table_dir, "gse73680_supplementary_download_log.tsv"), sep = "\t")
writeLines(c(
  paste("timestamp", Sys.time()),
  paste("manifest", manifest_path),
  paste("n_requested", nrow(manifest)),
  paste("n_ok", sum(log$gzip_valid & log$readable_header))
), file.path(log_dir, "download_gse73680_supplementary.log"))
message("wrote ", file.path(table_dir, "gse73680_supplementary_download_log.tsv"))
