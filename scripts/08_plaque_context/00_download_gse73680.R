suppressPackageStartupMessages(library(data.table))
options(timeout = 600)

dir.create("data/raw/GSE73680", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

download_one <- function(url, dest) {
  status <- "not_attempted"
  error <- NA_character_
  bytes <- NA_real_
  tmp <- paste0(dest, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  tryCatch({
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
    if (file.exists(tmp) && file.info(tmp)$size > 0) {
      ok <- TRUE
      if (grepl("\\.gz$", dest, ignore.case = TRUE)) {
        ok <- system2("gzip", c("-t", tmp), stdout = TRUE, stderr = TRUE) |> length() == 0
      }
      if (ok) {
        file.rename(tmp, dest)
        status <- "downloaded"
        bytes <- file.info(dest)$size
      } else {
        status <- "failed_integrity_check"
        error <- "gzip integrity check failed"
        unlink(tmp)
      }
    } else {
      status <- "empty_or_missing"
    }
  }, error = function(e) {
    status <<- "failed"
    error <<- conditionMessage(e)
    if (file.exists(tmp)) unlink(tmp)
  })
  data.table(url = url, local_path = dest, status = status, bytes = bytes, error = error)
}

resources <- data.table(
  resource_type = c("series_matrix", "supplementary_filelist", "series_soft"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73680/matrix/GSE73680_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73680/suppl/filelist.txt",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73680&targ=self&form=text&view=full"
  ),
  local_path = c(
    "data/raw/GSE73680/GSE73680_series_matrix.txt.gz",
    "data/raw/GSE73680/GSE73680_suppl_filelist.txt",
    "data/raw/GSE73680/GSE73680_series_full.soft.txt"
  )
)

status <- rbindlist(Map(download_one, resources$url, resources$local_path))
status <- merge(resources, status, by = c("url", "local_path"), all.x = TRUE)
fwrite(status, "results/tables/gse73680_download_status.tsv", sep = "\t")

files <- list.files("data/raw/GSE73680", recursive = TRUE, full.names = TRUE)
base <- basename(files)
has_local <- function(pattern) any(grepl(pattern, base, ignore.case = TRUE))
valid_gzip <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(TRUE)
  length(system2("gzip", c("-t", path), stdout = TRUE, stderr = TRUE)) == 0
}

file_decision <- data.table(
  resource_type = c("Series matrix", "processed expression matrix", "platform annotation", "sample metadata", "supplementary processed files", "raw CEL or FASTQ"),
  expected_file = c(
    "GSE73680_series_matrix.txt.gz",
    "processed expression matrix in GEO supplementary files",
    "platform annotation table or GPL reference",
    "sample metadata in series matrix or SOFT",
    "GEO supplementary processed files",
    "raw CEL or FASTQ"
  ),
  local_detected = c(
    valid_gzip("data/raw/GSE73680/GSE73680_series_matrix.txt.gz"),
    has_local("expression|matrix|normalized|processed|counts"),
    has_local("GPL|annot|platform"),
    has_local("series_matrix|soft|metadata|sample"),
    file.exists("data/raw/GSE73680/GSE73680_suppl_filelist.txt") && file.info("data/raw/GSE73680/GSE73680_suppl_filelist.txt")$size > 0,
    has_local("\\.(cel|fastq|fq)(\\.gz)?$")
  )
)
file_decision[, download_required := !local_detected]
file_decision[, usable_for_analysis := local_detected & resource_type %in% c("Series matrix", "processed expression matrix", "sample metadata", "supplementary processed files")]
file_decision[, notes := fifelse(
  local_detected,
  "Local resource detected; manual metadata and matrix-format audit required before disease-context analysis.",
  "Not detected locally; acquire processed expression and metadata before disease-context analysis."
)]
fwrite(file_decision, "results/tables/gse73680_file_decision.tsv", sep = "\t")

series_path <- "data/raw/GSE73680/GSE73680_series_matrix.txt.gz"
if (valid_gzip(series_path)) {
  lines <- readLines(gzfile(series_path), n = 2000, warn = FALSE)
  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  sample_ids <- unique(unlist(regmatches(lines, gregexpr("GSM[0-9]+", lines))))
  metadata_audit <- data.table(
    sample_id = sample_ids,
    geo_accession = sample_ids,
    group_raw = NA_character_,
    group_curated = NA_character_,
    tissue_context = NA_character_,
    disease_status = NA_character_,
    plaque_status = NA_character_,
    stone_status = NA_character_,
    sample_type = NA_character_,
    include_in_analysis = FALSE,
    exclude_reason = "metadata_group_not_curated",
    notes = "Metadata extracted only at resource-landing stage; manual curation required before differential analysis."
  )
  fwrite(metadata_audit, "results/tables/gse73680_metadata_audit.tsv", sep = "\t")
  fwrite(metadata_audit, "config/gse73680_sample_metadata_curated.tsv", sep = "\t")

  table_start <- grep("^!series_matrix_table_begin", lines)
  table_end <- grep("^!series_matrix_table_end", lines)
  expression_qc <- data.table(
    dataset = "GSE73680",
    series_matrix_available = TRUE,
    n_samples_detected_in_header = length(sample_ids),
    sample_metadata_lines = length(sample_lines),
    expression_table_detected_in_first_2000_lines = length(table_start) > 0,
    n_features = NA_integer_,
    missing_rate = NA_real_,
    log2_status = "not_audited",
    normalization_status = "not_audited",
    platform_annotation = "not_audited",
    notes = "Series matrix landed; full expression parsing and probe-to-gene annotation require metadata/platform audit."
  )
  fwrite(expression_qc, "results/tables/gse73680_expression_qc.tsv", sep = "\t")

  gene_sets <- data.table(
    gene_set = c("P1_core_TAL_candidates", "MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "TAL_marker_set", "injury_remodeling_marker_set"),
    n_input_genes = c(6L, 50L, 100L, NA_integer_, NA_integer_, NA_integer_, 12L),
    n_detected = NA_integer_,
    detected_fraction = NA_real_,
    usability = "pending_probe_to_gene_annotation",
    notes = "Gene availability cannot be finalized until platform/probe-to-gene mapping is audited."
  )
  fwrite(gene_sets, "results/tables/gse73680_gene_availability.tsv", sep = "\t")
} else {
  fwrite(data.table(
    dataset = "GSE73680",
    series_matrix_available = FALSE,
    n_samples_detected_in_header = 0L,
    sample_metadata_lines = 0L,
    expression_table_detected_in_first_2000_lines = FALSE,
    n_features = NA_integer_,
    missing_rate = NA_real_,
    log2_status = "not_available",
    normalization_status = "not_available",
    platform_annotation = "not_available",
    notes = "Series matrix not available locally; no expression QC or disease-context analysis was run."
  ), "results/tables/gse73680_expression_qc.tsv", sep = "\t")
}

message("wrote GSE73680 download/resource landing status")
