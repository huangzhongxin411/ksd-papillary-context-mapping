suppressPackageStartupMessages(library(data.table))

infile <- "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
outfile <- "data/processed/twas_input/ksd_2025_for_fusion.tsv"
qc_out <- "results/tables/fusion_sumstats_qc.tsv"

dt <- fread(infile)
required <- c("SNP", "A1", "A2", "Z", "P", "N")
missing <- setdiff(required, names(dt))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

raw_n <- nrow(dt)
out <- dt[, ..required]
out <- out[!is.na(SNP) & !is.na(A1) & !is.na(A2) & !is.na(Z) & !is.na(P) & !is.na(N)]
missing_required_removed <- raw_n - nrow(out)
duplicated_removed <- nrow(out) - uniqueN(out$SNP)
out <- unique(out, by = "SNP")

fwrite(out, outfile, sep = "\t")

qc <- data.table(
  metric = c(
    "input_file",
    "output_file",
    "input_rows",
    "output_rows",
    "duplicated_removed",
    "missing_required_removed",
    "min_p",
    "max_abs_z",
    "min_n",
    "median_n",
    "max_n"
  ),
  value = c(
    infile,
    outfile,
    as.character(raw_n),
    as.character(nrow(out)),
    as.character(duplicated_removed),
    as.character(missing_required_removed),
    format(min(out$P, na.rm = TRUE), scientific = TRUE),
    format(max(abs(out$Z), na.rm = TRUE), scientific = TRUE),
    as.character(min(out$N, na.rm = TRUE)),
    as.character(stats::median(out$N, na.rm = TRUE)),
    as.character(max(out$N, na.rm = TRUE))
  )
)

fwrite(qc, qc_out, sep = "\t")
message("wrote\t", outfile)
message("wrote\t", qc_out)
