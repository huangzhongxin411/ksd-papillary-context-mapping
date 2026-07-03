suppressPackageStartupMessages(library(data.table))

twas_input <- "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
resource_manifest <- "results/tables/twas_resource_manifest.tsv"
out_by_tissue <- "results/tables/twas_harmonization_by_tissue.tsv"
out_summary <- "results/tables/twas_harmonization_summary.tsv"

dt <- fread(twas_input, select = c("SNP", "A1", "A2"))
input_snps <- uniqueN(dt$SNP)
palindromic <- dt[paste0(A1, A2) %in% c("AT", "TA", "CG", "GC"), uniqueN(SNP)]

manifest <- fread(resource_manifest)

by_tissue <- manifest[, .(
  method,
  tissue,
  input_snps = input_snps,
  matched_snps = NA_integer_,
  unmatched_snps = NA_integer_,
  palindromic_snps = palindromic,
  flipped_snps = NA_integer_,
  ambiguous_snps = NA_integer_,
  retained_snps = NA_integer_,
  retention_rate = NA_real_,
  harmonization_status = ifelse(status == "ready", "pending_model_level_check", "not_evaluable_resource_missing"),
  notes = ifelse(status == "ready",
    "Resource appears ready; run model-specific SNP/allele matching before TWAS interpretation.",
    paste("Resource status:", status, "-", notes)
  )
)]

summary <- by_tissue[, .(
  n_tissues = .N,
  n_ready = sum(harmonization_status == "pending_model_level_check"),
  n_not_evaluable = sum(harmonization_status == "not_evaluable_resource_missing"),
  input_snps = input_snps[1],
  palindromic_snps_in_input = palindromic[1],
  min_retention_rate = suppressWarnings(min(retention_rate, na.rm = TRUE)),
  notes = "True harmonization requires FUSION weight SNPs or PredictDB covariance/model SNPs; current missing resources prevent allele matching."
), by = method]

summary[is.infinite(min_retention_rate), min_retention_rate := NA_real_]

fwrite(by_tissue, out_by_tissue, sep = "\t")
fwrite(summary, out_summary, sep = "\t")

message("wrote\t", out_by_tissue)
message("wrote\t", out_summary)
