suppressPackageStartupMessages(library(data.table))

manifest <- fread("results/tables/eqtl_resource_manifest_v0.2.tsv")
available <- manifest[status != "missing"]

if (!nrow(available)) {
  fwrite(data.table(
    status = "not_ready",
    notes = "No eQTL resources are available locally; coloc eQTL preparation not run."
  ), "results/tables/eqtl_for_coloc_preparation_status.tsv", sep = "\t")
  message("No eQTL resources available; wrote status table.")
  quit(save = "no", status = 0)
}

stop("eQTL resource preparation parser must be implemented for the specific resource format once files are available.")
