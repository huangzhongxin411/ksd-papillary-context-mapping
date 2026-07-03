suppressPackageStartupMessages(library(data.table))

expected_path <- "results/tables/spatial_expected_files.tsv"
detection_path <- "results/tables/spatial_file_detection.tsv"
out_path <- "results/tables/spatial_input_inventory_v0.1.tsv"

if (!file.exists(expected_path)) {
  stop("Missing spatial_expected_files.tsv.")
}

expected <- fread(expected_path)
if (file.exists(detection_path)) {
  detection <- fread(detection_path)
} else {
  detection <- data.table(dataset = unique(expected$dataset))
  detection[, `:=`(
    local_root = file.path("data/raw", dataset),
    n_files = 0L,
    has_matrix = FALSE,
    has_image = FALSE,
    has_coordinates = FALSE,
    status = "not_detected"
  )]
}

inventory <- merge(expected, detection, by = "dataset", all.x = TRUE, suffixes = c("_expected", "_detected"))
inventory[, available := fifelse(
  expected_file_type == "expression_matrix", has_matrix,
  fifelse(expected_file_type == "histology_image", has_image,
          fifelse(expected_file_type == "spatial_coordinates", has_coordinates, FALSE))
)]
inventory[, inventory_status := fifelse(available, "available", "missing")]
inventory[, notes := fifelse(
  available,
  "Detected a compatible local file category; manual format inspection still required before analysis.",
  "Required spatial file category not detected locally; no spatial analysis should be interpreted."
)]
inventory[, `:=`(
  local_root = fifelse(!is.na(local_root_detected), local_root_detected, local_root_expected),
  expected_status = status_expected,
  detected_status = status_detected
)]
setcolorder(inventory, c(
  "dataset", "expected_file_type", "required", "local_root", "accepted_patterns",
  "n_files", "has_matrix", "has_image", "has_coordinates", "available",
  "inventory_status", "expected_status", "detected_status", "notes"
))
inventory[, c("local_root_expected", "local_root_detected", "status_expected", "status_detected") := NULL]

fwrite(inventory, out_path, sep = "\t")
message("wrote\t", out_path)
