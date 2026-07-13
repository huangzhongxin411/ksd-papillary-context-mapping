#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- getwd()
dir.create(file.path(root, "config"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "notes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "codex_tasks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "data/raw/spatial/geo_metadata"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
rel <- function(path) sub(paste0("^", normalizePath(root, mustWork = TRUE), "/?"), "", normalizePath(path, mustWork = FALSE), fixed = TRUE)

manifest_path <- file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv")
section_audit_path <- file.path(root, "results/tables/phase3_step1_spatial_section_audit.tsv")
gene_overlap_path <- file.path(root, "results/tables/phase3_step1_spatial_scrna_gene_overlap.tsv")
roi_report_path <- file.path(root, "notes/phase3_step1_roi_annotation_search_report.md")
existing_config_path <- file.path(root, "config/spatial_sample_metadata.tsv")

manifest <- fread(manifest_path)
section_audit <- fread(section_audit_path)
gene_overlap <- fread(gene_overlap_path)
existing_config <- if (file.exists(existing_config_path)) fread(existing_config_path, fill = TRUE) else data.table()

soft_files <- c(
  GSE206306 = file.path(root, "data/raw/spatial/geo_metadata/GSE206306_family.soft.gz"),
  GSE231630 = file.path(root, "data/raw/spatial/geo_metadata/GSE231630_family.soft.gz")
)
soft_urls <- c(
  GSE206306 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206306/soft/GSE206306_family.soft.gz",
  GSE231630 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231630/soft/GSE231630_family.soft.gz"
)
for (nm in names(soft_files)) {
  if (!file.exists(soft_files[[nm]])) {
    tryCatch(
      utils::download.file(soft_urls[[nm]], soft_files[[nm]], mode = "wb", quiet = TRUE),
      error = function(e) message("GEO metadata download unavailable for ", nm, "; using local files only.")
    )
  }
}

parse_soft <- function(path, dataset_label) {
  if (!file.exists(path)) return(data.table())
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_starts <- grep("^\\^SAMPLE = ", lines)
  out <- list()
  for (i in seq_along(sample_starts)) {
    start <- sample_starts[i]
    end <- if (i < length(sample_starts)) sample_starts[i + 1] - 1 else length(lines)
    block <- lines[start:end]
    val <- function(prefix) {
      hit <- block[startsWith(block, prefix)]
      if (!length(hit)) return("")
      sub(prefix, "", hit[1], fixed = TRUE)
    }
    vals <- function(prefix) {
      hit <- block[startsWith(block, prefix)]
      if (!length(hit)) return(character())
      sub(prefix, "", hit, fixed = TRUE)
    }
    accession <- val("!Sample_geo_accession = ")
    characteristics <- vals("!Sample_characteristics_ch1 = ")
    supp <- vals("!Sample_supplementary_file")
    out[[length(out) + 1]] <- data.table(
      geo_accession = accession,
      dataset_label = dataset_label,
      sample_title = val("!Sample_title = "),
      source_name = val("!Sample_source_name_ch1 = "),
      characteristics_joined = paste(characteristics, collapse = " | "),
      supplementary_files = paste(supp, collapse = " | "),
      soft_file = rel(path)
    )
  }
  rbindlist(out, fill = TRUE)
}

soft_meta <- rbindlist(Map(parse_soft, soft_files, names(soft_files)), fill = TRUE)
soft_meta <- soft_meta[geo_accession %in% unique(manifest$geo_accession) |
                         geo_accession %in% c("GSM7166168", "GSM7166169")]

get_meta <- function(accession) {
  canonical <- sub("_inferred_from_archive$", "", accession)
  row <- soft_meta[geo_accession == canonical & dataset_label == "GSE206306"]
  if (!nrow(row)) row <- soft_meta[geo_accession == canonical]
  if (!nrow(row)) return(NULL)
  row[1]
}

extract_characteristic <- function(text, key) {
  parts <- trimws(unlist(strsplit(text, "\\|", fixed = FALSE)))
  hit <- parts[grepl(paste0("^", key, ":"), parts, ignore.case = TRUE)]
  if (!length(hit)) return("")
  trimws(sub("^[^:]+:\\s*", "", hit[1]))
}

section_subject <- function(section_id, title, description) {
  candidates <- unique(c(
    regmatches(section_id, gregexpr("KRP[0-9]+|F[0-9]+|M[0-9]+|20-[0-9]+|K[0-9]+", section_id, perl = TRUE))[[1]],
    regmatches(title, gregexpr("KRP[0-9]+|F[0-9]+|M[0-9]+|20-[0-9]+|K[0-9]+", title, perl = TRUE))[[1]]
  ))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) return(candidates[1])
  ""
}

curate_condition <- function(raw) {
  raw_l <- tolower(raw)
  if (grepl("reference", raw_l)) return(list(curated = "Reference", disease_status = "reference_control", stone_type = "none_reported", control = "control", confidence = "moderate"))
  if (grepl("calcium oxalate", raw_l)) return(list(curated = "Calcium Oxalate stone disease", disease_status = "stone disease", stone_type = "Calcium Oxalate", control = "stone_disease", confidence = "high"))
  if (grepl("brushite", raw_l)) return(list(curated = "Stone disease; mixed calcium oxalate/brushite metadata", disease_status = "stone disease", stone_type = "mixed_or_section_specific", control = "stone_disease", confidence = "moderate"))
  list(curated = "unclear", disease_status = "unclear", stone_type = "unknown", control = "unclear", confidence = "low")
}

curated_rows <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  m <- manifest[i]
  meta <- get_meta(m$geo_accession)
  if (is.null(meta)) {
    raw_condition <- ""
    title <- ""
    source <- ""
    description <- ""
    tissue_source <- ""
    metadata_source <- "No matching GEO SOFT sample metadata found; local manifest only."
    confidence <- "low"
  } else {
    raw_condition <- extract_characteristic(meta$characteristics_joined, "disease")
    title <- meta$sample_title
    description <- extract_characteristic(meta$characteristics_joined, "description")
    tissue_source <- extract_characteristic(meta$characteristics_joined, "tissue source")
    source <- meta$source_name
    metadata_source <- paste0(meta$soft_file, ":", meta$geo_accession)
    confidence <- "moderate"
  }
  cc <- curate_condition(raw_condition)
  confidence <- cc$confidence
  if (grepl("inferred_from_archive", m$geo_accession)) confidence <- ifelse(confidence == "high", "moderate", confidence)
  subj <- section_subject(m$section_id, title, description)
  has_image <- nzchar(m$hires_image_file) || nzchar(m$lowres_image_file)
  complete <- m$complete_visium_input == "complete"
  overlap <- gene_overlap[sample_id == m$sample_id & section_id == m$section_id]
  overlap_ok <- nrow(overlap) && overlap$suitable_for_rctd[1] == "yes"
  include <- ifelse(complete && overlap_ok, "yes", "no")
  exclusion <- ifelse(include == "yes", "", "Incomplete Visium input or insufficient gene overlap.")
  data.table(
    sample_id = m$sample_id,
    geo_accession = sub("_inferred_from_archive$", "", m$geo_accession),
    dataset_label = "GSE206306",
    section_id = m$section_id,
    original_section_name = m$section_id,
    donor_id_if_available = subj,
    subject_id_if_available = subj,
    condition_label_raw = raw_condition,
    condition_label_curated = cc$curated,
    disease_status = cc$disease_status,
    stone_type_if_available = cc$stone_type,
    control_or_disease = cc$control,
    anatomical_region = "kidney papilla",
    tissue_type = ifelse(nzchar(tissue_source), tissue_source, source),
    roi_annotation_available = "no",
    histology_image_available = ifelse(has_image, "yes", "no"),
    included_for_phase3_step2 = include,
    exclusion_reason = exclusion,
    metadata_confidence = confidence,
    metadata_source = metadata_source,
    notes = paste(na.omit(c(
      ifelse(grepl("inferred_from_archive", m$geo_accession), "GEO accession inferred from archive naming; sample-level GEO metadata supports parent sample but section mapping should be human-checked.", ""),
      ifelse(grepl("Two samples|Three samples", description), description, ""),
      "No ROI annotation available; include only for anatomical/tissue-context projection."
    )), collapse = " ")
  )
}), fill = TRUE)

write_tsv(curated_rows, file.path(root, "config/spatial_sample_metadata_curated_phase3.tsv"))

evidence_fields <- c(
  "condition_label_raw", "condition_label_curated", "control_or_disease",
  "stone_type_if_available", "anatomical_region", "roi_annotation_available",
  "histology_image_available", "included_for_phase3_step2"
)
evidence <- rbindlist(lapply(seq_len(nrow(curated_rows)), function(i) {
  row <- curated_rows[i]
  rbindlist(lapply(evidence_fields, function(field) {
    src <- row$metadata_source
    ev <- switch(field,
      condition_label_raw = "GEO SOFT !Sample_characteristics_ch1 disease field",
      condition_label_curated = "Curated from GEO disease field without inferring beyond stated terms",
      control_or_disease = "Mapped from curated condition; Reference -> control, Calcium Oxalate -> stone_disease",
      stone_type_if_available = "GEO disease field and/or sample description",
      anatomical_region = "GEO tissue source/source_name and project accession memo; kidney papilla context",
      roi_annotation_available = "Phase 3-Step 1 ROI report found no spot-level ROI table",
      histology_image_available = "Phase 3-Step 1 sample manifest hires/lowres image availability",
      included_for_phase3_step2 = "Complete Visium input plus sufficient spatial-snRNA gene overlap",
      ""
    )
    data.table(
      sample_id = row$sample_id,
      metadata_field = field,
      curated_value = as.character(row[[field]]),
      evidence_source_file_or_url = src,
      evidence_text_or_column = ev,
      confidence = row$metadata_confidence,
      notes = ifelse(field == "roi_annotation_available", "No plaque/mineral/lesion/fibrosis ROI annotation found; no ROI-level claims allowed.", "")
    )
  }))
}), fill = TRUE)
write_tsv(evidence, file.path(root, "results/tables/phase3_step1B_spatial_metadata_evidence.tsv"))

decision <- merge(manifest[, .(sample_id, section_id, complete_visium_input)],
                  curated_rows[, .(sample_id, section_id, metadata_confidence, roi_annotation_available, control_or_disease, included_for_phase3_step2)],
                  by = c("sample_id", "section_id"), all.x = TRUE)
decision <- merge(decision, gene_overlap[, .(sample_id, section_id, suitable_for_rctd)],
                  by = c("sample_id", "section_id"), all.x = TRUE)
decision[, metadata_locked := ifelse(metadata_confidence %in% c("high", "moderate"), "yes", "no")]
decision[, gene_overlap_sufficient := ifelse(suitable_for_rctd == "yes", "yes", "no")]
decision[, include_for_deconvolution := ifelse(complete_visium_input == "complete" & gene_overlap_sufficient == "yes", "yes", "no")]
decision[, include_for_module_scoring := ifelse(include_for_deconvolution == "yes", "yes", "no")]
decision[, include_for_group_comparison := ifelse(metadata_locked == "yes" & control_or_disease %in% c("control", "stone_disease"), "yes", "no")]
decision[, reason := fifelse(include_for_deconvolution == "yes",
                             "Complete Visium input and sufficient gene overlap; no ROI-specific claims.",
                             "Incomplete input or insufficient gene overlap.")]
decision[, notes := fifelse(include_for_group_comparison == "yes",
                            "Disease/control comparison is metadata-supported, but ROI/plaque-specific claims remain disallowed.",
                            "Use only for within-section anatomical/tissue-context projection unless metadata is fixed.")]
decision <- decision[, .(
  sample_id, section_id, complete_visium_input, metadata_locked,
  roi_annotation_available, gene_overlap_sufficient,
  include_for_deconvolution, include_for_module_scoring,
  include_for_group_comparison, reason, notes
)]
write_tsv(decision, file.path(root, "results/tables/phase3_step1B_section_inclusion_decision.tsv"))

claim_note <- c(
  "# Phase 3-Step 1B Spatial Claim Boundary",
  "",
  "- No plaque-, mineral-, lesion- or fibrosis-resolved spot-level ROI annotation is available in the locked spatial inputs.",
  "- Spatial deconvolution or label transfer can support anatomical/tissue-context projection.",
  "- Spatial analysis cannot support plaque-specific localization.",
  "- Spatial analysis cannot support lesion-stage enrichment.",
  "- Spatial analysis cannot validate a causal niche.",
  "- Disease/control comparison is allowed only for samples with high or moderate metadata confidence and must remain section/sample-contextual.",
  "",
  "Recommended wording:",
  "",
  "Because no plaque-, mineral-, lesion- or fibrosis-resolved spot-level ROI annotation was available, spatial transcriptomics was used as an anatomical tissue-context projection layer rather than a plaque-specific localization assay."
)
writeLines(claim_note, file.path(root, "notes/phase3_step1B_spatial_claim_boundary.md"))

pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)
rctd_available <- pkg_available("spacexr") || pkg_available("RCTD")
seurat_available <- pkg_available("Seurat")
included_sections <- curated_rows[included_for_phase3_step2 == "yes", .N]
group_sections <- decision[include_for_group_comparison == "yes", .N]
decision_option <- if (rctd_available) {
  "A. RCTD primary, Seurat label transfer fallback"
} else if (seurat_available) {
  "B. Seurat label transfer primary because RCTD environment unavailable"
} else {
  "D. human decision required"
}
method_note <- c(
  "# Phase 3-Step 1B Phase 3-Step 2 Method Decision",
  "",
  paste0("- RCTD/spacexr environment appears available: ", rctd_available, "."),
  paste0("- Seurat label transfer environment appears available: ", seurat_available, "."),
  paste0("- Selected primary method: ", ifelse(rctd_available, "RCTD", "Seurat label transfer fallback"), "."),
  paste0("- Selected fallback method: ", ifelse(seurat_available, "Seurat label transfer", "human environment setup required"), "."),
  paste0("- Sections to include for deconvolution/label transfer: ", included_sections, " of ", nrow(curated_rows), "."),
  paste0("- Disease/control comparison allowed: ", ifelse(group_sections == nrow(curated_rows), "yes, for all included sections with high/moderate metadata confidence", "partially or no; see inclusion table"), "."),
  "- ROI/plaque-specific localization is not allowed.",
  "",
  paste0("Decision: **", decision_option, "**."),
  "",
  "If RCTD is installed before Phase 3-Step 2, use RCTD as primary. Otherwise proceed with Seurat label transfer as the conservative and available fallback."
)
writeLines(method_note, file.path(root, "notes/phase3_step1B_phase3_step2_method_decision.md"))

conf_counts <- curated_rows[, .N, by = metadata_confidence]
conf_summary <- paste(paste(conf_counts$metadata_confidence, conf_counts$N, sep = "="), collapse = "; ")
report <- c(
  "# Phase 3-Step 1B Report",
  "",
  "## Metadata Files Inspected",
  "",
  paste0("- `", rel(manifest_path), "`"),
  paste0("- `", rel(section_audit_path), "`"),
  paste0("- `", rel(gene_overlap_path), "`"),
  paste0("- `", rel(existing_config_path), "`"),
  paste0("- `", rel(soft_files[["GSE206306"]]), "`"),
  paste0("- `", rel(soft_files[["GSE231630"]]), "`"),
  "",
  "## Metadata Locking Outcome",
  "",
  "- Curated metadata table created: `config/spatial_sample_metadata_curated_phase3.tsv`.",
  paste0("- Metadata confidence counts: ", conf_summary, "."),
  paste0("- Sections included for deconvolution/label transfer: ", included_sections, " of ", nrow(curated_rows), "."),
  paste0("- Sections eligible for disease/control comparison: ", group_sections, " of ", nrow(curated_rows), "."),
  "- ROI annotation status: no plaque/mineral/lesion/fibrosis ROI annotation available.",
  "",
  "## Final Claim Boundary",
  "",
  "Spatial analysis can support anatomical/tissue-context projection only. It cannot support plaque-specific localization, lesion-stage enrichment, or causal niche validation.",
  "",
  "## Recommended Next Action",
  "",
  paste0("Recommendation: **", ifelse(decision_option == "A. RCTD primary, Seurat label transfer fallback", "A. proceed to Phase 3-Step 2: RCTD / label transfer on included sections", "B. proceed with Seurat label transfer only"), "**.")
)
writeLines(report, file.path(root, "notes/phase3_step1B_report.md"))

checklist <- data.table(
  task_id = sprintf("P3S1B-%02d", 1:8),
  task_name = c(
    "Create spatial metadata locking script",
    "Create curated spatial sample metadata table",
    "Create metadata evidence table",
    "Create section inclusion decision table",
    "Create no-ROI claim boundary note",
    "Create Phase 3-Step 2 method decision",
    "Create Phase 3-Step 1B report",
    "Respect stop rule and analysis boundary"
  ),
  completed = c(
    file.exists(file.path(root, "scripts/08_spatial_analysis/phase3_step1B_spatial_metadata_locking.R")),
    file.exists(file.path(root, "config/spatial_sample_metadata_curated_phase3.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1B_spatial_metadata_evidence.tsv")),
    file.exists(file.path(root, "results/tables/phase3_step1B_section_inclusion_decision.tsv")),
    file.exists(file.path(root, "notes/phase3_step1B_spatial_claim_boundary.md")),
    file.exists(file.path(root, "notes/phase3_step1B_phase3_step2_method_decision.md")),
    file.exists(file.path(root, "notes/phase3_step1B_report.md")),
    TRUE
  ),
  output_file = c(
    "scripts/08_spatial_analysis/phase3_step1B_spatial_metadata_locking.R",
    "config/spatial_sample_metadata_curated_phase3.tsv",
    "results/tables/phase3_step1B_spatial_metadata_evidence.tsv",
    "results/tables/phase3_step1B_section_inclusion_decision.tsv",
    "notes/phase3_step1B_spatial_claim_boundary.md",
    "notes/phase3_step1B_phase3_step2_method_decision.md",
    "notes/phase3_step1B_report.md",
    "No RCTD, Seurat label transfer, Cell2location, spatial scoring, TWAS, bulk, or manuscript edit performed"
  ),
  blocking_issue = c("none", "none", "none", "none", "no ROI annotation", ifelse(rctd_available, "none", "RCTD environment unavailable"), "none", "none"),
  manual_review_needed = c("no", "yes", "yes", "yes", "yes", "yes", "yes", "no"),
  notes = c(
    "Script parses downloaded/local GEO SOFT and prior Step 1 inputs.",
    "No disease labels fabricated; all labels trace to GEO fields.",
    "Evidence table links curated fields to source fields/files.",
    "Group comparison allowed only where metadata confidence is high/moderate.",
    "Plaque-specific localization remains disallowed.",
    "Seurat is available; RCTD availability checked without running it.",
    "Prepared for human review before Phase 3-Step 2.",
    "Stop rule respected."
  )
)
write_tsv(checklist, file.path(root, "codex_tasks/phase3_step1B_completion_checklist.tsv"))

message("Phase 3-Step 1B spatial metadata locking complete. No deconvolution, spatial scoring, TWAS, bulk, or manuscript edit was run.")
