suppressPackageStartupMessages({
  library(data.table)
})

dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

src <- "manuscript/manuscript_draft_v0.6.md"
stopifnot(file.exists(src))
x <- readLines(src, warn = FALSE)

phase8_start <- grep("^## Phase 8 targeted figure-enhancement note", x)
if (length(phase8_start) == 1) {
  phase8_note <- x[phase8_start:length(x)]
  writeLines(phase8_note, "docs/internal_phase8_note.md", useBytes = TRUE)
  x <- x[seq_len(phase8_start - 1)]
}

x <- gsub("locked kidney stone GWAS", "primary public KSD GWAS", x, fixed = TRUE)
x <- gsub("Locked kidney stone disease GWAS", "Primary public kidney stone disease GWAS", x, fixed = TRUE)
x <- gsub("locked kidney stone disease GWAS", "primary public kidney stone disease GWAS", x, fixed = TRUE)
x <- gsub("locked KSD GWAS", "primary KSD GWAS", x, fixed = TRUE)
x <- gsub("locked kidney stone GWAS reconstruction", "primary public KSD GWAS reconstruction", x, fixed = TRUE)

intro_hit <- grep("^Here we integrated", x)
innovation <- paste(
  "Unlike analyses that stop at locus-to-gene mapping, this study links",
  "MAGMA-prioritized KSD genes to audited renal papillary single-nucleus",
  "contexts and an independent plaque/stone papilla expression dataset, while",
  "explicitly separating supported cellular-context associations from untested",
  "causal mechanisms."
)
if (length(intro_hit) == 1 && !any(grepl("Unlike analyses that stop at locus-to-gene mapping", x, fixed = TRUE))) {
  x <- append(x, innovation, after = intro_hit - 1)
}

lim_start <- grep("^## Limitations", x)
if (length(lim_start) == 1) {
  next_section <- grep("^## ", x)
  next_section <- next_section[next_section > lim_start][1]
  if (is.na(next_section)) next_section <- length(x) + 1
  limitation_body <- x[(lim_start + 1):(next_section - 1)]
  limitation_body <- limitation_body[nzchar(limitation_body)]
  limitation_body <- limitation_body[!duplicated(limitation_body)]
  x <- c(x[seq_len(lim_start)],
         "",
         "TWAS, SMR/coloc and spatial transcriptomic validation were not completed because required external expression-weight, eQTL and spatial matrix/image resources were unavailable in the current analysis environment. GSE231569 supports single-nucleus expression-context mapping but does not provide spatial validation or causal cell-type mediation. GSE73680 supports disease-context module association but not causality or cell-type-specific disease response. P1 single-gene disease differential expression was not FDR-supported, and PKD2 should be treated only as a nominal exploratory observation.",
         "",
         x[next_section:length(x)])
}

writeLines(x, "manuscript/manuscript_draft_v0.7.md", useBytes = TRUE)

qc <- data.table(
  item = c("phase8_internal_note", "locked_language", "innovation_sentence", "limitations_deduplicated"),
  status = c(
    ifelse(file.exists("docs/internal_phase8_note.md"), "moved_to_docs", "not_present"),
    ifelse(any(grepl("locked", x, ignore.case = TRUE)), "review_remaining_locked_terms", "clean"),
    ifelse(any(grepl("Unlike analyses that stop at locus-to-gene mapping", x, fixed = TRUE)), "added", "missing"),
    "completed"
  )
)
fwrite(qc, "results/tables/phase9_manuscript_cleanup_qc_v0.1.tsv", sep = "\t")

message("wrote manuscript_draft_v0.7.md")
