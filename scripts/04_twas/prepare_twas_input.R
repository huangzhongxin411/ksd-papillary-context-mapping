suppressPackageStartupMessages({
  library(data.table)
})

infile <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
out_dir <- "data/processed/twas_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("SNP", "CHR", "BP", "EA", "NEA", "P", "BETA", "SE", "EAF", "N")
dt <- fread(infile)
missing_cols <- setdiff(required, names(dt))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

raw_n <- nrow(dt)
dt <- dt[, ..required]

dt[, `:=`(
  SNP = as.character(SNP),
  CHR = as.integer(CHR),
  BP = as.integer(BP),
  A1 = toupper(as.character(EA)),
  A2 = toupper(as.character(NEA)),
  P = as.numeric(P),
  BETA = as.numeric(BETA),
  SE = as.numeric(SE),
  EAF = as.numeric(EAF),
  N = as.numeric(N)
)]

dt[, EA := NULL]
dt[, NEA := NULL]

dt[, missing_required := is.na(SNP) | is.na(CHR) | is.na(BP) | is.na(A1) | is.na(A2) |
  is.na(P) | is.na(BETA) | is.na(SE) | is.na(N)]
dt[, invalid_p := is.na(P) | P <= 0 | P > 1]
dt[, invalid_se := is.na(SE) | SE <= 0]
dt[, invalid_eaf := !is.na(EAF) & (EAF < 0 | EAF > 1)]
dt[, duplicated_snp := duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]
dt[, palindromic := paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")]
dt[, non_acgt_allele := !(A1 %in% c("A", "C", "G", "T")) | !(A2 %in% c("A", "C", "G", "T"))]
dt[, allele_same := A1 == A2]

flags <- dt[, .(
  SNP, CHR, BP, A1, A2, EAF,
  missing_required,
  invalid_p,
  invalid_se,
  invalid_eaf,
  duplicated_snp,
  palindromic,
  non_acgt_allele,
  allele_same
)]

filtered <- dt[
  missing_required == FALSE &
    invalid_p == FALSE &
    invalid_se == FALSE &
    invalid_eaf == FALSE &
    duplicated_snp == FALSE &
    non_acgt_allele == FALSE &
    allele_same == FALSE
]

filtered[, Z := BETA / SE]
filtered <- filtered[, .(SNP, CHR, BP, A1, A2, BETA, SE, Z, P, N, EAF)]
setorder(filtered, CHR, BP, SNP)

qc <- data.table(
  metric = c(
    "source_file",
    "raw_rows",
    "output_rows",
    "missing_required_rows",
    "invalid_p_rows",
    "invalid_se_rows",
    "invalid_eaf_rows",
    "duplicate_snp_rows",
    "non_acgt_allele_rows",
    "allele_same_rows",
    "palindromic_rows_retained",
    "min_n",
    "median_n",
    "max_n",
    "min_p",
    "max_abs_z",
    "a1_definition",
    "palindromic_policy"
  ),
  value = c(
    infile,
    as.character(raw_n),
    as.character(nrow(filtered)),
    as.character(sum(dt$missing_required, na.rm = TRUE)),
    as.character(sum(dt$invalid_p, na.rm = TRUE)),
    as.character(sum(dt$invalid_se, na.rm = TRUE)),
    as.character(sum(dt$invalid_eaf, na.rm = TRUE)),
    as.character(sum(dt$duplicated_snp, na.rm = TRUE)),
    as.character(sum(dt$non_acgt_allele, na.rm = TRUE)),
    as.character(sum(dt$allele_same, na.rm = TRUE)),
    as.character(sum(filtered[, paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")], na.rm = TRUE)),
    as.character(min(filtered$N, na.rm = TRUE)),
    as.character(stats::median(filtered$N, na.rm = TRUE)),
    as.character(max(filtered$N, na.rm = TRUE)),
    format(min(filtered$P, na.rm = TRUE), scientific = TRUE),
    format(max(abs(filtered$Z), na.rm = TRUE), scientific = TRUE),
    "A1 is the GWAS effect allele, copied from EA; A2 is copied from NEA.",
    "Palindromic SNPs are retained in the standardized input and flagged separately for weight-level harmonization."
  )
)

fwrite(filtered, file.path(out_dir, "ksd_2025_for_twas.tsv.gz"), sep = "\t")
fwrite(flags, file.path(out_dir, "ksd_2025_for_twas_variant_flags.tsv.gz"), sep = "\t")
fwrite(qc, file.path("results", "tables", "twas_input_qc_report.tsv"), sep = "\t")

message("wrote\t", file.path(out_dir, "ksd_2025_for_twas.tsv.gz"))
message("wrote\t", file.path(out_dir, "ksd_2025_for_twas_variant_flags.tsv.gz"))
message("wrote\t", file.path("results", "tables", "twas_input_qc_report.tsv"))
