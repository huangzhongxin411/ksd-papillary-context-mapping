suppressPackageStartupMessages(library(data.table))

clean_gwas <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
raw_gwas <- "data/raw/gwas/2025_trans_ancestry/meta_sumstats"
qc_report_path <- "results/tables/phase1_gwas_qc_report.tsv"
leads_path <- "results/tables/phase1_2025_lead_snps.tsv"
loci_path <- "results/tables/phase1_2025_loci.tsv"
candidates_path <- "results/tables/phase1_candidate_genes.tsv"
fig_dir <- "results/figures"
table_dir <- "results/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

qc <- fread(qc_report_path)
qc_value <- function(metric) as.numeric(qc[metric == !!metric, value][1])

dt <- fread(clean_gwas, select = c("SNP", "CHR", "BP", "P"))
dt <- dt[!is.na(P) & P > 0 & P <= 1]
leads <- fread(leads_path)
loci <- fread(loci_path)
cand <- fread(candidates_path)

lambda_from_p <- function(p) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  q <- qchisq(p, df = 1, lower.tail = FALSE)
  median(q, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
}

plot_qq <- function(p, path, title) {
  p <- sort(p[!is.na(p) & p > 0 & p <= 1])
  expected <- -log10(ppoints(length(p)))
  observed <- -log10(p)
  pdf(path, width = 5, height = 5)
  plot(expected, observed, pch = 16, cex = 0.25, col = adjustcolor("black", 0.25),
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = title)
  abline(0, 1, col = "red", lwd = 1)
  dev.off()
}

dt[, CHR_clean := as.character(CHR)]
dt[, CHR_clean := gsub("^chr", "", CHR_clean, ignore.case = TRUE)]
dt[, CHR_clean := fifelse(CHR_clean %in% c("X", "Y"), CHR_clean, as.character(as.integer(CHR_clean)))]

no_gws <- dt[P >= 5e-8]
plot_qq(no_gws$P, file.path(fig_dir, "phase1_gwas_qq_no_gws.pdf"), "KSD GWAS QQ, excluding P < 5e-8")

mask <- rep(FALSE, nrow(dt))
for (i in seq_len(nrow(loci))) {
  chr <- as.character(loci$CHR[i])
  start <- as.integer(loci$start[i])
  end <- as.integer(loci$end[i])
  mask <- mask | (dt$CHR_clean == chr & dt$BP >= start & dt$BP <= end)
}
outside_loci <- dt[!mask]
plot_qq(outside_loci$P, file.path(fig_dir, "phase1_gwas_qq_without_known_loci.pdf"), "KSD GWAS QQ, outside Phase 1 loci")

candidate_counts <- cand[, .(n_candidate_genes = uniqueN(gene)), by = evidence_level]
get_count <- function(level) {
  val <- candidate_counts[evidence_level == level, n_candidate_genes]
  if (length(val) == 0) 0 else as.integer(val[1])
}

sanity <- data.table(
  metric = c(
    "raw_snp_n",
    "clean_snp_n",
    "genome_wide_sig_snp_n",
    "lead_snp_n",
    "independent_locus_n",
    "published_2025_locus_n",
    "min_p",
    "lambda_gc",
    "lambda_gc_no_gws",
    "lambda_gc_without_phase1_loci",
    "n_chr_with_gws_hits",
    "n_candidate_genes_50kb",
    "n_candidate_genes_250kb",
    "n_candidate_genes_500kb"
  ),
  value = as.character(c(
    qc[metric == "raw_rows", value][1],
    qc[metric == "qc_pass_rows", value][1],
    qc[metric == "genome_wide_significant_p_lt_5e-8", value][1],
    nrow(leads),
    nrow(loci),
    59,
    signif(min(dt$P, na.rm = TRUE), 6),
    signif(as.numeric(qc[metric == "lambda_gc", value][1]), 6),
    signif(lambda_from_p(no_gws$P), 6),
    signif(lambda_from_p(outside_loci$P), 6),
    uniqueN(dt[P < 5e-8, CHR_clean]),
    get_count("within_50kb"),
    get_count("within_250kb"),
    get_count("within_500kb")
  )),
  note = c(
    "Rows in downloaded 2025 trans-ancestry GWAS summary statistics.",
    "Rows after Phase 1 SNP QC.",
    "Significant SNPs after QC.",
    "Distance-pruned lead SNPs at P < 5e-8, +/-1 Mb.",
    "Distance-merged loci from Phase 1 lead SNPs.",
    "Reported loci in Cao et al. 2025 Nat Commun.",
    "Minimum P value after QC.",
    "Lambda GC after QC.",
    "Lambda GC after excluding all P < 5e-8 SNPs.",
    "Lambda GC after excluding Phase 1 lead loci windows.",
    "Chromosomes with genome-wide significant SNPs.",
    "Unique genes within 50 kb of lead SNPs.",
    "Unique genes within 250 kb of lead SNPs.",
    "Unique genes within 500 kb of lead SNPs."
  )
)
fwrite(sanity, file.path(table_dir, "phase1_gwas_sanity_check.tsv"), sep = "\t")

cat("wrote\t", file.path(table_dir, "phase1_gwas_sanity_check.tsv"), "\n", sep = "")
cat("wrote\t", file.path(fig_dir, "phase1_gwas_qq_no_gws.pdf"), "\n", sep = "")
cat("wrote\t", file.path(fig_dir, "phase1_gwas_qq_without_known_loci.pdf"), "\n", sep = "")

