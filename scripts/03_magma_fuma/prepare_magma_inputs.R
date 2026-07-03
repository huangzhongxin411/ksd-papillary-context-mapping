suppressPackageStartupMessages(library(data.table))

clean_gwas <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
out_dir <- "data/processed/magma/2025_trans_ancestry"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(clean_gwas, select = c("SNP", "CHR", "BP", "P", "N"))
dt <- dt[!is.na(SNP) & !is.na(CHR) & !is.na(BP) & !is.na(P)]
dt <- unique(dt, by = "SNP")

snploc <- dt[, .(SNP, CHR, BP)]
pval <- dt[, .(SNP, P, N)]

fwrite(snploc, file.path(out_dir, "ksd_2025.snploc"), sep = "\t", col.names = FALSE)
fwrite(pval, file.path(out_dir, "ksd_2025.pval"), sep = "\t")

cat("wrote\t", file.path(out_dir, "ksd_2025.snploc"), "\n", sep = "")
cat("wrote\t", file.path(out_dir, "ksd_2025.pval"), "\n", sep = "")
cat("n_snps\t", nrow(dt), "\n", sep = "")

