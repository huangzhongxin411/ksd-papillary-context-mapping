suppressPackageStartupMessages(library(data.table))

input <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input, select = c("CHR", "BP", "P"))
dt <- dt[!is.na(P) & P > 0 & P <= 1]

qq_p <- sort(dt$P)
expected <- -log10(ppoints(length(qq_p)))
observed <- -log10(qq_p)
pdf(file.path(fig_dir, "phase1_gwas_2025_qq_plot.pdf"), width = 5, height = 5)
plot(expected, observed, pch = 16, cex = 0.25, col = adjustcolor("black", 0.25),
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
     main = "KSD GWAS QQ plot")
abline(0, 1, col = "red", lwd = 1)
dev.off()

dt[, CHR := as.integer(gsub("^chr", "", as.character(CHR), ignore.case = TRUE))]
dt <- dt[!is.na(CHR) & CHR >= 1 & CHR <= 22]
setorder(dt, CHR, BP)
chr_len <- dt[, .(chr_max = max(BP, na.rm = TRUE)), by = CHR][order(CHR)]
chr_len[, offset := c(0, cumsum(chr_max[-.N] + 1e6))]
dt <- merge(dt, chr_len[, .(CHR, offset)], by = "CHR")
dt[, x := BP + offset]
dt[, y := -log10(P)]
axis_dt <- dt[, .(center = (min(x) + max(x)) / 2), by = CHR][order(CHR)]

png(file.path(fig_dir, "phase1_gwas_2025_manhattan_plot.png"), width = 2400, height = 900, res = 180)
cols <- c("#2f4858", "#4f7cac")
plot(dt$x, dt$y, pch = 16, cex = 0.18, col = cols[(dt$CHR %% 2) + 1],
     xaxt = "n", xlab = "Chromosome", ylab = "-log10(P)",
     main = "KSD trans-ancestry GWAS")
axis(1, at = axis_dt$center, labels = axis_dt$CHR, cex.axis = 0.65)
abline(h = -log10(5e-8), col = "red", lwd = 1)
dev.off()

cat("wrote\t", file.path(fig_dir, "phase1_gwas_2025_qq_plot.pdf"), "\n", sep = "")
cat("wrote\t", file.path(fig_dir, "phase1_gwas_2025_manhattan_plot.png"), "\n", sep = "")

