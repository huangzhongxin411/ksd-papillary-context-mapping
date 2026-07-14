#!/usr/bin/env Rscript

source_file <- "source_data/figures/phase1_step4_gwas_plot_source_data.tsv.gz"
out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
  stop("Missing plot source data: run scripts/01_gwas_qc/phase1_step4_gwas_qc_audit.py first")
}

d <- read.delim(gzfile(source_file), stringsAsFactors = FALSE)
d <- d[is.finite(d$P) & d$P > 0 & d$P <= 1 & is.finite(d$CHR) & is.finite(d$BP), ]
d$CHR <- as.integer(d$CHR)
d$BP <- as.integer(d$BP)
d$neglog10p <- -log10(d$P)

pdf(file.path(out_dir, "phase1_step4_gwas_pvalue_distribution.pdf"), width = 7, height = 5)
hist(d$P, breaks = 80, col = "grey75", border = "white",
     main = "GWAS P-value distribution (QC sample)",
     xlab = "P value")
mtext("Diagnostic plot from cleaned GWAS downsample plus genome-wide significant variants", side = 3, line = 0.2, cex = 0.75)
dev.off()

qq_p <- sort(d$P)
n <- length(qq_p)
expected <- -log10(ppoints(n))
observed <- -log10(qq_p)
pdf(file.path(out_dir, "phase1_step4_gwas_qq_plot.pdf"), width = 5.5, height = 5.5)
plot(expected, observed, pch = 16, cex = 0.35, col = rgb(0, 0, 0, 0.35),
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
     main = "GWAS QQ plot (QC sample)")
abline(0, 1, col = "firebrick", lwd = 1.5)
mtext("Visualization only; no biological inference", side = 3, line = 0.2, cex = 0.75)
dev.off()

chr_sizes <- aggregate(BP ~ CHR, d, max)
chr_sizes <- chr_sizes[order(chr_sizes$CHR), ]
offsets <- c(0, cumsum(as.numeric(chr_sizes$BP))[-nrow(chr_sizes)])
names(offsets) <- chr_sizes$CHR
d$cum_pos <- d$BP + offsets[as.character(d$CHR)]
axis_pos <- tapply(d$cum_pos, d$CHR, function(x) mean(range(x)))
palette <- rep(c("#334155", "#0f766e"), length.out = length(unique(d$CHR)))
cols <- palette[match(d$CHR, sort(unique(d$CHR)))]

pdf(file.path(out_dir, "phase1_step4_gwas_manhattan_plot.pdf"), width = 10, height = 4.8)
plot(d$cum_pos, d$neglog10p, pch = 16, cex = 0.28, col = adjustcolor(cols, alpha.f = 0.55),
     xaxt = "n", xlab = "Chromosome", ylab = "-log10(P)",
     main = "GWAS Manhattan plot (QC sample)")
axis(1, at = axis_pos, labels = names(axis_pos), cex.axis = 0.7)
abline(h = -log10(5e-8), col = "firebrick", lty = 2, lwd = 1.2)
mtext("Downsampled cleaned GWAS plus all sampled genome-wide significant variants; QC visualization only", side = 3, line = 0.2, cex = 0.75)
dev.off()
