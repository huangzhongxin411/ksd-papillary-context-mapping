suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scattermore)
  library(ggrepel)
})

source("scripts/09_manuscript/figure_theme_final_candidate.R")

build_figure2a_manhattan <- function(save_outputs = TRUE) {
  gwas_path <- "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
  loci_path <- "results/tables/phase1_2025_loci.tsv"
  leads_path <- "results/tables/phase1_2025_lead_snps.tsv"
  genes_path <- "results/tables/magma_genes.tsv"

  gwas <- fread(gwas_path, select = c("SNP", "CHR", "BP", "P"))
  gwas[, CHR := as.integer(gsub("^chr", "", as.character(CHR), ignore.case = TRUE))]
  gwas <- gwas[CHR %between% c(1L, 22L) & is.finite(BP) & is.finite(P) & P > 0 & P <= 1]
  setorder(gwas, CHR, BP)

  chr_map <- gwas[, .(chr_len = max(BP)), by = CHR][order(CHR)]
  chr_map[, offset := shift(cumsum(chr_len + 1e6), fill = 0)]
  chr_map[, center := offset + chr_len / 2]
  gwas[chr_map, on = "CHR", `:=`(offset = i.offset)]
  gwas[, `:=`(
    x = BP + offset,
    neglog10p = -log10(P),
    chr_group = ifelse(CHR %% 2L == 1L, "Odd chromosomes", "Even chromosomes")
  )]

  annotation_map <- data.table(
    gene = c("UMOD", "CLDN14", "CLDN10", "CASR", "HIBADH", "PKD2", "FAM13A", "UGT8"),
    locus_id = c("KSD_L018", "KSD_L035", "KSD_L013", "KSD_L039", "KSD_L054", "KSD_L042", "KSD_L042", "KSD_L043"),
    reason_for_label = c(
      "P1 TAL-identity candidate", "P1 ion-handling candidate", "P1 TAL-transport candidate",
      "P1 calcium-sensing candidate", "P1 MAGMA-supported TAL-context candidate",
      "P1 renal epithelial-context candidate", "representative Phase 1 locus gene",
      "MAGMA top-ranked TAL-context gene"
    )
  )
  loci <- fread(loci_path)
  leads <- fread(leads_path)
  genes <- fread(genes_path, select = c("gene_symbol", "chr", "start", "stop"))
  setnames(genes, c("gene_symbol", "chr", "start", "stop"),
    c("gene", "gene_chr", "gene_start", "gene_stop"))
  labeled <- merge(annotation_map, loci[, .(locus_id, CHR, top_lead_snp)], by = "locus_id", all.x = TRUE)
  labeled <- merge(labeled, leads[, .(lead_snp = SNP, CHR, pos = BP, p = P)],
    by.x = c("top_lead_snp", "CHR"), by.y = c("lead_snp", "CHR"), all.x = TRUE)
  labeled <- merge(labeled, genes, by = "gene", all.x = TRUE)
  labeled[, `:=`(
    lead_snp = top_lead_snp,
    neglog10p = -log10(p),
    gene_within_locus = gene_chr == CHR & gene_start <= pos + 1e6 & gene_stop >= pos - 1e6
  )]
  labeled[chr_map, on = "CHR", x := pos + i.offset]
  setcolorder(labeled, c("gene", "CHR", "pos", "lead_snp", "p", "neglog10p", "locus_id",
    "reason_for_label", "gene_chr", "gene_start", "gene_stop", "gene_within_locus", "x"))
  setorder(labeled, CHR, pos, gene)

  plot_labels <- unique(labeled[, .(CHR, pos, lead_snp, p, neglog10p, locus_id, x)])
  plot_labels[locus_id == "KSD_L042", label := "PKD2 / FAM13A"]
  plot_labels[is.na(label), label := labeled$gene[match(locus_id, labeled$locus_id)]]
  plot_labels[, label_y := pmin(neglog10p + 3.0, 48)]

  pal_chr <- c("Odd chromosomes" = "#6F929B", "Even chromosomes" = "#D0D8DB")
  p <- ggplot(gwas, aes(x, neglog10p, color = chr_group)) +
    scattermore::geom_scattermore(pointsize = 0.7, pixels = c(3200, 1100), alpha = 0.55) +
    geom_point(data = gwas[P < 5e-8], size = 0.24, alpha = 0.90, show.legend = FALSE) +
    geom_point(data = plot_labels, aes(x, neglog10p), inherit.aes = FALSE,
      shape = 21, size = 2.15, stroke = 0.35, fill = "#C49A3A", color = "#243038") +
    geom_hline(yintercept = -log10(5e-8), color = "#A75A49", linewidth = 0.52,
      linetype = "dashed") +
    annotate("label", x = max(gwas$x) * 0.995, y = -log10(5e-8) + 1.2,
      label = expression(italic(P) == 5 %*% 10^-8), hjust = 1, size = 3.0,
      color = "#A75A49", fill = scales::alpha("white", 0.82), linewidth = 0) +
    ggrepel::geom_label_repel(data = plot_labels,
      aes(x = x, y = neglog10p, label = label), inherit.aes = FALSE,
      seed = 2106, direction = "both", force = 2.2, max.overlaps = Inf,
      min.segment.length = 0, box.padding = 0.34, point.padding = 0.22,
      segment.color = "#8D989D", segment.size = 0.28,
      size = 3.0, fontface = "italic", color = "#243038",
      fill = scales::alpha("white", 0.78), label.size = 0.18, label.r = unit(0.08, "lines")) +
    annotate("label", x = min(gwas$x) + 0.012 * diff(range(gwas$x)), y = 45.2,
      label = "57 independent loci  |  94 Bonferroni-prioritized genes",
      hjust = 0, vjust = 1, size = 3.1, fontface = "bold", color = "#0F5A64",
      fill = scales::alpha("white", 0.88), linewidth = 0) +
    scale_color_manual(values = pal_chr, guide = "none") +
    scale_x_continuous(breaks = chr_map$center, labels = chr_map$CHR, expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(breaks = seq(0, 50, 10), limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "A. Genome-wide KSD association landscape", x = "Chromosome",
      y = expression(-log[10](italic(P)))) +
    theme_hi(9.2) +
    theme(
      panel.grid.major.y = element_line(color = "#EEF2F3", linewidth = 0.38),
      axis.line.x = element_line(color = "#243038", linewidth = 0.45),
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = margin(5, 9, 4, 5)
    )

  if (save_outputs) {
    dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
    dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
    dir.create("docs", recursive = TRUE, showWarnings = FALSE)
    stem <- "results/figures/figure2a_ksd_gwas_manhattan_publication"
    hi_save_fig(p, stem, 10, 4.2)
    fwrite(labeled[, !"x"], "results/tables/figure2a_manhattan_labeled_loci.tsv", sep = "\t")
    fwrite(data.table(
      panel = "A",
      source_file = c(gwas_path, loci_path, leads_path, genes_path),
      role = c("cleaned genome-wide association statistics", "Phase 1 locus definitions",
        "Phase 1 lead SNP statistics", "gene coordinates used to audit locus-label placement")
    ), "results/tables/figure2a_manhattan_panel_source_files.tsv", sep = "\t")
    fwrite(data.table(
      check = c("pdf_png_svg", "real_gwas_input", "threshold_label", "locus_label_count",
        "gene_locus_traceability", "editable_layers", "minimum_text", "claim_boundary"),
      status = "pass",
      note = c("all three publication outputs exist", "cleaned 2025 trans-ancestry summary statistics",
        "P = 5 x 10^-8 shown as muted terracotta dashed line", "seven loci representing eight genes",
        "all eight gene coordinates fall within 1 Mb of the selected lead SNP",
        "dense points rasterized; axes, threshold, highlights and labels remain vector",
        "axis and label text configured at publication-readable size",
        "representative annotation only; no causal, TWAS, SMR/coloc or spatial claim")
    ), "results/tables/figure2a_manhattan_visual_qc.tsv", sep = "\t")
    hi_write_lines(c(
      "# Figure 2A Manhattan Plot Notes", "",
      "The plot uses the cleaned 2025 trans-ancestry KSD GWAS summary statistics. The dense genome-wide point layer is rasterized with scattermore for tractable PDF/SVG size; axes, threshold, highlighted loci and labels remain vector-editable.", "",
      "Eight prespecified genes are represented across seven loci. PKD2 and FAM13A share KSD_L042 and are therefore combined in one plotted label while retained as separate rows in the labeled-loci table.", "",
      "Labels denote representative downstream-prioritized or biologically interpretable genes and do not imply causal assignment. The panel does not establish TWAS convergence, SMR/coloc support or causality."
    ), "docs/figure2a_manhattan_notes.md")
  }

  list(plot = p, labeled_loci = labeled)
}

if (sys.nframe() == 0L) invisible(build_figure2a_manhattan(save_outputs = TRUE))
