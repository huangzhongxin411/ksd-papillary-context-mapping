suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

go <- fread("results/tables/magma_pathway_enrichment.tsv")
if (!all(c("gene_set", "Description", "p.adjust", "Count", "geneID") %in% names(go))) {
  stop("magma_pathway_enrichment.tsv lacks expected GO enrichment columns.")
}

background <- data.table(
  item = c("GO universe", "display_FDR_threshold", "display_min_gene_count", "term_redundancy_method", "main_figure_focus"),
  value = c("MAGMA-tested genes mapped to Entrez IDs",
            "FDR < 0.10",
            "Count >= 2",
            "greedy removal of terms with gene-overlap Jaccard >= 0.70 within each gene set",
            "nephron development, transport, mineral/ion handling and injury-context terms")
)
fwrite(background, "results/tables/pathway_enrichment_background_audit.tsv", sep = "\t")

go[, topic_priority := grepl("loop|Henle|tubule|nephron|kidney|renal|water|vitamin D|metal ion|phosphate|ion|transport|epithelial|urate", Description, ignore.case = TRUE)]
filtered <- go[p.adjust < 0.10 & Count >= 2]
filtered[, display_priority := fcase(
  topic_priority, "theme_relevant",
  grepl("MHC|antigen|immune", Description, ignore.case = TRUE), "immune_or_injury_context",
  default = "other_functional_context"
)]
fwrite(filtered, "results/tables/pathway_gene_count_filter.tsv", sep = "\t")

split_genes <- function(x) unique(unlist(strsplit(x, "/", fixed = TRUE)))
reduced <- rbindlist(lapply(unique(filtered$gene_set), function(gs) {
  d <- filtered[gene_set == gs][order(p.adjust, -Count)]
  keep <- rep(TRUE, nrow(d))
  kept_genes <- list()
  for (i in seq_len(nrow(d))) {
    gi <- split_genes(d$geneID[i])
    if (length(kept_genes) > 0) {
      jac <- vapply(kept_genes, function(gj) length(intersect(gi, gj)) / length(union(gi, gj)), numeric(1))
      if (any(jac >= 0.70, na.rm = TRUE)) keep[i] <- FALSE
    }
    if (keep[i]) kept_genes[[length(kept_genes) + 1]] <- gi
  }
  d[, redundancy_reduced_keep := keep]
  d
}), fill = TRUE)
fwrite(reduced, "results/tables/go_bp_redundancy_reduced_terms.tsv", sep = "\t")

plot_dt <- reduced[redundancy_reduced_keep == TRUE & display_priority != "other_functional_context"]
plot_dt <- plot_dt[order(display_priority, p.adjust)][, head(.SD, 8), by = gene_set]
plot_dt[, Description_short := ifelse(nchar(Description) > 48, paste0(substr(Description, 1, 45), "..."), Description)]
plot_dt[, Description_short := factor(Description_short, levels = rev(unique(Description_short)))]
plot_dt[, gene_set_label := factor(gene_set,
                                   levels = c("MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "Loop_TAL_influential", "P1_core"),
                                   labels = c("MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "Loop/TAL contributors", "P1 core"))]

p <- ggplot(plot_dt, aes(-log10(p.adjust), Description_short, size = Count, fill = gene_set_label)) +
  geom_point(shape = 21, color = "#555555", stroke = 0.22) +
  scale_fill_manual(values = c("MAGMA top 100" = "#3E6672", "MAGMA FDR" = "#6F8F98",
                               "MAGMA suggestive" = "#7F9AA3", "Loop/TAL contributors" = "#B59A5B",
                               "P1 core" = "#9A5F52")) +
  labs(title = "Redundancy-reduced GO BP terms",
       subtitle = "Filtered to FDR < 0.10 and gene count >= 2; theme-relevant and immune/injury-context terms retained",
       x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555555"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave("results/figures/figure5_go_redundancy_reduced_dotplot.pdf", p,
       width = 8.2, height = 6.2, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure5_go_redundancy_reduced_dotplot.png", p,
       width = 8.2, height = 6.2, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# Pathway Enrichment Audit v0.1",
  "",
  "GO enrichment is used for functional interpretation, not pathway-level validation.",
  "Displayed terms are filtered by FDR, gene count and redundancy to reduce overinterpretation from small or overlapping GO terms."
), "docs/pathway_enrichment_audit_v0.1.md", useBytes = TRUE)

message("wrote pathway enrichment audit outputs")
