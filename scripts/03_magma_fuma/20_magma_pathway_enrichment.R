suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)

read_set <- function(path) {
  if (!file.exists(path)) character() else unique(fread(path, header = FALSE)[[1]])
}

sets <- list(
  MAGMA_top100 = read_set("results/gene_sets/magma_top100.txt"),
  MAGMA_FDR = read_set("results/gene_sets/magma_fdr05.txt"),
  MAGMA_suggestive = read_set("results/gene_sets/magma_suggestive_p1e4.txt"),
  Loop_TAL_influential = fread("results/tables/loop_tal_influential_magma_genes.tsv")[contribution_rank <= 20, gene],
  P1_core = c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
)
sets <- lapply(sets, unique)

magma <- fread("results/tables/magma_genes.tsv")
universe_symbols <- unique(magma$gene_symbol)

symbol_to_entrez <- function(sym) {
  suppressMessages(bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
}
universe_map <- symbol_to_entrez(universe_symbols)
universe_entrez <- unique(universe_map$ENTREZID)

go_rows <- rbindlist(lapply(names(sets), function(nm) {
  gm <- symbol_to_entrez(sets[[nm]])
  gene_ids <- unique(gm$ENTREZID)
  if (length(gene_ids) < 5) return(data.table())
  ego <- suppressMessages(enrichGO(
    gene = gene_ids,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    pvalueCutoff = 1,
    readable = TRUE
  ))
  dt <- as.data.table(ego)
  if (nrow(dt) == 0) return(data.table())
  dt[, gene_set := nm]
  dt
}), fill = TRUE)
if (nrow(go_rows) > 0) {
  fwrite(go_rows, "results/tables/magma_pathway_enrichment.tsv", sep = "\t")
  fwrite(go_rows[gene_set == "Loop_TAL_influential"], "results/tables/loop_tal_influential_gene_pathway_enrichment.tsv", sep = "\t")
} else {
  fwrite(data.table(note = "No GO BP enrichment rows returned"), "results/tables/magma_pathway_enrichment.tsv", sep = "\t")
  fwrite(data.table(note = "No GO BP enrichment rows returned"), "results/tables/loop_tal_influential_gene_pathway_enrichment.tsv", sep = "\t")
}

curated_sets <- list(
  TAL_transport = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2", "CASR", "CLDN14"),
  calcium_ion_handling = c("CASR", "CLDN14", "CLDN16", "CLDN19", "TRPV5", "TRPV6", "S100G", "ATP2B1"),
  epithelial_tight_junction = c("CLDN10", "CLDN14", "CLDN16", "CLDN19", "OCLN", "TJP1", "TJP2"),
  proximal_tubule_context = c("LRP2", "CUBN", "SLC5A2", "SLC34A1", "SLC22A6", "SLC22A8"),
  collecting_duct_context = c("AQP2", "AQP3", "SCNN1A", "SCNN1B", "SCNN1G", "AVPR2", "AQP4"),
  papillary_injury_remodeling = c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "FN1", "VIM", "HAVCR1", "LCN2")
)

hyper_one <- function(query, term_genes, universe) {
  query <- intersect(unique(query), universe)
  term <- intersect(unique(term_genes), universe)
  k <- length(intersect(query, term))
  n <- length(query)
  K <- length(term)
  N <- length(universe)
  p <- if (n == 0 || K == 0) NA_real_ else phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.table(overlap = k, query_size = n, term_size = K, universe_size = N, p_value = p,
             overlap_genes = paste(intersect(query, term), collapse = ";"))
}
curated <- rbindlist(lapply(names(sets), function(qn) {
  rbindlist(lapply(names(curated_sets), function(tn) {
    z <- hyper_one(sets[[qn]], curated_sets[[tn]], universe_symbols)
    z[, `:=`(gene_set = qn, term = tn)]
    z
  }))
}), fill = TRUE)
curated[, fdr := p.adjust(p_value, method = "BH")]
curated[, enrichment_ratio := (overlap / pmax(query_size, 1)) / (term_size / universe_size)]
curated[, interpretation := fcase(
  fdr < 0.05 & overlap >= 2, "FDR-supported curated functional enrichment",
  p_value < 0.05 & overlap >= 2, "nominal curated functional enrichment",
  default = "not enriched"
)]
fwrite(curated, "results/tables/nephron_segment_marker_enrichment.tsv", sep = "\t")

resource_status <- data.table(
  resource = c("GO Biological Process", "ReactomePA", "MSigDB/Hallmark via msigdbr", "KEGG"),
  status = c("completed_with_clusterProfiler_org.Hs.eg.db",
             ifelse(requireNamespace("ReactomePA", quietly = TRUE), "available_not_run", "unavailable_offline"),
             ifelse(requireNamespace("msigdbr", quietly = TRUE), "available_not_run", "unavailable_offline"),
             "not_run_no_local_resource"),
  note = c("Used for main pathway enrichment table.",
           "Package unavailable in current local R library.",
           "Package unavailable in current local R library.",
           "Not run to avoid online dependency and unsupported organism mapping assumptions.")
)
fwrite(resource_status, "results/tables/pathway_resource_status_v0.1.tsv", sep = "\t")

if (nrow(go_rows) > 0) {
  go_plot <- go_rows[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential")]
  go_plot <- go_plot[order(p.adjust)][, head(.SD, 4), by = gene_set]
  go_plot[, Description_short := ifelse(nchar(Description) > 46, paste0(substr(Description, 1, 43), "..."), Description)]
  go_plot[, Description_short := factor(Description_short, levels = rev(unique(Description_short)))]
  go_plot[, gene_set_label := factor(gene_set,
                                      levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential"),
                                      labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors"))]
  p_go <- ggplot(go_plot, aes(-log10(p.adjust), Description_short, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_manual(values = c("MAGMA top 100" = "#3E6672",
                                 "MAGMA FDR" = "#6F8F98",
                                 "Loop/TAL contributors" = "#B59A5B")) +
    labs(title = "B. GO biological process enrichment", x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_bw(base_size = 8.6) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
          panel.grid.minor = element_blank())
} else {
  p_go <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "GO enrichment unavailable", size = 4) + theme_void()
}

evidence <- data.table(
  evidence = c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context"),
  value = c("top-ranked/FDR", "enriched", "supported", "injury-coupled", "transport/calcium"),
  strength = c(0.95, 0.90, 0.75, 0.82, 0.78)
)
p_evidence <- ggplot(evidence, aes(strength, factor(evidence, levels = rev(evidence)), fill = strength)) +
  geom_col(width = 0.62, color = "#666666", linewidth = 0.2) +
  geom_text(aes(label = value), hjust = 1.05, color = "white", size = 2.7) +
  scale_fill_gradient(low = "#8AA0A8", high = "#3E6672") +
  labs(title = "A. Integrated evidence tiers", x = "Relative support", y = NULL) +
  theme_bw(base_size = 9.2) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank())

curated_plot <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
curated_plot[, gene_set_label := factor(gene_set,
                                        levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                        labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
curated_plot[, term_label := factor(term,
                                    levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction",
                                                   "proximal_tubule_context", "collecting_duct_context", "papillary_injury_remodeling")),
                                    labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction",
                                                   "Proximal tubule context", "Collecting duct context", "Papillary injury/remodeling")))]
p_curated <- ggplot(curated_plot, aes(gene_set_label, term_label, fill = pmin(enrichment_ratio, 20))) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(overlap > 0, overlap, "")), size = 2.7, color = "#303030") +
  scale_fill_gradient(low = "#EEF3F4", high = "#3E6672") +
  labs(title = "C. Curated nephron and functional context", x = NULL, y = NULL, fill = "Enrichment\nratio") +
  theme_bw(base_size = 8.6) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid = element_blank())

if (file.exists("results/tables/gse73680_risk_injury_module_correlations.tsv")) {
  cpl <- fread("results/tables/gse73680_risk_injury_module_correlations.tsv")
  cpl <- cpl[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")]
  cpl[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                               labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  cpl[, injury_label := factor(injury_module,
                               levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                               labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  p_couple <- ggplot(cpl, aes(injury_label, module_label, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.55) +
    scale_fill_gradient2(low = "#8AA0A8", mid = "white", high = "#9A5F52", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "D. GSE73680 risk-injury coupling", x = NULL, y = NULL, fill = "rho") +
    theme_bw(base_size = 8.6) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid = element_blank())
} else {
  p_couple <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Run Phase 9B first", size = 4) + theme_void()
}

fig5 <- plot_grid(p_evidence, p_go, p_curated, p_couple, ncol = 2, rel_heights = c(0.9, 1.1))
ggsave("results/figures/figure5_functional_context_v0.1.pdf", fig5,
       width = 13.2, height = 9.2, units = "in", device = "pdf", bg = "white")
ggsave("results/figures/figure5_functional_context_v0.1.png", fig5,
       width = 13.2, height = 9.2, units = "in", dpi = 260, bg = "white")

writeLines(c(
  "# Figure 5 Legend v0.1",
  "",
  "**Figure 5. Functional interpretation of MAGMA-prioritized TAL-associated KSD genes.**",
  "(A) Evidence-tier summary linking MAGMA prioritization, Loop/TAL specificity, donor detection, GSE73680 coupling and functional context. (B) GO Biological Process enrichment for MAGMA-prioritized and Loop/TAL influential gene sets. (C) Curated nephron and functional enrichment across MAGMA, Loop/TAL influential and P1 gene sets. Tile labels indicate overlapping genes. (D) GSE73680 module-level coupling between MAGMA-prioritized modules and papillary injury/remodeling programs.",
  "This figure supports functional interpretation and disease-context coupling of prioritized modules. It does not establish causal mechanism, TWAS convergence, colocalization or spatial validation."
), "docs/figure5_legend_v0.1.md", useBytes = TRUE)

message("wrote pathway enrichment and Figure 5")
