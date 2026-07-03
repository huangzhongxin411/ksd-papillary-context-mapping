suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

obj_path <- "results/scrna/gse231569/objects/gse231569_annotated_audited.rds"
annot_path <- "config/gse231569_cluster_to_celltype_audited.tsv"

p1_genes <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
tal_celltype <- "Loop_of_Henle_TAL"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(obj_path)
annot <- fread(annot_path)
annot[, cluster := as.character(cluster)]

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, seurat_clusters := as.character(seurat_clusters)]
meta <- merge(
  meta,
  annot[, .(
    seurat_clusters = as.character(cluster),
    audited_broad_cell_type = broad_cell_type,
    audited_fine_cell_type = fine_cell_type,
    annotation_confidence,
    immune_review_flag,
    mural_review_flag
  )],
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE
)
setkey(meta, cell)
meta <- meta[colnames(obj)]

expr <- LayerData(obj, assay = DefaultAssay(obj), layer = "data")
present_genes <- intersect(p1_genes, rownames(expr))
missing_genes <- setdiff(p1_genes, rownames(expr))

p1_expr <- as.matrix(expr[present_genes, , drop = FALSE])

availability <- rbindlist(lapply(p1_genes, function(g) {
  if (g %in% rownames(p1_expr)) {
    v <- p1_expr[g, ]
    data.table(
      gene = g,
      present_in_object = TRUE,
      matched_symbol = g,
      n_cells_detected = sum(v > 0),
      global_pct_detected = mean(v > 0),
      global_avg_expression = mean(v),
      notes = "present"
    )
  } else {
    data.table(
      gene = g,
      present_in_object = FALSE,
      matched_symbol = NA_character_,
      n_cells_detected = 0L,
      global_pct_detected = 0,
      global_avg_expression = 0,
      notes = "not_detected_or_symbol_missing"
    )
  }
}))
fwrite(availability, "results/tables/p1_gene_availability.tsv", sep = "\t")

celltype_summary <- rbindlist(lapply(present_genes, function(g) {
  v <- p1_expr[g, ]
  tmp <- copy(meta)
  tmp[, gene := g]
  tmp[, expr := as.numeric(v[cell])]
  tmp[, .(
    n_cells = .N,
    n_donors = uniqueN(donor_id),
    avg_expression = mean(expr),
    median_expression = median(expr),
    pct_expressed = mean(expr > 0),
    annotation_confidence = paste(sort(unique(na.omit(annotation_confidence))), collapse = ";"),
    immune_review_flag = any(as.logical(immune_review_flag), na.rm = TRUE),
    mural_review_flag = any(as.logical(mural_review_flag), na.rm = TRUE)
  ), by = .(gene, cell_type = audited_broad_cell_type)]
}), fill = TRUE)

celltype_summary[, expression_rank_across_celltypes := frank(-avg_expression, ties.method = "min"), by = gene]
celltype_summary[, is_TAL := cell_type == tal_celltype]
setorder(celltype_summary, gene, expression_rank_across_celltypes, cell_type)
fwrite(celltype_summary, "results/tables/p1_tal_gene_celltype_summary.tsv", sep = "\t")

donor_summary <- rbindlist(lapply(present_genes, function(g) {
  v <- p1_expr[g, ]
  tmp <- copy(meta)
  tmp[, gene := g]
  tmp[, expr := as.numeric(v[cell])]
  tmp[, .(
    n_cells = .N,
    avg_expression = mean(expr),
    median_expression = median(expr),
    pct_expressed = mean(expr > 0),
    detected_flag = mean(expr > 0) > 0,
    low_cell_count = .N < 20
  ), by = .(
    gene,
    donor_id,
    disease_status,
    cell_type = audited_broad_cell_type
  )]
}), fill = TRUE)
donor_summary[, is_TAL := cell_type == tal_celltype]
setorder(donor_summary, gene, donor_id, cell_type)
fwrite(donor_summary, "results/tables/p1_tal_gene_by_donor.tsv", sep = "\t")

specificity <- rbindlist(lapply(present_genes, function(g) {
  x <- celltype_summary[gene == g]
  tal <- x[cell_type == tal_celltype]
  non <- x[cell_type != tal_celltype][order(-avg_expression)][1]
  ratio <- ifelse(non$avg_expression == 0, Inf, tal$avg_expression / non$avg_expression)
  cls <- if (tal$expression_rank_across_celltypes == 1 && ratio >= 2) {
    "strong_TAL_preferential"
  } else if (tal$expression_rank_across_celltypes <= 2 && ratio >= 1.2) {
    "moderate_TAL_preferential"
  } else if (tal$expression_rank_across_celltypes <= 3) {
    "broad_expression_with_TAL_component"
  } else {
    "not_TAL_preferential"
  }
  data.table(
    gene = g,
    TAL_avg_expression = tal$avg_expression,
    TAL_pct_expressed = tal$pct_expressed,
    max_nonTAL_celltype = non$cell_type,
    max_nonTAL_avg_expression = non$avg_expression,
    max_nonTAL_pct_expressed = non$pct_expressed,
    specificity_ratio_avg = ratio,
    specificity_delta_avg = tal$avg_expression - non$avg_expression,
    TAL_rank = tal$expression_rank_across_celltypes,
    specificity_class = cls,
    notes = fifelse(is.infinite(ratio), "nonTAL_avg_zero_check_pct_expression", "")
  )
}), fill = TRUE)
fwrite(specificity, "results/tables/p1_tal_gene_specificity.tsv", sep = "\t")

detection_summary <- merge(availability, specificity[, .(gene, TAL_pct_expressed, TAL_avg_expression, TAL_rank, specificity_class)], by = "gene", all.x = TRUE)
fwrite(detection_summary, "results/tables/p1_tal_gene_detection_summary.tsv", sep = "\t")

canonical_tal <- c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2")
canonical_present <- intersect(canonical_tal, rownames(expr))

marker_audit <- fread("results/tables/gse231569_marker_audit.tsv")
audited_markers <- unique(unlist(strsplit(marker_audit[candidate_cell_type == tal_celltype, marker_genes_present], ",")))
audited_markers <- intersect(audited_markers[nzchar(audited_markers)], rownames(expr))

tal_program_reference <- rbindlist(list(
  data.table(program_name = "canonical_TAL_markers", gene = canonical_tal, source = "canonical_marker_set", included = canonical_tal %in% canonical_present, notes = "P1 gene excluded from score when correlating same gene."),
  data.table(program_name = "audited_TAL_cluster_markers", gene = audited_markers, source = "GSE231569_marker_audit", included = TRUE, notes = "Marker set from existing audited marker table.")
), fill = TRUE)
fwrite(tal_program_reference, "results/tables/tal_program_reference.tsv", sep = "\t")

calc_program_score <- function(program_genes, exclude_gene = NULL) {
  genes <- setdiff(intersect(program_genes, rownames(expr)), exclude_gene)
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[genes, , drop = FALSE])
}

cor_rows <- list()
for (g in present_genes) {
  gv <- as.numeric(p1_expr[g, ])
  for (program_name in c("canonical_TAL_markers", "audited_TAL_cluster_markers")) {
    pg <- if (program_name == "canonical_TAL_markers") canonical_present else audited_markers
    score <- calc_program_score(pg, exclude_gene = g)
    for (scope in c("all_cells", "TAL_cells_only")) {
      idx <- if (scope == "all_cells") rep(TRUE, length(gv)) else meta$audited_broad_cell_type == tal_celltype
      ok <- idx & is.finite(gv) & is.finite(score)
      ct <- suppressWarnings(cor.test(gv[ok], score[ok], method = "spearman", exact = FALSE))
      cor_rows[[length(cor_rows) + 1]] <- data.table(
        gene = g,
        program_name,
        analysis_scope = scope,
        n_cells = sum(ok),
        correlation_method = "spearman",
        rho = unname(ct$estimate),
        p_value = ct$p.value
      )
    }
  }
}
cor_dt <- rbindlist(cor_rows)
cor_dt[, fdr := p.adjust(p_value, method = "BH")]
cor_dt[, interpretation := fifelse(rho > 0.30 & fdr < 0.05, "positive_consistency_with_TAL_program",
  fifelse(rho > 0.10, "weak_to_moderate_consistency", "little_TAL_program_consistency")
)]
fwrite(cor_dt, "results/tables/p1_gene_vs_tal_program_correlation.tsv", sep = "\t")

donor_tal <- donor_summary[cell_type == tal_celltype, .(
  n_TAL_donors_available = uniqueN(donor_id),
  n_TAL_donors_detected = sum(detected_flag & !low_cell_count),
  TAL_donor_detection_fraction = sum(detected_flag & !low_cell_count) / uniqueN(donor_id),
  TAL_mean_expression_across_donors = mean(avg_expression),
  TAL_expression_cv_across_donors = ifelse(mean(avg_expression) > 0, sd(avg_expression) / mean(avg_expression), NA_real_)
), by = gene]
donor_tal[, donor_support_class := fifelse(TAL_donor_detection_fraction >= 0.70, "donor_level_stable",
  fifelse(TAL_donor_detection_fraction >= 0.40, "moderate_donor_support", "weak_donor_support")
)]

tiers <- fread("results/tables/candidate_gene_tiers_v0.2.tsv")
p1_meta <- tiers[gene %in% p1_genes, .(
  gene,
  priority_class = current_tier,
  magma_rank,
  magma_p,
  magma_fdr,
  current_tier,
  scrna_top_celltype
)]

best_cor <- cor_dt[program_name == "canonical_TAL_markers" & analysis_scope == "all_cells", .(
  gene,
  TAL_program_rho = rho,
  TAL_program_fdr = fdr
)]

evidence <- Reduce(function(x, y) merge(x, y, by = "gene", all.x = TRUE), list(
  data.table(gene = p1_genes),
  p1_meta,
  specificity[, .(gene, TAL_rank, TAL_avg_expression, TAL_pct_expressed, specificity_ratio_avg, specificity_class)],
  donor_tal,
  best_cor
))

evidence[, overall_evidence_class := fifelse(
  !is.na(TAL_rank) & TAL_rank == 1 & donor_support_class == "donor_level_stable" & specificity_class %in% c("strong_TAL_preferential", "moderate_TAL_preferential"),
  "P1_strong_TAL_context",
  fifelse(
    !is.na(TAL_avg_expression) & TAL_avg_expression > 0 & donor_support_class %in% c("donor_level_stable", "moderate_donor_support") & grepl("TAL|broad", specificity_class),
    "P1_moderate_TAL_context",
    fifelse(!is.na(TAL_avg_expression) & TAL_avg_expression > 0, "P1_context_gene", "P1_uncertain")
  )
)]
evidence[, recommended_interpretation := fifelse(
  gene == "UMOD", "Strong TAL-preferential candidate; representative TAL-associated locus/gene.",
  fifelse(gene == "CLDN14", "TAL/transport-associated candidate; interpret in calcium/ion handling context.",
    fifelse(gene == "CLDN10", "TAL/epithelial transport candidate; useful for TAL program interpretation.",
      fifelse(gene == "CASR", "Calcium-sensing candidate; interpret as biologically important TAL-associated context.",
        fifelse(gene == "HIBADH", "MAGMA-supported TAL-associated candidate; requires TWAS/coloc support.",
          "Renal epithelial/context candidate; interpret cautiously if expression is broad."
        )
      )
    )
  )
)]
evidence[, claim_boundary := "Supports TAL-associated expression context; does not establish causality, TWAS support, or colocalized causal variants."]
fwrite(evidence, "results/tables/p1_tal_gene_evidence.tsv", sep = "\t")

umap <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "cell")
setnames(umap, old = names(umap)[2:3], new = c("UMAP_1", "UMAP_2"))
umap <- merge(umap, meta[, .(cell, audited_broad_cell_type)], by = "cell", all.x = TRUE)

pdf("results/figures/p1_tal_gene_featureplots.pdf", width = 7, height = 5)
for (g in present_genes) {
  d <- copy(umap)
  d[, expr := as.numeric(p1_expr[g, cell])]
  print(
    ggplot(d, aes(UMAP_1, UMAP_2, color = expr)) +
      geom_point(size = 0.12, alpha = 0.65) +
      scale_color_viridis_c(option = "magma") +
      labs(title = paste(g, "normalized expression"), color = "expr") +
      theme_void(base_size = 10)
  )
}
dev.off()

dot_data <- copy(celltype_summary)
dot_data[, cell_type := factor(cell_type, levels = rev(unique(celltype_summary[order(is_TAL, cell_type), cell_type])))]
pdf("results/figures/p1_tal_gene_dotplot.pdf", width = 8, height = 4.8)
print(
  ggplot(dot_data, aes(gene, cell_type)) +
    geom_point(aes(size = pct_expressed, color = avg_expression)) +
    scale_color_viridis_c(option = "plasma") +
    scale_size(range = c(1, 8)) +
    labs(x = NULL, y = NULL, color = "Avg expr", size = "Pct expr") +
    theme_bw(base_size = 10)
)
dev.off()

heat <- dcast(celltype_summary, gene ~ cell_type, value.var = "avg_expression", fill = 0)
heat_long <- melt(heat, id.vars = "gene", variable.name = "cell_type", value.name = "avg_expression")
heat_long[, scaled_avg_expression := as.numeric(scale(avg_expression)), by = gene]
pdf("results/figures/p1_tal_gene_heatmap_by_celltype.pdf", width = 8.5, height = 4.5)
print(
  ggplot(heat_long, aes(cell_type, gene, fill = scaled_avg_expression)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(x = NULL, y = NULL, fill = "Scaled\navg expr") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
)
dev.off()

pdf("results/figures/p1_tal_gene_by_donor_boxplot.pdf", width = 10, height = 5.5)
print(
  ggplot(donor_summary, aes(cell_type, avg_expression, color = cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_point(aes(shape = low_cell_count), position = position_jitter(width = 0.15, height = 0), size = 1.2, alpha = 0.8) +
    facet_wrap(~gene, scales = "free_y", nrow = 2) +
    labs(x = NULL, y = "Donor-level avg expression", shape = "n_cells < 20") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
)
dev.off()

pdf("results/figures/p1_tal_gene_specificity_barplot.pdf", width = 7, height = 4)
print(
  ggplot(specificity, aes(gene, specificity_ratio_avg, fill = specificity_class)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = c(1.2, 2), linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "TAL avg / max non-TAL avg", fill = "Specificity") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
)
dev.off()

pdf("results/figures/p1_gene_vs_tal_program_correlation.pdf", width = 8, height = 4.5)
print(
  ggplot(cor_dt, aes(gene, rho, fill = interpretation)) +
    geom_col(width = 0.7) +
    facet_grid(program_name ~ analysis_scope) +
    geom_hline(yintercept = c(0.1, 0.3), linetype = "dashed", color = "grey45") +
    labs(x = NULL, y = "Spearman rho", fill = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
)
dev.off()

message("wrote P1 TAL evidence pack tables and figures")
