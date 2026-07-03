suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

p1 <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
meta[, group_curated := factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))]
mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
mat <- mat[gene %in% p1]
sample_ids <- intersect(meta$sample_id, names(mat))
meta <- meta[match(sample_ids, sample_id)]
expr <- as.matrix(mat[, ..sample_ids])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
expr <- log2(expr + 1)

design <- model.matrix(~ group_curated, data = meta)
coef_name <- "group_curatedplaque_or_stone_papilla"
method <- "limma"
fit <- tryCatch({
  if (uniqueN(meta$patient_id) < nrow(meta)) {
    corfit <- duplicateCorrelation(expr, design, block = meta$patient_id)
    f <- lmFit(expr, design, block = meta$patient_id, correlation = corfit$consensus)
    attr(f, "consensus_correlation") <- corfit$consensus
    method <<- "limma_duplicateCorrelation_log2"
    f
  } else {
    method <<- "limma_sample_level_log2"
    lmFit(expr, design)
  }
}, error = function(e) {
  method <<- paste0("limma_sample_level_log2_fallback_after_duplicateCorrelation_error: ", conditionMessage(e))
  lmFit(expr, design)
})
fit <- eBayes(fit)
tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), keep.rownames = "gene")

summary_dt <- rbindlist(lapply(p1, function(g) {
  vals <- data.table(sample_id = sample_ids, value = as.numeric(expr[g, sample_ids]), group = meta$group_curated)
  stats <- vals[, .(N = .N, mean = mean(value, na.rm = TRUE)), by = group]
  row <- tt[gene == g]
  data.table(
    gene = g,
    n_detected = sum(!is.na(expr[g, sample_ids]) & expr[g, sample_ids] > 0),
    group_1 = "control_or_adjacent",
    group_2 = "plaque_or_stone_papilla",
    n_group_1 = stats[group == "control_or_adjacent", N],
    n_group_2 = stats[group == "plaque_or_stone_papilla", N],
    mean_group_1 = stats[group == "control_or_adjacent", mean],
    mean_group_2 = stats[group == "plaque_or_stone_papilla", mean],
    logFC = row$logFC,
    p_value = row$P.Value,
    fdr = row$adj.P.Val,
    direction = fifelse(row$logFC > 0, "higher_in_plaque_or_stone_papilla", "lower_in_plaque_or_stone_papilla"),
    method = method,
    interpretation = fcase(
      row$adj.P.Val < 0.05, "disease-context differential expression support",
      row$P.Value < 0.05, "nominal exploratory signal",
      default = "no detectable disease-context difference"
    )
  )
}), fill = TRUE)
fwrite(summary_dt, file.path(table_dir, "gse73680_p1_gene_response.tsv"), sep = "\t")

plot_dt <- melt(as.data.table(expr, keep.rownames = "gene"), id.vars = "gene", variable.name = "sample_id", value.name = "log2_expression")
plot_dt <- merge(plot_dt, meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
patient_dt <- plot_dt[, .(patient_level_log2_expression = mean(log2_expression, na.rm = TRUE)),
                      by = .(gene, patient_id, group_curated)]
patient_wide <- dcast(patient_dt, gene + patient_id ~ group_curated, value.var = "patient_level_log2_expression")
patient_resp <- rbindlist(lapply(p1, function(g) {
  x <- patient_wide[gene == g]
  paired <- x[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
  if (nrow(paired) >= 3) {
    tt_p <- t.test(paired$plaque_or_stone_papilla, paired$control_or_adjacent, paired = TRUE)
    pval <- tt_p$p.value
    method_s <- "patient_level_paired_t_test_log2"
  } else {
    pval <- NA_real_
    method_s <- "patient_level_paired_t_test_not_run_n_lt_3"
  }
  data.table(gene = g,
             n_paired_patients = nrow(paired),
             mean_control_or_adjacent = mean(paired$control_or_adjacent, na.rm = TRUE),
             mean_plaque_or_stone_papilla = mean(paired$plaque_or_stone_papilla, na.rm = TRUE),
             paired_delta = mean(paired$plaque_or_stone_papilla - paired$control_or_adjacent, na.rm = TRUE),
             p_value = pval,
             method = method_s)
}), fill = TRUE)
patient_resp[, fdr := p.adjust(p_value, method = "BH")]
patient_resp[, interpretation := fcase(
  fdr < 0.05, "patient-level paired disease-context expression support",
  p_value < 0.05, "patient-level nominal exploratory signal",
  default = "no detectable patient-level paired difference"
)]
fwrite(patient_resp, file.path(table_dir, "gse73680_patient_level_p1_gene_response.tsv"), sep = "\t")
plot_dt[, gene := factor(gene, levels = p1)]
p_box <- ggplot(plot_dt, aes(group_curated, log2_expression, fill = group_curated)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35, color = "#777777") +
  geom_point(aes(group = patient_id), position = position_jitter(width = 0.08, height = 0), size = 1.6, alpha = 0.75, color = "#3D4B53") +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(control_or_adjacent = "#8AA0A8", plaque_or_stone_papilla = "#B08A45")) +
  labs(x = NULL, y = "log2(expression + 1)", title = "GSE73680 P1 gene disease-context response") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic"), axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "gse73680_p1_gene_boxplot.pdf"), p_box, width = 8.5, height = 5.6)

z <- t(scale(t(expr)))
z[!is.finite(z)] <- 0
heat_dt <- melt(as.data.table(z, keep.rownames = "gene"), id.vars = "gene", variable.name = "sample_id", value.name = "z")
ord <- meta[order(group_curated, patient_id), sample_id]
heat_dt[, sample_id := factor(sample_id, levels = ord)]
heat_dt[, gene := factor(gene, levels = rev(p1))]
p_heat <- ggplot(heat_dt, aes(sample_id, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C72B0", mid = "white", high = "#B35C44", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "z", title = "GSE73680 P1 gene expression heatmap") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(face = "italic"), panel.grid = element_blank())
ggsave(file.path(fig_dir, "gse73680_p1_gene_heatmap.pdf"), p_heat, width = 8.5, height = 3.6)
message("wrote GSE73680 P1 gene response")
