suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_set <- function(path, fallback = character()) if (file.exists(path)) unique(fread(path, header = FALSE)[[1]]) else fallback
sets <- list(
  P1_core_TAL_candidates = c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2"),
  MAGMA_top50 = read_set("results/gene_sets/magma_top50.txt"),
  MAGMA_top100 = read_set("results/gene_sets/magma_top100.txt"),
  MAGMA_FDR = read_set("results/gene_sets/magma_fdr05.txt"),
  MAGMA_suggestive = read_set("results/gene_sets/magma_suggestive_p1e4.txt"),
  TAL_marker_set = c("UMOD", "SLC12A1", "KCNJ1", "CLDN10", "CLDN16", "CLDN19", "FXYD2"),
  injury_remodeling_marker_set = c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "HAVCR1", "CCL2", "CCL7", "VCAM1", "KRT8", "KRT18")
)

meta <- fread("config/gse73680_sample_metadata_curated.tsv")
meta <- meta[include_in_analysis == TRUE & group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla")]
meta[, group_curated := factor(group_curated, levels = c("control_or_adjacent", "plaque_or_stone_papilla"))]
mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
sample_ids <- intersect(meta$sample_id, names(mat))
meta <- meta[match(sample_ids, sample_id)]
expr <- as.matrix(mat[, ..sample_ids])
rownames(expr) <- mat$gene
mode(expr) <- "numeric"
expr <- log2(expr + 1)
zexpr <- t(scale(t(expr)))
zexpr[!is.finite(zexpr)] <- NA_real_

score_rows <- rbindlist(lapply(names(sets), function(nm) {
  genes <- unique(sets[[nm]])
  detected <- intersect(genes, rownames(zexpr))
  frac <- length(detected) / max(1, length(genes))
  status <- fcase(frac >= 0.70, "full_module", frac >= 0.40, "exploratory_module", default = "do_not_test")
  scores <- if (length(detected) > 0 && status != "do_not_test") colMeans(zexpr[detected, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(sample_ids))
  data.table(module_name = nm, sample_id = sample_ids, module_score = as.numeric(scores),
             n_genes_input = length(genes), n_genes_detected = length(detected),
             detected_fraction = frac, module_status = status)
}), fill = TRUE)
fwrite(score_rows, file.path(table_dir, "gse73680_module_score_matrix.tsv"), sep = "\t")

score_mat <- dcast(score_rows[module_status != "do_not_test"], module_name ~ sample_id, value.var = "module_score")
module_meta <- unique(score_rows[, .(module_name, n_genes_input, n_genes_detected, detected_fraction, module_status)])
score_expr <- as.matrix(score_mat[, ..sample_ids])
rownames(score_expr) <- score_mat$module_name
mode(score_expr) <- "numeric"
design <- model.matrix(~ group_curated, data = meta)
coef_name <- "group_curatedplaque_or_stone_papilla"
method <- "limma"
fit <- tryCatch({
  if (uniqueN(meta$patient_id) < nrow(meta)) {
    corfit <- duplicateCorrelation(score_expr, design, block = meta$patient_id)
    method <<- "limma_duplicateCorrelation_module_score"
    lmFit(score_expr, design, block = meta$patient_id, correlation = corfit$consensus)
  } else {
    method <<- "limma_sample_level_module_score"
    lmFit(score_expr, design)
  }
}, error = function(e) {
  method <<- paste0("limma_sample_level_module_score_fallback_after_duplicateCorrelation_error: ", conditionMessage(e))
  lmFit(score_expr, design)
})
fit <- eBayes(fit)
tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), keep.rownames = "module_name")
resp <- merge(module_meta, tt[, .(module_name, logFC, P.Value, adj.P.Val)], by = "module_name", all.x = TRUE)
means <- merge(score_rows, meta[, .(sample_id, group_curated)], by = "sample_id")[, .(mean = mean(module_score, na.rm = TRUE), N = .N), by = .(module_name, group_curated)]
resp[, `:=`(
  group_1 = "control_or_adjacent",
  group_2 = "plaque_or_stone_papilla",
  mean_group_1 = means[group_curated == "control_or_adjacent"][match(module_name, means[group_curated == "control_or_adjacent"]$module_name), mean],
  mean_group_2 = means[group_curated == "plaque_or_stone_papilla"][match(module_name, means[group_curated == "plaque_or_stone_papilla"]$module_name), mean],
  delta = logFC,
  p_value = P.Value,
  fdr = adj.P.Val,
  method = method
)]
resp[, interpretation := fcase(
  module_status == "do_not_test", "not tested due to low detected fraction",
  fdr < 0.05, "disease-context module association",
  p_value < 0.05, "nominal exploratory module signal",
  default = "no detectable module score difference"
)]
resp <- resp[, .(module_name, n_genes_input, n_genes_detected, detected_fraction, group_1, group_2,
                 mean_group_1, mean_group_2, delta, p_value, fdr, method, module_status, interpretation)]
fwrite(resp, file.path(table_dir, "gse73680_module_score_response.tsv"), sep = "\t")

plot_dt <- merge(score_rows[module_status != "do_not_test"], meta[, .(sample_id, patient_id, group_curated)], by = "sample_id")
patient_score <- plot_dt[, .(patient_level_module_score = mean(module_score, na.rm = TRUE),
                             n_samples = .N,
                             n_genes_input = n_genes_input[1],
                             n_genes_detected = n_genes_detected[1],
                             detected_fraction = detected_fraction[1],
                             module_status = module_status[1]),
                          by = .(module_name, patient_id, group_curated)]
fwrite(patient_score, file.path(table_dir, "gse73680_patient_level_module_score_matrix.tsv"), sep = "\t")
patient_wide <- dcast(patient_score, module_name + patient_id ~ group_curated, value.var = "patient_level_module_score")
patient_resp <- rbindlist(lapply(unique(patient_score$module_name), function(mn) {
  x <- patient_wide[module_name == mn]
  paired <- x[is.finite(control_or_adjacent) & is.finite(plaque_or_stone_papilla)]
  if (nrow(paired) >= 3) {
    tt_p <- t.test(paired$plaque_or_stone_papilla, paired$control_or_adjacent, paired = TRUE)
    pval <- tt_p$p.value
    method_s <- "patient_level_paired_t_test_module_score"
  } else {
    pval <- NA_real_
    method_s <- "patient_level_paired_t_test_not_run_n_lt_3"
  }
  mm <- module_meta[module_name == mn]
  data.table(module_name = mn,
             n_genes_input = mm$n_genes_input,
             n_genes_detected = mm$n_genes_detected,
             detected_fraction = mm$detected_fraction,
             module_status = mm$module_status,
             n_paired_patients = nrow(paired),
             mean_control_or_adjacent = mean(paired$control_or_adjacent, na.rm = TRUE),
             mean_plaque_or_stone_papilla = mean(paired$plaque_or_stone_papilla, na.rm = TRUE),
             paired_delta = mean(paired$plaque_or_stone_papilla - paired$control_or_adjacent, na.rm = TRUE),
             p_value = pval,
             method = method_s)
}), fill = TRUE)
patient_resp[, fdr := p.adjust(p_value, method = "BH")]
patient_resp[, interpretation := fcase(
  fdr < 0.05, "patient-level paired module support",
  p_value < 0.05, "patient-level nominal exploratory module signal",
  default = "no detectable patient-level paired module difference"
)]
fwrite(patient_resp, file.path(table_dir, "gse73680_patient_level_module_response.tsv"), sep = "\t")
p <- ggplot(plot_dt, aes(group_curated, module_score, fill = group_curated)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35, color = "#777777") +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.4, alpha = 0.75, color = "#3D4B53") +
  facet_wrap(~ module_name, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(control_or_adjacent = "#8AA0A8", plaque_or_stone_papilla = "#B08A45")) +
  labs(x = NULL, y = "Mean z-score module score", title = "GSE73680 disease-context module score response") +
  theme_bw(base_size = 8.5) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "gse73680_module_score_boxplot.pdf"), p, width = 8.5, height = 8.2)
message("wrote GSE73680 module score response")
