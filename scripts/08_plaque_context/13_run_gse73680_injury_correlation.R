suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- "results/gse73680/tables"
fig_dir <- "results/gse73680/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(file.path(table_dir, "gse73680_module_score_matrix.tsv"))
scores_w <- dcast(scores[module_status != "do_not_test"], sample_id ~ module_name, value.var = "module_score")
mat <- fread("data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz")
mat <- mat[!is.na(gene) & gene != ""]
injury <- c("SPP1", "MMP7", "MMP9", "GPNMB", "COL1A1", "COL1A2", "HAVCR1", "CCL2", "CCL7", "VCAM1", "KRT8", "KRT18")
sample_ids <- intersect(scores_w$sample_id, names(mat))
expr <- as.matrix(mat[gene %in% injury, ..sample_ids])
rownames(expr) <- mat[gene %in% injury]$gene
mode(expr) <- "numeric"
expr <- log2(expr + 1)
zexpr <- t(scale(t(expr)))
zexpr[!is.finite(zexpr)] <- NA_real_
inj_dt <- melt(as.data.table(zexpr, keep.rownames = "injury_marker_or_module"), id.vars = "injury_marker_or_module",
               variable.name = "sample_id", value.name = "injury_score")
inj_module <- scores[module_name == "injury_remodeling_marker_set", .(sample_id, injury_score = module_score)]
inj_module[, injury_marker_or_module := "injury_remodeling_module"]
inj_dt <- rbind(inj_dt, inj_module, fill = TRUE)

module_names <- setdiff(names(scores_w), "sample_id")
module_names <- setdiff(module_names, "injury_remodeling_marker_set")
res <- rbindlist(lapply(module_names, function(mn) {
  x <- scores_w[[mn]]
  rbindlist(lapply(unique(inj_dt$injury_marker_or_module), function(im) {
    y_dt <- inj_dt[injury_marker_or_module == im][match(scores_w$sample_id, sample_id)]
    ok <- is.finite(x) & is.finite(y_dt$injury_score)
    ct <- suppressWarnings(cor.test(x[ok], y_dt$injury_score[ok], method = "spearman", exact = FALSE))
    data.table(module_name = mn, injury_marker_or_module = im, n_samples = sum(ok), cor_method = "Spearman",
               rho = unname(ct$estimate), p_value = ct$p.value)
  }))
}), fill = TRUE)
res[, fdr := p.adjust(p_value, method = "BH")]
res[, interpretation := fcase(
  rho >= 0.50 & fdr < 0.05, "strong disease-context association",
  rho >= 0.30 & rho < 0.50 & p_value < 0.05, "moderate exploratory association",
  default = "weak or no association"
)]
fwrite(res, file.path(table_dir, "gse73680_module_injury_correlation.tsv"), sep = "\t")

heat <- res[injury_marker_or_module == "injury_remodeling_module"]
heat[, module_name := factor(module_name, levels = module_name[order(rho)])]
p <- ggplot(heat, aes("Injury remodeling module", module_name, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#4C72B0", mid = "white", high = "#B35C44", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Module association with injury-remodeling program") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "gse73680_module_injury_correlation_heatmap.pdf"), p, width = 5.4, height = 4.6)
message("wrote GSE73680 module-injury correlation")
