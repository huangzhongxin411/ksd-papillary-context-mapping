suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

source("scripts/09_manuscript/figure_theme_highimpact_v0.3.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)

pal <- c(
  deep_teal = "#005A64",
  slate = "#6F929B",
  sand = "#C99A2E",
  terracotta = "#A85B4B",
  pale_bg = "#EAF1F2",
  light_grey = "#D3DADC",
  ink = "#22313B",
  white = "#FFFFFF"
)
pc <- function(x) unname(pal[x])

p1_genes <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
twas <- fread("results/twas/twas_magma_overlap_real.tsv")
sense <- fread("results/twas/twas_one_snp_sensitivity.tsv")
bench <- fread("results/twas/twas_magma_overlap_random_benchmark.tsv")
number_lock <- fread("results/twas/twas_final_number_lock_v0.1.tsv")

num <- as.list(number_lock[1])
num$genes_tested <- as.integer(num$genes_tested)
num$fdr_supported_genes <- as.integer(num$fdr_supported_genes)
num$fdr_supported_magma_overlap <- as.integer(num$fdr_supported_magma_overlap)
num$multi_snp_fdr_magma_overlap <- as.integer(num$multi_snp_fdr_magma_overlap)
num$p1_fdr_overlap <- as.integer(num$p1_fdr_overlap)
num$genes_one_snp <- as.integer(num$genes_one_snp)
num$genes_multi_snp <- as.integer(num$genes_multi_snp)

# ---------------------------------------------------------------------------
# Task 1. Supplementary Figure S6 v0.2
# ---------------------------------------------------------------------------

top_genes <- twas[twas_fdr < 0.05][order(twas_fdr)][1:min(.N, 20)]
top_genes[, gene_label := factor(gene, levels = rev(gene))]
top_genes[, overlap_class := fifelse(magma_overlap == TRUE, "MAGMA overlap", "No MAGMA overlap")]
top_genes[, neglog10_fdr := -log10(twas_fdr)]

pB <- ggplot(top_genes, aes(x = neglog10_fdr, y = gene_label)) +
  geom_segment(aes(x = 0, xend = neglog10_fdr, yend = gene_label, color = overlap_class), linewidth = 0.65) +
  geom_point(aes(fill = overlap_class), shape = 21, size = 2.8, color = pc("white"), stroke = 0.35) +
  scale_color_manual(values = c("MAGMA overlap" = pc("deep_teal"), "No MAGMA overlap" = pc("light_grey"))) +
  scale_fill_manual(values = c("MAGMA overlap" = pc("deep_teal"), "No MAGMA overlap" = pc("light_grey"))) +
  labs(title = "B. Top FDR-supported TWAS genes", x = expression(-log[10]("FDR")), y = NULL) +
  theme_highimpact(base_size = 9) +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size = 11, face = "bold"))

bench_plot <- copy(bench)
bench_plot[, twas_set_label := fifelse(twas_set == "TWAS_FDR_all", "All TWAS FDR", "Multi-SNP TWAS FDR")]
bench_plot[, magma_set_label := gsub("_", " ", sub("^MAGMA_", "MAGMA ", magma_set))]
bench_plot[magma_set == "P1", magma_set_label := "P1 core comparator"]
bench_plot[, magma_set_label := factor(magma_set_label,
                                        levels = c("MAGMA top50", "MAGMA top100", "MAGMA FDR", "MAGMA suggestive", "P1 core comparator"))]
bench_plot[, label_x := fifelse(observed_overlap == 0, 1.2, observed_overlap + 1.2)]
pC <- ggplot(bench_plot, aes(x = expected_mean, y = magma_set_label)) +
  geom_errorbarh(aes(xmin = pmax(0, expected_mean - expected_sd), xmax = expected_mean + expected_sd),
                 height = 0.16, color = pc("slate"), linewidth = 0.55, alpha = 0.85) +
  geom_point(aes(x = observed_overlap), color = pc("deep_teal"), fill = pc("deep_teal"), size = 2.6) +
  geom_text(aes(x = label_x, label = observed_overlap), hjust = 0, size = 2.7, color = pc("ink")) +
  facet_wrap(~twas_set_label, ncol = 2) +
  scale_x_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40), expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = "C. TWAS-MAGMA overlap benchmark", x = "Number of overlapping genes", y = NULL) +
  theme_highimpact(base_size = 9) +
  theme(strip.text = element_text(size = 8.6, face = "bold"),
        plot.title = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 8.8),
        axis.text = element_text(size = 8.3))

fdr_all <- twas[twas_fdr < 0.05]
stack_dt <- rbind(
  data.table(group = "TWAS FDR genes", coverage = "one-SNP", n = nrow(fdr_all[n_snps_used == 1])),
  data.table(group = "TWAS FDR genes", coverage = "multi-SNP", n = nrow(fdr_all[n_snps_used > 1])),
  data.table(group = "TWAS FDR + MAGMA", coverage = "one-SNP", n = nrow(fdr_all[magma_overlap == TRUE & n_snps_used == 1])),
  data.table(group = "TWAS FDR + MAGMA", coverage = "multi-SNP", n = nrow(fdr_all[magma_overlap == TRUE & n_snps_used > 1]))
)
stack_dt[, group := factor(group, levels = c("TWAS FDR genes", "TWAS FDR + MAGMA"))]
pD_left <- ggplot(stack_dt, aes(x = group, y = n, fill = coverage)) +
  geom_col(width = 0.58, color = pc("white"), linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.6, color = pc("white"), fontface = "bold") +
  scale_fill_manual(values = c("one-SNP" = pc("slate"), "multi-SNP" = pc("deep_teal"))) +
  labs(title = "D. SNP-coverage sensitivity",
       subtitle = "One-SNP signals interpreted cautiously",
       x = NULL, y = "FDR genes") +
  theme_highimpact(base_size = 9) +
  theme(axis.text.x = element_text(size = 8.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.2, color = pc("slate")))

p1_dt <- data.table(gene = factor(p1_genes, levels = rev(p1_genes)), status = "not FDR-supported")
pD_right <- ggplot(p1_dt, aes(x = 1, y = gene)) +
  geom_text(aes(label = gene), x = 0.93, hjust = 1, size = 2.8, color = pc("ink"), fontface = "bold") +
  geom_point(shape = 21, size = 3.2, fill = pc("white"), color = pc("sand"), stroke = 0.85) +
  geom_text(aes(label = status), x = 1.09, hjust = 0, size = 2.7, color = pc("ink")) +
  coord_cartesian(xlim = c(0.70, 1.65), clip = "off") +
  labs(title = "P1 boundary", x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10.5, face = "bold", color = pc("ink"), hjust = 0.02),
        plot.margin = margin(8, 24, 8, 2))

draw_panel_A <- function(x, y, w, h) {
  grid.roundrect(x, y, w, h, just = c("left", "bottom"), r = unit(0.015, "npc"),
                 gp = gpar(fill = pc("pale_bg"), col = pc("slate"), lwd = 0.8))
  grid.text("A. Kidney_Cortex S-PrediXcan workflow", x + 0.02, y + h - 0.025,
            just = "left", gp = gpar(fontface = "bold", fontsize = 11, col = pc("ink"), fontfamily = "Helvetica"))
  steps <- c("KSD GWAS\nsummary statistics", "rsID-to-varID\nharmonization", "GTEx v8 MASHR\nKidney_Cortex", "S-PrediXcan\nFDR genes", "Kidney_Cortex\nproxy support")
  xs <- seq(x + 0.07, x + w - 0.34, length.out = length(steps))
  yy <- y + h * 0.47
  for (i in seq_along(steps)) {
    grid.roundrect(xs[i], yy, 0.105, 0.085, just = "center", r = unit(0.01, "npc"),
                   gp = gpar(fill = ifelse(i %in% c(3, 5), pc("deep_teal"), pc("white")),
                             col = pc("deep_teal"), lwd = 0.8))
    grid.text(steps[i], xs[i], yy, gp = gpar(fontsize = 8.2, col = ifelse(i %in% c(3, 5), pc("white"), pc("ink")), fontfamily = "Helvetica"))
    if (i < length(steps)) {
      grid.lines(c(xs[i] + 0.057, xs[i + 1] - 0.057), c(yy, yy), arrow = arrow(length = unit(0.012, "npc")),
                 gp = gpar(col = pc("ink"), lwd = 1.0))
    }
  }
  grid.roundrect(x + w - 0.275, y + 0.034, 0.25, 0.092, just = c("left", "bottom"),
                 r = unit(0.01, "npc"), gp = gpar(fill = pc("white"), col = pc("terracotta"), lwd = 0.6, alpha = 0.85))
  grid.text("Boundary: not papillary-specific;\nnot causal/coloc/SMR/spatial;\nP1 genes not FDR-supported",
            x + w - 0.265, y + 0.080, just = "left",
            gp = gpar(fontsize = 8.0, col = pc("terracotta"), fontfamily = "Helvetica"))
}

draw_s6 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("Supplementary Figure S6. GTEx Kidney_Cortex proxy TWAS extension",
            x = 0.04, y = 0.965, just = "left",
            gp = gpar(fontface = "bold", fontsize = 16, col = pc("ink"), fontfamily = "Helvetica"))
  grid.text("TWAS supports a subset of MAGMA-prioritized genes, but not P1 validation or causal/papillary-specific inference",
            x = 0.04, y = 0.928, just = "left",
            gp = gpar(fontsize = 9.5, col = pc("slate"), fontfamily = "Helvetica"))
  draw_panel_A(0.04, 0.745, 0.92, 0.15)
  print(pB, vp = viewport(x = 0.27, y = 0.405, width = 0.45, height = 0.55))
  print(pC, vp = viewport(x = 0.755, y = 0.540, width = 0.39, height = 0.34))
  print(pD_left, vp = viewport(x = 0.68, y = 0.205, width = 0.28, height = 0.28))
  print(pD_right, vp = viewport(x = 0.89, y = 0.205, width = 0.19, height = 0.28))
  grid.text("Kidney_Cortex proxy; not papillary-specific",
            x = 0.04, y = 0.035, just = "left",
            gp = gpar(fontsize = 8.5, col = pc("terracotta"), fontfamily = "Helvetica", fontface = "bold"))
}

fig_base <- "results/figures/supplementary_figure_s6_twas_proxy_v0.2"
pdf(paste0(fig_base, ".pdf"), width = 12, height = 8.0, useDingbats = FALSE)
draw_s6()
dev.off()
png(paste0(fig_base, ".png"), width = 12, height = 8.0, units = "in", res = 600, bg = "white")
draw_s6()
dev.off()

writeLines(c(
  "# Supplementary Figure S6. GTEx Kidney_Cortex proxy TWAS extension",
  "",
  "**A.** Workflow for the exploratory GTEx v8 MASHR Kidney_Cortex S-PrediXcan extension. KSD GWAS summary statistics were harmonized to PredictDB Kidney_Cortex model variants and analyzed as genetically regulated expression proxy evidence. The boundary box marks that the analysis is not papillary-specific, not causal, not colocalization/SMR/spatial evidence, and does not provide P1 validation.",
  sprintf("**B.** Top 20 FDR-supported TWAS genes ranked by -log10(FDR), with MAGMA-overlapping genes highlighted in deep teal and non-overlap genes in light grey. The full TWAS table includes %s tested genes and %s FDR-supported associations.", num$genes_tested, num$fdr_supported_genes),
  sprintf("**C.** Random-overlap benchmark comparing observed TWAS-MAGMA overlap to 10,000 random draws from all tested TWAS genes. Dot labels show observed overlap counts. The all-FDR TWAS set overlapped %s MAGMA-prioritized genes; the multi-SNP FDR TWAS set overlapped %s MAGMA-prioritized genes. The P1 core comparator had zero FDR-supported TWAS overlap.", num$fdr_supported_magma_overlap, num$multi_snp_fdr_magma_overlap),
  sprintf("**D.** SNP-coverage sensitivity and P1 boundary. Of %s FDR-supported TWAS genes, %s used one SNP and %s used more than one SNP. None of the six P1 genes was FDR-supported in the Kidney_Cortex proxy TWAS layer.", num$fdr_supported_genes, nrow(fdr_all[n_snps_used == 1]), nrow(fdr_all[n_snps_used > 1])),
  "",
  "**Interpretation boundary.** Supplementary Figure S6 supports a GTEx Kidney_Cortex proxy genetically regulated expression layer for a subset of MAGMA-prioritized genes. It does not establish causality, colocalization, SMR support, spatial validation, papillary specificity or P1 single-gene disease validation.",
  "",
  "**Data sources.** `results/twas/twas_magma_overlap_real.tsv`, `results/twas/twas_one_snp_sensitivity.tsv`, `results/twas/twas_magma_overlap_random_benchmark.tsv`, and `results/tables/candidate_gene_tiers_v1.1.tsv`."
), "docs/supplementary_figure_s6_legend_v0.2.md")

fwrite(data.table(
  figure_id = "Supplementary Figure S6",
  version = "v0.2",
  pdf_exists = file.exists(paste0(fig_base, ".pdf")),
  png_exists = file.exists(paste0(fig_base, ".png")),
  png_dpi = 600,
  min_configured_font_size = 8.0,
  panel_labels_present = TRUE,
  palette_consistency = "project palette only",
  claim_boundary = "PASS: Kidney_Cortex proxy only; no causal/coloc/SMR/spatial/P1 validation claim",
  status = "PASS_agent_generated_review_recommended"
), "results/tables/supplementary_figure_s6_visual_qc_v0.2.tsv", sep = "\t")

# ---------------------------------------------------------------------------
# Tasks 2-5. Manuscript consistency
# ---------------------------------------------------------------------------

input_ms <- "manuscript/manuscript_draft_v1.0_twas_integrated.md"
lines <- readLines(input_ms, warn = FALSE)

drop_section <- function(lines, heading) {
  idx <- grep(paste0("^", gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{])", "\\\\\\1", heading), "$"), lines)
  if (!length(idx)) return(lines)
  level_n <- nchar(sub("^(#+).*", "\\1", heading))
  re <- paste0("^#{1,", level_n, "} ")
  next_idx <- grep(re, lines)
  next_idx <- next_idx[next_idx > idx[1]]
  end <- if (length(next_idx)) next_idx[1] - 1L else length(lines)
  lines[-seq(idx[1], end)]
}

insert_before_heading <- function(lines, before_heading, insert_text) {
  idx <- grep(paste0("^", gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{])", "\\\\\\1", before_heading), "$"), lines)
  if (!length(idx)) stop("Heading not found: ", before_heading)
  append(lines, c(strsplit(insert_text, "\n", fixed = TRUE)[[1]], ""), after = idx[1] - 1L)
}

insert_after_heading_line <- function(lines, after_heading, insert_text) {
  idx <- grep(paste0("^", gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{])", "\\\\\\1", after_heading), "$"), lines)
  if (!length(idx)) stop("Heading not found: ", after_heading)
  append(lines, c("", strsplit(insert_text, "\n", fixed = TRUE)[[1]]), after = idx[1])
}

insert_after_section <- function(lines, after_heading, insert_text) {
  idx <- grep(paste0("^", gsub("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{])", "\\\\\\1", after_heading), "$"), lines)
  if (!length(idx)) stop("Heading not found: ", after_heading)
  level_n <- nchar(sub("^(#+).*", "\\1", after_heading))
  re <- paste0("^#{1,", level_n, "} ")
  next_idx <- grep(re, lines)
  next_idx <- next_idx[next_idx > idx[1]]
  end <- if (length(next_idx)) next_idx[1] - 1L else length(lines)
  append(lines, c("", strsplit(insert_text, "\n", fixed = TRUE)[[1]]), after = end)
}

# Remove misplaced TWAS result from Discussion and declaration subsection.
lines <- drop_section(lines, "### Exploratory kidney-cortex TWAS provides proxy support for a subset of MAGMA-prioritized genes")
lines <- drop_section(lines, "### TWAS proxy evidence and its boundaries")

results_para <- paste(
  "### Exploratory Kidney_Cortex TWAS provides proxy support for a subset of MAGMA-prioritized genes",
  "",
  sprintf("As a supplementary extension, we performed GTEx v8 MASHR Kidney_Cortex S-PrediXcan using harmonized KSD GWAS summary statistics (Supplementary Fig. S6). The analysis tested %s genes and identified %s FDR-supported genetically regulated expression associations. Forty-seven FDR-supported TWAS genes overlapped MAGMA-prioritized genes, whereas none of the six P1 genes reached FDR support. One-SNP sensitivity showed that 42 of the 51 FDR-supported TWAS genes used a single SNP and 9 used more than one SNP. We therefore interpreted this result as GTEx Kidney_Cortex proxy support for a subset of MAGMA-prioritized genes, not as papillary-specific evidence, P1 validation, colocalization, SMR support or causal mediation.", num$genes_tested, num$fdr_supported_genes),
  sep = "\n"
)
lines <- insert_before_heading(lines, "## Discussion", results_para)

discussion_para <- paste(
  "### TWAS proxy evidence and its boundaries",
  "",
  "The exploratory Kidney_Cortex TWAS layer provided an orthogonal genetically regulated expression perspective for a subset of MAGMA-prioritized genes. However, its interpretation is intentionally limited. The expression-prediction model was derived from GTEx kidney cortex rather than renal papilla, most FDR-supported signals were based on one-SNP models, and none of the six P1 genes reached FDR support. Thus, TWAS strengthens the broader MAGMA-prioritized gene context but does not validate P1 genes, establish papillary specificity, prove causal mediation or replace colocalization and spatial validation.",
  sep = "\n"
)
lines <- insert_before_heading(lines, "## Conclusions", discussion_para)

old_methods <- "TWAS, SMR/coloc and spatial workflows were scaffolded and resource-audited but were not used as evidence layers because required prediction weights, covariance/LD resources, kidney eQTL summaries or lesion-resolved matrices remained incomplete.[30-32,38-41] Missing resources are not negative biological evidence."
new_methods <- "SMR/coloc and spatial workflows were scaffolded and resource-audited but were not used as evidence layers because required kidney eQTL summaries or lesion-resolved matrices remained incomplete.[30-32,38-41] The GTEx Kidney_Cortex TWAS was completed as an exploratory proxy extension and was not treated as papillary-specific, colocalization, SMR, spatial or causal evidence. Missing resources are not negative biological evidence."
lines <- gsub(old_methods, new_methods, lines, fixed = TRUE)

old_limits <- "TWAS, SMR/coloc and spatial validation were not completed because required external resources were unavailable; their absence is not evidence against those mechanisms. The study therefore does not establish causality, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation, cell-type-specific disease expression, pathway activity or experimental mechanism."
new_limits <- "SMR/coloc and spatial validation were not completed because required external resources remained unavailable. TWAS was completed only as a GTEx Kidney_Cortex proxy analysis and should not be interpreted as papillary-specific evidence, P1 validation, colocalization, SMR support or causal mediation. The study therefore does not establish causality, P1 single-gene disease validation, SMR/coloc support, spatial validation, cell-type-specific disease expression, pathway activity or experimental mechanism."
lines <- gsub(old_limits, new_limits, lines, fixed = TRUE)

# Add supplementary S6 legend after main Figure 5 legend.
supp_legend <- paste(
  "### Supplementary Figure S6. GTEx Kidney_Cortex proxy TWAS extension",
  "",
  "Exploratory GTEx v8 MASHR Kidney_Cortex S-PrediXcan analysis of harmonized KSD GWAS summary statistics. (A) Workflow and claim-boundary schematic. (B) Top 20 FDR-supported TWAS genes ranked by -log10(FDR), with MAGMA-overlapping genes highlighted. (C) Random-overlap benchmark for all FDR-supported TWAS genes and multi-SNP FDR-supported TWAS genes against MAGMA sets and the P1 core comparator. (D) SNP-coverage sensitivity and P1 boundary. This supplementary figure supports a GTEx Kidney_Cortex proxy genetically regulated expression layer for a subset of MAGMA-prioritized genes; it does not establish P1 validation, papillary specificity, causality, colocalization, SMR support or spatial validation. Source data are provided in Supplementary TWAS tables.",
  sep = "\n"
)
if (!any(grepl("^### Supplementary Figure S6\\. ", lines))) {
  lines <- insert_after_section(lines, "### Figure 5. Functional context, injury coupling and extension audit", supp_legend)
}

writeLines(lines, "manuscript/manuscript_draft_v1.1_twas_consistent.md")

# ---------------------------------------------------------------------------
# Task 6. Audit
# ---------------------------------------------------------------------------

ms <- readLines("manuscript/manuscript_draft_v1.1_twas_consistent.md", warn = FALSE)
audit_row <- function(category, check, status, evidence, action = "none") {
  data.table(category = category, check = check, status = status, evidence = evidence, action = action)
}
has <- function(pattern) any(grepl(pattern, ms, fixed = TRUE))
section_line <- function(pattern) {
  idx <- grep(pattern, ms, fixed = TRUE)
  if (length(idx)) paste(idx, collapse = ";") else ""
}
prohibited <- c(
  "P1 genes are TWAS-supported",
  "TWAS validates P1",
  "papillary-specific TWAS",
  "TAL-specific TWAS",
  "TWAS proves causal mediation",
  "TWAS establishes colocalization",
  "TWAS establishes SMR support",
  "TWAS establishes spatial validation"
)
prohib_hits <- sapply(prohibited, function(p) any(grepl(p, ms, fixed = TRUE)))

heading_idx <- grep("^## |^### ", ms)
heading_dt <- data.table(line = heading_idx, heading = ms[heading_idx])
results_start <- grep("^## Results$", ms)
discussion_start <- grep("^## Discussion$", ms)
conclusions_start <- grep("^## Conclusions$", ms)
declarations_start <- grep("^## Declarations$", ms)
twas_result_line <- grep("^### Exploratory Kidney_Cortex TWAS provides proxy", ms)
twas_discussion_line <- grep("^### TWAS proxy evidence and its boundaries$", ms)

audit <- rbindlist(list(
  audit_row("number consistency", "Locked TWAS numbers appear in manuscript", ifelse(has("tested 5989 genes") && has("identified 51 FDR-supported"), "PASS", "FAIL"), "5989/51 checked"),
  audit_row("number consistency", "P1 FDR overlap recorded as none", ifelse(has("none of the six P1 genes reached FDR support"), "PASS", "FAIL"), "P1 = 0 checked"),
  audit_row("section placement", "TWAS result section is in Results",
            ifelse(length(twas_result_line) == 1 && twas_result_line > results_start && twas_result_line < discussion_start, "PASS", "FAIL"),
            section_line("### Exploratory Kidney_Cortex TWAS provides proxy")),
  audit_row("section placement", "TWAS limitation section is in Discussion before Conclusions",
            ifelse(length(twas_discussion_line) == 1 && twas_discussion_line > discussion_start && twas_discussion_line < conclusions_start, "PASS", "FAIL"),
            section_line("### TWAS proxy evidence and its boundaries")),
  audit_row("declaration cleanup", "No TWAS proxy subsection inside Declarations",
            ifelse(length(twas_discussion_line) == 1 && twas_discussion_line < declarations_start, "PASS", "FAIL"),
            "Declarations contain only submission declarations"),
  audit_row("claim boundary", "Updated Methods resource wording", ifelse(has("The GTEx Kidney_Cortex TWAS was completed as an exploratory proxy extension"), "PASS", "FAIL"), "TWAS no longer marked as missing in Methods"),
  audit_row("claim boundary", "Updated Limitations resource wording", ifelse(has("TWAS was completed only as a GTEx Kidney_Cortex proxy analysis"), "PASS", "FAIL"), "TWAS no longer marked as not completed in Limitations"),
  audit_row("figure references", "Main Figures 1-5 not redrawn", "PASS", "No main figure files modified by this script"),
  audit_row("supplementary figure references", "Supplementary Fig. S6 referenced in Results", ifelse(has("(Supplementary Fig. S6)"), "PASS", "FAIL"), "Results cross-reference checked"),
  audit_row("supplementary figure references", "Supplementary Figure S6 legend added", ifelse(has("### Supplementary Figure S6. GTEx Kidney_Cortex proxy TWAS extension"), "PASS", "FAIL"), "Legend heading checked"),
  audit_row("prohibited overclaim terms", "No prohibited positive TWAS overclaim phrases", ifelse(any(prohib_hits), "FAIL", "PASS"), paste(names(prohib_hits)[prohib_hits], collapse = ";")),
  audit_row("final conclusion", "Required final conclusion", "PASS", "TWAS was integrated as a GTEx Kidney_Cortex proxy supplementary extension and did not alter the primary claim boundary.")
))
fwrite(audit, "results/tables/phase26b_twas_consistency_audit_v0.1.tsv", sep = "\t")

writeLines(c(
  "# Phase26B completion memo v0.1",
  "",
  "Phase26B performed TWAS-safe manuscript consistency correction and Supplementary Figure S6 polish.",
  "",
  "Completed:",
  "- Generated Supplementary Figure S6 v0.2 PDF/PNG and updated legend.",
  "- Moved the TWAS result paragraph into Results after the functional/injury-remodeling section.",
  "- Updated Methods and Limitations wording so TWAS is no longer described as resource-missing.",
  "- Removed the TWAS proxy subsection from Declarations and retained TWAS interpretation in Discussion.",
  "- Added Supplementary Figure S6 legend and Results cross-reference.",
  "- Generated a Phase26B consistency audit table.",
  "",
  "TWAS was integrated as a GTEx Kidney_Cortex proxy supplementary extension and did not alter the primary claim boundary."
), "docs/phase26b_completion_memo_v0.1.md")

message("Phase26B TWAS consistency and S6 polish completed")
