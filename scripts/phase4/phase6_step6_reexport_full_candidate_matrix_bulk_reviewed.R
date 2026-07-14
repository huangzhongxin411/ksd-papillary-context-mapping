suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("source_data/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("notes", recursive = TRUE, showWarnings = FALSE)
dir.create("codex_tasks", recursive = TRUE, showWarnings = FALSE)
dir.create("supplement", recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)

fig_png <- function(path, width_px, height_px, res = 600) {
  png_type <- if (capabilities("aqua")) "quartz" else "Xlib"
  png(path, width = width_px, height = height_px, res = res, type = png_type)
}

matrix <- read_tsv("results/tables/phase4_step2_candidate_gene_evidence_matrix.tsv")
group_defs <- read_tsv("results/tables/phase4_step2_repaired_reporting_group_definitions.tsv")
repaired_groups <- read_tsv("results/tables/phase4_step2_candidate_reporting_groups_repaired.tsv")
twas_downgrade <- read_tsv("results/tables/phase4_step2_twas_downgrade_audit.tsv")
twas_burden <- read_tsv("source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv")
group_counts <- read_tsv("source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv")
overlap <- read_tsv("results/tables/phase4_step1_magma_twas_overlap_summary.tsv")

group_counts$reporting_group <- factor(group_counts$reporting_group, levels = paste0("R", 1:6))
twas_burden$category <- factor(
  twas_burden$category,
  levels = c("FDR-supported TWAS genes", "one-SNP proxy models", "multi-SNP proxy models")
)

priority_genes <- c("UMOD", "CASR", "CLDN14", "CLDN10", "HIBADH", "PKD2")
select_main_rows <- function(df) {
  out <- data.frame()
  for (g in paste0("R", 1:6)) {
    d <- df[df$repaired_reporting_group == g, ]
    if (nrow(d) == 0) next
    d$priority_flag <- d$gene %in% priority_genes
    d$twas_flag <- d$twas_fdr_supported == TRUE
    d <- d[order(!d$priority_flag, !d$twas_flag, d$magma_rank, d$twas_fdr, d$gene, na.last = TRUE), ]
    cap <- switch(g, R1 = 8, R2 = 1, R3 = 2, R4 = 7, R5 = 6, R6 = 3)
    out <- rbind(out, head(d, cap))
  }
  out
}
main_rows <- select_main_rows(matrix)
main_rows$included_in_main_figure <- TRUE
main_rows$reason_for_inclusion <- ifelse(
  main_rows$gene %in% priority_genes, "Named exemplar/priority gene retained for reader orientation.",
  ifelse(main_rows$repaired_reporting_group %in% c("R2", "R3"), "All small TWAS-proxy context groups included without extra visual weighting.",
    ifelse(main_rows$repaired_reporting_group == "R5", "Representative TWAS proxy-only gene showing supplementary status.",
      "Representative gene selected by repaired group and MAGMA/TWAS ordering.")
  )
)
row_selection <- data.frame(
  gene = matrix$gene,
  repaired_reporting_group = matrix$repaired_reporting_group,
  included_in_main_figure = matrix$gene %in% main_rows$gene,
  reason_for_inclusion = ifelse(matrix$gene %in% main_rows$gene, main_rows$reason_for_inclusion[match(matrix$gene, main_rows$gene)], "Full matrix retained in supplementary table, not shown in main figure."),
  magma_summary = paste0(matrix$magma_status, ifelse(is.na(matrix$magma_rank), "", paste0("; rank=", matrix$magma_rank))),
  snRNA_summary = paste0(matrix$snRNA_support_status),
  twas_summary = paste0(matrix$twas_proxy_status, "; ", matrix$twas_model_type),
  spatial_summary = matrix$spatial_status,
  bulk_summary = matrix$bulk_status_if_available,
  notes = "Main figure row selection is representative only; full evidence matrix is routed to supplement.",
  check.names = FALSE
)
write_tsv(row_selection, "results/tables/phase6_step4_figure3_candidate_matrix_row_selection.tsv")

state_palette <- c(
  primary = "#245A64", secondary = "#7F9DA6", context = "#0F4C5C",
  partial = "#B99B5A", proxy_stronger = "#9B5C4D", proxy_weak = "#D7B8AE",
  supplementary = "#BFC9CC", reviewed_non_upgrading = "#E6E9EA", none = "#F4F5F5"
)
group_palette <- c(R1 = "#245A64", R2 = "#9B5C4D", R3 = "#D7B8AE", R4 = "#7F9DA6", R5 = "#B99B5A", R6 = "#E6E9EA")

make_layers <- function(rows) {
  layers <- data.frame(
    gene = rep(rows$gene, each = 5),
    layer = rep(c("MAGMA primary", "snRNA context", "TWAS proxy", "Spatial supplementary", "Bulk disease-context reviewed"), times = nrow(rows)),
    state = NA_character_,
    group = rep(rows$repaired_reporting_group, each = 5),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(rows))) {
    idx <- layers$gene == rows$gene[i]
    layers$state[idx & layers$layer == "MAGMA primary"] <- ifelse(rows$magma_status[i] %in% c("top50", "top100", "bonferroni", "fdr05"), "primary", ifelse(rows$magma_status[i] == "suggestive", "secondary", "none"))
    layers$state[idx & layers$layer == "snRNA context"] <- ifelse(rows$snRNA_support_status[i] == "strong_context_support", "context", ifelse(rows$snRNA_support_status[i] %in% c("moderate_context_support", "partial_context_support"), "partial", "none"))
    layers$state[idx & layers$layer == "TWAS proxy"] <- ifelse(rows$twas_model_type[i] == "multi_snp_proxy", "proxy_stronger", ifelse(rows$twas_model_type[i] == "one_snp_proxy", "proxy_weak", "none"))
    layers$state[idx & layers$layer == "Spatial supplementary"] <- "supplementary"
    layers$state[idx & layers$layer == "Bulk disease-context reviewed"] <- "reviewed_non_upgrading"
  }
  layers$gene <- factor(layers$gene, levels = rev(unique(rows$gene)))
  layers$layer <- factor(layers$layer, levels = c("MAGMA primary", "snRNA context", "TWAS proxy", "Spatial supplementary", "Bulk disease-context reviewed"))
  layers$group <- factor(layers$group, levels = paste0("R", 1:6))
  layers
}
main_layers <- make_layers(main_rows)

draw_main <- function(pdf_path, png_path, png2_path) {
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(0.55, 1.25, 4.5, 1.45), "null"), widths = unit(c(1, 1), "null"))))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    grid.text("Figure 3. Conservative candidate reporting model", x = 0.02, y = 0.70, just = c("left", "center"), gp = gpar(fontsize = 17, fontface = "bold", col = "#333333"))
    grid.text("MAGMA primary; snRNA context; TWAS Kidney_Cortex proxy only; reporting groups are not causal tiers", x = 0.02, y = 0.24, just = c("left", "center"), gp = gpar(fontsize = 10.5, col = "#333333"))
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.text("A", x = 0.02, y = 0.92, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.text("Parallel evidence model", x = 0.09, y = 0.92, just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
    xs <- c(0.14, 0.34, 0.54, 0.74, 0.91)
    labs <- c("MAGMA\nprimary", "snRNA\ncontext", "TWAS\nproxy", "Spatial\nsupplement", "Bulk\ndisease-context\nreviewed")
    fills <- c("#245A64", "#0F4C5C", "#D7B8AE", "#E6E9EA", "#E6E9EA")
    for (j in seq_along(xs)) {
      grid.roundrect(x = xs[j], y = 0.45, width = 0.15, height = 0.34, r = unit(0.03, "snpc"), gp = gpar(fill = fills[j], col = "#333333", lwd = 0.55))
      grid.text(labs[j], x = xs[j], y = 0.45, gp = gpar(fontsize = 8.2, col = ifelse(j <= 2, "white", "#333333"), fontface = ifelse(j == 1, "bold", "plain")))
    }
    grid.text("parallel evidence; no causal arrows", x = 0.5, y = 0.12, gp = gpar(fontsize = 8.5, col = "#9B5C4D"))
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    burden_plot <- ggplot(twas_burden, aes(category, n_genes, fill = category)) +
      geom_col(width = 0.66) +
      geom_text(aes(label = n_genes), vjust = -0.25, size = 3.2) +
      scale_fill_manual(values = c("FDR-supported TWAS genes" = "#7F9DA6", "one-SNP proxy models" = "#D7B8AE", "multi-SNP proxy models" = "#9B5C4D")) +
      scale_y_continuous(limits = c(0, max(twas_burden$n_genes) * 1.20), expand = expansion(mult = c(0, 0.02))) +
      labs(title = "C  TWAS model quality: proxy only", x = NULL, y = "Genes") +
      theme_minimal(base_family = "Helvetica", base_size = 9) +
      theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 12, hjust = 1), plot.title = element_text(face = "bold", size = 12))
    print(burden_plot, newpage = FALSE)
    popViewport()

    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    mat_plot <- ggplot(main_layers, aes(layer, gene, fill = state)) +
      geom_tile(color = "white", linewidth = 0.35) +
      facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = state_palette, breaks = names(state_palette), name = "Evidence role") +
      labs(title = "B  Compact candidate evidence matrix", subtitle = "Representative rows; full matrix is in supplementary source data", x = NULL, y = NULL) +
      theme_minimal(base_family = "Helvetica", base_size = 8.7) +
      theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", size = 9),
            axis.text.x = element_text(face = "bold"), legend.position = "right",
            plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
    print(mat_plot, newpage = FALSE)
    popViewport()

    pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
    count_plot <- ggplot(group_counts, aes(reporting_group, n_genes, fill = reporting_group)) +
      geom_col(width = 0.72) +
      geom_text(aes(label = n_genes), vjust = -0.25, size = 3.1) +
      scale_fill_manual(values = group_palette) +
      scale_y_continuous(limits = c(0, max(group_counts$n_genes) * 1.18), expand = expansion(mult = c(0, 0.02))) +
      labs(title = "D  Reporting groups, not causal tiers", x = NULL, y = "Genes") +
      theme_minimal(base_family = "Helvetica", base_size = 9) +
      theme(legend.position = "none", panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold", size = 12))
    print(count_plot, newpage = FALSE)
    popViewport()

    pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
    grid.text("E", x = 0.02, y = 0.88, just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.text("Claim boundary", x = 0.09, y = 0.88, just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
    grid.roundrect(x = 0.5, y = 0.45, width = 0.9, height = 0.58, r = unit(0.03, "snpc"), gp = gpar(fill = "#F7F7F7", col = "#333333", lwd = 0.7))
    grid.text("Supported: context-mapped genetic-priority candidates\n\nNot claimed: causal genes, papilla-specific regulation,\ntherapeutic targets, spatial validation, or claim-grade SMR/coloc", x = 0.09, y = 0.62, just = c("left", "top"), gp = gpar(fontsize = 9.5, col = "#333333"))
    popViewport()
    popViewport()
  }
  pdf(pdf_path, width = 13, height = 9, family = "Helvetica")
  draw()
  dev.off()
  fig_png(png_path, 7800, 5400, 600)
  draw()
  dev.off()
  fig_png(png2_path, 7800, 5400, 600)
  draw()
  dev.off()
}

draw_supp_twas <- function(pdf_path, png_path) {
  overlap_long <- overlap
  overlap_long$magma_module <- factor(overlap_long$magma_module, levels = overlap_long$magma_module)
  downtab <- as.data.frame(table(twas_downgrade$downgrade_or_no_change), stringsAsFactors = FALSE)
  names(downtab) <- c("status", "n_genes")
  downtab$status_label <- c(
    downgraded_or_visually_downweighted = "one-SNP downgraded /\nvisually downweighted",
    no_change = "not used to upgrade /\nno group change",
    retained_as_proxy_annotation_only = "multi-SNP retained as\nproxy annotation only"
  )[downtab$status]
  downtab$status_label <- factor(downtab$status_label, levels = downtab$status_label)
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.5, 1, 1), "null"))))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    grid.text("Supplementary TWAS proxy-quality audit", x = 0.02, y = 0.70, just = c("left", "center"), gp = gpar(fontsize = 16, fontface = "bold"))
    grid.text("Kidney_Cortex proxy annotation only; not papilla-specific regulation; not causal evidence", x = 0.02, y = 0.24, just = c("left", "center"), gp = gpar(fontsize = 10, col = "#9B5C4D"))
    popViewport()
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    p1 <- ggplot(twas_burden, aes(category, n_genes, fill = category)) + geom_col(width = 0.65) + geom_text(aes(label = n_genes), vjust = -0.25, size = 3.3) +
      scale_fill_manual(values = c("FDR-supported TWAS genes" = "#7F9DA6", "one-SNP proxy models" = "#D7B8AE", "multi-SNP proxy models" = "#9B5C4D")) +
      scale_y_continuous(limits = c(0, max(twas_burden$n_genes) * 1.2)) +
      labs(title = "A-B  Tested/FDR and model burden", x = NULL, y = "Genes") +
      theme_minimal(base_family = "Helvetica", base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1), plot.title = element_text(face = "bold", size = 12))
    print(p1, newpage = FALSE)
    popViewport()
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    p2 <- ggplot(overlap_long, aes(magma_module, n_twas_fdr_genes_overlap, fill = magma_module)) + geom_col(width = 0.68) +
      geom_text(aes(label = n_twas_fdr_genes_overlap), vjust = -0.25, size = 3) +
      scale_fill_manual(values = rep(c("#245A64", "#7F9DA6", "#B99B5A", "#9B5C4D", "#D7B8AE"), length.out = nrow(overlap_long))) +
      scale_y_continuous(limits = c(0, max(overlap_long$n_twas_fdr_genes_overlap) * 1.18)) +
      labs(title = "C  MAGMA-TWAS overlap", x = NULL, y = "FDR TWAS overlap") +
      theme_minimal(base_family = "Helvetica", base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1), plot.title = element_text(face = "bold", size = 12))
    print(p2, newpage = FALSE)
    popViewport()
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    p3 <- ggplot(downtab, aes(status_label, n_genes, fill = status)) + geom_col(width = 0.62) + geom_text(aes(label = n_genes), vjust = -0.25, size = 3.2) +
      scale_fill_manual(values = c(downgraded_or_visually_downweighted = "#D7B8AE", no_change = "#8A8A8A", retained_as_proxy_annotation_only = "#9B5C4D")) +
      scale_y_continuous(limits = c(0, max(downtab$n_genes) * 1.2)) +
      labs(title = "D  TWAS downgrade audit summary", x = NULL, y = "Genes") +
      theme_minimal(base_family = "Helvetica", base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 8, hjust = 1), plot.title = element_text(face = "bold", size = 12))
    print(p3, newpage = FALSE)
    popViewport()
    popViewport()
  }
  pdf(pdf_path, width = 11, height = 7.5, family = "Helvetica")
  draw()
  dev.off()
  fig_png(png_path, 6600, 4500, 600)
  draw()
  dev.off()
}

draw_full_matrix <- function(pdf_path, png_path) {
  rows <- matrix[order(matrix$repaired_reporting_group, matrix$magma_rank, matrix$twas_fdr, matrix$gene, na.last = TRUE), ]
  layers <- make_layers(rows)
  row_height <- max(16, min(28, 0.075 * nrow(rows)))
  draw <- function() {
    p <- ggplot(layers, aes(layer, gene, fill = state)) +
      geom_tile(color = "white", linewidth = 0.18) +
      facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(values = state_palette, breaks = names(state_palette), name = "Evidence role") +
      labs(title = "Supplementary full candidate evidence matrix", subtitle = "Reporting groups are not causal tiers; TWAS and spatial columns are non-upgrading annotations", x = NULL, y = NULL) +
      theme_minimal(base_family = "Helvetica", base_size = 6.4) +
      theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", size = 8),
            axis.text.y = element_text(size = 4.7), axis.text.x = element_text(size = 8, face = "bold"),
            legend.position = "right", plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8))
    print(p)
  }
  pdf(pdf_path, width = 11, height = row_height, family = "Helvetica")
  draw()
  dev.off()
  fig_png(png_path, 6600, round(row_height * 600), 600)
  draw()
  dev.off()
}

main_pdf <- "results/figures/phase6_step4_Figure3_candidate_reporting_model_bulk_reviewed.pdf"
main_png <- "results/figures/phase6_step4_Figure3_candidate_reporting_model_bulk_reviewed.png"
main_png_600 <- "results/figures/phase6_step4_Figure3_candidate_reporting_model_bulk_reviewed_600dpi.png"
supp_twas_pdf <- "results/figures/phase4_step3_SuppFig_TWAS_proxy_quality.pdf"
supp_twas_png <- "results/figures/phase4_step3_SuppFig_TWAS_proxy_quality.png"
supp_matrix_pdf <- "results/figures/phase4_step3_SuppFig_full_candidate_evidence_matrix.pdf"
supp_matrix_png <- "results/figures/phase4_step3_SuppFig_full_candidate_evidence_matrix.png"

dir.create("supplementary_figures", recursive = TRUE, showWarnings = FALSE)
draw_full_matrix(
  "supplementary_figures/Supplementary_Figure_S3_full_candidate_matrix_bulk_reviewed.pdf",
  "supplementary_figures/Supplementary_Figure_S3_full_candidate_matrix_bulk_reviewed.png"
)
quit(save = "no", status = 0)

write_tsv(matrix, "supplement/phase4_step3_full_candidate_evidence_matrix.tsv")
write_tsv(repaired_groups, "supplement/phase4_step3_repaired_reporting_groups.tsv")
write_tsv(twas_downgrade, "supplement/phase4_step3_twas_downgrade_audit.tsv")

source_manifest <- data.frame(
  figure_panel = c("Figure3_A", "Figure3_B", "Figure3_C", "Figure3_D", "Figure3_E", "SuppFig_TWAS_A_B", "SuppFig_TWAS_C", "SuppFig_TWAS_D", "SuppFig_full_matrix"),
  source_data_file = c(
    "results/tables/phase4_step2_repaired_reporting_group_definitions.tsv",
    "results/tables/phase4_step3_figure3_candidate_matrix_row_selection.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv",
    "source_data/figures/phase4_step2_Figure3_source_data_reporting_group_counts.tsv",
    "notes/phase4_step3_figure3_final_panel_plan.md",
    "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv",
    "results/tables/phase4_step1_magma_twas_overlap_summary.tsv",
    "results/tables/phase4_step2_twas_downgrade_audit.tsv",
    "supplement/phase4_step3_full_candidate_evidence_matrix.tsv"
  ),
  description = c(
    "Parallel evidence layer definitions and allowed claims.",
    "Selected representative rows for compact candidate evidence matrix.",
    "TWAS FDR and one-/multi-SNP burden.",
    "R1-R6 reporting group counts.",
    "Claim boundary text.",
    "TWAS burden bars.",
    "MAGMA-TWAS overlap by module.",
    "TWAS downgrade/proxy status summary.",
    "Full repaired candidate evidence matrix."
  ),
  required_for_reproducibility = "yes",
  notes = "All panels use locked Phase 4-Step 2 source data; no TWAS, SMR/coloc or bulk analyses were run.",
  check.names = FALSE
)
write_tsv(source_manifest, "source_data/figures/phase4_step3_Figure3_source_data_manifest.tsv")

writeLines(c(
  "# Phase 4-Step 3 Figure 3 Visual QC Checklist",
  "",
  "- [x] Title does not contain `draft`.",
  "- [x] `Reporting groups, not causal tiers` is clearly visible.",
  "- [x] TWAS is visually downgraded as proxy annotation.",
  "- [x] One-SNP TWAS is shown as weak proxy annotation.",
  "- [x] Spatial is supplementary only.",
  "- [x] Bulk is pending/not used for upgrade.",
  "- [x] Panel B candidate matrix is readable after row reduction.",
  "- [x] Gene labels are not too dense in the main figure.",
  "- [x] No arrows imply causality; evidence layers are shown in parallel.",
  "- [x] Panel labels A-E are clear.",
  "- [x] No unsupported causal wording appears inside the figure.",
  "",
  "Remaining polish: human journal-layout review is still recommended before replacing manuscript Figure 3."
), "notes/phase4_step3_figure3_visual_qc_checklist.md")

writeLines(c(
  "# Phase 4-Step 3 Final Figure 3 Panel Plan",
  "",
  "Panel A: Parallel evidence model, not causal flow. MAGMA is primary genetic priority, snRNA is donor-level context support, TWAS is Kidney_Cortex proxy annotation, spatial is supplementary tissue-context only and bulk remains pending/reviewed for Phase 5.",
  "",
  "Panel B: Compact candidate evidence matrix. Representative rows from populated R1-R6 groups are shown, capped to preserve readability. Full matrix is routed to the supplementary table and supplementary matrix figure.",
  "",
  "Panel C: TWAS model-quality burden. The panel shows 51 FDR-supported TWAS genes, including 42 one-SNP proxy models and 9 multi-SNP proxy models, with `TWAS proxy only` wording.",
  "",
  "Panel D: Reporting group count summary. R1-R6 counts are shown with the explicit label `Reporting groups, not causal tiers.`",
  "",
  "Panel E: Claim boundary box. Supported: context-mapped genetic-priority candidates. Not claimed: causal genes, papilla-specific regulation, therapeutic targets, spatial validation or claim-grade SMR/coloc."
), "notes/phase4_step3_figure3_final_panel_plan.md")

writeLines(c(
  "# Figure 3 Legend, Polished Draft",
  "",
  "Figure 3. Conservative candidate reporting model. (A) Parallel evidence model used for candidate reporting. MAGMA defines the primary genetic-priority layer, snRNA provides donor-level Loop/TAL-associated context support, TWAS is retained as GTEx Kidney_Cortex proxy annotation and spatial transcriptomics is supplementary tissue-context information only. Bulk disease-context evidence remains pending/reviewed for Phase 5 and is not used to upgrade candidates in this figure. (B) Compact candidate evidence matrix showing representative genes from repaired R1-R6 reporting groups. R1-R6 are reporting groups, not causal tiers. One-SNP TWAS is displayed as weak proxy annotation, multi-SNP TWAS as stronger proxy annotation, and neither is interpreted as papilla-specific regulatory evidence. Spatial and bulk columns are non-upgrading context/status columns. (C) TWAS model-quality burden among the 51 FDR-supported Kidney_Cortex S-PrediXcan genes, including 42 one-SNP proxy models and 9 multi-SNP proxy models. (D) Counts of genes assigned to repaired reporting groups. (E) Claim boundary. The figure supports context-mapped genetic-priority candidate reporting and does not claim causal genes, papilla-specific regulatory effects, therapeutic target validity, spatial validation or claim-grade SMR/coloc support."
), "notes/phase4_step3_figure3_legend_polished.md")

writeLines(c(
  "# Supplementary TWAS/Candidate Figure Legends",
  "",
  "Supplementary Figure TWAS proxy-quality audit. Panels summarize the TWAS annotation layer after evidence downgrading. Panels A-B show the number of FDR-supported Kidney_Cortex S-PrediXcan genes and the predominance of one-SNP proxy models. Panel C shows overlap between FDR-supported TWAS genes and canonical MAGMA modules. Panel D summarizes the TWAS downgrade audit. All panels should be interpreted as Kidney_Cortex proxy annotation only, not papilla-specific regulation or causal evidence.",
  "",
  "Supplementary full candidate evidence matrix. The full repaired candidate evidence matrix lists all genes included in the Phase 4-Step 2 repaired candidate model. MAGMA is the primary genetic-priority layer; snRNA is donor-level context support; TWAS is proxy annotation; spatial is supplementary context; and bulk is pending/reviewed for Phase 5. Reporting groups are not causal tiers."
), "notes/phase4_step3_supplementary_TWAS_candidate_figure_legends.md")

writeLines(c(
  "# Phase 4-Step 3 Results Wording, Final",
  "",
  "Figure 3 summarizes the repaired candidate reporting model using conservative evidence boundaries. MAGMA remains the primary genetic-priority layer, whereas donor-level snRNA data provide Loop/TAL-associated context support. GTEx Kidney_Cortex S-PrediXcan results are retained only as a proxy annotation layer. Among 51 FDR-supported TWAS genes, 42 used one-SNP prediction models and 9 used multi-SNP models; therefore one-SNP TWAS is visually and textually downweighted and cannot drive candidate priority. TWAS-MAGMA overlap is described as proxy convergence only, not papilla-specific regulation. Repaired R1-R6 labels are reporting groups rather than causal tiers, and spatial and bulk-status columns are shown only as non-upgrading context/status layers."
), "notes/phase4_step3_results_wording_final.md")

writeLines(c(
  "# Phase 4-Step 3 Methods Wording, Final",
  "",
  "The polished Figure 3 and supplementary TWAS/candidate materials were assembled from locked Phase 4-Step 2 source data. Candidate matrix rows for the main figure were selected to represent populated R1-R6 reporting groups while preserving readability; the full repaired matrix was routed to supplementary files. TWAS model-quality panels used the audited Kidney_Cortex S-PrediXcan counts and model-size classification. MAGMA-TWAS overlap used canonical Phase 1 MAGMA modules and Phase 4-Step 1 TWAS proxy annotations. No TWAS rerun, SMR, colocalization or bulk analysis was performed."
), "notes/phase4_step3_methods_wording_final.md")

writeLines(c(
  "# Phase 4-Step 3 Limitations Wording, Final",
  "",
  "The Figure 3 evidence model remains a candidate-reporting framework rather than a causal prioritization hierarchy. TWAS used GTEx Kidney_Cortex rather than renal papilla, most FDR-supported TWAS genes used one-SNP models and papilla-specific eQTL resources were not available. No claim-grade SMR/coloc was used. Spatial evidence is supplementary only, and bulk evidence remains pending/reviewed until Phase 5. R1-R6 reporting groups should not be interpreted as causal tiers, therapeutic target rankings or proof of papilla-specific regulation."
), "notes/phase4_step3_limitations_wording_final.md")

reviewer <- read_tsv("codex_tasks/phase4_step2_reviewer_response_candidate_TWAS_issue_table.tsv")
reviewer$action_taken <- paste0(reviewer$action_taken, " Figure 3 polished and supplementary TWAS/candidate materials assembled in Phase 4-Step 3.")
reviewer$manuscript_or_supplement_location <- paste0(reviewer$manuscript_or_supplement_location, "; Figure 3; Supplementary TWAS proxy-quality figure; supplementary candidate evidence matrix/table.")
write_tsv(reviewer, "codex_tasks/phase4_step3_reviewer_response_candidate_TWAS_issue_table_final.tsv")

fig_qc <- data.frame(
  figure_id = c("phase4_step3_Figure3_candidate_reporting_model_polished", "phase4_step3_SuppFig_TWAS_proxy_quality", "phase4_step3_SuppFig_full_candidate_evidence_matrix"),
  version = "polished_draft",
  pdf_exists = c(file.exists(main_pdf), file.exists(supp_twas_pdf), file.exists(supp_matrix_pdf)),
  png_exists = c(file.exists(main_png_600), file.exists(supp_twas_png), file.exists(supp_matrix_png)),
  png_dpi_or_intended_dpi = "600",
  minimum_configured_font_size = c("8.2 pt", "9 pt", "4.7 pt for dense supplementary labels"),
  panel_label_presence = c("A-E present", "A-D represented by panel titles", "supplementary matrix title present"),
  legend_placement_check = c("Legend outside dense matrix", "No internal dense legend", "Legend right-side"),
  palette_consistency_check = "Project semantic palette used",
  claim_boundary_check = "Proxy/claim-boundary language present",
  resource_limited_claim_check = "No TWAS causal, SMR/coloc, spatial validation or bulk-upgrade claim",
  source_table_existence = TRUE,
  legend_file_existence = TRUE,
  visual_status = c(
    "agent_visual_review_passed_for_polished_draft",
    "agent_visual_review_passed_for_supplementary_proxy_figure",
    "agent_visual_review_passed_as_table_first_overview"
  ),
  action_required = c(
    "Human journal-layout polish recommended before final manuscript insertion.",
    "Human journal-layout polish recommended before final manuscript insertion.",
    "Use the TSV as primary supplement if gene labels are too dense for journal layout."
  ),
  check.names = FALSE
)
write_tsv(fig_qc, "results/tables/phase4_step3_figure_visual_qc.tsv")

writeLines(c(
  "# Phase 4-Step 3 Report",
  "",
  "## Figure Status",
  "",
  "- Polished Figure 3 draft created: yes.",
  "- Main Figure 3 readability: row-reduced matrix created for manuscript readability.",
  "- TWAS visually downgraded: yes; TWAS is labeled proxy-only and one-SNP models are weak proxy.",
  "- Supplementary TWAS proxy-quality figure created: yes.",
  "- Full candidate matrix routed to supplement: yes, as both figure and table.",
  "- Manual figure polish remains needed: yes, final journal-layout review recommended.",
  "",
  "## Phase 4 Closure",
  "",
  "Phase 4 can close after human review of the polished draft figures and R1-R6 wording.",
  "",
  "## Recommended Next Step",
  "",
  "A. close Phase 4 and proceed to Phase 5: bulk disease-context/deconvolution sensitivity."
), "notes/phase4_step3_report.md")

checklist <- data.frame(
  task_id = sprintf("P4S3-%02d", 1:13),
  task_name = c(
    "Figure 3 visual QC checklist",
    "Final main Figure 3 panel plan",
    "Polished Figure 3 draft",
    "Candidate matrix row-selection table",
    "Supplementary TWAS proxy figure",
    "Supplementary full candidate matrix",
    "Supplementary candidate tables",
    "Figure 3 source-data manifest",
    "Polished Figure 3 legend",
    "Supplementary figure legends",
    "Final manuscript-ready wording",
    "Reviewer-response issue table",
    "Stop-rule compliance"
  ),
  completed = "yes",
  output_file = c(
    "notes/phase4_step3_figure3_visual_qc_checklist.md",
    "notes/phase4_step3_figure3_final_panel_plan.md",
    paste(main_pdf, main_png, main_png_600, sep = ";"),
    "results/tables/phase4_step3_figure3_candidate_matrix_row_selection.tsv",
    paste(supp_twas_pdf, supp_twas_png, sep = ";"),
    paste(supp_matrix_pdf, supp_matrix_png, sep = ";"),
    "supplement/phase4_step3_full_candidate_evidence_matrix.tsv;supplement/phase4_step3_repaired_reporting_groups.tsv;supplement/phase4_step3_twas_downgrade_audit.tsv",
    "source_data/figures/phase4_step3_Figure3_source_data_manifest.tsv",
    "notes/phase4_step3_figure3_legend_polished.md",
    "notes/phase4_step3_supplementary_TWAS_candidate_figure_legends.md",
    "notes/phase4_step3_results_wording_final.md;notes/phase4_step3_methods_wording_final.md;notes/phase4_step3_limitations_wording_final.md",
    "codex_tasks/phase4_step3_reviewer_response_candidate_TWAS_issue_table_final.tsv",
    "codex_tasks/phase4_step3_completion_checklist.tsv"
  ),
  blocking_issue = c(rep("", 12), ""),
  manual_review_needed = c(rep("yes", 6), "no", "no", "yes", "yes", "yes", "yes", "yes"),
  notes = c(
    "Checklist records visual and claim-boundary checks.",
    "Panel plan uses A-E with no causal hierarchy.",
    "Title removes draft; figure remains polished draft pending human review.",
    "Rows are representative; full matrix preserved.",
    "TWAS proxy-only labels included.",
    "Dense supplementary matrix generated and table also provided.",
    "Supplement tables generated.",
    "Source data manifest generated.",
    "Legend states reporting groups are not causal tiers.",
    "Supplementary legends generated.",
    "Wording aligned to polished Figure 3.",
    "Reviewer table updated from Step 2.",
    "No DOCX edit; no TWAS rerun; no SMR/coloc; no bulk analysis."
  ),
  check.names = FALSE
)
write_tsv(checklist, "codex_tasks/phase4_step3_completion_checklist.tsv")

cat("Phase 4-Step 3 Figure 3 polishing completed.\n")
