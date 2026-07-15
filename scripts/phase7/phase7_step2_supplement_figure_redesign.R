#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

root <- normalizePath(getwd())
out_dir <- file.path(root, "supplementary_figures_v3.3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL <- c(
  teal = "#205B66", deep = "#123F4A", bluegrey = "#7897A1",
  gold = "#B5964A", terra = "#A55846", rose = "#D5B2A7",
  olive = "#697A5B", ink = "#273238", grey = "#68757B",
  pale = "#E8EDEE", light = "#F5F7F7", white = "#FFFFFF"
)

theme_journal <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(color = COL["ink"]),
      plot.title = element_text(face = "bold", size = base_size + 1.5),
      plot.subtitle = element_text(size = base_size, color = COL["grey"]),
      plot.caption = element_text(size = base_size - 0.5, color = COL["grey"], hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, color = COL["ink"]),
      strip.background = element_rect(fill = COL["pale"], color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.2),
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      panel.grid.major = element_line(color = "#E8ECEC", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 6, 5, 6)
    )
}

read_tsv <- function(path) fread(path, sep = "\t", na.strings = c("NA", ""), data.table = TRUE)
read_gz <- function(path) fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", na.strings = c("NA", ""))
pdf_device <- function(filename, ...) grDevices::pdf(filename, ..., family = "Helvetica", useDingbats = FALSE, onefile = TRUE)

save_page <- function(plot, pdf_path, png_path, width = 11, height = 8.5) {
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = pdf_device, bg = "white")
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

save_multipage <- function(plots, path, width = 11, height = 8.5) {
  grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
  for (p in plots) print(p)
  dev.off()
}

panel_title <- function(letter, title) paste0(letter, "  ", title)

# Supplementary Figure S2 ---------------------------------------------------
spots <- read_tsv(file.path(root, "source_data/figures/phase3_step2_label_transfer_qc_source_data.tsv"))
composition <- read_tsv(file.path(root, "results/tables/phase3_step2B_predicted_compartment_composition.tsv"))
loop <- read_tsv(file.path(root, "results/tables/phase3_step2B_loop_tal_projection_usability.tsv"))
spatial_manifest <- read_tsv(file.path(root, "results/tables/phase3_step1_spatial_sample_manifest.tsv"))
module_spots <- read_gz(file.path(root, "results/tables/phase3_step3_spatial_module_scores.tsv.gz"))
module_summary <- read_tsv(file.path(root, "results/tables/phase3_step3_predicted_compartment_module_summary.tsv"))

section_map <- unique(spatial_manifest[, .(sample_id, section_id, dataset_label)])
section_map[, dataset_label := factor(dataset_label, levels = c("GSE206306", "GSE231630"))]
setorder(section_map, dataset_label, sample_id, section_id)
section_map[, section_short := sprintf("S%02d", seq_len(.N))]
section_map[, section_facet := paste0(section_short, " | ", dataset_label)]

spots <- merge(spots, section_map, by = c("sample_id", "section_id"), all.x = TRUE, suffixes = c("", ".locked"))
composition <- merge(composition, section_map, by = c("sample_id", "section_id"), all.x = TRUE)
loop <- merge(loop, section_map, by = c("sample_id", "section_id"), all.x = TRUE)
module_spots <- merge(module_spots, section_map, by = c("sample_id", "section_id"), all.x = TRUE)
module_summary <- merge(module_summary, section_map, by = c("sample_id", "section_id"), all.x = TRUE)

section_levels <- section_map$section_short
facet_levels <- section_map$section_facet
spots[, `:=`(section_short = factor(section_short, levels = section_levels), section_facet = factor(section_facet, levels = facet_levels))]
composition[, section_short := factor(section_short, levels = rev(section_levels))]
loop[, section_short := factor(section_short, levels = rev(section_levels))]
module_spots[, section_facet := factor(section_facet, levels = facet_levels)]

resource <- data.table(
  metric = factor(c("Complete sections", "GSE206306", "GSE231630", "ROI annotation"),
                  levels = c("Complete sections", "GSE206306", "GSE231630", "ROI annotation")),
  value = c("10", "4", "6", "None")
)
p_s2a_a <- ggplot(resource, aes(metric, 1)) +
  geom_tile(width = 0.88, height = 0.66, fill = "white", color = COL["bluegrey"], linewidth = 0.5) +
  geom_text(aes(label = value), y = 1.12, fontface = "bold", size = 5.1, color = COL["deep"]) +
  geom_text(aes(label = metric), y = 0.86, size = 2.8, color = COL["ink"]) +
  coord_cartesian(ylim = c(0.62, 1.38), clip = "off") +
  labs(title = panel_title("A", "spatial resource"), x = NULL, y = NULL) +
  theme_void(base_family = "Helvetica") + theme(plot.title = element_text(face = "bold", size = 11))

score_data <- spots[in_tissue == 1 & is.finite(prediction_score_max)]
p_s2a_b <- ggplot(score_data, aes(prediction_score_max)) +
  geom_histogram(bins = 32, fill = COL["bluegrey"], color = "white", linewidth = 0.15) +
  facet_wrap(~section_facet, ncol = 5, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = panel_title("B", "maximum label-transfer score"), x = "Maximum prediction score", y = "In-tissue spots") +
  theme_journal(8) + theme(strip.text = element_text(size = 7.1), axis.text = element_text(size = 6.8))

comp_cols <- c(
  "Collecting duct" = "#6F98D9", "Injured/undifferentiated epithelial" = "#E57F73",
  "Fibroblast/stromal" = "#68B66C", "Endothelial" = COL["gold"],
  "Pericyte/smooth muscle" = COL["olive"], "Loop/TAL" = COL["deep"]
)
comp_cols <- setNames(unname(comp_cols), names(comp_cols))
p_s2a_c <- ggplot(composition, aes(percent_spots, section_short, fill = predicted_compartment)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = comp_cols, name = "Predicted broad compartment") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = label_percent(scale = 1)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = panel_title("C", "max-predicted broad-compartment composition"), x = "In-tissue spots", y = NULL) +
  theme_journal(8.3) + theme(legend.position = "bottom", legend.text = element_text(size = 7), legend.key.width = unit(0.38, "cm")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

p_s2a_d <- ggplot(loop, aes(percent_nonzero_Loop_TAL_score_spots, section_short)) +
  geom_vline(xintercept = 0, color = COL["grey"], linewidth = 0.4) +
  geom_point(size = 2.4, shape = 21, fill = COL["terra"], color = COL["ink"], stroke = 0.3) +
  geom_text(aes(label = "0%"), hjust = -0.45, size = 2.8, color = COL["terra"]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = function(x) paste0(x, "%")) +
  labs(title = panel_title("D", "Loop/TAL projection usability"), subtitle = "Nonzero prediction-score fraction was 0% in every section", x = "Nonzero Loop/TAL score spots", y = NULL) +
  theme_journal(8.3)

s2a <- ((p_s2a_a | p_s2a_b) / (p_s2a_c | p_s2a_d)) +
  plot_layout(widths = c(0.9, 1.35), heights = c(1, 1.15)) +
  plot_annotation(
    title = "Supplementary Figure S2A | Spatial projection quality control",
    subtitle = "Ten complete sections: four GSE206306 and six GSE231630",
    caption = "Broad-compartment tissue-context projection only; no lesion ROI was available.",
    theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 15, color = COL["ink"]),
                  plot.subtitle = element_text(family = "Helvetica", size = 10, color = COL["grey"]),
                  plot.caption = element_text(family = "Helvetica", size = 8.5, color = COL["grey"], hjust = 0))
  )

save_page(
  s2a,
  file.path(out_dir, "Supplementary_Figure_S2A_spatial_QC_polished.pdf"),
  file.path(out_dir, "Supplementary_Figure_S2A_spatial_QC_polished.png")
)

module_labels <- c(
  MAGMA_Bonferroni = "Bonferroni", MAGMA_top100 = "Top 100", MAGMA_suggestive_p1e4 = "Suggestive",
  MAGMA_top50 = "Top 50", MAGMA_FDR05 = "FDR 0.05"
)
overlay_modules <- c("MAGMA_Bonferroni", "MAGMA_top100", "MAGMA_suggestive_p1e4")
overlay_pages <- lapply(seq_along(overlay_modules), function(i) {
  mod <- overlay_modules[i]
  d <- module_spots[module_name == mod & in_tissue == 1 & is.finite(module_score)]
  lim <- range(d$module_score, finite = TRUE)
  ggplot(d, aes(x_coord, -y_coord, color = module_score)) +
    geom_point(size = 0.26, alpha = 0.92, stroke = 0) +
    facet_wrap(~section_facet, ncol = 5, scales = "free") +
    scale_color_gradient(low = "#EEF2F2", high = COL["deep"], limits = lim, name = "Module score") +
    labs(
      title = paste0("Supplementary Figure S2B", i, " | ", module_labels[mod], " spatial module overlay"),
      subtitle = "Ten sections: four GSE206306 and six GSE231630",
      caption = "Descriptive tissue-context projection; no plaque localization, disease-control test or Loop/TAL enrichment.",
      x = NULL, y = NULL
    ) +
    theme_void(base_family = "Helvetica") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold", size = 15, color = COL["ink"]),
      plot.subtitle = element_text(size = 10, color = COL["grey"]),
      plot.caption = element_text(size = 8.5, color = COL["grey"], hjust = 0),
      strip.background = element_rect(fill = COL["pale"], color = NA),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "bottom", legend.title = element_text(face = "bold", size = 8.5), legend.text = element_text(size = 8),
      plot.margin = margin(12, 14, 10, 14)
    )
})

module_summary[, module_display := factor(module_labels[module_name], levels = c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive"))]
module_summary[, predicted_compartment_short := factor(
  predicted_compartment,
  levels = c("Collecting duct", "Injured/undifferentiated epithelial", "Fibroblast/stromal", "Endothelial", "Pericyte/smooth muscle"),
  labels = c("Collecting duct", "Injured epithelial", "Fibroblast/stromal", "Endothelial", "Pericyte/smooth muscle")
)]
summary_page <- ggplot(module_summary, aes(predicted_compartment_short, median_module_score, color = dataset_label, group = dataset_label)) +
  geom_boxplot(aes(group = predicted_compartment_short), color = COL["grey"], fill = COL["light"], outlier.shape = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.45, alpha = 0.85) +
  facet_wrap(~module_display, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c(GSE206306 = unname(COL["teal"]), GSE231630 = unname(COL["terra"])), name = "Dataset") +
  labs(
    title = "Supplementary Figure S2B4 | Module scores by predicted broad compartment",
    subtitle = "Section-level descriptive summaries",
    caption = "Transferred labels provide broad context only and are not cell fractions or histological ROIs.",
    x = NULL, y = "Median spatial module score"
  ) +
  theme_journal(9) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 7.4), legend.position = "bottom",
        plot.title = element_text(size = 15), plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 8.5))

s2b_pages <- c(overlay_pages, list(summary_page))
save_multipage(s2b_pages, file.path(out_dir, "Supplementary_Figure_S2B_spatial_overlays_polished.pdf"))
for (i in seq_along(s2b_pages)) {
  ggsave(file.path(out_dir, sprintf("Supplementary_Figure_S2B_page%02d_polished.png", i)), s2b_pages[[i]], width = 11, height = 8.5, units = "in", dpi = 600, bg = "white")
}
save_multipage(c(list(s2a), s2b_pages), file.path(out_dir, "Supplementary_Figure_S2_spatial_projection_polished.pdf"))
ggsave(file.path(out_dir, "Supplementary_Figure_S2_spatial_projection_polished.png"), s2a, width = 11, height = 8.5, units = "in", dpi = 600, bg = "white")

# Supplementary Figure S3 ---------------------------------------------------
burden <- read_tsv(file.path(root, "source_data/figures/phase4_step2_Figure3_source_data_twas_burden.tsv"))
overlap <- read_tsv(file.path(root, "results/tables/phase4_step1_magma_twas_overlap_summary.tsv"))
downgrade <- read_tsv(file.path(root, "results/tables/phase4_step2_twas_downgrade_audit.tsv"))
candidate <- read_tsv(file.path(root, "supplementary_tables/Supplementary_Table_5_candidate_reporting_matrix.tsv"))

burden[, category_short := factor(
  category,
  levels = c("FDR-supported TWAS genes", "one-SNP proxy models", "multi-SNP proxy models"),
  labels = c("FDR supported", "One-SNP", "Multi-SNP")
)]
p_s3a_a <- ggplot(burden, aes(category_short, n_genes, fill = category_short)) +
  geom_col(width = 0.66) + geom_text(aes(label = n_genes), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("FDR supported" = unname(COL["bluegrey"]), "One-SNP" = unname(COL["rose"]), "Multi-SNP" = unname(COL["terra"])), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = panel_title("A", "Kidney_Cortex model burden"), x = NULL, y = "Genes") + theme_journal(9)

overlap[, module_display := factor(
  magma_module,
  levels = c("MAGMA_top50", "MAGMA_top100", "MAGMA_Bonferroni", "MAGMA_FDR05", "MAGMA_suggestive"),
  labels = c("Top 50", "Top 100", "Bonferroni", "FDR 0.05", "Suggestive")
)]
p_s3a_b <- ggplot(overlap, aes(module_display, n_twas_fdr_genes_overlap, fill = module_display)) +
  geom_col(width = 0.68) + geom_text(aes(label = n_twas_fdr_genes_overlap), vjust = -0.3, size = 3.0) +
  scale_fill_manual(values = c("Top 50" = unname(COL["teal"]), "Top 100" = unname(COL["bluegrey"]), "Bonferroni" = unname(COL["gold"]), "FDR 0.05" = unname(COL["terra"]), "Suggestive" = unname(COL["rose"])), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = panel_title("B", "MAGMA-TWAS overlap"), x = NULL, y = "FDR-supported TWAS genes") +
  theme_journal(9) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

down_counts <- downgrade[, .(n_genes = .N), by = downgrade_or_no_change]
down_counts[, status_label := factor(
  downgrade_or_no_change,
  levels = c("downgraded_or_visually_downweighted", "no_change", "retained_as_proxy_annotation_only"),
  labels = c("One-SNP\ndownweighted", "No group\nchange", "Multi-SNP\nproxy only")
)]
p_s3a_c <- ggplot(down_counts, aes(status_label, n_genes, fill = status_label)) +
  geom_col(width = 0.62) + geom_text(aes(label = n_genes), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("One-SNP\ndownweighted" = unname(COL["rose"]), "No group\nchange" = "#8A8A8A", "Multi-SNP\nproxy only" = unname(COL["terra"])), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = panel_title("C", "proxy-evidence handling"), x = NULL, y = "Genes") + theme_journal(9)

s3a <- ((p_s3a_a | p_s3a_b) / p_s3a_c) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title = "Supplementary Figure S3A | Kidney_Cortex TWAS proxy quality",
    subtitle = "Proxy annotation only; model support does not establish papilla-specific regulation",
    caption = "Counts and overlaps are reproduced from the locked Phase 4 tables.",
    theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 15, color = COL["ink"]),
                  plot.subtitle = element_text(family = "Helvetica", size = 10, color = COL["grey"]),
                  plot.caption = element_text(family = "Helvetica", size = 8.5, color = COL["grey"], hjust = 0))
  )
save_page(
  s3a,
  file.path(out_dir, "Supplementary_Figure_S3A_TWAS_proxy_quality_polished.pdf"),
  file.path(out_dir, "Supplementary_Figure_S3A_TWAS_proxy_quality_polished.png")
)

state_cols <- setNames(
  c(unname(COL["teal"]), unname(COL["deep"]), unname(COL["gold"]), unname(COL["terra"]), unname(COL["rose"]), "#B9C6CA", unname(COL["pale"]), "#F4F6F6"),
  c("genetic_priority", "context", "partial", "proxy_multi", "proxy_one", "supplementary", "reviewed", "none")
)

make_layers <- function(rows) {
  z <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i]
    data.table(
      gene = r$gene,
      group = r$repaired_reporting_group,
      layer = c("MAGMA", "snRNA", "TWAS proxy", "Spatial", "Bulk review"),
      state = c(
        ifelse(r$magma_status %in% c("top50", "top100", "bonferroni", "fdr05"), "genetic_priority", ifelse(r$magma_status == "suggestive", "partial", "none")),
        ifelse(r$snRNA_support_status == "strong_context_support", "context", ifelse(r$snRNA_support_status %in% c("moderate_context_support", "partial_context_support"), "partial", "none")),
        ifelse(r$twas_model_type == "multi_snp_proxy", "proxy_multi", ifelse(r$twas_model_type == "one_snp_proxy", "proxy_one", "none")),
        "supplementary", "reviewed"
      )
    )
  }))
  z[, `:=`(
    gene = factor(gene, levels = rev(rows$gene)),
    layer = factor(layer, levels = c("MAGMA", "snRNA", "TWAS proxy", "Spatial", "Bulk review"))
  )]
  z
}

matrix_pages <- list()
page_meta <- list()
max_rows <- 32L
page_index <- 0L
for (grp in paste0("R", 1:6)) {
  grp_rows <- candidate[repaired_reporting_group == grp]
  if (nrow(grp_rows) == 0) next
  chunks <- split(seq_len(nrow(grp_rows)), ceiling(seq_len(nrow(grp_rows)) / max_rows))
  for (chunk_id in seq_along(chunks)) {
    page_index <- page_index + 1L
    rows <- grp_rows[chunks[[chunk_id]]]
    layers <- make_layers(rows)
    subtitle <- sprintf("%s | rows %d-%d of %d | source order preserved", grp,
                        min(chunks[[chunk_id]]), max(chunks[[chunk_id]]), nrow(grp_rows))
    p <- ggplot(layers, aes(layer, gene, fill = state)) +
      geom_tile(color = "white", linewidth = 0.28) +
      scale_fill_manual(
        values = state_cols, breaks = names(state_cols), name = "Evidence role",
        labels = c("MAGMA priority", "Strong context", "Partial context", "Multi-SNP proxy", "One-SNP proxy", "Supplementary", "Reviewed, non-upgrading", "None")
      ) +
      labs(
        title = sprintf("Supplementary Figure S3B | Candidate evidence matrix, page %d", page_index),
        subtitle = subtitle,
        caption = "Direct rendering of Supplementary Table 5; R1-R6 are reporting groups, not causal tiers.",
        x = NULL, y = NULL
      ) +
      theme_journal(9.2) +
      theme(
        panel.grid = element_blank(), axis.text.y = element_text(size = 7.2, face = "italic"),
        axis.text.x = element_text(size = 9, face = "bold"), legend.position = "right",
        legend.text = element_text(size = 7.5), plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 9.5), plot.caption = element_text(size = 8.5)
      )
    matrix_pages[[page_index]] <- p
    page_meta[[page_index]] <- data.table(page = page_index, reporting_group = grp, chunk = chunk_id,
                                          first_gene = rows$gene[1], last_gene = rows$gene[nrow(rows)], n_rows = nrow(rows))
  }
}

save_multipage(matrix_pages, file.path(out_dir, "Supplementary_Figure_S3B_candidate_matrix_paginated_polished.pdf"))
for (i in seq_along(matrix_pages)) {
  ggsave(file.path(out_dir, sprintf("Supplementary_Figure_S3B_matrix_page%02d_polished.png", i)), matrix_pages[[i]], width = 11, height = 8.5, units = "in", dpi = 600, bg = "white")
}
fwrite(rbindlist(page_meta), file.path(out_dir, "Supplementary_Figure_S3B_pagination_index.tsv"), sep = "\t", quote = FALSE, na = "")

save_multipage(c(list(s3a), matrix_pages), file.path(out_dir, "Supplementary_Figure_S3_TWAS_candidate_proxy_polished.pdf"))
ggsave(file.path(out_dir, "Supplementary_Figure_S3_TWAS_candidate_proxy_polished.png"), s3a, width = 11, height = 8.5, units = "in", dpi = 600, bg = "white")

message("Phase 7-Step 2 Supplementary Figures S2 and S3 created from locked source tables.")
