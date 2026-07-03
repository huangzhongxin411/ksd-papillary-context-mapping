suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(scales)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript", recursive = TRUE, showWarnings = FALSE)
dir.create("results/supervisor_review_package_v1.0", recursive = TRUE, showWarnings = FALSE)

pal <- list(
  genetic = "#245B64",
  scrna = "#6F8F98",
  p1 = "#B59A5B",
  disease = "#9A5F52",
  muted = "#D8D8D8",
  pale = "#EEF3F4",
  text = "#303030",
  grid = "#A8A8A8"
)

theme_pub <- function(base_size = 8.8) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = pal$text),
      plot.subtitle = element_text(size = base_size - 0.6, color = "#555555"),
      axis.title = element_text(color = pal$text),
      axis.text = element_text(color = pal$text),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.6)
    )
}

write_lines <- function(x, path) writeLines(x, path, useBytes = TRUE)

label_cell <- function(x, short = FALSE) {
  full <- c(
    Collecting_duct_principal = "Collecting duct",
    Fibroblast_stromal = "Fibroblast/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epithelial",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivascular/mural-like",
    Pericyte_smooth_muscle = "Perivascular/mural-like"
  )
  short_map <- c(
    Collecting_duct_principal = "CD",
    Fibroblast_stromal = "Fib/stromal",
    Endothelial = "Endothelial",
    Injured_undifferentiated_epithelial = "Injured epi.",
    Loop_of_Henle_TAL = "Loop/TAL",
    Perivascular_mural_like = "Perivasc.",
    Pericyte_smooth_muscle = "Perivasc."
  )
  map <- if (short) short_map else full
  unname(ifelse(x %in% names(map), map[x], x))
}

metric_value <- function(dt, key) {
  v <- dt[metric == key, value]
  if (length(v) == 0) return(NA_character_)
  as.character(v[1])
}

make_manuscript_v10 <- function() {
  sanity <- fread("results/tables/phase1_gwas_sanity_check.tsv")
  qc <- fread("results/tables/phase1_gwas_qc_report.tsv")
  raw_n <- metric_value(sanity, "raw_snp_n")
  clean_n <- metric_value(sanity, "clean_snp_n")
  gws_n <- metric_value(sanity, "genome_wide_sig_snp_n")
  lead_n <- metric_value(sanity, "lead_snp_n")
  loci_n <- metric_value(sanity, "independent_locus_n")
  min_p <- metric_value(sanity, "min_p")
  removed_missing <- metric_value(qc, "removed_missing_required")
  removed_dup <- metric_value(qc, "removed_duplicate_snp")
  if (is.na(removed_dup)) removed_dup <- metric_value(qc, "duplicate_snp_n")
  if (is.na(removed_dup)) removed_dup <- "not separately detected in the final retained QC report"

  title <- "# Post-GWAS mapping of kidney stone risk identifies Loop/TAL-associated renal papillary cellular and disease-context signatures"
  abstract <- c(
    "## Abstract",
    "",
    "**Background:** Kidney stone disease (KSD) has a substantial genetic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.",
    "",
    "**Methods:** We integrated public KSD GWAS summary statistics, MAGMA gene-based prioritization, audited GSE231569 renal papillary single-nucleus annotations, six-gene P1 candidate evidence scoring, and patient-aware GSE73680 papillary plaque/stone disease-context expression analysis.",
    "",
    "**Results:** MAGMA-prioritized KSD genes localized to a Loop/TAL-associated single-nucleus context in GSE231569, with enrichment exceeding random gene-set expectations and remaining supported across robustness checks. The six P1 candidate genes formed an interpretable TAL, epithelial transport and calcium-handling expression spectrum, but did not behave as a uniform disease-validated gene panel. In GSE73680, MAGMA-prioritized modules showed module-level plaque/stone disease-context association, whereas P1 single-gene responses were heterogeneous and not FDR-supported. Functional-context analyses linked prioritized modules to nephron development, TAL transport, calcium/ion handling and papillary injury/remodeling programs.",
    "",
    "**Conclusion:** These findings provide a claim-bounded post-GWAS framework linking KSD genetic prioritization to a Loop/TAL-associated renal papillary cellular context and MAGMA module-level plaque/stone disease-context programs. They nominate biologically interpretable candidate genes and modules for future causal, TWAS, colocalization and spatial validation."
  )
  intro <- c(
    "## Introduction",
    "",
    "Kidney stone disease is a common and recurrent disorder with contributions from urinary chemistry, renal epithelial transport, papillary microenvironments and inherited susceptibility. Genome-wide association studies have identified multiple loci associated with kidney stone risk, but most association signals do not immediately nominate a causal gene, cell type or disease-stage context. Post-GWAS interpretation therefore requires a framework that connects statistical genetic prioritization to renal cell states and disease-relevant expression patterns while preserving clear boundaries around causal inference.",
    "",
    "The renal papilla is a plausible tissue context for stone formation because it contains epithelial, interstitial, vascular and immune compartments involved in urine concentration, mineral handling and papillary plaque biology. The Loop of Henle and thick ascending limb (Loop/TAL) are central to renal electrolyte transport, urinary concentration and calcium/ion handling. These functions make TAL-associated epithelial programs biologically plausible contexts for interpreting genetic signals related to urinary mineral balance and stone susceptibility, but whether KSD risk genes converge on this renal papillary epithelial compartment has not been systematically evaluated using audited single-nucleus and disease-context transcriptomic resources.",
    "",
    "Single-cell and single-nucleus transcriptomic resources provide an opportunity to localize genetic signals to renal cell states, but such analyses depend on careful annotation audit and should not be interpreted as causal mediation on their own. Similarly, disease-context expression datasets can provide support for expression shifts in plaque or stone-associated tissue, but they cannot substitute for colocalization, TWAS or spatial validation.",
    "",
    "Here we integrated MAGMA gene-level prioritization from primary public KSD GWAS summary statistics with audited GSE231569 renal papillary single-nucleus annotations and GSE73680 papillary plaque/stone expression data. The study was designed to ask whether kidney stone genetic risk genes converge on a renal papillary cellular context, whether a core set of candidate genes forms a coherent expression spectrum, and whether an independent disease-context dataset supports the prioritized modules. We explicitly distinguish MAGMA plus snRNA-supported cellular context from genetic causality, TWAS convergence, colocalization and spatial validation."
  )
  methods <- c(
    "## Methods",
    "",
    "### GWAS summary statistics and QC",
    "",
    sprintf("Public KSD GWAS summary statistics were inspected for required variant, genomic position, allele, effect-size and P-value fields before downstream analysis. The raw file contained %s rows, of which %s SNP records were retained after initial GWAS QC; records with missing required fields (%s in the retained QC report), invalid P values or invalid genomic positions were excluded before lead-locus construction. Duplicate SNP handling was audited during harmonization (%s). QQ and Manhattan plots were generated for visual QC, and genome-wide significant variants were distance-pruned into %s lead SNPs and %s independent loci from %s genome-wide significant SNPs. The minimum retained P value was %s. Harmonized TWAS-ready input files were generated as an analysis scaffold, but TWAS association testing was not used as a manuscript evidence layer because external expression weights and LD/model resources were incomplete.", raw_n, clean_n, removed_missing, removed_dup, lead_n, loci_n, gws_n, min_p),
    "",
    "### MAGMA gene-based prioritization",
    "",
    "MAGMA v1.10 was used for gene-based prioritization with a GRCh37/hg19-compatible reference. Gene-level results were ranked by MAGMA P value, and downstream gene sets were defined from fixed ranking and significance thresholds, including MAGMA top 50, MAGMA top 100, FDR-significant and suggestive gene sets. These MAGMA-derived sets formed the main genetic prioritization layer for all later analyses.",
    "",
    "### GSE231569 single-nucleus processing and annotation audit",
    "",
    "Renal papillary single-nucleus data from GSE231569 were processed with an audited broad cell-type annotation. Marker expression, cluster assignment and renal epithelial compartment labels were reviewed before gene-set projection. Loop of Henle/thick ascending limb cells were harmonized as the Loop/TAL compartment and treated as an audited epithelial transport context. Cell types with review flags or low abundance were retained with explicit caution labels rather than being used as unqualified primary localization claims.",
    "",
    "### Gene-set projection to single-nucleus cell types",
    "",
    "MAGMA-prioritized gene sets were projected onto audited GSE231569 renal papillary cell types using module-score summaries across detected genes. The main projected sets were MAGMA top 50, MAGMA top 100, MAGMA FDR and MAGMA suggestive genes. For each cell type and gene set, the observed expression-context score was compared with size-matched random gene sets drawn from the detected single-nucleus background. Benchmark percentiles were used to identify cell-type contexts exceeding random expectation.",
    "",
    "Robustness checks were summarized without introducing additional primary claims. Locus-balanced benchmarking evaluated whether the MAGMA top 50 TAL signal persisted after balancing locus contribution. Conservative annotation sensitivity excluded low-confidence or review-flagged contexts where appropriate. Leave-one-locus-out sensitivity recalculated the Loop/TAL benchmark after removing each locus group and used the minimum retained percentile as a robustness summary. These analyses were interpreted as single-nucleus cellular context mapping rather than causal mediation, TWAS convergence, colocalization or spatial validation.",
    "",
    "### P1 candidate gene evidence scoring",
    "",
    "Six P1 candidate genes, UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2, were evaluated across audited GSE231569 cell types. For each gene, the evidence summary included TAL expression rank, average expression, detection frequency, donor-level detection, TAL specificity ratio, TAL program correlation and manuscript role assignment. The purpose was to define a gene-centric evidence spectrum across TAL, epithelial transport and calcium-handling biology, not to force all P1 genes into a uniform TAL marker class.",
    "",
    "### GSE73680 reconstruction from Agilent feature files",
    "",
    "GSE73680 was reconstructed from the local GSE73680_RAW archive. The archive yielded 62 TXT.gz files, all of which passed gzip validation and Agilent text-structure checks. Agilent FEATURES blocks were parsed into a 44,661-feature by 62-sample expression matrix. Direct gene-symbol mapping assigned 32,055 mapped genes, and curated metadata linked expression profiles to patient and tissue-context labels. After inclusion filtering, the disease-context analysis used 55 samples from 29 patients, including 26 paired patients contributing both control/adjacent and plaque/stone papilla samples.",
    "",
    "### Patient-aware disease-context analysis",
    "",
    "Primary GSE73680 disease-context analyses used limma duplicateCorrelation to account for repeated patient-level sampling. Paired patient-level delta analyses were used as sensitivity checks because 26 patients contributed both analysis groups. P1 single-gene responses were analyzed separately from module-level responses. Module scores were computed as mean z-scores across detected genes for P1 core candidates, MAGMA top 50, MAGMA top 100, MAGMA FDR, MAGMA suggestive, TAL marker and injury/remodeling marker sets.",
    "",
    "### Random gene-set benchmarking",
    "",
    "Size-matched random gene-set benchmarking was used as the main Figure 4 benchmark to evaluate whether observed disease-context module shifts exceeded random expectations for gene sets of comparable size. Leave-one-gene-out analyses, MAGMA-without-P1 sensitivity and paired direction-consistency checks were used to evaluate robustness. Expression-matched random benchmarking was retained as a conservative sensitivity analysis; under the current implementation, it did not provide primary support and was not used as the main benchmark.",
    "",
    "### Cell-level MAGMA module score projection",
    "",
    "An audited GSE231569 Seurat object with UMAP coordinates was used for cell-level visualization. MAGMA top 50, top 100, FDR and suggestive gene sets were scored at the cell level as mean z-scored expression across detected module genes. Donor-cell-type summaries were then computed from per-cell scores and used as the analysis unit for Loop/TAL module-score support. Cell-level UMAP score maps were treated as visualization rather than independent-cell statistical evidence.",
    "",
    "### Functional enrichment analysis",
    "",
    "GO Biological Process enrichment was performed with clusterProfiler and org.Hs.eg.db using MAGMA-tested genes mapped to Entrez IDs as the background universe. Displayed GO terms were filtered by FDR, gene count and redundancy. Curated nephron and functional marker-set enrichment was analyzed by hypergeometric testing against the MAGMA-tested gene universe. These analyses were used for functional interpretation rather than pathway-level validation.",
    "",
    "### GSE73680 injury/remodeling coupling analysis",
    "",
    "GSE73680 injury/remodeling, inflammation/immune, fibrosis/ECM and epithelial-injury programs were scored as mean z-scored marker-set modules. Coupling between MAGMA modules and injury programs was evaluated using sample-level Spearman correlations, paired patient-level delta correlations and patient/group residual correlations. Paired-delta and residual correlations were prioritized for main-text interpretation because sample-level correlations are more robust to disease-group composition than sample-level correlations.",
    "",
    "### Resource-limited TWAS, SMR/coloc and spatial extensions",
    "",
    "TWAS, SMR/coloc and spatial transcriptomic pipelines were scaffolded and resource-audited but were not used as evidence layers because required external expression weights, matched eQTL resources or spatial matrix/image/coordinate files were incomplete."
  )
  results <- c(
    "## Results",
    "",
    "### MAGMA prioritization defines KSD genetic modules",
    "",
    sprintf("The primary public KSD GWAS reconstruction retained %s independent loci and was carried forward into MAGMA gene-level prioritization. MAGMA v1.10 tested 17,316 genes using a GRCh37/hg19-compatible reference and identified 94 Bonferroni-significant genes, 369 FDR-significant genes and 187 suggestive genes at P < 1e-4. These outputs defined the genetic prioritization layer for downstream single-nucleus and disease-context analyses.", loci_n),
    "",
    "### MAGMA-prioritized genes localize to a Loop/TAL-associated snRNA context",
    "",
    "We next projected MAGMA-prioritized KSD gene sets onto the audited GSE231569 single-nucleus atlas (Figure 2). The analysis used harmonized audited cell types and retained Loop/TAL as a high-confidence epithelial transport compartment. Across MAGMA top 50, top 100, FDR-significant and suggestive gene sets, expression-context scores were highest in Loop/TAL cells compared with other audited cell types.",
    "",
    "Random gene-set benchmarking supported the specificity of this localization pattern. Loop/TAL benchmark percentiles were 0.998 for MAGMA top 50, 1.000 for MAGMA top 100, 1.000 for MAGMA suggestive genes and 0.968 for MAGMA FDR genes. Locus-balanced benchmarking, conservative annotation sensitivity and leave-one-locus-out analyses supported robustness of the MAGMA top 50 Loop/TAL signal, with the minimum retained leave-one-locus-out percentile remaining above 0.99. These results support a Loop/TAL-associated renal papillary single-nucleus cellular context for MAGMA-prioritized KSD genes, without establishing causal mediation or spatial validation.",
    "",
    "### P1 genes form a TAL, transport and calcium-handling role spectrum",
    "",
    "The six P1 candidate genes did not behave as a uniform TAL-specific marker set. UMOD served as a representative TAL-associated gene, CLDN10 supported an epithelial transport pattern, CLDN14 and CASR linked the signal to ion-handling and calcium-sensing biology, HIBADH remained a supporting MAGMA-associated context gene and PKD2 represented a broader renal epithelial context. Figure 3 summarizes this evidence across expression, specificity, donor-level descriptive consistency and gene-role assignment. The supported claim is a MAGMA plus snRNA-based TAL-associated cellular context, not causal mediation or colocalized genetic validation.",
    "",
    "### GSE73680 supports MAGMA module-level plaque/stone disease-context association",
    "",
    "We next examined whether the MAGMA- and single-nucleus-prioritized signals were reflected in an independent papillary plaque/stone disease-context dataset. After reconstructing GSE73680 from Agilent supplementary feature files, 55 samples from 29 patients were included, including 26 paired patients with control/adjacent and plaque/stone papilla samples. To account for non-independence among repeated samples, we used a patient-aware limma duplicateCorrelation model and paired patient-level sensitivity analyses.",
    "",
    "At the single-gene level, none of the six P1 candidate genes reached FDR-supported differential expression, although PKD2 showed a nominal exploratory paired response. In contrast, MAGMA-prioritized modules showed consistent paired disease-context shifts. Patient-level MAGMA module responses reached q <= 0.05, retained directionality in leave-one-gene-out analyses, remained robust after removing the six P1 genes, and exceeded size-matched random gene-set expectations. Paired direction analyses showed positive shifts in approximately 69% to 73% of paired patients for MAGMA modules. A stricter expression-matched benchmark did not provide primary support under the current implementation and is treated as a conservative boundary check. These findings indicate that GSE73680 supports plaque/stone disease-context expression association at the MAGMA module level rather than uniform single-gene P1 differential expression.",
    "",
    "### Functional and injury-remodeling analyses link prioritized modules to papillary disease programs",
    "",
    "Cell-level GSE231569 score projection supported the Figure 2 localization pattern on the audited UMAP, with MAGMA top-ranked module scores visibly enriched in the Loop/TAL-associated region. Functional enrichment added a pathway-level interpretation layer. GO Biological Process analysis of MAGMA top 100 genes identified loop of Henle development, distal tubule development, response to vitamin D, response to metal ion and phosphate ion homeostasis among the leading terms. Curated nephron and functional marker-set analyses further linked prioritized and Loop/TAL-influential genes to TAL transport, calcium-ion handling and epithelial junction contexts.",
    "",
    "In GSE73680, the main injury-coupling interpretation used paired patient delta and patient/group residual correlations, with Spearman rho values FDR-corrected across module-program pairs. Sample-level correlations were retained as supplementary context because they are more vulnerable to disease-group composition. MAGMA modules showed positive module-level coupling with papillary injury/remodeling programs in both prioritized analyses. For example, MAGMA top 50 showed paired-delta coupling with injury/remodeling scores (Spearman rho = 0.82, FDR = 4.5e-6) and residual coupling after accounting for patient and group structure (rho = 0.79, FDR = 7.2e-12); MAGMA FDR genes showed the same direction in paired-delta analysis (rho = 0.63, FDR = 0.0015). These analyses support functional interpretation and disease-context coupling, but do not establish a causal injury mechanism."
  )
  discussion <- c(
    "## Discussion",
    "",
    "The main finding of this study is that MAGMA-prioritized kidney stone disease genes localize to a Loop/TAL-associated renal papillary single-nucleus context and show disease-context support at the module level. This supports a module-level post-GWAS interpretation rather than a single causal gene or uniform disease-gene panel model. The major contribution is not a single nominated causal gene, but an auditable evidence framework that separates genetic prioritization, renal papillary cell-context localization, gene-centric interpretation and disease-context expression association.",
    "",
    "The Figure 2 consolidation is central to this interpretation. MAGMA top-ranked, FDR-significant and suggestive gene sets showed Loop/TAL-associated expression-context scores in the audited GSE231569 atlas, exceeded random gene-set expectations and remained supported in locus-balanced and leave-one-locus-out robustness summaries. Because Loop/TAL cells represented a relatively small but manually audited compartment, the interpretation relies on gene-set-level enrichment and robustness checks rather than cell abundance. These results support a TAL-associated single-nucleus cellular context, but they do not identify TAL as a causal cell type, provide spatial validation or substitute for TWAS or colocalization.",
    "",
    "The P1 gene analysis should be interpreted as gene-centric context rather than disease validation. UMOD, CLDN10, CLDN14, CASR, HIBADH and PKD2 form a biologically interpretable TAL, epithelial transport and calcium-handling spectrum, but they do not behave as a uniform TAL marker set. The value of the P1 panel is therefore interpretive rather than confirmatory: it separates representative TAL genes, epithelial transport candidates, calcium-handling genes and broader epithelial-context genes.",
    "",
    "GSE73680 provided disease-context support at the MAGMA module level, while arguing against uniform P1 single-gene differential expression. The module-level signal remained supported in patient-level paired sensitivity analyses and random gene-set benchmarks, but should not be interpreted as cell-type-specific disease expression because GSE73680 is bulk/microarray-based. This distinction keeps Result 4 aligned with the broader post-GWAS mapping framework rather than turning it into single-gene validation.",
    "",
    "The functional-context analyses further help interpret why the Loop/TAL-associated signal is biologically plausible. GO and curated marker-set enrichment linked prioritized genes to nephron development, TAL transport, calcium/ion handling and epithelial junction contexts. In GSE73680, MAGMA modules showed positive paired-delta and patient/group residual coupling with papillary injury/remodeling programs. These findings connect genetic prioritization to renal epithelial transport and papillary injury contexts, but remain functional annotations and disease-context associations rather than causal mechanism tests.",
    "",
    "Future work should prioritize kidney papilla-specific eQTL resources, lesion-resolved spatial transcriptomics and experimental perturbation of top MAGMA-prioritized modules in TAL-relevant epithelial models. These extensions would test whether the TAL-associated cellular context identified here corresponds to causal regulatory mechanisms or spatially localized papillary disease programs.",
    "",
    "## Limitations",
    "",
    "This study is computational and hypothesis-generating. TWAS, SMR/coloc and spatial transcriptomic validation were not completed because required external expression-weight, eQTL and spatial matrix/image resources were unavailable in the current analysis environment. GSE231569 supports single-nucleus expression-context mapping but does not provide spatial validation or causal cell-type mediation. GSE73680 supports disease-context module association but not causality or cell-type-specific disease response. P1 single-gene disease differential expression was not FDR-supported, and PKD2 should be treated only as a nominal exploratory observation. Functional enrichment and risk-injury coupling analyses are gene-set- and module-level interpretations and may be influenced by gene-set size, background definition, expression detectability and disease-group composition.",
    "",
    "## Data availability",
    "",
    "This study used public GWAS summary statistics and public GEO transcriptomic resources. Processed analysis tables, QC summaries and figure source outputs generated for this manuscript are organized with the analysis code and can be deposited in a public repository at submission. Resource limitations for TWAS, SMR/coloc and spatial analyses are documented in the supplementary resource-status tables.",
    "",
    "## Code availability",
    "",
    "Analysis scripts used to reconstruct datasets, perform QC, generate module-level analyses and produce manuscript figures are organized by workflow stage and can be released with the processed outputs at submission."
  )
  legends <- c(
    "## Figure legends",
    "",
    "### Figure 1. Post-GWAS framework for KSD cellular and disease-context mapping",
    "",
    "The study is organized as four evidence layers: GWAS/MAGMA gene prioritization, audited GSE231569 single-nucleus localization, P1 candidate gene evidence and GSE73680 disease-context module analysis. TWAS, SMR/coloc and spatial transcriptomic analyses were prepared as resource-limited extensions but were not used for manuscript claims.",
    "",
    "### Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated single-nucleus expression context",
    "",
    "MAGMA-prioritized gene sets were evaluated across audited renal papillary single-nucleus cell types from GSE231569. Projection, random gene-set benchmarking and robustness analyses supported the Loop/TAL compartment as the strongest cellular expression context for prioritized kidney stone risk genes. Cell-level UMAP score maps are visualization; donor-cell-type summaries provide the statistical interpretation layer. These analyses do not establish causal mediation, TWAS convergence, colocalization or spatial validation.",
    "",
    "### Figure 3. P1 genes form a TAL, transport and calcium-handling role spectrum",
    "",
    "P1 candidate genes were evaluated across audited renal papillary cell types, TAL specificity, donor-level descriptive consistency and gene-role assignments. These analyses support a TAL/transport/calcium-handling evidence spectrum but do not establish causal mediation, TWAS convergence, colocalization or spatial validation.",
    "",
    "### Figure 4. GSE73680 disease-context analysis supports MAGMA-prioritized modules rather than uniform P1 single-gene differential expression",
    "",
    "Patient-aware GSE73680 analyses showed MAGMA module-level disease-context expression association in plaque/stone papilla samples. P1 single-gene responses were heterogeneous and no P1 gene reached q <= 0.05 after FDR correction. Size-matched random benchmarking supported MAGMA module shifts, whereas expression-matched benchmarking was retained as a conservative sensitivity analysis rather than primary support. The analysis supports MAGMA-prioritized module-level disease-context association and does not establish genetic causality, TWAS convergence, colocalization or spatial validation.",
    "",
    "### Figure 5. Functional interpretation of MAGMA-prioritized TAL-associated KSD genes",
    "",
    "Figure 5 summarizes integrated evidence tiers, GO Biological Process enrichment, curated nephron/functional enrichment and GSE73680 risk-injury module coupling. The figure supports functional interpretation and module-level disease-context coupling of prioritized KSD genes, but does not establish causal mechanism, TWAS convergence, colocalization or spatial validation."
  )
  out <- c(title, "", abstract, "", intro, "", methods, "", results, "", discussion, "", legends)
  write_lines(out, "manuscript/manuscript_draft_v1.0.md")
  write_lines(out, "manuscript/manuscript_clean_for_supervisor_v1.0.md")
}

make_figure2_v06 <- function() {
  scores <- fread("results/tables/gse231569_celllevel_magma_scores.tsv")
  top50 <- scores[module_name == "MAGMA_top50"]
  top50[, cell_label := label_cell(audited_broad_cell_type)]
  top50[, is_tal := audited_broad_cell_type == "Loop_of_Henle_TAL"]
  set.seed(12)
  other <- top50[is_tal == FALSE]
  tal <- top50[is_tal == TRUE]
  other <- if (nrow(other) > 24000) other[sample(.N, 24000)] else other
  plot_cells <- rbind(other, tal, fill = TRUE)
  cell_cols <- c(
    "Collecting duct" = "#D9DEE1",
    "Endothelial" = "#BAC7B1",
    "Fibroblast/stromal" = "#B5CDBA",
    "Injured epithelial" = "#ADDADB",
    "Loop/TAL" = pal$genetic,
    "Perivascular/mural-like" = "#D7C5D7"
  )
  center <- tal[, .(x = median(UMAP_1), y = median(UMAP_2))]
  p_a <- ggplot(plot_cells[is_tal == FALSE], aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_label), size = 0.16, alpha = 0.18) +
    geom_point(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2), color = pal$genetic, size = 0.42, alpha = 0.95) +
    stat_ellipse(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                 color = pal$p1, linewidth = 0.5, level = 0.85) +
    annotate("curve", x = center$x + 3.0, y = center$y + 3.5, xend = center$x + 0.45, yend = center$y + 0.45,
             curvature = -0.25, arrow = arrow(length = unit(0.11, "inches")), color = pal$p1, linewidth = 0.45) +
    annotate("label", x = center$x + 3.2, y = center$y + 3.7, label = "Loop/TAL\nn = 540", hjust = 0,
             size = 2.6, fill = "white", color = pal$genetic, fontface = "bold") +
    scale_color_manual(values = cell_cols, name = "Audited cell type") +
    labs(title = "A. Audited GSE231569 snRNA-seq atlas", x = "UMAP 1", y = "UMAP 2") +
    theme_pub(8.8) +
    theme(legend.position = "right", panel.grid = element_blank())

  qlim <- quantile(plot_cells$celllevel_module_score, c(0.02, 0.99), na.rm = TRUE)
  p_b <- ggplot(plot_cells, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celllevel_module_score), size = 0.18, alpha = 0.78) +
    stat_density_2d(data = plot_cells[is_tal == TRUE], aes(UMAP_1, UMAP_2), inherit.aes = FALSE,
                    color = pal$p1, linewidth = 0.48, bins = 4) +
    scale_color_gradient(low = "#F1F3F4", high = pal$genetic, limits = qlim, oob = squish,
                         name = "Relative\nmodule score") +
    labs(title = "B. MAGMA top 50 module score", x = "UMAP 1", y = "UMAP 2") +
    theme_pub(8.8) +
    theme(panel.grid = element_blank())

  random <- fread("results/tables/magma_scrna_random_benchmark.tsv")
  random <- random[gene_set %in% c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4")]
  random[, module_label := factor(gene_set, levels = c("magma_top50", "magma_top100", "magma_fdr05", "magma_suggestive_p1e4"),
                                  labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  y_order <- c("Loop/TAL", "Perivascular/mural-like", "Fibroblast/stromal", "Endothelial", "Injured epithelial", "Collecting duct")
  random[, cell_label := factor(label_cell(audited_broad_cell_type), levels = rev(y_order))]
  random[, support := fifelse(audited_broad_cell_type == "Loop_of_Henle_TAL" & benchmark_percentile >= 0.95, "Loop/TAL exceeded expectation", "Other")]
  p_c <- ggplot(random, aes(benchmark_percentile, cell_label, fill = support)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "#555555", linewidth = 0.35) +
    annotate("text", x = 0.94, y = Inf, label = "95th percentile", hjust = 1, vjust = 1.2,
             size = 2.8, fontface = "bold", color = "#555555") +
    geom_point(shape = 21, size = 2.8, color = "#555555", stroke = 0.22) +
    facet_wrap(~ module_label, ncol = 2) +
    scale_fill_manual(values = c("Loop/TAL exceeded expectation" = pal$genetic, Other = pal$muted)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1.0)) +
    labs(title = "C. Size-matched benchmark percentile", x = "Benchmark percentile", y = NULL, fill = NULL) +
    theme_pub(8.8) +
    theme(legend.position = "bottom")

  infl <- fread("results/tables/loop_tal_influential_magma_genes.tsv")[order(contribution_rank)][1:12]
  infl[, gene := factor(gene, levels = rev(gene))]
  infl[, group := fifelse(candidate_role == "P1_core", "P1 candidate", "Other MAGMA gene")]
  p_d <- ggplot(infl, aes(contribution_score, gene)) +
    geom_segment(aes(x = 0, xend = contribution_score, yend = gene), color = pal$grid, linewidth = 0.42) +
    geom_point(aes(fill = group, size = donor_detection), shape = 21, color = "#444444", stroke = 0.22) +
    scale_fill_manual(values = c("P1 candidate" = pal$p1, "Other MAGMA gene" = pal$genetic)) +
    scale_size_continuous(range = c(2.2, 4.6), labels = percent_format()) +
    labs(title = "D. Genes contributing to Loop/TAL-associated signal", x = "Contribution score", y = NULL,
         fill = NULL, size = "Detection") +
    theme_pub(9) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(1.05, 1.05))
  ggsave("results/figures/figure2_magma_scrna_localization_v0.6.pdf", fig, width = 13.2, height = 9.6, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure2_magma_scrna_localization_v0.6.png", fig, width = 13.2, height = 9.6, units = "in", dpi = 320, bg = "white")

  donor <- fread("results/tables/gse231569_donor_celltype_magma_scores_v0.2.tsv")
  donor2 <- donor[module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive")]
  donor2[, cell_label := label_cell(audited_broad_cell_type)]
  contrasts <- rbindlist(lapply(unique(donor2$module_name), function(m) {
    d <- donor2[module_name == m]
    rbindlist(list(
      data.table(module_name = m, contrast = "Loop/TAL vs all other audited contexts",
                 n_loop_tal = d[audited_broad_cell_type == "Loop_of_Henle_TAL", .N],
                 n_comparator = d[audited_broad_cell_type != "Loop_of_Henle_TAL", .N],
                 loop_tal_median = d[audited_broad_cell_type == "Loop_of_Henle_TAL", median(mean_score, na.rm = TRUE)],
                 comparator_median = d[audited_broad_cell_type != "Loop_of_Henle_TAL", median(mean_score, na.rm = TRUE)],
                 wilcox_p = suppressWarnings(wilcox.test(mean_score ~ (audited_broad_cell_type == "Loop_of_Henle_TAL"), data = d)$p.value),
                 interpretation = "descriptive donor-cell-type support; not independent-cell inference"),
      data.table(module_name = m, contrast = "Loop/TAL vs collecting duct",
                 n_loop_tal = d[audited_broad_cell_type == "Loop_of_Henle_TAL", .N],
                 n_comparator = d[audited_broad_cell_type == "Collecting_duct_principal", .N],
                 loop_tal_median = d[audited_broad_cell_type == "Loop_of_Henle_TAL", median(mean_score, na.rm = TRUE)],
                 comparator_median = d[audited_broad_cell_type == "Collecting_duct_principal", median(mean_score, na.rm = TRUE)],
                 wilcox_p = suppressWarnings(wilcox.test(d[audited_broad_cell_type == "Loop_of_Henle_TAL", mean_score],
                                                         d[audited_broad_cell_type == "Collecting_duct_principal", mean_score])$p.value),
                 interpretation = "descriptive donor-cell-type support; not independent-cell inference"),
      data.table(module_name = m, contrast = "Loop/TAL vs injured epithelial",
                 n_loop_tal = d[audited_broad_cell_type == "Loop_of_Henle_TAL", .N],
                 n_comparator = d[audited_broad_cell_type == "Injured_undifferentiated_epithelial", .N],
                 loop_tal_median = d[audited_broad_cell_type == "Loop_of_Henle_TAL", median(mean_score, na.rm = TRUE)],
                 comparator_median = d[audited_broad_cell_type == "Injured_undifferentiated_epithelial", median(mean_score, na.rm = TRUE)],
                 wilcox_p = suppressWarnings(wilcox.test(d[audited_broad_cell_type == "Loop_of_Henle_TAL", mean_score],
                                                         d[audited_broad_cell_type == "Injured_undifferentiated_epithelial", mean_score])$p.value),
                 interpretation = "descriptive donor-cell-type support; not independent-cell inference")
    ), fill = TRUE)
  }), fill = TRUE)
  contrasts[, fdr := p.adjust(wilcox_p, "BH")]
  fwrite(contrasts, "results/tables/gse231569_donor_celltype_magma_score_tests.tsv", sep = "\t")

  donor_plot <- donor2[audited_broad_cell_type %in% c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")]
  donor_plot[, cell_label := factor(label_cell(audited_broad_cell_type), levels = c("Loop/TAL", "Collecting duct", "Injured epithelial", "Endothelial", "Fibroblast/stromal", "Perivascular/mural-like"))]
  donor_plot[, module_label := factor(module_name, levels = c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive"),
                                      labels = c("Top 50", "Top 100", "FDR", "Suggestive"))]
  p_supp <- ggplot(donor_plot, aes(cell_label, mean_score, fill = cell_label)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, color = "#555555", linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.4, alpha = 0.8, color = "#303030") +
    facet_wrap(~ module_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("Loop/TAL" = pal$genetic, "Collecting duct" = "#D9DEE1", "Injured epithelial" = "#ADDADB",
                                 "Endothelial" = "#BAC7B1", "Fibroblast/stromal" = "#B5CDBA", "Perivascular/mural-like" = "#D7C5D7")) +
    labs(title = "Donor-cell-type MAGMA module score support", x = NULL, y = "Mean module score") +
    theme_pub(8.6) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  ggsave("results/figures/supp_gse231569_donor_celltype_scores.pdf", p_supp, width = 9.0, height = 6.2, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/supp_gse231569_donor_celltype_scores.png", p_supp, width = 9.0, height = 6.2, units = "in", dpi = 320, bg = "white")

  write_lines(c(
    "# Figure 2 Legend v0.6",
    "",
    "**Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated renal papillary single-nucleus context.**",
    "(A) Audited GSE231569 snRNA-seq UMAP with non-Loop/TAL cells shown in low-saturation colors and the Loop/TAL compartment highlighted by deep blue-green points, an outline and an arrow. (B) MAGMA top 50 module score projected onto the same UMAP. Relative module score represents mean z-scored expression across detected MAGMA top 50 genes; the gold contour marks the audited Loop/TAL compartment. (C) Size-matched random gene-set benchmark percentiles, with the dashed line marking the 95th percentile. (D) Genes contributing to the Loop/TAL-associated signal. Contribution scores summarize MAGMA prioritization strength, Loop/TAL expression preference and detection support for descriptive ranking, not causal driver inference. Cell-level score maps are visualization; donor-cell-type summaries provide the statistical interpretation layer. These analyses do not establish causal cell-type mediation, TWAS convergence, colocalization or spatial validation."
  ), "docs/figure2_legend_v0.6.md")
}

make_figure3_v08 <- function() {
  p1 <- fread("results/tables/p1_tal_gene_interpretation_summary.tsv")
  gene_order <- c("UMOD", "CLDN10", "CLDN14", "CASR", "HIBADH", "PKD2")
  p1[, gene := factor(gene, levels = gene_order)]
  cards <- data.table(
    gene = factor(gene_order, levels = gene_order),
    role = c("Representative TAL", "Transport", "Ion handling", "Calcium sensing", "Supporting", "Broad epithelial"),
    x = c(1, 2, 3, 4, 2, 3),
    y = c(2, 2, 2, 2, 1, 1)
  )
  role_cols <- c("Representative TAL" = pal$genetic, "Transport" = "#557F89",
                 "Ion handling" = pal$p1, "Calcium sensing" = "#9A7E43",
                 "Supporting" = pal$scrna, "Broad epithelial" = pal$disease)
  p_a <- ggplot(cards, aes(x, y)) +
    annotate("segment", x = 0.8, xend = 4.2, y = 2.55, yend = 2.55, color = pal$grid, linewidth = 0.4,
             arrow = arrow(length = unit(0.10, "inches"))) +
    annotate("text", x = 0.78, y = 2.72, label = "TAL identity", hjust = 0, size = 2.5, color = "#555555") +
    annotate("text", x = 4.2, y = 2.72, label = "broader epithelial context", hjust = 1, size = 2.5, color = "#555555") +
    geom_label(aes(label = paste0(as.character(gene), "\n", role), fill = role),
               color = "white", fontface = "bold", size = 2.75, lineheight = 0.9,
               label.padding = unit(0.22, "lines"), linewidth = 0) +
    scale_fill_manual(values = role_cols) +
    coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.55, 2.95), clip = "off") +
    labs(title = "A. P1 gene role cards") +
    theme_void(base_size = 9.2) +
    theme(plot.title = element_text(face = "bold", color = pal$text), legend.position = "none")

  cell <- fread("results/tables/p1_tal_gene_celltype_summary.tsv")
  keep_cells <- c("Loop_of_Henle_TAL", "Collecting_duct_principal", "Injured_undifferentiated_epithelial", "Endothelial", "Fibroblast_stromal", "Perivascular_mural_like")
  cell <- cell[cell_type %in% keep_cells]
  cell[, cell_label := factor(label_cell(cell_type, short = TRUE), levels = label_cell(keep_cells, short = TRUE))]
  cell[, gene := factor(gene, levels = gene_order)]
  p_b <- ggplot(cell, aes(cell_label, gene)) +
    geom_point(aes(size = pct_expressed, fill = avg_expression), shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    scale_size_continuous(range = c(0.8, 5.0), labels = percent_format()) +
    labs(title = "B. P1 expression across audited cell types", x = NULL, y = NULL, fill = "Average\nexpression", size = "Detected") +
    theme_pub(8.8) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_text(face = "italic"),
          legend.position = "bottom", panel.grid = element_blank())

  p_c <- ggplot(p1, aes(specificity_ratio_avg, gene)) +
    geom_segment(aes(x = 0, xend = specificity_ratio_avg, yend = gene), color = pal$grid, linewidth = 0.5) +
    geom_point(aes(fill = specificity_class), shape = 21, size = 3.8, color = "#555555", stroke = 0.25) +
    annotate("text", x = Inf, y = Inf, label = "All shown genes detected\nin 3/4 TAL donors", hjust = 1.05, vjust = 1.25,
             size = 2.55, color = "#555555") +
    scale_fill_manual(values = c(strong_TAL_preferential = pal$genetic, moderate_TAL_preferential = pal$p1), na.value = pal$muted) +
    labs(title = "C. TAL specificity ratio", x = "TAL specificity ratio", y = NULL, fill = "Specificity") +
    theme_pub(8.8) +
    theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

  ev <- p1[, .(
    gene,
    MAGMA = fifelse(magma_p < 5e-8, "+++", "++"),
    TAL_specificity = fifelse(specificity_class == "strong_TAL_preferential", "+++", "++"),
    Donor_detection = sprintf("%d/4", round(TAL_donor_detection_fraction * 4)),
    Bulk_response = fifelse(as.character(gene) == "PKD2", "nominal", "no FDR"),
    Role = fifelse(as.character(gene) == "UMOD", "Representative TAL",
                   fifelse(as.character(gene) == "CLDN10", "Transport",
                           fifelse(as.character(gene) == "CLDN14", "Ion handling",
                                   fifelse(as.character(gene) == "CASR", "Calcium sensing",
                                           fifelse(as.character(gene) == "HIBADH", "Supporting", "Broad epithelial")))))
  )]
  ev_long <- melt(ev, id.vars = "gene", variable.name = "evidence", value.name = "call")
  ev_long[, gene := factor(gene, levels = rev(gene_order))]
  ev_long[, evidence := factor(evidence, levels = c("MAGMA", "TAL_specificity", "Donor_detection", "Bulk_response", "Role"),
                               labels = c("MAGMA", "TAL specificity", "Donor detection", "Bulk response", "Role"))]
  ev_long[, fill_class := fcase(
    call == "+++", "strong",
    call == "++", "moderate",
    call == "nominal", "nominal",
    call == "no FDR", "no_fdr",
    default = "role"
  )]
  p_d <- ggplot(ev_long, aes(evidence, gene, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 2.25, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, nominal = pal$p1,
                                 no_fdr = pal$muted, role = "#F7F8F8")) +
    labs(title = "D. Discrete evidence matrix", x = NULL, y = NULL) +
    theme_pub(8.8) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.text.y = element_text(face = "italic"),
          legend.position = "none", panel.grid = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.86, 1.16))
  ggsave("results/figures/figure3_p1_gene_evidence_v0.8.pdf", fig, width = 13.2, height = 8.9, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure3_p1_gene_evidence_v0.8.png", fig, width = 13.2, height = 8.9, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 3 Legend v0.8",
    "",
    "**Figure 3. P1 genes form a TAL, transport and calcium-handling role spectrum.**",
    "(A) Gene role cards summarize the six P1 candidates from representative TAL biology to broader epithelial context. (B) P1 gene expression and detection across audited GSE231569 cell types. Short labels indicate Loop/TAL, collecting duct (CD), injured epithelial (Injured epi.), endothelial, fibroblast/stromal (Fib/stromal) and perivascular/mural-like (Perivasc.) compartments. (C) TAL specificity ratio, with donor support summarized as all shown genes being detected in 3/4 TAL donors. (D) Discrete evidence matrix summarizing MAGMA support, TAL specificity, donor detection, bulk single-gene response and assigned manuscript role. The Bulk response column refers to GSE73680 single-gene disease-context response. Discrete calls are used to avoid implying that heterogeneous evidence dimensions are directly additive."
  ), "docs/figure3_legend_v0.8.md")
}

make_figure5_v03 <- function() {
  checklist <- data.table(
    evidence = rep(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context"), each = 4),
    gene_set = rep(c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"), times = 5),
    call = c("+++", "+++", "++", "+",
             "++", "++", "+++", "++",
             "+", "+", "+", "+",
             "++", "++", "NA", "+",
             "++", "++", "++", "++")
  )
  checklist[, evidence := factor(evidence, levels = rev(c("MAGMA rank", "Loop/TAL specificity", "Donor detection", "GSE73680 coupling", "Functional context")))]
  checklist[, gene_set := factor(gene_set, levels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  checklist[, fill_class := fcase(call == "+++", "strong", call == "++", "moderate", call == "+", "support", call == "NA", "not_applicable")]
  p_a <- ggplot(checklist, aes(gene_set, evidence, fill = fill_class)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = call), size = 3.1, color = "#303030") +
    scale_fill_manual(values = c(strong = pal$genetic, moderate = pal$scrna, support = pal$p1, not_applicable = pal$muted)) +
    labs(title = "A. Evidence checklist", x = NULL, y = NULL) +
    theme_pub(8.8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none", panel.grid = element_blank())

  go <- fread("results/tables/go_bp_redundancy_reduced_terms.tsv")
  go <- go[redundancy_reduced_keep == TRUE & p.adjust < 0.10 & Count >= 2]
  go <- go[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
  selected_terms <- c(
    "distal tubule development",
    "loop of Henle development",
    "intracellular calcium ion homeostasis",
    "regulation of calcium ion import",
    "regulation of monoatomic ion transport",
    "urate metabolic process",
    "phosphate ion homeostasis",
    "cell-cell junction assembly",
    "cell-cell adhesion via plasma-membrane adhesion molecules"
  )
  go <- go[Description %in% selected_terms]
  go[, theme := fcase(
    Description %in% c("distal tubule development", "loop of Henle development"), "Nephron development",
    Description %in% c("intracellular calcium ion homeostasis", "regulation of calcium ion import",
                       "regulation of monoatomic ion transport", "urate metabolic process",
                       "phosphate ion homeostasis"), "Ion/mineral handling",
    Description %in% c("cell-cell junction assembly", "cell-cell adhesion via plasma-membrane adhesion molecules"), "Epithelial context"
  )]
  go[, theme_rank := match(theme, c("Nephron development", "Ion/mineral handling", "Epithelial context"))]
  go[, term_rank := match(Description, selected_terms)]
  setorder(go, theme_rank, term_rank, p.adjust)
  display_levels <- unique(go[order(theme_rank, term_rank), paste0(theme, " | ", Description)])
  go[, Description_short := paste0(theme, " | ", Description)]
  go[, Description_short := factor(Description_short, levels = rev(display_levels))]
  go[, gene_set_label := factor(gene_set,
                                levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  p_b <- ggplot(go, aes(-log10(p.adjust), Description_short, size = Count, fill = gene_set_label)) +
    geom_point(shape = 21, color = "#555555", stroke = 0.22) +
    scale_fill_manual(values = c("MAGMA top 100" = pal$genetic, "MAGMA FDR" = pal$scrna,
                                 "Loop/TAL contributors" = pal$p1, "P1 core" = pal$disease)) +
    labs(title = "B. Redundancy-reduced GO BP terms", subtitle = "Functional interpretation; not pathway validation",
         x = "-log10(FDR)", y = NULL, fill = "Gene set", size = "Genes") +
    theme_pub(8.2) +
    theme(legend.position = "bottom")

  curated <- fread("results/tables/nephron_segment_marker_enrichment.tsv")
  curated <- curated[gene_set %in% c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core")]
  curated <- curated[term != "papillary_injury_remodeling"]
  curated[, gene_set_label := factor(gene_set, levels = c("MAGMA_top100", "MAGMA_FDR", "Loop_TAL_influential", "P1_core"),
                                     labels = c("MAGMA top 100", "MAGMA FDR", "Loop/TAL contributors", "P1 core"))]
  curated[, term_label := factor(term,
                                 levels = rev(c("TAL_transport", "calcium_ion_handling", "epithelial_tight_junction",
                                                "proximal_tubule_context", "collecting_duct_context")),
                                 labels = rev(c("TAL transport", "Calcium ion handling", "Epithelial tight junction",
                                                "Proximal tubule context", "Collecting duct context")))]
  p_c <- ggplot(curated, aes(gene_set_label, term_label, fill = pmin(enrichment_ratio, 20))) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(overlap > 0, overlap, "")), size = 2.65, color = "#303030") +
    scale_fill_gradient(low = "#EEF3F4", high = pal$genetic) +
    labs(title = "C. Curated nephron and functional context", subtitle = "Color = enrichment ratio; number = overlapping genes",
         x = NULL, y = NULL, fill = "Enrichment\nratio") +
    theme_pub(8.5) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  robust <- fread("results/tables/gse73680_risk_injury_correlation_robustness.tsv")
  d <- robust[analysis %in% c("Paired delta", "Patient/group residual") &
                module_name %in% c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates") &
                injury_module %in% c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune")]
  d[, module_label := factor(module_name, levels = rev(c("MAGMA_top50", "MAGMA_top100", "MAGMA_FDR", "MAGMA_suggestive", "P1_core_TAL_candidates")),
                             labels = rev(c("MAGMA top 50", "MAGMA top 100", "MAGMA FDR", "MAGMA suggestive", "P1 core")))]
  d[, injury_label := factor(injury_module,
                             levels = c("injury_remodeling", "epithelial_injury", "fibrosis_ecm", "inflammation_immune"),
                             labels = c("Injury/remodeling", "Epithelial injury", "Fibrosis/ECM", "Inflammation/immune"))]
  d[, analysis := factor(analysis, levels = c("Paired delta", "Patient/group residual"),
                         labels = c("Paired patient delta", "Patient/group residual"))]
  d[, sig := fifelse(fdr <= 0.001, "***", fifelse(fdr <= 0.01, "**", fifelse(fdr <= 0.05, "*", "")))]
  d[, label := sprintf("%.2f%s", rho, sig)]
  p_d <- ggplot(d, aes(injury_label, module_label, fill = pmin(pmax(rho, 0), 1))) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.35, color = "#303030") +
    facet_wrap(~ analysis, ncol = 1) +
    scale_fill_gradient(low = "#EEF3F4", high = pal$disease, limits = c(0, 1), name = "Spearman\nrho") +
    labs(title = "D. Risk-injury coupling robustness", x = NULL, y = NULL) +
    theme_pub(8.4) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 2, rel_heights = c(0.9, 1.1))
  ggsave("results/figures/figure5_functional_context_v0.3.pdf", fig, width = 13.2, height = 9.5, units = "in", device = "pdf", bg = "white")
  ggsave("results/figures/figure5_functional_context_v0.3.png", fig, width = 13.2, height = 9.5, units = "in", dpi = 320, bg = "white")
  write_lines(c(
    "# Figure 5 Legend v0.3",
    "",
    "**Figure 5. Functional context and risk-injury coupling of MAGMA-prioritized TAL-associated KSD genes.**",
    "(A) Evidence checklist across MAGMA-prioritized, Loop/TAL-contributing and P1 gene sets. Symbols denote +++ strong support, ++ moderate support, + detectable or contextual support and NA not applicable. (B) Redundancy-reduced GO Biological Process terms grouped by nephron development, ion/mineral handling and epithelial context. GO enrichment is used for functional interpretation rather than pathway-level validation. (C) Curated nephron and functional marker-set enrichment; color indicates enrichment ratio and tile labels indicate overlapping genes. The papillary injury/remodeling marker row is not shown in this curated panel because injury/remodeling coupling is evaluated separately in panel D. (D) Risk-injury coupling robustness in GSE73680 using paired patient delta and patient/group residual correlations. Values represent Spearman rho; FDR correction was applied across tested module-program pairs, with *, ** and *** marking FDR <= 0.05, 0.01 and 0.001. These analyses support module-level disease-context coupling, not causal injury mechanisms."
  ), "docs/figure5_legend_v0.3.md")
}

make_supplement_and_package <- function() {
  supp_tables <- data.table(
    table_id = paste0("Table S", 1:13),
    title = c(
      "GWAS QC and lead loci",
      "MAGMA gene-based results",
      "MAGMA gene sets and benchmark results",
      "GSE231569 annotation audit and cell counts",
      "Donor-cell-type MAGMA module score statistics",
      "P1 gene evidence summary",
      "GSE73680 sample metadata and patient pairing",
      "GSE73680 P1 single-gene response",
      "GSE73680 MAGMA module response and robustness",
      "Risk-injury coupling robustness",
      "GO and curated functional enrichment",
      "TWAS/SMR/coloc/spatial resource status",
      "Integrated candidate gene tiers"
    ),
    source_file = c(
      "results/phase1_v0.1/tables/phase1_gwas_qc_report.tsv; results/phase1_v0.1/tables/phase1_2025_lead_snps.tsv; results/phase1_v0.1/tables/phase1_2025_loci.tsv",
      "results/tables/magma_genes.tsv",
      "results/tables/magma_scrna_random_benchmark.tsv; results/tables/loop_tal_influential_magma_genes.tsv",
      "results/tables/gse231569_celltype_counts.tsv or annotation audit source tables",
      "results/tables/gse231569_donor_celltype_magma_score_tests.tsv; results/tables/gse231569_donor_celltype_magma_score_stats.tsv",
      "results/tables/p1_tal_gene_interpretation_summary.tsv; results/tables/p1_tal_gene_celltype_summary.tsv",
      "results/gse73680/tables/gse73680_sample_metadata.tsv; results/gse73680/tables/gse73680_patient_pairing.tsv",
      "results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv",
      "results/gse73680/tables/gse73680_patient_level_module_response.tsv; results/tables/gse73680_risk_injury_leave_one_patient_out.tsv",
      "results/tables/gse73680_risk_injury_correlation_robustness.tsv",
      "results/tables/go_bp_redundancy_reduced_terms.tsv; results/tables/nephron_segment_marker_enrichment.tsv",
      "results/twas/twas_resource_status.tsv; results/smr_coloc/eqtl_resource_status.tsv; results/spatial/spatial_resource_status.tsv",
      "results/tables/candidate_gene_tiers_v1.0.tsv"
    ),
    status = c("ready", "ready", "ready", "source_audit_needed_for_final_table", "ready", "ready",
               "source_audit_needed_for_final_table", "ready", "ready", "ready", "ready", "ready", "ready")
  )
  fwrite(supp_tables, "docs/supplementary_table_plan_v1.0.tsv", sep = "\t")

  supp_figs <- data.table(
    figure_id = paste0("Figure S", 1:10),
    title = c(
      "GWAS QQ and Manhattan plots",
      "GSE231569 annotation marker audit",
      "MAGMA module donor-cell-type boxplots",
      "Leave-one-locus-out robustness",
      "P1 individual gene expression/detection plots",
      "GSE73680 QC and normalization",
      "Expression-matched conservative benchmark",
      "Risk-injury coupling robustness",
      "GO enrichment full/redundancy-reduced terms",
      "TWAS/SMR/spatial resource audit workflow"
    ),
    source_file = c(
      "results/figures/phase1_gwas_2025_manhattan_plot.png; results/figures/phase1_gwas_2025_qq_plot.pdf",
      "annotation marker audit figures/source tables",
      "results/figures/supp_gse231569_donor_celltype_scores.pdf",
      "leave-one-locus-out source figure/table",
      "P1 individual expression source plots",
      "GSE73680 reconstruction QC source plots",
      "expression-matched benchmark source figure",
      "results/figures/figure5_risk_injury_correlation_robustness.pdf; results/figures/figure5_functional_context_v0.3.pdf",
      "GO full/redundancy-reduced source figure/table",
      "resource status workflow schematic planned from resource_status tables"
    ),
    status = c("existing", "planned_or_existing_source", "ready", "planned_or_existing_source", "planned_or_existing_source",
               "planned_or_existing_source", "planned_or_existing_source", "existing", "planned_or_existing_source", "planned")
  )
  fwrite(supp_figs, "docs/supplementary_figure_plan_v1.0.tsv", sep = "\t")

  fig_qc <- data.table(
    figure_id = paste0("Figure ", 1:5),
    file = c(
      "results/figures/figure1_integrative_framework_v0.4.pdf",
      "results/figures/figure2_magma_scrna_localization_v0.6.pdf",
      "results/figures/figure3_p1_gene_evidence_v0.8.pdf",
      "results/figures/figure4_gse73680_disease_context_v1.0.pdf",
      "results/figures/figure5_functional_context_v0.3.pdf"
    ),
    claim_supported = c(
      "Claim-bounded four-layer post-GWAS evidence framework",
      "MAGMA-prioritized genes localize to a Loop/TAL-associated snRNA context",
      "P1 genes form a TAL, transport and calcium-handling role spectrum",
      "GSE73680 supports MAGMA module-level plaque/stone disease-context association",
      "Prioritized modules link to nephron/transport functional context and injury-remodeling coupling"
    ),
    claim_not_supported = c(
      "mechanism model or completed TWAS/SMR/spatial evidence",
      "causal cell type, TWAS convergence, colocalization or spatial validation",
      "P1 disease-gene validation or additive evidence score",
      "P1 single-gene disease validation or cell-type-specific bulk response",
      "causal injury mechanism or pathway validation"
    ),
    statistical_unit = c(
      "Evidence layer schematic",
      "Cell-level visualization; donor-cell-type summaries for interpretation",
      "Gene-level descriptive evidence and donor detection",
      "Patient-aware bulk samples and paired patients",
      "Gene sets/modules; paired patient delta and patient/group residual correlations"
    ),
    primary_statistic = c(
      "Fixed analysis counts and resource boundary",
      "Size-matched benchmark percentile; donor-cell-type score tests",
      "TAL specificity ratio; expression/detection summaries",
      "Patient-aware module response; size-matched random benchmark",
      "GO FDR/enrichment ratio; Spearman rho with FDR"
    ),
    boundary = c(
      "Supported and not-established claims explicitly separated",
      "Context mapping only",
      "Interpretive candidate spectrum only",
      "Module-level disease context only",
      "Functional interpretation and non-causal coupling only"
    ),
    visual_status = c("carried_forward_v0.4", "updated_v0.6", "updated_v0.8", "carried_forward_v1.0", "updated_v0.3"),
    legend_status = c("ready_v0.4", "ready_v0.6", "ready_v0.8", "ready_v1.0", "ready_v0.3")
  )
  fwrite(fig_qc, "results/tables/main_figure_qc_v1.0.tsv", sep = "\t")

  boundary_md <- c(
    "# Claim Boundary Audit v1.0",
    "",
    "## Supported Claims",
    "",
    "1. MAGMA-prioritized KSD genes support a Loop/TAL-associated renal papillary snRNA context.",
    "2. P1 genes provide an interpretable TAL, transport and calcium-handling role spectrum.",
    "3. GSE73680 supports MAGMA module-level plaque/stone disease-context association.",
    "4. Functional enrichment and risk-injury coupling support module-level functional interpretation and disease-context coupling.",
    "",
    "## Not Established",
    "",
    "1. Causal mediation by Loop/TAL cells.",
    "2. A single causal KSD gene or uniform P1 disease-gene validation.",
    "3. TWAS convergence, SMR/colocalization or spatial transcriptomic validation.",
    "4. Cell-type-specific disease response in GSE73680 bulk microarray data.",
    "5. Causal injury/remodeling mechanism or therapeutic target validation.",
    "",
    "## Required Wording",
    "",
    "Use: supports, maps to, localizes to, disease-context association, module-level coupling, functional interpretation.",
    "",
    "Avoid: proves, validates causality, causal cell type, driver gene, therapeutic target, colocalized, TWAS-supported, spatially validated."
  )
  write_lines(boundary_md, "docs/claim_boundary_audit_v1.0.md")

  package_dir <- "results/supervisor_review_package_v1.0"
  files_to_copy <- c(
    "manuscript/manuscript_clean_for_supervisor_v1.0.md",
    "results/figures/figure1_integrative_framework_v0.4.pdf",
    "results/figures/figure1_integrative_framework_v0.4.png",
    "results/figures/figure2_magma_scrna_localization_v0.6.pdf",
    "results/figures/figure2_magma_scrna_localization_v0.6.png",
    "results/figures/figure3_p1_gene_evidence_v0.8.pdf",
    "results/figures/figure3_p1_gene_evidence_v0.8.png",
    "results/figures/figure4_gse73680_disease_context_v1.0.pdf",
    "results/figures/figure4_gse73680_disease_context_v1.0.png",
    "results/figures/figure5_functional_context_v0.3.pdf",
    "results/figures/figure5_functional_context_v0.3.png",
    "docs/figure1_legend_v0.4.md",
    "docs/figure2_legend_v0.6.md",
    "docs/figure3_legend_v0.8.md",
    "docs/figure4_legend_v1.0.md",
    "docs/figure5_legend_v0.3.md",
    "docs/supplementary_table_plan_v1.0.tsv",
    "docs/supplementary_figure_plan_v1.0.tsv",
    "docs/claim_boundary_audit_v1.0.md",
    "results/tables/main_figure_qc_v1.0.tsv",
    "results/tables/candidate_gene_tiers_v1.0.tsv",
    "results/twas/twas_resource_status.tsv",
    "results/smr_coloc/eqtl_resource_status.tsv",
    "results/spatial/spatial_resource_status.tsv"
  )
  copied <- data.table(file = files_to_copy, exists = file.exists(files_to_copy), copied = FALSE)
  for (i in seq_len(nrow(copied))) {
    if (copied$exists[i]) {
      copied$copied[i] <- file.copy(copied$file[i], file.path(package_dir, basename(copied$file[i])), overwrite = TRUE)
    }
  }
  fwrite(copied, file.path(package_dir, "package_manifest_v1.0.tsv"), sep = "\t")

  write_lines(c(
    "# Supervisor Review Cover Memo v1.0",
    "",
    "This package contains the v1.0 supervisor-review manuscript, five main candidate figures, figure legends, supplementary material plans, main-figure QC and claim-boundary audit.",
    "",
    "Core claim: MAGMA-prioritized KSD genes support a Loop/TAL-associated renal papillary snRNA context and MAGMA module-level plaque/stone papilla disease-context association.",
    "",
    "Boundary: TWAS, SMR/coloc and spatial transcriptomic extensions remain resource-limited and are not used as manuscript evidence layers. P1 genes are interpreted as a TAL/transport/calcium-handling role spectrum, not as FDR-supported single-gene disease validation."
  ), file.path(package_dir, "supervisor_review_cover_memo_v1.0.md"))
}

make_manuscript_v10()
make_figure2_v06()
make_figure3_v08()
make_figure5_v03()
make_supplement_and_package()

message("Phase 12 v1.0 convergence outputs written")
