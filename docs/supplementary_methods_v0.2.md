# Supplementary Methods v0.2

## Scope and claim boundary

These methods document the completed evidence layers supporting Loop/TAL-associated single-nucleus expression context and module-level plaque/stone-associated papillary bulk expression association. TWAS, SMR/coloc and spatial workflows were audited but incomplete and were not used for claims. No analysis described here establishes causality, P1 single-gene disease validation, spatial validation or experimental mechanism.

## GWAS reconstruction and variant quality control

The 2025 trans-ancestry KSD summary statistics were treated as GRCh37/hg19-compatible. Required fields were variant identifier, chromosome, base-pair position, effect and non-effect alleles, P value, beta, standard error, effect-allele frequency and sample size. The file did not contain an INFO field, so no imputation-quality filter was possible. Rows with missing required values, invalid P values or coordinates, non-SNV/non-ACGT alleles, duplicate variant identifiers, strand-ambiguous A/T or C/G alleles, or minor-allele frequency <0.01 were removed. The sequential counts were 5,960,489 raw rows; 72,210 non-SNV/bad-allele rows, 895,763 ambiguous rows and 77,483 low-frequency rows removed; and 4,915,033 retained variants. The genomic-control lambda was 1.127 and 3,209 retained variants had P <5 x 10^-8.

Genome-wide significant variants were sorted by P value and greedily pruned using a 1-Mb distance around each retained lead. Each lead defined a +/-1-Mb window, and overlapping windows were merged. This distance-based procedure produced 60 leads and 57 merged loci. It is not LD clumping and is not intended to reproduce exactly the 59 loci reported in the source publication.

## MAGMA gene analysis

MAGMA v1.10 was run with NCBI37.3 gene locations and the 1000 Genomes European `g1000_eur` LD panel. The P-value input contained 4,915,033 unique variants; MAGMA used 4,882,602, mapped 1,931,992 to genes and tested 17,316 genes. Annotation used default gene coordinates with no added flanking window. MAGMA reported 133 synonymous SNP-pair warnings. Bonferroni significance was P <0.05/17,316 = 2.89 x 10^-6; BH FDR <0.05 and nominal P <1 x 10^-4 defined FDR and suggestive sets. This yielded 94 Bonferroni-significant, 369 FDR-significant and 187 suggestive genes. Top-50 and top-100 sets were fixed by MAGMA rank. Because the GWAS was trans-ancestry and the LD panel European, gene statistics require ancestry-reference caution.

## GSE231569 reconstruction, QC and annotation audit

Four papillary donors (two healthy and two stone-disease samples) were retained; the study medulla sample was excluded. Seurat objects were created with `min.cells=3` and `min.features=200`. Nuclei were filtered to 200-6,000 detected features and <=20% mitochondrial reads, log-normalized with `NormalizeData`, represented by 3,000 variable features, scaled, reduced with 30 PCs, clustered with dimensions 1-25 at resolution 0.4 and embedded by UMAP with dimensions 1-25.

Initial labels were obtained by broad marker-score comparison and then audited against cluster markers. The final analysis harmonized six compartments: collecting-duct principal (34,849 nuclei; four donors), fibroblast/stromal (3,422; four), endothelial (2,558; four), injured/undifferentiated epithelial (2,361; four), Loop/TAL (540; four) and perivascular/mural (148; three). Loop/TAL identity required concordant expression of established markers including UMOD, SLC12A1, CLDN10, KCNJ1 and CLDN16. Loop/TAL donor counts were 29, 244, 263 and 4 nuclei; donor-level results are therefore descriptive and the fourth donor is a sparse-compartment limitation.

## Single-nucleus MAGMA scoring and benchmarks

For a gene set G and nucleus i, the score was `S_i = mean(x_ig)` for detected genes g in G, where x is the Seurat RNA `data` layer after log normalization. This is not Seurat `AddModuleScore` and does not use per-gene z scores. Scores were averaged by audited cell type and donor-cell-type combination.

The expression-stratified size-matched random-set benchmark used a universe comprising genes with finite positive global mean expression. Genes were divided into 20 quantile bins of global mean expression. For each real gene, one replacement was sampled from the same bin; 1,000 sets were generated with seed 20260617. The empirical percentile was the proportion of random scores below or equal to the observed score. Loop/TAL percentiles were 0.998 for top 50, 1.000 for top 100, 1.000 for suggestive genes and 0.968 for FDR genes.

For the locus-balanced top-50 set, one highest-ranked gene was retained per locus group, with MAGMA-only genes assigned individual groups. Full and conservative analyses used 1,000 matched random sets; the conservative analysis excluded low/exploratory and immune-review labels. The conservative Loop/TAL percentile was 0.987. Leave-one-locus-out analysis removed each locus group and repeated 1,000-set benchmarking; retained Loop/TAL percentiles ranged from 0.997 to 1.000.

## P1 candidate evidence scoring

UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were assessed as priority-1 candidates. The evidence table records MAGMA membership, expression/detection in audited compartments, Loop/TAL rank and specificity, donor detection, correlation with a TAL reference program, GSE73680 availability and biological role. Evidence-badge summaries are descriptive integrations, not a combined probability or causal score.

## GSE73680 expression reconstruction and metadata

All 62 local supplementary feature files were parsed. Agilent feature expression was selected in priority order from `gProcessedSignal`, `gBGSubSignal`, `gMeanSignal` or `gMedianSignal`. Control and dark-corner features were excluded. Direct gene-symbol identifiers were retained; 32,055 of 44,661 features mapped and 12,606 remained unmapped because platform annotation was unavailable. The prespecified collapse rule retained the probe with maximum mean expression for duplicate symbols, although no multi-probe genes remained after direct-symbol filtering. Because direct gene-symbol parsing removed platform probes lacking unambiguous symbols, the resulting gene-level matrix should be interpreted as a conservative direct-symbol reconstruction rather than a full platform reannotation. The reconstructed scale was continuous and appeared normalized, but platform-wide normalization could not be independently verified. Values were transformed as log2(expression + 1). No reliable batch variable was available.

Curated filename and sample metadata assigned 27 control/adjacent and 28 plaque/stone-associated papillary samples. These 55 samples came from 29 patients. Twenty-six patients contributed both groups and formed the paired analysis; three patients contributed one included group and informed the primary duplicate-correlation model only.

## GSE73680 gene and module models

The primary gene-level model was `expression ~ group_curated`, with `patient_id` as the blocking factor in limma `duplicateCorrelation`; the reported coefficient was plaque/stone-associated papilla versus control/adjacent. The paired sensitivity first averaged repeated samples within patient and group, then calculated plaque/stone-associated minus control/adjacent values and used paired t tests. BH adjustment was applied across the six P1 genes separately from module analyses.

For each gene, expression was standardized across samples. A module score was the sample mean of available gene z scores. Full modules required >=70% gene detection; modules at 40-70% were exploratory and lower-coverage modules were not tested. Patient-level scores averaged samples within patient and group. Paired module shifts were tested among 26 paired patients, with BH adjustment across the tested module family. Top-50, top-100, FDR and suggestive modules each had q=0.0499; the P1 module had q=0.299 and was not FDR-supported.

The size-matched benchmark drew 1,000 random sets. Leave-one-gene-out tests, paired-direction consistency and removal of all six P1 genes assessed dependence on individual genes. After P1 removal, the MAGMA modules retained q=0.0251. Expression-matched random sets used 20 mean-expression bins and 1,000 iterations (seed 73682); because this stricter analysis did not provide primary support, it was retained only as a sensitivity boundary.

## Functional enrichment and curated contexts

MAGMA-tested genes were mapped to Entrez identifiers, yielding a 15,421-gene GO universe. `clusterProfiler::enrichGO` tested Biological Process terms and calculated BH-adjusted P values. Main-figure display required FDR <0.10 and at least two query genes. Terms were ordered by adjusted P value and Count; subsequent terms with gene-overlap Jaccard >=0.70 to an already retained term were removed.

Prespecified contextual sets covered TAL transport, calcium/ion handling, epithelial junction and injury/remodeling. They were assembled from canonical renal marker genes used in the annotation and disease-context workflows and are listed in the source scripts and supplementary source tables. Hypergeometric tests used unique MAGMA symbols as the universe and applied BH correction across all tested query-set by context-set combinations. These sets provide interpretation, not externally validated pathway activity.

## Injury/remodeling coupling

Marker programs were injury/remodeling (`SPP1,MMP7,MMP9,GPNMB,COL1A1,COL1A2,FN1,VIM,HAVCR1,LCN2`), inflammation/immune (`CCL2,CCL7,CXCL8,IL1B,TNF,CD68,LST1`), fibrosis/ECM (`COL1A1,COL3A1,DCN,LUM,ACTA2,TAGLN`) and epithelial injury (`HAVCR1,LCN2,VCAM1,SOX9,PROM1`). Scores were means of row-wise z-scored genes.

Spearman correlations were computed for all samples, for the 26 paired patient deltas, and for residuals from separate linear models of risk-module score and injury score on group plus patient factor. BH adjustment was performed across module-program pairs within each analysis. Paired-delta and residual results were prioritized. Leave-one-patient-out summaries and sample-level correlations were sensitivity analyses. Correlation denotes coupling, not mediation or mechanism.

## Resource-limited analyses

FUSION/S-PrediXcan TWAS required tissue weights and covariance/LD resources that were incomplete. SMR/coloc required kidney eQTL summaries or BESD/ESI/EPI inputs and suitable LD that were incomplete. The GTEx atlas was reviewed as a potential eQTL resource, but available kidney context did not complete the prespecified papillary analysis. Spatial candidate datasets lacked a complete audited combination of matrices, barcodes/features, coordinates, scalefactors and images. These workflows were not used for scientific claims.

## Software

Captured versions were R 4.4.3, MAGMA 1.10, Python 3.9.6, data.table 1.18.2.1, ggplot2 4.0.2, Seurat 5.4.0, clusterProfiler 4.14.6, cowplot 1.2.0 and svglite 2.2.2. Earlier preprocessing logs and scripts remain the authoritative record where package-version capture was unavailable.

## Reproducibility and source-file inventory

Supplementary Tables S1-S13 contain the frozen analytical outputs used for manuscript reporting and sensitivity boundaries. Figure source-data records F1-F5 map each main figure to its source tables and generating scripts through `results/tables/figure_source_data_manifest_v0.1.tsv`. The current supplementary table inventory is `results/tables/supplementary_tables_manifest_v0.2.tsv`. Checksums were generated for the local release scaffold before repository deposition. Repository URL to be added before submission; no public identifier is assigned in this internal precheck package.
