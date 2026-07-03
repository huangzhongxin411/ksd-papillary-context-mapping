# Reproducibility Summary v1

## Scope

This document records the Phase 18 frozen analysis and figure package supporting the manuscript claim that MAGMA-prioritized KSD genes map to a Loop/TAL-associated papillary cellular context and show module-level plaque/stone disease-context association.

## Input Datasets

| Evidence layer | Dataset/input | Local role |
|---|---|---|
| Genetic prioritization | 2025 trans-ancestry KSD GWAS summary statistics | GWAS QC, lead loci and MAGMA gene-level prioritization |
| Single-nucleus context | GEO GSE231569 renal papillary snRNA-seq | Audited cell annotation, module projection, cell-context benchmarking and P1 expression summaries |
| Disease context | GEO GSE73680 Agilent plaque/stone papilla bulk expression | Patient-aware duplicateCorrelation analysis, 26-patient paired sensitivity, module benchmarks and injury coupling |
| Gene annotation/reference | GRCh37/hg19-compatible MAGMA gene locations and 1000 Genomes EUR LD reference | SNP-to-gene aggregation and gene-level testing |
| Functional annotation | GO Biological Process via clusterProfiler/org.Hs.eg.db | Background-audited, redundancy-reduced functional interpretation |

## Key Scripts

- GWAS and MAGMA: `scripts/00_download_phase1_gwas.sh` and the frozen Phase 1/Phase 2 MAGMA workflow documented in `docs/magma_analysis_notes.md`
- GSE231569 reconstruction/audit: `scripts/06_reconstruct_gse231569_raw.sh` and scripts under `scripts/07_scrna_gene_mapping/`
- GSE73680 reconstruction/analysis: scripts under `scripts/08_plaque_context/`
- Figure 1: `scripts/09_manuscript/figure1_evidence_map_real_manhattan_restored.R`
- Figure 2: `scripts/09_manuscript/figure2_gwas_magma_scrna_localization_v1.0.3.R`
- Figure 3 Phase 18 micro-polish: `scripts/09_manuscript/phase18_figure3_micro_polish.R`
- Figure 4 Phase 18 micro-polish: `scripts/09_manuscript/phase18_figure4_micro_polish.R`
- Figure 5 Phase 18 micro-polish: `scripts/09_manuscript/phase18_figure5_micro_polish.R`
- Figure freeze/package: `scripts/09_manuscript/phase18_freeze_main_figures.R`

## Key Frozen Outputs

- Main figure package: `results/figures/final_main_figures_v1/`
- Package manifest: `results/figures/final_main_figures_v1/MANIFEST.tsv`
- Package QC: `results/figures/final_main_figures_v1/figure_qc_summary.tsv`
- Final legends: `docs/figure1_legend_final.md` through `docs/figure5_legend_final.md`
- Result/figure/table crosswalk: `results/tables/result_figure_table_crosswalk_v1.tsv`
- Claim-boundary audit: `results/tables/claim_boundary_audit_v1.tsv` and `docs/claim_boundary_audit_v1.md`
- Hardened manuscript: `manuscript/manuscript_draft_v1.0_pre.md`

## Software Versions Captured in the Phase 18 Runtime

- R 4.4.3
- MAGMA 1.10
- Python 3.9.6
- data.table 1.18.2.1
- ggplot2 4.0.2
- cowplot 1.2.0
- svglite 2.2.2
- Seurat 5.4.0
- clusterProfiler 4.14.6

Package versions reflect the runtime used for the Phase 18 hardening/export pass. Earlier processing environments should additionally be reconstructed from existing logs and lock files where available.

## Missing External Resources

- **TWAS:** FUSION weights and LD reference; S-PrediXcan model database and covariance resources were incomplete.
- **SMR/coloc:** Kidney eQTL resources, SMR BESD/ESI/EPI files and coloc-ready eQTL summaries were incomplete.
- **Spatial transcriptomics:** Local matrix, barcode/feature, coordinate, scalefactor and image resources were incomplete for audited candidate datasets.

## Analyses Not Used for Claims

- TWAS association testing was not used as an evidence layer.
- SMR/coloc was not used as an evidence layer.
- Spatial transcriptomic validation was not used as an evidence layer.
- Cell-level UMAP score maps were used for visualization; donor-cell-type and benchmark summaries carry the interpretation.
- Expression-matched random benchmarking in GSE73680 was retained as a conservative sensitivity boundary rather than primary support.
- P1 single-gene responses were not used as disease validation; PKD2 remained nominal only.
- Functional enrichment and injury/remodeling correlations support interpretation and coupling, not pathway activity validation or causal mechanism.

## Reproduction Boundary

All frozen main figures are exported as PDF, SVG and 600-dpi PNG. The SVG files are checked for embedded raster/base64 image elements, and source paths plus MD5 checksums are recorded in the package manifest. Re-running external-resource-limited extensions requires obtaining the missing resources above and must not be treated as reproduction of a completed evidence layer.
