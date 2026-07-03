# Reproducibility guide

## Software overview

The workflow uses R, Python, shell scripts, MAGMA v1.10, Seurat-based single-nucleus/spatial processing, S-PrediXcan with GTEx v8 Kidney_Cortex MASHR models, and tabular audit utilities. MAGMA used NCBI build 37 gene locations and the 1000 Genomes European LD reference. It is therefore an EUR-LD-reference-based prioritization layer, not ancestry-generalizable fine mapping.

Exact R `sessionInfo()` and a frozen Python dependency lock were not available for every historical stage. `environment/ENVIRONMENT_STATUS.md` records this release-candidate limitation. Users should capture fresh environment details when rerunning the public scripts.

## Staged workflow

1. GWAS QC and locus reconstruction (`scripts/01_gwas_qc`, `scripts/02_locus_mapping`).
2. MAGMA gene prioritization (`scripts/03_magma_fuma`).
3. Kidney_Cortex proxy TWAS and feasibility audits (`scripts/04_twas`, `scripts/05_smr_coloc`).
4. GSE231569 single-nucleus processing and donor-level context analyses (`scripts/06_scrna_processing`, `scripts/07_scrna_gene_mapping`).
5. GSE73680 paired bulk and adjustment analyses (`scripts/08_plaque_context`, `scripts/09_bulk_plaque_context`).
6. GSE206306 spatial processing and complexity-adjusted context analysis (`scripts/11_spatial`).
7. Revision-stage evidence tables, source data and figures (`scripts/revision_utils`).

Scripts with hard-coded local absolute paths were excluded from the release candidate and are listed in the external repository inventory. Run from the repository root and inspect each stage's configuration before execution.

## Key-number checks

| Check | Expected value | Primary evidence location |
|---|---:|---|
| GWAS QC-passed rows | 4,915,033 | `supplementary_tables/S1_gwas_qc_summary.tsv` |
| MAGMA-tested genes | 17,316 | `supplementary_tables/S3_magma_gene_results_and_gene_sets.tsv` |
| Bonferroni genes | 94 | `supplementary_tables/S3_magma_gene_results_and_gene_sets.tsv` |
| Papillary nuclei | 43,878 | `source_data/stage4C2R_draft_figures_v0.2/figure2_panelA_source.tsv` |
| Loop/TAL nuclei | 540 | `source_data/stage4C2R_draft_figures_v0.2/figure2_panelA_source.tsv` |
| Paired patients | 26 | `source_data/stage5C1_gse73680_figure4_draft/figure4_panelA_source.tsv` |
| Spatial sections / spots | 5 / 7,747 | `source_data/stage6C_spatial_twas_figure5_draft/figure5_panelA_source.tsv` |

## Boundaries

The workflow is observational and context-mapping focused. It does not demonstrate causal genes, a causal cell type, plaque-specific localization, papilla-specific regulatory effects, SMR/colocalization support or therapeutic targets. Raw public data and ancestry-specific LD resources are not bundled.
