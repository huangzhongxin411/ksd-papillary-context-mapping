# Post-GWAS mapping of kidney stone disease genetic risk to renal papillary contexts

## Repository status

This branch prepares the manuscript v3.3-synchronized reproducibility package for a future human-approved GitHub Release v1.0.1. GitHub Release v1.0.0 remains an immutable pre-Zenodo checkpoint. No repository archive DOI has yet been minted.

## Study purpose and claim boundary

The study provides post-GWAS functional interpretation and renal papillary context mapping for kidney stone disease genetic risk. MAGMA defines an EUR-LD-reference-based genetic-priority layer. Donor-level single-nucleus analysis supports a Loop/TAL-associated context; spatial transcriptomics provides supplementary broad-compartment tissue projection; GTEx Kidney_Cortex TWAS is proxy annotation; and paired bulk expression provides disease-context and tissue-state sensitivity.

The repository does **not** establish causal genes, causal cell types, plaque-specific localization, papilla-specific regulatory effects or therapeutic targets. R1-R6 are reporting groups, not causal tiers.

## Current v3.3 outputs

- Clean manuscript: Markdown, DOCX and rendered PDF under `manuscript/`.
- Main Figures: Figure 1 through Figure 4 under `figures/`.
- Supplementary Figures: Supplementary Figure S1 through S3 under `supplementary_figures/`.
- Supplementary Tables: Supplementary Table 1 through 6 under `supplementary_tables/`, with an XLSX workbook and old-to-final crosswalk.
- Spatial projection: ten complete sections (10 total), comprising four GSE206306 and six GSE231630 sections. Projection is descriptive and broad-compartment only; no lesion ROI or plaque-specific localization is claimed.

## Public data sources

- KSD GWAS summary statistics: https://doi.org/10.5281/zenodo.14790324
- GWAS Catalog: GCST90652506
- GEO: GSE231569, GSE73680, GSE206306 and GSE231630
- GTEx v8 / PredictDB Kidney_Cortex MASHR models from the original distribution resource

Raw third-party data, LD reference panels and prediction-model files are not redistributed here. Review their original terms before use.

## Reproducing major analyses

Follow `REPRODUCIBILITY.md` for the Phase 1-7 workflow. Exact script/output links are in `scripts/SCRIPT_MAP_v3.3.tsv`. Panel source files are indexed by `source_data/Source_Data_manifest_v3.3.tsv`. All repository files are listed in `MANIFEST.tsv` and hashed in `CHECKSUMS.sha256`.

## Directory map

- `manuscript/`: submission-clean v3.3 manuscript and revision log.
- `figures/`: current main Figure 1-4 files and figure manifest.
- `supplementary_figures/`: current Supplementary Figure S1-S3 files and figure manifest.
- `supplementary_tables/`: current Supplementary Table 1-6 files, workbook and crosswalk.
- `supplementary_materials/`: current legends and captions.
- `source_data/`: panel source data, locked supporting tables and v3.3 source-data manifest.
- `scripts/`: selected Phase 1-7 scripts and exact mapping.
- `manifests/`: submission-package provenance manifests.
- `archive/`: deprecation documentation and pre-v3.2 remote inventory only.

## Citation

Until a Zenodo archive is minted, cite the associated manuscript title and this repository URL: https://github.com/huangzhongxin411/ksd-papillary-context-mapping. After human approval of v1.0.1 and Zenodo publication, update `CITATION.cff`, `.zenodo.json` and the manuscript Data Availability statement with the minted DOI.

Repository archive DOI: pending; do not cite an unminted DOI.
