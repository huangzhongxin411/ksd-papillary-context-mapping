# Post-GWAS mapping of kidney stone disease genetic risk to renal papillary contexts

## Repository status

This branch synchronizes the manuscript v3.4 human-review package for a future human-approved GitHub Release v1.0.1. GitHub Release v1.0.0 remains an immutable pre-Zenodo checkpoint. GitHub Release v1.0.1 has not been created, and no repository archive DOI has been minted.

## Study purpose and claim boundary

The study provides post-GWAS functional interpretation and renal papillary context mapping for kidney stone disease genetic risk. MAGMA defines an EUR-LD-reference-based genetic-priority layer. Donor-level single-nucleus analysis supports a Loop/TAL-associated context; spatial transcriptomics provides supplementary broad-compartment tissue projection; GTEx Kidney_Cortex TWAS is proxy annotation; and paired bulk expression provides disease-context and tissue-state sensitivity.

The repository does **not** establish causal genes, causal cell types, plaque-specific localization, papilla-specific regulatory effects or therapeutic targets. R1-R6 are reporting groups, not causal tiers.

## Current v3.4 review package

- Manuscript: v3.4 Markdown, DOCX, rendered PDF and revision change log under `manuscript/`.
- Changed from v3.3: targeted Abstract/Discussion/Limitations wording, Figure 3 layout/title policy and Supplementary Figure S2 overlay contrast/title policy.
- Retained unchanged: Figure 1, Figure 2, Figure 4, Supplementary Figure S1, Supplementary Figure S3 and Supplementary Tables 1-6.
- Figure 1 remains an unlabeled GWAS quality-control Manhattan plot.
- Spatial projection remains five pages across ten complete sections: four GSE206306 and six GSE231630. It is descriptive and broad-compartment only; no lesion ROI, Loop/TAL enrichment or plaque localization is claimed.

## Public data sources

- KSD GWAS summary statistics: https://doi.org/10.5281/zenodo.14790324
- GWAS Catalog: GCST90652506
- GEO: GSE231569, GSE73680, GSE206306 and GSE231630
- GTEx v8 / PredictDB Kidney_Cortex MASHR models from the original distribution resource

Raw third-party data, LD reference panels and prediction-model files are not redistributed here. Review their original terms before use.

## Reproducing major analyses

Follow `REPRODUCIBILITY.md` for the Phase 1-8 workflow. No biological analysis was rerun for v3.4. The v3.3 script and source-data maps remain the authoritative locked-analysis provenance because the v3.4 changes are manuscript language and figure rendering only. All repository files are listed in `MANIFEST.tsv` and hashed in `CHECKSUMS.sha256`.

## Directory map

- `manuscript/`: v3.4 human-review manuscript plus retained v3.3 checkpoint files.
- `figures/`: current canonical Figure 1-4 files and versioned figure manifests.
- `supplementary_figures/`: current canonical Supplementary Figure S1-S3 files and versioned manifests.
- `supplementary_tables/`: unchanged Supplementary Table 1-6 files and workbook.
- `supplementary_materials/`: current legends and captions.
- `source_data/`: unchanged panel source data and locked supporting tables.
- `scripts/`: selected Phase 1-7 analysis/package scripts; no analysis script changed for v3.4.
- `submission_package_v3.4_draft/`: lightweight human-review index only.
- `archive/`: deprecation documentation and pre-v3.2 remote inventory only.

## Citation

Until a Zenodo archive is minted, cite the associated manuscript title and this repository URL: https://github.com/huangzhongxin411/ksd-papillary-context-mapping. After human approval of v1.0.1 and Zenodo publication, update `CITATION.cff`, `.zenodo.json` and the manuscript Data Availability statement with the real minted DOI.

Repository archive DOI: pending; do not cite an unminted DOI.
