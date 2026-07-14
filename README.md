# Post-GWAS renal papillary context mapping of kidney stone disease genetic risk

## Repository status

This is the v3.2 submission-aligned reproducibility package for a BMC Genomics Research article. It contains final figures, supplementary figures, supplementary tables, panel source data, selected reproducibility scripts, manifests and release metadata. The scientific manuscript is content-locked; the repository archive DOI remains pending until the approved v1.0.0 release.

## Study purpose and claim boundary

The study provides post-GWAS functional interpretation and renal papillary context mapping for kidney stone disease genetic risk. MAGMA defines an EUR-LD-reference-based genetic-priority layer. Donor-level single-nucleus analysis supports a Loop/TAL-associated context; spatial transcriptomics provides supplementary broad-compartment tissue projection; GTEx Kidney_Cortex TWAS is proxy annotation; and paired bulk expression provides disease-context and tissue-state sensitivity.

The repository does **not** establish causal genes, causal cell types, plaque-specific localization, papilla-specific regulatory effects, therapeutic targets or claim-grade SMR/colocalization. R1-R6 are reporting groups, not causal tiers.

## Final v3.2 outputs

- Main Figures: Figure 1, Figure 2, Figure 3 and Figure 4 under `figures/`.
- Supplementary Figures: Supplementary Figure S1, S2 and S3 under `supplementary_figures/`.
- Supplementary Tables: Supplementary Table 1 through 6 under `supplementary_tables/`, with an XLSX workbook and old-to-final crosswalk.
- Spatial projection: ten complete sections (10 total), comprising four GSE206306 and six GSE231630 sections. Projection is descriptive and broad-compartment only; no lesion ROI or plaque-specific localization is claimed.

## Public data sources

- KSD GWAS summary statistics: https://doi.org/10.5281/zenodo.14790324
- GWAS Catalog: GCST90652506
- GEO: GSE231569, GSE73680, GSE206306 and GSE231630
- GTEx v8 / PredictDB Kidney_Cortex MASHR models from the original distribution resource

Raw third-party data, LD reference panels and prediction-model files are not redistributed here. Review their original terms before use.

## Reproducing major analyses

1. Follow Phase 1 paths in `REPRODUCIBILITY.md` for GWAS QC, locus reconstruction and MAGMA freezing.
2. Follow Phase 2 for donor-level GSE231569 module scoring and sensitivity analyses.
3. Follow Phase 3 for broad-compartment spatial projection across the ten complete sections.
4. Follow Phase 4 for Kidney_Cortex proxy TWAS auditing and the 235-gene R1-R6 reporting model.
5. Follow Phase 5 for paired GSE73680 disease-context and tissue-state sensitivity.
6. Follow Phase 6 for final figure, supplement and manuscript-package materialization.

Exact script/output links are in `scripts/SCRIPT_MAP_v3.2.tsv`. Panel source files are indexed by `source_data/Source_Data_manifest_v3.2.tsv`. All release files are listed in `MANIFEST.tsv` and hashed in `CHECKSUMS.sha256`.

## Directory map

- `figures/`: final main Figure 1-4 files.
- `supplementary_figures/`: final Supplementary Figure S1-S3 files.
- `supplementary_tables/`: final Supplementary Table 1-6 files, workbook and crosswalk.
- `supplementary_materials/`: final legends and captions.
- `source_data/`: panel source data and v3.2 source-data manifest.
- `scripts/`: selected Phase 1-6 scripts and exact mapping.
- `manifests/`: submission-package provenance manifests.
- `archive/`: deprecation documentation and pre-v3.2 remote inventory only. Archived items are not current evidence outputs.

## Citation

Until the repository release is minted, cite the associated manuscript title and repository URL: https://github.com/huangzhongxin411/ksd-papillary-context-mapping. After the v1.0.0 Zenodo archive is created, replace this instruction with the minted DOI and use `CITATION.cff`.

Repository archive DOI to be inserted after v1.0.0 release.
