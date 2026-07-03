# Post-GWAS renal papillary context mapping of kidney stone disease genetic risk

This repository contains code, derived result tables, figure source data, supplementary tables, manifests and audit logs supporting a BMC Genomics manuscript on post-GWAS renal papillary context mapping. MAGMA-prioritized kidney stone disease genes map to a Loop/TAL-associated papillary context and an injury/remodeling-associated bulk papillary disease background across conservative single-nucleus, bulk and spatial analyses.

## Claim boundary

The repository supports context mapping, not causal-gene validation. It does not establish validated disease genes, a causal cell type or niche, plaque-specific localization, plaque nucleation, papilla-specific TWAS regulation, SMR/colocalization support or therapeutic targets.

## Repository structure

- `scripts/`: analysis and release-support scripts without local absolute paths
- `config/`: dataset and sample configuration files
- `results/tables/`: audited claim-support tables
- `source_data/`: figure panel source data and manifests
- `supplementary_tables/`: Supplementary Tables 1-13
- `figures/`: working PDF, PNG and SVG candidates pending final author artwork approval
- `logs/`: selected reproducibility logs
- `environment/`: software and environment documentation
- `docs/`: supplementary and release documentation

## Public data sources

- Cao et al. KSD GWAS summary statistics: https://doi.org/10.5281/zenodo.14790324
- GWAS Catalog study record: GCST90652506
- GEO: GSE231569, GSE73680 and GSE206306
- GTEx v8 / PredictDB Kidney_Cortex MASHR models

Raw third-party data are not redistributed here. Retrieve them from their original repositories and review their applicable terms.

## Reproduction outline

1. Download the public inputs listed in `DATA_AVAILABILITY.md`.
2. Review `config/` and run scripts in the staged order described in `REPRODUCIBILITY.md`.
3. Compare outputs against `MANIFEST.tsv` and the source-data manifests.
4. Regenerate the derived source tables and figures using the figure-stage scripts.

## Expected key outputs

The package documents 4,915,033 GWAS QC-passed rows, 17,316 MAGMA-tested genes, 94 Bonferroni genes, 43,878 nuclei including 540 Loop/TAL nuclei, 26 paired patients, and five spatial sections comprising 7,747 spots.

## Citation and licenses

Use `CITATION.cff` when citing this repository. Code is licensed under MIT (`LICENSE-CODE`); repository documentation and generated tabular data are licensed under CC BY 4.0 (`LICENSE-DATA`). Third-party inputs remain subject to their original terms.

## Corresponding authors and contacts

- Xiaolu Duan: 94302304@qq.com
- Guohua Zeng: gzgyzgh@vip.sina.com

Repository URL: https://github.com/huangzhongxin411/ksd-papillary-context-mapping

Archive DOI: [TO FILL: Zenodo DOI]
