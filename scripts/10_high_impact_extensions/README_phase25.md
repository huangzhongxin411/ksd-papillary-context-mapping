# Phase 25 high-journal enhancement modules

This folder prepares the three requested extension layers without changing the frozen manuscript claims.

## Priority order

1. TWAS first.
2. Targeted coloc second.
3. Spatial projection only if complete Visium resources are available.

## TWAS resource placement

Preferred S-PrediXcan paths:

- Model databases: `external/twas/predixcan/models/`
- Covariance files: `external/twas/predixcan/covariance/`
- Software: `external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py`

Required Python packages for S-PrediXcan:

- numpy
- scipy
- pandas

The current run found MetaXcan code but not the required Python packages, model `.db` files, or covariance files.

## FUSION resource placement

Optional FUSION paths:

- `external/twas/fusion/software/FUSION.assoc_test.R`
- `external/twas/fusion/weights/`
- `external/twas/fusion/ref_ld/`

The current run found a full PLINK 1000G EUR reference at `external/reference/1000G_EUR/g1000_eur.*`, but not FUSION chr-split `1000G.EUR.{1..22}` files.

## Generated status files

- `results/twas/twas_resource_status.tsv`
- `results/twas/twas_results.tsv`
- `results/twas/twas_magma_overlap.tsv`
- `results/twas/twas_download_manifest_v0.1.tsv`
- `results/coloc/priority_locus_coloc_plan_v0.1.tsv`
- `results/spatial/spatial_projection_resource_status_v0.2.tsv`
- `docs/phase25_high_journal_extension_status_v0.1.md`

## Claim boundary

Only FDR-supported TWAS genes may be promoted as enhanced genetic evidence. Missing TWAS, coloc or spatial resources are resource limitations, not negative biological evidence.
