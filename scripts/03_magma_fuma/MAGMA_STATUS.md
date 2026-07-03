# MAGMA Status

Checked on 2026-06-17.

## Status

Phase 2A MAGMA gene-based analysis is complete.

The previous blocker is resolved:

- `external/magma/magma` is present and executable.
- `external/reference/1000G_EUR/g1000_eur.{bed,bim,fam}` is present.
- `external/reference/gene_loc/NCBI37.3.gene.loc` is present.

MAGMA version is recorded in:

- `results/logs/magma_version.txt`

## Prepared Inputs

MAGMA input files:

- `data/processed/magma_input/ksd_2025.dedup.snploc`
- `data/processed/magma_input/ksd_2025.dedup.pval`

Input QC:

- `results/tables/magma_input_qc.tsv`
- Rows/unique SNPs retained for MAGMA input: 4,915,033
- Duplicate SNPs detected in snploc/pval: 0 / 0

## Gene-Based Run

Run command:

```bash
bash scripts/03_magma_fuma/run_magma_gene_analysis.sh
```

Primary outputs:

- `results/magma/2025_trans_ancestry/ksd_2025.genes.annot`
- `results/magma/2025_trans_ancestry/ksd_2025.genes.out`
- `results/magma/2025_trans_ancestry/ksd_2025.genes.raw`
- `results/magma/2025_trans_ancestry/ksd_2025.log`

Run summary:

- 19,427 gene locations read.
- 4,915,033 SNP locations read.
- 1,931,992 SNPs mapped to at least one gene.
- 17,316 genes tested with valid SNPs in reference genotype data.

## Post-Processing

Post-processing command:

```bash
Rscript scripts/03_magma_fuma/postprocess_magma_genes.R
```

Generated:

- `results/tables/magma_qc_summary.tsv`
- `results/tables/magma_gene_set_summary.tsv`
- `results/tables/magma_genes.tsv`
- `results/tables/magma_vs_locus_overlap.tsv`
- `results/gene_sets/magma_top50.txt`
- `results/gene_sets/magma_top100.txt`
- `results/gene_sets/magma_top200.txt`
- `results/gene_sets/magma_suggestive_p1e4.txt`
- `results/gene_sets/magma_fdr05.txt`

Freeze-summary command:

```bash
Rscript scripts/03_magma_fuma/make_magma_freeze_summaries.R
```

## Single-Cell Reprojection and Sensitivity

MAGMA-GSE231569 projection:

```bash
Rscript scripts/07_scrna_gene_mapping/magma_scrna_projection_benchmark.R
```

MAGMA top50 sensitivity:

```bash
Rscript scripts/07_scrna_gene_mapping/magma_driver_balanced_sensitivity.R
```

Candidate tiering:

```bash
Rscript scripts/03_magma_fuma/make_candidate_gene_tiers_v0.1.R
```

Key results:

- MAGMA top50 TAL benchmark percentile: 0.998.
- MAGMA top100 TAL benchmark percentile: 1.000.
- MAGMA top200 TAL benchmark percentile: 0.999.
- MAGMA FDR 0.05 TAL benchmark percentile: 0.968.
- MAGMA locus-balanced top50 TAL percentile: 0.988 full, 0.987 conservative.
- Leave-one-locus-out retained TAL percentile above 0.90 for all tested locus groups.

## Interpretation

The MAGMA result supports a TAL-associated KSD genetic-prioritization model. It supersedes the frozen locus-based projection as the current main line, while TWAS/SMR-colocalization remains required before causal gene-expression mediation claims.

## Caution

This run uses `NCBI37.3.gene.loc` and 1000G EUR reference files. Because the primary GWAS is trans-ancestry EUR + EAS, this should be described as a formal first-pass MAGMA prioritization, not a complete ancestry-matched LD sensitivity analysis.
