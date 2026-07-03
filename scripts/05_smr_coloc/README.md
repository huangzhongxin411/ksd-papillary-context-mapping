# 05 SMR / Colocalization

## Goal

Strengthen candidate gene prioritization by testing whether GWAS and eQTL/TWAS signals share compatible genetic evidence.

## Candidate Criteria

- SMR significant after correction.
- HEIDI P > 0.01 as evidence against obvious heterogeneity.
- coloc PP4 > 0.75 for support and PP4 > 0.9 for strong support.

## Output

- `results/tables/smr_heidi_results.tsv`
- `results/tables/coloc_results.tsv`
- High-confidence subset for Tier 1 candidate genes.

