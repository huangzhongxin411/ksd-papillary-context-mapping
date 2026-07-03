#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/processed/magma_input results/tables

cp data/processed/magma/2025_trans_ancestry/ksd_2025.snploc \
  data/processed/magma_input/ksd_2025.snploc
cp data/processed/magma/2025_trans_ancestry/ksd_2025.pval \
  data/processed/magma_input/ksd_2025.pval

python3 scripts/03_magma_fuma/deduplicate_magma_input.py \
  --snploc data/processed/magma_input/ksd_2025.snploc \
  --pval data/processed/magma_input/ksd_2025.pval \
  --out-snploc data/processed/magma_input/ksd_2025.dedup.snploc \
  --out-pval data/processed/magma_input/ksd_2025.dedup.pval \
  --qc-out results/tables/magma_input_qc.tsv

printf 'wrote\t%s\n' data/processed/magma_input/ksd_2025.snploc
printf 'wrote\t%s\n' data/processed/magma_input/ksd_2025.pval
printf 'wrote\t%s\n' data/processed/magma_input/ksd_2025.dedup.snploc
printf 'wrote\t%s\n' data/processed/magma_input/ksd_2025.dedup.pval
printf 'wrote\t%s\n' results/tables/magma_input_qc.tsv
