#!/usr/bin/env bash
set -euo pipefail

files=(
  "data/processed/twas_input/ksd_2025_for_twas.tsv.gz"
  "data/processed/twas_input/ksd_2025_for_twas_variant_flags.tsv.gz"
  "results/tables/twas_input_qc_report.tsv"
  "docs/twas_input_notes.md"
  "results/tables/magma_genes.tsv"
  "results/tables/candidate_gene_tiers_v0.1.tsv"
)

for f in "${files[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "[MISSING] $f"
    exit 1
  else
    echo "[OK] $f"
  fi
done

gzip -dc data/processed/twas_input/ksd_2025_for_twas.tsv.gz | sed -n '1,3p'
echo "Phase 2B input check passed."
