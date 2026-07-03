#!/usr/bin/env bash
set -euo pipefail

gwas='data/raw/gwas/2025_trans_ancestry/meta_sumstats'

if [[ ! -s "$gwas" ]]; then
  printf 'missing_or_empty\t%s\n' "$gwas" >&2
  printf 'Download from: https://zenodo.org/records/14790324/files/meta_sumstats?download=1\n' >&2
  exit 1
fi

python scripts/01_gwas_qc/inspect_gwas_columns.py "$gwas"

