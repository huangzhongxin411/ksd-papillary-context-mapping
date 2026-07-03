#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/gwas/2025_trans_ancestry

url='https://zenodo.org/records/14790324/files/meta_sumstats?download=1'
out='data/raw/gwas/2025_trans_ancestry/meta_sumstats'

curl -L --fail --show-error --retry 5 --retry-delay 5 \
  -C - \
  "$url" \
  -o "$out"

printf 'downloaded\t%s\n' "$out"
if command -v md5sum >/dev/null 2>&1; then
  md5sum "$out"
elif command -v md5 >/dev/null 2>&1; then
  md5 "$out"
fi
