#!/usr/bin/env bash
set -euo pipefail

TISSUE="${1:-kidney_cortex}"
CHR="${2:-22}"

bash scripts/04_twas/fusion/03_run_fusion_per_tissue.sh "$TISSUE" "$CHR"

out="results/twas/fusion/per_tissue/${TISSUE}/${TISSUE}.chr${CHR}.fusion.tsv"
log="results/twas/fusion/logs/${TISSUE}.chr${CHR}.log"

if [[ ! -s "$out" ]]; then
  echo "[FAILED] FUSION smoke test output missing: $out"
  if [[ -s "$log" ]]; then
    tail -n 80 "$log"
  fi
  exit 1
fi

echo "[OK] FUSION smoke test passed: $out"
sed -n '1,5p' "$out"
