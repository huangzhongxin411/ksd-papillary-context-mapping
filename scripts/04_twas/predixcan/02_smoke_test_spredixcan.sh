#!/usr/bin/env bash
set -euo pipefail

TISSUE="${1:-kidney_cortex}"

bash scripts/04_twas/predixcan/02_run_spredixcan_per_tissue.sh "$TISSUE"

out="results/twas/predixcan/per_tissue/${TISSUE}.spredixcan.tsv"
log="results/twas/predixcan/logs/${TISSUE}.spredixcan.log"

if [[ ! -s "$out" ]]; then
  echo "[FAILED] S-PrediXcan smoke test output missing: $out"
  if [[ -s "$log" ]]; then
    tail -n 80 "$log"
  fi
  exit 1
fi

echo "[OK] S-PrediXcan smoke test passed: $out"
sed -n '1,5p' "$out"
