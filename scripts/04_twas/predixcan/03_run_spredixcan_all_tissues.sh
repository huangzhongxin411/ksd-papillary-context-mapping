#!/usr/bin/env bash
set -euo pipefail

tissues=(
  kidney_cortex
  whole_blood
  artery_aorta
  artery_tibial
  adipose_subcutaneous
  liver
  colon_transverse
  small_intestine_terminal_ileum
)

for tissue in "${tissues[@]}"; do
  echo "[RUN] $tissue"
  bash scripts/04_twas/predixcan/02_run_spredixcan_per_tissue.sh "$tissue"
done
