#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
SPREDIXCAN="${SPREDIXCAN:-external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py}"
GWAS="${GWAS:-data/processed/twas_input/ksd_2025_for_twas.tsv.gz}"

TISSUES=(
  kidney_cortex
  whole_blood
  artery_aorta
  artery_tibial
  adipose_subcutaneous
  liver
  colon_transverse
  small_intestine_terminal_ileum
)

mkdir -p results/twas/predixcan/per_tissue results/twas/predixcan/logs

if [[ ! -s "$SPREDIXCAN" ]]; then
  echo "[ERROR] SPrediXcan script not found: $SPREDIXCAN" >&2
  exit 1
fi

if [[ ! -s "$GWAS" ]]; then
  echo "[ERROR] TWAS GWAS input not found: $GWAS" >&2
  exit 1
fi

"$PYTHON_BIN" - <<'PY'
import numpy, scipy, pandas
print("python_dependencies_ok")
PY

for tissue in "${TISSUES[@]}"; do
  model=$(find external/twas/predixcan/models -iname "*${tissue}*.db" | head -n 1 || true)
  cov=$(find external/twas/predixcan/covariance -iname "*${tissue}*" | head -n 1 || true)
  if [[ -z "$model" || -z "$cov" ]]; then
    echo "[SKIP] ${tissue}: model or covariance missing" >&2
    continue
  fi

  out="results/twas/predixcan/per_tissue/${tissue}.spredixcan.tsv"
  log="results/twas/predixcan/logs/${tissue}.spredixcan.log"
  echo "[RUN] ${tissue}"
  "$PYTHON_BIN" "$SPREDIXCAN" \
    --model_db_path "$model" \
    --covariance "$cov" \
    --gwas_file "$GWAS" \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --zscore_column Z \
    --pvalue_column P \
    --output_file "$out" > "$log" 2>&1
done

Rscript scripts/10_high_impact_extensions/02_integrate_twas_results.R
Rscript scripts/10_high_impact_extensions/phase25_twas_coloc_spatial_status.R
