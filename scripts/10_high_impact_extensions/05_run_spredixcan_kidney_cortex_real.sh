#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-external/twas/predixcan/metaxcan_venv/bin/python3}"
SPREDIXCAN="${SPREDIXCAN:-external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py}"
PREDICTDB_DIR="${PREDICTDB_DIR:-external/twas/predixcan/predictdb_gtex_v8_mashr}"
GWAS="${GWAS:-data/processed/twas_input/ksd_2025_for_twas.tsv.gz}"
MAPPED_GWAS="${MAPPED_GWAS:-data/processed/twas_input/ksd_2025_for_twas.Kidney_Cortex_varID.tsv.gz}"
MAPPED_GWAS_QC="${MAPPED_GWAS_QC:-results/twas/ksd_2025_for_twas.Kidney_Cortex_varID_qc.tsv}"
OUT="${OUT:-results/twas/spredixcan_Kidney_Cortex_KSD.csv}"
LOG="${LOG:-logs/twas/spredixcan_Kidney_Cortex_KSD.log}"
SNP_MAP="${SNP_MAP:-external/twas/predixcan/predictdb_gtex_v8_mashr/eqtl/mashr/mashr_Kidney_Cortex.rsid_to_varID.tsv}"

mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")" "$(dirname "$MAPPED_GWAS")"

MODEL=$(find "$PREDICTDB_DIR" -iname "*Kidney*Cortex*.db" -o -iname "mashr_Kidney_Cortex.db" | head -n 1 || true)
COV=$(find "$PREDICTDB_DIR" -iname "*Kidney*Cortex*.txt.gz" -o -iname "mashr_Kidney_Cortex.txt.gz" | head -n 1 || true)

if [[ ! -s "$PYTHON_BIN" ]]; then
  echo "[ERROR] Python environment missing: $PYTHON_BIN" >&2
  exit 1
fi
if [[ ! -s "$SPREDIXCAN" ]]; then
  echo "[ERROR] SPrediXcan missing: $SPREDIXCAN" >&2
  exit 1
fi
if [[ ! -s "$MODEL" ]]; then
  echo "[ERROR] Kidney_Cortex model db missing under $PREDICTDB_DIR" >&2
  exit 1
fi
if [[ ! -s "$COV" ]]; then
  echo "[ERROR] Kidney_Cortex covariance missing under $PREDICTDB_DIR" >&2
  exit 1
fi
if [[ ! -s "$GWAS" ]]; then
  echo "[ERROR] GWAS input missing: $GWAS" >&2
  exit 1
fi

"$PYTHON_BIN" - <<'PY'
import numpy, scipy, pandas, h5py
print("python_dependencies_ok")
PY

"$PYTHON_BIN" - "$GWAS" > logs/twas/ksd_2025_for_twas_header.txt <<'PY'
import gzip
import sys
with gzip.open(sys.argv[1], "rt") as handle:
    print(handle.readline().rstrip("\n"))
PY

if [[ ! -s "$SNP_MAP" ]]; then
  "$PYTHON_BIN" scripts/10_high_impact_extensions/07_make_predixcan_snp_map.py "$MODEL" "$SNP_MAP"
fi

"$PYTHON_BIN" scripts/10_high_impact_extensions/08_make_kidney_cortex_varid_gwas.py \
  "$MODEL" "$GWAS" "$MAPPED_GWAS" "$MAPPED_GWAS_QC"

"$PYTHON_BIN" "$SPREDIXCAN" \
  --model_db_path "$MODEL" \
  --model_db_snp_key varID \
  --covariance "$COV" \
  --gwas_file "$MAPPED_GWAS" \
  --snp_column SNP \
  --effect_allele_column A1 \
  --non_effect_allele_column A2 \
  --beta_column BETA \
  --se_column SE \
  --pvalue_column P \
  --keep_non_rsid \
  --gwas_N 796644 \
  --output_file "$OUT" \
  --additional_output \
  --verbosity 7 > "$LOG" 2>&1

Rscript scripts/10_high_impact_extensions/06_summarize_spredixcan_kidney_cortex_real.R
Rscript scripts/10_high_impact_extensions/03_phase25b_resource_landing_status.R
