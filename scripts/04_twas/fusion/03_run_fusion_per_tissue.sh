#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <tissue> [chr]" >&2
  exit 2
fi

TISSUE="$1"
CHR="${2:-}"

FUSION_R="external/twas/fusion/software/FUSION.assoc_test.R"
SUMSTATS="data/processed/twas_input/ksd_2025_for_fusion.tsv"
WEIGHT_DIR="external/twas/fusion/weights/${TISSUE}"
REF_LD_PREFIX="external/twas/fusion/ref_ld/1000G.EUR."

OUTDIR="results/twas/fusion/per_tissue/${TISSUE}"
LOGDIR="results/twas/fusion/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$FUSION_R" ]]; then
  echo "[ERROR] FUSION script missing: $FUSION_R" >&2
  exit 1
fi

if [[ ! -s "$SUMSTATS" ]]; then
  echo "[ERROR] FUSION sumstats missing: $SUMSTATS" >&2
  exit 1
fi

if [[ ! -d "$WEIGHT_DIR" ]]; then
  echo "[ERROR] weight directory missing for ${TISSUE}: $WEIGHT_DIR" >&2
  exit 1
fi

WEIGHT_POS=$(find "$WEIGHT_DIR" -name "*.pos" | head -n 1 || true)
if [[ -z "$WEIGHT_POS" || ! -s "$WEIGHT_POS" ]]; then
  echo "[ERROR] weight pos file not found for ${TISSUE}" >&2
  exit 1
fi

if [[ ! -s "${REF_LD_PREFIX}1.bim" && ! -s "${REF_LD_PREFIX}chr1.bim" ]]; then
  echo "[ERROR] FUSION LD reference missing for prefix: $REF_LD_PREFIX" >&2
  exit 1
fi

if [[ -n "$CHR" ]]; then
  chromosomes=("$CHR")
else
  chromosomes=($(seq 1 22))
fi

for chr in "${chromosomes[@]}"; do
  out="${OUTDIR}/${TISSUE}.chr${chr}.fusion.tsv"
  log="${LOGDIR}/${TISSUE}.chr${chr}.log"

  Rscript "$FUSION_R" \
    --sumstats "$SUMSTATS" \
    --weights "$WEIGHT_POS" \
    --weights_dir "$WEIGHT_DIR" \
    --ref_ld_chr "$REF_LD_PREFIX" \
    --chr "$chr" \
    --out "$out" > "$log" 2>&1
done
