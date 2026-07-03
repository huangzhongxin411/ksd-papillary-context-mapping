#!/usr/bin/env bash
set -euo pipefail

TARGET_DIR="external/twas/predixcan/predictdb_gtex_v8_mashr"
TAR_FILE="${TARGET_DIR}/mashr_eqtl.tar"
URL="https://zenodo.org/records/3518299/files/mashr_eqtl.tar?download=1"

mkdir -p "$TARGET_DIR"

if [[ ! -s "$TAR_FILE" ]]; then
  echo "[INFO] downloading mashr_eqtl.tar from Zenodo"
  curl -L --fail --show-error -C - "$URL" -o "$TAR_FILE"
else
  echo "[INFO] found existing $TAR_FILE"
fi

echo "[INFO] extracting $TAR_FILE"
tar -xvf "$TAR_FILE" -C "$TARGET_DIR"

echo "[INFO] Kidney_Cortex candidates:"
find "$TARGET_DIR" -iname "*Kidney*Cortex*" -print
