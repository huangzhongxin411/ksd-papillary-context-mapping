#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/gse73680 results/gse73680/tables results/gse73680/logs

RAW_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73680/suppl/GSE73680_RAW.tar"
OUT="data/raw/gse73680/GSE73680_RAW.tar"
LOG="results/gse73680/logs/download_gse73680_raw_tar.log"
STATUS="results/gse73680/tables/gse73680_raw_tar_status.tsv"

{
  printf 'field\tvalue\n'
  printf 'raw_url\t%s\n' "$RAW_URL"
  printf 'local_path\t%s\n' "$OUT"
  printf 'download_started\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} > "$STATUS"

curl -L -C - --retry 5 --retry-delay 20 "$RAW_URL" -o "$OUT" 2>&1 | tee "$LOG"

if [[ ! -s "$OUT" ]]; then
  printf 'download_status\tfailed_missing_or_empty\n' >> "$STATUS"
  echo "[FAILED] RAW tar missing or empty: $OUT"
  exit 1
fi

size=$(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT")
{
  printf 'download_finished\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf 'download_status\tdownloaded_or_resumed\n'
  printf 'tar_size_bytes\t%s\n' "$size"
} >> "$STATUS"

ls -lh "$OUT" | tee -a "$LOG"
echo "[OK] downloaded RAW tar"
