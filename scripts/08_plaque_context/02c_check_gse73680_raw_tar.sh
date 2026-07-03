#!/usr/bin/env bash
set -euo pipefail

TAR="data/raw/gse73680/GSE73680_RAW.tar"
OUTDIR="results/gse73680/tables"
LOGDIR="results/gse73680/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$TAR" ]]; then
  echo "[FAILED] missing tar"
  exit 1
fi

tar -tf "$TAR" > "$OUTDIR/gse73680_raw_tar_filelist.tsv" 2> "$LOGDIR/gse73680_raw_tar_check.log"

n_files=$(wc -l < "$OUTDIR/gse73680_raw_tar_filelist.tsv" | tr -d ' ')
n_txtgz=$(grep -Ei '\.txt\.gz$|\.TXT\.gz$|\.txt\.GZ$|\.TXT\.GZ$' "$OUTDIR/gse73680_raw_tar_filelist.tsv" | wc -l | tr -d ' ')
size=$(stat -f%z "$TAR" 2>/dev/null || stat -c%s "$TAR")

{
  printf 'check\tvalue\n'
  printf 'tar_exists\tTRUE\n'
  printf 'tar_size_bytes\t%s\n' "$size"
  printf 'n_files_in_tar\t%s\n' "$n_files"
  printf 'n_txtgz_in_tar\t%s\n' "$n_txtgz"
} > "$OUTDIR/gse73680_raw_tar_qc.tsv"

if [[ "$n_txtgz" -lt 1 ]]; then
  echo "[FAILED] no TXT.gz detected in tar"
  exit 1
fi

echo "[OK] RAW tar is readable"
