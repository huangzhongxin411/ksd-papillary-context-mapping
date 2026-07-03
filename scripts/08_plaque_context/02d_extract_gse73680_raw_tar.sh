#!/usr/bin/env bash
set -euo pipefail

TAR="data/raw/gse73680/GSE73680_RAW.tar"
EXTRACT_DIR="data/raw/gse73680/extracted"
OUTDIR="results/gse73680/tables"
LOGDIR="results/gse73680/logs"

mkdir -p "$EXTRACT_DIR" "$OUTDIR" "$LOGDIR"

if [[ ! -s "$TAR" ]]; then
  echo "[FAILED] missing tar"
  exit 1
fi

tar -xvf "$TAR" -C "$EXTRACT_DIR" > "$LOGDIR/gse73680_raw_tar_extract.log" 2>&1

find "$EXTRACT_DIR" -type f | sort > "$OUTDIR/gse73680_extracted_filelist.tsv"
find "$EXTRACT_DIR" -type f | grep -Ei '\.txt\.gz$|\.TXT\.gz$' | sort > "$OUTDIR/gse73680_extracted_txtgz_files.tsv" || true

n_txtgz=$(wc -l < "$OUTDIR/gse73680_extracted_txtgz_files.tsv" | tr -d ' ')

{
  printf 'metric\tvalue\n'
  printf 'extract_dir\t%s\n' "$EXTRACT_DIR"
  printf 'n_txtgz_extracted\t%s\n' "$n_txtgz"
} > "$OUTDIR/gse73680_extract_summary.tsv"

if [[ "$n_txtgz" -lt 1 ]]; then
  echo "[FAILED] no TXT.gz extracted"
  exit 1
fi

echo "[OK] extracted $n_txtgz TXT.gz files"
