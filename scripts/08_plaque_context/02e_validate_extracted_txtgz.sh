#!/usr/bin/env bash
set -euo pipefail

IN="results/gse73680/tables/gse73680_extracted_txtgz_files.tsv"
OUT="results/gse73680/tables/gse73680_extracted_txtgz_validation.tsv"

if [[ ! -s "$IN" ]]; then
  echo "[FAILED] missing extracted TXT.gz file list"
  exit 1
fi

printf 'file\tgzip_valid\tfile_size_bytes\n' > "$OUT"

while IFS= read -r f; do
  size=$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")
  if gzip -t "$f" 2>/dev/null; then
    valid="TRUE"
  else
    valid="FALSE"
  fi
  printf '%s\t%s\t%s\n' "$f" "$valid" "$size" >> "$OUT"
done < "$IN"

n_valid=$(awk -F'\t' 'NR>1 && $2=="TRUE"{n++} END{print n+0}' "$OUT")
n_total=$(awk 'NR>1{n++} END{print n+0}' "$OUT")

echo "[OK] valid gzip: $n_valid / $n_total"
