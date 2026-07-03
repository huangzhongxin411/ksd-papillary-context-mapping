#!/usr/bin/env bash
set -euo pipefail

gwas="${1:-data/raw/gwas/2025_trans_ancestry/meta_sumstats}"
expected_md5='83f1644eaf103c8a3e74f78be9c48b16'

if [[ ! -s "$gwas" ]]; then
  printf 'missing_or_empty\t%s\n' "$gwas" >&2
  printf 'Download from: https://zenodo.org/records/14790324/files/meta_sumstats?download=1\n' >&2
  exit 1
fi

printf 'path\t%s\n' "$gwas"
ls -lh "$gwas"
file "$gwas"

if command -v md5sum >/dev/null 2>&1; then
  actual="$(md5sum "$gwas" | awk '{print $1}')"
elif command -v md5 >/dev/null 2>&1; then
  actual="$(md5 -q "$gwas")"
else
  printf 'md5_tool\tNOT_FOUND\n'
  actual=''
fi

if [[ -n "$actual" ]]; then
  printf 'md5\t%s\n' "$actual"
  if [[ "$actual" != "$expected_md5" ]]; then
    printf 'md5_status\tMISMATCH expected=%s\n' "$expected_md5" >&2
    exit 2
  fi
  printf 'md5_status\tOK\n'
fi

printf 'format_preview\n'
case "$(file -b "$gwas")" in
  *gzip*)
    gzip -dc "$gwas" | head -n 5
    ;;
  *Zip*)
    unzip -l "$gwas" | head
    ;;
  *tar*)
    tar -tf "$gwas" | head
    ;;
  *)
    head -n 5 "$gwas"
    ;;
esac

