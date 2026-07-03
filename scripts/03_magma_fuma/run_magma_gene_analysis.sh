#!/usr/bin/env bash
set -euo pipefail

MAGMA_BIN="${MAGMA_BIN:-external/magma/magma}"
GENE_LOC="${GENE_LOC:-external/reference/gene_loc/NCBI37.3.gene.loc}"
LD_PREFIX="${LD_PREFIX:-external/reference/1000G_EUR/g1000_eur}"
SNPLOC="${SNPLOC:-data/processed/magma_input/ksd_2025.dedup.snploc}"
PVAL="${PVAL:-data/processed/magma_input/ksd_2025.dedup.pval}"
OUT_PREFIX="${OUT_PREFIX:-results/magma/2025_trans_ancestry/ksd_2025}"

mkdir -p "$(dirname "$OUT_PREFIX")" results/logs

if [[ ! -x "$MAGMA_BIN" ]]; then
  echo "MAGMA executable not found or not executable: $MAGMA_BIN" >&2
  echo "Place the MAGMA binary at external/magma/magma, then rerun this script." >&2
  exit 2
fi
if [[ ! -s "$GENE_LOC" ]]; then
  echo "Gene location file missing: $GENE_LOC" >&2
  exit 3
fi
if [[ ! -s "${LD_PREFIX}.bed" || ! -s "${LD_PREFIX}.bim" || ! -s "${LD_PREFIX}.fam" ]]; then
  echo "LD reference missing: ${LD_PREFIX}.{bed,bim,fam}" >&2
  exit 4
fi

"$MAGMA_BIN" --version > results/logs/magma_version.txt

"$MAGMA_BIN" \
  --annotate \
  --snp-loc "$SNPLOC" \
  --gene-loc "$GENE_LOC" \
  --out "$OUT_PREFIX"

"$MAGMA_BIN" \
  --bfile "$LD_PREFIX" \
  --pval "$PVAL" use=SNP,P ncol=N \
  --gene-annot "${OUT_PREFIX}.genes.annot" \
  --out "$OUT_PREFIX"

printf 'wrote\t%s\n' "${OUT_PREFIX}.genes.out"
printf 'wrote\t%s\n' "${OUT_PREFIX}.genes.raw"
