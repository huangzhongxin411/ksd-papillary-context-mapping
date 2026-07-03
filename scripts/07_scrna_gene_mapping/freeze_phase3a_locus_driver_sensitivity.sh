#!/usr/bin/env bash
set -euo pipefail

freeze_dir="results/phase3a_locus_driver_sensitivity_v0.1"
mkdir -p "$freeze_dir/tables" "$freeze_dir/figures" "$freeze_dir/gene_sets" "$freeze_dir/docs" "$freeze_dir/scripts"

cp results/tables/audited_locus_top50_tal_driver_genes.tsv "$freeze_dir/tables/"
cp results/tables/locus_balanced_top50_genes.tsv "$freeze_dir/tables/"
cp results/tables/audited_locus_balanced_scrna_benchmark.tsv "$freeze_dir/tables/"
cp results/tables/audited_locus_leave_one_locus_out.tsv "$freeze_dir/tables/"
cp results/tables/umod_locus_context_summary.tsv "$freeze_dir/tables/"
cp results/gene_sets/locus_balanced_top50.txt "$freeze_dir/gene_sets/"
cp results/figures/audited_locus_leave_one_locus_out_tal.pdf "$freeze_dir/figures/"
cp docs/phase3a_locus_driver_sensitivity_notes.md "$freeze_dir/docs/"
cp docs/phase3a_freeze_decision.md "$freeze_dir/docs/"
cp scripts/07_scrna_gene_mapping/audited_locus_driver_balanced_sensitivity.R "$freeze_dir/scripts/"
cp scripts/07_scrna_gene_mapping/make_umod_locus_context_summary.R "$freeze_dir/scripts/"

{
  printf 'file\tsize_bytes\tmd5\n'
  find "$freeze_dir" -type f ! -name 'MANIFEST.tsv' | sort | while read -r f; do
    size=$(wc -c < "$f" | tr -d ' ')
    md5=$(md5 -q "$f")
    printf '%s\t%s\t%s\n' "$f" "$size" "$md5"
  done
} > "$freeze_dir/MANIFEST.tsv"

cat > "$freeze_dir/README.md" <<'EOF'
# Phase 3A Locus Driver Sensitivity v0.1

This frozen bundle contains the exploratory locus-based single-cell robustness checks for the audited GSE231569 projection.

Interpretation boundary:

Top-ranked locus-mapped KSD genes suggest a TAL-associated localization pattern, partly driven by influential loci including the UMOD/HIBADH/FAM13A-related loci. Locus-balanced and leave-one-locus-out analyses do not support broad TAL enrichment. MAGMA/TWAS confirmation is required before using TAL as a formal manuscript claim.
EOF

printf 'wrote\t%s\n' "$freeze_dir/MANIFEST.tsv"
printf 'wrote\t%s\n' "$freeze_dir/README.md"
