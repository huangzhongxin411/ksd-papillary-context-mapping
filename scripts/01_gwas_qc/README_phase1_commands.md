# Phase 1 GWAS QC Commands

## 1. Inspect columns

```bash
scripts/01_gwas_qc/verify_phase1_gwas_file.sh
scripts/01_gwas_qc/run_phase1_column_inspection.sh
```

Use the printed `suggested_mapping` to choose arguments for the QC script.

## 2. Run QC

Example:

```bash
python3 scripts/01_gwas_qc/qc_gwas_sumstats.py \
  --input data/raw/gwas/<file> \
  --out data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz \
  --report results/tables/phase1_gwas_qc_report.tsv \
  --prefix results/figures/phase1_gwas_2025 \
  --snp SNP \
  --chr CHR \
  --bp BP \
  --ea EA \
  --nea NEA \
  --p P \
  --beta BETA \
  --se SE \
  --eaf EAF \
  --n N \
  --info INFO
```

Expected outputs:

- `data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz`
- `results/tables/phase1_gwas_qc_report.tsv`
- `results/figures/phase1_gwas_2025.qq_plot.pdf`
- `results/figures/phase1_gwas_2025.manhattan_plot.pdf`
