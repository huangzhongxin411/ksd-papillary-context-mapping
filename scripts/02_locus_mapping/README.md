# 02 Locus Mapping

## Goal

Define lead SNPs, independent loci, and mapped genes from cleaned KSD GWAS summary statistics.

## Planned Outputs

- `results/tables/gwas_lead_snps.tsv`
- `results/tables/gwas_independent_loci.tsv`
- `results/tables/locus_gene_mapping.tsv`

## Phase 1 Distance-Based Fallback

If LD reference is not ready, use distance-based pruning only for the minimal loop:

```bash
python3 scripts/02_locus_mapping/make_phase1_leads_loci.py \
  --input data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz \
  --p 5e-8 \
  --window 1000000 \
  --out-prefix results/tables/phase1_2025
```

These outputs must be described as exploratory `distance_pruned` results, not formal LD-clumped loci.

Then map genes using GENCODE:

```bash
python3 scripts/02_locus_mapping/map_phase1_candidate_genes.py \
  --leads results/tables/phase1_2025_lead_snps.tsv \
  --loci results/tables/phase1_2025_loci.tsv \
  --gtf data/references/gencode/gencode.v44.annotation.gtf.gz \
  --out results/tables/phase1_candidate_genes.tsv
```

## Mapping Evidence

- Positional mapping within predefined windows.
- eQTL mapping where suitable references exist.
- Chromatin or enhancer-promoter mapping if using FUMA or external annotations.
