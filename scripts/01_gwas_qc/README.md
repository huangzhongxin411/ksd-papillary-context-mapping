# 01 GWAS QC

## Input

KSD GWAS summary statistics，至少包含：

- SNP or rsID
- chromosome
- position
- effect allele
- non-effect allele
- beta or OR
- standard error
- P value
- effect allele frequency
- sample size
- imputation INFO, if available

## QC

- Harmonize genome build, preferably GRCh37/hg19 for compatibility with many TWAS weights and LD references.
- Remove INFO < 0.8 when INFO is available.
- Remove MAF < 0.01.
- Remove ambiguous A/T and C/G SNPs unless allele frequency can resolve strand.
- Remove duplicated variants.
- Remove variants without beta/SE/P.

## Output

- `data/processed/gwas/<dataset>.clean.tsv.gz`
- `results/tables/<dataset>.snp_qc_report.tsv`
- QQ plot
- Manhattan plot
- genomic inflation factor and LDSC intercept if available

