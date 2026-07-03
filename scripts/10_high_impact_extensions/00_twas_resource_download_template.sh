#!/usr/bin/env bash
set -euo pipefail

mkdir -p external/twas/predixcan/models
mkdir -p external/twas/predixcan/covariance
mkdir -p external/twas/fusion/software
mkdir -p external/twas/fusion/weights
mkdir -p external/twas/fusion/ref_ld

cat <<'MSG'
Manual resource targets:

1. PredictDB GTEx v8 model .db files
   Put files under:
   external/twas/predixcan/models/

2. PredictDB GTEx v8 covariance files
   Put files under:
   external/twas/predixcan/covariance/

3. Optional FUSION resources
   FUSION.assoc_test.R -> external/twas/fusion/software/
   GTEx weights       -> external/twas/fusion/weights/
   1000G EUR LD       -> external/twas/fusion/ref_ld/

This template intentionally does not auto-download large resources because
PredictDB covariance/model archives and FUSION weights can be multi-GB.
After placing resources, run:

Rscript scripts/10_high_impact_extensions/phase25_twas_coloc_spatial_status.R
bash scripts/10_high_impact_extensions/01_run_spredixcan_priority_tissues.sh
MSG
