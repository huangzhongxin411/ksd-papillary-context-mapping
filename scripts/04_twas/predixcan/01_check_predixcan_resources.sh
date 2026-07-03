#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tables

SPREDIXCAN="external/twas/predixcan/software/MetaXcan/software/SPrediXcan.py"
if [[ ! -s "$SPREDIXCAN" ]]; then
  SPREDIXCAN="external/twas/predixcan/software/SPrediXcan.py"
fi
software_status="missing"
if [[ -s "$SPREDIXCAN" ]]; then
  software_status="present"
fi

{
  printf "tissue\tmodel_file\tcovariance_file\tspredixcan_status\tstatus\tnotes\n"
  for tissue in kidney_cortex whole_blood artery_aorta artery_tibial adipose_subcutaneous liver colon_transverse small_intestine_terminal_ileum; do
    model=$(find external/twas/predixcan/models -iname "*${tissue}*.db" | head -n 1 || true)
    cov=$(find external/twas/predixcan/covariance -iname "*${tissue}*" | head -n 1 || true)
    [[ -z "$model" ]] && model="NA"
    [[ -z "$cov" ]] && cov="NA"

    status="missing"
    notes=""
    if [[ "$software_status" == "present" && "$model" != "NA" && "$cov" != "NA" ]]; then
      status="ready"
    fi
    if [[ "$software_status" != "present" ]]; then notes="${notes}SPrediXcan.py missing; "; fi
    if [[ "$model" == "NA" ]]; then notes="${notes}model db missing; "; fi
    if [[ "$cov" == "NA" ]]; then notes="${notes}covariance missing; "; fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tissue" "$model" "$cov" "$software_status" "$status" "$notes"
  done
} > results/tables/predixcan_resource_check.tsv

cat results/tables/predixcan_resource_check.tsv
