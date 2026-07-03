#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tables

FUSION_R="external/twas/fusion/software/FUSION.assoc_test.R"
REF_LD_PREFIX="external/twas/fusion/ref_ld/1000G.EUR."

software_status="missing"
if [[ -s "$FUSION_R" ]]; then
  software_status="present"
fi

ld_status="missing"
if [[ -s "${REF_LD_PREFIX}1.bim" || -s "${REF_LD_PREFIX}chr1.bim" ]]; then
  ld_status="present"
fi

{
  printf "tissue\tweight_dir\tweight_pos_file\tn_weight_files\tfusion_software_status\tld_reference_status\tstatus\tnotes\n"
  for tissue in kidney_cortex whole_blood artery_aorta artery_tibial adipose_subcutaneous liver colon_transverse small_intestine_terminal_ileum; do
    dir="external/twas/fusion/weights/${tissue}"
    pos="NA"
    n_weights=0
    status="missing"
    notes=""

    if [[ -d "$dir" ]]; then
      pos=$(find "$dir" -name "*.pos" | head -n 1 || true)
      [[ -z "$pos" ]] && pos="NA"
      n_weights=$(find "$dir" -name "*.RDat" | wc -l | tr -d ' ')
    fi

    if [[ "$software_status" == "present" && "$ld_status" == "present" && "$pos" != "NA" && "$n_weights" -gt 0 ]]; then
      status="ready"
    fi

    if [[ "$software_status" != "present" ]]; then notes="${notes}FUSION.assoc_test.R missing; "; fi
    if [[ "$ld_status" != "present" ]]; then notes="${notes}FUSION LD reference missing; "; fi
    if [[ "$pos" == "NA" ]]; then notes="${notes}weights.pos missing; "; fi
    if [[ "$n_weights" -eq 0 ]]; then notes="${notes}RDat weights missing; "; fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tissue" "$dir" "$pos" "$n_weights" "$software_status" "$ld_status" "$status" "$notes"
  done
} > results/tables/fusion_resource_check.tsv

cat results/tables/fusion_resource_check.tsv
