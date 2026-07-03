#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tables

{
  printf 'dataset\tlocal_root\tn_files\thas_matrix\thas_image\thas_coordinates\tstatus\n'
  for dataset in GSE206306 GSE231630; do
    root="data/raw/${dataset}"
    if [ ! -d "${root}" ]; then
      printf '%s\t%s\t0\tFALSE\tFALSE\tFALSE\tmissing_directory\n' "${dataset}" "${root}"
      continue
    fi
    n_files=$(find "${root}" -type f | wc -l | tr -d ' ')
    has_matrix=FALSE
    has_image=FALSE
    has_coordinates=FALSE
    if find "${root}" -type f \( -name '*.h5' -o -name 'matrix.mtx*' -o -name '*count*matrix*' -o -name '*filtered*feature*' \) | grep -q .; then
      has_matrix=TRUE
    fi
    if find "${root}" -type f \( -name '*.tif' -o -name '*.tiff' -o -name '*.png' -o -name '*.jpg' -o -name '*.jpeg' \) | grep -q .; then
      has_image=TRUE
    fi
    if find "${root}" -type f \( -name '*tissue_positions*' -o -name '*spatial*coordinates*' -o -name '*positions*csv' \) | grep -q .; then
      has_coordinates=TRUE
    fi
    status=missing_required_files
    if [ "${has_matrix}" = TRUE ] && [ "${has_image}" = TRUE ] && [ "${has_coordinates}" = TRUE ]; then
      status=ready_for_inventory
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${dataset}" "${root}" "${n_files}" "${has_matrix}" "${has_image}" "${has_coordinates}" "${status}"
  done
} > results/tables/spatial_file_detection.tsv

printf 'Wrote results/tables/spatial_file_detection.tsv\n'
