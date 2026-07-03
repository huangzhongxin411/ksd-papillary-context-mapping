#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/GSE206306 data/raw/GSE231630 results/tables

cat > results/tables/spatial_download_status.tsv <<'EOF'
dataset	status	notes
GSE206306	not_downloaded	Spatial resource download is intentionally not automatic until exact GEO supplementary file names are confirmed.
GSE231630	not_downloaded	Spatial resource download is intentionally not automatic until exact GEO supplementary file names are confirmed.
EOF

printf 'Wrote results/tables/spatial_download_status.tsv\n'
