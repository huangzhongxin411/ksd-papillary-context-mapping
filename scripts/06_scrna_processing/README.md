# 06 snRNA-seq Processing

## Dataset

Primary dataset: GSE231569 human kidney papilla single-nucleus atlas.

## Key Principles

- Treat nuclei, not cells; mitochondrial thresholds should be data-driven and not copied mechanically from scRNA-seq.
- Avoid cell-level pseudoreplication for inferential statistics.
- Prefer donor x cell type pseudobulk for statistical comparisons.

## Broad Cell Types

- Proximal tubule: `LRP2`, `SLC34A1`, `CUBN`
- Loop of Henle / thin limb: `UMOD`, `SLC12A1`, `CLDN10`
- Collecting duct principal cell: `AQP2`, `AQP3`, `SCNN1G`
- Intercalated cell: `ATP6V1B1`, `SLC4A1`, `FOXI1`
- Undifferentiated epithelial: `KRT8`, `KRT18`, `PROM1`
- Endothelial: `PECAM1`, `VWF`, `KDR`
- Fibroblast / stromal: `DCN`, `LUM`, `COL1A1`
- Macrophage / myeloid: `LST1`, `C1QA`, `CD68`
- T cell: `CD3D`, `TRAC`
- Smooth muscle / pericyte: `RGS5`, `ACTA2`

## Output

- Clean Seurat or AnnData object.
- Cell-type annotation table.
- Donor x cell type pseudobulk matrix.
- Marker gene table.

