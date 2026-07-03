# Environment status

- MAGMA: v1.10; NCBI build 37 gene locations; 1000 Genomes European LD reference.
- R: used across single-nucleus, bulk, spatial and figure stages; a single final cross-stage `sessionInfo()` is not available.
- Python: used for auditing, packaging and selected TWAS utilities; a fully frozen cross-stage lock file is not available.
- S-PrediXcan: GTEx v8 Kidney_Cortex MASHR models; external models are not redistributed.

Before a full rerun, record `R --version`, `sessionInfo()`, `python --version`, installed package versions, MAGMA executable checksum and the LD-reference prefix/checksums.
