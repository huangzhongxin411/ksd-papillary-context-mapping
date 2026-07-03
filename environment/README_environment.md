# Environment documentation

`R_sessionInfo.txt` and `python_requirements.txt` capture the environment available during the Stage 7J-R2 release-preparation pass. They improve transparency but do not reconstruct every package version used across all historical analysis stages.

The project used R and Python across multiple stages, MAGMA v1.10 with NCBI build 37 gene locations and the 1000 Genomes European LD reference, and S-PrediXcan with GTEx v8 Kidney_Cortex MASHR models. External LD and prediction-model resources are not redistributed.

A complete cross-stage `renv.lock`, Conda environment or original session snapshot was not available. Reusers should treat the supplied files as release-time references, inspect individual scripts for package calls, and record fresh environment details for any full rerun.
