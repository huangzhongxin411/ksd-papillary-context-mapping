# Reproducibility guide for manuscript v3.2

## Interpretation boundary

The workflow is observational and context-mapping focused. MAGMA is an EUR-LD-reference-based priority layer; snRNA provides donor-level cell context; spatial projection is supplementary broad-compartment context; Kidney_Cortex TWAS is proxy annotation; and paired bulk analysis is disease-context sensitivity. None of these layers establishes causal genes, plaque-specific localization or papilla-specific regulation.

## External inputs

Large third-party inputs are not bundled. Retrieve the KSD GWAS, GSE231569, GSE73680, GSE206306, GSE231630, GTEx/PredictDB models and the required reference resources from their original distributions. Configure paths locally without committing credentials or personal absolute paths.

## Phase 1: GWAS and MAGMA

- GWAS QC and diagnostic plotting: `scripts/phase1/phase1_step4_gwas_qc_audit.py`, `scripts/phase1/phase1_step4_gwas_qc_plots.R`.
- Lead/locus reconstruction: `scripts/phase1/phase1_step4_reconstruct_leads_loci.py`.
- MAGMA gene-set freeze and provenance audit: `scripts/phase1/phase1_step3_freeze_magma_gene_sets.py`, `scripts/phase1/phase1_step5_magma_provenance_audit.py`.
- Final outputs: Figure 1, Supplementary Figure S1, Supplementary Tables 1-2.

## Phase 2: snRNA

- Input/object audit, scoring and sensitivity: `scripts/phase2/phase2_step1_scrna_object_audit.R` through `phase2_step4_known_driver_removal_sensitivity.R`.
- Final Figure 2 assembly: `scripts/phase2/phase2_step5C_patch_figure2.R`.
- Final outputs: Figure 2 and Supplementary Table 3.

## Phase 3: spatial supplementary projection

- Input/metadata locking, label transfer, usability audit, module overlays and integration: `scripts/phase3/phase3_step1_spatial_input_audit.R` through `phase3_step4_spatial_evidence_integration.R`.
- Final data scope: ten complete sections (10 total) = four GSE206306 + six GSE231630.
- Final output: Supplementary Figure S2 only. Superseded spatial/Figure 5 workflows are not current.

## Phase 4: TWAS and candidate reporting

- TWAS proxy audit: `scripts/phase4/phase4_step1_twas_audit.R`.
- Candidate-model repair and figure assembly: `scripts/phase4/phase4_step2_candidate_evidence_model_repair.R`, `scripts/phase4/phase6_step4_reexport_figure3_bulk_reviewed.R`.
- Final outputs: Figure 3, Supplementary Figure S3, Supplementary Tables 4-5.

## Phase 5: bulk disease context

- Input audit, paired module models and tissue-state sensitivity: `scripts/phase5/phase5_step1_bulk_input_audit.R` through `phase5_step3_bulk_tissue_state_sensitivity.R`.
- Final Figure 4 assembly: `scripts/phase5/phase5_step4_assemble_figure4_bulk_context.R`.
- Final outputs: Figure 4 and Supplementary Table 6.

## Phase 6: manuscript and package generation

- v3.0/v3.1/v3.2 assembly and final supplement/package materialization: `scripts/phase6/phase6_step2_build_manuscript_v3.py`, `phase6_step4_build_manuscript_v3_1.py`, `phase6_step6_build_manuscript_v3_2.py`, `phase6_step6_materialize_supplement.py`, `phase6_step6_finalize_submission_package.py`.
- No new biological analyses are performed by these packaging scripts.

## Exact source-data and script mapping

Use `scripts/SCRIPT_MAP_v3.2.tsv` to map each final output to its generating or packaging script. Use `source_data/Source_Data_manifest_v3.2.tsv` for panel-level source files. `MANIFEST.tsv` lists every allowlisted release file and `CHECKSUMS.sha256` verifies byte identity.

## Clean-checkout verification

From the repository root, confirm that every `MANIFEST.tsv` path exists and validate hashes with `shasum -a 256 -c CHECKSUMS.sha256` on macOS or `sha256sum -c CHECKSUMS.sha256` on Linux after converting the command format if needed. Raw-input-dependent scripts require locally configured public inputs and are not expected to run without those external resources.
