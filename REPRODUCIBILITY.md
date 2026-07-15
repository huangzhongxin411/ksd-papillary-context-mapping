# Reproducibility guide for manuscript v3.4

## Interpretation boundary

The workflow is observational and context-mapping focused. MAGMA is an EUR-LD-reference-based priority layer; snRNA provides donor-level cell context; spatial projection is supplementary broad-compartment context; Kidney_Cortex TWAS is proxy annotation; and paired bulk analysis is disease-context sensitivity. None of these layers establishes causal genes, plaque-specific localization or papilla-specific regulation. R1-R6 are reporting groups, not causal tiers.

## External inputs

Large third-party inputs are not bundled. Retrieve the KSD GWAS, GSE231569, GSE73680, GSE206306, GSE231630, GTEx/PredictDB models and required reference resources from their original distributions. Configure paths locally without committing credentials or personal absolute paths.

## Phase 1: GWAS and MAGMA

Use `scripts/phase1/` for GWAS QC, diagnostic plotting, lead/locus reconstruction and MAGMA freezing. Current outputs are Figure 1, Supplementary Figure S1 and Supplementary Tables 1-2.

## Phase 2: snRNA

Use `scripts/phase2/` for object audit, donor-level module scoring, matched-random benchmarking, low-cell sensitivity, known-driver removal and Figure 2 assembly. Current outputs are Figure 2 and Supplementary Table 3.

## Phase 3: spatial supplementary projection

Use `scripts/phase3/` for input locking, broad-compartment label transfer, usability auditing, descriptive module overlays and evidence integration. The final scope is ten complete sections (10 total) = four GSE206306 + six GSE231630. The current output is Supplementary Figure S2 only; no lesion ROI or plaque-specific localization is claimed.

## Phase 4: TWAS and candidate reporting

Use `scripts/phase4/` for Kidney_Cortex proxy TWAS auditing, candidate-model repair and Figure 3 assembly. Current outputs are Figure 3, Supplementary Figure S3 and Supplementary Tables 4-5.

## Phase 5: bulk disease context

Use `scripts/phase5/` for input audit, paired module models, tissue-state sensitivity and Figure 4 assembly. Current outputs are Figure 4 and Supplementary Table 6.

## Phase 6: manuscript and package generation

Use `scripts/phase6/` for v3.0-v3.2 manuscript assembly and supplement/package materialization. These scripts do not perform new biological analyses.

## Phase 7: v3.3 polish and synchronization

Use `scripts/phase7/` for visual polishing, targeted language polishing, clean manuscript rendering, package auditing and repository synchronization. These steps use locked inputs and do not change scientific results, numerical values or claim boundaries.

## Phase 8: v3.4 targeted polish and repository synchronization

The v3.4 package reorders approved manuscript wording and re-renders Figure 3 and Supplementary Figure S2 from locked inputs. No biological analysis, numerical result, source table, candidate assignment, MAGMA module or R1-R6 membership was changed. The spatial scope remains ten sections: four GSE206306 and six GSE231630.

## Exact source-data and script mapping

Because no analysis was rerun, use `scripts/SCRIPT_MAP_v3.3.tsv` to map each locked analysis output to its generating or packaging script. Use `source_data/Source_Data_manifest_v3.3.tsv` for panel-level source files. `MANIFEST.tsv` lists the repository payload and `CHECKSUMS.sha256` verifies byte identity.

## Clean-checkout verification

From the repository root, confirm that every `MANIFEST.tsv` path exists and run `shasum -a 256 -c CHECKSUMS.sha256` on macOS. On Linux, use `sha256sum -c CHECKSUMS.sha256`. Raw-input-dependent scripts require locally configured public inputs and are not expected to run without those external resources.
