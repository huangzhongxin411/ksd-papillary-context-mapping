#!/usr/bin/env python3
"""Build evidence-locked Stage 7A blueprint tables from frozen Stage 2-6 outputs."""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/tables/revision/stage7A_manuscript_integration_blueprint"
LOG = ROOT / "logs/revision/stage7A_manuscript_integration_blueprint/stage7A_build_blueprint_tables.log"
MANUSCRIPT = ROOT / "docs/revision/stage1/manuscript_v1.4_clean_working.md"

REQUIRED_INPUTS = [
    "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv",
    "results/tables/revision/stage2_genetic/magma_output_audit.tsv",
    "results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv",
    "results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv",
    "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
    "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
    "results/tables/revision/stage4B2_scrna_robustness/loop_tal_claim_decision_table.tsv",
    "results/tables/revision/stage5C1_gse73680_figure4_draft/gse73680_final_claim_lock_stage5C1.tsv",
    "results/tables/revision/stage6C_spatial_twas_figure5_draft/spatial_twas_final_claim_lock_stage6C.tsv",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def evidence_rows() -> list[dict[str, str]]:
    return [
        {
            "evidence_layer": "GWAS_MAGMA_prioritization",
            "main_result": "GWAS QC retained 4,915,033 of 5,960,489 rows and reconstructed 57 distance-defined loci; EUR-LD-reference MAGMA tested 17,316 genes and identified 94 Bonferroni, 369 FDR and 187 suggestive genes.",
            "support_strength": "primary_prioritization",
            "allowed_claim": "Auditable EUR-LD-reference-based gene prioritization for downstream renal papillary context mapping.",
            "disallowed_claim": "Ancestry-generalizable fine mapping; causal gene identification; trans-ancestry precise gene mapping.",
            "main_or_supplementary": "main",
            "figure": "Figure 1; Figure 3",
            "manuscript_section": "Results 1; Methods: GWAS/MAGMA",
            "notes": "1000 Genomes EUR LD is a hard interpretation boundary.",
        },
        {
            "evidence_layer": "TWAS_Kidney_Cortex_proxy",
            "main_result": "GTEx v8 Kidney_Cortex S-PrediXcan tested 5,989 genes; 51 were FDR-supported, including 42 one-SNP and 9 multi-SNP models, with 47 overlapping MAGMA support.",
            "support_strength": "supplementary_proxy_only",
            "allowed_claim": "Kidney_Cortex genetically regulated expression proxy support for a subset of MAGMA-prioritized genes.",
            "disallowed_claim": "Papilla-specific TWAS regulation; causal expression; colocalization; SMR support; evidence upgrade based on one-SNP models.",
            "main_or_supplementary": "supplementary",
            "figure": "Figure 3; Figure 5E",
            "manuscript_section": "Results 2 and 5; Methods: TWAS audit",
            "notes": "The one-SNP majority and tissue mismatch prevent use as a main validation layer.",
        },
        {
            "evidence_layer": "candidate_gene_evidence_model",
            "main_result": "A two-axis model separates MAGMA genetic priority from Kidney_Cortex proxy TWAS status across mutually exclusive R1-R6 reporting groups.",
            "support_strength": "primary_prioritization",
            "allowed_claim": "Transparent reporting groups distinguish the genetic prioritization axis from supplementary TWAS proxy evidence.",
            "disallowed_claim": "Causal tiers; high-confidence causal genes; validated disease genes; ranking by biological plausibility alone.",
            "main_or_supplementary": "main",
            "figure": "Figure 3",
            "manuscript_section": "Results 2; Methods: evidence model",
            "notes": "R1=79, R2=1, R3=13, R4=33, R5=4, R6=0 in the frozen model.",
        },
        {
            "evidence_layer": "curated_exemplar_panel",
            "main_result": "UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 provide a biology-informed role spectrum independent of reporting-group assignment.",
            "support_strength": "interpretive_only",
            "allowed_claim": "Curated exemplars illustrate renal transport, calcium/ion-handling and epithelial biological context.",
            "disallowed_claim": "Evidence upgrade; validated disease genes; therapeutic targets; uniform Loop/TAL markers.",
            "main_or_supplementary": "supplementary_interpretive",
            "figure": "Figure 3C-D",
            "manuscript_section": "Results 2; Discussion 2",
            "notes": "The exemplar flag is orthogonal to the two evidence axes.",
        },
        {
            "evidence_layer": "snRNA_LoopTAL_context",
            "main_result": "Across 43,878 nuclei, including 540 Loop/TAL nuclei from four donors, primary MAGMA modules showed donor-level Loop/TAL-associated patterns with partial matched-random support and robustness to driver removal.",
            "support_strength": "moderate_context_support",
            "allowed_claim": "Moderate donor-level Loop/TAL-associated single-nucleus context for MAGMA-prioritized modules.",
            "disallowed_claim": "Strong enrichment; causal cell type; causal mediation; plaque nucleation site; cell-level inferential support.",
            "main_or_supplementary": "main",
            "figure": "Figure 2",
            "manuscript_section": "Results 3; Methods: snRNA",
            "notes": "Matched-random evidence is partial, not uniformly positive.",
        },
        {
            "evidence_layer": "GSE73680_bulk_disease_context",
            "main_result": "Among 55 samples from 29 patients, including 26 paired patients, primary-module direction was retained overall but attenuated in composition/injury/remodeling-aware sensitivity models.",
            "support_strength": "attenuated_bulk_context_support",
            "allowed_claim": "Injury/remodeling-associated paired bulk disease-context support for MAGMA-prioritized modules.",
            "disallowed_claim": "Independent validation; composition-independent effect; cell-type-specific response; genetic causality; plaque causation; mediation.",
            "main_or_supplementary": "main_context_layer",
            "figure": "Figure 4",
            "manuscript_section": "Results 4; Methods: GSE73680",
            "notes": "Marker signatures are expression proxies, not validated cell fractions.",
        },
        {
            "evidence_layer": "spatial_LoopTAL_tissue_context",
            "main_result": "Five GSE206306 sections comprising 7,747 spots provided moderate supplementary Loop/TAL co-distribution; injury/ECM was not retained and mineralization was mixed after complexity adjustment.",
            "support_strength": "moderate_supplementary_context",
            "allowed_claim": "Moderate supplementary Loop/TAL-associated papillary tissue-context projection for MAGMA-prioritized modules.",
            "disallowed_claim": "Plaque-specific validation; genetic risk localizes to plaque; causal papillary niche; spatial injury or mineralization mechanism.",
            "main_or_supplementary": "supplementary",
            "figure": "Figure 5",
            "manuscript_section": "Results 5; Methods: spatial",
            "notes": "No plaque, mineral, lesion, calcification or fibrosis ROI annotation was available.",
        },
        {
            "evidence_layer": "overall_integrated_model",
            "main_result": "The evidence layers place MAGMA-prioritized KSD genes in a Loop/TAL-associated renal papillary context and an injury/remodeling-associated bulk papillary disease background.",
            "support_strength": "moderate_context_support",
            "allowed_claim": "MAGMA-prioritized kidney stone risk genes map to a Loop/TAL-associated renal papillary context and an injury/remodeling-associated bulk papillary disease background, supported by conservative single-nucleus, bulk and spatial context analyses.",
            "disallowed_claim": "Causal or validated genes; causal cell type or niche; plaque localization; papilla-specific TWAS; therapeutic target validation; SMR/coloc-supported genes.",
            "main_or_supplementary": "main_synthesis",
            "figure": "Figure 1",
            "manuscript_section": "Abstract; Results 6; Discussion; Conclusions",
            "notes": "This is a post-GWAS context-mapping model, not a causal validation model.",
        },
    ]


def figure_rows() -> list[dict[str, str]]:
    return [
        {
            "figure_id": "Figure_1_integrated_study_design",
            "figure_role": "Integrated workflow and evidence hierarchy",
            "primary_message": "Separate primary genetic prioritization, main context layers and supplementary/interpretive layers in one bounded workflow.",
            "allowed_claim": "The study performs conservative post-GWAS renal papillary context mapping across genetic, snRNA, bulk and spatial resources.",
            "disallowed_claim": "Linear validation cascade; causal convergence; every layer independently validates the same mechanism.",
            "must_show_limitation": "EUR LD reference; TWAS proxy; no claim-grade SMR/coloc; no lesion ROI; context mapping rather than causality.",
            "current_status": "planned_not_generated",
            "source_legend_file": "docs/revision/stage7A_manuscript_integration_blueprint/figure_legends_integrated_v0.1.md",
            "source_figure_file": "not_generated_in_stage7A",
            "notes": "Generate only after the Stage 7B text architecture is accepted.",
        },
        {
            "figure_id": "Figure_2_snRNA_context",
            "figure_role": "Primary donor-level cell-context localization and robustness",
            "primary_message": "Primary MAGMA modules show a moderate donor-level Loop/TAL-associated pattern with partial matched-random support.",
            "allowed_claim": "Moderate donor-level Loop/TAL-associated snRNA context, robust to donor and driver-removal sensitivity analyses.",
            "disallowed_claim": "Strong cell-type enrichment; causal cell type; plaque nucleation site; cell-level inference.",
            "must_show_limitation": "Four donors; descriptive ranks; partial matched-random support; no causal or lesion inference.",
            "current_status": "source_data_backed_conservative_draft_v0.2",
            "source_legend_file": "docs/revision/stage4C2R_draft_figures_v0.2/draft_figure2_figure3_legends_v0.2.md",
            "source_figure_file": "results/figures/revision/stage4C2R_draft_figures_v0.2/figure2_snRNA_context_draft_v0.2.pdf",
            "notes": "Main context figure; retain the current claim boundary.",
        },
        {
            "figure_id": "Figure_3_candidate_evidence_model",
            "figure_role": "Two-axis gene evidence reporting and exemplar boundaries",
            "primary_message": "MAGMA priority and Kidney_Cortex TWAS proxy are separate axes; curated exemplars are interpretive flags.",
            "allowed_claim": "The evidence model transparently distinguishes genetic priority, proxy support and biological interpretation.",
            "disallowed_claim": "Causal tiering; high-confidence genes; validated exemplars; therapeutic targets.",
            "must_show_limitation": "One-SNP TWAS burden; no claim-grade SMR/coloc; exemplar status does not upgrade evidence.",
            "current_status": "source_data_backed_conservative_draft_v0.2",
            "source_legend_file": "docs/revision/stage4C2R_draft_figures_v0.2/draft_figure2_figure3_legends_v0.2.md",
            "source_figure_file": "results/figures/revision/stage4C2R_draft_figures_v0.2/figure3_candidate_evidence_draft_v0.2.pdf",
            "notes": "Replaces all P1 and single-axis tier language.",
        },
        {
            "figure_id": "Figure_4_GSE73680_bulk_context",
            "figure_role": "Paired bulk disease-context analysis and attenuation boundary",
            "primary_message": "Primary modules shift directionally in paired papillary bulk tissue, but support attenuates under composition/injury/remodeling-aware sensitivity models.",
            "allowed_claim": "Injury/remodeling-associated paired bulk disease-context support for MAGMA-prioritized modules.",
            "disallowed_claim": "Independent validation; composition-independent support; cell-type-specific disease response; mediation or mechanism.",
            "must_show_limitation": "Bulk tissue; proxy signatures; attenuation; no deconvolution or causal inference.",
            "current_status": "source_data_backed_conservative_draft_v0.1",
            "source_legend_file": "docs/revision/stage5C1_gse73680_figure4_draft/draft_figure4_legend_v0.1.md",
            "source_figure_file": "results/figures/revision/stage5C1_gse73680_figure4_draft/figure4_gse73680_bulk_context_draft_v0.1.pdf",
            "notes": "Main disease-background figure, not a validation cohort figure.",
        },
        {
            "figure_id": "Figure_5_spatial_TWAS_context",
            "figure_role": "Supplementary spatial context and TWAS proxy boundary",
            "primary_message": "Spatial projection provides moderate supplementary Loop/TAL context, while TWAS remains a Kidney_Cortex proxy.",
            "allowed_claim": "Supplementary Loop/TAL-associated papillary tissue-context projection and Kidney_Cortex TWAS proxy evidence.",
            "disallowed_claim": "Plaque-specific spatial validation; causal niche; papilla-specific TWAS; SMR/coloc support.",
            "must_show_limitation": "Five sections and 7,747 spots; no ROI; injury/ECM null after adjustment; mineralization mixed; 42/51 TWAS models one-SNP.",
            "current_status": "source_data_backed_conservative_draft_v0.1",
            "source_legend_file": "docs/revision/stage6C_spatial_twas_figure5_draft/draft_figure5_legend_v0.1.md",
            "source_figure_file": "results/figures/revision/stage6C_spatial_twas_figure5_draft/figure5_spatial_twas_context_draft_v0.1.pdf",
            "notes": "Supplementary-context figure with explicit null and boundary panels.",
        },
    ]


AUDIT_SPECS = [
    ("# Post-GWAS renal papillary", "claim_strength", "major", "replace", "Post-GWAS renal papillary context mapping of kidney stone disease genetic risk", "Current title implies prioritization of a specific program more strongly than the integrated evidence supports."),
    ("structured interpretation of P1", "obsolete_terminology", "major", "replace", "two-axis candidate-gene evidence modeling", "Remove P1 throughout."),
    ("remained supported in expression-stratified", "outdated_result", "critical", "replace", "showed moderate donor-level Loop/TAL-associated patterns with partial matched-random support and robustness to driver removal", "Stage 4B2 supersedes the old blanket robustness statement."),
    ("whether prioritized KSD genes converge", "claim_strength", "major", "replace", "whether MAGMA-prioritized genes map to a renal papillary cellular context", "Use mapping language, not convergence."),
    ("### P1 candidate evidence spectrum", "obsolete_structure", "critical", "replace", "### Candidate-gene two-axis evidence model", "Stage 3R replaced P1 with independent genetic-priority and TWAS-proxy axes."),
    ("did not test causal mediation or spatial validation", "boundary_statement", "minor", "preserve_with_revision", "did not test causal mediation, identify a causal cell type or establish lesion-resolved spatial localization", "Boundary is sound but can be aligned to the canonical vocabulary."),
    ("### MAGMA-prioritized genes converge", "claim_strength", "major", "replace", "### MAGMA-prioritized modules show moderate donor-level Loop/TAL-associated snRNA context", "Heading must encode the Stage 4 claim strength."),
    ("as the strongest cellular expression context", "claim_strength", "major", "replace", "as a moderate donor-level Loop/TAL-associated context with partial matched-random support", "Avoid strongest/strong support language."),
    ("exceeded expression-stratified random expectation", "outdated_result", "critical", "replace", "showed partial support beyond expression/detection-matched random expectations", "The conservative delta metric, not saturated rank metrics, carries support."),
    ("**Figure 3. Gene-centric single-nucleus evidence for P1", "obsolete_figure_role", "critical", "replace", "**Figure 3. Two-axis candidate-gene evidence model and curated exemplar boundaries**", "Current Figure 3 is no longer a P1 gene-centric cell-expression panel."),
    ("### GSE73680 supports module-level disease-context association", "claim_strength", "major", "replace", "### GSE73680 provides attenuated injury/remodeling-associated paired bulk disease context", "The subsection title must disclose attenuation."),
    ("remained robust after removing the six", "outdated_result", "critical", "replace", "retained overall direction but attenuated under composition/injury/remodeling-aware sensitivity models", "Stage 5B2/C1 supersedes old random-set and leave-one-gene emphasis."),
    ("### Functional context and injury coupling", "obsolete_structure", "critical", "remove_or_merge", "Integrate injury/remodeling coupling into the GSE73680 subsection", "Current Figure 5 has been reassigned to spatial/TWAS context."),
    ("**Figure 5. Functional context and injury coupling", "obsolete_figure_role", "critical", "replace", "**Figure 5. Spatial and TWAS analyses provide supplementary renal papillary context for MAGMA-prioritized modules**", "The old Figure 5 role is structurally obsolete."),
    ("The analysis tested 5963 genes", "incorrect_count", "critical", "replace", "The analysis tested 5,989 genes and identified 51 FDR-supported genes, including 42 one-SNP and nine multi-SNP models; 47 overlapped MAGMA support.", "Use the frozen Stage 2 audit counts."),
    ("Multi-SNP TWAS FDR plus MAGMA overlap was limited to eight genes", "incorrect_count", "critical", "replace", "Nine FDR-supported TWAS genes used multi-SNP models, and 47 of 51 FDR-supported genes overlapped MAGMA support.", "Do not conflate multi-SNP totals with overlap subsets."),
    ("**Supplementary Figure S7.", "obsolete_figure_role", "major", "replace", "Integrate the spatial and TWAS layers into main Figure 5 with supplementary-strength interpretation.", "Stage 6C produced the current Figure 5."),
    ("showed moderate co-distribution with Loop/TAL and ion/mineral-handling", "outdated_result", "critical", "replace", "supported moderate supplementary Loop/TAL-associated co-distribution; injury/ECM was not retained and mineralization/remodeling was mixed after adjustment", "The earlier ion/mineral summary is too positive."),
    ("KSD genetic risk converges", "claim_strength", "major", "replace", "MAGMA-prioritized KSD genes map", "Integrated interpretation should use context-mapping language."),
    ("The independent disease-context analysis", "claim_strength", "critical", "replace", "The paired bulk disease-context analysis", "GSE73680 is attenuated support, not independent validation."),
    ("Patient-aware paired shifts, random-set benchmarking", "outdated_result", "critical", "replace", "Paired and patient-aware estimates were directionally positive, with attenuation in composition/injury/remodeling-aware sensitivity analyses", "Stage 5 claim lock controls this sentence."),
    ("with Loop/TAL and ion/mineral-handling signatures", "outdated_result", "critical", "replace", "with moderate supplementary Loop/TAL context; mineralization/remodeling remained mixed and descriptive", "Retain adjusted negative results."),
    ("and P1 detection in spatial spots", "obsolete_terminology", "moderate", "replace", "and curated-exemplar detection in spatial spots", "Detection indicates assay resolvability only."),
    ("provides independent module-level disease-context support", "claim_strength", "critical", "replace", "provides attenuated injury/remodeling-associated paired bulk disease-context support", "Conclusion must match Stage 5C1."),
    ("P1: priority-1", "obsolete_terminology", "minor", "remove", "Remove P1 from the abbreviation list.", "P1 is no longer part of the manuscript architecture."),
]


def audit_rows(lines: list[str]) -> list[dict[str, str]]:
    rows = []
    for needle, problem, severity, action, replacement, notes in AUDIT_SPECS:
        matches = [(i + 1, text.strip()) for i, text in enumerate(lines) if needle.lower() in text.lower()]
        if not matches:
            raise RuntimeError(f"Audit anchor not found: {needle}")
        line_no, text = matches[0]
        rows.append({
            "line_or_section": f"line {line_no}",
            "current_text": text,
            "problem_type": problem,
            "severity": severity,
            "recommended_action": action,
            "replacement_claim": replacement,
            "notes": notes,
        })
    return rows


def supplement_rows() -> list[dict[str, str]]:
    entries = [
        ("Supplementary Table 1", "GWAS quality-control and locus reconstruction manifest", "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv", "Methods: GWAS QC; Results 1", "yes", "Supports all headline GWAS row and locus counts."),
        ("Supplementary Table 2", "MAGMA reproducibility audit", "results/tables/revision/stage2_genetic/magma_output_audit.tsv", "Methods: MAGMA; Results 1", "yes", "Record executable, version, gene locations, LD reference, tested genes and significant counts."),
        ("Supplementary Table 3", "Kidney_Cortex proxy TWAS one-SNP sensitivity audit", "results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv", "Methods: TWAS; Results 2 and 5", "yes", "Report 5,989 tested, 51 FDR, 42 one-SNP and 9 multi-SNP models."),
        ("Supplementary Table 4", "SMR and colocalization feasibility boundary", "results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv", "Methods: evidence boundary; Discussion 6", "yes", "Documents why SMR/coloc is not used as claim-grade evidence."),
        ("Supplementary Table 5", "Candidate-gene two-axis evidence model", "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv", "Results 2; Figure 3", "yes", "Mutually exclusive R1-R6 reporting groups."),
        ("Supplementary Table 6", "Curated biological exemplar panel", "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv", "Results 2; Figure 3", "no", "Interpretive flag only; does not upgrade evidence."),
        ("Supplementary Table 7", "Donor-by-compartment snRNA module scores", "results/tables/revision/stage4B1_scrna_donor_level/scrna_donor_compartment_module_scores.tsv", "Results 3; Figure 2", "yes", "Primary donor-level source table."),
        ("Supplementary Table 8", "Matched-random snRNA benchmark summary", "results/tables/revision/stage4B2_scrna_robustness/scrna_random_set_benchmark_summary.tsv", "Results 3; Figure 2", "yes", "Supports partial rather than strong matched-random evidence."),
        ("Supplementary Table 9", "Known-driver removal sensitivity", "results/tables/revision/stage4B2_scrna_robustness/scrna_known_driver_removal_sensitivity.tsv", "Results 3; Figure 2", "yes", "Documents panel-level driver-removal robustness."),
        ("Supplementary Table 10", "GSE73680 paired module-response summary", "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta_summary.tsv", "Results 4; Figure 4", "yes", "Paired response in 26 patients."),
        ("Supplementary Table 11", "GSE73680 composition-aware retention summary", "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_adjustment_retention_summary.tsv", "Results 4; Figure 4", "yes", "Primary attenuation classification across sensitivity models."),
        ("Supplementary Table 12", "Spatial section-level co-distribution consistency", "results/tables/revision/stage6B1_spatial_context_projection/spatial_section_level_consistency_summary.tsv", "Results 5; Figure 5", "yes", "Includes adjusted support classes across five sections."),
        ("Supplementary Table 13", "Spatial and TWAS final claim lock", "results/tables/revision/stage6C_spatial_twas_figure5_draft/spatial_twas_final_claim_lock_stage6C.tsv", "Results 5; Figure 5", "yes", "Preserves no-ROI, adjusted-null and TWAS-proxy boundaries."),
        ("Source Data Figure 2-3", "Figure 2 and Figure 3 source-data manifest", "results/tables/revision/stage4C2R_draft_figures_v0.2/figure_source_data_manifest_v0.2.tsv", "Figures 2-3", "yes", "Map every panel to its source extract."),
        ("Source Data Figure 4", "Figure 4 source-data manifest", "results/tables/revision/stage5C1_gse73680_figure4_draft/figure4_source_data_manifest_v0.1.tsv", "Figure 4", "yes", "Map every panel to its source extract."),
        ("Source Data Figure 5", "Figure 5 source-data manifest", "results/tables/revision/stage6C_spatial_twas_figure5_draft/figure5_source_data_manifest_v0.1.tsv", "Figure 5", "yes", "Map every panel to its source extract."),
    ]
    return [dict(zip(["supplementary_item", "proposed_title", "source_file", "main_text_reference", "required_for_reproducibility", "notes"], item)) for item in entries]


def main() -> None:
    missing = [path for path in REQUIRED_INPUTS if not (ROOT / path).is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))
    lines = MANUSCRIPT.read_text(encoding="utf-8").splitlines()

    outputs = [
        ("integrated_evidence_hierarchy_v0.1.tsv", ["evidence_layer", "main_result", "support_strength", "allowed_claim", "disallowed_claim", "main_or_supplementary", "figure", "manuscript_section", "notes"], evidence_rows()),
        ("figure_claim_map_v0.1.tsv", ["figure_id", "figure_role", "primary_message", "allowed_claim", "disallowed_claim", "must_show_limitation", "current_status", "source_legend_file", "source_figure_file", "notes"], figure_rows()),
        ("current_manuscript_claim_audit_v0.1.tsv", ["line_or_section", "current_text", "problem_type", "severity", "recommended_action", "replacement_claim", "notes"], audit_rows(lines)),
        ("supplementary_table_manifest_plan_v0.1.tsv", ["supplementary_item", "proposed_title", "source_file", "main_text_reference", "required_for_reproducibility", "notes"], supplement_rows()),
    ]
    for name, fields, rows in outputs:
        write_tsv(OUT / name, fields, rows)

    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_lines = [
        f"timestamp={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"root={ROOT}",
        f"manuscript={MANUSCRIPT.relative_to(ROOT)}",
        f"manuscript_lines={len(lines)}",
        f"required_inputs={len(REQUIRED_INPUTS)}",
        "missing_inputs=0",
    ]
    log_lines.extend(f"output={name}\trows={len(rows)}" for name, _, rows in outputs)
    LOG.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
