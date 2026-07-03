#!/usr/bin/env python3
"""Build Stage 7B traceability, source-map and overclaim-audit tables."""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT = ROOT / "docs/revision/stage7B_full_manuscript_draft/manuscript_v2.0_evidence_locked_draft.md"
OUT = ROOT / "results/tables/revision/stage7B_full_manuscript_draft"
LOG = ROOT / "logs/revision/stage7B_full_manuscript_draft/stage7B_build_supporting_tables.log"


def write_tsv(name: str, fields: list[str], rows: list[dict[str, str]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def claim(section: str, text: str, layer: str, source: str, strength: str, status: str = "allowed_with_boundary", notes: str = "") -> dict[str, str]:
    return {
        "manuscript_section": section,
        "claim_text": text,
        "evidence_layer": layer,
        "source_file": source,
        "support_strength": strength,
        "allowed_claim_status": status,
        "disallowed_claim_checked": "yes",
        "notes": notes,
    }


def claim_rows() -> list[dict[str, str]]:
    s2 = "results/tables/revision/stage2_genetic"
    s3 = "results/tables/revision/stage3R_gene_tiering"
    s4a = "results/tables/revision/stage4B1_scrna_donor_level"
    s4b = "results/tables/revision/stage4B2_scrna_robustness"
    s5a = "results/tables/revision/stage5B1_gse73680_module_response"
    s5b = "results/tables/revision/stage5B2_gse73680_composition_adjusted"
    s5c = "results/tables/revision/stage5C1_gse73680_figure4_draft"
    s6a = "results/tables/revision/stage6A_spatial_twas_audit"
    s6b = "results/tables/revision/stage6B1_spatial_context_projection"
    s6c = "results/tables/revision/stage6C_spatial_twas_figure5_draft"
    return [
        claim("Abstract", "GWAS QC retained 4,915,033 of 5,960,489 rows and reconstructed 57 loci.", "GWAS_MAGMA_prioritization", f"{s2}/gwas_qc_manifest.tsv", "primary_prioritization", "allowed"),
        claim("Abstract", "MAGMA tested 17,316 genes and identified 94 Bonferroni, 369 FDR and 187 suggestive genes.", "GWAS_MAGMA_prioritization", f"{s2}/magma_output_audit.tsv", "primary_prioritization", "allowed"),
        claim("Abstract", "Of 51 FDR-supported TWAS genes, 42 were one-SNP models.", "TWAS_Kidney_Cortex_proxy", f"{s2}/twas_one_snp_model_audit.tsv", "supplementary_proxy_only", notes="The abstract explicitly restricts TWAS to a supplementary role."),
        claim("Abstract", "Primary MAGMA modules showed moderate donor-level Loop/TAL-associated patterns.", "snRNA_LoopTAL_context", f"{s4b}/loop_tal_claim_decision_table.tsv", "moderate_context_support"),
        claim("Abstract", "GSE73680 support was attenuated in composition/injury/remodeling-aware sensitivity models.", "GSE73680_bulk_disease_context", f"{s5b}/gse73680_adjustment_retention_summary.tsv", "attenuated_bulk_context_support"),
        claim("Abstract", "Spatial data provided moderate supplementary Loop/TAL co-distribution with adjusted null or mixed findings for other contexts.", "spatial_LoopTAL_tissue_context", f"{s6c}/spatial_twas_final_claim_lock_stage6C.tsv", "moderate_supplementary_context"),
        claim("Results 2.1", "QC retained 4,915,033 rows and reproducibly reconstructed 57 distance-defined loci; the published 59-locus table was unavailable.", "GWAS_MAGMA_prioritization", f"{s2}/gwas_qc_manifest.tsv; {s2}/gwas_loci_reconciliation_input_audit.tsv", "primary_prioritization", notes="The 57 versus 59 discrepancy is left unresolved."),
        claim("Results 2.1", "EUR-LD-reference MAGMA identified 94 Bonferroni, 369 FDR and 187 suggestive genes among 17,316 tested genes.", "GWAS_MAGMA_prioritization", f"{s2}/magma_output_audit.tsv; docs/revision/stage2_genetic/manuscript_ready_genetic_wording.md", "primary_prioritization"),
        claim("Results 2.2", "The two-axis model assigned 79, 1, 13, 33, 4 and 0 genes to R1-R6.", "candidate_gene_evidence_model", f"{s3}/evidence_model_summary_counts_v0.2.tsv", "primary_prioritization", "allowed"),
        claim("Results 2.2", "Kidney_Cortex TWAS identified 51 FDR genes, including 42 one-SNP and nine multi-SNP models; 47 overlapped MAGMA support.", "TWAS_Kidney_Cortex_proxy", f"{s2}/twas_one_snp_model_audit.tsv", "supplementary_proxy_only"),
        claim("Results 2.2", "Six curated exemplars provide biological interpretation without changing evidence-group assignment.", "curated_exemplar_panel", f"{s3}/curated_exemplar_panel_v0.2.tsv", "interpretive_only"),
        claim("Results 2.3", "The audited snRNA resource contained 43,878 nuclei, including 540 Loop/TAL nuclei from four donors.", "snRNA_LoopTAL_context", f"{s4a}/stage4B1_metadata_validation.tsv", "moderate_context_support", "allowed"),
        claim("Results 2.3", "Primary modules showed donor-level Loop/TAL-associated profiles retained during single-donor exclusion.", "snRNA_LoopTAL_context", f"{s4a}/scrna_donor_compartment_module_scores.tsv; {s4a}/scrna_leave_one_donor_out_module_ranks.tsv", "moderate_context_support"),
        claim("Results 2.3", "Matched-random support was partial and the pattern was robust to panel-level known-driver removal.", "snRNA_LoopTAL_context", f"{s4b}/scrna_random_set_benchmark_summary.tsv; {s4b}/scrna_known_driver_removal_sensitivity.tsv", "moderate_context_support"),
        claim("Results 2.4", "GSE73680 included 55 samples from 29 patients, with 26 paired patients.", "GSE73680_bulk_disease_context", f"{s5a}/gse73680_stage5B1_pairing_summary.tsv", "attenuated_bulk_context_support", "allowed"),
        claim("Results 2.4", "Primary modules shifted positively before composition-aware adjustment.", "GSE73680_bulk_disease_context", f"{s5a}/gse73680_paired_module_delta_summary.tsv; {s5a}/gse73680_module_response_model.tsv", "attenuated_bulk_context_support"),
        claim("Results 2.4", "All four primary modules were classified as attenuated_but_retained and coupled to bulk injury/remodeling signatures.", "GSE73680_bulk_disease_context", f"{s5b}/gse73680_adjustment_retention_summary.tsv; {s5b}/gse73680_injury_remodeling_coupling_interpretation.tsv; {s5c}/gse73680_final_claim_lock_stage5C1.tsv", "attenuated_bulk_context_support"),
        claim("Results 2.5", "Five spatial sections comprising 7,747 spots lacked lesion-resolved ROI annotation.", "spatial_LoopTAL_tissue_context", f"{s6a}/spatial_section_qc_summary.tsv; docs/revision/stage6A_spatial_twas_audit/spatial_roi_plaque_annotation_audit.md", "moderate_supplementary_context"),
        claim("Results 2.5", "Two modules were consistently positive and two were positive but attenuated for adjusted Loop/TAL co-distribution.", "spatial_LoopTAL_tissue_context", f"{s6b}/spatial_section_level_consistency_summary.tsv", "moderate_supplementary_context"),
        claim("Results 2.5", "Adjusted injury/ECM co-distribution was not retained and mineralization/remodeling was mixed.", "spatial_LoopTAL_tissue_context", f"{s6b}/spatial_section_level_consistency_summary.tsv; {s6c}/spatial_twas_final_claim_lock_stage6C.tsv", "moderate_supplementary_context", "allowed", "Negative and mixed findings are preserved."),
        claim("Results 2.6", "The integrated hierarchy assigns unequal inferential roles to MAGMA, snRNA, bulk, spatial, TWAS and exemplar layers.", "overall_integrated_model", "results/tables/revision/stage7A_manuscript_integration_blueprint/integrated_evidence_hierarchy_v0.1.tsv", "moderate_context_support"),
        claim("Results 2.6", "MAGMA-prioritized KSD genes map to a Loop/TAL-associated papillary context and injury/remodeling-associated bulk background.", "overall_integrated_model", "results/tables/revision/stage7A_manuscript_integration_blueprint/integrated_evidence_hierarchy_v0.1.tsv", "moderate_context_support"),
        claim("Discussion paragraph 1", "The most consistent integrated pattern is Loop/TAL-associated context with attenuated bulk disease-context support.", "overall_integrated_model", "results/tables/revision/stage7A_manuscript_integration_blueprint/integrated_evidence_hierarchy_v0.1.tsv", "moderate_context_support"),
        claim("Discussion paragraph 2", "Curated transport and ion-handling exemplars illustrate biological compatibility without resolving locus mediators.", "curated_exemplar_panel", f"{s3}/curated_exemplar_panel_v0.2.tsv", "interpretive_only"),
        claim("Discussion paragraph 3", "The Loop/TAL snRNA interpretation is moderate because donor sensitivity was supportive but matched-random evidence was partial.", "snRNA_LoopTAL_context", f"{s4b}/loop_tal_claim_decision_table.tsv", "moderate_context_support"),
        claim("Discussion paragraph 4", "Bulk module shifts are associated with injury/remodeling and attenuate after composition-aware adjustment.", "GSE73680_bulk_disease_context", f"{s5c}/gse73680_final_claim_lock_stage5C1.tsv", "attenuated_bulk_context_support"),
        claim("Discussion paragraph 5", "Spatial evidence is supplementary because no lesion ROI was available and non-Loop/TAL contexts were null or mixed after adjustment.", "spatial_LoopTAL_tissue_context", f"{s6c}/spatial_twas_final_claim_lock_stage6C.tsv", "moderate_supplementary_context"),
        claim("Discussion paragraph 6", "Kidney_Cortex TWAS is a one-SNP-heavy proxy and claim-grade papilla eQTL/SMR/colocalization evidence is missing.", "TWAS_Kidney_Cortex_proxy", f"{s2}/twas_one_snp_model_audit.tsv; {s2}/smr_coloc_feasibility.tsv", "supplementary_proxy_only"),
        claim("Discussion paragraph 7", "Ancestry reference, locus reconciliation, replicate scale, bulk attenuation, no-ROI spatial data and no perturbation limit interpretation.", "overall_integrated_model", "docs/revision/stage7A_manuscript_integration_blueprint/stage7A_report.md", "limitation", "allowed"),
        claim("Discussion paragraph 8", "Ancestry-matched genetics, papilla regulatory resources, lesion-resolved spatial data and donor-aware perturbation are future needs.", "overall_integrated_model", "docs/revision/stage7A_manuscript_integration_blueprint/discussion_rewrite_plan_v0.1.md", "limitation", "allowed", "Future-work statement, not a current result."),
    ]


def section_map_rows() -> list[dict[str, str]]:
    rows = [
        ("Title and Abstract", "Stage 7A title/abstract options; integrated evidence hierarchy; Stage 2-6 claim locks", "Figures 1-5 referenced conceptually", "Supplementary Tables 1-13", "Abstract is 300 words or fewer and uses only locked counts."),
        ("Introduction", "docs/revision/stage1/manuscript_v1.4_clean_working.md; References 1-41", "Figure 1", "None", "Literature statements retain existing references; no new citation metadata added."),
        ("Results 2.1 GWAS/MAGMA", "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv; results/tables/revision/stage2_genetic/magma_output_audit.tsv", "Figure 1; Figure 3", "Supplementary Tables 1-2", "Includes unresolved 57 versus 59 boundary and EUR LD limitation."),
        ("Results 2.2 evidence model/TWAS", "results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv; results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv; results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv; results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv", "Figure 3; Figure 5E", "Supplementary Tables 3-6", "Separates genetic, proxy and interpretive layers."),
        ("Results 2.3 snRNA", "Stage 4B1 donor-level tables; Stage 4B2 random, driver-removal and claim-decision tables", "Figure 2", "Supplementary Tables 7-9; Figure 2 Source Data", "Donor-level support; no cell-level inference."),
        ("Results 2.4 GSE73680", "Stage 5B1 paired/model tables; Stage 5B2 adjustment/coupling tables; Stage 5C1 claim lock", "Figure 4", "Supplementary Tables 10-11; Figure 4 Source Data", "Attenuation and bulk-signature limits remain adjacent to positive results."),
        ("Results 2.5 spatial/TWAS", "Stage 6A ROI/QC audit; Stage 6B1 section consistency; Stage 6C claim lock", "Figure 5", "Supplementary Tables 12-13; Figure 5 Source Data", "Supplementary context in main Results; no ROI inference."),
        ("Results 2.6 integrated model", "Stage 7A integrated evidence hierarchy and figure claim map", "Figure 1", "Supplementary Tables 1-13", "Synthesis introduces no new quantitative evidence."),
        ("Discussion", "Stage 2-6 claim locks; Stage 7A discussion plan; existing References 1-41", "Figures 1-5", "Supplementary Tables 1-13", "Eight paragraphs match the prescribed evidence hierarchy and limitations."),
        ("Methods", "Stage 2-6 scripts, audits, tables, logs and source-data manifests listed in Stage 7A methods plan", "Figures 1-5", "Supplementary Tables 1-13; Figure 2-5 Source Data", "No analysis beyond the locked stages is described."),
        ("Figure legends", "Locked Stage 4C2R, Stage 5C1 and Stage 6C legends; Stage 7A integrated legends", "Figures 1-5", "Figure 2-5 Source Data", "Figure 1 remains text-only; Figures 2-5 are unchanged."),
        ("Data availability and Declarations", "Current manuscript placeholders; Stage 7B instruction", "None", "None", "Repository/source URLs and author declarations remain TO FILL/VERIFY."),
        ("References", "docs/revision/stage1/manuscript_v1.4_clean_working.md", "None", "None", "All 41 entries carried forward without adding unverified references."),
    ]
    fields = ["manuscript_section", "primary_source_files", "figures_used", "supplementary_tables_needed", "notes"]
    return [dict(zip(fields, row)) for row in rows]


TERMS = [
    "causal",
    "validated",
    "validation",
    "high-confidence",
    "plaque-specific",
    "spatial validation",
    "genetic risk localizes",
    "causal cell type",
    "plaque nucleation",
    "therapeutic target",
    "SMR-supported",
    "coloc-supported",
    "papilla-specific TWAS",
    "independent validation",
    "strong enrichment",
]


def sectioned_lines(text: str) -> list[tuple[str, str]]:
    section = "Title"
    output = []
    for line in text.splitlines():
        if line.startswith("## ") or line.startswith("### "):
            section = line.lstrip("# ")
        if line.strip():
            output.append((section, line.strip()))
    return output


def acceptable_context(line: str) -> bool:
    lower = line.lower()
    boundary_markers = [
        "not ", "no ", "does not", "did not", "cannot", "rather than",
        "without", "unavailable", "missing", "future", "needed", "required",
        "do not", "doesn't", "absence", "was not", "were not", "limits",
    ]
    return any(marker in lower for marker in boundary_markers)


def overclaim_rows(text: str) -> list[dict[str, str]]:
    rows = []
    lines = sectioned_lines(text)
    for term in TERMS:
        hits = [(section, line) for section, line in lines if term.lower() in line.lower()]
        if not hits:
            rows.append({
                "term": term,
                "section": "not_applicable",
                "exact_text": "",
                "problem_status": "not_present",
                "recommended_fix": "none",
                "notes": "No occurrence in the manuscript.",
            })
            continue
        for section, line in hits:
            acceptable = acceptable_context(line) or section == "References"
            rows.append({
                "term": term,
                "section": section,
                "exact_text": line,
                "problem_status": "acceptable_context" if acceptable else "needs_revision",
                "recommended_fix": "retain as an explicit boundary statement" if acceptable else "remove or rewrite as an explicit negative boundary",
                "notes": "Reference title" if section == "References" else "Negative, limiting or future-work context" if acceptable else "Potential positive overclaim",
            })
    return rows


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    claims = claim_rows()
    source_map = section_map_rows()
    overclaims = overclaim_rows(text)

    write_tsv("claim_traceability_matrix_v2.0.tsv", ["manuscript_section", "claim_text", "evidence_layer", "source_file", "support_strength", "allowed_claim_status", "disallowed_claim_checked", "notes"], claims)
    write_tsv("section_source_file_map_v2.0.tsv", ["manuscript_section", "primary_source_files", "figures_used", "supplementary_tables_needed", "notes"], source_map)
    write_tsv("manuscript_overclaim_audit_v2.0.tsv", ["term", "section", "exact_text", "problem_status", "recommended_fix", "notes"], overclaims)

    status_counts: dict[str, int] = {}
    for row in overclaims:
        status_counts[row["problem_status"]] = status_counts.get(row["problem_status"], 0) + 1
    missing_sources = []
    for row in claims:
        for source in row["source_file"].split("; "):
            if not (ROOT / source).exists():
                missing_sources.append(source)

    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("\n".join([
        f"timestamp={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"manuscript={MANUSCRIPT.relative_to(ROOT)}",
        f"claim_rows={len(claims)}",
        f"section_map_rows={len(source_map)}",
        f"overclaim_rows={len(overclaims)}",
        f"overclaim_statuses={status_counts}",
        f"missing_claim_sources={len(set(missing_sources))}",
        *[f"missing_source={source}" for source in sorted(set(missing_sources))],
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
