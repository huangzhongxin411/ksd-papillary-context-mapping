#!/usr/bin/env python3
"""Build Stage 7D assembly manifests from frozen revision outputs."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/tables/revision/stage7D_figures_supplementary_planning"
MANUSCRIPT = ROOT / "docs/revision/stage7C_scientific_review/manuscript_v2.1_scientific_reviewed.md"


SUPPLEMENTARY_TABLES = [
    ("Supplementary Table 1", "GWAS quality control and locus reconstruction", "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv;results/tables/revision/stage2_genetic/gwas_loci_reconciliation_input_audit.tsv", "Results: Auditable GWAS/MAGMA prioritization", "The published 59-locus table remains unavailable; retain the documented 57-locus reconstruction boundary."),
    ("Supplementary Table 2", "MAGMA output and reproducibility audit", "results/tables/revision/stage2_genetic/magma_output_audit.tsv;results/tables/revision/stage2_genetic/magma_reproducibility_file_evidence.tsv", "Results: Auditable GWAS/MAGMA prioritization", "Preserve executable, version, gene-location and 1000G EUR LD provenance."),
    ("Supplementary Table 3", "TWAS output and one-SNP model audit", "results/tables/revision/stage2_genetic/twas_output_audit.tsv;results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv", "Results: Two-axis evidence modeling", "TWAS is a GTEx Kidney_Cortex proxy; one-SNP model status must remain explicit."),
    ("Supplementary Table 4", "SMR and colocalization feasibility audit", "results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv", "Results: Two-axis evidence modeling", "Resource feasibility only; no claim-grade SMR/colocalization support."),
    ("Supplementary Table 5", "Two-axis candidate evidence model", "results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv;results/tables/revision/stage3R_gene_tiering/candidate_gene_decision_log_v0.2.tsv", "Results: Two-axis evidence modeling", "Reporting groups are evidence roles, not causal tiers."),
    ("Supplementary Table 6", "Curated biological exemplar panel", "results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv", "Results: Two-axis evidence modeling", "Exemplar status is an interpretive flag and does not alter reporting-group assignment."),
    ("Supplementary Table 7", "Single-nucleus donor-level module scores", "results/tables/revision/stage4B1_scrna_donor_level/scrna_donor_compartment_module_scores.tsv", "Results: Donor-level single-nucleus context", "Donors, not nuclei, define independent biological replication."),
    ("Supplementary Table 8", "Single-nucleus matched-random benchmark", "results/tables/revision/stage4B2_scrna_robustness/scrna_random_set_benchmark_summary.tsv;results/tables/revision/stage4B2_scrna_robustness/scrna_random_set_benchmark_donor_detail.tsv", "Results: Donor-level single-nucleus context", "Retain the partial matched-random support classification."),
    ("Supplementary Table 9", "Single-nucleus known-driver removal sensitivity", "results/tables/revision/stage4B2_scrna_robustness/scrna_known_driver_removal_sensitivity.tsv;results/tables/revision/stage4B2_scrna_robustness/scrna_gene_contribution_to_loop_tal_signal.tsv", "Results: Donor-level single-nucleus context", "Sensitivity analysis does not upgrade the evidence classification."),
    ("Supplementary Table 10", "GSE73680 paired module response", "results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta.tsv;results/tables/revision/stage5B1_gse73680_module_response/gse73680_paired_module_delta_summary.tsv;results/tables/revision/stage5B1_gse73680_module_response/gse73680_module_response_model.tsv", "Results: Paired bulk papillary context", "Twenty-six paired patients; preserve patient pairing and unadjusted model provenance."),
    ("Supplementary Table 11", "GSE73680 composition- and context-adjusted models", "results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_adjustment_retention_summary.tsv;results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_delta_adjusted_single_covariate_models.tsv;results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_delta_adjusted_compact_models.tsv;results/tables/revision/stage5B2_gse73680_composition_adjusted/gse73680_patient_fixed_effect_adjusted_models.tsv", "Results: Paired bulk papillary context", "Attenuation after adjustment is the locked interpretation."),
    ("Supplementary Table 12", "Spatial section-level co-distribution summary", "results/tables/revision/stage6B1_spatial_context_projection/spatial_section_level_consistency_summary.tsv;results/tables/revision/stage6B1_spatial_context_projection/spatial_within_section_codistribution.tsv", "Results: Section-level spatial context", "Five sections and no plaque, mineral or lesion ROI annotation."),
    ("Supplementary Table 13", "Spatial and TWAS claim-lock summary", "results/tables/revision/stage6C_spatial_twas_figure5_draft/spatial_twas_final_claim_lock_stage6C.tsv;results/tables/revision/stage6B1_spatial_context_projection/spatial_stage6B1_claim_decision_table.tsv", "Results: Section-level spatial context; Integrated evidence", "Keep spatial evidence supplementary and TWAS as a Kidney_Cortex proxy."),
]


FIGURE_MANIFESTS = [
    (ROOT / "results/tables/revision/stage7D_figures_supplementary_planning/figure1_panel_source_manifest_v0.1.tsv", "source_file"),
    (ROOT / "results/tables/revision/stage4C2R_draft_figures_v0.2/figure_source_data_manifest_v0.2.tsv", "source_table"),
    (ROOT / "results/tables/revision/stage5C1_gse73680_figure4_draft/figure4_source_data_manifest_v0.1.tsv", "source_table"),
    (ROOT / "results/tables/revision/stage6C_spatial_twas_figure5_draft/figure5_source_data_manifest_v0.1.tsv", "source_table"),
]


def paths_exist(value: str) -> bool:
    return all((ROOT / item.strip()).exists() for item in value.split(";") if item.strip())


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_supplementary_manifest() -> None:
    rows = []
    for number, title, sources, reference, notes in SUPPLEMENTARY_TABLES:
        exists = paths_exist(sources)
        rows.append({
            "supplementary_table": number,
            "proposed_title": title,
            "source_file": sources,
            "source_exists": "yes" if exists else "no",
            "ready_for_assembly": "yes_with_boundary_review" if exists else "no_missing_source",
            "main_text_reference": reference,
            "notes": notes,
        })
    write_tsv(
        OUT / "supplementary_table_assembly_manifest_v0.1.tsv",
        ["supplementary_table", "proposed_title", "source_file", "source_exists", "ready_for_assembly", "main_text_reference", "notes"],
        rows,
    )


def normalize_figure(value: str) -> str:
    match = re.search(r"figure\s*([1-5])", value, flags=re.IGNORECASE)
    return f"Figure {match.group(1)}" if match else value


def build_source_data_manifest() -> None:
    rows = []
    for manifest, source_column in FIGURE_MANIFESTS:
        with manifest.open(encoding="utf-8", newline="") as handle:
            for record in csv.DictReader(handle, delimiter="\t"):
                sources = record[source_column]
                exists = paths_exist(sources)
                readiness = record.get("ready_for_publication_source_data", "")
                note_parts = [record.get("notes", "").strip()]
                if readiness:
                    note_parts.append(f"Upstream status: {readiness}.")
                if normalize_figure(record["figure"]) == "Figure 1":
                    note_parts.append("Schematic source/provenance only; no plotted quantitative source data.")
                rows.append({
                    "figure": normalize_figure(record["figure"]),
                    "panel": record["panel"],
                    "source_data_file": sources,
                    "source_exists": "yes" if exists else "no",
                    "ready_for_source_data_package": "yes_draft_requires_final_packaging" if exists else "no_missing_source",
                    "notes": " ".join(part for part in note_parts if part),
                })
    write_tsv(
        OUT / "source_data_assembly_manifest_v0.1.tsv",
        ["figure", "panel", "source_data_file", "source_exists", "ready_for_source_data_package", "notes"],
        rows,
    )


def cited_reference_numbers(text: str) -> set[int]:
    body = text.split("## References", 1)[0]
    cited: set[int] = set()
    for block in re.findall(r"\[([0-9,\-\s]+)\]", body):
        for token in block.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                start, end = (int(x.strip()) for x in token.split("-", 1))
                cited.update(range(start, end + 1))
            else:
                cited.add(int(token))
    return cited


def build_reference_checklist() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    references = text.split("## References", 1)[1]
    cited = cited_reference_numbers(text)
    rows = []
    for match in re.finditer(r"(?m)^(\d+)\.\s+(.+)$", references):
        number = int(match.group(1))
        entry = match.group(2).strip()
        year_match = re.search(r"\b(19|20)\d{2}\b", entry)
        year = year_match.group(0) if year_match else ""
        first_author = entry.split(",", 1)[0].split(".", 1)[0].strip()
        title_start = entry.find(". ") + 2
        title_end = entry.find(". ", title_start)
        title = entry[title_start:title_end if title_end >= 0 else len(entry)].strip()
        title_short = title if len(title) <= 100 else title[:97].rstrip() + "..."
        has_doi = "doi.org/" in entry.lower()
        has_url = bool(re.search(r"https?://", entry))
        issues = []
        if not year:
            issues.append("year not parsed")
        if not (has_doi or has_url):
            issues.append("DOI/URL absent")
        if has_url and not has_doi:
            issues.append("URL present but DOI absent")
        if number not in cited:
            issues.append("not detected in main-text citations")
        rows.append({
            "reference_number": str(number),
            "first_author": first_author,
            "year": year,
            "title_short": title_short,
            "doi_or_url_present": "yes" if has_doi or has_url else "no",
            "needs_external_verification": "yes",
            "possible_issue": "; ".join(issues) if issues else "none detected by structural audit",
            "notes": "Carried forward unchanged from manuscript v2.1; verify authors, title, journal, year, volume, issue, pages/article number, DOI and citation use against an external authoritative record.",
        })
    write_tsv(
        OUT / "reference_verification_checklist_v0.1.tsv",
        ["reference_number", "first_author", "year", "title_short", "doi_or_url_present", "needs_external_verification", "possible_issue", "notes"],
        rows,
    )


def main() -> None:
    build_supplementary_manifest()
    build_source_data_manifest()
    build_reference_checklist()
    print("Stage 7D manifests generated")


if __name__ == "__main__":
    main()
