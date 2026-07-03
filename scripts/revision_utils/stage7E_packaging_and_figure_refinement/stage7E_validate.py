#!/usr/bin/env python3
"""Validate Stage 7E deliverables without modifying scientific outputs."""

from __future__ import annotations

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs/revision/stage7E_packaging_and_figure_refinement"
FIGS = ROOT / "results/figures/revision/stage7E_packaging_and_figure_refinement"
TABLES = ROOT / "results/tables/revision/stage7E_packaging_and_figure_refinement"
LOG = ROOT / "logs/revision/stage7E_packaging_and_figure_refinement/stage7E_validation.log"


REQUIRED = [
    FIGS / "figure1_evidence_hierarchy_draft_v0.2.pdf",
    FIGS / "figure1_evidence_hierarchy_draft_v0.2.png",
    FIGS / "figure1_evidence_hierarchy_draft_v0.2.svg",
    DOCS / "figure1_qc_audit_v0.2.md",
    TABLES / "figures_1_to_5_visual_harmonization_checklist_v0.2.tsv",
    TABLES / "figure_harmonization_action_log_v0.1.tsv",
    TABLES / "supplementary_tables_draft_package_manifest_v0.1.tsv",
    TABLES / "source_data_draft_package_manifest_v0.1.tsv",
    TABLES / "reference_verification_checklist_v0.2.tsv",
    TABLES / "metadata_placeholder_action_table_v0.1.tsv",
    DOCS / "stage7F_next_command_draft.md",
    DOCS / "stage7E_report.md",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def add(checks: list[tuple[str, bool, str]], name: str, ok: bool, detail: str) -> None:
    checks.append((name, ok, detail))


def main() -> None:
    checks: list[tuple[str, bool, str]] = []
    for path in REQUIRED:
        exists = path.exists() and path.stat().st_size > 0
        add(checks, f"required_file:{path.relative_to(ROOT)}", exists, str(path.stat().st_size if path.exists() else 0))

    schemas = {
        "figures_1_to_5_visual_harmonization_checklist_v0.2.tsv": ["figure", "current_version", "visual_issue", "severity", "recommended_action", "requires_data_change", "requires_claim_change", "priority", "notes"],
        "figure_harmonization_action_log_v0.1.tsv": ["figure", "panel", "action_category", "action_description", "data_changed", "claim_changed", "source_data_refresh_needed", "human_approval_required", "notes"],
        "supplementary_tables_draft_package_manifest_v0.1.tsv": ["supplementary_table", "final_proposed_filename", "title", "source_file", "source_exists", "row_count_if_readable", "column_count_if_readable", "ready_for_packaging", "needs_human_review", "main_text_reference", "notes"],
        "source_data_draft_package_manifest_v0.1.tsv": ["figure", "panel", "final_proposed_filename", "source_data_file", "source_exists", "row_count_if_readable", "column_count_if_readable", "ready_for_packaging", "needs_human_review", "notes"],
        "reference_verification_checklist_v0.2.tsv": ["reference_number", "first_author", "year", "title_short", "doi_or_url_present", "needs_external_verification", "possible_issue", "reference_status", "action_needed", "priority", "notes"],
        "metadata_placeholder_action_table_v0.1.tsv": ["item", "current_placeholder", "required_human_input", "where_used", "priority", "blocking_for_submission", "notes"],
    }
    tables: dict[str, list[dict[str, str]]] = {}
    for filename, expected in schemas.items():
        fields, rows = read_tsv(TABLES / filename)
        tables[filename] = rows
        add(checks, f"schema:{filename}", fields == expected, f"columns={len(fields)}")

    visual = tables["figures_1_to_5_visual_harmonization_checklist_v0.2.tsv"]
    actions = tables["figure_harmonization_action_log_v0.1.tsv"]
    supplements = tables["supplementary_tables_draft_package_manifest_v0.1.tsv"]
    source = tables["source_data_draft_package_manifest_v0.1.tsv"]
    refs = tables["reference_verification_checklist_v0.2.tsv"]
    metadata = tables["metadata_placeholder_action_table_v0.1.tsv"]

    add(checks, "visual_checklist_count", len(visual) == 18, str(len(visual)))
    add(checks, "visual_checklist_no_data_change", all(row["requires_data_change"] == "no" for row in visual), f"{sum(row['requires_data_change'] == 'no' for row in visual)}/{len(visual)}")
    add(checks, "visual_checklist_no_claim_change", all(row["requires_claim_change"] == "no" for row in visual), f"{sum(row['requires_claim_change'] == 'no' for row in visual)}/{len(visual)}")
    add(checks, "action_log_no_data_change", all(row["data_changed"] == "no" for row in actions), f"{sum(row['data_changed'] == 'no' for row in actions)}/{len(actions)}")
    add(checks, "action_log_no_claim_change", all(row["claim_changed"] == "no" for row in actions), f"{sum(row['claim_changed'] == 'no' for row in actions)}/{len(actions)}")

    add(checks, "supplementary_count", len(supplements) == 13, str(len(supplements)))
    add(checks, "supplementary_sources_exist", all(row["source_exists"] == "yes" for row in supplements), f"{sum(row['source_exists'] == 'yes' for row in supplements)}/13")
    add(checks, "supplementary_counts_recorded", all(row["row_count_if_readable"] and row["column_count_if_readable"] for row in supplements), "13/13")
    add(checks, "supplementary_pending_final_format", all(row["ready_for_packaging"] == "ready_for_packaging_pending_final_format" for row in supplements), "13/13")

    panel_counts = Counter(row["figure"] for row in source)
    add(checks, "source_data_count", len(source) == 31, str(len(source)))
    add(checks, "source_data_panel_distribution", panel_counts == {"Figure 1": 6, "Figure 2": 7, "Figure 3": 5, "Figure 4": 7, "Figure 5": 6}, str(dict(panel_counts)))
    add(checks, "source_data_sources_exist", all(row["source_exists"] == "yes" for row in source), f"{sum(row['source_exists'] == 'yes' for row in source)}/31")
    add(checks, "source_data_counts_recorded", all(row["row_count_if_readable"] and row["column_count_if_readable"] for row in source), "31/31")

    statuses = Counter(row["reference_status"] for row in refs)
    add(checks, "reference_count", len(refs) == 41, str(len(refs)))
    add(checks, "reference_status_distribution", statuses == {"pending_external_verification": 37, "DOI_missing_or_not_available": 1, "citation_use_needs_human_review": 3}, str(dict(statuses)))
    add(checks, "all_references_preserved_for_verification", all(row["needs_external_verification"] == "yes" for row in refs), "41/41")
    flagged = {row["reference_number"]: row["reference_status"] for row in refs if row["reference_number"] in {"14", "29", "35", "36"}}
    add(checks, "required_reference_flags", flagged == {"14": "DOI_missing_or_not_available", "29": "citation_use_needs_human_review", "35": "citation_use_needs_human_review", "36": "citation_use_needs_human_review"}, str(flagged))

    expected_metadata = {
        "GWAS accession/source URL", "code/source-data repository URL", "ethics approval", "consent for publication",
        "competing interests", "funding", "authors' contributions", "acknowledgements", "target journal", "primary backup journal",
    }
    add(checks, "metadata_count_and_items", len(metadata) == 10 and {row["item"] for row in metadata} == expected_metadata, str(len(metadata)))
    add(checks, "metadata_submission_blockers", all(row["blocking_for_submission"] == "yes" for row in metadata), "10/10")

    svg = (FIGS / "figure1_evidence_hierarchy_draft_v0.2.svg").read_text(encoding="utf-8")
    required_svg = [
        "Evidence-stratified renal papillary context mapping of KSD genetic risk",
        "MAGMA output", "17,316 genes tested", "94 Bonferroni genes", "1000G EUR LD",
        "MAGMA only", "MAGMA +", "proxy TWAS", "TWAS proxy", "4 donors",
        "partial matched-random support", "26 paired patients", "attenuated after adjustment",
        "5 sections / 7,747 spots", "No lesion ROI", "Kidney_Cortex TWAS proxy",
        "Not causal-gene", "Not a causal", "Not plaque-specific", "identification", "cell type", "localization",
    ]
    for label in required_svg:
        add(checks, f"figure1_text:{label}", label in svg, "present" if label in svg else "missing")
    add(checks, "figure1_panels_A_to_F", all(f">{letter}<" in svg for letter in "ABCDEF"), "A-F")

    tracker = (ROOT / "docs/revision/STAGE_TRACKER.tsv").read_text(encoding="utf-8")
    add(checks, "tracker_stage7E_complete", "stage7E_completed_stage7F_not_started" in tracker, "Stage 7F not started")

    passed = all(ok for _, ok, _ in checks)
    lines = [
        f"timestamp={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"status={'PASS' if passed else 'FAIL'}",
    ]
    lines.extend(f"{'PASS' if ok else 'FAIL'}\t{name}\t{detail}" for name, ok, detail in checks)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
