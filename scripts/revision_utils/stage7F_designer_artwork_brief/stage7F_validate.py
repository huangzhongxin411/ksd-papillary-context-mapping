#!/usr/bin/env python3
"""Validate Stage 7F designer-handoff deliverables."""

from __future__ import annotations

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs/revision/stage7F_designer_artwork_brief"
TABLES = ROOT / "results/tables/revision/stage7F_designer_artwork_brief"
LOG = ROOT / "logs/revision/stage7F_designer_artwork_brief/stage7F_validation.log"

REQUIRED = [
    DOCS / "figure1_biorender_illustrator_brief_v1.0.md",
    DOCS / "figures_2_to_5_harmonization_brief_v1.0.md",
    TABLES / "figure_text_minimization_table_v1.0.tsv",
    TABLES / "figure_panel_redesign_priority_table_v1.0.tsv",
    DOCS / "biorender_prompt_figure1_v1.0.md",
    DOCS / "illustrator_layer_manifest_figures_1_to_5_v1.0.md",
    DOCS / "figure_claim_boundary_qc_checklist_v1.0.md",
    DOCS / "stage7G_next_command_draft.md",
    DOCS / "stage7F_report.md",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def add(checks: list[tuple[str, bool, str]], name: str, condition: bool, detail: str) -> None:
    checks.append((name, condition, detail))


def contains_all(text: str, terms: list[str]) -> tuple[bool, list[str]]:
    missing = [term for term in terms if term not in text]
    return not missing, missing


def main() -> None:
    checks: list[tuple[str, bool, str]] = []
    for path in REQUIRED:
        exists = path.exists() and path.stat().st_size > 0
        add(checks, f"required_file:{path.relative_to(ROOT)}", exists, str(path.stat().st_size if path.exists() else 0))

    text_fields, text_rows = read_tsv(TABLES / "figure_text_minimization_table_v1.0.tsv")
    priority_fields, priority_rows = read_tsv(TABLES / "figure_panel_redesign_priority_table_v1.0.tsv")
    expected_text_fields = ["figure", "panel", "current_text_or_label", "problem", "recommended_short_text", "move_detail_to_legend", "claim_boundary_preserved", "notes"]
    expected_priority_fields = ["figure", "panel", "priority", "reason", "recommended_action", "requires_source_data_change", "requires_claim_change", "human_designer_needed", "notes"]
    add(checks, "text_table_schema", text_fields == expected_text_fields, f"columns={len(text_fields)}")
    add(checks, "priority_table_schema", priority_fields == expected_priority_fields, f"columns={len(priority_fields)}")
    add(checks, "text_table_count", len(text_rows) == 46, str(len(text_rows)))
    add(checks, "text_claim_boundaries_preserved", all(row["claim_boundary_preserved"] == "yes" for row in text_rows), f"{sum(row['claim_boundary_preserved'] == 'yes' for row in text_rows)}/{len(text_rows)}")
    add(checks, "text_move_values_valid", all(row["move_detail_to_legend"] in {"yes", "no", "partial"} for row in text_rows), "yes/no/partial")

    panel_counts = Counter(row["figure"] for row in priority_rows)
    expected_counts = {"Figure 1": 6, "Figure 2": 7, "Figure 3": 5, "Figure 4": 7, "Figure 5": 6}
    add(checks, "priority_table_count", len(priority_rows) == 31, str(len(priority_rows)))
    add(checks, "all_current_panels_covered", panel_counts == expected_counts, str(dict(panel_counts)))
    add(checks, "priority_values_valid", all(row["priority"] in {"high", "medium", "low"} for row in priority_rows), str(dict(Counter(row['priority'] for row in priority_rows))))
    add(checks, "no_source_data_change", all(row["requires_source_data_change"] == "no" for row in priority_rows), "31/31")
    add(checks, "no_claim_change", all(row["requires_claim_change"] == "no" for row in priority_rows), "31/31")
    add(checks, "human_designer_required", all(row["human_designer_needed"] == "yes" for row in priority_rows), "31/31")

    brief = (DOCS / "figure1_biorender_illustrator_brief_v1.0.md").read_text(encoding="utf-8")
    brief_terms = [
        "Final purpose", "One-sentence message", "Preferred composition", "Visual metaphor",
        "Genetic prioritization", "Evidence model", "Parallel context-mapping branches",
        "Integrated interpretation boundary", "Text to remove from artwork", "Palette", "Icon guidance",
        "Spacing and hierarchy", "Claim-boundary requirements", "4,915,033 QC rows", "17,316 genes tested",
        "94 Bonferroni genes", "1000G EUR LD", "4 donors", "26 patients", "7,747 spots",
        "No lesion ROI", "Kidney_Cortex proxy",
    ]
    ok, missing = contains_all(brief, brief_terms)
    add(checks, "figure1_brief_required_content", ok, f"missing={missing}")

    prompt = (DOCS / "biorender_prompt_figure1_v1.0.md").read_text(encoding="utf-8")
    forbidden = [
        "causal gene", "validated gene", "high-confidence causal gene", "independent validation",
        "strong enrichment", "plaque-specific validation", "genetic risk localizes to plaque",
        "causal cell type", "causal niche", "plaque nucleation mechanism", "papilla-specific TWAS",
        "SMR/coloc-supported", "therapeutic target",
    ]
    prompt_terms = ["16:9 white canvas", "colorblind-safe", "kidney silhouette", "Loop/TAL", "GWAS", "MAGMA", "evidence hierarchy", "three equal parallel branches", "No lesion ROI", "Kidney_Cortex proxy", *forbidden]
    ok, missing = contains_all(prompt, prompt_terms)
    add(checks, "biorender_prompt_required_and_forbidden_terms", ok, f"missing={missing}")

    harmonization = (DOCS / "figures_2_to_5_harmonization_brief_v1.0.md").read_text(encoding="utf-8")
    harmonization_terms = [
        "Figure 2: donor-level snRNA context", "Figure 3: two-axis candidate evidence model",
        "Figure 4: paired bulk disease context", "Figure 5: spatial and TWAS supplementary context",
        "Main versus supplementary recommendation", "Source-data lock", "Manual editing recommendation",
    ]
    ok, missing = contains_all(harmonization, harmonization_terms)
    add(checks, "figures_2_to_5_brief_required_content", ok, f"missing={missing}")

    layers = (DOCS / "illustrator_layer_manifest_figures_1_to_5_v1.0.md").read_text(encoding="utf-8")
    layer_terms = [
        "03_DATA_PLOTS_LINKED_LOCKED", "07_CLAIM_BOUNDARIES_EDITABLE", "F1_A_genetic_prioritization",
        "Figure 2 groups", "Figure 3 groups", "Figure 4 groups", "Figure 5 groups",
        "editable master", "vector PDF", "editable SVG", "600 dpi",
    ]
    ok, missing = contains_all(layers, layer_terms)
    add(checks, "layer_manifest_required_content", ok, f"missing={missing}")

    qc = (DOCS / "figure_claim_boundary_qc_checklist_v1.0.md").read_text(encoding="utf-8")
    qc_terms = [
        "validation cascade", "causal-gene", "plaque-specific spatial validation", "Kidney_Cortex proxy",
        "strong enrichment", "independent validation", "lesion localization", "source-data",
        "180-183 mm", "Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5",
    ]
    ok, missing = contains_all(qc, qc_terms)
    add(checks, "claim_qc_required_content", ok, f"missing={missing}")

    next_command = (DOCS / "stage7G_next_command_draft.md").read_text(encoding="utf-8")
    next_terms = [
        "Do not execute Stage 7G until", "Revised editable Figure 1 artwork", "Import the revised artwork",
        "claim-boundary checklist", "Supplementary Tables 1-13", "41 references", "manuscript v2.2",
        "Stop before DOCX",
    ]
    ok, missing = contains_all(next_command, next_terms)
    add(checks, "stage7G_command_required_content", ok, f"missing={missing}")

    report = (DOCS / "stage7F_report.md").read_text(encoding="utf-8")
    report_terms = [
        "Why the current code-generated drafts should not be final artwork", "Recommended tool route",
        "Figure 1 redesign", "Figures 2-5 priorities", "Human decisions required before Stage 7G",
        "Stage 7G has not started",
    ]
    ok, missing = contains_all(report, report_terms)
    add(checks, "stage7F_report_required_content", ok, f"missing={missing}")

    prohibited_artifacts = []
    for directory in [DOCS, TABLES]:
        for path in directory.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".docx", ".pdf", ".png", ".svg"}:
                prohibited_artifacts.append(str(path.relative_to(ROOT)))
    add(checks, "no_stage7F_artwork_or_docx_generated", not prohibited_artifacts, str(prohibited_artifacts))

    tracker = (ROOT / "docs/revision/STAGE_TRACKER.tsv").read_text(encoding="utf-8")
    add(checks, "tracker_stage7F_complete", "stage7F_completed_stage7G_waiting_for_designer_artwork" in tracker, "Stage 7G waiting for revised artwork")

    passed = all(condition for _, condition, _ in checks)
    lines = [
        f"timestamp={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"status={'PASS' if passed else 'FAIL'}",
    ]
    lines.extend(f"{'PASS' if condition else 'FAIL'}\t{name}\t{detail}" for name, condition, detail in checks)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
