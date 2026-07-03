#!/usr/bin/env python3
"""Validate Stage 7D deliverables without changing scientific outputs."""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs/revision/stage7D_figures_supplementary_planning"
FIGS = ROOT / "results/figures/revision/stage7D_figures_supplementary_planning"
TABLES = ROOT / "results/tables/revision/stage7D_figures_supplementary_planning"
LOG = ROOT / "logs/revision/stage7D_figures_supplementary_planning/stage7D_validation.log"


REQUIRED = [
    FIGS / "figure1_evidence_hierarchy_draft_v0.1.pdf",
    FIGS / "figure1_evidence_hierarchy_draft_v0.1.png",
    FIGS / "figure1_evidence_hierarchy_draft_v0.1.svg",
    DOCS / "figure1_qc_audit_v0.1.md",
    DOCS / "figures_1_to_5_harmonization_plan_v0.1.md",
    TABLES / "supplementary_table_assembly_manifest_v0.1.tsv",
    TABLES / "source_data_assembly_manifest_v0.1.tsv",
    TABLES / "reference_verification_checklist_v0.1.tsv",
    DOCS / "journal_target_adaptation_plan_v0.1.md",
    DOCS / "stage7D_report.md",
    DOCS / "stage7E_next_command_draft.md",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    checks: list[tuple[str, bool, str]] = []
    for path in REQUIRED:
        checks.append((f"required_file:{path.relative_to(ROOT)}", path.exists() and path.stat().st_size > 0, str(path.stat().st_size if path.exists() else 0)))

    supplemental = read_tsv(TABLES / "supplementary_table_assembly_manifest_v0.1.tsv")
    source_data = read_tsv(TABLES / "source_data_assembly_manifest_v0.1.tsv")
    references = read_tsv(TABLES / "reference_verification_checklist_v0.1.tsv")
    checks.extend([
        ("supplementary_table_count", len(supplemental) == 13, str(len(supplemental))),
        ("supplementary_sources_exist", all(row["source_exists"] == "yes" for row in supplemental), "13/13"),
        ("source_data_panel_count", len(source_data) == 31, str(len(source_data))),
        ("source_data_sources_exist", all(row["source_exists"] == "yes" for row in source_data), "31/31"),
        ("reference_count", len(references) == 41, str(len(references))),
        ("references_flagged_for_external_verification", all(row["needs_external_verification"] == "yes" for row in references), "41/41"),
    ])

    svg = (FIGS / "figure1_evidence_hierarchy_draft_v0.1.svg").read_text(encoding="utf-8")
    required_svg_text = [
        "4,915,033", "17,316", "94 Bonferroni", "1000G EUR LD", "PRIMARY PRIORITIZATION",
        "42/51", "4 donors", "MODERATE CONTEXT SUPPORT", "26 paired patients",
        "attenuated after adjustment", "ATTENUATED BULK CONTEXT", "5 sections / 7,747 spots",
        "No lesion ROI", "SUPPLEMENTARY CONTEXT", "Kidney_Cortex proxy",
    ]
    for label in required_svg_text:
        checks.append((f"figure1_text:{label}", label in svg, "present" if label in svg else "missing"))
    split_boundaries = {
        "Not causal-gene identification": ("Not causal-gene", "identification"),
        "Not a causal cell type": ("Not a causal", "cell type"),
        "Not plaque-specific localization": ("Not plaque-specific", "localization"),
    }
    for label, fragments in split_boundaries.items():
        present = all(fragment in svg for fragment in fragments)
        checks.append((f"figure1_text:{label}", present, "present as wrapped SVG text" if present else "missing"))
    checks.extend([
        ("figure1_panel_labels_A_to_F", all(re.search(fr">\s*{letter}\s*<", svg) for letter in "ABCDEF"), "A-F"),
        ("tracker_stage7D_complete", "stage7D_completed_stage7E_not_started" in (ROOT / "docs/revision/STAGE_TRACKER.tsv").read_text(encoding="utf-8"), "Stage 7E not started"),
    ])

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
