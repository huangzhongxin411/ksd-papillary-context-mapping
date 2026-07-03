#!/usr/bin/env python3
"""Validate Stage 7A deliverables, schemas, source paths and abstract lengths."""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs/revision/stage7A_manuscript_integration_blueprint"
TAB = ROOT / "results/tables/revision/stage7A_manuscript_integration_blueprint"
LOG = ROOT / "logs/revision/stage7A_manuscript_integration_blueprint/stage7A_validation.log"

DOCS = [
    "stage7A_report.md",
    "manuscript_architecture_v0.1.md",
    "results_rewrite_plan_v0.1.md",
    "methods_rewrite_plan_v0.1.md",
    "discussion_rewrite_plan_v0.1.md",
    "title_abstract_options_v0.1.md",
    "figure_legends_integrated_v0.1.md",
    "stage7B_full_rewrite_command_draft.md",
    "stage7A_simulated_reviewer_check.md",
]

TABLES = {
    "integrated_evidence_hierarchy_v0.1.tsv": (8, ["evidence_layer", "support_strength", "allowed_claim", "disallowed_claim"]),
    "figure_claim_map_v0.1.tsv": (5, ["figure_id", "figure_role", "current_status", "source_figure_file"]),
    "current_manuscript_claim_audit_v0.1.tsv": (25, ["line_or_section", "severity", "recommended_action", "replacement_claim"]),
    "supplementary_table_manifest_plan_v0.1.tsv": (16, ["supplementary_item", "source_file", "required_for_reproducibility"]),
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def word_count(text: str) -> int:
    text = re.sub(r"\*\*[^*]+:\*\*", "", text)
    return len(re.findall(r"\b[\w+/-]+(?:['’-][\w]+)?\b", text))


def main() -> None:
    checks: list[str] = []
    failures: list[str] = []

    for name in DOCS:
        ok = (DOC / name).is_file()
        checks.append(f"document\t{name}\t{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"missing document: {name}")

    tables: dict[str, list[dict[str, str]]] = {}
    for name, (expected_rows, required_fields) in TABLES.items():
        path = TAB / name
        if not path.is_file():
            failures.append(f"missing table: {name}")
            continue
        fields, rows = read_tsv(path)
        tables[name] = rows
        rows_ok = len(rows) == expected_rows
        fields_ok = all(field in fields for field in required_fields)
        checks.append(f"table\t{name}\trows={len(rows)}/{expected_rows}\tschema={'PASS' if fields_ok else 'FAIL'}")
        if not rows_ok:
            failures.append(f"row count mismatch: {name}")
        if not fields_ok:
            failures.append(f"schema mismatch: {name}")

    manifest = tables.get("supplementary_table_manifest_plan_v0.1.tsv", [])
    missing_sources = [row["source_file"] for row in manifest if not (ROOT / row["source_file"]).is_file()]
    checks.append(f"supplementary_sources\tpresent={len(manifest) - len(missing_sources)}/{len(manifest)}")
    failures.extend(f"missing source: {path}" for path in missing_sources)

    abstract_text = (DOC / "title_abstract_options_v0.1.md").read_text(encoding="utf-8")
    abstract_a = abstract_text.split("## Abstract A: compact journal format\n", 1)[1].split("\n## Abstract B:", 1)[0]
    abstract_b = abstract_text.split("## Abstract B: structured format\n", 1)[1]
    count_a, count_b = word_count(abstract_a), word_count(abstract_b)
    checks.append(f"abstract_A\twords={count_a}\ttarget=250-300")
    checks.append(f"abstract_B\twords={count_b}\ttarget=350-400")
    if not 250 <= count_a <= 300:
        failures.append("Abstract A outside target range")
    if not 350 <= count_b <= 400:
        failures.append("Abstract B outside target range")

    tracker = (ROOT / "docs/revision/STAGE_TRACKER.tsv").read_text(encoding="utf-8")
    tracker_ok = "stage7A_completed_stage7B_not_started" in tracker
    checks.append(f"tracker\tstage7A_completed_stage7B_not_started\t{'PASS' if tracker_ok else 'FAIL'}")
    if not tracker_ok:
        failures.append("tracker status not updated")

    content = [
        f"timestamp={datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"status={'FAIL' if failures else 'PASS'}",
        *checks,
    ]
    if failures:
        content.extend(f"failure\t{failure}" for failure in failures)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("\n".join(content) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit("; ".join(failures))


if __name__ == "__main__":
    main()
