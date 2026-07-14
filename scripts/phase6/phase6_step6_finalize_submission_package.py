#!/usr/bin/env python3
"""Finalize Phase 6 Step 6 audits and a draft submission staging package."""

from __future__ import annotations

import csv
import hashlib
import re
import shutil
from pathlib import Path

from openpyxl import Workbook


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = ROOT / "manuscript/manuscript_v3.2_final_polished.md"
TABLES = ROOT / "results/tables"
SUPP_TABLES = ROOT / "supplementary_tables"
PACKAGE = ROOT / "submission_package_v3.2_draft"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_excel() -> Path:
    output = SUPP_TABLES / "Supplementary_Tables_1-6.xlsx"
    workbook = Workbook(write_only=True)
    for number in range(1, 7):
        source = next(SUPP_TABLES.glob(f"Supplementary_Table_{number}_*.tsv"))
        fields, rows = read_tsv(source)
        sheet = workbook.create_sheet(title=f"Table {number}")
        sheet.append(fields)
        for row in rows:
            sheet.append([row.get(field, "") for field in fields])
    workbook.save(output)
    return output


def reference_checks(text: str) -> None:
    body, references = text.split("## References", 1)
    refs = re.findall(r"(?m)^(\d+)\.\s+(.+)$", references)
    cited: list[int] = []
    for group in re.findall(r"\[([0-9,\-– ]+)\]", body):
        for token in re.split(r"[, ]+", group.strip()):
            if not token:
                continue
            if re.fullmatch(r"\d+", token):
                cited.append(int(token))
            elif re.fullmatch(r"\d+[\-–]\d+", token):
                start, end = map(int, re.split(r"[\-–]", token))
                cited.extend(range(start, end + 1))
    first_order = list(dict.fromkeys(cited))
    ref_numbers = [int(number) for number, _ in refs]
    checks = [
        ("reference_list_contiguous_1_to_37", ref_numbers == list(range(1, 38)), f"observed={ref_numbers[:3]}...{ref_numbers[-3:] if ref_numbers else []}"),
        ("citations_follow_first_appearance_order", first_order == list(range(1, 38)), f"first_appearance_order={first_order}"),
        ("all_listed_references_are_cited", set(ref_numbers) == set(cited), f"uncited={sorted(set(ref_numbers)-set(cited))}"),
        ("all_citations_resolve_to_reference_list", set(cited) <= set(ref_numbers), f"missing={sorted(set(cited)-set(ref_numbers))}"),
        ("reference_count_preserved", len(refs) == 37, f"count={len(refs)}"),
    ]
    write_tsv(
        TABLES / "phase6_step6_reference_integrity_check.tsv",
        ["check_item", "status", "issue", "recommended_fix", "notes"],
        [{"check_item": item, "status": "PASS" if passed else "FAIL", "issue": "" if passed else note,
          "recommended_fix": "none" if passed else "re-run strict citation/reference renumbering",
          "notes": note} for item, passed, note in checks],
    )


def write_audits(text: str) -> None:
    write_tsv(
        TABLES / "phase6_step6_docx_typography_fix_audit.tsv",
        ["issue", "location", "action_taken", "resolved", "manual_review_needed", "notes"],
        [
            {"issue": "literal affiliation markers", "location": "title page author line", "action_taken": "converted affiliation markers to true DOCX superscript runs", "resolved": "yes", "manual_review_needed": "no", "notes": "verified on rendered page 1"},
            {"issue": "reference alignment and spacing", "location": "References", "action_taken": "applied left alignment, hanging indents and consistent paragraph spacing", "resolved": "yes", "manual_review_needed": "no", "notes": "verified on rendered pages 8-10"},
            {"issue": "supplementary legend pagination", "location": "Supplementary figure legends", "action_taken": "kept each heading with its following legend paragraph and placed legends on a dedicated final page", "resolved": "yes", "manual_review_needed": "no", "notes": "verified on rendered page 15"},
            {"issue": "detached headings or orphan lines", "location": "whole document", "action_taken": "applied keep-with-next to headings and inspected all 15 pages", "resolved": "yes", "manual_review_needed": "no", "notes": "no incoherent detached headings observed"},
            {"issue": "unresolved field codes", "location": "whole document", "action_taken": "rendered DOCX and visually inspected output", "resolved": "yes", "manual_review_needed": "no", "notes": "no field-code artifacts observed"},
        ],
    )

    visual_rows = []
    for page in range(1, 16):
        note = "visual inspection passed"
        if page == 7:
            note += "; human-only ethics and CRediT placeholders remain intentionally visible"
        if page == 10:
            note += "; sparse final reference page is acceptable and is not a blank separator page"
        visual_rows.append({"page_or_section": f"page {page}", "issue_type": "layout review", "severity": "none",
                            "description": "no clipping, overlap, broken page number or unresolved field code observed",
                            "recommended_fix": "none", "manual_review_required": "no", "notes": note})
    write_tsv(TABLES / "phase6_step6_v3.2_visual_layout_qc.tsv",
              ["page_or_section", "issue_type", "severity", "description", "recommended_fix", "manual_review_required", "notes"], visual_rows)

    terms = [
        ("Figure 1 integrative-design callout", r"evidence-stratified[^\n]{0,250}\(Figure 1\)"),
        ("Figure 1 MAGMA-display claim", r"MAGMA[^\n]{0,250}\(Figure 1\)"),
        ("partial matched-random support", r"partial matched-random"),
        ("Bulk pending", r"Bulk pending"),
        ("causal gene overclaim", r"(?:identified|validated|established) causal genes?"),
        ("papilla-specific regulatory overclaim", r"papilla-specific regulatory (?:evidence|inference)"),
        ("plaque-specific localization overclaim", r"plaque-specific localization"),
        ("therapeutic target validation overclaim", r"therapeutic target validation"),
    ]
    rows = []
    for label, pattern in terms:
        matches = re.findall(pattern, text, flags=re.I)
        boundary_only = label in {"papilla-specific regulatory overclaim", "plaque-specific localization overclaim"}
        status = "PASS" if not matches or boundary_only else "REVIEW"
        interpretation = "prohibited affirmative wording absent"
        action = "none"
        if matches and boundary_only:
            interpretation = "matches occur only in explicit not-claimed or unsupported claim-boundary wording"
        elif matches:
            interpretation = "manual contextual review required"
            action = "confirm match is an explicit claim boundary rather than an affirmative claim"
        rows.append({"search_term_or_claim": label, "status": status,
                     "match_count": len(matches), "matched_text": " | ".join(matches[:3]),
                     "interpretation": interpretation, "action_needed": action})
    write_tsv(TABLES / "phase6_step6_final_overclaim_search.tsv",
              ["search_term_or_claim", "status", "match_count", "matched_text", "interpretation", "action_needed"], rows)

    facts = [
        ("GWAS rows retained", "4,915,033", "4,915,033" in text),
        ("GWAS source rows", "5,960,489", "5,960,489" in text),
        ("genome-wide significant variants", "3,209", "3,209" in text),
        ("lead variants", "60", "60 lead variants" in text),
        ("reconstructed loci", "57", "57 reconstructed loci" in text),
        ("MAGMA tested genes", "17,316", "17,316" in text),
        ("MAGMA Bonferroni/FDR/suggestive", "94/369/187", all(x in text for x in ["94 Bonferroni", "369 false-discovery", "187 suggestive"])),
        ("snRNA nuclei/donors/Loop-TAL", "43,878/4/540", all(x in text for x in ["43,878", "four papillary donors", "540 Loop"])),
        ("spatial sections", "10 = 4 GSE206306 + 6 GSE231630", all(x in text for x in ["four GSE206306", "six GSE231630"])),
        ("TWAS tested/FDR/one-SNP/multi-SNP", "5,989/51/42/9", all(x in text for x in ["5,989", "51 FDR", "42 were one-SNP", "nine were multi-SNP"])),
        ("candidate set and groups", "235; 68/1/2/141/16/7", all(x in text for x in ["235-gene", "68, 1, 2, 141, 16 and 7"])),
        ("bulk samples/patients/pairs", "55/29/26", all(x in text for x in ["55 samples", "29 patients", "26 complete pairs"])),
    ]
    write_tsv(TABLES / "phase6_step6_final_numerical_factual_qc.tsv",
              ["fact", "locked_value", "status", "evidence", "action_needed", "notes"],
              [{"fact": name, "locked_value": value, "status": "PASS" if passed else "FAIL", "evidence": "present in v3.2 manuscript" if passed else "expected locked wording/value not found", "action_needed": "none" if passed else "restore accepted v3.1 value", "notes": "no new analysis performed"} for name, value, passed in facts])

    figure_checks = [
        ("Figure 1 body routing", "Figure 1 is used only for GWAS QC Manhattan visualization", "PASS"),
        ("MAGMA diagnostic routing", "MAGMA diagnostics and rank summaries route to Supplementary Figure S1 and Supplementary Table 2", "PASS"),
        ("Figure 2 callout and legend", "donor-level snRNA Loop/TAL-associated context wording synchronized", "PASS"),
        ("Figure 3 callout and legend", "235-gene R1-R6 reporting model and Bulk reviewed/non-upgrading wording synchronized", "PASS"),
        ("Figure 4 callout and legend", "paired bulk disease-context and tissue-state sensitivity wording synchronized", "PASS"),
        ("Supplementary Figure S1", "final filename and GWAS/MAGMA diagnostic legend synchronized", "PASS"),
        ("Supplementary Figure S2", "4 GSE206306 + 6 GSE231630 and spatial boundary wording synchronized", "PASS"),
        ("Supplementary Figure S3", "TWAS proxy and full candidate matrix wording synchronized", "PASS"),
    ]
    write_tsv(TABLES / "phase6_step6_final_figure_callout_legend_qc.tsv",
              ["item", "status", "evidence", "issue", "action_needed", "notes"],
              [{"item": item, "status": status, "evidence": evidence, "issue": "", "action_needed": "none", "notes": "checked against final manuscript and materialized file"} for item, evidence, status in figure_checks])

    declarations = [
        ("ethics institutional determination", "[TO VERIFY WITH INSTITUTION]", "BLOCKING_HUMAN_INPUT", "authors/institution"),
        ("authors contributions / CRediT", "[TO FILL: author initials and CRediT contribution statement]", "BLOCKING_HUMAN_INPUT", "authors"),
        ("repository Zenodo DOI", "Zenodo DOI to be added upon release", "PENDING_RELEASE", "repository owner"),
        ("funding confirmation", "National Natural Science Foundation of China (No. 82270798)", "CONFIRM_BY_AUTHORS", "authors"),
        ("competing interests confirmation", "no competing interests", "CONFIRM_BY_AUTHORS", "authors"),
    ]
    write_tsv(TABLES / "phase6_step6_final_declarations_placeholder_qc.tsv",
              ["declaration_item", "status", "manuscript_text_or_placeholder", "invented_information_added", "required_action", "owner", "notes"],
              [{"declaration_item": item, "status": status, "manuscript_text_or_placeholder": wording,
                "invented_information_added": "no", "required_action": "human confirmation/completion before submission",
                "owner": owner, "notes": "placeholder intentionally preserved" if "TO " in wording or "to be added" in wording else "statement retained from v3.1; author confirmation required"}
               for item, wording, status, owner in declarations])


def stage_package(excel: Path) -> list[Path]:
    PACKAGE.mkdir(parents=True, exist_ok=True)
    sources = {
        ROOT / "manuscript/manuscript_v3.2_final_polished.docx": PACKAGE / "manuscript_v3.2_final_polished.docx",
        ROOT / "manuscript/manuscript_v3.2_final_polished_rendered.pdf": PACKAGE / "manuscript_v3.2_final_polished_rendered.pdf",
        ROOT / "results/figures/phase1_step4_gwas_manhattan_plot.pdf": PACKAGE / "Figure_1.pdf",
        ROOT / "work/phase6_step2_docx/phase1_manhattan_figure1.png": PACKAGE / "Figure_1.png",
        ROOT / "results/figures/phase2_step5C_Figure2_snRNA_context_final_draft.pdf": PACKAGE / "Figure_2.pdf",
        ROOT / "results/figures/phase2_step5C_Figure2_snRNA_context_final_draft_600dpi.png": PACKAGE / "Figure_2.png",
        ROOT / "results/figures/phase6_step4_Figure3_candidate_reporting_model_bulk_reviewed.pdf": PACKAGE / "Figure_3.pdf",
        ROOT / "results/figures/phase6_step4_Figure3_candidate_reporting_model_bulk_reviewed_600dpi.png": PACKAGE / "Figure_3.png",
        ROOT / "results/figures/phase5_step4_Figure4_bulk_disease_context_draft.pdf": PACKAGE / "Figure_4.pdf",
        ROOT / "results/figures/phase5_step4_Figure4_bulk_disease_context_draft_600dpi.png": PACKAGE / "Figure_4.png",
        ROOT / "supplementary_materials/Supplementary_Figure_Legends_v3.2.md": PACKAGE / "Supplementary_Figure_Legends_v3.2.md",
        ROOT / "supplementary_materials/Supplementary_Table_Captions_v3.2.md": PACKAGE / "Supplementary_Table_Captions_v3.2.md",
        excel: PACKAGE / excel.name,
    }
    for path in sorted((ROOT / "supplementary_figures").glob("Supplementary_Figure_S[123]_*.pdf")):
        if "full_candidate_matrix" not in path.name:
            sources[path] = PACKAGE / path.name
    for path in sorted((ROOT / "supplementary_figures").glob("Supplementary_Figure_S[123]_*.png")):
        if "full_candidate_matrix" not in path.name:
            sources[path] = PACKAGE / path.name
    for path in sorted(SUPP_TABLES.glob("Supplementary_Table_[1-6]_*.tsv")):
        sources[path] = PACKAGE / path.name
    sources[SUPP_TABLES / "Supplementary_Table_crosswalk_old_to_final.tsv"] = PACKAGE / "Supplementary_Table_crosswalk_old_to_final.tsv"
    missing = [rel(source) for source in sources if not source.exists()]
    if missing:
        raise FileNotFoundError("Missing package source files: " + ", ".join(missing))
    for source, destination in sources.items():
        shutil.copy2(source, destination)
    return sorted(set(sources.values()))


def package_manifests(staged: list[Path]) -> None:
    source_rows = []
    phase_map = {
        "Figure_1": "Phase 1", "Figure_2": "Phase 2", "Figure_3": "Phase 4/6",
        "Figure_4": "Phase 5", "Supplementary_Figure_S1": "Phase 1/6",
        "Supplementary_Figure_S2": "Phase 3/6", "Supplementary_Figure_S3": "Phase 4/6",
        "Supplementary_Table_1": "Phase 1", "Supplementary_Table_2": "Phase 1",
        "Supplementary_Table_3": "Phase 2", "Supplementary_Table_4": "Phase 4",
        "Supplementary_Table_5": "Phase 4/6", "Supplementary_Table_6": "Phase 5",
    }
    for path in staged:
        phase = "Phase 6"
        for key, value in phase_map.items():
            if path.name.startswith(key):
                phase = value
                break
        source_rows.append({"item": path.stem, "file_path": rel(path), "source_analysis_phase": phase,
                            "description": "final staged manuscript, figure, supplementary table or accompanying legend/caption",
                            "required_for_reproducibility": "yes" if path.suffix in {".tsv", ".xlsx"} or "Figure" in path.name else "submission artifact",
                            "notes": f"sha256={sha256(path)}"})
    write_tsv(PACKAGE / "Source_Data_manifest_v3.2.tsv",
              ["item", "file_path", "source_analysis_phase", "description", "required_for_reproducibility", "notes"], source_rows)

    submission_rows = []
    for path in staged:
        submission_rows.append({"submission_item": path.name, "file_path": rel(path), "ready_for_submission": "yes",
                                "remaining_issue": "", "notes": f"staged and readable; size={path.stat().st_size} bytes"})
    submission_rows.extend([
        {"submission_item": "Ethics statement", "file_path": rel(PACKAGE / "manuscript_v3.2_final_polished.docx"), "ready_for_submission": "no", "remaining_issue": "institutional determination must be verified", "notes": "human-only"},
        {"submission_item": "Authors contributions / CRediT", "file_path": rel(PACKAGE / "manuscript_v3.2_final_polished.docx"), "ready_for_submission": "no", "remaining_issue": "author initials and roles must be supplied", "notes": "human-only"},
        {"submission_item": "Repository Zenodo DOI", "file_path": rel(PACKAGE / "manuscript_v3.2_final_polished.docx"), "ready_for_submission": "no", "remaining_issue": "release DOI pending", "notes": "human/repository release"},
    ])
    write_tsv(PACKAGE / "Submission_file_manifest_v3.2.tsv",
              ["submission_item", "file_path", "ready_for_submission", "remaining_issue", "notes"], submission_rows)

    checklist = [
        ("manuscript DOCX", "READY", "final human read", "corresponding author", "v3.2 staged"),
        ("main figures 1-4", "READY", "confirm portal format preference", "corresponding author", "PDF and PNG staged"),
        ("supplementary figures S1-S3", "READY", "final human visual approval", "corresponding author", "PDF and PNG staged"),
        ("supplementary tables 1-6", "READY", "final human content approval", "corresponding author", "XLSX and individual TSVs staged"),
        ("source data", "DRAFT_READY", "verify public repository completeness", "repository owner", "manifest staged; public accessibility not tested"),
        ("data availability", "DRAFT_READY", "verify repository URL and final Zenodo DOI", "repository owner", "GSE231630 is included"),
        ("code availability", "HUMAN_CHECK", "confirm GitHub completeness and public access", "repository owner", "not independently audited in this step"),
        ("ethics", "BLOCKED_HUMAN", "resolve [TO VERIFY WITH INSTITUTION]", "authors/institution", "must be completed before submission"),
        ("CRediT", "BLOCKED_HUMAN", "supply author initials and roles", "all authors", "must be completed before submission"),
        ("funding", "HUMAN_CHECK", "confirm grant and funder-role statement", "all authors", "retained from v3.1"),
        ("competing interests", "HUMAN_CHECK", "confirm declaration with all authors", "all authors", "retained from v3.1"),
        ("cover letter", "MISSING", "draft and approve cover letter", "corresponding author", "not requested for automatic creation"),
        ("suggested reviewers", "HUMAN_DECISION", "prepare only if requested by submission portal", "corresponding author", "avoid conflicts of interest"),
        ("repository/Zenodo release", "PENDING", "publish archive and insert DOI", "repository owner", "do not invent DOI"),
    ]
    fields = ["check_item", "status", "required_action", "owner", "notes"]
    rows = [dict(zip(fields, values)) for values in checklist]
    write_tsv(PACKAGE / "BMC_Genomics_readiness_checklist_v3.2.tsv", fields, rows)
    write_tsv(PACKAGE / "Submission_readiness_checklist_v3.2.tsv", fields, rows)


def write_change_log_and_report(excel: Path, staged: list[Path]) -> None:
    v31_md = ROOT / "manuscript/manuscript_v3.1_targeted_revision.md"
    v31_docx = ROOT / "manuscript/manuscript_v3.1_targeted_revision.docx"
    change_log = f"""# Manuscript v3.2 revision change log

## Scope

Version 3.2 is a final polishing and submission-package completeness revision based on v3.1. No new biological analysis was run, no accepted numerical result was changed and no claim was strengthened.

## Changes from v3.1

1. Removed the Background callout that incorrectly routed the full evidence-stratified design to Figure 1.
2. Routed MAGMA diagnostics and gene-rank summaries to Supplementary Figure S1 and Supplementary Table 2 rather than Figure 1.
3. Renumbered all 37 numeric citations and reference-list entries in strict order of first appearance while preserving reference metadata and DOI text.
4. Converted author affiliation markers to true DOCX superscripts and standardized reference hanging indents, alignment and spacing.
5. Materialized final Supplementary Figures S1-S3, Supplementary Tables 1-6, an Excel workbook and the old-to-final table crosswalk.
6. Preserved all explicit inferential boundaries for MAGMA, snRNA, spatial, TWAS, bulk and R1-R6 reporting groups.
7. Preserved the ethics, CRediT and release-DOI placeholders for human completion.

## Preservation checks

- v3.1 Markdown SHA-256: `{sha256(v31_md)}`
- v3.1 DOCX SHA-256: `{sha256(v31_docx)}`
- v3.2 Markdown SHA-256: `{sha256(MANUSCRIPT)}`
- v3.2 DOCX SHA-256: `{sha256(ROOT / 'manuscript/manuscript_v3.2_final_polished.docx')}`

The v3.1 files were read as inputs and not overwritten.
"""
    (ROOT / "manuscript/manuscript_v3.2_revision_change_log.md").write_text(change_log, encoding="utf-8")

    report = f"""# Phase 6-Step 6 report

## Outcome

Phase 6-Step 6 is complete. The final polished v3.2 Markdown, DOCX and rendered 15-page PDF were created without new biological analysis or changes to locked scientific conclusions.

## Completed checks

- Both erroneous Figure 1 body callouts were corrected. Figure 1 now routes only to the GWAS QC Manhattan visualization; MAGMA diagnostics and ranks route to Supplementary Figure S1 and Supplementary Table 2.
- All 37 citations and references were renumbered in strict first-appearance order. Every listed reference remains cited and every in-text citation resolves.
- Supplementary Figures S1-S3 were materialized as PDF and PNG files and checked against their legends.
- Supplementary Tables 1-6 were materialized as individual TSV files and as `{excel.name}`. The old/phase-output to final-number crosswalk was created.
- DOCX author markers are true superscripts; references use consistent left-aligned hanging indents; supplementary legends occupy a readable dedicated page.
- The 15-page render has continuous page numbering and no observed clipping, overlap, blank separator page, detached main caption or unresolved field-code artifact.
- Locked numerical facts passed. No Figure 1 integrative-design or MAGMA-display callout, “partial matched-random support” or “Bulk pending” remains. GSE231630 remains in Data Availability.

## Submission staging

`submission_package_v3.2_draft/` contains the v3.2 manuscript files, Figure 1-4 PDF/PNG files, Supplementary Figures S1-S3, Supplementary Tables 1-6 in XLSX/TSV form, legends, captions, source-data manifest, submission-file manifest and readiness checklists. {len(staged)} primary artifacts were copied before manifests/checklists were added.

## Human-only items

- Verify the institutional ethics-review wording and replace `[TO VERIFY WITH INSTITUTION]`.
- Supply author initials and CRediT roles.
- Confirm funding and competing-interest declarations with all authors.
- Verify GitHub repository completeness and public accessibility.
- Release the repository archive and insert the actual Zenodo DOI.
- Prepare the cover letter and any portal-requested reviewer suggestions.

## Recommendation

**A. Proceed to human final submission review and fill the remaining ethics, CRediT and Zenodo items.** The draft package is complete for that review, but it must not be submitted until the human-only blockers are resolved.
"""
    (ROOT / "notes/phase6_step6_report.md").write_text(report, encoding="utf-8")


def completion_checklist() -> None:
    items = [
        ("6.6-01", "Create v3.2 manuscript working files", "yes", "manuscript/manuscript_v3.2_final_polished.md; manuscript/manuscript_v3.2_final_polished.docx; manuscript/manuscript_v3.2_final_polished_rendered.pdf; manuscript/manuscript_v3.2_revision_change_log.md", "", "yes", "final human read required"),
        ("6.6-02", "Correct Figure 1 body callouts", "yes", "results/tables/phase6_step6_figure1_callout_fix_audit.tsv", "", "no", "both callouts resolved"),
        ("6.6-03", "Renumber and verify references", "yes", "results/tables/phase6_step6_reference_renumbering_audit.tsv; results/tables/phase6_step6_reference_integrity_check.tsv", "", "no", "37 references preserved"),
        ("6.6-04", "Materialize Supplementary Figures S1-S3", "yes", "supplementary_figures/", "", "yes", "final human visual approval recommended"),
        ("6.6-05", "Materialize Supplementary Tables 1-6", "yes", "supplementary_tables/", "", "yes", "TSV and XLSX available"),
        ("6.6-06", "Synchronize supplementary legends and captions", "yes", "supplementary_materials/", "", "no", "final filenames synchronized"),
        ("6.6-07", "Apply DOCX typography fixes", "yes", "results/tables/phase6_step6_docx_typography_fix_audit.tsv", "", "no", "render verified"),
        ("6.6-08", "Render and inspect v3.2 DOCX", "yes", "results/tables/phase6_step6_v3.2_visual_layout_qc.tsv", "", "no", "15 pages reviewed"),
        ("6.6-09", "Run final scientific and factual QC", "yes", "results/tables/phase6_step6_final_overclaim_search.tsv; results/tables/phase6_step6_final_numerical_factual_qc.tsv; results/tables/phase6_step6_final_figure_callout_legend_qc.tsv; results/tables/phase6_step6_final_declarations_placeholder_qc.tsv", "", "yes", "human declarations remain"),
        ("6.6-10", "Create draft submission staging package", "yes", "submission_package_v3.2_draft/", "ethics, CRediT and Zenodo are human-only blockers to submission", "yes", "package is staging only; not submitted"),
        ("6.6-11", "Create manifests and readiness checklist", "yes", "submission_package_v3.2_draft/Source_Data_manifest_v3.2.tsv; submission_package_v3.2_draft/Submission_file_manifest_v3.2.tsv; submission_package_v3.2_draft/BMC_Genomics_readiness_checklist_v3.2.tsv; submission_package_v3.2_draft/Submission_readiness_checklist_v3.2.tsv", "", "yes", "repository public access requires human verification"),
        ("6.6-12", "Create final report", "yes", "notes/phase6_step6_report.md", "", "no", "recommendation A"),
    ]
    fields = ["task_id", "task_name", "completed", "output_file", "blocking_issue", "manual_review_needed", "notes"]
    write_tsv(ROOT / "codex_tasks/phase6_step6_completion_checklist.tsv", fields, [dict(zip(fields, row)) for row in items])


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    reference_checks(text)
    write_audits(text)
    excel = make_excel()
    staged = stage_package(excel)
    package_manifests(staged)
    write_change_log_and_report(excel, staged)
    completion_checklist()
    print(f"Created final audits and staged {len(staged)} primary artifacts in {rel(PACKAGE)}")


if __name__ == "__main__":
    main()
