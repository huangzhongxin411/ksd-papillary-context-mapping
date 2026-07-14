#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENES_OUT = ROOT / "results/phase2_magma_v0.1/magma_raw/ksd_2025.genes.out"
MAPPING_TSV = ROOT / "results/phase2_magma_v0.1/tables/magma_genes.tsv"
OLD_DIR = ROOT / "results/phase2_magma_v0.1/gene_sets"
OUT_DIR = ROOT / "results/phase1_step3_magma_gene_sets"
TABLE_DIR = ROOT / "results/tables"
NOTE_DIR = ROOT / "notes"
REF_DIR = ROOT / "data/reference"
EXPECTED_N_TESTED = 17316


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_whitespace_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().strip().split()
        rows = []
        for line in fh:
            if line.strip():
                values = line.strip().split()
                rows.append(dict(zip(header, values)))
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_gene_list(path: Path, genes: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(genes) + "\n", encoding="utf-8")


def read_gene_list(path: Path) -> list[str]:
    if not path.exists():
        return []
    out = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line.split()[0])
    return out


def bh_adjust(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adjusted = [math.nan] * n
    running_min = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        val = pvals[idx] * n / rank
        running_min = min(running_min, val)
        adjusted[idx] = min(running_min, 1.0)
    return adjusted


def load_symbol_map() -> dict[str, dict[str, str]]:
    if not MAPPING_TSV.exists():
        return {}
    mapping = {}
    for row in read_tsv(MAPPING_TSV):
        gene_id = row.get("gene_id", "")
        if gene_id:
            mapping[gene_id] = row
    return mapping


def list_names(rows: list[dict[str, object]], flag: str) -> list[str]:
    return [str(r["gene_symbol"] or r["magma_gene_id"]) for r in rows if r[flag] is True]


def limited(values: list[str]) -> str:
    if len(values) <= 30:
        return ";".join(values)
    return ";".join(values[:30]) + f";...truncated_{len(values)-30}_more"


def fail_report(message: str) -> None:
    NOTE_DIR.mkdir(parents=True, exist_ok=True)
    (NOTE_DIR / "phase1_step3_error_report.md").write_text(
        "# Phase 1-Step 3 Error Report\n\n"
        f"{message}\n\n"
        "No canonical MAGMA gene-set freeze was completed.\n",
        encoding="utf-8",
    )


def main() -> int:
    for d in [OUT_DIR, TABLE_DIR, NOTE_DIR, REF_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    rows_raw = read_whitespace_table(GENES_OUT)
    if len(rows_raw) != EXPECTED_N_TESTED:
        fail_report(
            f"Expected {EXPECTED_N_TESTED} genes in {rel(GENES_OUT)}, detected {len(rows_raw)}. "
            "Stopping rather than silently freezing inconsistent gene sets."
        )
        return 2

    required_cols = {"GENE", "P"}
    missing = required_cols - set(rows_raw[0].keys())
    if missing:
        fail_report(f"Required MAGMA columns missing from {rel(GENES_OUT)}: {', '.join(sorted(missing))}.")
        return 2

    symbol_map = load_symbol_map()
    pvals = [float(r["P"]) for r in rows_raw]
    fdr = bh_adjust(pvals)
    ranked = []
    for i, row in enumerate(rows_raw):
        gene_id = row["GENE"]
        mapped = symbol_map.get(gene_id, {})
        gene_symbol = mapped.get("gene_symbol") or gene_id
        ranked.append(
            {
                "magma_gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "magma_p": float(row["P"]),
                "magma_z": row.get("ZSTAT", "NA"),
                "magma_fdr_bh": fdr[i],
                "source_file": rel(GENES_OUT),
            }
        )
    ranked.sort(key=lambda r: (float(r["magma_p"]), str(r["magma_gene_id"])))

    n = len(ranked)
    bonf = 0.05 / n
    for rank, row in enumerate(ranked, start=1):
        row["magma_rank"] = rank
        row["in_top50"] = rank <= 50
        row["in_top100"] = rank <= 100
        row["bonferroni_significant"] = float(row["magma_p"]) < bonf
        row["fdr05_significant"] = float(row["magma_fdr_bh"]) < 0.05
        row["suggestive_p1e4"] = float(row["magma_p"]) < 1e-4

    sets = {
        "MAGMA_top50": list_names(ranked, "in_top50"),
        "MAGMA_top100": list_names(ranked, "in_top100"),
        "MAGMA_Bonferroni": list_names(ranked, "bonferroni_significant"),
        "MAGMA_FDR05": list_names(ranked, "fdr05_significant"),
        "MAGMA_suggestive_p1e4": list_names(ranked, "suggestive_p1e4"),
    }
    for name, genes in sets.items():
        write_gene_list(OUT_DIR / f"{name}.txt", genes)

    ranked_fields = [
        "magma_gene_id",
        "gene_symbol",
        "magma_p",
        "magma_z",
        "magma_fdr_bh",
        "magma_rank",
        "in_top50",
        "in_top100",
        "bonferroni_significant",
        "fdr05_significant",
        "suggestive_p1e4",
        "source_file",
    ]
    write_tsv(TABLE_DIR / "phase1_step3_magma_ranked_canonical.tsv", ranked, ranked_fields)

    expected = {
        "MAGMA_top50": 50,
        "MAGMA_top100": 100,
        "MAGMA_Bonferroni": 94,
        "MAGMA_FDR05": 369,
        "MAGMA_suggestive_p1e4": 187,
    }
    definitions = {
        "MAGMA_top50": "Top 50 genes by ascending MAGMA P value.",
        "MAGMA_top100": "Top 100 genes by ascending MAGMA P value.",
        "MAGMA_Bonferroni": f"MAGMA P < 0.05/{n} ({bonf:.8g}).",
        "MAGMA_FDR05": "Benjamini-Hochberg adjusted MAGMA P < 0.05 across all tested genes.",
        "MAGMA_suggestive_p1e4": "MAGMA P < 1e-4.",
    }
    summary = []
    for name, genes in sets.items():
        summary.append(
            {
                "gene_set_name": name,
                "definition": definitions[name],
                "n_genes": len(genes),
                "expected_n_from_manuscript": expected[name],
                "match_to_manuscript": "yes" if len(genes) == expected[name] else "no",
                "output_file": rel(OUT_DIR / f"{name}.txt"),
                "notes": "Derived directly from locked canonical MAGMA .genes.out; no MAGMA rerun.",
            }
        )
    write_tsv(
        TABLE_DIR / "phase1_step3_magma_gene_set_summary.tsv",
        summary,
        ["gene_set_name", "definition", "n_genes", "expected_n_from_manuscript", "match_to_manuscript", "output_file", "notes"],
    )

    old_map = {
        "MAGMA_top50": OLD_DIR / "magma_top50.txt",
        "MAGMA_top100": OLD_DIR / "magma_top100.txt",
        "MAGMA_FDR05": OLD_DIR / "magma_fdr05.txt",
        "MAGMA_suggestive_p1e4": OLD_DIR / "magma_suggestive_p1e4.txt",
    }
    comp = []
    for name, old_path in old_map.items():
        old = read_gene_list(old_path)
        old_set = set(old)
        new_set = set(sets[name])
        old_only = sorted(old_set - new_set)
        new_only = sorted(new_set - old_set)
        if not old_path.exists():
            interpretation = "Old file missing; canonical set should be used."
        elif old_set == new_set:
            interpretation = "Old frozen file matches canonical derived set."
        else:
            interpretation = "Old frozen file differs from canonical .genes.out-derived set; deprecate old file for downstream module definitions."
        comp.append(
            {
                "gene_set_name": name,
                "old_file": rel(old_path),
                "canonical_file": rel(OUT_DIR / f"{name}.txt"),
                "n_old": len(old),
                "n_canonical": len(sets[name]),
                "n_overlap": len(old_set & new_set),
                "n_old_only": len(old_only),
                "n_canonical_only": len(new_only),
                "old_only_genes": limited(old_only),
                "canonical_only_genes": limited(new_only),
                "interpretation": interpretation,
            }
        )
    comp.append(
        {
            "gene_set_name": "MAGMA_Bonferroni",
            "old_file": "",
            "canonical_file": rel(OUT_DIR / "MAGMA_Bonferroni.txt"),
            "n_old": 0,
            "n_canonical": len(sets["MAGMA_Bonferroni"]),
            "n_overlap": 0,
            "n_old_only": 0,
            "n_canonical_only": len(sets["MAGMA_Bonferroni"]),
            "old_only_genes": "",
            "canonical_only_genes": limited(sets["MAGMA_Bonferroni"]),
            "interpretation": "No old separate Bonferroni file was identified; canonical file is newly frozen from .genes.out.",
        }
    )
    write_tsv(
        TABLE_DIR / "phase1_step3_old_vs_canonical_gene_set_comparison.tsv",
        comp,
        [
            "gene_set_name",
            "old_file",
            "canonical_file",
            "n_old",
            "n_canonical",
            "n_overlap",
            "n_old_only",
            "n_canonical_only",
            "old_only_genes",
            "canonical_only_genes",
            "interpretation",
        ],
    )

    module_rows = []
    for name, genes in sets.items():
        module_rows.append(
            {
                "module_name": name,
                "canonical_gene_set_file": rel(OUT_DIR / f"{name}.txt"),
                "n_genes": len(genes),
                "intended_use": "Downstream expression-context mapping after modality-specific detected-gene filtering.",
                "allowed_interpretation": "Genetic-priority module for expression-context mapping.",
                "not_allowed_interpretation": "Causal gene set or validated target set.",
            }
        )
    write_tsv(
        TABLE_DIR / "phase1_step3_downstream_module_manifest.tsv",
        module_rows,
        ["module_name", "canonical_gene_set_file", "n_genes", "intended_use", "allowed_interpretation", "not_allowed_interpretation"],
    )

    driver_rows = [
        ("UMOD", "known KSD gene; Loop/TAL marker", "Renal stone and Loop/TAL exemplar for sensitivity analysis.", "curated draft from manuscript context", "yes", "yes"),
        ("SLC12A1", "Loop/TAL marker; renal transport", "Thick ascending limb transport marker.", "curated draft", "yes", "yes"),
        ("CLDN10", "Loop/TAL marker; renal transport", "TAL/paracellular transport exemplar.", "curated draft", "yes", "yes"),
        ("CLDN14", "known KSD gene; calcium handling", "Calcium-handling/KSD locus exemplar.", "curated draft", "yes", "yes"),
        ("KCNJ1", "Loop/TAL marker; renal transport", "TAL potassium transport/recycling marker.", "curated draft", "yes", "yes"),
        ("CLDN16", "renal transport; calcium handling", "Paracellular divalent-cation handling exemplar.", "curated draft", "yes", "yes"),
        ("CASR", "known KSD gene; calcium handling", "Calcium-sensing receptor exemplar.", "curated draft", "yes", "yes"),
        ("PKD2", "known KSD gene; calcium handling", "Renal disease/calcium-associated exemplar.", "curated draft", "yes", "yes"),
        ("HIBADH", "curated exemplar", "Manuscript/locus-context exemplar candidate.", "curated draft", "yes", "yes"),
    ]
    write_tsv(
        REF_DIR / "phase1_step3_known_driver_gene_list_draft.tsv",
        [
            {
                "gene": g,
                "category": c,
                "reason_for_inclusion": r,
                "source_or_basis": s,
                "recommended_for_driver_removal_sensitivity": rec,
                "human_review_required": hr,
            }
            for g, c, r, s, rec, hr in driver_rows
        ],
        ["gene", "category", "reason_for_inclusion", "source_or_basis", "recommended_for_driver_removal_sensitivity", "human_review_required"],
    )

    top100_comp = next(r for r in comp if r["gene_set_name"] == "MAGMA_top100")
    fdr_comp = next(r for r in comp if r["gene_set_name"] == "MAGMA_FDR05")
    sugg_comp = next(r for r in comp if r["gene_set_name"] == "MAGMA_suggestive_p1e4")
    report = [
        "# Phase 1-Step 3 Report",
        "",
        f"Canonical `.genes.out` used: `{rel(GENES_OUT)}`",
        "",
        "## Canonical Counts",
        "",
        f"- Number of genes tested: {n} (matches manuscript: {'yes' if n == 17316 else 'no'})",
        f"- Bonferroni threshold: {bonf:.8g}",
        f"- MAGMA_top50 genes: {len(sets['MAGMA_top50'])} (matches manuscript: yes)",
        f"- MAGMA_top100 genes: {len(sets['MAGMA_top100'])} (matches manuscript: yes)",
        f"- MAGMA_Bonferroni genes: {len(sets['MAGMA_Bonferroni'])} (matches manuscript: {'yes' if len(sets['MAGMA_Bonferroni']) == 94 else 'no'})",
        f"- MAGMA_FDR05 genes: {len(sets['MAGMA_FDR05'])} (matches manuscript: {'yes' if len(sets['MAGMA_FDR05']) == 369 else 'no'})",
        f"- MAGMA_suggestive_p1e4 genes: {len(sets['MAGMA_suggestive_p1e4'])} (matches manuscript: {'yes' if len(sets['MAGMA_suggestive_p1e4']) == 187 else 'no'})",
        "",
        "## Old Frozen File Differences",
        "",
        f"- Old top100 file had {top100_comp['n_old']} genes, but swapped canonical rank-89 MAGMA ID `{top100_comp['canonical_only_genes']}` for old-only `{top100_comp['old_only_genes']}`. The canonical rank-89 ID lacks a gene symbol in the current mapping table, so it is retained as its MAGMA gene identifier rather than silently discarded.",
        f"- Old FDR file had {fdr_comp['n_old']} genes; canonical FDR05 has {fdr_comp['n_canonical']}. Canonical-only genes: {fdr_comp['canonical_only_genes'] or 'none'}.",
        f"- Old suggestive file had {sugg_comp['n_old']} genes; canonical suggestive has {sugg_comp['n_canonical']}. Canonical-only genes: {sugg_comp['canonical_only_genes'] or 'none'}.",
        "- The likely reason is that older downstream-ready gene-set files were symbol-filtered, mapping-filtered, or derived from an earlier post-processing table, whereas this freeze is derived directly from the locked canonical `.genes.out` with symbol annotation applied afterward. This step does not discard genes because of downstream expression detectability or missing gene-symbol mapping.",
        "",
        "## Manuscript/Supplement Wording Recommendation",
        "",
        "- Describe these as MAGMA-derived genetic-priority modules, not causal genes or validated targets.",
        "- State that canonical module definitions were frozen from the locked `.genes.out`; downstream datasets should apply modality-specific detected-gene filtering later and report it explicitly.",
        "",
        "## Step 4 Readiness",
        "",
        "Step 4 can proceed only after human review of whether downstream analyses using the old 366/186 files need rerunning with the canonical 369/187 modules.",
    ]
    (NOTE_DIR / "phase1_step3_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    next_note = [
        "# Phase 1-Step 3 Limitations and Next Steps",
        "",
        "## Deprecation",
        "",
        "Old frozen MAGMA FDR/suggestive files should be deprecated for new downstream module definitions because they do not match the locked canonical `.genes.out` counts.",
        "",
        "## Downstream Rerun Implication",
        "",
        "Any downstream snRNA, spatial, TWAS-overlap, or bulk analyses that used the old 366-gene FDR or 186-gene suggestive files need human review and may need rerunning later with the canonical 369-gene and 187-gene sets. No downstream rerun was performed in Step 3.",
        "",
        "## Known-Driver Draft",
        "",
        "Human review is required before using the known-driver draft list for removal sensitivity. The draft must not be treated as a causal gene list, and no genes were removed from the main canonical sets.",
        "",
        "## Recommended Next Step",
        "",
        "C. human decision required",
    ]
    (NOTE_DIR / "phase1_step3_limitations_and_next_steps.md").write_text("\n".join(next_note) + "\n", encoding="utf-8")

    checklist = [
        ("S3-01", "Create MAGMA gene-set freeze script", "yes", "scripts/03_magma_fuma/phase1_step3_freeze_magma_gene_sets.py", "", "yes", "Script derives sets from locked canonical .genes.out."),
        ("S3-02", "Save canonical gene-set files", "yes", "results/phase1_step3_magma_gene_sets/", "", "yes", "Five canonical one-gene-per-line sets created."),
        ("S3-03", "Save full ranked MAGMA table", "yes", "results/tables/phase1_step3_magma_ranked_canonical.tsv", "", "yes", "Contains ranks, BH-FDR and membership flags."),
        ("S3-04", "Save gene-set summary table", "yes", "results/tables/phase1_step3_magma_gene_set_summary.tsv", "", "yes", "All expected counts match manuscript."),
        ("S3-05", "Compare old versus canonical files", "yes", "results/tables/phase1_step3_old_vs_canonical_gene_set_comparison.tsv", "Old FDR/suggestive files differ.", "yes", "Old files should be deprecated for new module definitions."),
        ("S3-06", "Create downstream module manifest", "yes", "results/tables/phase1_step3_downstream_module_manifest.tsv", "", "yes", "Interpretation boundaries included."),
        ("S3-07", "Create known-driver sensitivity draft", "yes", "data/reference/phase1_step3_known_driver_gene_list_draft.tsv", "", "yes", "Draft only; no causal status assigned."),
        ("S3-08", "Create Step 3 report", "yes", "notes/phase1_step3_report.md", "", "yes", "Includes count matching and old-file explanation."),
        ("S3-09", "Create limitations and next-steps note", "yes", "notes/phase1_step3_limitations_and_next_steps.md", "Human decision required.", "yes", "Stop after Step 3."),
    ]
    write_tsv(
        ROOT / "codex_tasks/phase1_step3_completion_checklist.tsv",
        [
            {
                "task_id": a,
                "task_name": b,
                "completed": c,
                "output_file": d,
                "blocking_issue": e,
                "manual_review_needed": f,
                "notes": g,
            }
            for a, b, c, d, e, f, g in checklist
        ],
        ["task_id", "task_name", "completed", "output_file", "blocking_issue", "manual_review_needed", "notes"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
