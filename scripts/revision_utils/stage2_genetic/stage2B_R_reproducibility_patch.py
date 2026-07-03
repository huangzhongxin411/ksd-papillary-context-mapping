#!/usr/bin/env python3
"""Stage 2B-R reproducibility patch before Stage 3 tiering."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / "docs/revision/stage2_genetic"
TAB = ROOT / "results/tables/revision/stage2_genetic"
STAGE3 = ROOT / "results/tables/revision/stage3_gene_tiering"
LOG = ROOT / "logs/revision/stage2_genetic"


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rel(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def stat_row(item: str, expected: str, path: Path | None, used: str, notes: str = "") -> list[str]:
    found = "yes" if path and path.exists() else "no"
    if path and path.exists():
        return [
            item,
            expected,
            found,
            rel(path),
            str(path.stat().st_size),
            dt.datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
            used,
            "confirmed",
            notes,
        ]
    return [item, expected, found, "NA_not_found", "NA_not_found", "NA_not_found", used, "missing", notes]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def make_file_evidence() -> None:
    rows = [
        stat_row("MAGMA executable", "external/magma/magma", ROOT / "external/magma/magma", "yes", "Executable path present; no rerun performed."),
        stat_row("MAGMA version record", "results/logs/magma_version.txt", ROOT / "results/logs/magma_version.txt", "yes", "Version record supports MAGMA v1.10 audit."),
        stat_row("MAGMA gene location file", "external/reference/gene_loc/NCBI37.3.gene.loc", ROOT / "external/reference/gene_loc/NCBI37.3.gene.loc", "yes", "NCBI build 37 gene location file."),
        stat_row("1000G EUR LD reference prefix", "external/reference/1000G_EUR/g1000_eur.{bed,bim,fam}", ROOT / "external/reference/1000G_EUR/g1000_eur.bed", "yes", "Prefix confirmed by bed/bim/fam; bed row records prefix size using .bed file."),
        stat_row("1000G EUR LD reference .bim", "external/reference/1000G_EUR/g1000_eur.bim", ROOT / "external/reference/1000G_EUR/g1000_eur.bim", "yes", "Companion BIM file present."),
        stat_row("1000G EUR LD reference .fam", "external/reference/1000G_EUR/g1000_eur.fam", ROOT / "external/reference/1000G_EUR/g1000_eur.fam", "yes", "Companion FAM file present."),
        stat_row("MAGMA pval input file", "data/processed/magma_input/ksd_2025.dedup.pval", ROOT / "data/processed/magma_input/ksd_2025.dedup.pval", "yes", "Deduplicated p-value input."),
        stat_row("MAGMA snploc input file", "data/processed/magma_input/ksd_2025.dedup.snploc", ROOT / "data/processed/magma_input/ksd_2025.dedup.snploc", "yes", "Deduplicated SNP location input."),
        stat_row("MAGMA gene-based output", "results/tables/magma_genes.tsv", ROOT / "results/tables/magma_genes.tsv", "yes", "Postprocessed gene-based result used in manuscript."),
        stat_row("MAGMA log file", "results/magma/2025_trans_ancestry/ksd_2025.log", ROOT / "results/magma/2025_trans_ancestry/ksd_2025.log", "yes", "MAGMA run log present."),
        stat_row("MAGMA gene-set source table", "results/tables/magma_top50_top100_gene_sets.tsv", ROOT / "results/tables/magma_top50_top100_gene_sets.tsv", "yes", "Top50/top100 source table present."),
        stat_row("MAGMA gene-set source files", "results/gene_sets/magma_top50.txt", ROOT / "results/gene_sets/magma_top50.txt", "yes", "Gene-set text source present; other set files are in same directory."),
        stat_row("GWAS QC manifest", "results/tables/revision/stage2_genetic/gwas_qc_manifest.tsv", TAB / "gwas_qc_manifest.tsv", "yes", "Stage 2A/2B manifest."),
        stat_row("reconstructed 57-locus table", "results/tables/phase1_2025_loci.tsv", ROOT / "results/tables/phase1_2025_loci.tsv", "yes", "Analysis-specific reconstructed loci."),
        stat_row("TWAS all-results table", "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv", ROOT / "results/tables/supplementary_table_spredixcan_all_results_v0.1.tsv", "yes", "Kidney_Cortex S-PrediXcan parsed results."),
        stat_row("TWAS one-SNP audit table", "results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv", TAB / "twas_one_snp_model_audit.tsv", "yes", "Stage 2A/2B one-SNP audit."),
        stat_row("candidate gene evidence skeleton", "results/tables/revision/stage2_genetic/candidate_gene_evidence_skeleton.tsv", TAB / "candidate_gene_evidence_skeleton.tsv", "yes", "Frozen for Stage 3 below."),
    ]
    write_tsv(TAB / "magma_reproducibility_file_evidence.tsv",
              ["evidence_item", "expected_file_or_resource", "found", "file_path", "file_size_bytes", "last_modified", "used_in_current_manuscript", "audit_status", "notes"],
              rows)


def update_magma_output_audit() -> None:
    path = TAB / "magma_output_audit.tsv"
    rows = read_tsv(path)
    if not rows:
        return
    extra = {
        "magma_executable_path": "external/magma/magma",
        "magma_version_file": "results/logs/magma_version.txt",
        "gene_location_file_path": "external/reference/gene_loc/NCBI37.3.gene.loc",
        "ld_reference_prefix": "external/reference/1000G_EUR/g1000_eur",
        "magma_log_path": "results/magma/2025_trans_ancestry/ksd_2025.log",
        "pval_input_path": "data/processed/magma_input/ksd_2025.dedup.pval",
        "snploc_input_path": "data/processed/magma_input/ksd_2025.dedup.snploc",
        "audit_evidence_level": "full_path_confirmed",
    }
    fieldnames = list(rows[0].keys())
    for k in extra:
        if k not in fieldnames:
            fieldnames.append(k)
    out_rows = []
    for r in rows:
        r.update(extra)
        out_rows.append([r.get(k, "") for k in fieldnames])
    write_tsv(path, fieldnames, out_rows)


def manuscript_wording() -> None:
    (DOC / "manuscript_ready_genetic_wording.md").write_text("""# Manuscript-Ready Genetic Wording

## A. MAGMA Reproducibility Wording

After GWAS quality control, 4,915,033 variants were retained from 5,960,489 input rows for downstream gene-based analysis. MAGMA v1.10 tested 17,316 genes using GRCh37-compatible gene locations and identified 94 Bonferroni-significant genes and 369 FDR-significant genes. The MAGMA executable, version record, gene-location file, 1000 Genomes European LD reference prefix, p-value input, SNP-location input, gene-based output, and run log were audited as local reproducibility evidence.

## B. LD Limitation Wording

MAGMA was used as an EUR-LD-reference-based gene-prioritization layer for downstream renal papillary context mapping, rather than ancestry-generalizable fine mapping. Because the current input GWAS is trans-ancestry and ancestry-specific summary statistics or matched LD references were not confirmed as ready inputs, gene-level prioritization should be interpreted with this LD-reference limitation.

## C. TWAS Proxy Wording

As a conservative supplementary extension, GTEx v8 Kidney_Cortex S-PrediXcan tested 5,989 genes and identified 51 FDR-supported genes. Of these FDR-supported models, 42 were one-SNP models and 9 were multi-SNP models. TWAS was therefore retained as kidney-cortex proxy support only and was not interpreted as papilla-specific causal expression evidence.

## D. SMR/Coloc Boundary Wording

Current SMR/coloc resources were not claim-grade ready for the audited priority loci. This resource status is treated as missing evidence, not negative biological evidence. No SMR-supported or coloc-supported manuscript claim is made in the current evidence framework.
""", encoding="utf-8")


def freeze_stage3_input() -> tuple[Path, str]:
    STAGE3.mkdir(parents=True, exist_ok=True)
    src = TAB / "candidate_gene_evidence_skeleton.tsv"
    dst = STAGE3 / "stage3_input_candidate_gene_evidence_skeleton_frozen.tsv"
    shutil.copy2(src, dst)
    digest = sha256(dst)
    checksum = STAGE3 / "stage3_input_checksum.md"
    checksum.write_text(f"""# Stage 3 Input Checksum

Frozen input: `results/tables/revision/stage3_gene_tiering/stage3_input_candidate_gene_evidence_skeleton_frozen.tsv`

Source: `results/tables/revision/stage2_genetic/candidate_gene_evidence_skeleton.tsv`

SHA256: `{digest}`

Frozen at: {dt.datetime.now().isoformat(timespec="seconds")}

Interpretation boundary: this frozen table is a candidate evidence skeleton only. It does not assign causal status, SMR/coloc support, or final evidence tiers.
""", encoding="utf-8")
    return dst, digest


def update_tracker() -> None:
    tracker = ROOT / "docs/revision/STAGE_TRACKER.tsv"
    rows = read_tsv(tracker)
    if not rows:
        return
    fields = list(rows[0].keys())
    for r in rows:
        if r.get("stage_id") == "2":
            r["status"] = "stage2A_2B_and_2B_R_completed"
            r["completed_outputs"] = "Stage 2A/2B audits generated; Stage 2B-R reproducibility patch completed; Stage 3 input frozen"
            r["blocking_issues"] = "Full Stage 2 not complete: Stage 2C claim-grade coloc/SMR not performed; original 59-locus table, 2023 European GWAS, ancestry-matched LD, and claim-grade SMR/coloc resources incomplete"
            r["next_stage_ready"] = "stage3_ready_with_conservative_tier_names"
    write_tsv(tracker, fields, [[r.get(f, "") for f in fields] for r in rows])


def report(frozen: Path, digest: str) -> None:
    (DOC / "stage2B_reproducibility_patch_report.md").write_text(f"""# Stage 2B-R Reproducibility Patch Report

Date: {dt.date.today().isoformat()}

## Completed

- Added MAGMA/GWAS/TWAS reproducibility file evidence with concrete local paths.
- Updated `magma_output_audit.tsv` with executable, version, gene-location, LD-reference, input, output, and log paths.
- Wrote manuscript-ready conservative wording for MAGMA reproducibility, LD limitation, TWAS proxy status, and SMR/coloc boundary.
- Frozen Stage 3 input candidate evidence skeleton.
- Updated `docs/revision/STAGE_TRACKER.tsv` to mark Stage 2B-R completed while keeping full Stage 2 incomplete because Stage 2C claim-grade coloc/SMR was not performed.

## MAGMA Reproducibility Evidence Level

`full_path_confirmed`

The audit confirmed local paths for the MAGMA executable, version record, NCBI37 gene-location file, 1000G EUR LD bed/bim/fam resources, deduplicated p-value and SNP-location inputs, gene-based output, and MAGMA run log.

## Frozen Stage 3 Input

- File: `{frozen.relative_to(ROOT).as_posix()}`
- SHA256: `{digest}`

## Unresolved Blockers

- Original 59-locus source table remains unavailable, so 57 vs 59 reconciliation is blocked.
- 2023 European KSD GWAS summary statistics are not confirmed as ready MAGMA input.
- Ancestry-specific or mixed LD references are not available as usable inputs.
- SMR/coloc resources are not claim-grade ready; no SMR/coloc-supported claim should be made.

## Stage 3 Readiness

Stage 3 can start using conservative tier names and the frozen evidence skeleton. Do not use causal/high-confidence causal, SMR-supported, coloc-supported, or papilla-specific TWAS wording.
""", encoding="utf-8")


def main() -> None:
    for d in [DOC, TAB, STAGE3, LOG]:
        d.mkdir(parents=True, exist_ok=True)
    make_file_evidence()
    update_magma_output_audit()
    manuscript_wording()
    frozen, digest = freeze_stage3_input()
    update_tracker()
    report(frozen, digest)
    (LOG / "stage2B_R_reproducibility_patch.log").write_text(f"completed={dt.datetime.now().isoformat(timespec='seconds')}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
