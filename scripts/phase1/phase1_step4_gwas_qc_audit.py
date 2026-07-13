#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import random
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/gwas/2025_trans_ancestry/meta_sumstats"
CLEAN = ROOT / "data/processed/gwas/2025_trans_ancestry/meta_sumstats.cleaned.tsv.gz"
ALIASES = ROOT / "config/gwas_column_aliases.tsv"
TABLE_DIR = ROOT / "results/tables"
NOTE_DIR = ROOT / "notes"
SOURCE_DIR = ROOT / "source_data/figures"
EXPECTED_RAW = 5_960_489
EXPECTED_CLEAN = 4_915_033
PLOT_SAMPLE_N = 200_000


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore", newline="")
    return path.open("r", encoding="utf-8", errors="ignore", newline="")


def read_aliases() -> dict[str, list[str]]:
    out = {}
    with ALIASES.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            std = row["standard"]
            out[std] = [x.strip() for x in row["aliases"].split(",") if x.strip()]
    return out


def map_columns(header: list[str], aliases: dict[str, list[str]]) -> dict[str, str]:
    lower = {h.lower(): h for h in header}
    mapping = {}
    for std, names in aliases.items():
        hit = ""
        for name in [std] + names:
            if name.lower() in lower:
                hit = lower[name.lower()]
                break
        mapping[std] = hit
    return mapping


def parse_float(value: str) -> float | None:
    try:
        if value is None or value == "":
            return None
        x = float(value)
        if math.isnan(x) or math.isinf(x):
            return None
        return x
    except Exception:
        return None


def is_acgt(a: str) -> bool:
    return a.upper() in {"A", "C", "G", "T"}


def is_ambiguous(a1: str, a2: str) -> bool:
    pair = {a1.upper(), a2.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def safe_int_chr(value: str) -> int | None:
    try:
        v = str(value).strip().replace("chr", "")
        if v.upper() == "X":
            return 23
        if v.upper() == "Y":
            return 24
        return int(float(v))
    except Exception:
        return None


def audit_file(path: Path, label: str, aliases: dict[str, list[str]], make_plot_source: bool = False) -> tuple[dict[str, object], dict[str, str]]:
    rng = random.Random(20260709)
    reservoir: list[dict[str, object]] = []
    gws_points: list[dict[str, object]] = []
    seen: set[str] = set()
    duplicated = 0
    total = 0
    missing_p = invalid_p = missing_chr_pos = 0
    non_acgt = ambiguous = maf_lt_001 = info_lt_08 = 0
    min_p = 1.0
    gws = 0
    chr_counter: Counter[str] = Counter()

    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        mapping = map_columns(header, aliases)
        for row in reader:
            total += 1
            snp = row.get(mapping.get("SNP", ""), "")
            if snp:
                if snp in seen:
                    duplicated += 1
                else:
                    seen.add(snp)
            p = parse_float(row.get(mapping.get("P", ""), ""))
            if p is None:
                missing_p += 1
            elif p <= 0 or p > 1:
                invalid_p += 1
            else:
                min_p = min(min_p, p)
                if p < 5e-8:
                    gws += 1
            chrom = row.get(mapping.get("CHR", ""), "")
            pos = row.get(mapping.get("BP", ""), "")
            chr_i = safe_int_chr(chrom)
            pos_f = parse_float(pos)
            if chr_i is None or pos_f is None:
                missing_chr_pos += 1
            else:
                chr_counter[str(chr_i)] += 1
            ea = row.get(mapping.get("EA", ""), "")
            nea = row.get(mapping.get("NEA", ""), "")
            if ea and nea:
                if not (is_acgt(ea) and is_acgt(nea)):
                    non_acgt += 1
                elif is_ambiguous(ea, nea):
                    ambiguous += 1
            eaf = parse_float(row.get(mapping.get("EAF", ""), ""))
            if eaf is not None:
                maf = min(eaf, 1 - eaf)
                if maf < 0.01:
                    maf_lt_001 += 1
            info = parse_float(row.get(mapping.get("INFO", ""), ""))
            if info is not None and info < 0.8:
                info_lt_08 += 1

            if make_plot_source and p is not None and p > 0 and p <= 1 and chr_i is not None and pos_f is not None:
                point = {
                    "SNP": snp,
                    "CHR": chr_i,
                    "BP": int(pos_f),
                    "P": p,
                    "neglog10p": -math.log10(p),
                    "is_gws": p < 5e-8,
                }
                if p < 5e-8:
                    gws_points.append(point)
                elif len(reservoir) < PLOT_SAMPLE_N:
                    reservoir.append(point)
                else:
                    j = rng.randrange(total)
                    if j < PLOT_SAMPLE_N:
                        reservoir[j] = point

    stats = {
        "file_label": label,
        "file_path": rel(path),
        "total_rows": total,
        "expected_rows": EXPECTED_RAW if label == "raw_gwas_sumstats" else EXPECTED_CLEAN,
        "matches_expected_rows": "yes" if total == (EXPECTED_RAW if label == "raw_gwas_sumstats" else EXPECTED_CLEAN) else "no",
        "missing_p_values": missing_p,
        "invalid_p_values_outside_0_1": invalid_p,
        "missing_chr_or_position": missing_chr_pos,
        "duplicated_variant_ids": duplicated,
        "ambiguous_at_cg_variants": ambiguous,
        "non_acgt_alleles": non_acgt,
        "maf_lt_0_01": maf_lt_001 if mapping.get("EAF") else "not_available",
        "info_lt_0_8": info_lt_08 if mapping.get("INFO") else "not_available",
        "genomewide_significant_variants_p_lt_5e_8": gws,
        "min_p": f"{min_p:.6g}",
        "chromosomes_observed": ";".join(sorted(chr_counter.keys(), key=lambda x: int(x))),
        "notes": "Streaming audit only; no GWAS re-analysis beyond requested QC counts.",
    }
    if make_plot_source:
        SOURCE_DIR.mkdir(parents=True, exist_ok=True)
        points = gws_points + reservoir
        points.sort(key=lambda x: (int(x["CHR"]), int(x["BP"]), float(x["P"])))
        with gzip.open(SOURCE_DIR / "phase1_step4_gwas_plot_source_data.tsv.gz", "wt", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=["SNP", "CHR", "BP", "P", "neglog10p", "is_gws"], delimiter="\t")
            writer.writeheader()
            writer.writerows(points)
        stats["plot_source_rows"] = len(points)
        stats["plot_source_file"] = rel(SOURCE_DIR / "phase1_step4_gwas_plot_source_data.tsv.gz")
    return stats, mapping


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    aliases = read_aliases()
    raw_stats, raw_map = audit_file(RAW, "raw_gwas_sumstats", aliases, make_plot_source=False)
    clean_stats, clean_map = audit_file(CLEAN, "cleaned_gwas_sumstats", aliases, make_plot_source=True)

    summary_fields = [
        "file_label",
        "file_path",
        "total_rows",
        "expected_rows",
        "matches_expected_rows",
        "missing_p_values",
        "invalid_p_values_outside_0_1",
        "missing_chr_or_position",
        "duplicated_variant_ids",
        "ambiguous_at_cg_variants",
        "non_acgt_alleles",
        "maf_lt_0_01",
        "info_lt_0_8",
        "genomewide_significant_variants_p_lt_5e_8",
        "min_p",
        "chromosomes_observed",
        "plot_source_rows",
        "plot_source_file",
        "notes",
    ]
    write_tsv(TABLE_DIR / "phase1_step4_gwas_qc_summary.tsv", [raw_stats, clean_stats], summary_fields)

    map_rows = []
    for std in aliases:
        map_rows.append(
            {
                "standard_field": std,
                "raw_column": raw_map.get(std, ""),
                "cleaned_column": clean_map.get(std, ""),
                "aliases_considered": ",".join(aliases[std]),
                "required_for_phase1": "true" if std in {"SNP", "CHR", "BP", "EA", "NEA", "P"} else "false",
                "notes": "INFO absent in current raw/cleaned GWAS." if std == "INFO" and not clean_map.get(std) else "",
            }
        )
    write_tsv(TABLE_DIR / "phase1_step4_gwas_column_mapping.tsv", map_rows, ["standard_field", "raw_column", "cleaned_column", "aliases_considered", "required_for_phase1", "notes"])

    NOTE_DIR.mkdir(parents=True, exist_ok=True)
    report = [
        "# Phase 1-Step 4 GWAS QC Audit",
        "",
        f"Raw GWAS file: `{rel(RAW)}`",
        f"Cleaned GWAS file: `{rel(CLEAN)}`",
        f"Column alias file: `{rel(ALIASES)}`",
        "",
        "## Row Counts",
        "",
        f"- Raw rows: {raw_stats['total_rows']} (expected manuscript denominator: {EXPECTED_RAW}; match: {raw_stats['matches_expected_rows']})",
        f"- Cleaned rows: {clean_stats['total_rows']} (expected retained rows: {EXPECTED_CLEAN}; match: {clean_stats['matches_expected_rows']})",
        "",
        "## Key QC Counts",
        "",
        f"- Cleaned genome-wide significant variants at P < 5e-8: {clean_stats['genomewide_significant_variants_p_lt_5e_8']}",
        f"- Cleaned ambiguous A/T or C/G variants: {clean_stats['ambiguous_at_cg_variants']}",
        f"- Cleaned non-A/C/G/T allele rows: {clean_stats['non_acgt_alleles']}",
        f"- Cleaned MAF < 0.01 rows from EAF: {clean_stats['maf_lt_0_01']}",
        f"- INFO < 0.8: {clean_stats['info_lt_0_8']}",
        "",
        "The audit used streaming row-by-row inspection and did not load the full GWAS into memory. QC plots should be interpreted as diagnostic visualizations only.",
    ]
    (NOTE_DIR / "phase1_step4_gwas_qc_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
