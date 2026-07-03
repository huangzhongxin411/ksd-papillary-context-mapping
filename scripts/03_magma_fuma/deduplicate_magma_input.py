#!/usr/bin/env python3
"""Deduplicate MAGMA snploc/pval inputs by retaining the smallest P per SNP."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_pval(path: Path) -> tuple[list[str], dict[str, list[str]], dict[str, int]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if "SNP" not in header or "P" not in header:
            raise SystemExit(f"{path} must contain SNP and P columns")
        snp_i = header.index("SNP")
        p_i = header.index("P")
        best: dict[str, list[str]] = {}
        seen: dict[str, int] = {}
        for row in reader:
            if not row or len(row) <= max(snp_i, p_i):
                continue
            snp = row[snp_i]
            seen[snp] = seen.get(snp, 0) + 1
            try:
                p = float(row[p_i])
            except ValueError:
                p = math.inf
            if snp not in best:
                best[snp] = row
                continue
            try:
                old_p = float(best[snp][p_i])
            except ValueError:
                old_p = math.inf
            if p < old_p:
                best[snp] = row
    return header, best, seen


def read_snploc(path: Path, keep: set[str]) -> tuple[list[str], dict[str, list[str]], dict[str, int], bool]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader)
        has_header = len(first) >= 3 and first[0].upper() == "SNP"
        rows = reader if has_header else iter([first, *reader])
        best: dict[str, list[str]] = {}
        seen: dict[str, int] = {}
        for row in rows:
            if len(row) < 3:
                continue
            snp = row[0]
            seen[snp] = seen.get(snp, 0) + 1
            if snp in keep and snp not in best:
                best[snp] = row[:3]
    return ["SNP", "CHR", "BP"], best, seen, has_header


def write_rows(path: Path, header: list[str], rows: list[list[str]], include_header: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if include_header:
            writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snploc", required=True, type=Path)
    parser.add_argument("--pval", required=True, type=Path)
    parser.add_argument("--out-snploc", required=True, type=Path)
    parser.add_argument("--out-pval", required=True, type=Path)
    parser.add_argument("--qc-out", type=Path, default=Path("results/tables/magma_input_qc.tsv"))
    args = parser.parse_args()

    pval_header, pval_best, pval_seen = read_pval(args.pval)
    snploc_header, snploc_best, snploc_seen, snploc_had_header = read_snploc(args.snploc, set(pval_best))
    shared_snps = sorted(set(pval_best).intersection(snploc_best))

    write_rows(args.out_snploc, snploc_header, [snploc_best[s] for s in shared_snps], snploc_had_header)
    write_rows(args.out_pval, pval_header, [pval_best[s] for s in shared_snps], include_header=True)

    args.qc_out.parent.mkdir(parents=True, exist_ok=True)
    metrics = [
        ("input_snploc_rows", sum(snploc_seen.values())),
        ("input_pval_rows", sum(pval_seen.values())),
        ("input_snploc_unique_snps", len(snploc_seen)),
        ("input_pval_unique_snps", len(pval_seen)),
        ("input_snploc_duplicate_snps", sum(v > 1 for v in snploc_seen.values())),
        ("input_pval_duplicate_snps", sum(v > 1 for v in pval_seen.values())),
        ("output_shared_unique_snps", len(shared_snps)),
        ("output_snploc_has_header", str(snploc_had_header).lower()),
    ]
    with args.qc_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)


if __name__ == "__main__":
    main()
