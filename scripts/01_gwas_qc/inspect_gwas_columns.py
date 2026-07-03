#!/usr/bin/env python3
"""Inspect GWAS summary statistic columns before QC.

This script is intentionally lightweight: it reads only the first rows and
prints column names plus suggested mappings from config/gwas_column_aliases.tsv.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def sniff_delimiter(sample: str) -> str:
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t, ;").delimiter
    except csv.Error:
        return "\t"


def load_aliases(path: Path) -> dict[str, list[str]]:
    aliases: dict[str, list[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            aliases[row["standard"]] = [x.strip() for x in row["aliases"].split(",")]
    return aliases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_file", type=Path, nargs="?")
    parser.add_argument("--input", dest="input_file", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--aliases",
        type=Path,
        default=Path("config/gwas_column_aliases.tsv"),
    )
    args = parser.parse_args()
    gwas_file = args.input_file or args.gwas_file
    if gwas_file is None:
        raise SystemExit("Provide a GWAS file as positional argument or --input")

    with open_text(gwas_file) as fh:
        sample = fh.read(65536)
    delimiter = sniff_delimiter(sample)
    header = sample.splitlines()[0].split(delimiter)

    aliases = load_aliases(args.aliases)
    lowered = {col.lower(): col for col in header}

    lines: list[str] = []
    lines.append(f"file\t{gwas_file}")
    lines.append(f"detected_delimiter\t{repr(delimiter)}")
    lines.append("columns")
    for col in header:
        lines.append(f"- {col}")

    lines.append("")
    lines.append("suggested_mapping")
    for standard, names in aliases.items():
        hit = ""
        for name in names:
            if name.lower() in lowered:
                hit = lowered[name.lower()]
                break
        lines.append(f"{standard}\t{hit or 'NOT_FOUND'}")

    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)


if __name__ == "__main__":
    main()
