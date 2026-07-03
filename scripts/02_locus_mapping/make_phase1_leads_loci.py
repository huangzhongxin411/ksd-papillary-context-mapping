#!/usr/bin/env python3
"""Create exploratory Phase 1 lead SNPs and loci from cleaned GWAS.

This is distance-based pruning for the minimal loop. It is not a replacement
for LD clumping; the output is explicitly marked `distance_pruned`.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep="\t", compression=compression)


def prune_leads(df: pd.DataFrame, p_threshold: float, window_bp: int) -> pd.DataFrame:
    sig = df[df["P"] < p_threshold].copy()
    sig["CHR"] = sig["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
    sig["BP"] = pd.to_numeric(sig["BP"], errors="coerce")
    sig = sig.dropna(subset=["CHR", "BP", "P"]).sort_values("P")

    selected = []
    occupied: dict[str, list[tuple[int, int]]] = {}
    for _, row in sig.iterrows():
        chrom = str(row["CHR"])
        bp = int(row["BP"])
        blocked = any(start <= bp <= end for start, end in occupied.get(chrom, []))
        if blocked:
            continue
        selected.append(row)
        occupied.setdefault(chrom, []).append((bp - window_bp, bp + window_bp))

    if not selected:
        return sig.iloc[0:0].copy()
    leads = pd.DataFrame(selected)
    leads = leads.sort_values(["CHR", "BP", "P"])
    leads["lead_method"] = "distance_pruned"
    return leads


def make_loci(leads: pd.DataFrame, locus_window_bp: int) -> pd.DataFrame:
    records = []
    for chrom, sub in leads.groupby("CHR", sort=True):
        sub = sub.sort_values("BP")
        current = None
        for _, row in sub.iterrows():
            start = int(row["BP"]) - locus_window_bp
            end = int(row["BP"]) + locus_window_bp
            if current is None or start > current["end"]:
                if current is not None:
                    records.append(current)
                current = {
                    "CHR": chrom,
                    "start": start,
                    "end": end,
                    "lead_snps": [row["SNP"]],
                    "min_p": float(row["P"]),
                    "top_lead_snp": row["SNP"],
                }
            else:
                current["end"] = max(current["end"], end)
                current["lead_snps"].append(row["SNP"])
                if float(row["P"]) < current["min_p"]:
                    current["min_p"] = float(row["P"])
                    current["top_lead_snp"] = row["SNP"]
        if current is not None:
            records.append(current)

    loci = pd.DataFrame(records)
    if loci.empty:
        return loci
    loci.insert(0, "locus_id", [f"KSD_L{idx:03d}" for idx in range(1, len(loci) + 1)])
    loci["lead_snps"] = loci["lead_snps"].map(lambda xs: ",".join(xs))
    loci["mapping_stage"] = "phase1_distance_locus"
    return loci


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean-gwas", type=Path)
    parser.add_argument("--input", dest="input_file", type=Path)
    parser.add_argument("--out-leads", type=Path)
    parser.add_argument("--out-loci", type=Path)
    parser.add_argument("--out-prefix", type=Path)
    parser.add_argument("--p-threshold", "--p", type=float, default=5e-8)
    parser.add_argument("--lead-window-bp", "--window", type=int, default=1_000_000)
    parser.add_argument("--locus-window-bp", type=int, default=1_000_000)
    args = parser.parse_args()
    clean_gwas = args.clean_gwas or args.input_file
    if clean_gwas is None:
        raise SystemExit("Provide --clean-gwas or --input")
    if args.out_prefix is not None:
        out_leads = Path(f"{args.out_prefix}_lead_snps.tsv")
        out_loci = Path(f"{args.out_prefix}_loci.tsv")
    else:
        if args.out_leads is None or args.out_loci is None:
            raise SystemExit("Provide --out-prefix or both --out-leads and --out-loci")
        out_leads = args.out_leads
        out_loci = args.out_loci

    df = read_table(clean_gwas)
    required = {"SNP", "CHR", "BP", "P"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in cleaned GWAS: {sorted(missing)}")

    leads = prune_leads(df, args.p_threshold, args.lead_window_bp)
    loci = make_loci(leads, args.locus_window_bp)

    out_leads.parent.mkdir(parents=True, exist_ok=True)
    out_loci.parent.mkdir(parents=True, exist_ok=True)
    leads.to_csv(out_leads, sep="\t", index=False)
    loci.to_csv(out_loci, sep="\t", index=False)
    print(f"lead_snps\t{len(leads)}\t{out_leads}")
    print(f"loci\t{len(loci)}\t{out_loci}")


if __name__ == "__main__":
    main()
