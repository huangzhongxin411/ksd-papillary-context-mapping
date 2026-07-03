#!/usr/bin/env python3
"""Phase 1 GWAS summary-statistics QC.

The script expects a TSV/CSV file and explicit column mappings. It writes a
cleaned TSV.GZ and a compact QC report. It also generates QQ and Manhattan
plots when matplotlib is available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd


AMBIGUOUS = {"AT", "TA", "CG", "GC"}


def read_table(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    opener = gzip.open if compression == "gzip" else open
    with opener(path, "rt", newline="") as fh:
        sample = fh.read(65536)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters="\t, ;").delimiter
    except csv.Error:
        sep = "\t"
    return pd.read_csv(path, sep=sep, compression=compression, low_memory=False)


def require_columns(df: pd.DataFrame, mapping: dict[str, str]) -> None:
    missing = [std for std, col in mapping.items() if col and col not in df.columns]
    if missing:
        raise SystemExit(f"Missing mapped columns: {missing}")


def lambda_gc(pvalues: pd.Series) -> float:
    try:
        from scipy.stats import chi2
        p = pvalues.dropna().clip(lower=np.nextafter(0, 1), upper=1)
        chisq = chi2.isf(p, 1)
        return float(np.nanmedian(chisq) / 0.454936423119572)
    except Exception:
        from statistics import NormalDist

        p = pvalues.dropna().clip(lower=np.nextafter(0, 1), upper=1).to_numpy()
        nd = NormalDist()
        upper = np.nextafter(1.0, 0.0)
        chisq = np.array([nd.inv_cdf(min(upper, 1 - float(x) / 2)) ** 2 for x in p])
    return float(np.nanmedian(chisq) / 0.454936423119572)


def write_plots(df: pd.DataFrame, out_prefix: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    p = df["P"].dropna().sort_values()
    exp = -np.log10(np.arange(1, len(p) + 1) / (len(p) + 1))
    obs = -np.log10(p.to_numpy())
    plt.figure(figsize=(4.5, 4.5))
    plt.scatter(exp, obs, s=2, alpha=0.5)
    lim = max(float(np.nanmax(exp)), float(np.nanmax(obs)))
    plt.plot([0, lim], [0, lim], color="black", linewidth=0.8)
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.qq_plot.pdf")
    plt.close()

    man = df[["CHR", "BP", "P"]].copy()
    man = man.dropna()
    man["CHR"] = man["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
    man = man[man["CHR"].str.match(r"^\d+$")]
    man["CHR"] = man["CHR"].astype(int)
    man = man.sort_values(["CHR", "BP"])
    offsets = {}
    offset = 0
    ticks = []
    labels = []
    xs = []
    for chrom, sub in man.groupby("CHR", sort=True):
        offsets[chrom] = offset
        x = sub["BP"].to_numpy() + offset
        xs.append(pd.Series(x, index=sub.index))
        ticks.append(float(np.nanmedian(x)))
        labels.append(str(chrom))
        offset += int(sub["BP"].max()) + 1_000_000
    if xs:
        man["x"] = pd.concat(xs).sort_index()
        man["mlog10p"] = -np.log10(man["P"].clip(lower=np.nextafter(0, 1)))
        plt.figure(figsize=(10, 3.8))
        colors = np.where(man["CHR"] % 2 == 0, "#4c78a8", "#222222")
        plt.scatter(man["x"], man["mlog10p"], c=colors, s=2, alpha=0.65)
        plt.axhline(-math.log10(5e-8), color="#d62728", linewidth=0.8)
        plt.xticks(ticks, labels, fontsize=7)
        plt.xlabel("Chromosome")
        plt.ylabel("-log10(P)")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.manhattan_plot.pdf")
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out-prefix", type=Path)
    parser.add_argument("--out", type=Path, help="Cleaned TSV.GZ output path.")
    parser.add_argument("--report", type=Path, help="QC report TSV output path.")
    parser.add_argument("--prefix", type=Path, help="Figure prefix; defaults to --out-prefix or --out stem.")
    parser.add_argument("--snp", required=True)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--bp", required=True)
    parser.add_argument("--ea", required=True)
    parser.add_argument("--nea", required=True)
    parser.add_argument("--p", required=True)
    parser.add_argument("--beta")
    parser.add_argument("--or-col")
    parser.add_argument("--se")
    parser.add_argument("--eaf")
    parser.add_argument("--n")
    parser.add_argument("--info")
    parser.add_argument("--info-min", type=float, default=0.8)
    parser.add_argument("--maf-min", type=float, default=0.01)
    args = parser.parse_args()
    if args.out_prefix is None and args.out is None:
        raise SystemExit("Provide --out-prefix or --out")

    raw = read_table(args.input)
    mapping = {
        "SNP": args.snp,
        "CHR": args.chr,
        "BP": args.bp,
        "EA": args.ea,
        "NEA": args.nea,
        "P": args.p,
        "BETA": args.beta,
        "OR": args.or_col,
        "SE": args.se,
        "EAF": args.eaf,
        "N": args.n,
        "INFO": args.info,
    }
    require_columns(raw, {k: v for k, v in mapping.items() if v})

    df = pd.DataFrame()
    for std, col in mapping.items():
        if col:
            df[std] = raw[col]

    report: list[tuple[str, int]] = [("raw_rows", len(df))]

    required = ["SNP", "CHR", "BP", "EA", "NEA", "P"]
    before = len(df)
    df = df.dropna(subset=required)
    report.append(("removed_missing_required", before - len(df)))

    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    before = len(df)
    df = df[df["P"].gt(0) & df["P"].le(1)]
    report.append(("removed_invalid_p", before - len(df)))

    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    before = len(df)
    df = df.dropna(subset=["BP"])
    report.append(("removed_invalid_bp", before - len(df)))
    df["BP"] = df["BP"].astype(int)

    for col in ["EA", "NEA"]:
        df[col] = df[col].astype(str).str.upper().str.strip()
    before = len(df)
    df = df[df["EA"].str.fullmatch("[ACGT]") & df["NEA"].str.fullmatch("[ACGT]")]
    report.append(("removed_non_snv_or_bad_allele", before - len(df)))

    before = len(df)
    df = df.drop_duplicates(subset=["SNP"], keep="first")
    report.append(("removed_duplicate_snp", before - len(df)))

    pair = df["EA"] + df["NEA"]
    before = len(df)
    df = df[~pair.isin(AMBIGUOUS)]
    report.append(("removed_ambiguous_at_cg", before - len(df)))

    if "INFO" in df.columns:
        df["INFO"] = pd.to_numeric(df["INFO"], errors="coerce")
        before = len(df)
        df = df[df["INFO"].isna() | df["INFO"].ge(args.info_min)]
        report.append((f"removed_info_lt_{args.info_min}", before - len(df)))

    if "EAF" in df.columns:
        df["EAF"] = pd.to_numeric(df["EAF"], errors="coerce")
        maf = np.minimum(df["EAF"], 1 - df["EAF"])
        before = len(df)
        df = df[maf.isna() | maf.ge(args.maf_min)]
        report.append((f"removed_maf_lt_{args.maf_min}", before - len(df)))

    report.append(("qc_pass_rows", len(df)))
    report.append(("genome_wide_significant_p_lt_5e-8", int(df["P"].lt(5e-8).sum())))
    report.append(("lambda_gc", lambda_gc(df["P"])))

    if args.out is not None:
        clean_path = args.out
    else:
        clean_path = Path(f"{args.out_prefix}.clean.tsv.gz")
    clean_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(clean_path, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)

    if args.report is not None:
        report_path = args.report
    elif args.out_prefix is not None:
        report_path = Path(f"{args.out_prefix}.qc_report.tsv")
    else:
        report_path = clean_path.with_suffix("").with_suffix(".qc_report.tsv")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(report, columns=["metric", "value"]).to_csv(report_path, sep="\t", index=False)
    fig_prefix = args.prefix or args.out_prefix or clean_path.with_suffix("").with_suffix("")
    fig_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_plots(df, fig_prefix)

    print(f"wrote\t{clean_path}")
    print(f"wrote\t{report_path}")


if __name__ == "__main__":
    main()
