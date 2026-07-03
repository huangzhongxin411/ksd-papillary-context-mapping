#!/usr/bin/env python3
"""Map Phase 1 lead SNPs and loci to nearby GENCODE genes."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


ATTR_RE = re.compile(r'(\S+) "([^"]+)"')


def parse_attrs(text: str) -> dict[str, str]:
    return {k: v for k, v in ATTR_RE.findall(text)}


def load_genes(gtf: Path, biotypes: set[str] | None = None) -> pd.DataFrame:
    opener = gzip.open if gtf.suffix == ".gz" else open
    records = []
    with opener(gtf, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attrs(fields[8])
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype") or ""
            if biotypes and gene_type not in biotypes:
                continue
            chrom = fields[0].removeprefix("chr")
            if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}:
                continue
            records.append(
                {
                    "chr": chrom,
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene": attrs.get("gene_name", attrs.get("gene_id", "")),
                    "gene_type": gene_type,
                }
            )
    genes = pd.DataFrame(records).drop_duplicates(["chr", "start", "end", "gene"])
    return genes.sort_values(["chr", "start", "end"])


def distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    return min(abs(pos - start), abs(pos - end))


def map_genes(
    leads: pd.DataFrame,
    loci: pd.DataFrame,
    genes: pd.DataFrame,
    windows: list[int],
) -> pd.DataFrame:
    records = []
    loci_by_lead = {}
    for _, locus in loci.iterrows():
        for snp in str(locus["lead_snps"]).split(","):
            loci_by_lead[snp] = locus

    for _, lead in leads.iterrows():
        chrom = str(lead["CHR"]).replace("chr", "")
        pos = int(lead["BP"])
        sub = genes[genes["chr"].astype(str) == chrom].copy()
        if sub.empty:
            continue
        sub["distance_to_lead_snp"] = [
            distance_to_interval(pos, int(s), int(e)) for s, e in zip(sub["start"], sub["end"])
        ]
        nearest = sub.sort_values(["distance_to_lead_snp", "start"]).head(1)
        locus = loci_by_lead.get(lead["SNP"])
        locus_id = locus["locus_id"] if locus is not None else ""

        for _, gene in nearest.iterrows():
            records.append(
                {
                    "gene": gene["gene"],
                    "gene_id": gene["gene_id"],
                    "gene_type": gene["gene_type"],
                    "chr": chrom,
                    "start": gene["start"],
                    "end": gene["end"],
                    "nearest_lead_snp": lead["SNP"],
                    "lead_snp_bp": pos,
                    "lead_snp_p": lead["P"],
                    "locus_id": locus_id,
                    "mapping_method": "nearest_gene",
                    "distance_to_lead_snp": int(gene["distance_to_lead_snp"]),
                    "evidence_level": "lead_locus_nearest",
                }
            )

        for window in windows:
            hits = sub[sub["distance_to_lead_snp"] <= window].copy()
            for _, gene in hits.iterrows():
                records.append(
                    {
                        "gene": gene["gene"],
                        "gene_id": gene["gene_id"],
                        "gene_type": gene["gene_type"],
                        "chr": chrom,
                        "start": gene["start"],
                        "end": gene["end"],
                        "nearest_lead_snp": lead["SNP"],
                        "lead_snp_bp": pos,
                        "lead_snp_p": lead["P"],
                        "locus_id": locus_id,
                        "mapping_method": f"lead_snp_within_{window // 1000}kb",
                        "distance_to_lead_snp": int(gene["distance_to_lead_snp"]),
                        "evidence_level": f"within_{window // 1000}kb",
                    }
                )

    out = pd.DataFrame(records)
    if out.empty:
        return out
    priority = {
        "lead_locus_nearest": 0,
        "within_50kb": 1,
        "within_250kb": 2,
        "within_500kb": 3,
    }
    out["evidence_priority"] = out["evidence_level"].map(priority).fillna(9).astype(int)
    out = out.sort_values(
        ["evidence_priority", "lead_snp_p", "distance_to_lead_snp", "gene"]
    )
    out = out.drop_duplicates(["gene", "nearest_lead_snp", "evidence_level"], keep="first")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--leads", required=True, type=Path)
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--windows", default="50000,250000,500000")
    parser.add_argument(
        "--all-biotypes",
        action="store_true",
        help="Keep all GENCODE biotypes; default keeps protein_coding and lncRNA.",
    )
    args = parser.parse_args()

    biotypes = None if args.all_biotypes else {"protein_coding", "lncRNA"}
    windows = [int(x) for x in args.windows.split(",") if x]
    leads = pd.read_csv(args.leads, sep="\t")
    loci = pd.read_csv(args.loci, sep="\t")
    genes = load_genes(args.gtf, biotypes=biotypes)
    mapped = map_genes(leads, loci, genes, windows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(args.out, sep="\t", index=False)
    print(f"genes_loaded\t{len(genes)}")
    print(f"candidate_rows\t{len(mapped)}")
    if not mapped.empty:
        print(f"candidate_genes\t{mapped['gene'].nunique()}")
    print(f"wrote\t{args.out}")


if __name__ == "__main__":
    main()

