#!/usr/bin/env python3
"""Create compact Phase 1 gene sets for single-nucleus projection."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--sizes", default="50,100")
    args = parser.parse_args()

    x = pd.read_csv(args.candidates, sep="\t")
    x = x[x["gene_type"].eq("protein_coding")].copy()
    x = x.sort_values(["evidence_priority", "lead_snp_p", "distance_to_lead_snp", "gene"])
    genes = x.drop_duplicates("gene", keep="first").copy()
    genes.insert(0, "rank_in_candidate_list", range(1, len(genes) + 1))

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    ranked_path = Path(f"{args.out_prefix}_ranked_protein_coding.tsv")
    genes.to_csv(ranked_path, sep="\t", index=False)
    print(f"ranked_genes\t{len(genes)}\t{ranked_path}")

    for size in [int(s) for s in args.sizes.split(",") if s]:
        sub = genes.head(size)
        out = Path(f"{args.out_prefix}_top{size}.txt")
        out.write_text("\n".join(sub["gene"].astype(str)) + "\n")
        print(f"top{size}\t{len(sub)}\t{out}")


if __name__ == "__main__":
    main()

