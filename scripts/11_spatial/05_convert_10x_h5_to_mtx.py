#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
import traceback

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, io


SAMPLES = ["GSM6250307", "GSM6250308", "GSM6250309", "GSM6250310", "GSM7166170"]
BASE = Path("data/processed/spatial/gse206306")
OUT_QC = Path("results/spatial/phase28b_loading_rescue/h5_to_mtx_conversion_qc_v0.1.tsv")


def decode_array(values):
    return [v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v) for v in values]


def convert_one(sample_id: str, force: bool = False):
    sample_dir = BASE / sample_id
    h5_path = sample_dir / "filtered_feature_bc_matrix.h5"
    out_dir = sample_dir / "filtered_feature_bc_matrix"
    matrix_out = out_dir / "matrix.mtx.gz"
    barcodes_out = out_dir / "barcodes.tsv.gz"
    features_out = out_dir / "features.tsv.gz"
    row = {
        "sample_id": sample_id,
        "input_h5": str(h5_path),
        "h5_exists": h5_path.exists(),
        "conversion_attempted": False,
        "conversion_success": False,
        "n_genes": "",
        "n_barcodes": "",
        "n_nonzero": "",
        "output_matrix": str(matrix_out),
        "output_barcodes": str(barcodes_out),
        "output_features": str(features_out),
        "error_message": "",
        "notes": "",
    }
    if not h5_path.exists():
        row["error_message"] = "input_h5_missing"
        return row
    if matrix_out.exists() and barcodes_out.exists() and features_out.exists() and not force:
        row["conversion_success"] = True
        row["notes"] = "Existing mtx outputs detected; skipped conversion because --force was not supplied."
        return row
    row["conversion_attempted"] = True
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        with h5py.File(h5_path, "r") as h5:
            grp = h5["matrix"]
            data = grp["data"][:]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]
            shape = tuple(grp["shape"][:])
            barcodes = decode_array(grp["barcodes"][:])
            feat = grp["features"]
            gene_ids = decode_array(feat["id"][:])
            gene_names = decode_array(feat["name"][:])
            if "feature_type" in feat:
                feature_type = decode_array(feat["feature_type"][:])
            else:
                feature_type = ["Gene Expression"] * len(gene_ids)

        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
        n_genes = len(gene_ids)
        n_barcodes = len(barcodes)
        if mat.shape == (n_barcodes, n_genes):
            mat = mat.T
        if mat.shape != (n_genes, n_barcodes):
            raise ValueError(f"Unexpected matrix shape {mat.shape}; features={n_genes}, barcodes={n_barcodes}")

        with gzip.open(matrix_out, "wb") as fh:
            io.mmwrite(fh, mat)
        with gzip.open(barcodes_out, "wt") as fh:
            for bc in barcodes:
                fh.write(f"{bc}\n")
        features = pd.DataFrame({
            "gene_id": gene_ids,
            "gene_name": gene_names,
            "feature_type": feature_type,
        })
        features.to_csv(features_out, sep="\t", index=False, header=False, compression="gzip")

        row.update({
            "conversion_success": True,
            "n_genes": n_genes,
            "n_barcodes": n_barcodes,
            "n_nonzero": int(mat.nnz),
            "notes": "Converted 10x HDF5 matrix to Matrix Market fallback format.",
        })
    except Exception as exc:
        row["error_message"] = f"{type(exc).__name__}: {exc}"
        row["notes"] = traceback.format_exc(limit=1).replace("\n", " | ")
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="overwrite existing mtx fallback files")
    args = parser.parse_args()
    OUT_QC.parent.mkdir(parents=True, exist_ok=True)
    rows = [convert_one(s, force=args.force) for s in SAMPLES]
    pd.DataFrame(rows).to_csv(OUT_QC, sep="\t", index=False)
    print(f"wrote {OUT_QC}")


if __name__ == "__main__":
    main()
