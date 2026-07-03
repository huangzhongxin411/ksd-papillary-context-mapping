#!/usr/bin/env python3
import gzip
import sqlite3
import sys

model_db, gwas_in, gwas_out, qc_out = sys.argv[1:5]

con = sqlite3.connect(model_db)
cur = con.cursor()
by_rsid = {}
for rsid, varid, ref, eff in cur.execute(
    """
    select distinct rsid, varID, ref_allele, eff_allele
    from weights
    where rsid is not null and rsid != ''
      and varID is not null and varID != ''
      and ref_allele is not null and eff_allele is not null
    """
):
    by_rsid.setdefault(rsid, []).append((varid, ref, eff))

stats = {
    "gwas_rows": 0,
    "rsid_overlap": 0,
    "allele_compatible": 0,
    "written_rows": 0,
    "flipped_rows": 0,
    "ambiguous_multi_mapping_skipped": 0,
}

seen_varids = set()
with gzip.open(gwas_in, "rt") as fin, gzip.open(gwas_out, "wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    out_header = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "Z", "P", "N", "EAF", "SOURCE_RSID"]
    fout.write("\t".join(out_header) + "\n")
    for line in fin:
        stats["gwas_rows"] += 1
        parts = line.rstrip("\n").split("\t")
        rsid = parts[idx["SNP"]]
        candidates = by_rsid.get(rsid)
        if not candidates:
            continue
        stats["rsid_overlap"] += 1
        a1 = parts[idx["A1"]].upper()
        a2 = parts[idx["A2"]].upper()
        compatible = []
        for varid, ref, eff in candidates:
            ref = ref.upper()
            eff = eff.upper()
            if {a1, a2} == {ref, eff}:
                compatible.append((varid, ref, eff))
        if not compatible:
            continue
        stats["allele_compatible"] += 1
        # Prefer the first compatible unique varID; skip if the same rsid maps to multiple compatible varIDs.
        unique = []
        for item in compatible:
            if item[0] not in [u[0] for u in unique]:
                unique.append(item)
        if len(unique) != 1:
            stats["ambiguous_multi_mapping_skipped"] += 1
            continue
        varid, ref, eff = unique[0]
        if varid in seen_varids:
            continue
        seen_varids.add(varid)
        beta = float(parts[idx["BETA"]])
        z = float(parts[idx["Z"]])
        eaf_raw = parts[idx["EAF"]]
        eaf = None
        try:
            eaf = float(eaf_raw)
        except Exception:
            pass
        flipped = a1 != eff
        if flipped:
            beta = -beta
            z = -z
            stats["flipped_rows"] += 1
            if eaf is not None:
                eaf = 1.0 - eaf
        out = [
            varid,
            parts[idx["CHR"]],
            parts[idx["BP"]],
            eff,
            ref,
            f"{beta:.12g}",
            parts[idx["SE"]],
            f"{z:.12g}",
            parts[idx["P"]],
            parts[idx["N"]],
            f"{eaf:.12g}" if eaf is not None else eaf_raw,
            rsid,
        ]
        fout.write("\t".join(out) + "\n")
        stats["written_rows"] += 1

with open(qc_out, "w") as out:
    out.write("metric\tvalue\n")
    for k, v in stats.items():
        out.write(f"{k}\t{v}\n")

print(f"wrote {stats['written_rows']} mapped rows to {gwas_out}")
