#!/usr/bin/env python3
import sqlite3
import sys

db_path = sys.argv[1]
out_path = sys.argv[2]

con = sqlite3.connect(db_path)
cur = con.cursor()
rows = cur.execute(
    """
    select distinct rsid, ref_allele, eff_allele, varID, ref_allele, eff_allele
    from weights
    where rsid is not null and rsid != ''
      and varID is not null and varID != ''
      and ref_allele is not null and eff_allele is not null
    """
).fetchall()

with open(out_path, "w") as out:
    out.write("rsid\ta0\ta1\tpanel_variant_id\tpanel_variant_a0\tpanel_variant_a1\tswap\n")
    seen = set()
    for rsid, a0, a1, panel_id, panel_a0, panel_a1 in rows:
        key = (rsid, panel_id, a0, a1)
        if key in seen:
            continue
        seen.add(key)
        out.write(f"{rsid}\t{a0}\t{a1}\t{panel_id}\t{panel_a0}\t{panel_a1}\t1\n")

print(f"wrote {len(seen)} rows to {out_path}")
