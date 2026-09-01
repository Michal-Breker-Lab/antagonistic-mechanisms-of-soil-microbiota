#!/usr/bin/env python
"""Build tables/MF6_ani_raw.tsv from the sparse all-vs-all ANI table.

skani was run with MF6 as the REFERENCE and every genome as the QUERY, but the
figures (fig10, fig12) key on Ref_file being the *genome* and Query_file being
MF6 -- the orientation of the retained original. So the two file columns, the
two name columns and the two aligned-fraction columns are all swapped here.
Swapping the fractions matters: they are not symmetric (MF7 is ref 99.25 /
query 89.07 against MF6, an asymmetry that is the ~940 kb duplication).
"""
import csv
from pathlib import Path

D = Path(__file__).resolve().parent.parent
SRC = D / "results" / "ani_sparse.tsv"
OUT = D / "tables" / "MF6_ani_raw.tsv"
COLS = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref",
        "Align_fraction_query", "Ref_name", "Query_name"]

out = []
with open(SRC, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if Path(r["Ref_file"]).name != "MF6.fna":
            continue
        out.append({"Ref_file": "genomes/" + Path(r["Query_file"]).name,
                    "Query_file": "MF6.fna",
                    "ANI": r["ANI"],
                    "Align_fraction_ref": r["Align_fraction_query"],
                    "Align_fraction_query": r["Align_fraction_ref"],
                    "Ref_name": r["Query_name"],
                    "Query_name": r["Ref_name"]})

out.sort(key=lambda r: -float(r["ANI"]))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
    w.writeheader(); w.writerows(out)

near = [r for r in out if float(r["ANI"]) >= 94.0]
print(f"wrote {OUT.name}: {len(out)} rows")
print(f"  ANI >= 94%: {len(near)} genomes")
for r in near:
    print(f"    {Path(r['Ref_file']).stem:20s} ANI {r['ANI']:>6s}  "
          f"AF_ref {r['Align_fraction_ref']:>6s}  AF_query {r['Align_fraction_query']:>6s}")
