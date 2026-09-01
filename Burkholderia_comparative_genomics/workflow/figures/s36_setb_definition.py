#!/usr/bin/env python
"""Define Set B -- MF6's ANI >=95% clade -- and name each member's pC3 contig.

Set B is the basis of the setB_* functional tables and the fig7/fig8 screens.
D15: MF7 is INCLUDED (ANI 100.00). It was isolated separately, so it is an
independent sampling of what a B. sola pC3 carries, not a re-sequencing of MF6.
Its assembly is a 63-contig draft, so any family present in the other five and
absent only in MF7 must be checked against MF7's contig boundaries before being
called a loss -- see D15.

Writes tables/setB_members.tsv: accession, ani_to_mf6, pc3_contig, pc3_bp.
"""
import csv
from pathlib import Path

D = Path(__file__).resolve().parent.parent
TAB = D / "tables"
# The report says "ANI >=95%", but that text is WRONG about what was run. At 95%
# the rule admits GCF_003966315.1 (95.32%), which is B. cenocepacia -- a different
# species -- and that genome appears ZERO times in the retained screen outputs.
# The report's own prose gives the real rule: "the four closest genomes sit at
# >=96% -- the B. sola clade, MF6's own species -- and the fifth drops to 95.3%,
# into assorted B. cenocepacia". So the operational cut is 96%, and the "95%"
# in the report must be corrected.
ANI_MIN = 96.0

ani = {"MF6": 100.0}
for r in csv.DictReader(open(TAB / "MF6_ani_raw.tsv", newline=""), delimiter="\t"):
    ani[Path(r["Ref_file"]).name.replace(".fna", "")] = float(r["ANI"])

members = sorted((a for a, v in ani.items() if v >= ANI_MIN),
                 key=lambda a: (-ani[a], a))

# pC3 contig: the replicon the frozen classifier called c3.
pc3 = {}
for r in csv.DictReader(open(TAB / "secondary_replicon_clusters.tsv", newline=""),
                        delimiter="\t"):
    if r["accession"] in ani and r["is_c3"] == "True":
        prev = pc3.get(r["accession"])
        if prev is None or int(r["length"]) > prev[1]:
            pc3[r["accession"]] = (r["contig"], int(r["length"]))

rows = []
for a in members:
    if a not in pc3:
        raise SystemExit(f"FAIL: no c3-classified contig for Set B member {a}")
    c, bp = pc3[a]
    rows.append(dict(accession=a, ani_to_mf6=f"{ani[a]:.2f}", pc3_contig=c, pc3_bp=bp))

with open(TAB / "setB_members.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)

print(f"Set B (ANI >= {ANI_MIN}%): {len(rows)} genomes, "
      f"{sum(r['pc3_bp'] for r in rows):,} bp of pC3")
for r in rows:
    print(f"  {r['accession']:<20} ANI {r['ani_to_mf6']:>6}  "
          f"{r['pc3_contig'][:42]:<44} {r['pc3_bp']:>10,} bp")
nxt = max((v for a, v in ani.items() if v < ANI_MIN), default=0)
print(f"  (next genome below the cut: ANI {nxt:.2f}%)")
