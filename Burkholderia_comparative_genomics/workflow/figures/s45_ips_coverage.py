#!/usr/bin/env python3
"""Family-level InterProScan rescue table -> tables/screen_ips_coverage.tsv (fig8 panel A).

The screen is per PROTEIN; section 6.4's annotation unit is the gene FAMILY. This joins
the two: for every pC3 clade family, does InterProScan assign an InterPro entry
where COG assigned nothing?

A family counts as InterPro-annotated when a STRICT MAJORITY of its members carry
an InterPro entry -- the same rule s43 applies to COG and KO, and for the same
reason: one member out of six matching a signature is not the family's
annotation.

Emits one row per (gene_set, cog_status): n_families, n_interpro_entry.
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--families", required=True, type=Path,
                help="setB_functional_families_*.tsv (gene_set, family, cog_categories)")
ap.add_argument("--fam-members", required=True, type=Path,
                help="fam_pC3_*.tsv from s42, which carries the member list")
ap.add_argument("--signatures", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
a = ap.parse_args()


def tsv(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


ipr = {}
for r in tsv(a.signatures):
    _, _, lt = r["protein"].partition("|")
    ipr[lt or r["protein"]] = r["interpro_entry"] == "True"

members = {r["family"]: [m.split("|", 1)[1] for m in r["members"].split(";")]
           for r in tsv(a.fam_members)}

rows = defaultdict(lambda: [0, 0])
for f in tsv(a.families):
    gs = f["gene_set"].split("/")[0]
    if not gs.startswith("pC3_clade_"):
        continue
    key = (gs.replace("pC3_clade_", ""), "COG" if f["cog_categories"] else "noCOG")
    mem = members.get(f["family"], [])
    hit = sum(1 for t in mem if ipr.get(t))
    rows[key][0] += 1
    rows[key][1] += int(hit * 2 > len(mem))

with open(a.out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene_set", "cog_status", "n_families", "n_interpro_entry"])
    for k in sorted(rows):
        w.writerow([k[0], k[1], rows[k][0], rows[k][1]])

for gs in ("core", "accessory"):
    n = sum(v[0] for k, v in rows.items() if k[0] == gs)
    ncog = sum(v[0] for k, v in rows.items() if k == (gs, "COG"))
    nipr = sum(v[1] for k, v in rows.items() if k == (gs, "noCOG"))
    print(f"pC3 clade {gs:<9} n={n:>4}  COG {ncog:>4} ({100*ncog/n:.1f}%)  "
          f"+InterPro {nipr:>4}  -> {100*(ncog+nipr)/n:.1f}% annotated")
print(f"wrote {a.out}")
