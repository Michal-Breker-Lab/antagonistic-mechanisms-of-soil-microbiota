#!/usr/bin/env python3
"""Resolve, per Set B genome, EVERY contig belonging to pC3 and to chromosome 1.

Why this exists: s36 wrote setB_members.tsv with a single `pc3_contig` column and
s37 carved that one contig. For the five closed genomes one contig IS the
replicon, so nothing was lost -- but MF7 is a draft assembly whose pC3 is split
across TWO contigs (NODE_3 857,810 bp + NODE_8 318,422 bp = 1,176,232 bp, i.e.
99.7% of MF6's pC3) and whose chromosome 1 is split across EIGHT. Carving one
contig silently dropped a quarter of MF7's pC3 and made a complete replicon look
like a truncation.

Closed genomes: chromosome 1 is the contig carrying the full core-marker set
(dnaA/dnaN/ftsZ/gyrB/recA/rpoB/rpoC) and the ribosomal superoperon; pC3 is the
contig named in setB_members.tsv. Both come from replicon_types.tsv.

MF7: single-contig marker typing cannot work on a draft assembly -- the markers
scatter across contigs -- so MF7's assignment comes from the skani-based
mf7_replicon_assignment.tsv, taking every contig with call == "assigned".

Emits setB_contigs.tsv: accession, replicon, contig, length, cds
"""
import argparse
import csv
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--members", required=True, type=Path)
ap.add_argument("--replicon-types", required=True, type=Path)
ap.add_argument("--mf7-assignment", required=True, type=Path)
ap.add_argument("--out", required=True, type=Path)
a = ap.parse_args()


def tsv(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


members = tsv(a.members)
rtypes = tsv(a.replicon_types)
by_acc = {}
for r in rtypes:
    by_acc.setdefault(r["accession"], []).append(r)

CORE = {"dnaA", "dnaN", "ftsZ", "gyrB", "recA", "rpoB", "rpoC"}
rows = []
for m in members:
    acc = m["accession"]
    if acc == "MF7":
        for r in tsv(a.mf7_assignment):
            if r["call"] != "assigned" or r["replicon"] not in ("pC3", "chr1"):
                continue
            rep = "pC3" if r["replicon"] == "pC3" else "chr1"
            rows.append(dict(accession=acc, replicon=rep, contig=r["contig"],
                             length=r["length"], cds=""))
        continue

    ctgs = by_acc[acc]
    chr1 = [c for c in ctgs
            if CORE.issubset(set(filter(None, c["core_markers"].split(";"))))]
    if len(chr1) != 1:
        raise SystemExit(f"FAIL: {acc} has {len(chr1)} full-marker contigs")
    rows.append(dict(accession=acc, replicon="chr1", contig=chr1[0]["contig"],
                     length=chr1[0]["length"], cds=chr1[0]["cds"]))

    want = m["pc3_contig"]
    pc3 = [c for c in ctgs if c["contig"] == want]
    if len(pc3) != 1:
        raise SystemExit(f"FAIL: {acc} pC3 contig {want} not found in replicon_types")
    rows.append(dict(accession=acc, replicon="pC3", contig=pc3[0]["contig"],
                     length=pc3[0]["length"], cds=pc3[0]["cds"]))

with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, ["accession", "replicon", "contig", "length", "cds"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"{'genome':<18}{'replicon':<8}{'contigs':>8}{'total bp':>14}")
seen = {}
for r in rows:
    seen.setdefault((r["accession"], r["replicon"]), []).append(int(r["length"]))
for (acc, rep), lens in sorted(seen.items()):
    print(f"{acc:<18}{rep:<8}{len(lens):>8}{sum(lens):>14,}")
print(f"\nwrote {a.out}  ({len(rows)} contigs)")
