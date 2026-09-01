#!/usr/bin/env python
"""Carry host_categories.tsv forward from the original 772-genome curation to
the rebuild's 773-genome set.

The delta is four rows and nothing else, so the curated evidence strings are
carried over verbatim rather than re-derived:

  GCA_050858985.1 -> GCF_050858985.1   RefSeq mirror of the same GenBank
  GCA_050955445.1 -> GCF_050955445.1   assembly; identical BioSample.
  GCF_057005415.1 -> GCF_057005415.2   version bump of the same RefSeq
                                       assembly; identical BioSample.
  MF7                (new)             same lab isolate as MF6 -- ANI 100.00,
                                       clone cluster 22 -- so it inherits MF6's
                                       curated row.

Genus/species were checked against rebuild/genus_species.tsv for all three
renames before carrying the rows over.
"""
import csv
import sys
from pathlib import Path

D = Path(__file__).resolve().parent.parent
SRC = D.parent / "tables" / "host_categories.tsv"
OUT = D / "tables" / "host_categories.tsv"
WANT = D / "genome_list_full.txt"

RENAME = {"GCA_050858985.1": "GCF_050858985.1",
          "GCA_050955445.1": "GCF_050955445.1",
          "GCF_057005415.1": "GCF_057005415.2"}

with open(SRC) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    fields = rd.fieldnames
    rows = {r["accession"]: r for r in rd}

for old, new in RENAME.items():
    if old not in rows:
        sys.exit(f"FAIL: expected {old} in {SRC.name}")
    r = rows.pop(old)
    r["accession"] = new
    r["evidence"] = (r["evidence"] + "; accession updated "
                     f"{old} -> {new} (same BioSample)").lstrip("; ")
    rows[new] = r

mf6 = rows["MF6"]
rows["MF7"] = dict(mf6, accession="MF7",
                   organism_name="Burkholderia sola (MF7)",
                   evidence="manual curation: same lab isolate as MF6 "
                            "(ANI 100.00, clone cluster 22); inherits MF6 row")

want = [l.strip() for l in open(WANT) if l.strip()]
missing = [a for a in want if a not in rows]
extra = [a for a in rows if a not in set(want)]
if missing or extra:
    sys.exit(f"FAIL: {len(missing)} missing {missing[:5]}, "
             f"{len(extra)} extra {extra[:5]}")

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for a in want:
        w.writerow(rows[a])

cats = {}
for a in want:
    cats[rows[a]["host_category"]] = cats.get(rows[a]["host_category"], 0) + 1
print(f"wrote {OUT.relative_to(D.parent)}  ({len(want)} genomes)")
for c, n in sorted(cats.items(), key=lambda kv: -kv[1]):
    print(f"  {c:20s} {n:4d}")
