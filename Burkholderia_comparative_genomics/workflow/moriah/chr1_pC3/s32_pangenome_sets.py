#!/usr/bin/env python3
"""Emit the contig lists (s27 label format) for the pC3 and chromosome-1 pangenomes.

pC3 set  : every replicon classified `c3` by Stage 4 (D11).
chr1 set : the chromosome-1 replicon of every genome, including fused ones.

Written in the same two-column-plus-role layout s27_subset_gff.py consumes, so
one subsetting path serves every pangenome arm.
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--types", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    c3 = [(r["accession"], r["contig"], r["length"])
          for r in csv.DictReader(open(a.clusters), delimiter="\t")
          if r["replicon_class"] == "c3"]
    ch1 = [(r["accession"], r["contig"], r["length"])
           for r in csv.DictReader(open(a.types), delimiter="\t")
           if r["replicon_type"] in ("chromosome1", "chromosome1_fused")]

    for name, rows in (("pc3", c3), ("chr1", ch1)):
        p = a.outdir / f"pangenome_{name}_labels.tsv"
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["accession", "contig", "role", "length"])
            for acc, ctg, ln in rows:
                w.writerow([acc, ctg, name, ln])
        print(f"{name}: {len({r[0] for r in rows})} genomes, {len(rows)} replicons -> {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
