#!/usr/bin/env python3
"""Collect Bakta proteins from every large secondary replicon, for MMseqs2 clustering.

Reproduces Stage 4 of the original pipeline, which the methods draft records as:

    "Proteins from all large secondary replicons were clustered into orthologous
     groups with MMseqs2 v18.8cc5c (easy-cluster, 50% identity, 80% coverage)."

Note this is NOT the Panaroo c3 pangenome (Stage 7) -- that is a separate
deliverable over the pC3-carrying subset. Stage 4 clusters ALL secondary
replicons at once, which is why it can classify every replicon in the set
without any training/application split.

Sequence ids are written as `accession|contig|locus_tag`, matching the retained
orthogroup identifiers (e.g. `GCF_964265025.1|NZ_OZ183803.1|BIENCE_06412`), so
residence can be computed by parsing the id alone.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def cds_by_contig(gff: Path, wanted: set[str]) -> dict[str, str]:
    """locus_tag -> contig, for CDS features on the wanted contigs."""
    out = {}
    for ln in open(gff):
        if ln.startswith("#"):
            if ln.startswith("##FASTA"):
                break
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "CDS" or f[0] not in wanted:
            continue
        attrs = f[8]
        tag = None
        for kv in attrs.strip().split(";"):
            if kv.startswith("locus_tag="):
                tag = kv[len("locus_tag="):].strip()
                break
        if tag:
            out[tag] = f[0]
    return out


def read_faa(path: Path):
    name, buf = None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name:
                yield name, "".join(buf)
            name = ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if name:
        yield name, "".join(buf)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types", required=True, type=Path)
    ap.add_argument("--annot", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--kind", default="secondary_large")
    a = ap.parse_args(argv)

    want = collections.defaultdict(set)
    for r in csv.DictReader(open(a.types), delimiter="\t"):
        if r["replicon_type"] == a.kind:
            want[r["accession"]].add(r["contig"])
    n_rep = sum(len(v) for v in want.values())
    print(f"genomes {len(want)}   {a.kind} replicons {n_rep}", flush=True)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    missing_gff = missing_faa = 0
    empty = []
    with open(a.out, "w") as oh:
        for i, (acc, contigs) in enumerate(sorted(want.items()), 1):
            gff = a.annot / acc / f"{acc}.gff3"
            faa = a.annot / acc / f"{acc}.faa"
            if not gff.exists():
                missing_gff += 1; continue
            if not faa.exists():
                missing_faa += 1; continue
            loc2ctg = cds_by_contig(gff, contigs)
            if not loc2ctg:
                empty.append(acc); continue
            n_here = 0
            for tag, seq in read_faa(faa):
                ctg = loc2ctg.get(tag)
                if ctg is None:
                    continue
                oh.write(f">{acc}|{ctg}|{tag}\n")
                for j in range(0, len(seq), 60):
                    oh.write(seq[j:j + 60] + "\n")
                n_here += 1
            written += n_here
            if i % 100 == 0:
                print(f"  {i}/{len(want)} genomes, {written:,d} proteins", flush=True)

    print(f"wrote {written:,d} proteins to {a.out}")
    if missing_gff or missing_faa:
        print(f"WARNING missing gff3 {missing_gff}, missing faa {missing_faa}",
              file=sys.stderr)
    if empty:
        print(f"WARNING {len(empty)} genomes had no CDS on their {a.kind} "
              f"contigs: {empty[:6]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
