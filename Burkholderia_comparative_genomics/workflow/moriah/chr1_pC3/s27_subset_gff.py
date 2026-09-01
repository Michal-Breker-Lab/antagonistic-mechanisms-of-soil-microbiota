#!/usr/bin/env python3
"""Subset each training genome's Bakta GFF3 to its two large secondary replicons.

The pC3 diagnostics are derived from a pangenome over the SECONDARY replicons
only (chromosome 2 + pC3), not whole genomes -- see the scope note in the D9
approval. This writes one Panaroo-ready GFF3 per training genome containing just
the two labelled contigs, with the matching FASTA records preserved in the
##FASTA block (Panaroo requires the embedded sequence).
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def split_gff(path: Path):
    """Return (header_lines, feature_lines, {seqid: fasta_record_lines})."""
    head, feats, fasta = [], [], collections.OrderedDict()
    cur = None
    in_fasta = False
    for ln in open(path):
        if ln.startswith("##FASTA"):
            in_fasta = True
            continue
        if in_fasta:
            if ln.startswith(">"):
                cur = ln[1:].split()[0]
                fasta[cur] = [ln]
            elif cur:
                fasta[cur].append(ln)
        elif ln.startswith("#"):
            head.append(ln)
        elif ln.strip():
            feats.append(ln)
    return head, feats, fasta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels", required=True, type=Path)
    ap.add_argument("--annot", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    keep = collections.defaultdict(set)
    for r in csv.DictReader(open(a.labels), delimiter="\t"):
        keep[r["accession"]].add(r["contig"])

    ok = bad = 0
    for acc, contigs in sorted(keep.items()):
        src = a.annot / acc / f"{acc}.gff3"
        if not src.exists():
            print(f"MISSING {src}", file=sys.stderr); bad += 1; continue
        head, feats, fasta = split_gff(src)
        f_keep = [l for l in feats if l.split("\t", 1)[0] in contigs]
        s_keep = {k: v for k, v in fasta.items() if k in contigs}
        if len(s_keep) != len(contigs):
            print(f"WARNING {acc}: wanted {sorted(contigs)}, "
                  f"found sequences {sorted(s_keep)}", file=sys.stderr)
            bad += 1
            continue
        if not f_keep:
            print(f"WARNING {acc}: no features on {sorted(contigs)}", file=sys.stderr)
            bad += 1
            continue
        out = a.outdir / f"{acc}.gff3"
        with open(out, "w") as fh:
            fh.write("##gff-version 3\n")
            for c in sorted(contigs):
                n = sum(len(l.strip()) for l in s_keep[c][1:])
                fh.write(f"##sequence-region {c} 1 {n}\n")
            fh.writelines(f_keep)
            fh.write("##FASTA\n")
            for c in sorted(contigs):
                fh.writelines(s_keep[c])
        ok += 1
    print(f"wrote {ok} subset GFF3s to {a.outdir}   (problems: {bad})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
