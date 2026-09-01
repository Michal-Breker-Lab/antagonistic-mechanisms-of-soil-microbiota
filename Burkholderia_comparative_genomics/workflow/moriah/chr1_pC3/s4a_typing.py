#!/usr/bin/env python3
"""s4a_typing.py -- per-contig replicon typing from Bakta GFF3.

REWRITE of a script that was stranded on Shannon (2026-08-20 outage).  The
rules below were reverse-engineered from the retained output
`tables/replicon_types.tsv` (1,107 rows / 306 genomes) and every threshold
sits inside a clean, unambiguous gap in that table -- see `--validate`.

Rules
-----
chromosome1        the contig carrying the ribosomal-protein cluster, i.e.
                   argmax(n ribosomal-protein genes).  Verified: 0 violations
                   in 306 retained genomes.
chromosome1_fused  that contig when it is >= 4.5 Mb.  Burkholderia c1 runs
                   ~3.3-4.0 Mb, so a longer one has c1 and c2 fused into a
                   single replicon.  Retained data: normal c1 tops out at
                   4,485,038 bp, fused starts at 4,531,000 bp -- a 46 kb gap
                   straddling 4.5 Mb.
secondary_large    any other contig >= 300 kb.  Retained gap: largest
                   small_plasmid 299,296 bp, smallest secondary_large
                   301,592 bp.
small_plasmid      any other contig < 300 kb.

Outputs `replicon_types.tsv` with the retained column order:
    accession contig size_rank length cds ribosomal_proteins
    core_markers n_core_markers replicon_type
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from urllib.parse import unquote

FUSED_MIN = 4_500_000
LARGE_MIN = 300_000

CORE_MARKERS = ["dnaA", "dnaN", "ftsZ", "gyrB", "recA", "rpoB", "rpoC"]

# rpsA / rplB / rpmC ... -- the 30S, 50S and 50S-L7/L12 gene families.  Bakta
# names these from UniRef, so the gene= attribute is the reliable handle; the
# product string is not (it also matches "ribosomal protein methyltransferase").
RIBO_RE = re.compile(r"^rp[sml][A-Z]\d*$")

GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")


def parse_gff3(path: Path):
    """-> (lengths {seqid: int} in file order, cds {seqid: int},
             ribo {seqid: int}, markers {seqid: set})"""
    lengths: dict[str, int] = {}
    cds: dict[str, int] = {}
    ribo: dict[str, int] = {}
    markers: dict[str, set] = {}
    core = set(CORE_MARKERS)
    with open(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    sid = p[1]
                    lengths[sid] = int(p[3])
                    cds.setdefault(sid, 0)
                    ribo.setdefault(sid, 0)
                    markers.setdefault(sid, set())
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            sid = f[0]
            cds[sid] = cds.get(sid, 0) + 1
            m = GENE_RE.search(f[8])
            if not m:
                continue
            gene = unquote(m.group(1)).strip()
            if RIBO_RE.match(gene):
                ribo[sid] = ribo.get(sid, 0) + 1
            elif gene in core:
                markers.setdefault(sid, set()).add(gene)
    return lengths, cds, ribo, markers


def type_genome(acc: str, lengths, cds, ribo, markers):
    """-> list of output rows, ordered by descending contig length."""
    if not lengths:
        return []
    order = sorted(lengths, key=lambda s: (-lengths[s], s))
    # chromosome 1 = the ribosomal-protein-cluster contig.  Ties (all-zero
    # genomes with no annotated r-proteins) fall to the longest contig, which
    # is what argmax over the length-sorted order gives.
    c1 = max(order, key=lambda s: (ribo.get(s, 0), lengths[s]))
    rows = []
    for rank, sid in enumerate(order, start=1):
        if sid == c1:
            rtype = ("chromosome1_fused" if lengths[sid] >= FUSED_MIN
                     else "chromosome1")
        elif lengths[sid] >= LARGE_MIN:
            rtype = "secondary_large"
        else:
            rtype = "small_plasmid"
        mk = sorted(markers.get(sid, ()), key=CORE_MARKERS.index)
        rows.append({
            "accession": acc,
            "contig": sid,
            "size_rank": rank,
            "length": lengths[sid],
            "cds": cds.get(sid, 0),
            "ribosomal_proteins": ribo.get(sid, 0),
            "core_markers": ";".join(mk),
            "n_core_markers": len(mk),
            "replicon_type": rtype,
        })
    return rows


def validate(new_rows, ref_path: Path):
    """Compare against the retained table, genome by genome."""
    ref_rows = list(csv.DictReader(open(ref_path), delimiter="\t"))
    ref = {(r["accession"], r["contig"]): r for r in ref_rows}
    new = {(r["accession"], str(r["contig"])): r for r in new_rows}
    ref_acc = {a for a, _ in ref}
    new_acc = {a for a, _ in new}
    shared_acc = ref_acc & new_acc
    print(f"  reference genomes : {len(ref_acc)}")
    print(f"  new genomes       : {len(new_acc)}")
    print(f"  shared genomes    : {len(shared_acc)}")
    missing = [k for k in ref if k[0] in shared_acc and k not in new]
    if missing:
        print(f"  contigs in reference but not new: {len(missing)} {missing[:4]}")
    ok = not missing
    for field in ["size_rank", "length", "cds", "ribosomal_proteins",
                  "core_markers", "n_core_markers", "replicon_type"]:
        diffs = []
        for k, rr in ref.items():
            if k[0] not in shared_acc or k not in new:
                continue
            rv, nv = rr[field], str(new[k][field])
            if rv != nv:
                diffs.append((k, rv, nv))
        if diffs:
            ok = False
            print(f"  {field:<20} {len(diffs)} MISMATCH")
            for k, rv, nv in diffs[:6]:
                print(f"      {k[0]} {k[1]}: ref={rv!r} new={nv!r}")
        else:
            print(f"  {field:<20} MATCH")
    print(f"  VERDICT: {'reproduces the retained table' if ok else 'DIVERGENCE'}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True, type=Path,
                    help="dir of <acc>/<acc>.gff3 Bakta outputs")
    ap.add_argument("--list", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--validate", type=Path, default=None)
    a = ap.parse_args(argv)

    accs = [l.strip() for l in open(a.list) if l.strip()]
    rows = []
    missing = []
    for acc in accs:
        gff = a.annot / acc / f"{acc}.gff3"
        if not gff.exists():
            missing.append(acc)
            continue
        rows.extend(type_genome(acc, *parse_gff3(gff)))
    if missing:
        print(f"WARNING: {len(missing)} genomes without GFF3: {missing[:5]}",
              file=sys.stderr)

    a.outdir.mkdir(parents=True, exist_ok=True)
    cols = ["accession", "contig", "size_rank", "length", "cds",
            "ribosomal_proteins", "core_markers", "n_core_markers",
            "replicon_type"]
    with open(a.outdir / "replicon_types.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_gen = len({r["accession"] for r in rows})
    print(f"wrote {len(rows)} contig rows for {n_gen} genomes -> {a.outdir}")

    if a.validate:
        print("\n=== validation against retained table ===")
        return 0 if validate(rows, a.validate) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
