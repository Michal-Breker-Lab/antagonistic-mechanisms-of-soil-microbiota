#!/usr/bin/env python3
"""Assign MF7's draft contigs to replicons by alignment to the closed MF6 genome.

MF7 (IMG taxon 2546825540) is a 63-contig, 8.70 Mb SPAdes draft with N50 858 kb;
8 contigs carry 7.0 Mb. Its replicon membership therefore cannot be read off the
within-genome size rank the way it is for the 771 complete genomes. But MF6 from
the same collection is a closed 5-replicon 7.80 Mb assembly, so MF7's contigs can
be placed by alignment instead.

This is a DIFFERENT method from the one used for the rest of the set and must be
footnoted as such wherever an MF7 replicon call is used.

Input : a PAF from `minimap2 -cx asm10 MF6.fna MF7.fna` (MF6 = reference).
Output: mf7_replicon_assignment.tsv -- one row per MF7 contig with the winning
        MF6 replicon, the aligned fraction, the margin over the runner-up, and a
        call of assigned / ambiguous / unplaced.

A contig is assigned only when the winner covers >= --min-frac of the contig AND
beats the runner-up by >= --min-margin fold. Everything else is reported as
ambiguous or unplaced rather than forced -- a wrong replicon call on MF7 would
propagate straight into the pC3 presence table.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def _merge(ivs: list[tuple[int, int]]) -> int:
    """Total length covered by a set of intervals, counting overlaps once."""
    if not ivs:
        return 0
    ivs = sorted(ivs)
    tot = 0
    cs, ce = ivs[0]
    for s, e in ivs[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        elif e > ce:
            ce = e
    return tot + (ce - cs)


def parse_alignments(path: Path):
    """query -> {target: covered_bases}, and query -> query_length.

    Accepts minimap2 PAF or blastn -outfmt "6 qseqid qlen qstart qend sstrand
    sseqid slen sstart send nident length", whose columns coincide with PAF's
    first eleven.

    Coverage is computed by MERGING query intervals, not by summing hit lengths.
    BLAST reports many HSPs per contig pair and they overlap, so summing would
    push aligned_frac above 1 and could hand a contig to the wrong replicon.
    Coordinates are normalised because BLAST reverses qstart/qend on minus-strand
    hits.
    """
    ivs: dict[str, dict[str, list[tuple[int, int]]]] = collections.defaultdict(
        lambda: collections.defaultdict(list))
    qlen: dict[str, int] = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            except ValueError:
                continue
            t = f[5]
            if qs > qe:
                qs, qe = qe, qs
            qlen[q] = ql
            ivs[q][t].append((qs, qe))

    aln: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for q, per_t in ivs.items():
        for t, lst in per_t.items():
            aln[q][t] = _merge(lst)
    return aln, qlen


# backwards-compatible alias
parse_paf = parse_alignments


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True, type=Path,
                    help="minimap2 PAF, or blastn -outfmt 6 with PAF-compatible columns")
    ap.add_argument("--replicon-map", type=Path, default=None,
                    help="TSV: MF6 contig name -> replicon label (chr1/chr2/pC3/...). "
                         "Without it, MF6 contig names are used as the labels.")
    ap.add_argument("--min-frac", type=float, default=0.50)
    ap.add_argument("--min-margin", type=float, default=2.0)
    ap.add_argument("--min-contig", type=int, default=1000,
                    help="contigs shorter than this are reported but never assigned; "
                         "the MF7 draft carries SPAdes fragments down to 147 bp")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    label = {}
    if args.replicon_map and args.replicon_map.exists():
        with open(args.replicon_map) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2 and not row[0].startswith("#"):
                    label[row[0]] = row[1]

    aln, qlen = parse_alignments(args.paf)
    rows = []
    tally = collections.Counter()
    bp = collections.Counter()

    for q in sorted(qlen, key=lambda x: -qlen[x]):
        ql = qlen[q]
        hits = aln[q].most_common()
        best_t, best_b = (hits[0] if hits else ("", 0))
        second_b = hits[1][1] if len(hits) > 1 else 0
        frac = best_b / ql if ql else 0.0
        margin = (best_b / second_b) if second_b else float("inf")
        rep = label.get(best_t, best_t)

        if ql < args.min_contig:
            call = "too_short"
        elif not hits or frac < args.min_frac:
            call = "unplaced"
        elif margin < args.min_margin:
            call = "ambiguous"
        else:
            call = "assigned"

        if call == "assigned":
            tally[rep] += 1
            bp[rep] += ql
        rows.append(dict(contig=q, length=ql, replicon=(rep if call == "assigned" else ""),
                         best_target=best_t, aligned_bp=best_b,
                         aligned_frac=f"{frac:.4f}",
                         margin=("inf" if margin == float("inf") else f"{margin:.2f}"),
                         call=call))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    calls = collections.Counter(r["call"] for r in rows)
    total_bp = sum(qlen.values())
    placed_bp = sum(bp.values())
    print(f"MF7 contigs           : {len(rows)}  ({total_bp:,} bp)")
    print(f"calls                 : {dict(calls)}")
    print(f"assigned bp           : {placed_bp:,} ({placed_bp/total_bp:.1%})")
    for rep, n in tally.most_common():
        print(f"  {rep:<24} {n:>3} contigs  {bp[rep]:>10,} bp")
    unl = [r for r in rows if r["call"] in ("unplaced", "ambiguous")
           and int(r["length"]) >= 50000]
    if unl:
        print(f"\nlarge contigs NOT confidently placed ({len(unl)}):")
        for r in unl:
            print(f"  {r['contig'][:48]:<50} {int(r['length']):>9,} bp  "
                  f"frac={r['aligned_frac']} margin={r['margin']} [{r['call']}]")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
