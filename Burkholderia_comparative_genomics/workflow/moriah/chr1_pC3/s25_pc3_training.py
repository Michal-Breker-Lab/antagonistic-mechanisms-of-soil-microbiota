#!/usr/bin/env python3
"""Build the clone-collapsed training set for pC3 diagnostic derivation (D9a).

Training genomes are those with EXACTLY two large secondary replicons, where the
positional assignment is unambiguous: the larger is chromosome 2, the smaller is
pC3. One genome is taken per 99% ANI clone cluster, because without
dereplication the training set is clone-dominated and a family carried by 200
near-identical Bcc strains would score as "diagnostic" on what is effectively a
single observation.

The representative for each cluster is chosen deterministically: among that
cluster's eligible (two-secondary) genomes, the one whose two secondary
replicons are most clearly separated in size, tie-broken by accession. A clean
size gap is exactly what makes the positional label trustworthy, so this picks
the cluster member whose training label is least likely to be wrong.

Outputs
  pc3_training_set.txt      one accession per line
  pc3_training_labels.tsv   accession, contig, role (chromosome2|pC3), length
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types", required=True, type=Path)
    ap.add_argument("--clones", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args(argv)

    sec = collections.defaultdict(list)          # acc -> [(len, contig)]
    for r in csv.DictReader(open(a.types), delimiter="\t"):
        if r["replicon_type"] == "secondary_large":
            sec[r["accession"]].append((int(r["length"]), r["contig"]))

    eligible = {acc: sorted(v, reverse=True) for acc, v in sec.items() if len(v) == 2}
    print(f"genomes with exactly two secondary_large replicons: {len(eligible)}")

    clone = {r["accession"]: r["cluster_id"]
             for r in csv.DictReader(open(a.clones), delimiter="\t")}

    # group eligible genomes by clone cluster, pick the clearest size separation
    by_cluster = collections.defaultdict(list)
    for acc in eligible:
        by_cluster[clone.get(acc, f"_unlabelled_{acc}")].append(acc)

    chosen = []
    for cid, accs in sorted(by_cluster.items()):
        def sep(acc):
            big, small = eligible[acc][0][0], eligible[acc][1][0]
            return (-(big - small), acc)      # widest gap first, then accession
        chosen.append(sorted(accs, key=sep)[0])
    chosen.sort()

    a.outdir.mkdir(parents=True, exist_ok=True)
    (a.outdir / "pc3_training_set.txt").write_text("\n".join(chosen) + "\n")

    with open(a.outdir / "pc3_training_labels.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "contig", "role", "length", "clone_cluster"])
        for acc in chosen:
            (lb, cb), (ls, cs) = eligible[acc]
            w.writerow([acc, cb, "chromosome2", lb, clone.get(acc, "")])
            w.writerow([acc, cs, "pC3", ls, clone.get(acc, "")])

    gaps = [eligible[acc][0][0] - eligible[acc][1][0] for acc in chosen]
    gaps.sort()
    print(f"clone clusters represented : {len(chosen)}")
    print(f"collapsed from             : {len(eligible)} eligible genomes")
    print(f"size gap chr2-pC3 (bp)     : min {gaps[0]:,d}  median "
          f"{gaps[len(gaps)//2]:,d}  max {gaps[-1]:,d}")
    tight = [acc for acc in chosen
             if eligible[acc][0][0] - eligible[acc][1][0] < 100_000]
    print(f"training genomes whose two replicons differ by <100 kb: {len(tight)}")
    if tight:
        print("  (positional label is least certain for these)")
        for acc in tight[:10]:
            (lb, _), (ls, _) = eligible[acc]
            print(f"    {acc}  chr2={lb:,d}  pC3={ls:,d}  gap={lb-ls:,d}")
    print(f"wrote {a.outdir/'pc3_training_set.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
