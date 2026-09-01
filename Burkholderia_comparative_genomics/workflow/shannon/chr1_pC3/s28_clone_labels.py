#!/usr/bin/env python3
"""Assign every genome in the full (non-dereplicated) set to a 99% ANI clone cluster.

The 99% ANI dereplication is NOT applied to the rebuild -- but its cluster
assignments are retained and demoted from a *filter* to a *covariate*. One
cluster holds 233 of 771 genomes (the B. pseudomallei / B. mallei complex), so
any statistic that is sensitive to how many times a strain happened to be
sequenced needs a clone-aware companion estimate. This script produces the
labels every such estimate consumes.

Inputs
  --membership  derep_cluster_membership.tsv  (accession, cluster_id, cluster_size)
  --genomes     genome_list_full.txt          one accession per line, the full set
  --ani         optional skani edge list (query, ref, ANI, ...) used to place
                genomes that post-date the original dereplication (newly
                downloaded accessions, MF6, MF7) into an existing cluster at
                >= --threshold ANI, or as new singletons.

Output
  clone_cluster.tsv : accession, cluster_id, cluster_size, is_representative, source

`source` records how the label was obtained -- "derep" (inherited),
"ani_join" (placed by ANI), or "singleton_new" -- so the report can state
exactly how many genomes were labelled by inference rather than inheritance.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def read_membership(path: Path) -> dict[str, str]:
    """accession -> cluster_id from the retained dereplication output."""
    out: dict[str, str] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["accession"]] = row["cluster_id"]
    return out


def read_genomes(path: Path) -> list[str]:
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def read_ani(path: Path, threshold: float) -> dict[str, list[tuple[str, float]]]:
    """query -> [(ref, ani), ...] above threshold, best first.

    skani's `dist`/`triangle -E` output is a TSV whose first two columns are file
    paths and whose third is ANI. Accessions are recovered from the basenames, so
    this tolerates both `genomes/GCF_x.fna` and a bare `GCF_x`.
    """
    hits: dict[str, list[tuple[str, float]]] = collections.defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("Ref_file") or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ref, qry = parts[0], parts[1]
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            if ani < threshold:
                continue
            q = Path(qry).name.replace(".fna", "").replace(".fasta", "")
            r = Path(ref).name.replace(".fna", "").replace(".fasta", "")
            if q == r:
                continue
            hits[q].append((r, ani))
    for q in hits:
        hits[q].sort(key=lambda t: -t[1])
    return hits


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--ani", type=Path, default=None)
    ap.add_argument("--threshold", type=float, default=99.0)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    inherited = read_membership(args.membership)
    genomes = read_genomes(args.genomes)
    ani = read_ani(args.ani, args.threshold) if args.ani else {}

    assign: dict[str, tuple[str, str]] = {}   # acc -> (cluster_id, source)
    unplaced: list[str] = []
    for acc in genomes:
        if acc in inherited:
            assign[acc] = (inherited[acc], "derep")
        else:
            unplaced.append(acc)

    # Place post-dereplication genomes by ANI to an already-labelled genome.
    next_id = max((int(c) for c in inherited.values() if c.isdigit()), default=-1) + 1
    for acc in unplaced:
        placed = False
        for ref, _a in ani.get(acc, []):
            if ref in assign:
                assign[acc] = (assign[ref][0], "ani_join")
                placed = True
                break
        if not placed:
            assign[acc] = (str(next_id), "singleton_new")
            next_id += 1

    sizes = collections.Counter(cid for cid, _ in assign.values())

    # The representative is the genome the original dereplication kept, where one
    # exists; otherwise the lexicographically first member, so the choice is
    # deterministic and reproducible rather than dependent on dict ordering.
    by_cluster: dict[str, list[str]] = collections.defaultdict(list)
    for acc, (cid, _s) in assign.items():
        by_cluster[cid].append(acc)
    rep: dict[str, str] = {cid: sorted(m)[0] for cid, m in by_cluster.items()}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "cluster_id", "cluster_size",
                    "is_representative", "source"])
        for acc in genomes:
            cid, src = assign[acc]
            w.writerow([acc, cid, sizes[cid], str(acc == rep[cid]), src])

    src_counts = collections.Counter(s for _c, s in assign.values())
    big = sizes.most_common(8)
    print(f"genomes labelled : {len(genomes)}")
    print(f"clusters         : {len(sizes)}")
    print(f"sources          : {dict(src_counts)}")
    print(f"largest clusters : {big}")
    print(f"singletons       : {sum(1 for v in sizes.values() if v == 1)}")
    print(f"in clusters >=5  : {sum(v for v in sizes.values() if v >= 5)}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
