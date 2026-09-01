#!/usr/bin/env python3
"""Cluster the full (non-dereplicated) genome set into 99% ANI clone groups.

The rebuild does NOT dereplicate. But the original 99% ANI dereplication is not
discarded either -- it is demoted from a *filter* to a *covariate*. One cluster
holds 233 of 771 genomes (the B. pseudomallei complex), so any statistic that is
sensitive to how many times a strain happened to be sequenced needs a clone-aware
companion estimate. This script produces the labels those estimates consume.

Method is reproduced verbatim from the original run, as recorded in
logs/s5_derep.log:

    skani triangle -E   ->  keep edges with ANI >= 99.0 AND AF >= 50.0
                        ->  single-linkage connected components

`AF` is min(align_fraction_ref, align_fraction_query) -- i.e. BOTH genomes
must be >=50% covered. The original log line ("edges >= 99.0% ANI and >= 50.0%
AF") is ambiguous between min and max; min is what reproduces its 256 clusters,
and is the defensible reading: max would let a 3.5 Mb genome be called a clone
of an 8 Mb genome that merely contains it. Running this over the original
771 accessions must return the original 256 clusters; --validate asserts that
partition equality before the full 773-genome labels are trusted.

Note on --min-af: the original skani call passed `--min-af 3`, this rebuild uses
skani's default 15. Immaterial here -- both are far below the 50% clustering
cut, so no edge that could join a cluster is lost -- but it means the raw edge
lists are not byte-comparable and validation is done on the partition instead.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

ANI_MIN_DEFAULT = 99.0
AF_MIN_DEFAULT = 50.0


def acc_of(path_or_name: str) -> str:
    """Recover an accession from a skani file field or a bare name."""
    n = Path(path_or_name.strip()).name
    for suf in (".fna.gz", ".fasta.gz", ".fna", ".fasta", ".fa"):
        if n.endswith(suf):
            return n[: -len(suf)]
    return n


def read_edges(path: Path, ani_min: float, af_min: float, af_mode: str = "min"):
    """Yield (acc_a, acc_b) for every skani pair passing both thresholds.

    Columns are located by header name so the parser survives skani's
    --short-header and any future column reordering.
    """
    kept = seen = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_ref = header.index("Ref_file")
            i_qry = header.index("Query_file")
            i_ani = header.index("ANI")
            i_afr = header.index("Align_fraction_ref")
            i_afq = header.index("Align_fraction_query")
        except ValueError:
            sys.exit(f"unexpected skani header in {path}: {header}")
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) <= i_afq:
                continue
            seen += 1
            try:
                ani = float(p[i_ani])
                afr, afq = float(p[i_afr]), float(p[i_afq])
                af = min(afr, afq) if af_mode == "min" else max(afr, afq)
            except ValueError:
                continue
            if ani < ani_min or af < af_min:
                continue
            a, b = acc_of(p[i_ref]), acc_of(p[i_qry])
            if a == b:
                continue
            kept += 1
            yield a, b
    print(f"  edges read {seen}, kept (ANI>={ani_min}, {af_mode}(AF)>={af_min}) {kept}",
          file=sys.stderr)


class Union:
    """Single-linkage connected components."""

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, x: str) -> None:
        self.parent.setdefault(x, x)

    def find(self, x: str) -> str:
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:      # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def partition(assign: dict[str, str]) -> frozenset:
    """Cluster-label-independent representation, for comparing two clusterings."""
    groups: dict[str, set] = collections.defaultdict(set)
    for acc, cid in assign.items():
        groups[cid].add(acc)
    return frozenset(frozenset(g) for g in groups.values())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ani", required=True, type=Path,
                    help="skani triangle -E edge list over the full set")
    ap.add_argument("--genomes", required=True, type=Path,
                    help="one accession per line; the full non-dereplicated set")
    ap.add_argument("--ani-min", type=float, default=ANI_MIN_DEFAULT)
    ap.add_argument("--af-min", type=float, default=AF_MIN_DEFAULT)
    ap.add_argument("--af-mode", choices=["min", "max"], default="min",
                    help="min = require RECIPROCAL coverage of both genomes "
                         "(reproduces the original run); max = either genome. "
                         "min prevents a large genome being called a clone of a "
                         "small one merely contained within it.")
    ap.add_argument("--validate", type=Path, default=None,
                    help="original derep_cluster_membership.tsv to reproduce")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    genomes = [ln.strip() for ln in open(args.genomes)
               if ln.strip() and not ln.startswith("#")]
    gset = set(genomes)
    if len(gset) != len(genomes):
        dup = [g for g, n in collections.Counter(genomes).items() if n > 1]
        sys.exit(f"duplicate accessions in --genomes: {dup[:5]}")
    print(f"genomes: {len(genomes)}")

    uf = Union()
    for g in genomes:
        uf.add(g)
    n_edge = 0
    for a, b in read_edges(args.ani, args.ani_min, args.af_min, args.af_mode):
        if a in gset and b in gset:        # ignore anything not in the study set
            uf.union(a, b)
            n_edge += 1
    print(f"edges joining study genomes: {n_edge}")

    assign = {g: uf.find(g) for g in genomes}
    sizes = collections.Counter(assign.values())
    print(f"clusters: {len(sizes)}  (from {len(genomes)} genomes)")

    # ---- validation: reproduce the original partition on the shared accessions
    if args.validate:
        orig = {}
        with open(args.validate) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                orig[row["accession"]] = row["cluster_id"]
        shared = sorted(set(orig) & gset)
        print(f"\n=== validation vs {args.validate.name} ===")
        print(f"original accessions {len(orig)}, shared with rebuild {len(shared)}")
        missing = sorted(set(orig) - gset)
        added = sorted(gset - set(orig))
        if missing:
            print(f"in original but NOT in rebuild ({len(missing)}): {missing[:8]}")
        if added:
            print(f"new in rebuild ({len(added)}): {added[:8]}")

        # Re-cluster using ONLY the shared accessions, so the comparison is
        # like-for-like: new genomes must not be allowed to bridge two original
        # clusters and then be scored as a disagreement.
        uf2 = Union()
        for g in shared:
            uf2.add(g)
        sh = set(shared)
        for a, b in read_edges(args.ani, args.ani_min, args.af_min, args.af_mode):
            if a in sh and b in sh:
                uf2.union(a, b)
        mine = {g: uf2.find(g) for g in shared}
        theirs = {g: orig[g] for g in shared}
        pm, pt = partition(mine), partition(theirs)
        print(f"clusters  rebuild {len(pm)}   original {len(pt)}")
        if pm == pt:
            print("PARTITION MATCH -- clone clustering reproduces the original exactly")
        else:
            only_mine = pm - pt
            only_theirs = pt - pm
            print(f"PARTITION MISMATCH: {len(only_mine)} rebuild-only groups, "
                  f"{len(only_theirs)} original-only groups")
            for g in sorted(only_mine, key=len, reverse=True)[:5]:
                print(f"  rebuild-only  n={len(g)}: {sorted(g)[:6]}")
            for g in sorted(only_theirs, key=len, reverse=True)[:5]:
                print(f"  original-only n={len(g)}: {sorted(g)[:6]}")

    # ---- stable, human-meaningful cluster ids: 0 = largest, ties by first member
    order = sorted(sizes, key=lambda r: (-sizes[r], min(
        a for a in genomes if assign[a] == r)))
    renum = {root: str(i) for i, root in enumerate(order)}
    members: dict[str, list[str]] = collections.defaultdict(list)
    for g in genomes:
        members[renum[assign[g]]].append(g)
    rep = {cid: sorted(m)[0] for cid, m in members.items()}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "cluster_id", "cluster_size", "is_representative"])
        for g in genomes:
            cid = renum[assign[g]]
            w.writerow([g, cid, len(members[cid]), str(g == rep[cid])])

    big = sorted(((len(m), c) for c, m in members.items()), reverse=True)[:8]
    print(f"\nsingletons        : {sum(1 for m in members.values() if len(m) == 1)}")
    print(f"largest clusters  : {[(c, n) for n, c in big]}")
    print(f"genomes in n>=2   : {sum(len(m) for m in members.values() if len(m) >= 2)}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
