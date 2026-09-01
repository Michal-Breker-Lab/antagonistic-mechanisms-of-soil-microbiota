#!/usr/bin/env python3
"""Render-time collapsing of clone clades for the 771-tip tree figures.

Nothing is removed from the ANALYSIS -- every genome stays in the alignment, the
tree and every statistic. This module only makes the tree DRAWABLE: a 233-tip
near-identical B. pseudomallei comb is 30% of the tips and none of the
information, so in the main figure panel it becomes one triangle labelled with
its member count, while the fully expanded tree ships as a supplementary file.

A clone cluster is collapsed only where it is monophyletic on the tree. If a
cluster is scattered, each maximal single-cluster clade collapses separately and
the rest stay as ordinary tips -- the figure never implies a monophyly the tree
does not support. `collapse_clone_clades` reports both, so the caption can state
how many clusters were non-monophyletic rather than hiding it.
"""
from __future__ import annotations

import collections
import copy


def _tipnames(clade) -> list[str]:
    return [t.name for t in clade.get_terminals()]


def collapse_clone_clades(tree, clusters: dict[str, str], min_n: int = 3,
                          label_fmt: str = "{cluster} (n={n})"):
    """Collapse maximal monophyletic single-cluster clades into one tip each.

    Args:
        tree: a Bio.Phylo tree (modified on a deep copy, not in place).
        clusters: accession -> clone cluster id, covering every tip.
        min_n: clades smaller than this are left expanded.
        label_fmt: label for a collapsed tip.

    Returns:
        (new_tree, info) where info holds:
          collapsed  : list of dicts (cluster, n, label, members, depth)
          kept_tips  : tips left as ordinary tips
          scattered  : cluster ids that are NOT monophyletic on this tree
    """
    t = copy.deepcopy(tree)
    tips = _tipnames(t)
    missing = [x for x in tips if x not in clusters]
    if missing:
        raise ValueError(f"{len(missing)} tips lack a cluster label, "
                         f"e.g. {missing[:5]}")

    total = collections.Counter(clusters[x] for x in tips)

    # Walk pre-order; the first clade whose tips are all one cluster is maximal.
    collapsed: list[dict] = []
    seen_members: set[str] = set()

    def walk(clade):
        names = _tipnames(clade)
        if not names:
            return
        cl = {clusters[n] for n in names}
        if len(cl) == 1 and len(names) >= min_n:
            cid = cl.pop()
            collapsed.append({
                "cluster": cid,
                "n": len(names),
                "label": label_fmt.format(cluster=cid, n=len(names)),
                "members": sorted(names),
                "clade": clade,
            })
            seen_members.update(names)
            return                       # maximal: do not descend further
        for k in clade.clades:
            walk(k)

    walk(t.root)

    # Rewrite each collapsed clade into a single terminal.
    for rec in collapsed:
        clade = rec.pop("clade")
        # keep the clade's own depth by giving the new tip the mean tip length
        lens = []
        for tip in clade.get_terminals():
            lens.append(t.distance(clade, tip))
        clade.clades = []
        clade.name = rec["label"]
        clade.branch_length = (clade.branch_length or 0.0)
        rec["mean_tip_depth"] = sum(lens) / len(lens) if lens else 0.0

    # A cluster is scattered if its members did not all end up in one collapse.
    collapsed_by_cluster = collections.Counter(r["cluster"] for r in collapsed)
    covered = collections.Counter()
    for r in collapsed:
        covered[r["cluster"]] += r["n"]
    scattered = sorted(
        cid for cid, n in total.items()
        if n >= min_n and (collapsed_by_cluster[cid] > 1 or covered[cid] < n)
    )

    kept = [x for x in _tipnames(t) if x not in {r["label"] for r in collapsed}]
    info = {
        "collapsed": collapsed,
        "n_collapsed_clades": len(collapsed),
        "n_tips_hidden": sum(r["n"] for r in collapsed),
        "kept_tips": kept,
        "n_tips_drawn": len(_tipnames(t)),
        "scattered": scattered,
    }
    return t, info


def simulate_full_tree(tree, membership: dict[str, list[str]],
                       clone_len: float = 1e-4):
    """Preview helper: expand each representative tip into its clone cluster.

    Used to rehearse the 771-tip figures against the existing 304-tip tree while
    the real full-size tree has not been built yet. Not used for any published
    number -- the branch lengths inside an expanded cluster are invented.
    """
    t = copy.deepcopy(tree)
    for tip in t.get_terminals():
        members = membership.get(tip.name)
        if not members or len(members) < 2:
            continue
        from Bio.Phylo.BaseTree import Clade
        kids = [Clade(branch_length=clone_len, name=m) for m in sorted(members)]
        tip.clades = kids
        tip.name = None
    return t


if __name__ == "__main__":
    import csv
    import sys
    from pathlib import Path

    from Bio import Phylo

    D = Path(__file__).resolve().parent.parent
    TAB = D / "tables"

    print("=== rehearsal: expand the 304-tip tree to full size, then collapse ===")
    tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
    print(f"published tree tips              : {len(tree.get_terminals())}")

    clusters, by_cluster = {}, collections.defaultdict(list)
    for r in csv.DictReader(open(TAB / "derep_cluster_membership.tsv"),
                            delimiter="\t"):
        clusters[r["accession"]] = r["cluster_id"]
        by_cluster[r["cluster_id"]].append(r["accession"])

    tip_names = {t.name for t in tree.get_terminals()}
    membership = {}
    for cid, members in by_cluster.items():
        reps = [m for m in members if m in tip_names]
        if reps:
            membership[reps[0]] = members

    full = simulate_full_tree(tree, membership)
    ftips = [t.name for t in full.get_terminals()]
    print(f"simulated full tree tips         : {len(ftips)}")

    cl = dict(clusters)
    for n in ftips:
        cl.setdefault(n, n)
    small, info = collapse_clone_clades(full, cl, min_n=3)
    print(f"collapsed clades                 : {info['n_collapsed_clades']}")
    print(f"tips hidden inside triangles     : {info['n_tips_hidden']}")
    print(f"tips actually drawn              : {info['n_tips_drawn']}")
    print(f"non-monophyletic clusters        : {len(info['scattered'])}")
    big = sorted(info["collapsed"], key=lambda r: -r["n"])[:6]
    print("largest collapsed clades         : "
          + ", ".join(f"{r['label']}" for r in big))
    sys.exit(0)
