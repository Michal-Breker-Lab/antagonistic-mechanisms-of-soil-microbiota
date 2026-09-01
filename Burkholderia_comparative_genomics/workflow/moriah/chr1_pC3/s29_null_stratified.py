#!/usr/bin/env python3
"""Randomisation tests for phylogenetic clustering, computed under TWO nulls.

Why two. The rebuild keeps every genome, so 233 of ~773 tips (30%) are one
99% ANI clone cluster -- the B. pseudomallei / B. mallei complex -- which is
sequencing effort on a biothreat pathogen, not biology. That distorts the
unrestricted null, and the distortion is CONSERVATIVE, not anti-conservative:

  p is defined as P(a random tip set is at least as tightly clustered as the
  carriers), so only the null's LEFT tail matters. A large clone cluster fills
  that tail with all-clone draws whose mean pairwise distance is ~0 -- tighter
  than any real carrier set can be. The floor this puts under p is close to
  P(every drawn tip comes from the clone cluster) = C(233,8)/C(771,8) ~= 7e-5
  for an 8-carrier set on the full tree. The null's *mean* also inflates
  (cross-clade pairs become common), but the mean is not what p is read from.

  Empirically, over clone fractions from 10% to 80% of the tree, p_unrestricted
  was >= p_stratified at every point, and the gap widened with clone fraction
  (see selftest 2). So the failure mode to guard against is a real clustering
  result being erased by an artefactual null floor -- NOT a false positive.

  1. unrestricted -- tips drawn uniformly from the pool, as in the original
     analysis. Comparable to the published numbers.
  2. clone-stratified -- draws match the observed carrier set's cluster
     composition: the same number of distinct clusters, contributing the same
     number of tips each. A clone cluster can therefore never contribute more
     tips to a null draw than the carriers actually drew from it, which removes
     the all-clone floor.

Both p-values are reported side by side. Where they disagree the stratified one
is the interpretable test, and the report says so explicitly rather than
printing whichever number is smaller.

Distances are computed once into a full patristic matrix (depth_a + depth_b -
2 * depth_lca) rather than per-call, because 10,000 replicates x 28 pairs would
otherwise dominate runtime.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

import numpy as np


# --------------------------------------------------------------- tree distances
def patristic_matrix(tree) -> tuple[list[str], np.ndarray]:
    """Full tip-by-tip patristic distance matrix from a Bio.Phylo tree.

    Uses d(a,b) = depth(a) + depth(b) - 2*depth(lca(a,b)). Each pair's LCA depth
    is written exactly once, when the post-order walk reaches the node that first
    joins the two tips, so the whole matrix costs O(n^2) writes rather than
    O(n^2) tree traversals.
    """
    tips = [t.name for t in tree.get_terminals()]
    idx = {n: i for i, n in enumerate(tips)}
    n = len(tips)

    depth: dict[int, float] = {}
    for clade, d in _depths(tree.root, 0.0):
        depth[id(clade)] = d

    tipdepth = np.zeros(n)
    for t in tree.get_terminals():
        tipdepth[idx[t.name]] = depth[id(t)]

    lca = np.zeros((n, n))
    for clade in tree.find_clades(order="postorder"):
        kids = clade.clades
        if len(kids) < 2:
            continue
        sets = [np.array([idx[t.name] for t in k.get_terminals()], dtype=int)
                for k in kids]
        d = depth[id(clade)]
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                lca[np.ix_(sets[i], sets[j])] = d
                lca[np.ix_(sets[j], sets[i])] = d

    D = tipdepth[:, None] + tipdepth[None, :] - 2.0 * lca
    np.fill_diagonal(D, 0.0)
    return tips, D


def _depths(clade, acc):
    """Yield (clade, root-to-clade distance) for every clade, iteratively."""
    stack = [(clade, acc)]
    while stack:
        c, d = stack.pop()
        yield c, d
        for k in c.clades:
            stack.append((k, d + (k.branch_length or 0.0)))


def mean_pairwise(D: np.ndarray, ix: np.ndarray) -> float:
    """Mean pairwise distance among the tips indexed by ix."""
    sub = D[np.ix_(ix, ix)]
    k = len(ix)
    return float(sub.sum() / (k * (k - 1)))


# ------------------------------------------------------------------------ nulls
def unrestricted_null(D, pool_ix, k, nreps, rng) -> np.ndarray:
    out = np.empty(nreps)
    for r in range(nreps):
        out[r] = mean_pairwise(D, rng.choice(pool_ix, size=k, replace=False))
    return out


def stratified_null(D, pool_by_cluster, composition, nreps, rng) -> np.ndarray:
    """Null draws matching the observed carrier set's cluster composition.

    `composition` is the multiset of per-cluster carrier counts, e.g. [2,1,1,1,1,1]
    means the carriers came from six distinct clusters, one of which contributed
    two genomes. Each replicate draws that many distinct clusters -- each large
    enough to supply its assigned count -- and samples that many tips from each.
    Multiplicities are filled largest-first so the scarce large clusters are
    placed while the choice is still unconstrained.
    """
    cluster_ids = list(pool_by_cluster)
    sizes = np.array([len(pool_by_cluster[c]) for c in cluster_ids])
    comp = sorted(composition, reverse=True)
    out = np.empty(nreps)

    for r in range(nreps):
        chosen: list[int] = []
        used: set[int] = set()
        for need in comp:
            cand = np.flatnonzero((sizes >= need))
            cand = np.array([c for c in cand if c not in used])
            if len(cand) == 0:
                raise ValueError(
                    f"no cluster left with >= {need} tips; the carrier set's "
                    "cluster composition cannot be matched from this pool")
            pick = int(rng.choice(cand))
            used.add(pick)
            members = pool_by_cluster[cluster_ids[pick]]
            chosen.extend(rng.choice(members, size=need, replace=False).tolist())
        out[r] = mean_pairwise(D, np.array(chosen, dtype=int))
    return out


def dual_test(D, tips, carriers, pool, clusters, nreps=10000, seed=20260819):
    """Run both nulls for one carrier set. Returns a dict of results."""
    idx = {n: i for i, n in enumerate(tips)}
    car_ix = np.array([idx[c] for c in carriers if c in idx], dtype=int)
    if len(car_ix) < 3:
        raise ValueError(f"need >=3 carriers on the tree, got {len(car_ix)}")
    pool_ix = np.array([idx[p] for p in pool if p in idx], dtype=int)
    obs = mean_pairwise(D, car_ix)

    rng = np.random.default_rng(seed)
    un = unrestricted_null(D, pool_ix, len(car_ix), nreps, rng)

    # composition of the OBSERVED carriers across clone clusters
    car_names = [tips[i] for i in car_ix]
    comp = collections.Counter(clusters[c] for c in car_names)
    pool_by_cluster: dict[str, np.ndarray] = collections.defaultdict(list)
    for i in pool_ix:
        pool_by_cluster[clusters[tips[i]]].append(i)
    pool_by_cluster = {c: np.array(v, dtype=int) for c, v in pool_by_cluster.items()}

    rng2 = np.random.default_rng(seed)
    st = stratified_null(D, pool_by_cluster, list(comp.values()), nreps, rng2)

    def _p(null):
        return (int((null <= obs).sum()) + 1) / (nreps + 1)

    return {
        "n_carriers": len(car_ix),
        "n_pool": len(pool_ix),
        "n_carrier_clusters": len(comp),
        "carrier_composition": sorted(comp.values(), reverse=True),
        "observed": obs,
        "p_unrestricted": _p(un),
        "null_unrestricted_mean": float(un.mean()),
        "null_unrestricted_sd": float(un.std()),
        "p_stratified": _p(st),
        "null_stratified_mean": float(st.mean()),
        "null_stratified_sd": float(st.std()),
    }


def verdict(res: dict, alpha: float = 0.05) -> str:
    """Interpret the pair of p-values.

    The stratified test is the interpretable one; the unrestricted test is kept
    for comparability with the published, dereplicated analysis.
    """
    pu, ps = res["p_unrestricted"], res["p_stratified"]
    if pu < alpha and ps < alpha:
        return "clustered: significant under both nulls"
    if ps < alpha <= pu:
        return ("clustered after clone correction; the unrestricted null is "
                "floored by all-clone draws and understates the result")
    if pu < alpha <= ps:
        return ("NOT supported once clone structure is respected - the "
                "unrestricted result is carried by clone sampling")
    return "not clustered under either null"


# ------------------------------------------------------------------- self-tests
def selftest() -> int:
    """Three checks: degeneracy, direction of the clone bias, and its mechanism."""
    from io import StringIO
    from math import comb

    from Bio import Phylo

    def comb_tree(n_clone: int, n_div: int = None):
        """Clade A: diverse shallow tips. Clade B: a deep near-identical comb."""
        n_div = n_div if n_div is not None else 100 - n_clone
        a = ",".join(f"a{i}:{0.02 + 0.001 * i:.4f}" for i in range(n_div))
        b = ",".join(f"b{i}:0.0005" for i in range(n_clone))
        tree = Phylo.read(StringIO(f"(({a}):0.05,({b}):1.20);"), "newick")
        tips, D = patristic_matrix(tree)
        clusters = {t: ("cloneB" if t.startswith("b") else t) for t in tips}
        return tips, D, clusters, n_div

    ok = []

    print("=== self-test 1: no clone structure -> the two nulls must agree ===")
    nw = "(" + ",".join(f"(t{i}:0.1,t{i+1}:0.1):0.1" for i in range(0, 40, 2)) + ");"
    tree = Phylo.read(StringIO(nw), "newick")
    tips, D = patristic_matrix(tree)
    clusters = {t: t for t in tips}
    res = dual_test(D, tips, tips[:8], tips, clusters, nreps=4000, seed=1)
    diff = abs(res["p_unrestricted"] - res["p_stratified"])
    print(f"  p_unrestricted={res['p_unrestricted']:.4f}  "
          f"p_stratified={res['p_stratified']:.4f}  |diff|={diff:.4f}")
    ok.append(diff < 0.05)
    print(f"  {'PASS' if ok[-1] else 'FAIL'} (agreement within Monte-Carlo error)")

    print("\n=== self-test 2: the unrestricted null is CONSERVATIVE, and more so "
          "as the clone cluster grows ===")
    print(f"  {'cloneN':>7} {'p_unrestricted':>15} {'p_stratified':>13}  ordering")
    gaps = []
    for nb in (10, 25, 40, 60, 80):
        tips2, D2, cl2, ndiv = comb_tree(nb)
        r = dual_test(D2, tips2, [f"a{i}" for i in range(8)], tips2, cl2,
                      nreps=4000, seed=7)
        good = r["p_unrestricted"] >= r["p_stratified"]
        gaps.append(r["p_unrestricted"] - r["p_stratified"])
        print(f"  {nb:>7} {r['p_unrestricted']:>15.4f} {r['p_stratified']:>13.4f}"
              f"  {'p_un >= p_st  OK' if good else 'UNEXPECTED'}")
        ok.append(good)
    widening = gaps[-1] > gaps[0]
    print(f"  gap widens with clone fraction: {widening}")
    ok.append(widening)
    print(f"  {'PASS' if all(ok[-6:]) else 'FAIL'}")

    print("\n=== self-test 3: the mechanism is the all-clone draw ===")
    # p_unrestricted should sit at or above P(all 8 tips drawn from the clone
    # cluster), because every such draw has ~zero mean pairwise distance and so
    # always counts toward p.
    nb = 60
    tips3, D3, cl3, ndiv = comb_tree(nb)
    r3 = dual_test(D3, tips3, [f"a{i}" for i in range(8)], tips3, cl3,
                   nreps=8000, seed=11)
    floor = comb(nb, 8) / comb(nb + ndiv, 8)
    print(f"  P(all 8 drawn from the clone cluster) = {floor:.5f}")
    print(f"  p_unrestricted                        = {r3['p_unrestricted']:.5f}")
    print(f"  p_stratified                          = {r3['p_stratified']:.5f}")
    near = r3["p_unrestricted"] >= floor * 0.5
    ok.append(near)
    print(f"  {'PASS' if near else 'FAIL'} (p_unrestricted is floored by that "
          "probability)")
    print(f"  verdict string: {verdict(r3)}")

    print(f"\n{'ALL PASS' if all(ok) else 'FAILURES PRESENT'}")
    return 0 if all(ok) else 1


# -------------------------------------------------------------------------- CLI
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tree", type=Path)
    ap.add_argument("--clusters", type=Path, help="clone_cluster.tsv from s28")
    ap.add_argument("--carriers", type=Path, help="one accession per line")
    ap.add_argument("--pool", type=Path, default=None,
                    help="restrict the null pool (e.g. pC3-positive tips only)")
    ap.add_argument("--label", default="carriers")
    ap.add_argument("--nreps", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260819)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()

    for req in ("tree", "clusters", "carriers"):
        if getattr(args, req) is None:
            ap.error(f"--{req} is required unless --selftest")

    from Bio import Phylo
    tree = Phylo.read(str(args.tree), "newick")
    tips, D = patristic_matrix(tree)

    clusters = {}
    with open(args.clusters) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            clusters[row["accession"]] = row["cluster_id"]
    missing = [t for t in tips if t not in clusters]
    if missing:
        raise SystemExit(f"{len(missing)} tips have no clone-cluster label, "
                         f"first few: {missing[:5]}")

    carriers = [ln.strip() for ln in open(args.carriers) if ln.strip()]
    pool = ([ln.strip() for ln in open(args.pool) if ln.strip()]
            if args.pool else tips)

    res = dual_test(D, tips, carriers, pool, clusters,
                    nreps=args.nreps, seed=args.seed)
    res["label"] = args.label
    res["verdict"] = verdict(res)

    for k, v in res.items():
        print(f"  {k:26s} {v}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        newfile = not args.out.exists()
        with open(args.out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(res), delimiter="\t")
            if newfile:
                w.writeheader()
            w.writerow(res)
        print(f"\nappended to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
