#!/usr/bin/env python3
"""Stage 9b - gene-content dendrogram of c3 replicons.

This is the pre-declared fallback (plan D7) for the case where the c3 core is too
small to support a core-gene phylogeny. It is NOT a phylogeny and must never be
presented as one: branch lengths are gene-content (Jaccard) distances, not
substitutions, and no evolutionary model is fitted. It answers "which c3s carry
similar gene repertoires", not "how are c3s related by descent".

Emits Newick so the same plotting code can consume it, with tip labels =
accessions, plus a summary of how much of the c3 pangenome is actually shared.
"""
import collections
import csv
import os
import sys

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist, squareform

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES = f"{W}/results"


def to_newick(node, labels):
    if node.is_leaf():
        return f"{labels[node.id]}:{max(node.dist, 0.0):.6f}"
    left = to_newick(node.get_left(), labels)
    right = to_newick(node.get_right(), labels)
    bl = max(node.dist - max(node.get_left().dist, node.get_right().dist), 0.0)
    return f"({left},{right}):{bl:.6f}"


def main():
    cls = {}
    for r in csv.DictReader(open(f"{RES}/secondary_replicon_clusters.tsv"),
                            delimiter="\t"):
        if r.get("replicon_class") == "c3" or r.get("is_c3") == "True":
            cls[(r["accession"], r["contig"])] = r

    rep_of = collections.defaultdict(set)
    with open(f"{W}/replicons/clu_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            acc, ctg, _ = mem.split("|", 2)
            if (acc, ctg) in cls:
                rep_of[(acc, ctg)].add(rep)

    keys = sorted(rep_of)
    print(f"confident c3 replicons: {len(keys)}")
    if len(keys) < 4:
        sys.exit("too few c3 replicons for a dendrogram")

    ogs = sorted({o for s in rep_of.values() for o in s})
    oidx = {o: i for i, o in enumerate(ogs)}
    M = np.zeros((len(keys), len(ogs)), dtype=bool)
    for i, k in enumerate(keys):
        for o in rep_of[k]:
            M[i, oidx[o]] = True
    print(f"c3 gene families (union): {len(ogs)}")

    n = len(keys)
    occ = M.sum(0)
    print("\nc3 pangenome structure (from confident c3 replicons only):")
    for thr, lab in ((1.00, "in ALL"), (0.95, ">=95%"), (0.90, ">=90%"),
                     (0.50, ">=50%"), (0.15, ">=15%")):
        print(f"  families {lab:>7} of {n} replicons: {(occ >= thr*n).sum():>6}")
    print(f"  singletons (1 replicon only):        {(occ == 1).sum():>6}")
    core = int((occ >= 0.95 * n).sum())
    print(f"\n-> core (>=95%) = {core} families")
    print("   A core-gene phylogeny needs roughly 50+ families; below that the"
          "\n   dendrogram below is the honest representation." if core < 50 else
          "   Core is large enough for a proper core-gene phylogeny; prefer that.")

    D = pdist(M, metric="jaccard")
    sq = squareform(D)
    iu = np.triu_indices(n, 1)
    print(f"\npairwise Jaccard distance: median={np.median(sq[iu]):.3f} "
          f"min={sq[iu].min():.3f} max={sq[iu].max():.3f}")

    Z = linkage(D, method="average")
    tree = to_tree(Z)
    labels = [k[0] for k in keys]
    nwk = to_newick(tree, labels) + ";"
    out = f"{W}/trees/c3_genecontent.treefile"
    os.makedirs(f"{W}/trees", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(nwk + "\n")
    print(f"\nwrote {out}")
    print("NOTE: gene-content dendrogram, NOT a phylogeny. Branch lengths are")
    print("      Jaccard distances; label it as such in any figure.")
    print("STAGE9B_DONE")


if __name__ == "__main__":
    sys.exit(main())
