#!/usr/bin/env python3
"""Does each of the 10 toxin loci track the chromosome-1 phylogeny?

Generalises warhead_phylo_test.py from one locus to all of them, with two
pre-declared guards that the single-locus version did not need:

1. TESTABILITY WINDOW. A locus is tested only if 3 <= carriers <= n_tips - 3.
   Below 3 carriers there is no pairwise distance distribution to speak of;
   within 3 of every tip there is no null distribution either, because almost
   every tip is a carrier. Excluded loci are listed with their reason rather
   than dropped silently.

   The bound is expressed in CARRIERS AND NON-CARRIERS, not as a literal tip
   count: it was 3-301 when the tree had 304 tips, and that constant did not
   move when the rebuilt tree grew to 763 (no dereplication). Left stale it
   discarded four testable loci -- MF6_003843 at 444 carriers is 58% of tips,
   nowhere near "almost every tip" -- and shrank the BH denominator from 10 to
   6, which is a correction applied to the wrong family.

2. MULTIPLE TESTING. One randomisation per qualifying locus, Benjamini-Hochberg
   across them. The previous session ran a single locus and needed no correction.

Both tests are ROOTING-INDEPENDENT by design. The chromosome-1 tree is unrooted
(IQ-TREE writes a trifurcating root) and its midpoint rooting is a display
convention only, so a monophyly test would import an assumption the data do not
support. The clan test and patristic distances do not.

Nulls: all 304 tips, and pC3-positive tips only. The second is sharper because
pC3 presence is itself clade-structured (D = -0.153, section 7.1), so clustering
against all tips could simply be re-detecting pC3's own distribution.
"""
import csv
import random
from collections import defaultdict

from pathlib import Path

import numpy as np
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

D = Path(__file__).resolve().parent.parent
TAB = D / "tables"
N_PERM = 10000
# Scales with the tree; see the TESTABILITY WINDOW note above. MAX_CARRIERS is
# set once `tips` is known.
MIN_CARRIERS = 3

# ------------------------------------------------------------------ inputs
tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
tips = [t.name for t in tree.get_terminals()]
MAX_CARRIERS = len(tips) - MIN_CARRIERS
tipset = set(tips)
print(f"tree tips: {len(tips)}")

c3, org = {}, {}
with open(TAB / "c3_calls_all_genomes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c3[r["accession"]] = (r["c3_present"] == "True")
# Organism names moved out of c3_calls_all_genomes.tsv in the rebuild (it now
# carries only the call, its evidence and n_secondary_large); they live in
# host_categories.tsv. Fall back so the script runs on either layout.
for src in ("host_categories.tsv", "c3_calls_all_genomes.tsv"):
    f = TAB / src
    if not f.exists():
        continue
    with open(f, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if "organism_name" not in (rd.fieldnames or []):
            continue
        for r in rd:
            org.setdefault(r["accession"], r["organism_name"])
    break
c3["MF6"] = True
org["MF6"] = "Burkholderia sola MF6"
pc3_tips = [t for t in tips if c3.get(t)]
print(f"pC3-positive tips: {len(pc3_tips)}\n")

carriers_all = defaultdict(set)
with open(TAB / "toxin12_carriers.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["tier1_carrier"] == "True":
            carriers_all[r["locus"]].add(r["accession"])

# Dereplication clusters, reconstructed from the retained skani all-vs-all edges
# with the same 99% ANI / 50% AF single-linkage rule s5_derep.py used. s5 wrote
# only the representatives, not the membership, so this is rebuilt rather than
# read -- it is what makes "was this carrier's lineage dropped from the tree?"
# answerable exactly instead of by inspecting species names.
cluster_of = {}
with open(TAB / "derep_cluster_membership.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cluster_of[r["accession"]] = r["cluster_id"]
tip_clusters = {cluster_of.get(t) for t in tips}

# ------------------------------------------------------------------ distances
# All-pairs patristic distances in one pass. Bio.Phylo's tree.distance() walks
# the tree per call, which is far too slow for 10,000 permutations across ten
# loci (the first attempt timed out). Instead: every pair of tips has its MRCA at
# exactly one internal node, so walking each internal node once and filling in
# all cross-subtree pairs there gives the full matrix in O(n^2).
tip_idx = {t: i for i, t in enumerate(tips)}
n_t = len(tips)
DIST = np.zeros((n_t, n_t), dtype=float)

depth = {id(tree.root): 0.0}
for cl in tree.find_clades(order="preorder"):
    for c in cl.clades:
        depth[id(c)] = depth[id(cl)] + max(c.branch_length or 0.0, 0.0)
TIP_DEPTH = np.array([depth[id(t)] for t in tree.get_terminals()])

below = {}
for cl in tree.find_clades(order="postorder"):
    if cl.is_terminal():
        below[id(cl)] = [tip_idx[cl.name]]
        continue
    groups = [np.array(below[id(c)]) for c in cl.clades]
    dn = depth[id(cl)]
    for gi in range(len(groups)):
        for gj in range(gi + 1, len(groups)):
            a, b = groups[gi], groups[gj]
            block = (TIP_DEPTH[a][:, None] + TIP_DEPTH[b][None, :]) - 2 * dn
            DIST[np.ix_(a, b)] = block
            DIST[np.ix_(b, a)] = block.T
    below[id(cl)] = [int(i) for g in groups for i in g]


def meanpd(names):
    ix = np.array([tip_idx[x] for x in names])
    sub = DIST[np.ix_(ix, ix)]
    k = len(ix)
    return sub.sum() / (k * (k - 1))


# The ten LOCI. The two "_Cterm" entries in toxin12_carriers.tsv are the
# C-terminal tips of MF6_003684 and MF6_004284 and are NOT separate loci -- they
# are reported below as a within-locus contrast and are excluded from the
# multiple-testing correction, otherwise those two loci would be counted twice.
LOCI = [k for k in sorted(carriers_all) if not k.endswith("_Cterm")]
CTERM = [k for k in sorted(carriers_all) if k.endswith("_Cterm")]

# ------------------------------------------------------------------ per locus
random.seed(20260812)
rows, tested_p, tested_idx = [], [], []
for locus in LOCI + CTERM:
    call = carriers_all[locus]
    on_tree = sorted(call & tipset)
    off_tree = sorted(call - tipset)
    # A carrier removed by dereplication is only genuinely MISSING from the tree
    # if no carrier from its own 99%-ANI cluster survived onto the tree. If a
    # cluster-mate did survive and is itself a carrier, the lineage is still
    # represented and the tree does not understate the distribution.
    rescued = [a for a in off_tree
               if any(cluster_of.get(x) == cluster_of.get(a) for x in on_tree)]
    lost = [a for a in off_tree if a not in rescued]

    row = dict(locus=locus, n_carriers_total=len(call),
               n_carriers_on_tree=len(on_tree),
               n_off_tree=len(off_tree),
               n_off_tree_rescued=len(rescued),
               n_off_tree_lost=len(lost),
               off_tree_lost=";".join(lost))

    if not (MIN_CARRIERS <= len(on_tree) <= MAX_CARRIERS):
        row.update(tested=False,
                   reason=(f"{len(on_tree)} carriers on tree, outside the "
                           f"pre-declared {MIN_CARRIERS}-{MAX_CARRIERS} window"),
                   exact_clan="", clade_size="", clade_ufboot="",
                   mean_patristic="", z_all="", p_all="", z_pc3="", p_pc3="")
        rows.append(row)
        print(f"{locus:<22} SKIPPED - {row['reason']}")
        continue

    cset = set(on_tree)
    clan = any({t.name for t in cl.get_terminals()} in (cset,)
               or (tipset - {t.name for t in cl.get_terminals()}) == cset
               for cl in tree.get_nonterminals())
    mrca = tree.common_ancestor(on_tree)
    inside = {t.name for t in mrca.get_terminals()}
    obs = meanpd(on_tree)

    res, note = {}, ""
    for key, pool in (("all", tips), ("pc3", pc3_tips)):
        # The pC3-only null cannot be built when a locus has more carriers than
        # there are pC3-positive tips -- you cannot draw 204 tips from a pool of
        # 140. That is reported as "undefined", not as a missing value, and such
        # loci are held out of the BH correction rather than entering it as NaN.
        if len(pool) < len(on_tree):
            res[key] = (None, None)
            note = (f"pC3 null undefined: {len(on_tree)} carriers exceed the "
                    f"{len(pc3_tips)} pC3-positive tips")
            continue
        null = [meanpd(random.sample(pool, len(on_tree))) for _ in range(N_PERM)]
        mu = sum(null) / N_PERM
        sd = (sum((v - mu) ** 2 for v in null) / N_PERM) ** 0.5
        p = (sum(1 for v in null if v <= obs) + 1) / (N_PERM + 1)
        res[key] = ((obs - mu) / sd if sd else float("nan"), p)

    def fmt(v, spec):
        return "" if v is None else format(v, spec)

    row.update(tested=True, reason=note,
               exact_clan=clan, clade_size=len(inside),
               clade_ufboot=mrca.confidence,
               n_noncarriers_in_clade=len(inside - cset),
               noncarriers_in_clade=";".join(sorted(inside - cset)),
               mean_patristic=f"{obs:.4f}",
               z_all=fmt(res["all"][0], "+.2f"), p_all=fmt(res["all"][1], ".4f"),
               z_pc3=fmt(res["pc3"][0], "+.2f"), p_pc3=fmt(res["pc3"][1], ".4f"))
    # only the ten loci with a DEFINED pC3-null p-value enter the correction
    if locus in LOCI and res["pc3"][1] is not None:
        tested_idx.append(len(rows))
        tested_p.append(res["pc3"][1])
    rows.append(row)
    print(f"{locus:<22} n={len(on_tree):>3} clan={str(clan):<5} "
          f"clade={len(inside):>3} (UFBoot {mrca.confidence}) "
          f"mean_pd={obs:.4f}  z_all={fmt(res['all'][0], '+.2f'):>6} "
          f"p={fmt(res['all'][1], '.4f'):>6}  "
          f"z_pC3={fmt(res['pc3'][0], '+.2f') or 'undef':>6} "
          f"p={fmt(res['pc3'][1], '.4f') or 'undef':>6}"
          + ("   [C-term contrast, not in BH]" if locus in CTERM else ""))

# ------------------------------------------------------------------ BH
for r in rows:
    r["p_pc3_bh"] = ""
    r["significant_bh_0.05"] = ""
if tested_p:
    rej, q, _, _ = multipletests(tested_p, method="fdr_bh")
    for i, qq, rj in zip(tested_idx, q, rej):
        rows[i]["p_pc3_bh"] = f"{qq:.4f}"
        rows[i]["significant_bh_0.05"] = bool(rj)

cols = ["locus", "n_carriers_total", "n_carriers_on_tree", "n_off_tree",
        "n_off_tree_rescued", "n_off_tree_lost", "off_tree_lost", "tested",
        "reason", "exact_clan", "clade_size", "clade_ufboot",
        "n_noncarriers_in_clade", "noncarriers_in_clade", "mean_patristic",
        "z_all", "p_all", "z_pc3", "p_pc3", "p_pc3_bh", "significant_bh_0.05"]
with open(TAB / "toxin12_phylo.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})

print(f"\nBH-corrected across the {len(tested_p)} of {len(LOCI)} loci with a "
      f"defined pC3-null p-value\n")
for i, r in enumerate(rows):
    if r["locus"] in CTERM:
        continue
    if i in tested_idx:
        print(f"  {r['locus']:<22} p={r['p_pc3']}  p_BH={r['p_pc3_bh']}  "
              f"significant={r['significant_bh_0.05']}")
    else:
        print(f"  {r['locus']:<22} not corrected - {r['reason'] or 'not tested'}")
print("\nOne-sided: p = P(a random set is at least as tightly clustered).")
print("A large p therefore means over-dispersed relative to that null, not "
      "'no signal'.")
print(f"\nwrote {TAB / 'toxin12_phylo.tsv'}")
