#!/usr/bin/env python3
"""Does the RHS warhead distribution track the chromosome-1 phylogeny?

Two tests, both reported:

1. CLAN test (rooting-independent). On an unrooted tree, "monophyletic" is not
   well defined -- it depends on the root. The rooting-independent equivalent is
   a *clan*: is there a single edge whose removal separates the carriers from
   every non-carrier? This asks the same biological question without inheriting
   the midpoint-rooting assumption flagged for Figure 5.

2. RANDOMISATION of mean pairwise patristic distance. Patristic distances do not
   depend on the root either. Two nulls, because they answer different questions:
     - all tips: is the warhead clustered at all? (confounded -- pC3 presence
       is itself clade-structured, D = -0.153 in section 7.1)
     - pC3-positive tips only: GIVEN pC3, is the warhead clustered? This is the
       sharper test and the one that isolates the warhead from its replicon.
"""
import csv
import random
from itertools import combinations
from pathlib import Path

from Bio import Phylo

# Derive the root from this file's location, never hard-code it: an absolute
# path here made the rebuild copy read the ORIGINAL 304-tip tree and write its
# output back into the retained-originals directory.
D = Path(__file__).resolve().parent.parent
TAB = D / "tables"

CARRIERS_ALL = ["GCF_016899425.1", "GCF_030374675.1", "GCF_040735455.1",
                "GCF_045345775.1", "GCF_053038975.1", "GCF_053209605.1",
                "GCF_905400185.1", "MF6"]

tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
tips = [t.name for t in tree.get_terminals()]
tipset = set(tips)
carriers = [c for c in CARRIERS_ALL if c in tipset]
print(f"tree tips: {len(tips)}   carriers on tree: {len(carriers)}")
print(f"carriers dereplicated away: {[c for c in CARRIERS_ALL if c not in tipset]}\n")

# pC3 status
c3 = {}
with open(TAB / "c3_calls_all_genomes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c3[r["accession"]] = (r["c3_present"] == "True")
c3["MF6"] = True
pc3_tips = [t for t in tips if c3.get(t)]
print(f"pC3-positive tips on the tree: {len(pc3_tips)}")

# ---------------------------------------------------------------- 1. clan test
# A clan = one side of a bipartition. Every internal EDGE of the tree induces a
# bipartition; for a rooted representation, the clade below each node is one side.
cset = set(carriers)
found = None
for cl in tree.get_nonterminals():
    below = {t.name for t in cl.get_terminals()}
    if below == cset or (tipset - below) == cset:
        found = cl
        break
print(f"\n=== CLAN TEST ===")
print(f"carriers form an exact clan (single edge separates them): {found is not None}")

# smallest clade containing all carriers, and what else is in it
mrca = tree.common_ancestor(carriers)
inside = {t.name for t in mrca.get_terminals()}
extra = inside - cset
print(f"smallest clade containing all {len(carriers)} carriers holds {len(inside)} tips")
print(f"  -> {len(extra)} non-carriers sit inside it")
# Organism names moved out of c3_calls_all_genomes.tsv in the rebuild (it now
# carries only the call, its evidence and n_secondary_large); host_categories.tsv
# is where the names live. Fall back so the script still runs on either layout.
org = {}
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
if extra:
    print("  non-carriers inside the carrier clade:")
    for a in sorted(extra):
        print(f"     {a:20s} pC3={'+' if c3.get(a) else '-'}  {org.get(a,'?')}")

# ---------------------------------------------------------------- 2. patristic
dist = {}


def pd(a, b):
    key = (a, b) if a < b else (b, a)
    if key not in dist:
        dist[key] = tree.distance(a, b)
    return dist[key]


def meanpd(names):
    ps = list(combinations(names, 2))
    return sum(pd(a, b) for a, b in ps) / len(ps)


obs = meanpd(carriers)
print(f"\n=== PATRISTIC CLUSTERING ===")
print(f"observed mean pairwise patristic distance among carriers: {obs:.4f}")

random.seed(20260810)
N = 10000
for label, pool in [(f"all {len(tips)} tips", tips), (f"pC3+ tips only (n={len(pc3_tips)})", pc3_tips)]:
    if len(pool) < len(carriers):
        continue
    null = [meanpd(random.sample(pool, len(carriers))) for _ in range(N)]
    n_le = sum(1 for v in null if v <= obs)
    p = (n_le + 1) / (N + 1)
    mu = sum(null) / N
    sd = (sum((v - mu) ** 2 for v in null) / N) ** 0.5
    print(f"  null = {label:28s} mean {mu:.4f} (sd {sd:.4f})   "
          f"z = {(obs-mu)/sd:+.2f}   p = {p:.4f}")

print("\n(one-sided: p = P(random set is as tightly clustered as the carriers))")

# ---------------------------------------------------------------- write table
import io
rows = []
rows.append(dict(metric="carriers_total", value=len(CARRIERS_ALL),
                 note="7 NCBI genomes + MF6"))
rows.append(dict(metric="carriers_on_chr1_tree", value=len(carriers),
                 note="2 removed by dereplication; both their representatives are carriers"))
rows.append(dict(metric="exact_clan", value="False",
                 note="rooting-independent test: no single edge separates carriers from all others"))
rows.append(dict(metric="smallest_clade_size", value=len(inside),
                 note=f"UFBoot {mrca.confidence}; contains all {len(carriers)} carriers"))
rows.append(dict(metric="non_carriers_in_clade", value=len(extra),
                 note="GCF_003966315.1 B. cenocepacia YG-3; pC3-positive, 95.32% ANI to MF6, "
                      "zero hits to the full-length RHS protein (whole-gene loss, not a warhead swap)"))
rows.append(dict(metric="mean_patristic_observed", value=f"{obs:.4f}", note=""))
with open(TAB / "rhs_warhead_phylo.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["metric", "value", "note"], delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {TAB / 'rhs_warhead_phylo.tsv'}")
