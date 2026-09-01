#!/usr/bin/env python3
"""Stage 4c - is c3 one evolutionary entity, or merely a size class?

WHY NOT NAIVE CLUSTERING. Clustering all secondary replicons by raw gene-content
Jaccard answers the wrong question: gene content is dominated by which ORGANISM a
replicon came from, not by which replicon it is. Chromosome 2 and c3 of the same
strain share far more genes with each other than two chromosome 2s from different
genera do, so a global clustering recovers taxonomy and says nothing about
replicon identity. An earlier version of this script did exactly that and
concluded "c3 does not separate from chr2" -- an artifact of the statistic.

Two taxonomy-controlled tests are used instead, both restricted to genomes that
carry exactly TWO large secondary replicons, where c2 (larger) and c3 (smaller)
are unambiguous by position:

  TEST A - orthogroup residence consistency. For each orthogroup seen in many
  such genomes, what fraction of the time does it sit on c3? If c3 is a coherent
  replicon, this is strongly BIMODAL (genes are reliably c3-resident or reliably
  c2-resident). If "chromosome 3" were merely a size class, residence would be
  near-random and the distribution unimodal around the base rate. Compared
  against a within-genome label-shuffled null.

  TEST B - paired cross-genome similarity. For genome pairs (A, B), is
  J(c3_A, c3_B) > J(c3_A, c2_B)? Taxonomy is held constant because both
  comparisons involve the same pair of genomes.
"""
import collections
import csv
import os
import random
import sys

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES = f"{W}/results"
MF6_C3 = ("MF6", "cluster_003_consensus")
MIN_GENOMES_PER_OG = 10
random.seed(0)
np.random.seed(0)


def main():
    rep_of = {}
    with open(f"{W}/replicons/clu_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            acc, ctg, _ = mem.split("|", 2)
            rep_of.setdefault((acc, ctg), set()).add(rep)
    replicons = sorted(rep_of)
    print(f"secondary replicons: {len(replicons)}")

    ogs = sorted({o for s in rep_of.values() for o in s})
    oidx = {o: i for i, o in enumerate(ogs)}
    print(f"orthogroups: {len(ogs)}")
    M = np.zeros((len(replicons), len(ogs)), dtype=bool)
    for i, r in enumerate(replicons):
        for o in rep_of[r]:
            M[i, oidx[o]] = True
    ridx = {r: i for i, r in enumerate(replicons)}

    types = {(r["accession"], r["contig"]): r for r in
             csv.DictReader(open(f"{RES}/replicon_types.tsv"), delimiter="\t")}
    cen = {r["accession"]: r for r in
           csv.DictReader(open(f"{RES}/replicon_census.tsv"), delimiter="\t")}

    bygen = collections.defaultdict(list)
    for r in replicons:
        bygen[r[0]].append(r)
    for acc in bygen:
        bygen[acc].sort(key=lambda x: -int(types[x]["length"]))
    sec_rank = {x: n for acc, rl in bygen.items() for n, x in enumerate(rl, 1)}

    # genomes with exactly two large secondary replicons: c2 = larger, c3 = smaller
    pairs = {acc: rl for acc, rl in bygen.items() if len(rl) == 2}
    print(f"genomes with exactly 2 large secondary replicons (training set): {len(pairs)}")
    if len(pairs) < 10:
        sys.exit("too few 2-secondary genomes to test coherence")

    # ------------------------------------------------------------------
    # TEST A - orthogroup residence consistency
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("TEST A - orthogroup residence consistency (c3 vs c2)")
    print("=" * 72)
    on_c3 = collections.Counter()
    on_any = collections.Counter()
    for acc, (c2, c3) in pairs.items():
        s2, s3 = rep_of[c2], rep_of[c3]
        for o in s2 | s3:
            on_any[o] += 1
            if o in s3:
                on_c3[o] += 1
    common = [o for o in on_any if on_any[o] >= MIN_GENOMES_PER_OG]
    frac = np.array([on_c3[o] / on_any[o] for o in common])
    print(f"orthogroups seen in >={MIN_GENOMES_PER_OG} of the {len(pairs)} "
          f"2-secondary genomes: {len(common)}")
    base = sum(on_c3.values()) / max(sum(on_any.values()), 1)
    print(f"base rate (overall fraction of gene occurrences on c3): {base:.3f}")

    def summarise(f, label):
        ext = np.mean((f < 0.1) | (f > 0.9))
        mid = np.mean((f > 0.35) & (f < 0.65))
        print(f"  {label:<22} extreme(<0.1 or >0.9)={100*ext:>5.1f}%   "
              f"intermediate(0.35-0.65)={100*mid:>5.1f}%")
        return ext, mid

    obs_ext, obs_mid = summarise(frac, "observed")

    # null: shuffle the c2/c3 label within each genome
    null_ext = []
    for _ in range(20):
        n_c3 = collections.Counter()
        n_any = collections.Counter()
        for acc, rl in pairs.items():
            a, b = rl if random.random() < 0.5 else rl[::-1]
            sa, sb = rep_of[a], rep_of[b]
            for o in sa | sb:
                n_any[o] += 1
                if o in sb:
                    n_c3[o] += 1
        f = np.array([n_c3[o] / n_any[o] for o in common if n_any[o]])
        null_ext.append(np.mean((f < 0.1) | (f > 0.9)))
    print(f"  {'shuffled null':<22} extreme={100*np.mean(null_ext):>5.1f}% "
          f"(sd {100*np.std(null_ext):.1f})")
    print(f"\n  VERDICT: {100*obs_ext:.1f}% of shared orthogroups are near-exclusively "
          f"resident on one replicon\n           vs {100*np.mean(null_ext):.1f}% expected "
          f"if residence were random.")
    coherent = obs_ext > np.mean(null_ext) + 5 * (np.std(null_ext) + 1e-9)
    print(f"           -> c3 IS a coherent replicon with its own gene complement"
          if coherent else
          "           -> residence looks random; c3 may be only a size class")

    # ------------------------------------------------------------------
    # TEST B - paired cross-genome similarity, taxonomy held constant
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("TEST B - paired cross-genome similarity")
    print("=" * 72)

    def jac(a, b):
        A, B = rep_of[a], rep_of[b]
        u = len(A | B)
        return len(A & B) / u if u else 0.0

    accs = sorted(pairs)
    sample = accs if len(accs) <= 120 else random.sample(accs, 120)
    same, cross, wins = [], [], 0
    n = 0
    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            a2, a3 = pairs[sample[i]]
            b2, b3 = pairs[sample[j]]
            s, c = jac(a3, b3), jac(a3, b2)
            same.append(s); cross.append(c)
            wins += (s > c)
            n += 1
    same, cross = np.array(same), np.array(cross)
    print(f"genome pairs compared: {n}")
    print(f"  mean J(c3_A, c3_B) = {same.mean():.4f}")
    print(f"  mean J(c3_A, c2_B) = {cross.mean():.4f}")
    print(f"  c3-c3 more similar than c3-c2 in {100*wins/n:.1f}% of pairs")
    from scipy.stats import wilcoxon
    try:
        st, p = wilcoxon(same, cross)
        print(f"  Wilcoxon signed-rank: statistic={st:.0f}, p={p:.3g}")
    except ValueError as e:
        print(f"  Wilcoxon failed: {e}")

    # ------------------------------------------------------------------
    # c3-diagnostic orthogroups and the c3 call
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("c3-DIAGNOSTIC ORTHOGROUPS AND c3 CALLS")
    print("=" * 72)
    diag = [o for o in common if on_c3[o] / on_any[o] >= 0.9]
    anti = [o for o in common if on_c3[o] / on_any[o] <= 0.1]
    print(f"c3-diagnostic  (>=90% of occurrences on c3): {len(diag)}")
    print(f"c2-diagnostic  (<=10% of occurrences on c3): {len(anti)}")
    didx = np.array([oidx[o] for o in diag], dtype=int)
    aidx = np.array([oidx[o] for o in anti], dtype=int)

    def content(r):
        """(c3-diagnostic content, c2-diagnostic content), each in [0, 1]."""
        d = M[ridx[r], didx].mean() if len(didx) else 0.0
        a = M[ridx[r], aidx].mean() if len(aidx) else 0.0
        return d, a

    def score(r):
        """+1 = looks like c3, -1 = looks like c2. Retained for the output table,
        but NOT used to call c3: a replicon carrying NEITHER signature scores ~0,
        and a midpoint boundary would silently call it c3. Positive evidence is
        required instead - see the three-class rule below."""
        d, a = content(r)
        return d - a

    ref_d3 = np.array([content(rl[1])[0] for rl in pairs.values()])
    ref_a3 = np.array([content(rl[1])[1] for rl in pairs.values()])
    ref_d2 = np.array([content(rl[0])[0] for rl in pairs.values()])
    ref_a2 = np.array([content(rl[0])[1] for rl in pairs.values()])
    print(f"\nreference smaller-secondary ('c3' by position): "
          f"c3-content mean={ref_d3.mean():.3f} sd={ref_d3.std():.3f} "
          f"min={ref_d3.min():.3f}")
    print(f"reference larger-secondary  ('c2' by position): "
          f"c3-content mean={ref_d2.mean():.3f} p95={np.percentile(ref_d2,95):.3f}  "
          f"c2-content mean={ref_a2.mean():.3f}")

    # The reference "c3" set is itself heterogeneous: high mean c3-content but a
    # near-zero minimum. Position alone does not identify pC3 -- in some lineages
    # the smaller secondary replicon is an unrelated megaplasmid. So classify by
    # POSITIVE evidence into three classes, with the c3 threshold set well above
    # the 95th percentile of chromosome-2 c3-content.
    D_CUT = max(0.10, float(np.percentile(ref_d2, 95)) * 2)
    A_CUT = 0.30
    print(f"\nthresholds: c3-content >= {D_CUT:.3f} -> c3;  "
          f"c2-content >= {A_CUT:.2f} -> chromosome 2;  otherwise other megaplasmid")
    frac_ref_c3_pass = float(np.mean(ref_d3 >= D_CUT))
    print(f"  of the {len(ref_d3)} position-defined c3 replicons, "
          f"{100*frac_ref_c3_pass:.0f}% carry enough c3-diagnostic content to qualify")
    print(f"  of the {len(ref_d2)} position-defined c2 replicons, "
          f"{100*np.mean(ref_d2 >= D_CUT):.0f}% would be misassigned (specificity check)")

    def classify(r):
        d, a = content(r)
        if d >= D_CUT and d > a:
            return "c3"
        if a >= A_CUT:
            return "chromosome2"
        return "other_megaplasmid"

    cls = {r: classify(r) for r in replicons}
    counts = collections.Counter(cls.values())
    print("\n=== REPLICON CLASSES ===")
    for k, v in counts.most_common():
        print(f"  {v:>5}  {k}")

    is_c3 = {r: cls[r] == "c3" for r in replicons}
    n_c3 = sum(is_c3.values())
    print(f"\nreplicons called c3: {n_c3}")
    print(f"genomes with >=1 c3: {len({r[0] for r in replicons if is_c3[r]})}")
    per = collections.Counter(r[0] for r in replicons if is_c3[r])
    print(f"c3 replicons per genome: {dict(collections.Counter(per.values()))}")
    c3_sizes = sorted(int(types[r]["length"]) for r in replicons if is_c3[r])
    if c3_sizes:
        print(f"called-c3 size: median={c3_sizes[len(c3_sizes)//2]:,} "
              f"min={c3_sizes[0]:,} max={c3_sizes[-1]:,}")
    om = [r for r in replicons if cls[r] == "other_megaplasmid"]
    if om:
        oms = sorted(int(types[r]["length"]) for r in om)
        print(f"\nother megaplasmids: n={len(om)} "
              f"median={oms[len(oms)//2]:,} max={oms[-1]:,}")
        print("  these carry neither a c3 nor a chromosome-2 gene signature and are")
        print("  EXCLUDED from the c3 pangenome; including them was inflating it with")
        print("  unrelated accessory content.")
    if MF6_C3 in rep_of:
        d, a = content(MF6_C3)
        print(f"\nMF6 cluster_003 ({int(types[MF6_C3]['length']):,} bp): "
              f"c3-content={d:.3f} c2-content={a:.3f} -> {cls[MF6_C3]}")
    else:
        print("\nMF6 not yet annotated - anchor unavailable")

    # ---- distances for the record ----
    D = pdist(M, metric="jaccard")
    Z = linkage(D, method="average")
    sq = squareform(D)
    lab6 = fcluster(Z, 6, criterion="maxclust")
    i_mf6 = ridx.get(MF6_C3)

    with open(f"{RES}/secondary_replicon_clusters.tsv", "w") as fh:
        fh.write("accession\tcontig\tlength\tsecondary_rank\tn_orthogroups\t"
                 "c3_content\tc2_content\tc3_score\treplicon_class\tcluster_k6\t"
                 "jaccard_to_MF6c3\tis_c3\torganism\n")
        for r in replicons:
            d = f"{sq[ridx[r], i_mf6]:.4f}" if i_mf6 is not None else "NA"
            fh.write(f"{r[0]}\t{r[1]}\t{types[r]['length']}\t{sec_rank[r]}\t"
                     f"{len(rep_of[r])}\t{content(r)[0]:.4f}\t{content(r)[1]:.4f}\t"
                     f"{score(r):.4f}\t{cls[r]}\t{lab6[ridx[r]]}\t{d}\t"
                     f"{is_c3[r]}\t{cen.get(r[0],{}).get('organism_name','')}\n")
    with open(f"{RES}/c3_diagnostic_orthogroups.txt", "w") as fh:
        fh.write("\n".join(diag) + "\n")
    with open(f"{RES}/orthogroup_residence.tsv", "w") as fh:
        fh.write("orthogroup\tn_genomes\tfrac_on_c3\n")
        for o in common:
            fh.write(f"{o}\t{on_any[o]}\t{on_c3[o]/on_any[o]:.4f}\n")
    print(f"\nwrote {RES}/secondary_replicon_clusters.tsv")
    print("STAGE4C_DONE")


if __name__ == "__main__":
    sys.exit(main())
