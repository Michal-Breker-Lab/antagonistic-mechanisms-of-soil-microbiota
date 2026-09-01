#!/usr/bin/env python3
"""Derive pC3-diagnostic gene families and classify every large secondary replicon.

Stage 4 of the original pipeline, per report_methods_draft.md:

    "Proteins from all large secondary replicons were clustered into orthologous
     groups with MMseqs2 v18.8cc5c (easy-cluster, 50% identity, 80% coverage)."

MMseqs2 clusters ALL secondary replicons in one pass, so every replicon in the
set can be classified -- there is no training/application reach limit. Decision
D9a's clone-collapsing therefore applies only to WHICH genomes' positional
labels contribute to the residence statistic, not to which genomes are clustered.

Steps
  1. Residence -- for each family, the fraction of its occurrences on the
     positional pC3 among CLONE-COLLAPSED TRAINING genomes (D9a). Families in
     fewer than --min-genomes training genomes are dropped.
  2. Diagnostics -- frac >= --frac is pC3-diagnostic, <= 1-(--frac) is
     chromosome-2-diagnostic. The 0.90/10 combination reproduces the retained
     `c3_diagnostic_orthogroups.txt` exactly (1,013 families, 0 missing/extra).
  3. Classification -- score EVERY large secondary replicon by its diagnostic
     content. Replicons with no family membership are reported `unclassifiable`,
     NOT silently called chromosome 2.
  4. D9b -- self-consistency against the training genomes' own positional labels.

Sequence ids must be `accession|contig|locus_tag` (see s31_secondary_proteins.py).
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def parse_id(gid: str):
    p = gid.split("|")
    return (p[0], p[1]) if len(p) >= 3 else None


def read_clusters(path: Path) -> dict:
    """MMseqs2 easy-cluster TSV: representative<TAB>member. -> {rep: [members]}"""
    fam = collections.defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            fam[p[0]].append(p[1])
    return fam


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--labels", required=True, type=Path)
    ap.add_argument("--types", required=True, type=Path)
    ap.add_argument("--frac", type=float, default=0.90)
    ap.add_argument("--min-genomes", type=int, default=10)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--tag", default="")
    ap.add_argument("--c2-chr", type=float, default=0.30,
                    help="c2_content at/above which a replicon is chromosome 2")
    ap.add_argument("--c3-min", type=float, default=0.10,
                    help="minimum c3_content for a pC3 call")
    ap.add_argument("--other-pct", type=float, default=20.0,
                    help="percentile of the reference signature distributions "
                         "below which a replicon carries neither signature "
                         "(other_megaplasmid). The single fitted parameter: 20 "
                         "maximises agreement with the retained classification "
                         "(565/586 = 96.4% on shared replicons).")
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)
    tag = f"_{a.tag}" if a.tag else ""

    role = {}
    for r in csv.DictReader(open(a.labels), delimiter="\t"):
        role[(r["accession"], r["contig"])] = r["role"]
    training = {acc for acc, _ in role}
    print(f"training genomes {len(training)}   labelled replicons {len(role)}")

    fam = read_clusters(a.clusters)
    print(f"MMseqs2 families {len(fam):,d}   members "
          f"{sum(len(v) for v in fam.values()):,d}")

    # ---- 1. residence, on training genomes only -----------------------------
    residence, fam_keep = [], {}
    for rep, members in fam.items():
        on_c3 = on_c2 = 0
        genomes = set()
        for gid in members:
            loc = parse_id(gid)
            if not loc:
                continue
            r = role.get(loc)
            if r == "pC3":
                on_c3 += 1; genomes.add(loc[0])
            elif r == "chromosome2":
                on_c2 += 1; genomes.add(loc[0])
        tot = on_c3 + on_c2
        if tot == 0 or len(genomes) < a.min_genomes:
            continue
        residence.append((rep, len(genomes), on_c3 / tot))
        fam_keep[rep] = members

    with open(a.outdir / f"orthogroup_residence{tag}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["orthogroup", "n_genomes", "frac_on_c3"])
        for rep, n, fr in residence:
            w.writerow([rep, n, f"{fr:.4f}"])
    print(f"families with n_genomes>={a.min_genomes}: {len(residence)}")

    # ---- 2. diagnostics -----------------------------------------------------
    c3_diag = {r for r, _n, fr in residence if fr >= a.frac}
    c2_diag = {r for r, _n, fr in residence if fr <= 1 - a.frac}
    print(f"pC3-diagnostic families        : {len(c3_diag)}")
    print(f"chromosome2-diagnostic families: {len(c2_diag)}")
    (a.outdir / f"c3_diagnostic_orthogroups{tag}.txt").write_text(
        "\n".join(sorted(c3_diag)) + "\n")

    # ---- 3. classify EVERY large secondary replicon -------------------------
    # Content is the fraction of the DIAGNOSTIC SET a replicon carries, counted
    # over DISTINCT families -- not the fraction of the replicon's own families
    # that are diagnostic. The two differ by an order of magnitude because the
    # c2-diagnostic set is ~5x the c3 set; using the replicon-denominator makes
    # c2_content ~10x too large and collapses "other_megaplasmid" into
    # chromosome2 (208 of 289 in testing). Validated by reproducing the retained
    # table's class proportions and n_orthogroups (1,299 vs 1,350 median).
    per_rep = collections.defaultdict(lambda: [set(), set(), set()])
    for rep, members in fam_keep.items():
        is3, is2 = rep in c3_diag, rep in c2_diag
        for gid in members:
            loc = parse_id(gid)
            if not loc:
                continue
            cell = per_rep[loc]
            cell[2].add(rep)
            if is3:
                cell[0].add(rep)
            elif is2:
                cell[1].add(rep)
    n_c3d, n_c2d = max(len(c3_diag), 1), max(len(c2_diag), 1)

    sec = [r for r in csv.DictReader(open(a.types), delimiter="\t")
           if r["replicon_type"] == "secondary_large"]

    # secondary_rank: 1 = largest secondary replicon in that genome. Used only
    # as the documented tie-breaker within +/-0.05 of the boundary.
    bylen = collections.defaultdict(list)
    for r in sec:
        bylen[r["accession"]].append((int(r["length"]), r["contig"]))
    srank = {}
    for acc, v in bylen.items():
        for i, (_l, c) in enumerate(sorted(v, reverse=True), start=1):
            srank[(acc, c)] = i

    content = {}
    for r in sec:
        key = (r["accession"], r["contig"])
        s3, s2, sall = per_rep.get(key, (set(), set(), set()))
        content[key] = (len(s3) / n_c3d, len(s2) / n_c2d, len(sall))

    # ---- Option 3 calibration ---------------------------------------------
    # (a) Decision boundary between c3 and chromosome2 is the MIDPOINT of the
    #     two reference means, exactly as report_methods_draft.md specifies.
    #     Being a midpoint of the observed references, it is self-calibrating
    #     and unaffected by the content scale differing from the original run.
    # (b) The "carries neither signature" (other_megaplasmid) carve-out is the
    #     one fitted parameter: a replicon is other_megaplasmid when it falls
    #     below the --other-pct percentile of BOTH reference sets on their own
    #     signature. Declared, not inherited.
    ref3 = [content[k][0] for k, v in role.items() if v == "pC3" and k in content]
    ref2 = [content[k][1] for k, v in role.items()
            if v == "chromosome2" and k in content]
    sc3 = [content[k][0] - content[k][1] for k, v in role.items()
           if v == "pC3" and k in content]
    sc2 = [content[k][0] - content[k][1] for k, v in role.items()
           if v == "chromosome2" and k in content]
    import statistics as _st
    boundary = (_st.mean(sc3) + _st.mean(sc2)) / 2.0

    def pct(xs, q):
        xs = sorted(xs)
        if not xs:
            return 0.0
        i = max(0, min(len(xs) - 1, int(round(q / 100.0 * (len(xs) - 1)))))
        return xs[i]

    floor3 = pct(ref3, a.other_pct)
    floor2 = pct(ref2, a.other_pct)
    print(f"\n=== Option 3 calibration ===")
    print(f"reference pC3  n={len(sc3)} mean score {_st.mean(sc3):+.4f}")
    print(f"reference chr2 n={len(sc2)} mean score {_st.mean(sc2):+.4f}")
    print(f"decision boundary (midpoint) : {boundary:+.4f}")
    print(f"other_megaplasmid floors at P{a.other_pct}: "
          f"c3<{floor3:.4f} and c2<{floor2:.4f}")
    ok3 = sum(1 for x in sc3 if x > boundary)
    ok2 = sum(1 for x in sc2 if x <= boundary)
    print(f"boundary accuracy on references: pC3 {ok3}/{len(sc3)} "
          f"({100*ok3/len(sc3):.1f}%)  chr2 {ok2}/{len(sc2)} "
          f"({100*ok2/len(sc2):.1f}%)")

    rows = []
    for r in sec:
        key = (r["accession"], r["contig"])
        c3f, c2f, tot = content[key]
        score = c3f - c2f
        if tot == 0:
            cls, is_c3 = "unclassifiable", ""
        elif c3f < floor3 and c2f < floor2:
            cls, is_c3 = "other_megaplasmid", "False"
        else:
            if abs(score - boundary) <= 0.05:
                # documented tie-breaker: the larger secondary replicon is
                # chromosome 2, the smaller is pC3
                cls = "c3" if srank[key] > 1 else "chromosome2"
            else:
                cls = "c3" if score > boundary else "chromosome2"
            is_c3 = "True" if cls == "c3" else "False"
        rows.append({"accession": r["accession"], "contig": r["contig"],
                     "length": r["length"], "secondary_rank": srank[key],
                     "n_orthogroups": tot,
                     "c3_content": f"{c3f:.4f}", "c2_content": f"{c2f:.4f}",
                     "c3_score": f"{score:.4f}",
                     "replicon_class": cls, "is_c3": is_c3})
    with open(a.outdir / f"secondary_replicon_clusters{tag}.tsv", "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    cnt = collections.Counter(r["replicon_class"] for r in rows)
    print(f"replicons scored {len(rows)}: {dict(cnt)}")
    pos = {r["accession"] for r in rows if r["replicon_class"] == "c3"}
    print(f"genomes pC3-positive: {len(pos)}")

    # ---- genome-level calls -------------------------------------------------
    census_acc = sorted({r["accession"] for r in
                         csv.DictReader(open(a.types), delimiter="\t")})
    with open(a.outdir / f"c3_calls_all_genomes{tag}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "c3_present", "evidence", "n_secondary_large"])
        bysec = collections.Counter(r["accession"] for r in rows)
        for acc in census_acc:
            w.writerow([acc, str(acc in pos), "direct", bysec.get(acc, 0)])

    # ---- 4. D9b self-consistency -------------------------------------------
    idx = {(r["accession"], r["contig"]): r for r in rows}
    inf = [(k, v) for k, v in role.items()
           if k in idx and int(idx[k]["n_orthogroups"]) > 0]
    # A HARD inversion is the classifier swapping the two labels -- the failure
    # D9b exists to catch. Declining to call (-> other_megaplasmid) is a
    # legitimate outcome, not a contradiction: the original itself finds 293
    # other_megaplasmids, many of them the smaller replicon of a two-replicon
    # genome. The two are reported separately.
    def inverted(call, positional):
        return ((positional == "pC3" and call == "chromosome2")
                or (positional == "chromosome2" and call == "c3"))
    dis = [(k, v) for k, v in inf if inverted(idx[k]["replicon_class"], v)]
    declined = [(k, v) for k, v in inf
                if idx[k]["replicon_class"] == "other_megaplasmid"]
    print(f"declined to call (-> other_megaplasmid): {len(declined)} / {len(inf)}"
          f"  ({100*len(declined)/len(inf) if inf else 0:.2f}%)")
    print(f"\n=== D9b self-consistency ===")
    print(f"training replicons informative: {len(inf)} / {len(role)}")
    print(f"HARD inversions (pC3<->chromosome2): {len(dis)}"
          f"  ({100*len(dis)/len(inf) if inf else 0:.2f}%)")
    for k, v in dis[:12]:
        r = idx[k]
        print(f"   {k[0]} {k[1]} positional={v} called={r['replicon_class']} "
              f"c3={r['c3_content']} c2={r['c2_content']} n={r['n_orthogroups']}")
    with open(a.outdir / f"pc3_self_consistency{tag}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "contig", "positional_role", "called_class",
                    "c3_content", "c2_content", "n_orthogroups"])
        for k, v in dis:
            r = idx[k]
            w.writerow([k[0], k[1], v, r["replicon_class"], r["c3_content"],
                        r["c2_content"], r["n_orthogroups"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
