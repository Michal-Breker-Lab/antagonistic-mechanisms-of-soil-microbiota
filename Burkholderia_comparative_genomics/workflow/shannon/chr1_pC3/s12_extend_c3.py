#!/usr/bin/env python3
"""Extend gene-content c3 calls from the 306 annotated genomes to all 771.

Why this is defensible: dereplication clustered at 99% ANI, which is
strain-level identity, and any genome whose replicon ARCHITECTURE differed from
its cluster representative was force-retained and annotated in its own right. So
a non-representative genome shares both ~99% ANI and the same architecture as its
representative, and inheriting the c3 call adds little risk.

Where it is NOT safe, the genome is marked `inferred_uncertain` rather than given
a call: singleton clusters with no annotated member, or architecture mismatches
that somehow survived. Both the direct and extended counts are reported so the
reader can use whichever they trust.
"""
import collections
import csv
import os
import sys

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES = f"{W}/results"
ANI_CUT, AF_CUT = 99.0, 50.0


def main():
    cen = {r["accession"]: r for r in
           csv.DictReader(open(f"{RES}/replicon_census.tsv"), delimiter="\t")}
    direct = {}
    for r in csv.DictReader(open(f"{RES}/secondary_replicon_clusters.tsv"),
                            delimiter="\t"):
        acc = r["accession"]
        if r.get("replicon_class") == "c3":
            direct[acc] = True
        else:
            direct.setdefault(acc, False)
    # genomes annotated but with NO large secondary replicon at all -> no c3
    annotated = set()
    for r in csv.DictReader(open(f"{RES}/replicon_types.tsv"), delimiter="\t"):
        annotated.add(r["accession"])
    for a in annotated:
        direct.setdefault(a, False)
    direct.pop("MF6", None)
    print(f"genomes with a direct gene-content c3 call: {len(direct)}")
    print(f"  of these, c3 present: {sum(direct.values())}")

    # rebuild the 99% ANI clusters over all 771
    accs = sorted(cen)
    idx = {a: i for i, a in enumerate(accs)}
    parent = list(range(len(accs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    with open(f"{RES}/skani_triangle_edges.tsv") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                ani = float(r["ANI"])
                af = min(float(r["Align_fraction_ref"]),
                         float(r["Align_fraction_query"]))
            except (KeyError, ValueError):
                continue
            if ani < ANI_CUT or af < AF_CUT:
                continue
            a = r["Ref_file"].split("/")[-1].replace(".fna", "")
            b = r["Query_file"].split("/")[-1].replace(".fna", "")
            if a in idx and b in idx:
                ra, rb = find(idx[a]), find(idx[b])
                if ra != rb:
                    parent[rb] = ra

    clusters = collections.defaultdict(list)
    for a in accs:
        clusters[find(idx[a])].append(a)

    out, stats = {}, collections.Counter()
    for members in clusters.values():
        called = [m for m in members if m in direct]
        for m in members:
            if m in direct:
                out[m] = (direct[m], "direct")
                stats["direct"] += 1
                continue
            # inherit only from a cluster-mate with the SAME architecture
            same = [c for c in called
                    if cen[c]["architecture"] == cen[m]["architecture"]]
            if same:
                v = any(direct[c] for c in same)
                out[m] = (v, "inherited")
                stats["inherited"] += 1
            elif called:
                out[m] = (None, "uncertain_arch_mismatch")
                stats["uncertain_arch_mismatch"] += 1
            else:
                out[m] = (None, "uncertain_no_annotated_member")
                stats["uncertain_no_annotated_member"] += 1

    print("\n=== extension of c3 calls to all genomes ===")
    for k, v in stats.most_common():
        print(f"  {v:>4}  {k}")
    n_yes = sum(1 for v, _ in out.values() if v is True)
    n_no = sum(1 for v, _ in out.values() if v is False)
    n_unk = sum(1 for v, _ in out.values() if v is None)
    tot = len(out)
    print(f"\nc3 PRESENT: {n_yes}/{tot} ({100*n_yes/tot:.1f}%)")
    print(f"c3 ABSENT : {n_no}/{tot} ({100*n_no/tot:.1f}%)")
    print(f"uncertain : {n_unk}/{tot} ({100*n_unk/tot:.1f}%)")

    with open(f"{RES}/c3_calls_all_genomes.tsv", "w") as fh:
        fh.write("accession\torganism_name\tarchitecture\tc3_present\tevidence\t"
                 "qc_pass\n")
        for a in sorted(out):
            v, ev = out[a]
            fh.write(f"{a}\t{cen[a]['organism_name']}\t{cen[a]['architecture']}\t"
                     f"{'' if v is None else v}\t{ev}\t{cen[a].get('qc_pass','')}\n")
    print(f"\nwrote {RES}/c3_calls_all_genomes.tsv")

    # architecture vs c3 identity -- the key comparison
    print("\n=== ARCHITECTURE vs c3 IDENTITY ===")
    print("(does 'has 3 large replicons' mean 'has pC3'?)")
    ct = collections.defaultdict(collections.Counter)
    for a, (v, ev) in out.items():
        ct[cen[a]["architecture"]][v] += 1
    print(f"{'architecture':<12} {'c3 present':>11} {'c3 absent':>10} {'unknown':>8}")
    for arch in ["1_large", "2_large", "3_large", "4+_large", "0_large"]:
        if arch not in ct:
            continue
        c = ct[arch]
        print(f"{arch:<12} {c[True]:>11} {c[False]:>10} {c[None]:>8}")

    # per-species table for the genomes with a confident call
    sp = collections.defaultdict(collections.Counter)
    for a, (v, ev) in out.items():
        if v is not None:
            sp[cen[a]["organism_name"]][v] += 1
    var = [(n, c) for n, c in sp.items() if c[True] and c[False]]
    var.sort(key=lambda x: -sum(x[1].values()))
    print(f"\n=== SPECIES WHERE c3 IS PRESENT IN SOME STRAINS AND ABSENT IN OTHERS ===")
    print(f"{'species':<42} {'with c3':>8} {'without':>8}")
    for n, c in var[:25]:
        print(f"{n[:42]:<42} {c[True]:>8} {c[False]:>8}")
    print(f"\n{len(var)} species show within-species variation in c3 presence")
    print("STAGE12_DONE")


if __name__ == "__main__":
    sys.exit(main())
