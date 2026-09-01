#!/usr/bin/env python3
"""Stage 5 - dereplicate at 99% ANI (single linkage) for the pangenome and trees.
The Stage 3 census deliberately keeps all 771 genomes; only downstream stages
use this reduced set.

Refinement over plain dereplication: within a cluster, genomes whose replicon
ARCHITECTURE differs from the chosen representative are force-retained. Without
this, a pair of 99.5%-identical strains where one has lost c3 would collapse to
one genome and the most interesting signal in the dataset would be discarded by
a QC step.
"""
import collections
import csv
import os
import subprocess
import sys

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES = f"{W}/results"
ANI_CUT = 99.0
AF_CUT = 50.0          # require decent alignment before calling two genomes redundant
TRI = f"{RES}/skani_triangle_edges.tsv"


def run_triangle():
    if os.path.exists(TRI) and os.path.getsize(TRI) > 0:
        print(f"reusing {TRI}")
        return
    print("running skani triangle (all-vs-all, 771 genomes)...", flush=True)
    subprocess.run(
        [f"{W}/envs/burk/bin/skani", "triangle", "-l", f"{W}/genome_list.txt",
         "-t", "48", "-E", "-o", TRI, "--min-af", "30"],
        check=True)


def main():
    run_triangle()
    cen = {r["accession"]: r for r in
           csv.DictReader(open(f"{RES}/replicon_census.tsv"), delimiter="\t")}
    accs = list(cen)
    idx = {a: i for i, a in enumerate(accs)}

    # union-find over edges above the ANI cut
    parent = list(range(len(accs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    n_edge = 0
    with open(TRI) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                ani = float(r["ANI"])
                af = min(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"]))
            except (KeyError, ValueError):
                continue
            if ani < ANI_CUT or af < AF_CUT:
                continue
            a = r["Ref_file"].split("/")[-1].replace(".fna", "")
            b = r["Query_file"].split("/")[-1].replace(".fna", "")
            if a in idx and b in idx:
                union(idx[a], idx[b])
                n_edge += 1
    print(f"edges >= {ANI_CUT}% ANI and >= {AF_CUT}% AF: {n_edge}")

    clusters = collections.defaultdict(list)
    for a in accs:
        clusters[find(idx[a])].append(a)
    print(f"clusters: {len(clusters)}  (from {len(accs)} genomes)")

    def quality(a):
        """higher is better: completeness, low contamination, long-read, N50"""
        r = cen[a]
        try:
            comp = float(r["checkm_completeness"] or 0)
        except ValueError:
            comp = 0.0
        try:
            cont = float(r["checkm_contamination"] or 100)
        except ValueError:
            cont = 100.0
        tech = (r["sequencing_tech"] or "").lower()
        longread = any(k in tech for k in ("pacbio", "nanopore", "minion", "hifi", "pabio"))
        ref = 1 if r.get("ncbi_n_chromosomes") else 0
        try:
            n50 = int(cen[a].get("largest") or 0)
        except ValueError:
            n50 = 0
        return (comp, -cont, longread, ref, n50)

    keep, forced = [], []
    for members in clusters.values():
        members.sort(key=quality, reverse=True)
        rep = members[0]
        keep.append(rep)
        rep_arch = cen[rep]["architecture"]
        for m in members[1:]:
            if cen[m]["architecture"] != rep_arch:
                keep.append(m)
                forced.append((m, rep, cen[m]["architecture"], rep_arch))

    keep = sorted(set(keep))
    print(f"representatives: {len(clusters)}")
    print(f"force-retained for differing architecture: {len(forced)}")
    print(f"FINAL DEREPLICATED SET: {len(keep)}")

    with open(f"{RES}/derep_representatives.tsv", "w") as fh:
        fh.write("accession\torganism_name\tarchitecture\tcluster_size\treason\n")
        cl_of = {}
        for cid, members in clusters.items():
            for m in members:
                cl_of[m] = (cid, len(members))
        forced_set = {f[0] for f in forced}
        for a in keep:
            cid, sz = cl_of[a]
            fh.write(f"{a}\t{cen[a]['organism_name']}\t{cen[a]['architecture']}\t{sz}\t"
                     f"{'architecture_variant' if a in forced_set else 'representative'}\n")
    with open(f"{RES}/derep_accessions.txt", "w") as fh:
        fh.write("\n".join(keep) + "\n")

    if forced:
        print("\nexamples of architecture variants rescued from dereplication:")
        for m, rep, am, ar in forced[:12]:
            print(f"  {m} ({am}) kept alongside {rep} ({ar})  [{cen[m]['organism_name']}]")

    big = sorted(clusters.values(), key=len, reverse=True)[:10]
    print("\nlargest redundancy clusters collapsed:")
    for members in big:
        names = collections.Counter(cen[m]["organism_name"] for m in members)
        top = names.most_common(1)[0]
        print(f"  n={len(members):<4} {top[0]}")
    print("STAGE5_DONE")


if __name__ == "__main__":
    sys.exit(main())
