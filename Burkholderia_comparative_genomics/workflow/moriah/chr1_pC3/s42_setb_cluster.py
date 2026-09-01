#!/usr/bin/env python3
"""Cluster the Set B replicon proteins into gene families, and recover the
identity/coverage the original (lost) generator used.

The retained tables/setB_functional_families.tsv is a THIRD clustering of these
genomes -- neither the Panaroo nor the PPanGGOLiN run reported in report section 6.3
(which give 1,346 families / 843 core). Its family ids are Bakta locus tags, i.e.
MMseqs2-style representatives. Its settings were never written down, so --sweep
re-runs easy-cluster across a grid and reports which cell reproduces the retained
counts on the ORIGINAL five closed genomes:

    pC3   1,382 families   833 core (5/5)   549 accessory (1-4/5)
    chr1                 2,847 core (5/5)

Once a cell reproduces those, --emit rebuilds the family table for real, over
whichever genomes are asked for.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--faa", required=True, type=Path)
ap.add_argument("--mmseqs", required=True, type=Path)
ap.add_argument("--tmp", required=True, type=Path)
ap.add_argument("--exclude", default="", help="comma-separated accessions to drop")
ap.add_argument("--sweep", action="store_true")
ap.add_argument("--min-seq-id", type=float, default=0.5)
ap.add_argument("--cov", type=float, default=0.8)
ap.add_argument("--cov-mode", type=int, default=0)
ap.add_argument("--threads", type=int, default=0,
                help="mmseqs --threads; 0 = leave to mmseqs. Cluster boundaries can\n                     shift slightly with thread count in tie cases, so pin it to\n                     make a run reproducible.")
ap.add_argument("--out", type=Path)
a = ap.parse_args()

EXCL = {x for x in a.exclude.split(",") if x}


def cluster(faa, ident, cov, covmode, workdir):
    pref = workdir / f"c_{ident}_{cov}_{covmode}"
    cmd = [str(a.mmseqs), "easy-cluster", str(faa), str(pref), str(workdir / "mmtmp"),
           "--min-seq-id", str(ident), "-c", str(cov), "--cov-mode", str(covmode),
           "-v", "0"]
    if a.threads:
        cmd += ["--threads", str(a.threads)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    fams = defaultdict(set)          # representative -> set of accessions
    members = defaultdict(list)      # representative -> [(acc, locus_tag)]
    for line in open(f"{pref}_cluster.tsv"):
        rep, mem = line.rstrip("\n").split("\t")
        acc, tag = mem.split("|", 1)
        if acc in EXCL:
            continue
        fams[rep].add(acc)
        members[rep].append((acc, tag))
    fams = {r: g for r, g in fams.items() if g}       # reps whose members were all excluded
    return fams, {r: members[r] for r in fams}


a.tmp.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(dir=str(a.tmp)) as wd:
    wd = Path(wd)
    if a.sweep:
        genomes = set()
        for line in open(a.faa):
            if line.startswith(">"):
                genomes.add(line[1:].split("|", 1)[0])
        n = len(genomes - EXCL)
        print(f"sweeping over {n} genomes (excluded: {sorted(EXCL) or 'none'})")
        print(f"{'min-seq-id':>11}{'cov':>6}{'cov-mode':>10}{'families':>10}{'core':>8}{'accessory':>11}")
        for covmode in (0, 1):
            for cov in (0.5, 0.8):
                for ident in (0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98):
                    fams, _ = cluster(a.faa, ident, cov, covmode, wd)
                    core = sum(1 for g in fams.values() if len(g) == n)
                    print(f"{ident:>11}{cov:>6}{covmode:>10}{len(fams):>10,}"
                          f"{core:>8,}{len(fams)-core:>11,}")
        sys.exit(0)

    fams, members = cluster(a.faa, a.min_seq_id, a.cov, a.cov_mode, wd)
    genomes = sorted({acc for g in fams.values() for acc in g})
    n = len(genomes)
    if not a.out:
        sys.exit("FAIL: --out required without --sweep")
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "n_genomes", "n_members", "is_core", "is_single_copy_core", "genomes", "members"])
        for rep in sorted(fams):
            g = sorted(fams[rep])
            mm = sorted(members[rep])
            w.writerow([rep, len(g), len(mm), int(len(g) == n),
                        int(len(g) == n and len(mm) == n), ";".join(g),
                        ";".join(f"{acc}|{tag}" for acc, tag in mm)])
    core = sum(1 for g in fams.values() if len(g) == n)
    # "single-copy core" = one member per genome and no paralogs. The retained
    # setB_functional_families.tsv carries n_members == 5 for every core family,
    # i.e. exactly one gene per genome, so this -- not the 843 all-genomes count
    # -- is what the lost generator called the clade core.
    sc = sum(1 for rep, g in fams.items()
             if len(g) == n and len(members[rep]) == n)
    print(f"genomes {n}: {', '.join(genomes)}")
    print(f"families {len(fams):,}   core ({n}/{n}) {core:,}   "
          f"single-copy core {sc:,}   accessory {len(fams)-core:,}")
    print(f"wrote {a.out}")
