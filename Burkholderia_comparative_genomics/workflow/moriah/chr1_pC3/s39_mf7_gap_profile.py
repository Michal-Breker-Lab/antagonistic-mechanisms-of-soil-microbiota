#!/usr/bin/env python3
"""Where, along MF6's pC3, do the genes missing from MF7's assembly sit?

MF6 and MF7 are the same strain (ANI 100.00) isolated and assembled separately.
MF6's pC3 is a closed 1.18 Mb Nanopore replicon; MF7's is a single 858 kb draft
contig. The gene-count difference (1012 vs 759) therefore has to be assembly, not
biology -- but WHERE the loss sits decides whether a strict 6/6 core would merely
shrink or would preferentially delete one functional block. A contiguous loss
biases the COG profile; a scattered one only costs power.

Usage: s39_mf7_gap_profile.py <setb_pc3_dir> <blast_bin_dir>
"""
import subprocess
import sys
import tempfile
from pathlib import Path

d = Path(sys.argv[1])
bindir = Path(sys.argv[2])
ID_MIN, COV_MIN = 90.0, 80.0

with tempfile.TemporaryDirectory(dir=str(d)) as tmp:
    tmp = Path(tmp)
    subprocess.run([str(bindir / "makeblastdb"), "-in", str(d / "MF7_pC3.faa"),
                    "-dbtype", "prot", "-out", str(tmp / "mf7")],
                   check=True, stdout=subprocess.DEVNULL)
    out = subprocess.run(
        [str(bindir / "blastp"), "-query", str(d / "MF6_pC3.faa"),
         "-db", str(tmp / "mf7"), "-evalue", "1e-5", "-max_target_seqs", "1",
         "-num_threads", "4", "-outfmt", "6 qseqid pident qcovhsp"],
        check=True, capture_output=True, text=True).stdout

best = {}
for line in out.splitlines():
    q, pid, cov = line.split("\t")
    if q not in best:                      # blast emits best HSP first
        best[q] = (float(pid), float(cov))
found = {q for q, (p, c) in best.items() if p >= ID_MIN and c >= COV_MIN}

# CDS start coordinates along the MF6 pC3, in replicon order
pos = []
for line in open(d / "MF6_pC3.gff3"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS":
        continue
    lt = next((a.split("=", 1)[1] for a in f[8].split(";")
               if a.startswith("locus_tag=")), None)
    if lt:
        pos.append((int(f[3]), lt))
pos.sort()

n = len(pos)
print(f"MF6 pC3 CDS: {n}   with a confident MF7 ortholog: {len(found & {l for _, l in pos})}")
print("\ndecile  start_bp     n_CDS  n_missing  pct_missing")
for dec in range(10):
    lo, hi = dec * n // 10, (dec + 1) * n // 10
    chunk = pos[lo:hi]
    miss = sum(1 for _, lt in chunk if lt not in found)
    print(f"{dec+1:>4}  {chunk[0][0]:>10,}  {len(chunk):>7}  {miss:>8}  {100*miss/len(chunk):>9.1f}%")

# longest run of consecutive missing CDS -- the signature of a truncation
run = best_run = 0
for _, lt in pos:
    run = run + 1 if lt not in found else 0
    best_run = max(best_run, run)
print(f"\nlongest consecutive run of MF6 pC3 CDS absent from MF7: {best_run}")

# ---- does the truncation delete a functionally biased slice? ---------------
# A contiguous loss is only a problem for the COG profile if the lost window
# differs in content from the rest. Compare COG-category composition of the
# retained vs deleted MF6 pC3 genes, using the Bakta DbXrefs already on disk.
from collections import Counter                                    # noqa: E402

ANNOT = Path(sys.argv[3]) if len(sys.argv) > 3 else None
if ANNOT:
    cats = {}
    for line in open(ANNOT):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[1] != "cds":
            continue
        # Bakta writes one "COG:<letters>" xref alongside "COG:COG#####"
        cc = [x.split(":", 1)[1] for x in f[8].split(", ")
              if x.startswith("COG:") and not x.startswith("COG:COG")]
        cats[f[5]] = cc[0] if cc else "-"

    keep = Counter(); lost = Counter()
    for _, lt in pos:
        for ch in (cats.get(lt, "-") or "-"):
            (keep if lt in found else lost)[ch] += 1
    nk = sum(1 for _, lt in pos if lt in found)
    nl = len(pos) - nk
    print(f"\nCOG composition, retained (n={nk}) vs deleted-by-truncation (n={nl}):")
    print("cat   retained%   deleted%   diff")
    for ch in sorted(set(keep) | set(lost)):
        pk, pl = 100 * keep[ch] / nk, 100 * lost[ch] / nl
        if max(pk, pl) >= 1.0:
            print(f"  {ch}  {pk:>8.1f}  {pl:>9.1f}  {pl-pk:>+6.1f}")
