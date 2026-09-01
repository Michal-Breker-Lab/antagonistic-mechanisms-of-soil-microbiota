#!/usr/bin/env python3
"""Core of MF6's closest relatives, as a SUBSET of the full pangenomes.

Complements s16_neighbour_pangenomes.sh. That script re-clusters on the 5/11
genomes alone, which is the direct answer to "what is their core". This one
subsets the 140/306-genome matrices instead, so the family definitions are
identical to the ones used everywhere else in the report and the two numbers can
be read side by side.

Caveat that must be stated wherever these numbers appear: families here were
clustered across the FULL set, so a family split into two by a distant genome's
paralogue stays split even when the 11 genomes carry only one copy. That biases
this version's core DOWN relative to the re-clustered version.

At n=11 any threshold from 0.91 to 1.00 selects the same families
(ceil(0.95*11) == ceil(1.00*11) == 11), so core_ge_95pct == core_ge_100pct is an
arithmetic identity here, not a coincidence -- it is used as a self-check.
"""
import csv
import sys
from pathlib import Path

W = Path("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")

SET_A = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1",
         "GCF_003966315.1", "GCF_001718515.1", "GCF_014211915.1", "GCF_003076415.1",
         "GCF_000019505.1", "GCF_000203955.1", "MF6"]
SET_B = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1",
         "GCF_053209605.1", "MF6"]
THRESHOLDS = [1.00, 0.95, 0.90, 0.80, 0.60]


def core_of(rtab, members):
    with open(rtab) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = {g.strip('"'): i for i, g in enumerate(hdr[1:], start=1)}
        present = [g for g in members if g in cols]
        idx = [cols[g] for g in present]
        n = len(present)
        if n == 0:
            return 0, 0, {}, []
        counts = []
        for line in fh:
            p = line.rstrip("\n").split("\t")
            counts.append(sum(1 for i in idx if i < len(p) and p[i] == "1"))
    any_present = sum(1 for c in counts if c > 0)
    missing = [g for g in members if g not in cols]
    return n, any_present, {t: sum(1 for c in counts if c >= t * n) for t in THRESHOLDS}, missing


targets = []
for lbl, p in (("panaroo_c3_moderate", W / "pangenome" / "c3_moderate" / "gene_presence_absence.Rtab"),
               ("panaroo_c3_sensitive", W / "pangenome" / "c3_sensitive" / "gene_presence_absence.Rtab"),
               ("panaroo_chr1_strict", W / "pangenome" / "chr1_strict" / "gene_presence_absence.Rtab")):
    if p.exists():
        targets.append((lbl, p))
for d in sorted((W / "ppanggolin").glob("*_id*")):
    r = d / "out" / "gene_presence_absence.Rtab"
    if r.exists():
        targets.append((f"ppanggolin_{d.name}", r))

rows = []
for label, rtab in targets:
    for sname, members in (("setA_10closest", SET_A), ("setB_sola_clade", SET_B)):
        n, anyp, cores, missing = core_of(rtab, members)
        row = {"pangenome": label, "genome_set": sname,
               "n_requested": len(members), "n_found": n,
               "missing": ",".join(missing) if missing else "",
               "n_families_in_any": anyp}
        for t in THRESHOLDS:
            row[f"core_ge_{int(t*100)}pct"] = cores.get(t, "")
        rows.append(row)
        print(row)

if not rows:
    sys.exit("no matrices found")

out = W / "results" / "mf6_neighbour_core_subset.tsv"
with open(out, "w", newline="") as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    wtr.writeheader()
    wtr.writerows(rows)
print(f"\nwrote {out}")

# self-check: at n<=20, ceil(0.95n) == n, so the 95% and 100% cores must agree
for r in rows:
    if r["n_found"] and r["n_found"] <= 20:
        if r["core_ge_95pct"] != r["core_ge_100pct"]:
            print(f"WARN {r['pangenome']}/{r['genome_set']}: 95% and 100% cores "
                  f"differ ({r['core_ge_95pct']} vs {r['core_ge_100pct']}) at "
                  f"n={r['n_found']} -- counting is wrong", file=sys.stderr)
