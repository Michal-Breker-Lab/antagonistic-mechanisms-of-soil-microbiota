#!/usr/bin/env python3
"""Paralogue-collapse diagnostic for the PPanGGOLiN identity arms.

At 60% identity with single-linkage clustering, paralogues can fuse into their
orthologue family. A fused family is trivially present in every genome, so the
core inflates for an ARTIFACTUAL reason, and a core alignment built from it
concatenates non-orthologues.

Two independent readouts, both from PPanGGOLiN's own outputs rather than a
reimplementation:

  1. mean_persistent_duplication.tsv -- the tool's own per-persistent-family
     duplication ratio and `is_single_copy_marker` flag.
  2. matrix.csv -- "No. isolates" and "Avg sequences per isolate" for every
     family, from which the >=95%-prevalence core and its copy number follow.

Pre-declared interpretation rule (written before the numbers were seen):
  if the single-copy fraction for c3_id060 falls MORE THAN 0.15 below c3_id080,
  the 0.60 core is substantially fused and the 0.60 pC3 tree is a sensitivity
  analysis only -- never the primary phylogeny.
"""
import csv
import sys
from pathlib import Path

W = Path("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
CORE_FRAC = 0.95
csv.field_size_limit(10_000_000)


def persistent_dup(path):
    """-> (n_persistent, frac_single_copy_marker, mean_duplication_ratio)"""
    n = n_sc = 0
    tot = 0.0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            tot += float(r["duplication_ratio"])
            if r["is_single_copy_marker"].strip().lower() == "true":
                n_sc += 1
    return n, (round(n_sc / n, 4) if n else ""), (round(tot / n, 4) if n else "")


def core_from_matrix(path):
    """-> (n_genomes, n_core_ge_95pct, frac_core_single_copy, mean_copies)"""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh)
        hdr = [h.strip().strip('"') for h in next(rd)]
        i_iso = hdr.index("No. isolates")
        i_avg = hdr.index("Avg sequences per isolate")
        # genome columns are everything after the fixed Roary metadata block
        n_gen = len(hdr) - (hdr.index("Avg group size nuc") + 1)
        n_core = n_single = 0
        tot = 0.0
        for rec in rd:
            if len(rec) <= i_avg:
                continue
            try:
                iso = int(float(rec[i_iso]))
                avg = float(rec[i_avg])
            except ValueError:
                continue
            if iso < CORE_FRAC * n_gen:
                continue
            n_core += 1
            tot += avg
            if abs(avg - 1.0) < 1e-9:
                n_single += 1
    if not n_core:
        return n_gen, 0, "", ""
    return n_gen, n_core, round(n_single / n_core, 4), round(tot / n_core, 3)


rows = []
for d in sorted((W / "ppanggolin").glob("*_id*")):
    mat, dup = d / "out" / "matrix.csv", d / "out" / "mean_persistent_duplication.tsv"
    if not mat.exists():
        print(f"SKIP {d.name}: no matrix.csv", file=sys.stderr)
        continue
    n_gen, n_core, frac_core_sc, mean_c = core_from_matrix(mat)
    n_pers, frac_pers_sc, mean_dup = persistent_dup(dup) if dup.exists() else ("", "", "")
    rows.append({
        "run": d.name, "n_genomes": n_gen,
        "n_core_ge_95pct": n_core,
        "frac_core_single_copy": frac_core_sc,
        "mean_copies_per_core_carrier": mean_c,
        "n_persistent": n_pers,
        "frac_persistent_single_copy_marker": frac_pers_sc,
        "mean_persistent_duplication_ratio": mean_dup,
    })
    print(rows[-1])

if not rows:
    raise SystemExit("no PPanGGOLiN runs found")

out = W / "results" / "collapse_diagnostic.tsv"
out.parent.mkdir(exist_ok=True)
with open(out, "w", newline="") as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    wtr.writeheader()
    wtr.writerows(rows)
print(f"\nwrote {out}\n")

by = {r["run"]: r for r in rows}
for setname in ("c3", "chr1"):
    lo, hi = by.get(f"{setname}_id060"), by.get(f"{setname}_id080")
    if not (lo and hi):
        continue
    for key, label in (("frac_core_single_copy", "core"),
                       ("frac_persistent_single_copy_marker", "persistent")):
        if lo[key] == "" or hi[key] == "":
            continue
        drop = hi[key] - lo[key]
        verdict = "FUSED -> sensitivity only" if drop > 0.15 else "acceptable"
        print(f"{setname:5s} {label:11s} single-copy  0.80={hi[key]:.4f}  "
              f"0.60={lo[key]:.4f}  drop={drop:+.4f}  -> {verdict}")
