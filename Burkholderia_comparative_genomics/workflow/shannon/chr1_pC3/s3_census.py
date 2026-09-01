#!/usr/bin/env python3
"""Stage 3 - replicon census. Answers "3 large replicons vs 2" from replicon
SIZES only, with no reliance on how NCBI labelled each molecule.

Two outputs matter:
  1. the census itself, at four size thresholds (sensitivity analysis, D1);
  2. the disagreement rate between the size-based call and NCBI's own
     chromosome/plasmid labelling, which is a reportable finding about database
     curation rather than a nuisance.
"""
import collections
import csv
import os
import sys

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
MD, RES = f"{W}/metadata", f"{W}/results"
THRESHOLDS = [200_000, 300_000, 500_000, 1_000_000]
PRIMARY = 300_000


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    os.makedirs(RES, exist_ok=True)
    asm = {r["accession"]: r for r in read_tsv(f"{MD}/assemblies.tsv")}
    reps = read_tsv(f"{MD}/replicons.tsv")
    print(f"assemblies: {len(asm)}   replicon rows: {len(reps)}")

    by_asm = collections.defaultdict(list)
    skipped_role = 0
    for r in reps:
        if r.get("role") and r["role"] != "assembled-molecule":
            skipped_role += 1
            continue
        try:
            r["length"] = int(r["length"])
        except (ValueError, KeyError):
            continue
        by_asm[r["accession"]].append(r)
    print(f"skipped non-assembled-molecule rows: {skipped_role}")
    print(f"assemblies with replicon data: {len(by_asm)}")

    missing = sorted(set(asm) - set(by_asm))
    if missing:
        print(f"WARNING: {len(missing)} assemblies have NO replicon data")

    rows = []
    for acc, rl in by_asm.items():
        a = asm.get(acc, {})
        rl.sort(key=lambda x: -x["length"])
        sizes = [x["length"] for x in rl]
        row = {
            "accession": acc,
            "organism_name": a.get("organism_name", ""),
            "genus": (a.get("organism_name", "") or " ").split()[0],
            "strain": a.get("strain", ""),
            "query_genus": a.get("query_genus", ""),
            "total_length": a.get("total_length", ""),
            "gc_percent": a.get("gc_percent", ""),
            "n_replicons_total": len(rl),
            "replicon_sizes": ";".join(str(s) for s in sizes),
            "largest": sizes[0] if sizes else 0,
            "second": sizes[1] if len(sizes) > 1 else 0,
            "third": sizes[2] if len(sizes) > 2 else 0,
            "sequencing_tech": a.get("sequencing_tech", ""),
            "checkm_completeness": a.get("checkm_completeness", ""),
            "checkm_contamination": a.get("checkm_contamination", ""),
            "ncbi_n_chromosomes": a.get("ncbi_n_chromosomes", ""),
        }
        for t in THRESHOLDS:
            row[f"n_large_{t//1000}kb"] = sum(1 for s in sizes if s >= t)
        n = row[f"n_large_{PRIMARY//1000}kb"]
        row["architecture"] = (f"{n}_large" if n <= 3 else "4+_large")
        # NCBI's own labelling, for the disagreement cross-check
        n_chr_lab = sum(1 for x in rl if x.get("molecule_type") == "Chromosome")
        n_pl_lab = sum(1 for x in rl if x.get("molecule_type") == "Plasmid")
        row["ncbi_labelled_chromosomes"] = n_chr_lab
        row["ncbi_labelled_plasmids"] = n_pl_lab
        # what does NCBI call the third-largest replicon?
        row["third_ncbi_label"] = rl[2].get("molecule_type", "") if len(rl) > 2 else ""
        row["third_chr_name"] = rl[2].get("chr_name", "") if len(rl) > 2 else ""
        row["size_vs_label_agree"] = str(n == n_chr_lab)
        # QC: "Complete Genome" is a submitter assertion NCBI does not enforce.
        # One genome in this set is labelled complete at 8.09% CheckM
        # completeness -- a single chromosome-2-sized fragment. Such genomes
        # would register as false c3 absences, so flag rather than trust.
        try:
            comp = float(a.get("checkm_completeness") or "nan")
        except ValueError:
            comp = float("nan")
        try:
            cont = float(a.get("checkm_contamination") or "nan")
        except ValueError:
            cont = float("nan")
        row["qc_pass"] = str(not (comp < 95.0 or cont > 5.0))
        row["qc_reason"] = ("completeness<95" if comp < 95.0 else
                            "contamination>5" if cont > 5.0 else "")
        rows.append(row)

    rows.sort(key=lambda r: (r["organism_name"], r["accession"]))
    cols = list(rows[0].keys())
    with open(f"{RES}/replicon_census.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {RES}/replicon_census.tsv  ({len(rows)} genomes)")

    # ---- threshold sensitivity ----
    print("\n=== THRESHOLD SENSITIVITY: n genomes by count of large replicons ===")
    print(f"{'threshold':>10} " + " ".join(f"{k:>8}" for k in
                                           ["1_large", "2_large", "3_large", "4+_large"]))
    for t in THRESHOLDS:
        c = collections.Counter()
        for r in rows:
            n = r[f"n_large_{t//1000}kb"]
            c[f"{n}_large" if n <= 3 else "4+_large"] += 1
        print(f"{t//1000:>8}kb " + " ".join(f"{c.get(k,0):>8}" for k in
                                            ["1_large", "2_large", "3_large", "4+_large"]))

    # ---- label agreement ----
    agree = sum(1 for r in rows if r["size_vs_label_agree"] == "True")
    print(f"\n=== SIZE-BASED vs NCBI LABEL (at {PRIMARY//1000} kb) ===")
    print(f"agree: {agree}/{len(rows)} ({100*agree/len(rows):.1f}%)")
    lab = collections.Counter(r["third_ncbi_label"] for r in rows if r["third"] >= PRIMARY)
    print("NCBI's label for the third-largest replicon, where it is >=300 kb:")
    for k, v in lab.most_common():
        print(f"   {k or '(none)':<14} {v}")

    # ---- census with and without the QC filter ----
    passing = [r for r in rows if r["qc_pass"] == "True"]
    print(f"\n=== QC FILTER (CheckM completeness >=95%, contamination <=5%) ===")
    print(f"pass: {len(passing)}/{len(rows)}   fail: {len(rows)-len(passing)}")
    reasons = collections.Counter(r["qc_reason"] for r in rows if r["qc_reason"])
    for k, v in reasons.most_common():
        print(f"   {v:>4}  {k}")
    print(f"\n{'set':<18}" + "".join(f"{k:>10}" for k in
          ["1_large", "2_large", "3_large", "4+_large"]))
    for label, rs in (("all genomes", rows), ("QC-passing", passing)):
        c = collections.Counter(r["architecture"] for r in rs)
        print(f"{label+f' (n={len(rs)})':<18}" + "".join(
            f"{c.get(k,0):>10}" for k in ["1_large", "2_large", "3_large", "4+_large"]))
    pc = collections.Counter(r["architecture"] for r in passing)
    ac = collections.Counter(r["architecture"] for r in rows)
    if ac["3_large"]:
        print(f"\n3-replicon share: all={100*ac['3_large']/len(rows):.1f}%  "
              f"QC-passing={100*pc['3_large']/len(passing):.1f}%")

    # ---- by genus and by species ----
    print("\n=== ARCHITECTURE BY GENUS (300 kb) ===")
    g = collections.defaultdict(collections.Counter)
    for r in rows:
        g[r["genus"]][r["architecture"]] += 1
    print(f"{'genus':<20}" + "".join(f"{k:>10}" for k in
                                     ["1_large", "2_large", "3_large", "4+_large"]) + "   n")
    for gen in sorted(g, key=lambda x: -sum(g[x].values())):
        c = g[gen]
        print(f"{gen:<20}" + "".join(f"{c.get(k,0):>10}" for k in
                                     ["1_large", "2_large", "3_large", "4+_large"])
              + f"   {sum(c.values())}")

    sp = collections.defaultdict(collections.Counter)
    for r in rows:
        sp[r["organism_name"]][r["architecture"]] += 1
    with open(f"{RES}/census_by_species.tsv", "w") as fh:
        fh.write("organism_name\tn_genomes\t1_large\t2_large\t3_large\t4+_large\tvariable\n")
        for name in sorted(sp):
            c = sp[name]
            tot = sum(c.values())
            var = sum(1 for k in c if c[k] > 0) > 1
            fh.write(f"{name}\t{tot}\t" + "\t".join(str(c.get(k, 0)) for k in
                     ["1_large", "2_large", "3_large", "4+_large"]) + f"\t{var}\n")

    # The headline: species where c3 presence VARIES between conspecific genomes.
    print("\n=== SPECIES WITH WITHIN-SPECIES ARCHITECTURE VARIATION (n>=3) ===")
    hits = [(n, sp[n]) for n in sp
            if sum(sp[n].values()) >= 3 and sum(1 for k in sp[n] if sp[n][k] > 0) > 1]
    hits.sort(key=lambda x: -sum(x[1].values()))
    if not hits:
        print("  none")
    for name, c in hits[:30]:
        print(f"  {name:<45} n={sum(c.values()):<4} " +
              " ".join(f"{k}={c[k]}" for k in sorted(c) if c[k]))
    print(f"\n({len(hits)} species show variation)")
    print("STAGE3_DONE")


if __name__ == "__main__":
    sys.exit(main())
