#!/usr/bin/env python3
"""
s22 - tabulate antiSMASH regions on the five B. sola pC3 replicons.

Reports region type, coordinates, size, and the best knownclusterblast hit with
its similarity. A knownclusterblast hit is a SIMILARITY score to a characterised
cluster, not an identification - wording in the report must reflect that.
"""
import csv, glob, json, os, re
from collections import Counter, defaultdict

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
OUT = f"{W}/results"
AS = f"{W}/screen/antismash"
SETB = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1", "MF6"]

rows = []
for g in SETB:
    js = glob.glob(f"{AS}/{g}/*_pC3.json")
    if not js:
        print(f"WARNING: no JSON for {g}")
        continue
    data = json.load(open(js[0]))
    for rec in data.get("records", []):
        reclen = len(rec.get("seq", {}).get("data", "")) or rec.get("length", 0)
        # regions live in 'areas'
        for i, area in enumerate(rec.get("areas", []), 1):
            prods = area.get("products", [])
            start, end = area.get("start", 0), area.get("end", 0)
            rows.append(dict(genome=g, region=i, products=";".join(prods),
                             start=start, end=end, length_bp=end - start,
                             best_known_hit="", n_proteins_hit="", cumulative_blast=""))

# knownclusterblast: parse the per-region txt reports for the top hit
for g in SETB:
    for f in sorted(glob.glob(f"{AS}/{g}/knownclusterblast/*_c*.txt")):
        m = re.search(r"_c(\d+)\.txt$", f)
        if not m:
            continue
        reg = int(m.group(1))
        txt = open(f, errors="replace").read()
        # "Significant hits:\n1. BGC0000001\tname" then later "Cumulative BLAST score"
        hit, nprot, score = "", "", ""
        sec = re.search(r"Significant hits:\s*\n(.*?)\n\s*\n", txt, re.S)
        if sec:
            first = sec.group(1).strip().split("\n")[0]
            hit = re.sub(r"^\d+\.\s*", "", first).strip().replace("\t", " ")
        # first entry of the Details block corresponds to the top hit
        det = re.search(r"Details:.*?Number of proteins with BLAST hits to this cluster:\s*(\d+)"
                        r".*?Cumulative BLAST score:\s*([\d.]+)", txt, re.S)
        if det:
            nprot, score = det.group(1), det.group(2)
        for r in rows:
            if r["genome"] == g and r["region"] == reg:
                r["best_known_hit"] = hit
                r["n_proteins_hit"] = nprot
                r["cumulative_blast"] = score

with open(f"{OUT}/screen_antismash_regions.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

print("=== antiSMASH regions per pC3 replicon ===", flush=True)
per = Counter(r["genome"] for r in rows)
for g in SETB:
    print(f"  {g:18s} {per[g]} regions", flush=True)

print("\n=== region types across the clade ===", flush=True)
tc = Counter()
for r in rows:
    for p in r["products"].split(";"):
        if p:
            tc[p] += 1
for p, n in tc.most_common():
    ing = len({r["genome"] for r in rows if p in r["products"].split(";")})
    print(f"  {p:24s} {n:3d} regions, in {ing}/5 genomes", flush=True)

print("\n=== per-genome detail ===", flush=True)
for g in SETB:
    print(f"-- {g}", flush=True)
    for r in [x for x in rows if x["genome"] == g]:
        kb = f"  ~ {r['best_known_hit'][:46]} [{r['n_proteins_hit']} prot, score {r['cumulative_blast']}]" if r['best_known_hit'] else ''
        print(f"   region {r['region']:>2} {r['products']:<22} "
              f"{r['length_bp']:>7,} bp{kb}", flush=True)

tot = sum(r["length_bp"] for r in rows)
print(f"\ntotal BGC bp across 5 replicons: {tot:,}", flush=True)
for g in SETB:
    b = sum(r["length_bp"] for r in rows if r["genome"] == g)
    print(f"  {g:18s} {b:>9,} bp in BGC regions", flush=True)
