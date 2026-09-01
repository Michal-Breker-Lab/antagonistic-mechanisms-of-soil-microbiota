#!/usr/bin/env python3
"""
s17b - follow-ups to s17_functional.py:
  (a) what ARE the 71% of pC3-core families with no COG assignment?
  (b) KEGG module completeness at WHOLE-REPLICON level (fairer than core-only)
  (c) modules carried on pC3 but not chr1, and vice versa
"""
import csv, os, re, json
from collections import Counter, defaultdict

_SRC = open("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/s17_functional.py").read()
assert "# 9. Redundancy" in _SRC, "marker moved - fix the split"
exec(_SRC.split("# 9. Redundancy")[0])   # everything up to the DIAMOND step

# ---------------------------------------------------------------- (a)
print("\n=== (a) the un-COG'd majority of the pC3 clade core ===", flush=True)
nocog = [k for k in SETS["pC3_clade_core"] if not LAB[k]["cats"]]
print(f"families with no COG: {len(nocog)} / {len(SETS['pC3_clade_core'])}", flush=True)

def bucket(p):
    pl = p.lower()
    if "hypothetical" in pl:
        return "hypothetical protein"
    if re.search(r"\bduf\d+", pl) or "uncharacterized" in pl or "uncharacterised" in pl:
        return "DUF / uncharacterised domain"
    if "domain-containing protein" in pl:
        return "domain-containing, no specific function"
    return "named product"

bk = Counter(bucket(LAB[k]["product"]) for k in nocog)
tot = len(nocog)
for b, n in bk.most_common():
    print(f"  {b:44s} {n:5d}  ({100*n/tot:5.1f}%)", flush=True)

print("\n  most frequent products among no-COG pC3-core families:", flush=True)
for p, n in Counter(LAB[k]["product"] for k in nocog).most_common(25):
    print(f"    {n:4d}  {p[:88]}", flush=True)

# same buckets for the COG-annotated part and for chr1, as reference
print("\n  reference - product buckets by set:", flush=True)
bucket_rows = []
for sname in ("pC3_clade_core", "pC3_clade_accessory", "chr1_clade_core"):
    c = Counter(bucket(LAB[k]["product"]) for k in SETS[sname])
    n = len(SETS[sname])
    print(f"    {sname}: " + ", ".join(f"{b} {100*v/n:.1f}%" for b, v in c.most_common()), flush=True)
    for b, v in c.items():
        bucket_rows.append(dict(gene_set=sname, product_bucket=b, n_families=v,
                                total=n, pct=round(100*v/n, 2)))
with open(f"{OUT}/setB_product_buckets.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(bucket_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(bucket_rows)

# ---------------------------------------------------------------- (b)
print("\n=== (b) whole-replicon KEGG module completeness (MF6) ===", flush=True)

def gff_loci(path):
    s = set()
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS":
            m = re.search(r"locus_tag=([^;]+)", f[8])
            if m:
                s.add(m.group(1))
    return s

mf6_c3 = gff_loci(f"{W}/pangenome/gff_c3/MF6.gff3")
mf6_chr1 = gff_loci(f"{W}/pangenome/gff_chr1/MF6.gff3")

def kos_of(loci):
    s = set()
    for lt in loci:
        a = ANN.get(lt)
        if a:
            s |= a["kos"]
    return s

WHOLE = {
    "MF6_pC3_whole":  kos_of(mf6_c3),
    "MF6_chr1_whole": kos_of(mf6_chr1),
    "MF6_genome":     kos_of(mf6_c3 | mf6_chr1),
    "pC3_clade_core": KOSETS["pC3_clade_core"],
    "chr1_clade_core": KOSETS["chr1_clade_core"],
}
for k, v in WHOLE.items():
    print(f"  {k:18s} {len(v):5d} distinct KOs", flush=True)

COMPW = {k: stepwise_completeness(v) for k, v in WHOLE.items()}
for k in WHOLE:
    n = sum(1 for m in COMPW[k] if COMPW[k][m][2] >= MODULE_THRESHOLD)
    n50 = sum(1 for m in COMPW[k] if COMPW[k][m][2] >= 0.5)
    print(f"  {k:18s} modules >=75%: {n:4d}   >=50%: {n50:4d}", flush=True)

with open(f"{OUT}/setB_module_completeness_whole.tsv", "w") as fh:
    cols = list(WHOLE.keys())
    fh.write("module\tname\tclass\tn_steps\t" +
             "\t".join(f"{c}_completeness" for c in cols) + "\n")
    for m in sorted(MODULES):
        if m not in COMPW[cols[0]]:
            continue
        vals = [f"{COMPW[c][m][2]:.4f}" for c in cols]
        fh.write(f"{m}\t{MODULES[m]['name']}\t{MODULES[m]['klass']}\t"
                 f"{COMPW[cols[0]][m][1]}\t" + "\t".join(vals) + "\n")

# ---------------------------------------------------------------- (c)
print("\n=== (c) modules complete on pC3 but NOT on chr1 (MF6, whole replicons) ===",
      flush=True)
only_c3, only_chr1 = [], []
for m in sorted(MODULES):
    if m not in COMPW["MF6_pC3_whole"]:
        continue
    a = COMPW["MF6_pC3_whole"][m][2]
    b = COMPW["MF6_chr1_whole"][m][2]
    if a >= MODULE_THRESHOLD and b < MODULE_THRESHOLD:
        only_c3.append((m, MODULES[m]["name"], a, b, MODULES[m]["klass"]))
    if b >= MODULE_THRESHOLD and a < MODULE_THRESHOLD:
        only_chr1.append((m, MODULES[m]["name"], a, b, MODULES[m]["klass"]))

print(f"  pC3-only complete modules: {len(only_c3)}", flush=True)
for m, nm, a, b, cl in only_c3:
    print(f"    {m}  {nm[:70]:70s}  pC3 {a:.2f} / chr1 {b:.2f}", flush=True)
print(f"  chr1-only complete modules: {len(only_chr1)} (showing 15)", flush=True)
for m, nm, a, b, cl in only_chr1[:15]:
    print(f"    {m}  {nm[:70]:70s}  pC3 {a:.2f} / chr1 {b:.2f}", flush=True)

with open(f"{OUT}/setB_modules_replicon_specific.tsv", "w") as fh:
    fh.write("module\tname\tclass\tpC3_completeness\tchr1_completeness\tcomplete_on\n")
    for m, nm, a, b, cl in only_c3:
        fh.write(f"{m}\t{nm}\t{cl}\t{a:.4f}\t{b:.4f}\tpC3\n")
    for m, nm, a, b, cl in only_chr1:
        fh.write(f"{m}\t{nm}\t{cl}\t{a:.4f}\t{b:.4f}\tchr1\n")

# top near-complete modules on the pC3 clade core specifically
print("\n=== pC3 CLADE CORE: modules >=50% complete ===", flush=True)
rows = sorted(((COMPW["pC3_clade_core"][m][2], m) for m in COMPW["pC3_clade_core"]
               if COMPW["pC3_clade_core"][m][2] >= 0.5), reverse=True)
for v, m in rows[:25]:
    print(f"    {v:.2f}  {m}  {MODULES[m]['name'][:74]}", flush=True)

json.dump(dict(
    n_nocog=len(nocog),
    product_buckets_pC3_core=dict(bk),
    modules_75={k: sum(1 for m in COMPW[k] if COMPW[k][m][2] >= MODULE_THRESHOLD)
                for k in WHOLE},
    modules_50={k: sum(1 for m in COMPW[k] if COMPW[k][m][2] >= 0.5) for k in WHOLE},
    n_modules_pC3_only=len(only_c3), n_modules_chr1_only=len(only_chr1),
), open(f"{OUT}/setB_functional_detail.json", "w"), indent=2)
print("\ndone", flush=True)
