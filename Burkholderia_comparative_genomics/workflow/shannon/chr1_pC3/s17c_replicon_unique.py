#!/usr/bin/env python3
"""
s17c - correct the replicon-uniqueness test.

s17b compared pC3 against chromosome 1 only, but MF6 has 7,075 CDS while
chr1 + pC3 account for 4,262 - the second replicon (and any plasmids) were
excluded. "Complete on pC3, absent from chr1" therefore does NOT establish
that a module is unique to pC3.

Here the comparison set is the WHOLE REST OF THE GENOME (all MF6 CDS minus
those on pC3), which is the test that actually supports a uniqueness claim.
"""
import csv, re, json
from collections import Counter

_SRC = open("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/s17_functional.py").read()
assert "# 9. Redundancy" in _SRC
exec(_SRC.split("# 9. Redundancy")[0])

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
mf6_all = {lt for lt in ANN if lt.startswith(next(iter(mf6_c3)).split("_")[0])}
print(f"MF6 CDS: total {len(mf6_all)}, pC3 {len(mf6_c3)}, chr1 {len(mf6_chr1)}, "
      f"other replicons {len(mf6_all - mf6_c3 - mf6_chr1)}", flush=True)

def kos_of(loci):
    s = set()
    for lt in loci:
        a = ANN.get(lt)
        if a:
            s |= a["kos"]
    return s

ko_c3 = kos_of(mf6_c3)
ko_rest = kos_of(mf6_all - mf6_c3)
ko_all = kos_of(mf6_all)
print(f"KOs: pC3 {len(ko_c3)}, rest-of-genome {len(ko_rest)}, whole genome {len(ko_all)}",
      flush=True)

ko_only_c3 = ko_c3 - ko_rest
print(f"KOs found ONLY on pC3 (absent from the rest of the genome): {len(ko_only_c3)}",
      flush=True)

C_c3 = stepwise_completeness(ko_c3)
C_rest = stepwise_completeness(ko_rest)
C_all = stepwise_completeness(ko_all)

n75 = lambda C: sum(1 for m in C if C[m][2] >= MODULE_THRESHOLD)
print(f"modules >=75% complete: pC3 {n75(C_c3)}, rest-of-genome {n75(C_rest)}, "
      f"whole genome {n75(C_all)}", flush=True)

print("\n=== modules complete on pC3 AND NOT complete on the rest of the genome ===",
      flush=True)
uniq = []
for m in sorted(MODULES):
    if m not in C_c3:
        continue
    a, b = C_c3[m][2], C_rest[m][2]
    if a >= MODULE_THRESHOLD and b < MODULE_THRESHOLD:
        uniq.append((m, MODULES[m]["name"], a, b))
        print(f"  {m}  {MODULES[m]['name'][:66]:66s} pC3 {a:.2f} / rest {b:.2f}", flush=True)
if not uniq:
    print("  (none)", flush=True)

print("\n=== modules where pC3 raises whole-genome completeness (contributes steps) ===",
      flush=True)
contrib = []
for m in sorted(MODULES):
    if m not in C_c3:
        continue
    if C_all[m][2] > C_rest[m][2] + 1e-9:
        contrib.append((m, MODULES[m]["name"], C_rest[m][2], C_all[m][2], C_c3[m][2]))
for m, nm, r, a, c in sorted(contrib, key=lambda x: x[3] - x[2], reverse=True)[:20]:
    print(f"  {m}  {nm[:60]:60s} rest {r:.2f} -> genome {a:.2f}  (pC3 alone {c:.2f})",
          flush=True)
print(f"  total modules to which pC3 contributes: {len(contrib)}", flush=True)

with open(f"{OUT}/setB_pC3_module_contribution.tsv", "w") as fh:
    fh.write("module\tname\tclass\tn_steps\tpC3_alone\trest_of_genome\twhole_genome\t"
             "pC3_unique_complete\tpC3_raises_completeness\n")
    for m in sorted(MODULES):
        if m not in C_c3:
            continue
        a, r, g = C_c3[m][2], C_rest[m][2], C_all[m][2]
        fh.write(f"{m}\t{MODULES[m]['name']}\t{MODULES[m]['klass']}\t{C_c3[m][1]}\t"
                 f"{a:.4f}\t{r:.4f}\t{g:.4f}\t"
                 f"{int(a >= MODULE_THRESHOLD and r < MODULE_THRESHOLD)}\t"
                 f"{int(g > r + 1e-9)}\n")

# KOs unique to pC3, with their products - the concrete gene-level answer
print("\n=== the pC3-only KOs (gene products) ===", flush=True)
rows = []
for lt in sorted(mf6_c3):
    a = ANN.get(lt)
    if not a:
        continue
    hits = a["kos"] & ko_only_c3
    if hits:
        rows.append((lt, ",".join(sorted(hits)), a["gene"], a["product"]))
for lt, k, g, p in rows[:40]:
    print(f"  {lt}  {k:10s} {g:12s} {p[:66]}", flush=True)
print(f"  total pC3 CDS carrying a pC3-only KO: {len(rows)}", flush=True)

with open(f"{OUT}/setB_pC3_unique_kos.tsv", "w") as fh:
    fh.write("locus_tag\tko\tgene\tproduct\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

json.dump(dict(mf6_cds_total=len(mf6_all), mf6_pC3=len(mf6_c3), mf6_chr1=len(mf6_chr1),
               mf6_other_replicons=len(mf6_all - mf6_c3 - mf6_chr1),
               kos_pC3=len(ko_c3), kos_rest=len(ko_rest), kos_only_pC3=len(ko_only_c3),
               modules75_pC3=n75(C_c3), modules75_rest=n75(C_rest),
               modules75_genome=n75(C_all),
               modules_unique_complete_pC3=[u[0] for u in uniq],
               n_modules_pC3_contributes=len(contrib)),
          open(f"{OUT}/setB_pC3_uniqueness.json", "w"), indent=2)
print("\ndone", flush=True)
