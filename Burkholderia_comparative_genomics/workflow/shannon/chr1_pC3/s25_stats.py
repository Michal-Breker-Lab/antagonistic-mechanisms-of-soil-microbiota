#!/usr/bin/env python3
"""s25 - association statistics and the per-residue coverage profile."""
import csv
import glob
import os

from scipy.stats import fisher_exact

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
c3 = {r["accession"]: (r["c3_present"] == "True")
      for r in csv.DictReader(open(f"{W}/results/c3_calls_all_genomes.tsv"), delimiter="\t")}
accs = sorted(os.path.basename(f)[:-4] for f in glob.glob(f"{W}/pyrodigal_prot/*.faa"))

FMT = ("qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send "
       "evalue bitscore qcovhsp qcovs slen qlen").split()
rows = [dict(zip(FMT, l.rstrip("\n").split("\t")))
        for l in open(f"{W}/search_rhs/blastp_hits.tsv")]

ct = {r["sseqid"].split("|")[0] for r in rows
      if r["qseqid"] == "CT354_CFFIHE_03684"
      and float(r["pident"]) >= 60 and float(r["qcovs"]) >= 70}
full = {r["sseqid"].split("|")[0] for r in rows if r["qseqid"] == "FULL_CFFIHE_03684"}


def tab(S, label):
    a = sum(1 for x in accs if x in S and c3.get(x))
    b = sum(1 for x in accs if x in S and not c3.get(x))
    c = sum(1 for x in accs if x not in S and c3.get(x))
    d = sum(1 for x in accs if x not in S and not c3.get(x))
    orr, p = fisher_exact([[a, b], [c, d]])
    orr_s = "inf" if orr == float("inf") else f"{orr:.2f}"
    print(f"{label}:  present on {a} pC3+ and {b} pC3- genomes "
          f"(of {a+c} pC3+ / {b+d} pC3-)  OR={orr_s}  Fisher p={p:.3g}")
    return dict(feature=label.strip(), n_pC3pos_with=a, n_pC3neg_with=b,
                n_pC3pos_total=a + c, n_pC3neg_total=b + d,
                odds_ratio=orr_s, fisher_p=f"{p:.3g}")


print(f"genomes searched: {len(accs)}   pC3-positive: {sum(1 for a in accs if c3.get(a))}\n")
s1 = tab(ct, "CT354 warhead")
s2 = tab(full, "FULL RHS protein")
with open(f"{W}/results/rhs_pC3_association.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(s1.keys()), delimiter="\t")
    w.writeheader(); w.writerows([s1, s2])

# ---- per-residue coverage profile of the query by CT-NEGATIVE genomes
L = 2867
cov = [0] * (L + 1)
for r in rows:
    if r["qseqid"] != "FULL_CFFIHE_03684":
        continue
    if r["sseqid"].split("|")[0] in ct:
        continue
    for i in range(int(r["qstart"]), int(r["qend"]) + 1):
        cov[i] += 1
with open(f"{W}/results/rhs_query_coverage_profile.tsv", "w") as fh:
    fh.write("query_residue\tn_hsps_ct_negative\n")
    for i in range(1, L + 1):
        fh.write(f"{i}\t{cov[i]}\n")
last = max(i for i in range(1, L + 1) if cov[i] > 0)
print(f"\nlast query residue covered by any CT-negative genome: {last}")
print(f"residues {last+1}-{L} ({L-last} aa) are exclusive to the {len(ct)} carriers")

# ---- DEGQ control audit
degq = {}
for r in rows:
    if r["qseqid"] != "DEGQ_CFFIHE_02734":
        continue
    a = r["sseqid"].split("|")[0]
    if a not in degq or float(r["bitscore"]) > float(degq[a]["bitscore"]):
        degq[a] = r
miss = [(a, float(v["pident"]), float(v["qcovs"])) for a, v in degq.items()
        if not (float(v["pident"]) >= 60 and float(v["qcovs"]) >= 70)]
print(f"\nDEGQ control: {len(degq)}/771 genomes hit; {len(miss)} below tier1:")
for a, p, q in sorted(miss):
    print(f"   {a}  id={p:.1f}  cov={q:.1f}")
