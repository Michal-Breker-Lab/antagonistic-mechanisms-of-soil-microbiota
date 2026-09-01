#!/usr/bin/env python3
"""
s20 - apply the PRE-DECLARED signature list to the InterProScan output.

Two things this script deliberately does:

  1. Reports EVERY declared category, including the ones that score zero, so the
     denominator of the screen is visible and a null is legible as a null.
  2. Separates accession-based evidence (exact InterPro/Pfam/NCBIfam accession,
     high specificity) from free-text description matches (low specificity, prone
     to substring false positives). This split was declared before results were
     seen; only accession evidence is treated as a confident call.

Anything interesting found OUTSIDE the declared list is written to a separate
"post_hoc" table and must be labelled as such wherever it is reported.
"""
import csv, os, re, sys
from collections import defaultdict, Counter

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
OUT = f"{W}/results"
IPS_TSV = f"{W}/screen/ips/pC3_ips.tsv"
SIGS = f"{W}/effector_toxin_signatures.tsv"

# ------------------------------------------------------------------ load rules
rules = []          # (category, match_type, pattern_lower, note)
for line in open(SIGS):
    line = line.rstrip("\n")
    if not line.strip() or line.lstrip().startswith("#"):
        continue
    f = line.split("\t")
    if len(f) < 3:
        continue
    rules.append((f[0].strip(), f[1].strip(), f[2].strip().lower(),
                  f[3].strip() if len(f) > 3 else ""))
CATS = []
for c, *_ in rules:
    if c not in CATS:
        CATS.append(c)
print(f"loaded {len(rules)} declared rules across {len(CATS)} categories", flush=True)

# ------------------------------------------------------------------ load IPS
# columns: acc, md5, len, analysis, sig_acc, sig_desc, start, stop, score,
#          status, date, ipr_acc, ipr_desc, GO, pathways
fam_hits = defaultdict(list)     # family -> list of hit dicts
fam_meta = {}
n_rows = 0
if not os.path.exists(IPS_TSV):
    sys.exit(f"missing {IPS_TSV}")
with open(IPS_TSV) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        n_rows += 1
        acc = f[0]
        parts = acc.split("|")
        fam = parts[0]
        fam_meta.setdefault(fam, dict(locus=parts[1] if len(parts) > 1 else "",
                                      set=parts[2] if len(parts) > 2 else "",
                                      cog=parts[3].replace("COG=", "") if len(parts) > 3 else ""))
        fam_hits[fam].append(dict(
            analysis=f[3], sig_acc=f[4], sig_desc=f[5],
            ipr_acc=f[11] if len(f) > 11 else "", ipr_desc=f[12] if len(f) > 12 else "",
            go=f[13] if len(f) > 13 else ""))
print(f"InterProScan: {n_rows} rows over {len(fam_hits)} families", flush=True)

# Denominator must be EVERY submitted family, not just those with a hit - otherwise
# proteins that matched nothing vanish and coverage is overstated.
ALL = {}
for line in open(f"{W}/screen/ips/pC3_clade_pangenome.faa"):
    if line.startswith(">"):
        pr = line[1:].split()[0].split("|")
        ALL[pr[0]] = dict(locus=pr[1] if len(pr) > 1 else "",
                          set=pr[2] if len(pr) > 2 else "",
                          cog=pr[3].replace("COG=", "") if len(pr) > 3 else "")
for f, m in ALL.items():
    fam_meta.setdefault(f, m)
print(f"submitted families: {len(ALL)}; with >=1 signature: {len(fam_hits)}; "
      f"with NO signature at all: {len(ALL) - len(fam_hits)}", flush=True)

# ------------------------------------------------------------------ coverage
# The headline question for the unknowns: did InterProScan assign anything at all?
UNINFORMATIVE = {"MobiDBLite", "Coils", "Phobius", "TMHMM", "SignalP"}
cov = Counter()
by_set = defaultdict(lambda: Counter())
for fam in ALL:
    hits = fam_hits.get(fam, [])
    s = fam_meta[fam]["set"]
    cog = fam_meta[fam]["cog"]
    informative = [h for h in hits if h["analysis"] not in UNINFORMATIVE]
    with_ipr = [h for h in informative if h["ipr_acc"].startswith("IPR")]
    key = (s, "noCOG" if cog in ("-", "") else "COG")
    by_set[key]["families"] += 1
    by_set[key]["any_signature"] += bool(informative)
    by_set[key]["interpro_entry"] += bool(with_ipr)

print("\n=== InterProScan coverage ===", flush=True)
print(f"{'set':12s} {'COG?':6s} {'families':>9s} {'any sig':>9s} {'IPR entry':>10s}", flush=True)
rows = []
for (s, c), v in sorted(by_set.items()):
    print(f"{s:12s} {c:6s} {v['families']:9d} {v['any_signature']:9d} "
          f"({100*v['any_signature']/v['families']:.1f}%) {v['interpro_entry']:6d} "
          f"({100*v['interpro_entry']/v['families']:.1f}%)", flush=True)
    rows.append(dict(gene_set=s, cog_status=c, n_families=v["families"],
                     n_any_signature=v["any_signature"],
                     pct_any_signature=round(100*v["any_signature"]/v["families"], 1),
                     n_interpro_entry=v["interpro_entry"],
                     pct_interpro_entry=round(100*v["interpro_entry"]/v["families"], 1)))
with open(f"{OUT}/screen_ips_coverage.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ apply rules
hits_by_cat = defaultdict(lambda: {"acc": set(), "desc": set()})
detail = []
for fam, hits in fam_hits.items():
    for h in hits:
        hay_acc = {h["sig_acc"].lower(), h["ipr_acc"].lower()}
        hay_desc = f"{h['sig_desc']} {h['ipr_desc']}".lower()
        for cat, mtype, pat, note in rules:
            matched = False
            if mtype == "acc" and pat in hay_acc:
                matched, ev = True, "accession"
            elif mtype == "desc" and pat and pat in hay_desc:
                matched, ev = True, "description"
            if matched:
                hits_by_cat[cat]["acc" if ev == "accession" else "desc"].add(fam)
                detail.append(dict(category=cat, evidence=ev, pattern=pat, family=fam,
                                   gene_set=fam_meta[fam]["set"], locus=fam_meta[fam]["locus"],
                                   cog=fam_meta[fam]["cog"], analysis=h["analysis"],
                                   sig_acc=h["sig_acc"], sig_desc=h["sig_desc"],
                                   ipr_acc=h["ipr_acc"], ipr_desc=h["ipr_desc"]))

print("\n=== declared-category results (zero rows included by design) ===", flush=True)
print(f"{'category':20s} {'acc-based':>10s} {'desc-only':>10s} {'core':>6s} {'acc':>5s}", flush=True)
summ = []
for cat in CATS:
    a = hits_by_cat[cat]["acc"]
    d = hits_by_cat[cat]["desc"] - a
    allf = a | d
    n_core = sum(1 for f in allf if fam_meta[f]["set"] == "core")
    print(f"{cat:20s} {len(a):10d} {len(d):10d} {n_core:6d} {len(allf)-n_core:5d}", flush=True)
    summ.append(dict(category=cat, n_families_accession_evidence=len(a),
                     n_families_description_only=len(d), n_total=len(allf),
                     n_in_core=n_core, n_in_accessory=len(allf)-n_core))
with open(f"{OUT}/screen_category_summary.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(summ)

with open(f"{OUT}/screen_category_hits.tsv", "w") as fh:
    if detail:
        w = csv.DictWriter(fh, fieldnames=list(detail[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(detail)
    else:
        fh.write("category\tevidence\tpattern\tfamily\tgene_set\tlocus\tcog\t"
                 "analysis\tsig_acc\tsig_desc\tipr_acc\tipr_desc\n")
print(f"\nwrote {len(detail)} hit rows", flush=True)

# ------------------------------------------------------------------ post hoc
# Families that InterProScan characterised but which the declared list never
# touched - reported separately and explicitly labelled post hoc.
declared_fams = set()
for cat in CATS:
    declared_fams |= hits_by_cat[cat]["acc"] | hits_by_cat[cat]["desc"]
nocog_rescued = []
for fam, hits in fam_hits.items():
    if fam_meta[fam]["cog"] not in ("-", ""):
        continue
    if fam in declared_fams:
        continue
    best = [h for h in hits if h["ipr_acc"].startswith("IPR")]
    if best:
        nocog_rescued.append((fam, fam_meta[fam]["set"], best[0]["ipr_acc"],
                              best[0]["ipr_desc"]))
with open(f"{OUT}/screen_posthoc_rescued.tsv", "w") as fh:
    fh.write("family\tgene_set\tinterpro\tinterpro_description\n")
    for r in sorted(nocog_rescued, key=lambda x: x[3]):
        fh.write("\t".join(r) + "\n")
print(f"un-COG'd families given an InterPro entry but outside the declared list: "
      f"{len(nocog_rescued)}  (post hoc, reported separately)", flush=True)
