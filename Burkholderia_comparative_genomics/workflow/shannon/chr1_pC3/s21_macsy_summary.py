#!/usr/bin/env python3
"""
s21 - tabulate the MacSyFinder / TXSScan results across the five B. sola pC3 replicons.

Reports BOTH:
  - complete systems that satisfied the model quorum (best_solution.tsv)
  - rejected candidates (rejected_candidates.tsv), which is where a partial or
    orphan module shows up. Scope here is pC3 alone, so a system whose mandatory
    genes live on the chromosome CANNOT satisfy quorum; reporting only the
    complete systems would turn that scoping choice into a false negative.
"""
import csv, os, re
from collections import defaultdict, Counter

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
OUT = f"{W}/results"
MS = f"{W}/screen/macsy"
SETB = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1", "MF6"]


def read_tsv(path):
    if not os.path.exists(path):
        return []
    rows, hdr = [], None
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            continue
        if len(f) < 2:
            continue
        rows.append(dict(zip(hdr, f)))
    return rows


# ---------------------------------------------------------------- complete systems
sys_rows, comp_rows = [], []
for g in SETB:
    best = read_tsv(f"{MS}/{g}/best_solution.tsv")
    bysys = defaultdict(list)
    for r in best:
        sid = r.get("sys_id") or r.get("model_fqn", "")
        bysys[(r.get("model_fqn", ""), sid)].append(r)
    for (model, sid), rs in sorted(bysys.items()):
        genes = [r.get("gene_name", "") for r in rs]
        pos = [int(r["hit_pos"]) for r in rs if r.get("hit_pos", "").isdigit()]
        sys_rows.append(dict(
            genome=g, model=model.split("/")[-1], system_id=sid,
            n_genes=len(rs), wholeness=rs[0].get("sys_wholeness", ""),
            score=rs[0].get("sys_score", ""),
            first_pos=min(pos) if pos else "", last_pos=max(pos) if pos else "",
            genes=",".join(sorted(set(genes))),
            loci=",".join(sorted({r.get("hit_id", "") for r in rs})),
        ))

with open(f"{OUT}/screen_secretion_systems.tsv", "w") as fh:
    if sys_rows:
        w = csv.DictWriter(fh, fieldnames=list(sys_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(sys_rows)

print("=== complete systems per genome (quorum satisfied) ===", flush=True)
cnt = defaultdict(Counter)
for r in sys_rows:
    cnt[r["genome"]][r["model"]] += 1
models = sorted({r["model"] for r in sys_rows})
print(f"{'genome':18s} " + " ".join(f"{m:>11s}" for m in models), flush=True)
for g in SETB:
    print(f"{g:18s} " + " ".join(f"{cnt[g][m]:>11d}" for m in models), flush=True)

# ---------------------------------------------------------------- rejected / partial
print("\n=== rejected candidates: partial or orphan modules ===", flush=True)
rej_rows = []
for g in SETB:
    rej = read_tsv(f"{MS}/{g}/rejected_candidates.tsv")
    bycand = defaultdict(list)
    for r in rej:
        bycand[(r.get("model_fqn", ""), r.get("candidate_id") or r.get("replicon", ""))].append(r)
    for (model, cid), rs in sorted(bycand.items()):
        genes = sorted({r.get("gene_name", "") for r in rs})
        loci = [r.get("hit_id", "") for r in rs]
        reason = rs[0].get("rejected_reason", "")
        rej_rows.append(dict(genome=g, model=model.split("/")[-1], candidate=cid,
                             n_genes=len(rs), genes=",".join(genes),
                             loci=",".join(loci), reason=reason))

with open(f"{OUT}/screen_secretion_rejected.tsv", "w") as fh:
    if rej_rows:
        w = csv.DictWriter(fh, fieldnames=list(rej_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rej_rows)

for r in rej_rows:
    if "T6SS" in r["model"]:
        print(f"  {r['genome']:18s} {r['model']:10s} genes={r['genes']}  loci={r['loci']}",
              flush=True)
t6 = [r for r in rej_rows if "T6SS" in r["model"]]
print(f"\npC3 replicons carrying a partial T6SS module: "
      f"{len({r['genome'] for r in t6})} / {len(SETB)}", flush=True)
print(f"pC3 replicons with a COMPLETE T6SS: "
      f"{len({r['genome'] for r in sys_rows if 'T6SS' in r['model']})} / {len(SETB)}", flush=True)
print(f"\nwrote {len(sys_rows)} system rows, {len(rej_rows)} rejected-candidate rows", flush=True)
