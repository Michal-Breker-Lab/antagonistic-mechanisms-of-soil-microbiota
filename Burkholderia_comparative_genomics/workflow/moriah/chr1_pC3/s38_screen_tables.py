#!/usr/bin/env python3
"""Summarise the three Set B screens into the small tables fig8 consumes.

Runs ON MORIAH: the InterProScan TSV is ~136 MB and the antiSMASH JSONs are
large, so they are reduced here and only the summaries travel to the Drive.

Outputs (schemas reproduce the retained originals):
  screen_antismash_regions.tsv    genome, region, products, start, end, ...
  screen_secretion_systems.tsv    one row per detected system
  screen_secretion_rejected.tsv   one row per rejected candidate
  screen_ips_signatures.tsv       per-protein signature rollup (feeds coverage)
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--screens", required=True, type=Path)
ap.add_argument("--members", required=True, type=Path)
ap.add_argument("--outdir", required=True, type=Path)
a = ap.parse_args()
a.outdir.mkdir(parents=True, exist_ok=True)
MEM = [r["accession"] for r in csv.DictReader(open(a.members, newline=""), delimiter="\t")]

def rows_of(p):
    """MacSyFinder tables carry '#' comment banners before the real header."""
    if not p.exists():
        return []
    lines = [l for l in open(p) if not l.startswith("#") and l.strip()]
    if not lines:
        return []
    hdr = lines[0].rstrip("\n").split("\t")
    return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in lines[1:]]

# ------------------------------------------------------------- antiSMASH
reg_cols = ["genome", "region", "products", "start", "end", "length_bp",
            "best_known_hit", "n_proteins_hit", "cumulative_blast"]
regions = []
for g in MEM:
    js = sorted((a.screens / f"antismash_{g}").glob("*.json"))
    if not js:
        print(f"  WARNING: no antiSMASH json for {g}")
        continue
    doc = json.load(open(js[0]))
    for rec in doc.get("records", []):
        for i, feat in enumerate(rec.get("areas", []), 1):
            start, end = int(feat["start"]), int(feat["end"])
            prods = ";".join(feat.get("products", []))
            best, nprot, cum = "", "", ""
            kc = (rec.get("modules", {})
                     .get("antismash.modules.clusterblast", {})
                     .get("knowncluster", {}).get("results", []))
            for kr in kc:
                if kr.get("region_number") == i and kr.get("ranking"):
                    top = kr["ranking"][0]
                    best = f'{top[0].get("accession","")} {top[0].get("description","")}'.strip()
                    nprot = top[1].get("hits", "")
                    cum = top[1].get("blast_score", "")
                    break
            regions.append(dict(genome=g, region=i, products=prods, start=start,
                                end=end, length_bp=end - start,
                                best_known_hit=best, n_proteins_hit=nprot,
                                cumulative_blast=cum))
with open(a.outdir / "screen_antismash_regions.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=reg_cols, delimiter="\t")
    w.writeheader(); w.writerows(regions)
print(f"screen_antismash_regions.tsv : {len(regions)} regions")

# ------------------------------------------------------- secretion systems
sys_cols = ["genome", "model", "system_id", "n_genes", "wholeness", "score",
            "first_pos", "last_pos", "genes", "loci"]
systems = []
for g in MEM:
    by = defaultdict(list)
    for r in rows_of(a.screens / f"macsy_{g}" / "best_solution.tsv"):
        by[r["sys_id"]].append(r)
    for sid, hs in by.items():
        hs.sort(key=lambda h: int(h["hit_pos"]))
        systems.append(dict(genome=g, model=hs[0]["model_fqn"].rsplit("/", 1)[-1],
                            system_id=sid, n_genes=len(hs),
                            wholeness=hs[0]["sys_wholeness"], score=hs[0]["sys_score"],
                            first_pos=hs[0]["hit_pos"], last_pos=hs[-1]["hit_pos"],
                            genes=",".join(h["gene_name"] for h in hs),
                            loci=",".join(h["hit_id"] for h in hs)))
systems.sort(key=lambda r: (r["genome"], r["model"], r["system_id"]))
with open(a.outdir / "screen_secretion_systems.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=sys_cols, delimiter="\t")
    w.writeheader(); w.writerows(systems)
print(f"screen_secretion_systems.tsv : {len(systems)} systems")

# ------------------------------------------------------ rejected candidates
rej_cols = ["genome", "model", "candidate", "n_genes", "genes", "loci", "reason"]
rej = []
for g in MEM:
    by = defaultdict(list)
    for r in rows_of(a.screens / f"macsy_{g}" / "rejected_candidates.tsv"):
        by[(r.get("candidate_id", ""), r.get("model_fqn", ""))].append(r)
    for (cid, mod), hs in by.items():
        rej.append(dict(genome=g, model=mod.rsplit("/", 1)[-1], candidate=cid,
                        n_genes=len(hs),
                        genes=",".join(h.get("gene_name", "") for h in hs),
                        loci=",".join(h.get("hit_id", "") for h in hs),
                        reason=hs[0].get("reasons", "")))
rej.sort(key=lambda r: (r["genome"], r["model"], r["candidate"]))
with open(a.outdir / "screen_secretion_rejected.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rej_cols, delimiter="\t")
    w.writeheader(); w.writerows(rej)
print(f"screen_secretion_rejected.tsv: {len(rej)} rejected candidates")

# -------------------------------------------------------- InterProScan roll-up
# Non-predictive signature types are EXCLUDED from "has a signature": these fire
# on almost every protein and would make coverage meaningless (TOOLS.md caveat).
NON_PREDICTIVE = {"MobiDBLite", "Coils", "Phobius", "TMHMM", "SignalP"}
sig = defaultdict(lambda: {"any": False, "ipr": False, "dbs": set()})
ips = a.screens / "setB_pC3_ips.tsv"
with open(ips) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        prot, db = p[0], p[3]
        s = sig[prot]
        s["dbs"].add(db)
        if db not in NON_PREDICTIVE:
            s["any"] = True
            if len(p) > 11 and p[11].startswith("IPR"):
                s["ipr"] = True
with open(a.outdir / "screen_ips_signatures.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["protein", "genome", "locus_tag", "any_signature",
                "interpro_entry", "n_databases"])
    for prot in sorted(sig):
        g, _, lt = prot.partition("|")
        s = sig[prot]
        w.writerow([prot, g, lt, s["any"], s["ipr"], len(s["dbs"])])
n_any = sum(1 for s in sig.values() if s["any"])
n_ipr = sum(1 for s in sig.values() if s["ipr"])
print(f"screen_ips_signatures.tsv    : {len(sig)} proteins with any row; "
      f"{n_any} predictive ({100*n_any/max(len(sig),1):.1f}%), {n_ipr} with an InterPro entry")
