#!/usr/bin/env python3
"""Proposed plant-compartment vocabulary for the plant/rhizosphere genomes.

Vocabulary approved 2026-08-05 (root+rhizosphere merged). Feeds the compartment
ring on Figure 4; it is DESCRIPTIVE ONLY and is not a term in any model.

Scope: the 106 genomes whose top-level host_category is `plant` or
`rhizosphere`. Compartment is parsed from the free-text BioSample fields
`isolation_source` and `host`, in that order of preference.

Rule order matters and is deliberate:
  nodule      before root_rhizosphere  -- "root nodule" is a nodule
  seed/grain  before plant             -- reproductive/processed tissue is not "plant"

`root_rhizosphere` merges rhizosphere and root on Moshe's instruction: free-text
records mix the two ("rhizosphere soil of X root") and neither stands alone at
n=18 / n=10.

`plant_unspecified` means the record names a host plant but no compartment
(e.g. host="Oryza sativa" with an empty isolation_source). It is an honest
"unknown", NOT a compartment, and must never be treated as one in a model.
"""
import collections
import csv
import re
import sys
from pathlib import Path

D = Path(__file__).resolve().parent.parent
TAB = D / "tables"

# (label, regex) -- evaluated in order, first match wins
RULES = [
    ("nodule",           r"nodule"),
    # rhizosphere and root are ONE category: the two are not cleanly separable
    # in free-text BioSample records ("rhizosphere soil of X root"), and split
    # they were n=18 and n=10, too small to stand alone.
    ("root_rhizosphere", r"rhizosph|root exudate|\broot\b|\broots\b"),
    ("seed",        r"\bseed"),
    ("grain_fruit", r"\bgrain|\bpanicle|pericarpium|broken rice|meal|food"),
    ("leaf",        r"\bleaf\b|\bleaves\b|phyllosph"),
    ("stem_shoot",  r"\bshoot|\bstem\b|\bseedling"),
    ("flower",      r"\bflower"),
]

BCC = {"cepacia", "multivorans", "cenocepacia", "stabilis", "vietnamiensis",
       "dolosa", "ambifaria", "anthina", "pyrrocinia", "ubonensis", "latens",
       "diffusa", "arboris", "seminalis", "metallica", "contaminans", "lata",
       "pseudomultivorans", "stagnalis", "territorii", "puraquae", "orbicola",
       "aenigmatica", "sola", "semiarida", "catudaia", "alpina"}


def clade(organism):
    p = organism.split()
    g, s = (p[0] if p else ""), (p[1] if len(p) > 1 else "")
    if g != "Burkholderia":
        return "other genus"
    return "Bcc" if s in BCC else "other Burkholderia"


def compartment(row):
    """-> (label, matched_text). isolation_source wins over host."""
    for field in ("raw_isolation_source", "raw_host"):
        t = row.get(field, "").strip().lower()
        if not t:
            continue
        for label, pat in RULES:
            if re.search(pat, t):
                return label, f"{field}={t}"
    joined = (row.get("raw_isolation_source", "") + row.get("raw_host", "")).strip()
    return ("plant_unspecified" if joined else "no_metadata"), ""


def main():
    hosts = list(csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"))
    calls = {r["accession"]: r["c3_present"] for r in
             csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t")}
    pl = [r for r in hosts if r["host_category"] in ("plant", "rhizosphere")]
    print(f"plant + rhizosphere genomes: {len(pl)}")

    out = TAB / "plant_compartment_for_review.tsv"
    n_c3, n_clade = collections.defaultdict(collections.Counter), \
        collections.defaultdict(collections.Counter)
    with open(out, "w") as fh:
        fh.write("accession\torganism_name\tclade\thost_category\tcompartment\t"
                 "matched_evidence\traw_isolation_source\traw_host\tc3_present\n")
        for r in sorted(pl, key=lambda x: x["accession"]):
            comp, ev = compartment(r)
            cl = clade(r["organism_name"])
            c3 = calls.get(r["accession"], "")
            n_c3[comp][c3] += 1
            n_clade[comp][cl] += 1
            fh.write(f"{r['accession']}\t{r['organism_name']}\t{cl}\t"
                     f"{r['host_category']}\t{comp}\t{ev}\t"
                     f"{r['raw_isolation_source']}\t{r['raw_host']}\t{c3}\n")
    print(f"wrote {out}")

    print(f"\n{'compartment':<20}{'n':>5}{'c3+':>6}{'c3-':>6}{'%c3':>6}"
          f"{'Bcc':>6}{'othBurk':>9}{'othGen':>8}")
    for comp in sorted(n_c3, key=lambda k: -sum(n_c3[k].values())):
        c, k = n_c3[comp], n_clade[comp]
        tot = c["True"] + c["False"]
        pct = f"{100*c['True']/tot:.0f}%" if tot else "-"
        print(f"{comp:<20}{tot:>5}{c['True']:>6}{c['False']:>6}{pct:>6}"
              f"{k['Bcc']:>6}{k['other Burkholderia']:>9}{k['other genus']:>8}")
    print("\nNOTE: compartment tracks clade closely -- see the last three columns.")


if __name__ == "__main__":
    sys.exit(main())
