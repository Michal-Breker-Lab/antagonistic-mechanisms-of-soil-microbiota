#!/usr/bin/env python3
"""Stage 10 - map BioSample free text to a controlled host/habitat vocabulary.

Design notes:
  * Evidence is combined in PRIORITY order, most specific first. host_disease is
    the strongest clinical signal; bs_host names an actual organism; only then
    isolation_source, which is the messiest field.
  * NULL_TOKENS matter more than they look: "not collected", "not applicable",
    "missing", "unknown", "na", "none" are all populated-but-empty and must not
    be mistaken for data. 100+ genomes carry one of these.
  * attr_sample_type is deliberately IGNORED. Its top value is "cell culture"
    (120 genomes) which describes sample preparation, not provenance.
  * Every rule is recorded per genome in `evidence` so the assignment is
    auditable, and a full raw-string -> category table is emitted for review.
"""
import collections
import csv
import os
import re
import sys

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES = f"{W}/results"

NULL_TOKENS = {"", "na", "n/a", "none", "null", "unknown", "missing",
               "not collected", "not applicable", "not available",
               "not determined", "not provided", "unspecified", "-"}

# (category, regex) - order within a field matters; first match wins.
HOST_RULES = [
    ("human_clinical", r"\bhomo sapiens\b|\bhuman\b"),
    ("fungus",         r"rhizopus|fungus|fungal|mycelium|hypha|aspergillus|mucor"),
    ("amoeba_protist", r"dictyostelium|acanthamoeba|amoeba|paramecium"),
    ("insect",         r"physopelta|insect|drosophila|beetle|aphid|termite|larva|"
                       r"anasa tristis|leptoglossus|lagria|plectrocnemia|"
                       r"hemiptera|coleoptera|bug\b"),
    ("animal",         r"equus|caballus|horse|goat|sheep|ovis|capra|iguana|bovine|"
                       r"cattle|\bcow\b|\bpig\b|\bsus scrofa\b|dog|canis|feline|"
                       r"felis|bird|chicken|gallus|donkey|camel|alpaca|mouse|mus "
                       r"musculus|\brat\b|rattus|monkey|primate|macaca|dolphin|"
                       r"reptile|parrot"),
    ("plant",          r"oryza|\brice\b|sugarcane|saccharum|arabidopsis|allium|"
                       r"onion|nepenthes|fir\b|cunninghamia|zea mays|\bmaize\b|"
                       r"\bcorn\b|glycine max|soybean|triticum|wheat|hordeum|"
                       r"solanum|tomato|potato|nicotiana|tobacco|vitis|grape|"
                       r"musa|banana|coffea|coffee|gladiolus|iris|orchid|"
                       r"phaseolus|bean|medicago|lotus|mimosa|populus|eucalyptus|"
                       r"pine|pinus|bamboo|plant|tree|moss|lichen|fern|"
                       r"pelargonium|freesia|gladiolus|dendrobium|cinnamomum|"
                       r"poplar|citri|camphora|orchid"),
    # a host field naming a habitat rather than an organism
    ("soil",           r"^soil$|^sand$"),
    ("water",          r"^water$"),
    ("environmental_other", r"^environment$"),
]

SOURCE_RULES = [
    ("human_clinical", r"^human$|human |\bblood\b|sputum|\burine\b|"
                       r"cystic fibrosis|\bcf\b|\bliver\b|brain tissue|"
                       r"wet tissue|"
                       r"clinical|patient|catheter|wound|abscess|nasal|"
                       r"bronch|tracheal|respiratory|melioidosis|sepsis|"
                       r"cepacia syndrome|swab|pus|cerebrospinal|csf|"
                       r"bloodstream|blood stream|throat|lung|"
                       r"hospital|infection|surgical|milk|mastitis"),
    ("rhizosphere",    r"rhizosphere|rhizoplane|root nodule|nodule|rhizospheric"),
    ("plant",          r"soybean|gladiolus|dendrobium|pericarpium|"
                       r"\bleaf\b|leaves|\broot\b|stem|seed|grain|panicle|"
                       r"\bplant\b|rice|maize|corn|onion|sugarcane|"
                       r"phyllosphere|endophyt|fruit|flower|bark|tuber|"
                       r"vegetable|crop|paddy|orchard|diseased"),
    ("fungus",         r"fungus|fungal|hypha|mycelium|rhizopus|mushroom|"
                       r"sporangi|lichen"),
    ("insect",         r"insect|beetle|larva|gut of|termite|aphid"),
    ("animal",         r"\bhorse\b|equine|goat|sheep|bovine|cattle|swine|"
                       r"veterinar|animal|glanders"),
    ("soil",           r"\bsoil\b|sediment|rhizospheric soil|compost|mud|"
                       r"\bfield\b|farmland|\bdust\b|sand"),
    ("water",          r"\bwater\b|\briver\b|\blake\b|\bpond\b|marine|"
                       r"seawater|groundwater|aquatic|stream|well water|"
                       r"drinking|effluent|wastewater|sewage"),
    ("industrial",     r"industrial|pharmaceutic|disinfectant|detergent|"
                       r"contaminated product|manufactur|antiseptic|"
                       r"mouthwash|cosmetic|solution|saline|reagent|"
                       r"food|dairy|beverage|juice|fermentation|baijiu|"
                       r"activated sludge|chemical product|sludge"),
    ("environmental_other", r"environment|aerosol|\bair\b|atmospher|"
                            r"biofilm|bioreactor|enrich"),
]

DISEASE_CLINICAL = re.compile(
    r"melioidosis|meliodosis|cystic fibrosis|sepsis|cepacia syndrome|"
    r"blood stream|bloodstream|mastitis|abscess|wound|infection|"
    r"nasal inflammation|surgical")


def clean(s):
    s = (s or "").strip().lower()
    return "" if s in NULL_TOKENS else s


def match(rules, text):
    for cat, pat in rules:
        if re.search(pat, text):
            return cat
    return None


def main():
    rows = list(csv.DictReader(open(f"{W}/metadata/assemblies.tsv"), delimiter="\t"))
    out, audit = [], collections.defaultdict(collections.Counter)

    for r in rows:
        host = clean(r.get("bs_host"))
        src = clean(r.get("bs_isolation_source"))
        dis = clean(r.get("attr_host_disease"))
        envb = clean(r.get("attr_env_broad_scale"))
        envl = clean(r.get("attr_env_local_scale"))
        envm = clean(r.get("attr_env_medium"))

        cat, ev = None, ""
        # 1. a named disease is unambiguous clinical evidence
        if dis and DISEASE_CLINICAL.search(dis):
            cat, ev = "human_clinical", f"host_disease={dis}"
        # 2. glanders is equine unless the host says human
        if cat == "human_clinical" and "glanders" in dis and "homo" not in host:
            cat, ev = "animal", f"host_disease={dis} (glanders, equine)"
        # 3. the host organism
        if cat is None and host:
            m = match(HOST_RULES, host)
            if m:
                cat, ev = m, f"host={host}"
        # 4. isolation source
        if cat is None and src:
            m = match(SOURCE_RULES, src)
            if m:
                cat, ev = m, f"isolation_source={src}"
        # 5. structured environment fields, last resort
        if cat is None:
            for fld, val in (("env_medium", envm), ("env_local", envl),
                             ("env_broad", envb)):
                if val:
                    m = match(SOURCE_RULES, val)
                    if m:
                        cat, ev = m, f"{fld}={val}"
                        break
        if cat is None:
            cat, ev = "unknown", "no usable metadata"

        out.append({
            "accession": r["accession"],
            "organism_name": r.get("organism_name", ""),
            "host_category": cat,
            "evidence": ev,
            "raw_host": r.get("bs_host", ""),
            "raw_isolation_source": r.get("bs_isolation_source", ""),
            "raw_host_disease": r.get("attr_host_disease", ""),
            "geo_loc_name": r.get("bs_geo_loc_name", ""),
        })
        if host:
            audit["bs_host"][(host, cat)] += 1
        if src:
            audit["bs_isolation_source"][(src, cat)] += 1

    cols = list(out[0].keys())
    with open(f"{RES}/host_categories.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    # full raw-string -> category table for manual review
    with open(f"{RES}/host_mapping_for_review.tsv", "w") as fh:
        fh.write("field\traw_value\tassigned_category\tn_genomes\n")
        for fld in audit:
            for (raw, cat), n in sorted(audit[fld].items(), key=lambda x: -x[1]):
                fh.write(f"{fld}\t{raw}\t{cat}\t{n}\n")

    c = collections.Counter(r["host_category"] for r in out)
    print(f"=== HOST CATEGORY ASSIGNMENT ({len(out)} genomes) ===")
    for k, v in c.most_common():
        print(f"  {v:>4}  {100*v/len(out):>5.1f}%  {k}")
    unk = c["unknown"]
    print(f"\nunknown: {unk}/{len(out)} = {100*unk/len(out):.1f}%")

    print("\n=== assignment by evidence field ===")
    ef = collections.Counter(r["evidence"].split("=")[0] for r in out)
    for k, v in ef.most_common():
        print(f"  {v:>4}  {k}")

    # what fell through despite having text
    ft = [r for r in out if r["host_category"] == "unknown"
          and (clean(r["raw_host"]) or clean(r["raw_isolation_source"]))]
    print(f"\n=== UNMATCHED DESPITE HAVING TEXT: {len(ft)} ===")
    seen = collections.Counter()
    for r in ft:
        seen[(clean(r["raw_host"]), clean(r["raw_isolation_source"]))] += 1
    for (h, s), n in seen.most_common(25):
        print(f"  {n:>3}  host={h[:30]!r:<32} source={s[:40]!r}")
    print("STAGE10_DONE")


if __name__ == "__main__":
    sys.exit(main())
