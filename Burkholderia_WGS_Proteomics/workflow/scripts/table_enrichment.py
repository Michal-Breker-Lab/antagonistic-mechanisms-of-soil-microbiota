#!/usr/bin/env python3
"""Published enrichment table: the day-2 contrast block plus the core-cluster block."""
import csv
import re
import sys
from pathlib import Path

import openpyxl
from supp_xlsx_style import style_body, write_header

sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821

S4 = Path(snakemake.input.s4)  # noqa: F821
ENRICH = Path(snakemake.params.enrich_dir)  # noqa: F821
OUT_XLSX = Path(snakemake.output.xlsx)  # noqa: F821
OUT_TSV = Path(snakemake.output.tsv)  # noqa: F821
VARIANT = snakemake.params.variant  # noqa: F821

SHEET = "enrichment"
TITLE = ("COG and KEGG enrichment for the day-2 MF6 co-culture proteomic response "
         "and the eight core expression clusters, Related to Figure 2")
DE_ANALYSIS = "day-2 co-culture DE"
CONTRAST_RENAME = {"coculture_d2": "MF6_C+_D2_vs_MF6_C-_D2"}
UP_GROUP = CONTRAST_RENAME["coculture_d2"].split("_vs_")[0]
DIR_LABEL = {"up": f"Up in {UP_GROUP}", "down": f"Down in {UP_GROUP}"}
COLS = ["analysis", "set", "direction", "ontology", "term", "term_name", "cog_group",
        "k", "n", "K", "N", "fold_enrichment", "p_value", "p_adj", "proteins"]


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def fix(v):
    """Punctuation commas become semicolons; a comma with no space after it belongs
    to a chemical name and is left alone.  Same rule build_supp_tables.py applies, so
    the two tables spell a term the same way."""
    return re.sub(r",\s+", "; ", v) if isinstance(v, str) and "," in v else v


s4 = rd(S4)
cluster = [r for r in s4 if r["analysis"] != DE_ANALYSIS]
cog_group = {r["term"]: r["cog_group"] for r in cluster
             if r["ontology"] == "COG" and r["cog_group"]}

de = []
for ont in ("COG", "KEGG"):
    src = ENRICH / f"{ont}_enrichment_{VARIANT.replace('de_plus_', '')}.tsv"
    for r in rd(src):
        if r["variant"] != VARIANT:
            continue
        de.append({
            "analysis": DE_ANALYSIS,
            "set": CONTRAST_RENAME.get(r["contrast"], r["contrast"]),
            "direction": DIR_LABEL.get(r["direction"], r["direction"]),
            "ontology": r["ontology"],
            "term": re.sub(r"^map", "ko", r["term"]) if ont == "KEGG" else r["term"],
            "term_name": fix(r["term_name"]),
            "cog_group": cog_group.get(r["term"], "") if ont == "COG" else "",
            "k": r["k_fg"], "n": r["n_fg"], "K": r["K_bg"], "N": r["N_bg"],
            "fold_enrichment": r["fold_enrichment"],
            "p_value": r["p_value"], "p_adj": r["p_adj"],
            "proteins": r["proteins"].replace(",", "|"),
        })
de.sort(key=lambda r: (r["ontology"], r["direction"], float(r["p_value"])))
rows = de + cluster

with OUT_TSV.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(COLS)
    for r in rows:
        w.writerow([r.get(c, "") for c in COLS])

wb = openpyxl.Workbook()
ws = wb.active
ws.title = SHEET
write_header(ws, TITLE, COLS)


def cast(v):
    if v in ("", None):
        return None
    try:
        f = float(v)
    except (ValueError, TypeError):
        return v
    return int(f) if f.is_integer() and "." not in str(v) else f


for r in rows:
    ws.append([cast(r.get(c, "")) for c in COLS])
style_body(ws)
ws.freeze_panes = "A4"
OUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
wb.save(OUT_XLSX)

print(f"{len(de)} DE rows (variant {VARIANT}) + {len(cluster)} cluster rows "
      f"= {len(rows)}", file=sys.stderr)
