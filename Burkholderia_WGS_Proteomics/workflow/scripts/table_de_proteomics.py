"""Published DE table: every contrast, with ON/OFF calls and antiSMASH annotation."""
import csv
import sys
from pathlib import Path

import openpyxl
from supp_xlsx_style import style_body, write_header

sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821

S2 = Path(snakemake.input.s2)  # noqa: F821
OUT_XLSX = Path(snakemake.output.xlsx)  # noqa: F821
OUT_TSV = Path(snakemake.output.tsv)  # noqa: F821
ON_GROUP_MIN_PCT = float(snakemake.params.gate)  # noqa: F821
LOCUS_PREFIX = snakemake.params.locus_prefix  # noqa: F821

SHEET = "DE and ON-OFF"
TITLE = ("Differential protein abundance and presence/absence calls across seven "
         "day-2 Bso MF6 proteomic contrasts, Related to Figure 2 and Figure S3")
OLD_NAME, NEW_NAME, PCT_NAME = "on_meanLFQ", "ON_group", "ON_percentile"
PCT_CONTRAST = "MF6_C+_D2_vs_MF6_C-_D2"
PCT_SRC = Path(snakemake.input.onoff)  # noqa: F821
GENES = Path(snakemake.input.genes)  # noqa: F821

ROLE = {"biosynthetic": "core biosynthetic",
        "biosynthetic-additional": "additional biosynthetic",
        "transport": "transport-related",
        "regulatory": "regulatory",
        "resistance": "resistance",
        "other": "other"}
NEW_COLS = ["antiSMASH region", "BGC class", "antiSMASH role",
            "MIBiG accession", "MIBiG similarity"]


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))




def region_label(region):
    """antismash.tsv's "chr1.region001" -> the published "Chr1.1".

    The gene table and the antiSMASH JSON name the same region two different ways -
    "chr1.region001" against record id "chr1" plus region key "1" - so both are
    normalised to one display label here.  Getting this wrong does not raise: the
    MIBiG lookup simply misses and the accession/similarity columns come back empty,
    which is exactly how the first rebuild lost 1,368 of them."""
    contig, _, num = region.partition(".region")
    if not num:
        return region
    return f"{contig.capitalize()}.{int(num)}"


def protein_annotation():
    """Region, BGC class, role and the MIBiG match, per protein."""
    ann, unmapped = {}, set()
    for r in rd(GENES):
        lab = region_label(r.get("region_label") or r.get("region") or "")
        kind = (r.get("gene_kind") or "").strip()
        if kind and kind not in ROLE:
            unmapped.add(kind)
        acc = r.get("mibig_accession", "")
        sim = r.get("mibig_similarity", "")
        ann[r["protein_id"]] = (lab, r.get("region_product", ""),
                                ROLE.get(kind, ""), acc, sim)
    if unmapped:
        sys.exit(f"unrecognised gene_kind value(s): {sorted(unmapped)}")
    return ann


rows = rd(S2)
if not rows:
    sys.exit(f"no rows in {S2}")
cols = [c for c in rows[0]]
if OLD_NAME not in cols:
    sys.exit(f"{OLD_NAME!r} not in {S2}: {cols}")

pct = {r["Protein"]: float(r["pct_in_opposite_group"]) for r in rd(PCT_SRC)}
kept = cleared = filled = missing = 0
for r in rows:
    on = r.get(OLD_NAME) or ""
    r[PCT_NAME] = ""
    if not on:
        r[NEW_NAME] = ""
        continue
    va, vb = int(r["n_detected_A"] or 0), int(r["n_detected_B"] or 0)
    if va > 0 and vb == 0:
        group = r["group_A"]
    elif vb > 0 and va == 0:
        group = r["group_B"]
    else:
        sys.exit(f"{r['Protein']} / {r['contrast']}: cannot resolve ON group "
                 f"({va}, {vb})")
    if r["contrast"] == PCT_CONTRAST:
        if r["Protein"] in pct:
            r[PCT_NAME] = pct[r["Protein"]]
            filled += 1
        else:
            missing += 1
    q = r[PCT_NAME]
    if q != "" and float(q) > ON_GROUP_MIN_PCT:
        r[NEW_NAME] = group
        kept += 1
    else:
        r[NEW_NAME] = ""
        cleared += 1

ann = protein_annotation()
hit = 0
for r in rows:
    a = ann.get(r["Protein"]) if r["Protein"].startswith(LOCUS_PREFIX) else None
    if a:
        hit += 1
    for name, v in zip(NEW_COLS, a or ("", "", "", "", "")):
        r[name] = v if v not in (None,) else ""

i = cols.index(OLD_NAME)
out_cols = cols[:i] + [NEW_NAME, PCT_NAME] + cols[i + 1:] + NEW_COLS

with open(OUT_TSV, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(out_cols)
    for r in rows:
        w.writerow(["" if r.get(c) in (None,) else r.get(c, "") for c in out_cols])

wb = openpyxl.Workbook()
ws = wb.active
ws.title = SHEET
write_header(ws, TITLE, out_cols)
for r in rows:
    ws.append([r.get(c, "") if r.get(c, "") != "" else None for c in out_cols])
style_body(ws)
ws.freeze_panes = "A4"
OUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
wb.save(OUT_XLSX)

print(f"{len(rows)} rows from {S2.name}; {NEW_NAME} kept on {kept} row(s) with "
      f"{PCT_NAME} > {ON_GROUP_MIN_PCT:g}, cleared on {cleared}; {PCT_NAME} filled "
      f"{filled}" + (f", {missing} NOT FOUND" if missing else "")
      + f"; antiSMASH on {hit} rows -> {len(out_cols)} columns", file=sys.stderr)
