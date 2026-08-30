"""Published coverage-breadth table."""
import csv
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from supp_xlsx_style import render_xlsx  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

COLS = [
    ("clone", "Clone"), ("platform", "Platform"), ("contig", "Contig"),
    ("length_bp", "Length_bp"), ("numreads_q10", "Numreads_q10"),
    ("covbases_q10", "Covbases_q10"), ("breadth_pct_q10", "Breadth_pct_q10"),
    ("meandepth_q10", "Meandepth_q10"), ("meanmapq_q10", "Meanmapq_q10"),
]

with open(snakemake.input.table) as fh:  # noqa: F821
    rows = list(csv.DictReader(fh, delimiter="\t"))

out = [[r[k] for k, _ in COLS] for r in rows]

with open(snakemake.output.tsv, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([h for _, h in COLS])
    w.writerows(out)

render_xlsx(snakemake.output.tsv, snakemake.output.xlsx,  # noqa: F821
            sheet='coverage breadth', title='Whole-genome sequencing coverage breadth')

clones = {r["clone"] for r in rows}
print(f"{len(rows)} rows: {len(clones)} clones x "
      f"{len(rows) // max(len(clones), 1)} contigs", file=sys.stderr)
