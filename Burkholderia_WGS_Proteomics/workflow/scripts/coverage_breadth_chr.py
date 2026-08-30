"""Reduce the raw coverage table to the MAPQ>=10 view the coverage table publishes."""
import csv
import sys
from decimal import Decimal

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

COLS = ["numreads_q10", "covbases_q10", "breadth_pct_q10", "meandepth_q10",
        "meanmapq_q10"]
CASSETTE = snakemake.params.cassette  # noqa: F821


def plain(v):
    """samtools prints %g; the published table carries no exponents."""
    return format(Decimal(v), "f") if "e" in v.lower() else v


with open(snakemake.input.table) as fh:  # noqa: F821
    rows = list(csv.DictReader(fh, delimiter="\t"))

by = {}
for r in rows:
    by.setdefault(r["clone"], {})[r["contig"]] = r

contigs = sorted({r["contig"] for r in rows} - {CASSETTE}) + [CASSETTE]

with open(snakemake.output.table, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["clone", "platform", "contig", "length_bp"] + COLS)
    for clone in sorted(by):
        for contig in contigs:
            r = by[clone].get(contig)
            if r is None:
                sys.exit(f"{clone}: no coverage row for {contig}")
            w.writerow([clone, r["platform"], contig, r["length_bp"]]
                       + [plain(r[c]) for c in COLS])

print(f"{len(by)} clones x {len(contigs)} contigs", file=sys.stderr)
