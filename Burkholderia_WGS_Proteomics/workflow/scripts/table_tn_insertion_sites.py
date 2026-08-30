"""Published Tn insertion site table, WGS clones plus the Sanger-only ones."""
import csv
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from supp_xlsx_style import render_xlsx  # noqa: E402
from tn_annotate import parse_gff  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

COLS = [
    ("clone", "Clone"), ("platform", "Platform"), ("chromosome", "Chromosome"),
    ("insertion_pos", "Insertion_pos"), ("tsd_len_junction", "tsd_len_junction"),
    ("tsd_len_insert", "tsd_len_insert"),
    ("cassette_orientation", "Cassette_orientation"),
    ("locus_tag", "Locus_tag"), ("gene", "Gene"), ("product", "Product"),
    ("cds_start", "Cds_start"), ("cds_end", "Cds_end"),
    ("cds_strand", "Cds_strand"), ("cds_len_bp", "Cds_len_bp"),
    ("bp_into_cds", "bp_into_cds"), ("codon", "Codon"), ("aa_pos_of", "aa_pos_of"),
]
PLATFORM_ORDER = {"illumina": 0, "nanopore": 1, "Sanger": 2}


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def chrom_number(contig):
    """chr1 -> 1.  Anything unexpected is passed through rather than guessed."""
    m = re.fullmatch(r"chr(\d+)", contig)
    return m.group(1) if m else contig


PRODUCT = {g["locus_tag"]: g["product"]
           for genes in parse_gff(snakemake.input.gff).values()  # noqa: F821
           for g in genes}


def flanking_tags(row):
    """The two CDS either side of a site with no CDS of its own.

    Read off find_insertions' `context`, which looks like
        AC1V0C_13885 ends 412 bp upstream (+); AC1V0C_13890 starts 88 bp downstream (-)
    """
    return re.findall(r"(AC1V0C_\d+)", row.get("context", ""))


rows = rd(snakemake.input.wgs) + rd(snakemake.input.sanger)  # noqa: F821

out = []
for r in rows:
    rec = {k: r.get(k, "") for k, _ in COLS}
    rec["chromosome"] = chrom_number(r["contig"])
    if not r.get("locus_tag"):
        tags = flanking_tags(r)
        rec["locus_tag"] = " - ".join(tags)
        rec["product"] = " - ".join(PRODUCT.get(t, "") for t in tags)
    out.append(rec)

out.sort(key=lambda r: (PLATFORM_ORDER.get(r["platform"], 9), r["clone"]))

with open(snakemake.output.tsv, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow([h for _, h in COLS])
    for r in out:
        w.writerow([r[k] for k, _ in COLS])

render_xlsx(snakemake.output.tsv, snakemake.output.xlsx,  # noqa: F821
            sheet='insertion info', title='Transposon insertion site')

n_cds = sum(1 for r in out if r["cds_start"])
print(f"{len(out)} clones ({n_cds} in a CDS, {len(out) - n_cds} intergenic)",
      file=sys.stderr)
for r in out:
    print(f"  {r['clone']:<8} {r['platform']:<9} chr{r['chromosome']}:"
          f"{r['insertion_pos']:<10} {r['locus_tag']}", file=sys.stderr)
