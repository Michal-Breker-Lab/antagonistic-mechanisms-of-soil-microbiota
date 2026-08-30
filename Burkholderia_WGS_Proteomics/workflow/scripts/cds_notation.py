"""Add HGVS-style CDS notation and nearest-gene columns to the insertion table."""
import csv
import sys
from collections import defaultdict

sys.stderr = open(snakemake.log[0], "w")

CASSETTE_LEN = snakemake.params.cassette_len
END = "nearest_gene_end"
NEW = ["cds_notation", "notation_c", "notation_g", "nearest_gene",
       "nearest_gene_symbol", "nearest_gene_distance_bp",
       "nearest_gene_relation", "nearest_gene_product"]


def read_gff(path):
    feats = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if "locus_tag" in a:
                feats[f[0]].append({
                    "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                    "locus_tag": a["locus_tag"], "gene": a.get("gene", ""),
                    "product": a.get("product", "").replace("%2C", ","),
                })
    for c in feats:
        feats[c].sort(key=lambda g: g["start"])
    return feats


def nearest(feats, contig, pos):
    """Closest CDS, with the relation named in the gene's own orientation."""
    best = None
    for g in feats.get(contig, []):
        if g["start"] <= pos <= g["end"]:
            return g, 0, "inside"
        if pos > g["end"]:
            d, side = pos - g["end"], ("downstream" if g["strand"] == "+"
                                       else "upstream")
        else:
            d, side = g["start"] - pos, ("upstream" if g["strand"] == "+"
                                         else "downstream")
        if best is None or d < best[1]:
            best = (g, d, side)
    return best if best else (None, "", "")


feats = read_gff(snakemake.input.gff)
with open(snakemake.input.sites) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

for r in rows:
    for c in [END] + NEW:
        r[c] = ""
    r["notation_g"] = f"g.{int(r['insertion_pos'])}"
    if r["locus_tag"] and r["bp_into_cds"]:
        nt = int(r["bp_into_cds"])
        r["cds_notation"] = (f"{r['locus_tag']}:c.{nt}_{nt + 1}"
                             f"insTn5[{CASSETTE_LEN}]")
        r["notation_c"] = f"{r['gene'] or r['locus_tag']}:c.{nt}"
        continue
    g, d, side = nearest(feats, r["contig"], int(r["insertion_pos"]))
    if g:
        r["nearest_gene"] = g["locus_tag"]
        r["nearest_gene_symbol"] = g["gene"]
        r["nearest_gene_distance_bp"] = d
        r["nearest_gene_relation"] = side
        r["nearest_gene_product"] = g["product"]
        r[END] = "5'" if side == "upstream" else "3'"
        off = f"-{d}" if side == "upstream" else f"*{d}"
        r["notation_c"] = f"{g['gene'] or g['locus_tag']}:c.{off}"

cols = [c for c in rows[0] if c not in NEW] + NEW
with open(snakemake.output.sites, "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

n_cds = sum(1 for r in rows if r["cds_notation"])
print(f"{n_cds} in a CDS, {len(rows) - n_cds} intergenic", file=sys.stderr)
