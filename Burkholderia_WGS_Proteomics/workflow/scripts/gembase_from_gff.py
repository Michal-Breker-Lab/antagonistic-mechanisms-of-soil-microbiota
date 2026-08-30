"""Write a genome's proteome in MacSyFinder gembase form, plus the map from those
ids back to the annotation's locus tags and coordinates.
"""
import csv
import sys

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

MAP_COLS = ["msf_id", "replicon", "gene_oid", "locus_tag", "contig", "start",
            "end", "strand", "product"]


def read_proteome(path):
    """gene_oid -> sequence.  Headers look like `>2546952507 BMF7_00001 product`."""
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.rstrip("\n"))
    return {k: "".join(v) for k, v in seqs.items()}


prot = read_proteome(snakemake.input.faa)  # noqa: F821
print(f"proteome: {len(prot)} sequences", file=sys.stderr)

cds = []
for line in open(snakemake.input.gff):  # noqa: F821
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS":
        continue
    a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    cds.append({
        "contig": f[0], "start": int(f[3]), "end": f[4], "strand": f[6],
        "gene_oid": a.get("ID", ""), "locus_tag": a.get("locus_tag", ""),
        "product": a.get("product", "").replace("%2C", ","),
    })
print(f"annotation: {len(cds)} CDS", file=sys.stderr)

seen = []
for c in cds:
    if c["contig"] not in seen:
        seen.append(c["contig"])
rank = {c: i for i, c in enumerate(seen)}
cds.sort(key=lambda c: (rank[c["contig"]], c["start"]))

names, counts = {}, {}
rows, missing = [], []
for c in cds:
    if c["contig"] not in names:
        names[c["contig"]] = f"ctg{len(names) + 1:04d}"
    rep = names[c["contig"]]
    counts[rep] = counts.get(rep, 0) + 1
    c["replicon"] = rep
    c["msf_id"] = f"{rep}_g{counts[rep]:05d}"
    if c["gene_oid"] not in prot:
        missing.append(c["gene_oid"])
    rows.append(c)

if missing:
    sys.exit(f"{len(missing)} CDS with no protein sequence, e.g. "
             f"{', '.join(missing[:5])}")

with open(snakemake.output.faa, "w") as fh:  # noqa: F821
    for c in rows:
        fh.write(f">{c['msf_id']}\n{prot[c['gene_oid']]}\n")

with open(snakemake.output.map, "w", newline="") as fh:  # noqa: F821
    w = csv.DictWriter(fh, MAP_COLS, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

print(f"\n{len(rows)} protein(s) over {len(names)} replicon(s)", file=sys.stderr)
