"""MF6's TXSScan hits, one row per protein per system, joined to the annotation for
coordinates and product.
"""
import csv
import re
import sys

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

COLS = ["MF6_ID", "locus", "component", "hit_status", "system_id", "system",
        "replicon", "i_evalue", "hmm_score", "PGAP_product"]

ann = {}
for line in open(snakemake.input.gff):  # noqa: F821
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS":
        continue
    m = re.search(r"locus_tag=([^;]+)", f[8])
    if not m or m.group(1) in ann:
        continue
    p = re.search(r"product=([^;]*)", f[8])
    ann[m.group(1)] = (f[0], f[3], f[4], f[6],
                       (p.group(1) if p else "").replace("%2C", ","))
print(f"annotation: {len(ann)} CDS", file=sys.stderr)

with open(snakemake.input.best) as fh:  # noqa: F821
    rows = list(csv.DictReader(
        (l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t"))
print(f"macsyfinder hits: {len(rows)}", file=sys.stderr)

out, missing = [], []
for r in rows:
    tag = r["hit_id"]
    a = ann.get(tag)
    if a is None:
        missing.append(tag)
    contig, start, end, strand, product = a if a else ("", "", "", "", "")
    out.append({
        "MF6_ID": tag,
        "locus": f"{contig}:{start}-{end}({strand})" if a else "",
        "component": r["gene_name"],
        "hit_status": r["hit_status"],
        "system_id": r["sys_id"],
        "system": r["model_fqn"].rsplit("/", 1)[-1],
        "replicon": r["replicon"],
        "i_evalue": r["hit_i_eval"],
        "hmm_score": r["hit_score"],
        "PGAP_product": product,
    })
if missing:
    sys.exit(f"{len(missing)} hit(s) absent from the annotation: "
             f"{', '.join(sorted(set(missing))[:5])}")

start_of = {}
for r in out:
    m = re.match(r"[^:]+:(\d+)-", r["locus"])
    start_of[r["MF6_ID"]] = int(m.group(1)) if m else 0
sys_start = {}
for r in out:
    s0 = start_of[r["MF6_ID"]]
    sid = r["system_id"]
    sys_start[sid] = min(sys_start.get(sid, s0), s0)
out.sort(key=lambda r: (r["replicon"], sys_start[r["system_id"]], r["system_id"],
                        start_of[r["MF6_ID"]]))

with open(snakemake.output.genes, "w", newline="") as fh:  # noqa: F821
    w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
print(f"\n{len(out)} row(s) over {len({r['system_id'] for r in out})} system(s)",
      file=sys.stderr)
