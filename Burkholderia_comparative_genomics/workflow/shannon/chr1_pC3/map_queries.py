#!/usr/bin/env python3
"""Map the 12 pasted toxin queries onto MF6's Bakta annotation by exact sequence
match (whole protein) or substring (fragment of a larger protein)."""
import csv
from pathlib import Path

W = Path("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")


def read_fasta(p):
    d, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                d[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        d[name] = "".join(buf)
    return d


q = read_fasta(W / "queries_toxins12.faa")
mf6 = read_fasta(W / "annot/MF6/MF6.faa")

# Bakta TSV: locus tag -> contig, start, stop, strand, gene, product
info = {}
with open(W / "annot/MF6/MF6.tsv") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        # cols: Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product, DbXrefs
        info[f[5]] = dict(contig=f[0], start=f[2], stop=f[3], strand=f[4],
                          gene=f[6], product=f[7],
                          dbxref=f[8] if len(f) > 8 else "")

rows = []
for qid, qs in q.items():
    hit, kind, off = None, None, ""
    for lt, ps in mf6.items():
        if ps == qs:
            hit, kind = lt, "exact"
            break
    if hit is None:
        for lt, ps in mf6.items():
            i = ps.find(qs)
            if i >= 0:
                hit, kind = lt, "substring"
                off = f"{i+1}-{i+len(qs)} of {len(ps)}"
                break
    d = info.get(hit, {})
    rows.append(dict(query=qid, qlen=len(qs), locus_tag=hit or "NO_MATCH",
                     match=kind or "none", span=off, contig=d.get("contig", ""),
                     start=d.get("start", ""), stop=d.get("stop", ""),
                     strand=d.get("strand", ""), gene=d.get("gene", ""),
                     product=d.get("product", ""), dbxref=d.get("dbxref", "")))

out = W / "queries_toxins12_map.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

for r in rows:
    print(f"{r['query']:<26} {r['qlen']:>5}aa  {r['locus_tag']:<16} {r['match']:<10} "
          f"{r['contig']:<10} {r['start']:>8}-{r['stop']:<8} {r['strand']}  {r['product'][:60]}")
    if r["span"]:
        print(f"{'':<26}        -> spans {r['span']}")
