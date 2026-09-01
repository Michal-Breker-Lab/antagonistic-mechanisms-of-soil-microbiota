#!/usr/bin/env python3
"""
s26 - call genes on MF6 with the SAME pyrodigal settings used for the 771 NCBI
assemblies in s23.

MF6 is not in search_rhs/all_pyrodigal.faa (that database is the 771 downloaded
genomes only), so without this step MF6 - the source of every query - cannot be
scored on the same footing as the genomes it is being compared against.

Settings are copied from s23 deliberately: meta=False (single/self-training),
translation_table=11, and the same header format, so a hit in MF6 and a hit in
any other genome are produced by an identical gene caller. MF6 is 7.8 Mb, far
above the 20 kb single-mode floor, so the meta fallback must not trigger; if it
does, that is an error, not a fallback.
"""
import os
import sys

import pyrodigal

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
FNA = f"{W}/annot/MF6/MF6.fna"
OUT = f"{W}/pyrodigal_mf6"
os.makedirs(OUT, exist_ok=True)
ACC = "MF6"


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    return recs


recs = read_fasta(FNA)
seqs = [s for _, s in recs]
total = sum(len(s) for s in seqs)
print(f"contigs: {len(recs)}  total: {total:,} bp", flush=True)
if total < 20000:
    sys.exit(f"ERROR: only {total} bp - s23 would have used meta mode here")

gf = pyrodigal.GeneFinder(meta=False)
gf.train(*seqs, translation_table=11)
print("mode=single", flush=True)

dest = f"{OUT}/{ACC}.faa"
n = 0
with open(dest + ".tmp", "w") as out:
    for contig, seq in recs:
        for i, g in enumerate(gf.find_genes(seq), 1):
            n += 1
            out.write(f">{ACC}|{contig}|{i}|{g.begin}|{g.end}|"
                      f"{'+' if g.strand > 0 else '-'}\n")
            prot = g.translate()
            for j in range(0, len(prot), 60):
                out.write(prot[j:j + 60] + "\n")
os.replace(dest + ".tmp", dest)
print(f"proteins: {n:,} -> {dest}", flush=True)
