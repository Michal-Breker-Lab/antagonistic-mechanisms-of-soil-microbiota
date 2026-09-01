#!/usr/bin/env python3
"""
s23 - uniform gene calls across ALL 771 assemblies with pyrodigal.

Why this exists: only 309 of the 771 genomes carry a Bakta proteome (the
dereplicated set used for the trees). Searching only those would answer the
question for 40% of the data and would additionally confound "absent" with
"not annotated". Pyrodigal re-calls genes on every assembly with one gene
caller and one parameter set, so presence/absence is not an artefact of which
annotation pipeline happened to touch a genome.

Mode: SINGLE (self-training), not meta. Prodigal's own guidance is that single
mode is more accurate whenever >=20 kb of sequence is available to train on;
these are 7-9 Mb bacterial genomes. train() merges contigs with TTAATTAATTAA
linkers internally, so passing all contigs is correct. Falls back to meta mode
only if training raises (recorded in the log, not silently).

Headers: >{accession}|{contig}|{gene_index}|{start}|{end}|{strand}
so every hit can be mapped back to a replicon by coordinate.
"""
import glob
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import pyrodigal

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
GEN = f"{W}/genomes"
OUT = f"{W}/pyrodigal_prot"
os.makedirs(OUT, exist_ok=True)


def read_fasta(path):
    recs = []
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    recs.append((name, "".join(seq)))
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    return recs


def call_one(fna):
    acc = os.path.basename(fna).replace(".fna", "")
    dest = f"{OUT}/{acc}.faa"
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return (acc, "cached", 0, "")
    try:
        recs = read_fasta(fna)
        seqs = [s for _, s in recs]
        total = sum(len(s) for s in seqs)
        mode = "single"
        try:
            if total < 20000:
                raise ValueError(f"only {total} bp, below single-mode minimum")
            gf = pyrodigal.GeneFinder(meta=False)
            gf.train(*seqs, translation_table=11)
        except ValueError as e:
            mode = f"meta({e})"
            gf = pyrodigal.GeneFinder(meta=True)

        n = 0
        with open(dest + ".tmp", "w") as out:
            for contig, seq in recs:
                for i, g in enumerate(gf.find_genes(seq), 1):
                    n += 1
                    out.write(f">{acc}|{contig}|{i}|{g.begin}|{g.end}|"
                              f"{'+' if g.strand > 0 else '-'}\n")
                    prot = g.translate()
                    for j in range(0, len(prot), 60):
                        out.write(prot[j:j + 60] + "\n")
        os.replace(dest + ".tmp", dest)
        return (acc, mode, n, "")
    except Exception:
        return (acc, "FAILED", 0, traceback.format_exc().splitlines()[-1])


def main():
    fnas = sorted(glob.glob(f"{GEN}/*.fna"))
    print(f"assemblies: {len(fnas)}", flush=True)
    ok = fail = 0
    meta_used = []
    with ProcessPoolExecutor(max_workers=48) as ex:
        futs = {ex.submit(call_one, f): f for f in fnas}
        for k, fut in enumerate(as_completed(futs), 1):
            acc, mode, n, err = fut.result()
            if mode == "FAILED":
                fail += 1
                print(f"  FAIL {acc}: {err}", flush=True)
            else:
                ok += 1
                if mode.startswith("meta"):
                    meta_used.append((acc, mode))
            if k % 100 == 0:
                print(f"  {k}/{len(fnas)} done", flush=True)
    print(f"\nok={ok} failed={fail}", flush=True)
    print(f"meta-mode fallbacks: {len(meta_used)}", flush=True)
    for a, m in meta_used:
        print(f"   {a} {m}", flush=True)
    if fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
