#!/usr/bin/env python3
"""s23_pyrodigal.py -- uniform gene calls across the whole genome set.

Reconstructed from the settings recorded in TOOLS.md for the original run:

    gf = pyrodigal.GeneFinder(meta=False)      # SINGLE mode, self-training
    gf.train(*contigs, translation_table=11)

Protein IDs are emitted as ``{accession}|{contig}|{n}`` where n is pyrodigal's
per-contig gene ordinal (1-based).  That is the format the retained
`toxin12_search_all_hits.tsv` splits into its acc / contig / gene columns.

Single mode (not meta) is deliberate: Prodigal's own guidance is that
self-training on a whole genome beats the metagenomic model whenever the input
really is one organism, and every input here is an isolate assembly.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import traceback
from pathlib import Path

import pyrodigal


def read_fasta(path: Path):
    seqs, sid, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if sid is not None:
                    seqs.append((sid, "".join(buf)))
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if sid is not None:
        seqs.append((sid, "".join(buf)))
    return seqs


def call_one(args):
    acc, fna, outdir = args
    try:
        contigs = read_fasta(Path(fna))
        if not contigs:
            return acc, 0, "empty fasta"
        gf = pyrodigal.GeneFinder(meta=False)
        # train() needs bytes-like sequence(s); short assemblies can be too
        # small to self-train, in which case fall back to the metagenomic model
        # rather than dropping the genome.
        try:
            gf.train(*[s for _, s in contigs], translation_table=11)
        except ValueError:
            gf = pyrodigal.GeneFinder(meta=True)
        n = 0
        out = Path(outdir) / f"{acc}.faa"
        # Coordinates are written in the SAME pass as the proteins so the gene
        # index in the FASTA header and the row in the coords table cannot drift
        # apart. Downstream searches report hit positions from this table.
        crd = Path(outdir) / f"{acc}.coords.tsv"
        with open(out, "w") as fh, open(crd, "w") as cf:
            cf.write("acc\tcontig\tgene\tstart\tend\tstrand\tcontig_len\n")
            for cid, seq in contigs:
                genes = gf.find_genes(seq)
                clen = len(seq)
                for i, g in enumerate(genes, start=1):
                    fh.write(f">{acc}|{cid}|{i}\n")
                    prot = g.translate()
                    if prot.endswith("*"):
                        prot = prot[:-1]
                    for j in range(0, len(prot), 60):
                        fh.write(prot[j:j + 60] + "\n")
                    cf.write(f"{acc}\t{cid}\t{i}\t{g.begin}\t{g.end}\t"
                             f"{'+' if g.strand > 0 else '-'}\t{clen}\n")
                    n += 1
        return acc, n, ""
    except Exception:
        return acc, 0, traceback.format_exc().splitlines()[-1]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--list", required=True, type=Path)
    ap.add_argument("--mf6", type=Path, default=None)
    ap.add_argument("--mf7", type=Path, default=None)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args(argv)

    a.outdir.mkdir(parents=True, exist_ok=True)
    extra = {}
    if a.mf6:
        extra["MF6"] = a.mf6
    if a.mf7:
        extra["MF7"] = a.mf7

    tasks, missing = [], []
    for line in open(a.list):
        acc = line.strip()
        if not acc:
            continue
        fna = extra.get(acc) or (a.genomes / f"{acc}.fna")
        if not Path(fna).exists():
            missing.append(acc)
            continue
        tasks.append((acc, str(fna), str(a.outdir)))
    if missing:
        print(f"WARNING: {len(missing)} missing FASTA: {missing[:5]}", file=sys.stderr)

    print(f"calling genes on {len(tasks)} genomes with {a.jobs} workers", flush=True)
    ok = fail = total = 0
    with mp.Pool(a.jobs) as pool:
        for i, (acc, n, err) in enumerate(pool.imap_unordered(call_one, tasks, chunksize=1), 1):
            if err:
                fail += 1
                print(f"  FAIL {acc}: {err}", file=sys.stderr, flush=True)
            else:
                ok += 1
                total += n
            if i % 50 == 0:
                print(f"  {i}/{len(tasks)} done", flush=True)
    print(f"done: {ok} ok, {fail} failed, {total} proteins")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
