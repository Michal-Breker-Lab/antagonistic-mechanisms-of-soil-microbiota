#!/usr/bin/env python3
"""Carve the Set B pC3 inputs the three screens need, out of the Bakta annotation.

Per Set B genome, from annot/<acc>/<acc>.{fna,gff3,faa}:
  setb_pc3/<acc>_pC3.fna   the pC3 contig only        -> antiSMASH
  setb_pc3/<acc>_pC3.gff3  its CDS features           -> antiSMASH (--genefinding-tool none)
  setb_pc3/<acc>_pC3.faa   proteins on that contig    -> MacSyFinder, InterProScan
  setb_pc3/<acc>_pC3.gbk   GenBank subset of the contig -> antiSMASH (it needs the
                           CDS features in GenBank form, not FASTA + GFF3)
and one concatenated setb_pc3/setB_pC3_all.faa for the InterProScan pass.

Bakta writes the locus tag as the FASTA id and the contig in the GFF3 seqid, so
proteins are selected by walking the GFF3 rather than by parsing FASTA headers.
"""
import argparse
import csv
from pathlib import Path


def read_fasta(p):
    name, buf, out = None, [], {}
    for line in open(p):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--members", required=True, type=Path)
ap.add_argument("--annot", required=True, type=Path)
ap.add_argument("--outdir", required=True, type=Path)
a = ap.parse_args()
a.outdir.mkdir(parents=True, exist_ok=True)

members = list(csv.DictReader(open(a.members, newline=""), delimiter="\t"))
all_faa = []
print(f"{'genome':<20}{'contig':<46}{'bp':>10}{'CDS':>7}")
for m in members:
    acc, contig = m["accession"], m["pc3_contig"]
    d = a.annot / acc
    fna, gff, faa = d / f"{acc}.fna", d / f"{acc}.gff3", d / f"{acc}.faa"
    for f in (fna, gff, faa):
        if not f.exists():
            raise SystemExit(f"FAIL: missing {f}")

    seqs = read_fasta(fna)
    if contig not in seqs:
        cand = [k for k in seqs if k.startswith(contig[:40])]
        if len(cand) != 1:
            raise SystemExit(f"FAIL: contig {contig} not in {fna.name} "
                             f"({len(seqs)} contigs, {len(cand)} prefix matches)")
        contig = cand[0]
    (a.outdir / f"{acc}_pC3.fna").write_text(f">{contig}\n{seqs[contig]}\n")

    keep, tags = [], set()
    for line in open(gff):
        if line.startswith("#"):
            if line.startswith("##FASTA"):
                break
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[0] != contig:
            continue
        keep.append(line)
        if p[2] == "CDS":
            for kv in p[8].split(";"):
                if kv.startswith("locus_tag="):
                    tags.add(kv.split("=", 1)[1])
    (a.outdir / f"{acc}_pC3.gff3").write_text("##gff-version 3\n" + "".join(keep))

    # antiSMASH takes a GenBank carrying the existing CDS features.
    gbff = d / f"{acc}.gbff"
    if gbff.exists():
        from Bio import SeqIO
        recs = [r for r in SeqIO.parse(str(gbff), "genbank") if r.id == contig
                or r.name == contig or contig.startswith(r.id)]
        if len(recs) != 1:
            raise SystemExit(f"FAIL: {len(recs)} GenBank records match {contig} in {gbff.name}")
        SeqIO.write(recs, str(a.outdir / f"{acc}_pC3.gbk"), "genbank")
    else:
        raise SystemExit(f"FAIL: missing {gbff}")

    prot = read_fasta(faa)
    sel = {t: prot[t] for t in sorted(tags) if t in prot}
    with open(a.outdir / f"{acc}_pC3.faa", "w") as fh:
        for t, s in sel.items():
            fh.write(f">{t}\n{s}\n")
            all_faa.append((f"{acc}|{t}", s))
    print(f"{acc:<20}{contig[:44]:<46}{len(seqs[contig]):>10,}{len(sel):>7}")

with open(a.outdir / "setB_pC3_all.faa", "w") as fh:
    for n, s in all_faa:
        fh.write(f">{n}\n{s}\n")
print(f"\nsetB_pC3_all.faa: {len(all_faa):,} proteins from {len(members)} genomes")
