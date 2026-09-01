#!/usr/bin/env python3
"""Carve the Set B replicon inputs -- ALL contigs per replicon, not just one.

Supersedes s37_setb_inputs.py, which took a single `pc3_contig` from
setB_members.tsv. That is correct for the five closed genomes but wrong for MF7,
whose pC3 spans two contigs and whose chromosome 1 spans eight; s37 therefore
carved 73% of MF7's pC3 and made a complete replicon look truncated. Contig
membership now comes from setB_contigs.tsv (see s40_setb_contigs.py).

Per genome and replicon, from annot/<acc>/<acc>.{fna,gff3,faa,gbff}:
  <out>/<acc>_<rep>.fna    the replicon's contigs        -> antiSMASH
  <out>/<acc>_<rep>.gff3   their features                -> antiSMASH
  <out>/<acc>_<rep>.faa    proteins on those contigs     -> MacSyFinder, InterProScan, clustering
  <out>/<acc>_<rep>.gbk    GenBank subset (pC3 only)     -> antiSMASH
plus a concatenated <out>/setB_<rep>_all.faa, ids "<acc>|<locus_tag>", for the
clade-core clustering.

Bakta writes the locus tag as the FASTA id and the contig as the GFF3 seqid, so
proteins are selected by walking the GFF3, not by parsing FASTA headers.
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
ap.add_argument("--contigs", required=True, type=Path)
ap.add_argument("--annot", required=True, type=Path)
ap.add_argument("--outdir", required=True, type=Path)
ap.add_argument("--replicon", required=True, choices=["pC3", "chr1"])
ap.add_argument("--genbank", action="store_true",
                help="also emit a GenBank subset (antiSMASH needs it; clustering does not)")
a = ap.parse_args()
a.outdir.mkdir(parents=True, exist_ok=True)

want = {}
for r in csv.DictReader(open(a.contigs, newline=""), delimiter="\t"):
    if r["replicon"] == a.replicon:
        want.setdefault(r["accession"], []).append(r["contig"])

all_faa = []
print(f"{'genome':<18}{'contigs':>8}{'bp':>13}{'CDS':>7}")
for acc, ctgs in sorted(want.items()):
    d = a.annot / acc
    fna, gff, faa = d / f"{acc}.fna", d / f"{acc}.gff3", d / f"{acc}.faa"
    for f in (fna, gff, faa):
        if not f.exists():
            raise SystemExit(f"FAIL: missing {f}")

    seqs = read_fasta(fna)
    # Bakta may truncate a long SPAdes contig name; fall back to a prefix match
    resolved = []
    for c in ctgs:
        if c in seqs:
            resolved.append(c)
            continue
        cand = [k for k in seqs if k.startswith(c[:40]) or c.startswith(k)]
        if len(cand) != 1:
            raise SystemExit(f"FAIL: contig {c} of {acc} matches {len(cand)} records")
        resolved.append(cand[0])
    keepset = set(resolved)

    with open(a.outdir / f"{acc}_{a.replicon}.fna", "w") as fh:
        for c in resolved:
            fh.write(f">{c}\n{seqs[c]}\n")

    keep, tags = [], set()
    for line in open(gff):
        if line.startswith("#"):
            if line.startswith("##FASTA"):
                break
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[0] not in keepset:
            continue
        keep.append(line)
        if p[2] == "CDS":
            for kv in p[8].split(";"):
                if kv.startswith("locus_tag="):
                    tags.add(kv.split("=", 1)[1])
    (a.outdir / f"{acc}_{a.replicon}.gff3").write_text("##gff-version 3\n" + "".join(keep))

    if a.genbank:
        from Bio import SeqIO
        recs = [r for r in SeqIO.parse(str(d / f"{acc}.gbff"), "genbank")
                if r.id in keepset or r.name in keepset
                or any(c.startswith(r.id) for c in keepset)]
        if len(recs) != len(keepset):
            raise SystemExit(f"FAIL: {len(recs)} GenBank records for {len(keepset)} contigs of {acc}")
        SeqIO.write(recs, str(a.outdir / f"{acc}_{a.replicon}.gbk"), "genbank")

    prot = read_fasta(faa)
    sel = {t: prot[t] for t in sorted(tags) if t in prot}
    with open(a.outdir / f"{acc}_{a.replicon}.faa", "w") as fh:
        for t, s in sel.items():
            fh.write(f">{t}\n{s}\n")
            all_faa.append((f"{acc}|{t}", s))
    print(f"{acc:<18}{len(resolved):>8}{sum(len(seqs[c]) for c in resolved):>13,}{len(sel):>7}")

with open(a.outdir / f"setB_{a.replicon}_all.faa", "w") as fh:
    for n, s in all_faa:
        fh.write(f">{n}\n{s}\n")
print(f"\nsetB_{a.replicon}_all.faa: {len(all_faa):,} proteins from {len(want)} genomes")
