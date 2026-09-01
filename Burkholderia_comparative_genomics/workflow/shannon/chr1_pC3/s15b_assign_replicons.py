#!/usr/bin/env python3
"""Assign chromosome-1 and c3 contigs for the 3 newly annotated MF6 neighbours.

These genomes were dereplicated away before Stage 6, so their c3 calls in
c3_calls_all_genomes.tsv read "inherited" -- propagated from their dereplication
representative, never measured. They are near-clones of retained genomes (ANI
~96-99.7%), so their c3 replicon is identified by skani against the KNOWN c3
contigs of the retained set, rather than by re-running the c3 classifier.
Retraining the classifier on a larger set would perturb every published c3 call;
this way the frozen calls stay frozen.

  chr1 := the longest contig.
  c3   := the secondary contig (>=300 kb, not chr1) whose best skani hit against
          the known c3 set clears ANI >= 95% and aligned fraction >= 60%.

Any genome failing that bar is reported as UNASSIGNED and must be dropped from
Set A rather than forced.
"""
import subprocess
import sys
from pathlib import Path

W = Path("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
SKANI = W / "envs" / "burk" / "bin" / "skani"
WORK = W / "tmp" / "assign_replicons"
WORK.mkdir(parents=True, exist_ok=True)

NEW = ["GCF_016899425.1", "GCF_053209605.1", "GCF_000203955.1"]
# retained relatives of the three, all confirmed c3 carriers
REF_ACC = ["GCF_905400185.1", "GCF_053038975.1", "GCF_003966315.1",
           "GCF_014211915.1", "GCF_003076415.1", "GCF_000019505.1"]
MIN_SECONDARY = 300_000
MIN_ANI, MIN_AF = 95.0, 60.0


def read_fasta(path, after_fasta_marker=False):
    """-> list of (name, [lines]). If after_fasta_marker, skip to ##FASTA first."""
    seqs, name, buf, started = [], None, [], not after_fasta_marker
    with open(path) as fh:
        for line in fh:
            if not started:
                if line.startswith("##FASTA"):
                    started = True
                continue
            if line.startswith(">"):
                if name:
                    seqs.append((name, buf))
                name, buf = line[1:].split()[0], []
            elif name:
                buf.append(line)
    if name:
        seqs.append((name, buf))
    return seqs


# ---- reference c3 contigs, taken straight from the gff_c3 FASTA blocks ----
ref_files = []
for acc in REF_ACC:
    g = W / "pangenome" / "gff_c3" / f"{acc}.gff3"
    if not g.exists():
        print(f"WARN: no gff_c3 for {acc}", file=sys.stderr)
        continue
    for name, buf in read_fasta(g, after_fasta_marker=True):
        p = WORK / f"ref__{acc}__{name}.fna"
        p.write_text(f">{name}\n" + "".join(buf))
        ref_files.append(str(p))
print(f"reference c3 contigs: {len(ref_files)}")
if not ref_files:
    sys.exit("no reference c3 contigs found")

rows = []
for acc in NEW:
    contigs = read_fasta(W / "genomes" / f"{acc}.fna")
    lengths = {n: sum(len(x.strip()) for x in b) for n, b in contigs}
    order = sorted(lengths, key=lambda k: -lengths[k])
    chr1 = order[0]
    secondary = [c for c in order[1:] if lengths[c] >= MIN_SECONDARY]
    print(f"\n{acc}: {len(contigs)} contigs; chr1={chr1} ({lengths[chr1]/1e6:.2f} Mb); "
          f"secondary>=300kb={[(c, round(lengths[c]/1e6, 2)) for c in secondary]}")

    qfiles = []
    for name, buf in contigs:
        if name in secondary:
            p = WORK / f"q__{acc}__{name}.fna"
            p.write_text(f">{name}\n" + "".join(buf))
            qfiles.append((name, str(p)))

    best = {}
    if qfiles:
        outp = WORK / f"dist__{acc}.tsv"
        cmd = ([str(SKANI), "dist", "-t", "16", "-s", "80", "-o", str(outp),
                "-q"] + [f for _, f in qfiles] + ["-r"] + ref_files)
        subprocess.run(cmd, check=True, capture_output=True)
        for line in outp.read_text().splitlines()[1:]:
            f = line.split("\t")
            if len(f) < 5:
                continue
            qname = Path(f[1]).stem.split("__")[-1]
            ani, af_r, af_q = float(f[2]), float(f[3]), float(f[4])
            af = max(af_r, af_q)
            if qname not in best or ani > best[qname][0]:
                best[qname] = (ani, af, Path(f[0]).stem)

    c3 = None
    for name in secondary:
        ani, af, ref = best.get(name, (0.0, 0.0, ""))
        ok = ani >= MIN_ANI and af >= MIN_AF
        print(f"   {name} {lengths[name]/1e6:.2f} Mb -> ANI {ani:.2f} AF {af:.1f} "
              f"({ref}) {'C3' if ok else '-'}")
        if ok and (c3 is None or ani > best[c3][0]):
            c3 = name
    rows.append({"accession": acc, "chr1_contig": chr1,
                 "chr1_len": lengths[chr1],
                 "c3_contig": c3 or "UNASSIGNED",
                 "c3_len": lengths.get(c3, 0),
                 "c3_ani_to_ref": round(best.get(c3, (0,))[0], 2) if c3 else "",
                 "c3_aligned_fraction": round(best.get(c3, (0, 0))[1], 1) if c3 else "",
                 "c3_best_ref": best.get(c3, (0, 0, ""))[2] if c3 else ""})

out = W / "results" / "neighbour_replicon_assignment.tsv"
with open(out, "w") as fh:
    fh.write("\t".join(rows[0].keys()) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[k]) for k in rows[0].keys()) + "\n")
print(f"\nwrote {out}")

# ---- build the GFF3 subsets into NEW directories ----
# New dirs, so the 140/306-genome inputs used for the main runs stay byte-identical.
for setname, key in (("gff_c3_plus", "c3_contig"), ("gff_chr1_plus", "chr1_contig")):
    d = W / "ppanggolin" / setname
    d.mkdir(parents=True, exist_ok=True)
    for r in rows:
        if r[key] == "UNASSIGNED":
            print(f"SKIP {r['accession']} for {setname}: unassigned")
            continue
        subprocess.run([sys.executable, str(W / "subset_gff.py"),
                        str(W / "annot" / r["accession"] / f"{r['accession']}.gff3"),
                        str(d / f"{r['accession']}.gff3"), r[key]], check=True)
print("done")
