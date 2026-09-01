#!/usr/bin/env python3
"""
s19 - build the inputs for the InterProScan / antiSMASH / MacSyFinder screen.

  ips/pC3_clade_pangenome.faa   1,382 family representatives (833 core + 549 accessory)
  antismash/<genome>_pC3.gbk    pC3 contigs only, carrying Bakta's annotation
  macsy/<genome>_pC3.faa        every pC3 protein, in genomic order (ordered_replicon)

Scope per Moshe 2026-08-10: the five B. sola pC3 replicons, no chromosomal baseline.
"""
import csv, os, re, sys

_SRC = open("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/s17_functional.py").read()
assert "# 9. Redundancy" in _SRC
exec(_SRC.split("# 9. Redundancy")[0])

from Bio import SeqIO

SETB = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1", "MF6"]
IPS = f"{W}/screen/ips"; AS = f"{W}/screen/antismash_in"; MS = f"{W}/screen/macsy_in"
for d in (IPS, AS, MS):
    os.makedirs(d, exist_ok=True)

# gff path per genome, straight from the list PPanGGOLiN was given
GFF = {}
for line in open(f"{W}/ppanggolin/neighbours/list_setB_c3.tsv"):
    g, p = line.rstrip("\n").split("\t")[:2]
    GFF[g] = p


def load_faa(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


SEQ = {}
for g in SETB:
    SEQ.update(load_faa(f"{W}/annot/{g}/{g}.faa"))


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


# ---------------------------------------------------------------- 1. IPS input
# Representative = MF6 where the family has one (core always does), otherwise the
# first available member. Header keeps family id and set so hits join back.
n_core = n_acc = 0
with open(f"{IPS}/pC3_clade_pangenome.faa", "w") as fh:
    for fam in sorted(c3_fams):
        cells = c3_fams[fam]
        rep = None
        for g in ["MF6"] + [x for x in SETB if x != "MF6"]:
            for lt in cells.get(g, []):
                if lt in SEQ:
                    rep = lt; break
            if rep:
                break
        if not rep:
            continue
        is_core = len(cells) == len(SETB)
        L = LAB[("c3", fam)]
        cats = "".join(sorted(L["cats"])) or "-"
        tag = "core" if is_core else "accessory"
        n_core += is_core; n_acc += not is_core
        fh.write(f">{fam}|{rep}|{tag}|COG={cats}|{L['product'].replace(chr(9),' ')}\n"
                 f"{wrap(SEQ[rep])}\n")
print(f"[1] IPS input: {n_core} core + {n_acc} accessory = {n_core + n_acc} proteins",
      flush=True)

# ---------------------------------------------------------------- 2+3. per-genome pC3
summary = []
for g in SETB:
    # contigs that make up this genome's pC3, taken from the subset GFF PPanGGOLiN used
    contigs, order = set(), []
    for line in open(GFF[g]):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS":
            m = re.search(r"locus_tag=([^;]+)", f[8])
            if m:
                contigs.add(f[0])
                order.append((f[0], int(f[3]), m.group(1)))
    order.sort()

    # MacSyFinder wants the proteome in genomic order for ordered_replicon mode
    with open(f"{MS}/{g}_pC3.faa", "w") as fh:
        n = 0
        for _c, _s, lt in order:
            if lt in SEQ:
                fh.write(f">{lt}\n{wrap(SEQ[lt])}\n"); n += 1

    # antiSMASH: subset the Bakta GenBank so its own gene calls are reused
    gb = f"{W}/annot/{g}/{g}.gbff"
    kept = 0
    if os.path.exists(gb):
        recs = [r for r in SeqIO.parse(gb, "genbank") if r.id in contigs or r.name in contigs]
        if recs:
            SeqIO.write(recs, f"{AS}/{g}_pC3.gbk", "genbank")
            kept = len(recs)
    summary.append((g, len(contigs), n, kept, sum(len(r) for r in recs) if kept else 0))
    print(f"[2] {g}: {len(contigs)} pC3 contigs, {n} proteins, {kept} GenBank records",
          flush=True)

with open(f"{W}/screen/input_summary.tsv", "w") as fh:
    fh.write("genome\tn_pC3_contigs\tn_pC3_proteins\tn_genbank_records\tpC3_bp\n")
    for r in summary:
        fh.write("\t".join(str(x) for x in r) + "\n")
print("[3] done", flush=True)
