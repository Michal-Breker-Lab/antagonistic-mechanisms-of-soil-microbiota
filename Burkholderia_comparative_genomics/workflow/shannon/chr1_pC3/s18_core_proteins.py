#!/usr/bin/env python3
"""
s18 - write FASTA files of the B. sola (Set B) pC3 core proteins.

Three files, all one-sequence-per-gene-family so that downstream per-protein
tools are not fed five near-identical copies of everything:

  setB_pC3_core_proteins.faa        833 families, MF6 representative
  setB_pC3_core_noCOG.faa           the 594 with no COG category (the "unknowns")
  setB_pC3_core_allmembers.faa      all 5 members of every core family (4,165)

MF6 is used as the representative because it is present in every core family by
definition, and because it is the strain the project is about. Headers carry the
family id, the locus tag, the COG categories and the product so a hit can be
traced back without a join.
"""
import csv, os, re

_SRC = open("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/s17_functional.py").read()
assert "# 9. Redundancy" in _SRC
exec(_SRC.split("# 9. Redundancy")[0])

SETB = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1", "MF6"]
FAA = {g: f"{W}/annot/{g}/{g}.faa" for g in SETB}


def load_faa(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


SEQ = {}
for g, p in FAA.items():
    SEQ.update(load_faa(p))
print(f"loaded {len(SEQ)} protein sequences from 5 genomes", flush=True)


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def hdr(fam, lt, genome, L):
    cats = "".join(sorted(L["cats"])) or "-"
    prod = L["product"].replace("\t", " ") or "hypothetical protein"
    return f">{fam}|{lt}|{genome}|COG={cats}|{prod}"


n_rep = n_nocog = n_all = 0
missing = []
with open(f"{OUT}/setB_pC3_core_proteins.faa", "w") as f_rep, \
     open(f"{OUT}/setB_pC3_core_noCOG.faa", "w") as f_no, \
     open(f"{OUT}/setB_pC3_core_allmembers.faa", "w") as f_all:
    for fam in sorted(c3_core):
        L = LAB[("c3", fam)]
        cells = c3_fams[fam]
        rep = cells.get("MF6", [None])[0]
        if rep and rep in SEQ:
            f_rep.write(hdr(fam, rep, "MF6", L) + "\n" + wrap(SEQ[rep]) + "\n")
            n_rep += 1
            if not L["cats"]:
                f_no.write(hdr(fam, rep, "MF6", L) + "\n" + wrap(SEQ[rep]) + "\n")
                n_nocog += 1
        else:
            missing.append(fam)
        for g in SETB:
            for lt in cells.get(g, []):
                if lt in SEQ:
                    f_all.write(hdr(fam, lt, g, L) + "\n" + wrap(SEQ[lt]) + "\n")
                    n_all += 1

print(f"setB_pC3_core_proteins.faa    {n_rep} sequences (1 per core family, MF6)", flush=True)
print(f"setB_pC3_core_noCOG.faa       {n_nocog} sequences (no COG category)", flush=True)
print(f"setB_pC3_core_allmembers.faa  {n_all} sequences (all members)", flush=True)
if missing:
    print(f"WARNING: {len(missing)} core families had no usable MF6 protein: "
          f"{missing[:5]}", flush=True)

# InterProScan rejects '*' and non-standard residues; check before we hand it over
bad = 0
for fam in sorted(c3_core):
    lt = c3_fams[fam].get("MF6", [None])[0]
    if lt and lt in SEQ and re.search(r"[^ACDEFGHIKLMNPQRSTVWYX]", SEQ[lt]):
        bad += 1
print(f"sequences containing non-standard residues (incl. '*'): {bad}", flush=True)
