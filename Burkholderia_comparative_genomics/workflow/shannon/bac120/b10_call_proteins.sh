#!/bin/bash
# Gene calls for the bac120 genomes that the original pipeline never annotated.
#
# WHY: s23 called genes for the 771 DEREPLICATED genomes. The bac120 set is not
# dereplicated, so 142 of its 792 genomes have no proteins and therefore were
# never searched for the 10 toxin/immunity loci -- which is why 38 of Figure 12b's
# 45 tips are cross-hatched "never assessed" rather than coloured.
#
# Settings are copied from s23_pyrodigal.py, NOT re-chosen: single mode with
# translation_table=11, falling back to meta mode below 20 kb, and the same
# ">acc|contig|i|begin|end|strand" header format the blast parser expects. A
# different gene caller or different settings would make these genomes' hit
# counts incomparable with the 771 already in the table.
set -euo pipefail
C=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
PY=$C/envs/bakta/bin/python          # pyrodigal 3.7.1, the env s23 used
OUT=$C/pyrodigal_prot                # SAME directory, so the set stays unified
mkdir -p "$OUT" "$W/logs"
cd "$W"

echo "=== preflight ==="
"$PY" -c "import pyrodigal; print('pyrodigal', pyrodigal.__version__)"
ls genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort > "$W/all_bac120.txt"
ls "$OUT"/*.faa 2>/dev/null | xargs -n1 basename | sed 's/\.faa$//' | sort > "$W/have_prot.txt"
comm -23 "$W/all_bac120.txt" "$W/have_prot.txt" > "$W/need_prot.txt"
echo "  bac120 genomes: $(wc -l < "$W/all_bac120.txt")"
echo "  already called: $(wc -l < "$W/have_prot.txt")"
echo "  to call:        $(wc -l < "$W/need_prot.txt")"

"$PY" - <<'PY'
import os, sys, glob
from concurrent.futures import ProcessPoolExecutor
import pyrodigal

C = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
W = "/mnt/LargeStorageNoBackup/Moshea/burk_bac120"
OUT = f"{C}/pyrodigal_prot"


def read_fasta(p):
    recs, name, seq = [], None, []
    with open(p) as fh:
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


def call_one(acc):
    fna = f"{W}/genomes/{acc}.fna"
    dest = f"{OUT}/{acc}.faa"
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return (acc, "cached", 0)
    recs = read_fasta(fna)
    seqs = [s for _, s in recs]
    total = sum(len(s) for s in seqs)
    mode = "single"
    try:
        if total < 20000:
            raise ValueError("below single-mode minimum")
        gf = pyrodigal.GeneFinder(meta=False)
        gf.train(*seqs, translation_table=11)
    except ValueError:
        mode = "meta"
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
    return (acc, mode, n)


need = [l.strip() for l in open(f"{W}/need_prot.txt") if l.strip()]
print(f"calling genes for {len(need)} genomes", flush=True)
done = 0
with ProcessPoolExecutor(max_workers=32) as ex:
    for acc, mode, n in ex.map(call_one, need):
        done += 1
        if done % 25 == 0 or done == len(need):
            print(f"  {done}/{len(need)}  last {acc} {mode} {n} proteins", flush=True)
print("gene calling done", flush=True)
PY

echo "=== verify ==="
MISS=0
while read -r a; do
    [ -s "$OUT/$a.faa" ] || { echo "  MISSING proteins: $a" >&2; MISS=$((MISS+1)); }
done < "$W/all_bac120.txt"
[ "$MISS" -eq 0 ] || { echo "FAIL: $MISS genomes without proteins" >&2; exit 1; }
echo "  all $(wc -l < "$W/all_bac120.txt") bac120 genomes have protein calls"
echo "=== PROTEINS OK ==="
