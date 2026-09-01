#!/bin/bash
# GTDB-Tk identify + align on the 792-genome set -> concatenated bac120 MSA.
#
# WHY NOT de_novo_wf: that wrapper also runs root/decorate, which read the
# taxonomy and RED files and make the output depend on the GTDB release in ways
# that are hard to audit. identify+align read only the marker HMMs and the
# release's column mask. This tree is rooted on its own topology, not by gtdbtk.
#
# --skip_gtdb_refs IS LOAD-BEARING: without it, thousands of GTDB reference
# genomes join the MSA and the tree stops being about our 792.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
P=/mnt/LargeStorageNoBackup/Moshea/envs/gtdbtk-2.7.2
cd "$W"
export GTDBTK_DATA_PATH=/mnt/scratch/gtdbtk_database/release232
export PATH="$P/bin:$PATH"
# Shannon's / is 95% full and /tmp lives on it. gtdbtk writes large per-genome
# temp files; keep them on LargeStorage or the run dies with an opaque error.
export TMPDIR="$W/tmp"; mkdir -p "$TMPDIR" logs
CPUS=${CPUS:-64}

echo "=== preflight ==="
gtdbtk --version | head -1
N=$(ls genomes/*.fna 2>/dev/null | wc -l)
echo "  genomes: $N   cpus: $CPUS   TMPDIR: $TMPDIR"
[ "$N" -eq 792 ] || { echo "FAIL: $N genomes, expected 792" >&2; exit 1; }
# MF7 is the methodological point of this tree -- fail here, not 3 h in.
for g in MF6 MF7; do
    [ -s "genomes/$g.fna" ] || { echo "FAIL: genomes/$g.fna absent" >&2; exit 1; }
done

# --- namespace the inputs -----------------------------------------------------
# gtdbtk REFUSES any input whose id matches a GTDB reference id, and 68 of our
# accessions ARE GTDB representative genomes (expected: this is a Burkholderia
# set and GTDB picks representatives from the same assemblies). The reference id
# set is the raw accession, so any prefix clears it.
# Prefix ALL 792, not just the 68: a two-class naming scheme means every
# downstream join has to know which class a tip is in, and that asymmetry is
# exactly the kind of thing that silently mislabels a tree. The prefix is
# stripped from the MSA headers at the end, so tip labels are unchanged.
echo "=== namespacing inputs (u_ prefix) ==="
mkdir -p genomes_u
for f in genomes/*.fna; do
    b=$(basename "$f" .fna)
    [ -e "genomes_u/u_$b.fna" ] || ln -s "$(readlink -f "$f")" "genomes_u/u_$b.fna"
done
NU=$(ls genomes_u/*.fna | wc -l)
[ "$NU" -eq 792 ] || { echo "FAIL: $NU namespaced inputs != 792" >&2; exit 1; }
echo "  namespaced: $NU"

if [ ! -s identify/identify/gtdbtk.bac120.markers_summary.tsv ]; then
    echo "=== gtdbtk identify ==="
    gtdbtk identify --genome_dir genomes_u --out_dir identify --extension fna --cpus "$CPUS"
else
    echo "=== identify already done, skipping ==="
fi
S=identify/identify/gtdbtk.bac120.markers_summary.tsv
[ -s "$S" ] || { echo "FAIL: no marker summary" >&2; exit 1; }
ROWS=$(( $(wc -l < "$S") - 1 ))
echo "  marker summary rows: $ROWS (expected 792)"
[ "$ROWS" -eq 792 ] || { echo "FAIL: $ROWS rows != 792" >&2; exit 1; }
for g in MF6 MF7; do grep -qP "^u_$g\t" "$S" || { echo "FAIL: $g absent from summary" >&2; exit 1; }; done

echo "=== gtdbtk align (--skip_gtdb_refs) ==="
gtdbtk align --identify_dir identify --out_dir align --skip_gtdb_refs --cpus "$CPUS"

echo "=== verify ==="
A=align/align/gtdbtk.bac120.user_msa.fasta.gz
[ -s "$A" ] || { echo "FAIL: no user MSA at $A" >&2; exit 1; }
# Strip the u_ prefix so tip labels are the real accessions / MF6 / MF7.
gunzip -c "$A" | sed 's/^>u_/>/' > msa_bac120.faa
grep -q "^>u_" msa_bac120.faa && { echo "FAIL: u_ prefix survived into the MSA" >&2; exit 1; }
TAXA=$(grep -c '^>' msa_bac120.faa)
COLS=$(awk '/^>/{if(n){print length(s); exit} n=1; next}{s=s $0}' msa_bac120.faa)
echo "  taxa: $TAXA   columns: $COLS"
for g in MF6 MF7; do
    grep -q "^>$g" msa_bac120.faa || { echo "FAIL: $g absent from the MSA" >&2; exit 1; }
done
echo "  MF6 and MF7 both in the alignment"
# TAXA < 792 is NOT an error: align drops genomes with too few markers. Report the
# shortfall by name -- a silent drop reads as "everything passed".
if [ "$TAXA" -ne 792 ]; then
    echo "  NOTE: $((792-TAXA)) genome(s) dropped for insufficient markers:"
    comm -23 <(ls genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort) \
             <(grep '^>' msa_bac120.faa | sed 's/^>//' | awk '{print $1}' | sort) | sed 's/^/    /'
    comm -23 <(ls genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort) \
             <(grep '^>' msa_bac120.faa | sed 's/^>//' | awk '{print $1}' | sort) > dropped_genomes.txt
fi
echo "=== ALIGN OK ==="
