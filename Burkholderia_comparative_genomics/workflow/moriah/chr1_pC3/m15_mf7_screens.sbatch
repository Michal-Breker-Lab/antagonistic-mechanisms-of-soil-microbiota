#!/bin/bash
#SBATCH --job-name=c3_mf7scr
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=0-04:00:00
#SBATCH --output=logs/m15_mf7scr.%j.out
#SBATCH --error=logs/m15_mf7scr.%j.err
# NO --partition (the -k variants are account-restricted).
#
# Re-run the three fig7/fig8 screens for MF7 ONLY, against its COMPLETE pC3.
#
# Job 45911143 screened MF7 using setb_pc3/MF7_pC3.* , which s37 carved from a
# single contig (NODE_3, 857,810 bp). MF7's pC3 is actually TWO contigs --
# NODE_3 + NODE_8 (318,422 bp) = 1,176,232 bp, 1,011 CDS, i.e. 99.7% of MF6's
# pC3 -- so that run saw 73% of the replicon and MF7's counts (5 antiSMASH
# regions vs MF6's 7; 21 MacSyFinder rows vs 37) are artefacts of the missing
# file, not of MF7's assembly. Inputs now come from setb_pc3_v2 (s41).
#
# The other five Set B genomes are closed single-contig replicons, so their
# 45911143 results stand unchanged and are NOT recomputed here.
#
# MacSyFinder: --db-type ordered_replicon assumes ONE ordered replicon. MF7's
# pC3 is two contigs, so the two are screened SEPARATELY and merged, rather than
# concatenated (which would invent co-localisation across the junction) or run
# --db-type unordered (which would drop the co-localisation criterion entirely
# and make MF7 incomparable to the other five).
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
IPS=/sci/backup/ofinkel/moshea/Tools/interproscan-5.76-107.0
IN="$W/setb_pc3_v2"; OUT="$W/screens_mf7"; DB="$W/databases/antismash_db"
MODELS="$W/macsy_models"
T=${SLURM_CPUS_PER_TASK:-16}
cd "$W"; mkdir -p "$OUT"

export TMPDIR=$(mktemp -d "$W/tmp/mf7scr.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight ==="
for e in fna gbk faa gff3; do
    [ -s "$IN/MF7_pC3.$e" ] || { echo "PREFLIGHT FAIL: missing $IN/MF7_pC3.$e" >&2; exit 1; }
done
NC=$(grep -c '^>' "$IN/MF7_pC3.fna")
[ "$NC" -eq 2 ] || { echo "PREFLIGHT FAIL: $NC contigs in MF7_pC3.fna, expected 2" >&2; exit 1; }
NP=$(grep -c '^>' "$IN/MF7_pC3.faa")
[ "$NP" -eq 1011 ] || { echo "PREFLIGHT FAIL: $NP proteins, expected 1011" >&2; exit 1; }
[ -x "$IPS/interproscan.sh" ] || { echo "PREFLIGHT FAIL: no interproscan.sh" >&2; exit 1; }
[ -d "$DB" ] || { echo "PREFLIGHT FAIL: no antiSMASH db at $DB" >&2; exit 1; }
[ -d "$MODELS/TXSScan/definitions/bacteria/diderm" ] || { echo "PREFLIGHT FAIL: no TXSScan models" >&2; exit 1; }
echo "  MF7 pC3: $NC contigs, $NP proteins; all inputs ok"

source "$M/etc/profile.d/conda.sh"

# ------------------------------------------------------------- InterProScan
# Same flags as m13 so the rows merge with the other five genomes' output.
echo "=== [1/3] InterProScan ($(date +%H:%M)) ==="
"$IPS/interproscan.sh" -i "$IN/MF7_pC3.faa" -b "$OUT/MF7_pC3_ips" \
    -f TSV,GFF3 -goterms -pa -iprlookup -dp -cpu "$T" -T "$TMPDIR"
echo "  InterProScan rows: $(wc -l < "$OUT/MF7_pC3_ips.tsv")"

# ---------------------------------------------------------------- antiSMASH
# antiSMASH handles a multi-record GenBank natively (one record per contig).
echo "=== [2/3] antiSMASH ($(date +%H:%M)) ==="
conda activate "$W/envs/antismash7"
antismash --taxon bacteria --genefinding-tool none --cpus "$T" \
    --databases "$DB" --cb-general --cb-knownclusters --asf --pfam2go \
    --output-dir "$OUT/antismash_MF7" "$IN/MF7_pC3.gbk"
conda deactivate

# -------------------------------------------------------------- MacSyFinder
echo "=== [3/3] MacSyFinder, per contig ($(date +%H:%M)) ==="
conda activate "$W/envs/macsy"
CTGS=$(grep '^>' "$IN/MF7_pC3.fna" | cut -c2- | awk '{print $1}')
i=0
for c in $CTGS; do
    i=$((i+1))
    echo "  --- contig $i: $c ---"
    # proteins on this contig only, in annotation order (ordered_replicon needs order)
    "$M/bin/python" - "$IN/MF7_pC3.gff3" "$IN/MF7_pC3.faa" "$c" \
        > "$TMPDIR/MF7_c$i.faa" <<'PY'
import sys
gff, faa, contig = sys.argv[1], sys.argv[2], sys.argv[3]
tags = []
for line in open(gff):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != contig or f[2] != "CDS":
        continue
    lt = next((a.split("=", 1)[1] for a in f[8].split(";") if a.startswith("locus_tag=")), None)
    if lt:
        tags.append((int(f[3]), lt))
tags.sort()
seqs, name, buf = {}, None, []
for line in open(faa):
    if line.startswith(">"):
        if name:
            seqs[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
if name:
    seqs[name] = "".join(buf)
for _, t in tags:
    if t in seqs:
        sys.stdout.write(f">{t}\n{seqs[t]}\n")
PY
    echo "     $(grep -c '^>' "$TMPDIR/MF7_c$i.faa") proteins"
    macsyfinder --db-type ordered_replicon --sequence-db "$TMPDIR/MF7_c$i.faa" \
        --models TXSScan/bacteria/diderm all --models-dir "$MODELS" \
        -o "$OUT/macsy_MF7_c$i" -w "$T"
done
conda deactivate

echo "=== MF7 SCREENS DONE ($(date +%H:%M)) ==="
ls -la "$OUT"
exit
