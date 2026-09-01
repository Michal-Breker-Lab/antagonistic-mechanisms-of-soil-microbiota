#!/bin/bash
# ANI of MF6 (and MF7) against ALL 792 bac120 genomes.
#
# WHY THIS IS NEEDED: rebuild/tables/MF6_ani_raw.tsv covers only 553 genomes,
# because it was computed against the DEREPLICATED genome set. The bac120 set
# does not dereplicate, so 238 of its genomes have NO ANI measurement at all --
# and 53 of those fall inside the MRCA of MF6's >=94% relatives, 37 of them
# B. sola, MF6's own species. Selecting "genomes >=94% ANI to MF6" out of a set
# where 30% are unmeasured excludes close relatives for want of a number, which
# is exactly the selection artefact that makes a close-relatives panel mislead.
#
# EXACTLY the original invocation, recovered from run_full_rebuild.sh:139-142:
#     skani dist -q MF6.fna MF7.fna -r genomes/*.fna      (plain defaults)
# Two traps were hit getting here, both worth recording:
#   1. `--medium` shifts ANI by a median +1.41 vs the default mode. ANY extra
#      skani flag makes these numbers incomparable with the values already
#      printed in Figure 12 and ranked in 6.3.
#   2. The original query is the WHOLE MF6 genome, not chromosome 1. The
#      `cluster_001_consensus` in the table's Query_name column is merely the
#      contig that anchored the best alignment -- it is output, not input.
# Change ONLY the reference set, which is the thing that was actually incomplete.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
SK=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/envs/burk/bin/skani
cd "$W"
mkdir -p logs
CPUS=${CPUS:-32}

echo "=== preflight ==="
[ -x "$SK" ] || { echo "FAIL: project skani not at $SK" >&2; exit 1; }
"$SK" --version 2>&1 | head -1
N=$(ls genomes/*.fna | wc -l)
for g in MF6 MF7; do
    [ -s "genomes/$g.fna" ] || { echo "FAIL: no genomes/$g.fna" >&2; exit 1; }
done
echo "  references: $N   threads: $CPUS"

"$SK" dist -t "$CPUS" -q genomes/MF6.fna genomes/MF7.fna -r genomes/*.fna \
      > mf6_ani_full.tsv 2> logs/skani_full.log

echo "=== verify ==="
[ -s mf6_ani_full.tsv ] || { echo "FAIL: no output" >&2; tail -5 logs/skani_full.log >&2; exit 1; }
ROWS=$(( $(wc -l < mf6_ani_full.tsv) - 1 ))
echo "  rows: $ROWS (queries MF6 + MF7 x $N references)"
head -1 mf6_ani_full.tsv
# skani reports only pairs above its screen, so ROWS < 2N is normal and means
# "below the screen", not "failed". Say so rather than asserting equality.
echo "  (skani reports only pairs above its screen; absent = below it, not an error)"
awk -F'\t' 'NR>1 && $2 ~ /MF6\.fna/ && $3+0 >= 94 {n++} END {print "  MF6 >=94% ANI: " n+0}' mf6_ani_full.tsv
echo "=== ANI OK ==="
