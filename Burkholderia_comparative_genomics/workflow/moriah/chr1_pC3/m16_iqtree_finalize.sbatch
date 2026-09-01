#!/bin/bash
#SBATCH --job-name=c3_iqfin
#SBATCH --partition=glacier
# --partition is set explicitly for the same reason as m10: the site job_submit
# plugin stamps Features=cpuonly but can leave the GPU partitions in the list,
# and SLURM 20.11 then rejects the request as BadConstraints instead of falling
# back to glacier.
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-06:00:00
#SBATCH --output=logs/m16_iqfin.%j.out
#SBATCH --error=logs/m16_iqfin.%j.err
set -euo pipefail
#
# FINALISE the chr1 tree by RESUMING job 45910701 from its checkpoint.
#
# 45910701 ran 1d23h58m and was cancelled deliberately. Its tree search had
# converged -- no better tree since iteration 27, 103 unsuccessful iterations --
# but UFBoot had not:
#
#   NOTE: Bootstrap correlation coefficient of split occurrence frequencies: 0.983
#   NOTE: UFBoot does not converge, continue at least 100 more iterations
#
# so at iteration ~127 it granted itself another 100 iterations and its own ETA
# started RISING (9h46m at iteration 100 -> 24h37m at iteration 130). The likely
# cause is in the alignment and was flagged at startup: 21 sequences carry >50%
# gaps/ambiguity, so their splits move between replicates and the correlation
# stalls just under the 0.99 default.
#
# Moshe's decision (2026-08-25): accept 0.983. A plain cancel does NOT deliver
# that -- IQ-TREE writes .contree/.splits.nex/.iqtree only on completion, so the
# support values existed solely in the checkpoint. This run lowers --bcor below
# where the correlation already sits and shortens the check interval, so the
# next check declares convergence and the finalised files get written.
#
# NOTE THE ABSENCE OF -redo: that flag would discard the checkpoint and restart
# from scratch, throwing away two days of compute. It is why m10 is not simply
# re-run.
W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
IQ="$M/envs/anvio-9/bin/iqtree3"
T=${SLURM_CPUS_PER_TASK:-64}
ALN="$W/results/chr1_core_763.aln"
PREFIX="$W/results/tree_chr1"
cd "$W"
export TMPDIR=$(mktemp -d "$W/tmp/iqfin.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight ==="
[ -s "$ALN" ]            || { echo "PREFLIGHT FAIL: missing $ALN" >&2; exit 1; }
[ -s "$PREFIX.ckp.gz" ]  || { echo "PREFLIGHT FAIL: no checkpoint at $PREFIX.ckp.gz -- resuming is the whole point" >&2; exit 1; }
squeue -j 45910701 -h -o "%T" | grep -q RUNNING && { echo "PREFLIGHT FAIL: 45910701 still running; two writers on one prefix would corrupt it" >&2; exit 1; }
echo "  alignment: $(grep -c '^>' "$ALN") taxa"
echo "  checkpoint: $(stat -c%s "$PREFIX.ckp.gz") bytes, $(stat -c%y "$PREFIX.ckp.gz")"
"$IQ" --version | head -1

# --bcor 0.98 : accept the 0.983 already reached (default is 0.99)
# --nstep 10  : check convergence every 10 iterations, not every 100
"$IQ" -s "$ALN" -m GTR+F+I+R5 -B 1000 -T "$T" \
      --bcor 0.98 --nstep 10 --prefix "$PREFIX"

echo "=== outputs ==="
for f in treefile contree splits.nex iqtree; do
    [ -s "$PREFIX.$f" ] || { echo "FAIL: $PREFIX.$f not written" >&2; exit 1; }
    echo "  $PREFIX.$f  $(stat -c%s "$PREFIX.$f") bytes"
done
grep -E "correlation coefficient|converge" "$PREFIX.log" | tail -3
echo "=== FINALISED ($(date +%H:%M)) ==="
exit
