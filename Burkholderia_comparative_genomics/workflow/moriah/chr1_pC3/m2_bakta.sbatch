#!/bin/bash
#SBATCH --job-name=burk_bakta
#SBATCH --array=1-773%40
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/bakta_%A_%a.out
#SBATCH --error=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/bakta_%A_%a.err
# No --partition (account-restricted -k variants; SLURM picks by account).
# 8 cpus, not 16: on Shannon diamond used ~1 core regardless of --threads, so
# wider tasks bought nothing and just idled cores. 8 is the configuration that
# actually annotated 300+ genomes there.
# mem 64G / 8 cpus = 8 GB per core, the TOP of the $0.008/core-h base tier.
set -uo pipefail
W=/sci/backup/ofinkel/moshea/burkholderia_c3
DB=/sci/backup/aerez/aerez/moshea/bakta_db/db
source /sci/backup/ofinkel/moshea/miniforge3/etc/profile.d/conda.sh
conda activate bakta

# tRNAscan-SE in this env has an install-time path baked into BOTH its `use lib`
# line and its default config location, pointing at /sci/home/... where the env
# was built before being relocated to /sci/backup. Two earlier attempts failed:
#   1. conda activate alone      -> every task died at "predict tRNAs"
#   2. PERL5LIB + /sci/home symlink -> worked, but ~5% of tasks still died on
#      transient mount latency resolving the symlink (8 of the first 170).
# $W/bin/tRNAscan-SE is a wrapper that passes -c explicitly, so only
# /sci/backup is ever consulted. It must come FIRST on PATH.
export PATH="$W/bin:$PATH"

# Per-task scratch on the LAB filesystem, not node-local /tmp.
# Typed memory records a previous bakta array on this cluster losing 38% of its
# tasks to EDQUOT because node /tmp is per-user quota'd. Observed here: glacier
# nodes carry a 234 GB /tmp at ~1% use, so the hazard did not bite this run --
# but "it happened to fit" is not a design.
# NOTE: this block was accidentally deleted by an earlier patch of mine, which
# left `--tmp-dir "$TMPDIR"` pointing at SLURM's default /tmp for ~759 genomes.
# Those annotations are unaffected (ample space); restored so it cannot recur.
export TMPDIR=$(mktemp -d "$W/tmp/bakta.${SLURM_ARRAY_TASK_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT


# Preflight, with retries: a single transient filesystem stall must not cost a
# whole genome. Only give up if it fails repeatedly.
ok=0
for attempt in 1 2 3; do
  if tRNAscan-SE -h 2>&1 | grep -qi "tRNAscan-SE"; then ok=1; break; fi
  echo "preflight attempt $attempt failed; retrying"
  sleep $((attempt * 5))
done
if [ "$ok" -ne 1 ]; then
  echo "FATAL: tRNAscan-SE not runnable after 3 attempts on $(hostname)"
  exit 1
fi

read -r acc genus sp < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$W/genus_species.tsv" | tr '\t' ' ')
[ -n "${acc:-}" ] || { echo "no accession for task $SLURM_ARRAY_TASK_ID"; exit 1; }

case "$acc" in
  MF6) fna="$W/MF6.fna"; cflag="--complete" ;;
  # MF7 is a 63-contig draft: --complete would assert every sequence is a
  # finished replicon and mis-handle the contig ends.
  MF7) fna="$W/MF7.fna"; cflag="" ;;
  *)   fna="$W/genomes/$acc.fna"; cflag="--complete" ;;
esac

out="$W/annot/$acc"
if [ -s "$out/$acc.gff3" ]; then echo "skip $acc (already annotated)"; exit 0; fi
[ -s "$fna" ] || { echo "MISSING FASTA $acc"; exit 1; }
mkdir -p "$out"

echo "task $SLURM_ARRAY_TASK_ID  acc=$acc  genus=$genus  species=$sp  complete=${cflag:-no}"
bakta --db "$DB" --threads "${SLURM_CPUS_PER_TASK:-8}" \
      --genus "$genus" --species "$sp" \
      --prefix "$acc" --output "$out" $cflag --keep-contig-headers \
      --skip-plot --force --tmp-dir "$TMPDIR" "$fna"
rc=$?

if [ -s "$out/$acc.gff3" ]; then
  echo "ok $acc CDS=$(grep -c $'\tCDS\t' "$out/$acc.gff3")"
else
  echo "FAIL $acc (bakta rc=$rc)"; exit 1
fi
# Explicit exit: a bare echo as the last line returns 0 and masks failures
# from any afterok dependency.
exit 0
