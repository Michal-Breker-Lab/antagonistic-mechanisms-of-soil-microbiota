#!/bin/bash
#SBATCH --job-name=bac120_dl
#SBATCH --partition=glacier
#SBATCH --exclude=glacier-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/b3_download.%j.out
#SBATCH --error=logs/b3_download.%j.err
#
# --partition is set explicitly, against the usual house rule, for the reason
# documented in m10/m16: the site job_submit plugin stamps Features=cpuonly but
# sometimes leaves the GPU partitions in the list, and SLURM 20.11 then rejects
# the whole request as BadConstraints instead of falling back to glacier.
# glacier-11 cannot see /sci/backup and fails with "couldn't chdir".
#
# --- WHY THIS SCRIPT IS RESUMABLE (job 45943757 post-mortem) ------------------
# The first attempt died after 3 min on chunk.ag with
#     Error: Download error: stream error: stream ID 5; INTERNAL_ERROR;
#            received from peer
# a TRANSIENT HTTP/2 error from NCBI, not a bad request. Two defects showed up:
#   1. no retry -- one transient network blip killed the whole job;
#   2. staging lived in TMPDIR, which the EXIT trap deleted, so the six chunks
#      that HAD downloaded were thrown away and the retry started from zero.
#
# The SECOND attempt (45943765) then failed differently, and this is the important
# one: `datasets` exited 0, wrote nothing alarming to stderr, and produced a
# CORRUPT 209 MB zip. unzip died with "invalid compressed data to inflate" and
# took the job with it. Verified afterwards: chunk.aa.zip passes `unzip -t`,
# chunk.ab.zip fails it -- so the archive was bad, not the extractor.
#
# So exit code + stderr are NOT sufficient acceptance criteria for a download.
# This is the same silent-corruption-with-clean-exit class as NCBI's truncated
# streams. Fixed here: every archive must pass `unzip -t` before it is accepted,
# and a failing archive is deleted and re-fetched rather than killing the job.
# Chunks are also halved to 50 accessions, so a corrupt archive costs half as
# much to re-fetch and the blast radius is smaller.
# The name-by-name assertions at the end are unchanged and remain the final guard.
set -euo pipefail
W=/sci/backup/ofinkel/moshea/burk_bac120
M=/sci/backup/ofinkel/moshea/miniforge3
C=/sci/backup/ofinkel/moshea/burkholderia_c3      # the existing c3 work dir
STAGE="$W/dl_stage"                               # PERSISTENT: survives the job
cd "$W"
mkdir -p "$STAGE/zips" "$STAGE/x" genomes logs
# TMPDIR still goes on /sci/backup -- Moriah's / is ~100% full -- but nothing we
# need to keep is written there any more.
export TMPDIR=$(mktemp -d "$W/tmp/dl.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT

source "$M/etc/profile.d/conda.sh"
conda activate datasets

echo "=== preflight ==="
[ -s accessions.txt ] || { echo "PREFLIGHT FAIL: no accessions.txt" >&2; exit 1; }
N=$(wc -l < accessions.txt)
[ "$N" -eq 790 ] || { echo "PREFLIGHT FAIL: $N accessions, expected 790" >&2; exit 1; }
for f in "$C/MF6.fna" "$C/MF7.fna"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing lab isolate $f" >&2; exit 1; }
done
datasets --version
echo "  accessions: $N   stage: $STAGE"

# --- download in chunks, with retry and resume ------------------------------
# CHUNKED ON PURPOSE: a single stream of 790 truncates at round counts and still
# exits 0 -- this project lost 40% of an Acidobacteria download that way.
split -l 50 accessions.txt "$STAGE/chunk."
RETRIES=6
for c in "$STAGE"/chunk.*; do
    case "$c" in *.done|*.err) continue;; esac
    b=$(basename "$c")
    if [ -f "$STAGE/$b.done" ]; then
        echo "  $b already done, skipping"
        continue
    fi
    ok=0
    for try in $(seq 1 $RETRIES); do
        rm -f "$STAGE/zips/$b.zip"
        if ! datasets download genome accession --inputfile "$c" --include genome \
                 --no-progressbar --filename "$STAGE/zips/$b.zip" \
                 2> "$STAGE/$b.err"; then
            echo "  $b attempt $try: datasets exited non-zero, retrying" >&2
            sed -n '1,3p' "$STAGE/$b.err" >&2
            sleep $(( try * 20 )); continue
        fi
        # A clean exit is NOT sufficient. Three independent checks, because this
        # download path has now failed two different ways with exit 0:
        #   (a) stderr can carry a stream/truncation error despite exit 0
        if grep -qiE "INTERNAL_ERROR|internal error|truncat" "$STAGE/$b.err"; then
            echo "  $b attempt $try: stderr reports a stream error, retrying" >&2
            sleep $(( try * 20 )); continue
        fi
        #   (b) the archive itself can be CORRUPT despite exit 0 and clean stderr
        if ! unzip -t "$STAGE/zips/$b.zip" >/dev/null 2>&1; then
            echo "  $b attempt $try: archive fails unzip -t (corrupt), retrying" >&2
            sleep $(( try * 20 )); continue
        fi
        #   (c) extraction must actually succeed before the chunk is accepted
        if ! unzip -qo "$STAGE/zips/$b.zip" -d "$STAGE/x" 2> "$STAGE/$b.unzip.err"; then
            echo "  $b attempt $try: extraction failed, retrying" >&2
            sed -n '1,3p' "$STAGE/$b.unzip.err" >&2
            sleep $(( try * 20 )); continue
        fi
        ok=1; break
    done
    [ "$ok" -eq 1 ] || { echo "FAIL: $b still failing after $RETRIES attempts" >&2; exit 1; }
    touch "$STAGE/$b.done"
    echo "  $b ok"
done

# --- flatten to one <accession>.fna per genome ------------------------------
# datasets lays out ncbi_dataset/data/<ACCESSION>/<file>.fna, so the accession is
# the PARENT DIRECTORY name.
echo "=== flattening ==="
find "$STAGE/x" -name "*.fna" -print0 | while IFS= read -r -d '' f; do
    acc=$(basename "$(dirname "$f")")
    cp -n "$f" "genomes/$acc.fna"
done
cp "$C/MF6.fna" "$C/MF7.fna" genomes/

echo "=== verify ==="
GOT=$(ls genomes/*.fna 2>/dev/null | wc -l)
echo "  genomes present: $GOT (expected 792 = 790 + MF6 + MF7)"
# Name-by-name, not just a count: a duplicate plus a miss would pass a bare count.
MISS=0
while read -r a; do
    [ -s "genomes/$a.fna" ] || { echo "  MISSING: $a" >&2; MISS=$((MISS+1)); }
done < accessions.txt
[ "$MISS" -eq 0 ] || { echo "FAIL: $MISS requested accessions absent" >&2; exit 1; }
[ "$GOT" -eq 792 ] || { echo "FAIL: $GOT != 792 -- incomplete set" >&2; exit 1; }
for f in genomes/MF6.fna genomes/MF7.fna; do
    printf "  %s  %s contigs\n" "$f" "$(grep -c '^>' "$f")"
done
echo "=== DOWNLOAD OK ==="
exit
