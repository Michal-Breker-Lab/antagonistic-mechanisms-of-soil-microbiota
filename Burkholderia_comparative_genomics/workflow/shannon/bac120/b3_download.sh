#!/bin/bash
# Fetch the bac120 genomes MISSING from the c3 dereplicated set, and assemble the
# full 790 by symlinking the 653 already on disk.
#
# WHY THE 137 ARE NOT OPTIONAL: the existing c3 genomes/ dir is the DEREPLICATED
# set. The 137 absent accessions are exactly what dereplication removed -- 81% of
# them sit in species with >=10 genomes in the target set (B. sola 37, B.
# multivorans 30, B. pseudomallei 17) against 61% for the 653 already present.
# Building on the 653 alone would silently re-impose the dereplication this tree
# exists to avoid. They are the non-redundancy, not padding.
#
# ACCEPTANCE CRITERIA (learned from Moriah jobs 45943757 and 45943765): `datasets`
# has now failed TWICE with exit code 0 -- once with a truncated HTTP/2 stream,
# once writing a CORRUPT 209 MB zip with clean stderr. Exit code and stderr are
# therefore NOT sufficient. Every chunk must pass, in order: clean stderr, a
# successful `unzip -t`, and a successful extraction. Any failure re-fetches
# rather than killing the run.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
C=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
DL=/home/moshea/miniconda/envs/ncbi_datasets/bin
STAGE="$W/dl_stage"
cd "$W"; mkdir -p "$STAGE/zips" "$STAGE/x" genomes logs
export PATH="$DL:$PATH"

echo "=== preflight ==="
[ -s accessions.txt ] || { echo "FAIL: no accessions.txt" >&2; exit 1; }
[ -s need.txt ]       || { echo "FAIL: no need.txt" >&2; exit 1; }
N=$(wc -l < accessions.txt); NEED=$(wc -l < need.txt)
[ "$N" -eq 790 ] || { echo "FAIL: $N accessions, expected 790" >&2; exit 1; }
for f in "$C/MF6.fna" "$C/MF7.fna"; do
    [ -s "$f" ] || { echo "FAIL: missing lab isolate $f" >&2; exit 1; }
done
datasets --version
echo "  target $N   already on disk $((N-NEED))   to download $NEED"

# --- link the 653 that already exist (verified identical to fresh NCBI copies) --
echo "=== linking existing genomes ==="
LINKED=0
while read -r a; do
    if [ -s "$C/genomes/$a.fna" ] && [ ! -e "genomes/$a.fna" ]; then
        ln -s "$C/genomes/$a.fna" "genomes/$a.fna"; LINKED=$((LINKED+1))
    fi
done < accessions.txt
echo "  linked: $LINKED"

# --- download the rest, in small chunks, with retry + resume ------------------
split -l 25 need.txt "$STAGE/chunk."
RETRIES=6
for c in "$STAGE"/chunk.*; do
    case "$c" in *.done|*.err) continue;; esac
    b=$(basename "$c")
    [ -f "$STAGE/$b.done" ] && { echo "  $b already done"; continue; }
    ok=0
    for try in $(seq 1 $RETRIES); do
        rm -f "$STAGE/zips/$b.zip"
        if ! datasets download genome accession --inputfile "$c" --include genome \
                 --no-progressbar --filename "$STAGE/zips/$b.zip" 2> "$STAGE/$b.err"; then
            echo "  $b try $try: datasets exited non-zero" >&2
            sed -n '1,3p' "$STAGE/$b.err" >&2; sleep $((try*15)); continue
        fi
        if grep -qiE "INTERNAL_ERROR|internal error|truncat" "$STAGE/$b.err"; then
            echo "  $b try $try: stream error in stderr despite exit 0" >&2
            sleep $((try*15)); continue
        fi
        if ! unzip -t "$STAGE/zips/$b.zip" >/dev/null 2>&1; then
            echo "  $b try $try: archive CORRUPT despite exit 0 + clean stderr" >&2
            sleep $((try*15)); continue
        fi
        if ! unzip -qo "$STAGE/zips/$b.zip" -d "$STAGE/x" 2> "$STAGE/$b.unzip.err"; then
            echo "  $b try $try: extraction failed" >&2
            sed -n '1,3p' "$STAGE/$b.unzip.err" >&2; sleep $((try*15)); continue
        fi
        ok=1; break
    done
    [ "$ok" -eq 1 ] || { echo "FAIL: $b failed $RETRIES times" >&2; exit 1; }
    touch "$STAGE/$b.done"; echo "  $b ok"
done

# datasets lays out ncbi_dataset/data/<ACCESSION>/<file>.fna -- the accession is
# the PARENT DIRECTORY, not the filename.
echo "=== flattening ==="
find "$STAGE/x" -name "*.fna" -print0 | while IFS= read -r -d '' f; do
    acc=$(basename "$(dirname "$f")")
    [ -e "genomes/$acc.fna" ] || cp "$f" "genomes/$acc.fna"
done
cp -n "$C/MF6.fna" "$C/MF7.fna" genomes/ 2>/dev/null || true

echo "=== verify ==="
# Name by name, not a bare count: one duplicate plus one miss passes a count check.
MISS=0
while read -r a; do
    [ -s "genomes/$a.fna" ] || { echo "  MISSING: $a" >&2; MISS=$((MISS+1)); }
done < accessions.txt
GOT=$(ls genomes/*.fna 2>/dev/null | wc -l)
echo "  genomes present: $GOT (expected 792 = 790 + MF6 + MF7)"
[ "$MISS" -eq 0 ] || { echo "FAIL: $MISS requested accessions absent" >&2; exit 1; }
[ "$GOT" -eq 792 ] || { echo "FAIL: $GOT != 792" >&2; exit 1; }
for g in MF6 MF7; do
    printf "  %-4s %s contigs\n" "$g" "$(grep -c '^>' "genomes/$g.fna")"
done
echo "=== DOWNLOAD OK ==="
