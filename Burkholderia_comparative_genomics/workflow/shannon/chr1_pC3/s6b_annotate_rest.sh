#!/usr/bin/env bash
# Annotate every genome in the full set that the dereplicated run never touched.
#
# Replicates s6_annotate.sh's invocation EXACTLY -- same Bakta binary, same DB,
# same flags -- because Panaroo clusters on gene calls and mixing annotation
# settings inside one run manufactures differences that look like biology.
# s6_annotate.sh itself cannot be reused directly: it takes no arguments and
# rebuilds its job list from results/derep_accessions.txt, i.e. the dereplicated
# set we are deliberately abandoning.
#
# ONE DELIBERATE DIFFERENCE: --complete is dropped for MF7. That flag asserts
# every sequence is a finished replicon; MF7 is a 63-contig SPAdes draft
# (8.70 Mb, N50 858 kb), so asserting completeness would mis-handle contig ends.
# Every other genome is a closed NCBI assembly, or MF6's closed 5-replicon
# assembly, and keeps --complete exactly as the original run had it.
set -uo pipefail

W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
export TMPDIR=$W/tmp
DB=/mnt/LargeStorageNoBackup/Datasets/Moshea/Databases/bakta_db/db
BAKTA=$W/envs/bakta/bin/bakta
# Bakta resolves tRNAscan-SE/aragorn/cmscan/pilercr/diamond from PATH, not from
# its own prefix. Absolute path alone fails with "tRNAscan-SE not found".
export PATH=$W/envs/bakta/bin:$PATH
JOBS=${JOBS:-14}          # concurrent genomes
THREADS=${THREADS:-8}     # 14 x 8 = 112 of 128 cores
LIST=${LIST:-$W/genome_list_full.txt}
CONFIRM=0
[ "${1:-}" = "--confirm" ] && CONFIRM=1
mkdir -p "$W/annot" "$W/logs/bakta" "$TMPDIR"

echo "=== environment ==="
for b in bakta diamond tRNAscan-SE aragorn cmscan; do
  p="$W/envs/bakta/bin/$b"
  if [ -x "$p" ]; then
    v=$("$p" --version 2>&1 | head -1)
    printf "  %-14s OK   %s\n" "$b" "$v"
  else
    printf "  %-14s MISSING at %s\n" "$b" "$p"
  fi
done
DV=$("$W/envs/bakta/bin/diamond" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "  diamond version parsed: ${DV:-unknown}"
case "$DV" in
  2.1.*) echo "  diamond 2.1.x - OK";;
  *)     echo "  WARNING: diamond is not 2.1.x. 2.2.x SIGSEGVs inside Bakta at the"
         echo "           sORF stage only, i.e. hours in, not at startup.";;
esac
[ -d "$DB" ] || { echo "FATAL: Bakta DB missing at $DB" >&2; exit 1; }
echo "  DB           : $DB"

echo
echo "=== job list ==="
[ -s "$LIST" ] || { echo "FATAL: $LIST missing (run stage 1 first)" >&2; exit 1; }
python3 - "$W" "$LIST" <<'PY'
import csv, os, sys
W, LIST = sys.argv[1], sys.argv[2]
want = [l.strip() for l in open(LIST) if l.strip()]
asm = {}
p = f"{W}/metadata/assemblies.tsv"
if os.path.exists(p):
    asm = {r["accession"]: r for r in csv.DictReader(open(p), delimiter="\t")}

rows, skipped = [], 0
for a in want:
    if os.path.getsize(f"{W}/annot/{a}/{a}.gff3") if os.path.exists(f"{W}/annot/{a}/{a}.gff3") else 0:
        skipped += 1
        continue
    if a == "MF6":
        # closed 5-replicon reassembly; B. sola by ANI
        rows.append((a, f"{W}/MF6.fna", "Burkholderia", "sola", "1"))
    elif a == "MF7":
        # 63-contig IMG draft -> NOT --complete. Species left as sp.: calling it
        # sola here would presume the ANI result this rebuild is meant to produce.
        rows.append((a, f"{W}/MF7.fna", "Burkholderia", "sp.", "0"))
    else:
        name = (asm.get(a, {}).get("organism_name") or "Burkholderia sp.").split()
        genus = name[0] if name else "Burkholderia"
        sp = name[1] if len(name) > 1 and not name[1].startswith("sp.") else "sp."
        rows.append((a, f"{W}/genomes/{a}.fna", genus, sp, "1"))

missing = [r[0] for r in rows if not os.path.exists(r[1])]
with open(f"{W}/annot_jobs_full.tsv", "w") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"  in full set        : {len(want)}")
print(f"  already annotated  : {skipped}")
print(f"  to annotate        : {len(rows)}")
print(f"  drafts (no --complete): {sum(1 for r in rows if r[4]=='0')}")
if missing:
    print(f"  FATAL: {len(missing)} have no FASTA: {missing[:5]}")
    sys.exit(1)
PY
[ $? -eq 0 ] || exit 1

N=$(wc -l < "$W/annot_jobs_full.tsv")
AVAIL=$(df -BG --output=avail /mnt/LargeStorageNoBackup | tail -1 | tr -dc '0-9')
NEED=$(( N / 20 + 10 ))
echo "  disk free          : ${AVAIL} GB (rough need ~${NEED} GB)"
[ "$AVAIL" -ge "$NEED" ] || { echo "FATAL: not enough disk" >&2; exit 1; }

echo
echo "=== invocation that will run (one genome shown) ==="
head -1 "$W/annot_jobs_full.tsv" | while IFS=$'\t' read -r acc fna genus sp comp; do
  cflag=""; [ "$comp" = "1" ] && cflag="--complete"
  echo "  $BAKTA --db \$DB --threads $THREADS --genus $genus --species $sp \\"
  echo "         --prefix $acc --output \$W/annot/$acc $cflag --keep-contig-headers \\"
  echo "         --skip-plot --force --tmp-dir \$TMPDIR $fna"
done
echo "  concurrency: $JOBS jobs x $THREADS threads = $((JOBS*THREADS)) of $(nproc) cores"

if [ "$CONFIRM" -ne 1 ]; then
  echo
  echo "Dry run. Re-run with --confirm to start."
  exit 0
fi

echo
echo "=== annotating ==="
date
run_one() {
  IFS=$'\t' read -r acc fna genus sp comp <<<"$1"
  out=$W/annot/$acc
  if [ -s "$out/$acc.gff3" ]; then echo "skip $acc"; return 0; fi
  mkdir -p "$out"
  cflag=""; [ "$comp" = "1" ] && cflag="--complete"
  $BAKTA --db $DB --threads $THREADS --genus "$genus" --species "$sp" \
         --prefix "$acc" --output "$out" $cflag --keep-contig-headers \
         --skip-plot --force --tmp-dir $TMPDIR "$fna" \
         > $W/logs/bakta/$acc.log 2>&1
  if [ -s "$out/$acc.gff3" ]; then echo "ok $acc"; else echo "FAIL $acc"; fi
}
export -f run_one; export W BAKTA DB THREADS TMPDIR PATH
cat "$W/annot_jobs_full.tsv" | xargs -d '\n' -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
date

echo
echo "=== verification ==="
bad=0; n=0
while read -r acc; do
  f="$W/annot/$acc/$acc.gff3"; n=$((n+1))
  if [ ! -s "$f" ]; then echo "  MISSING $acc"; bad=$((bad+1)); continue; fi
  c=$(grep -c $'\tCDS\t' "$f" 2>/dev/null); c=${c:-0}
  if [ "$c" -lt 2000 ]; then echo "  LOW CDS $acc ($c)"; bad=$((bad+1)); fi
done < "$LIST"
echo "  checked $n genomes, problems: $bad"
[ "$bad" -eq 0 ] || { echo "  FIX BEFORE CONTINUING" >&2; exit 1; }
echo "ANNOT_FULL_DONE"
