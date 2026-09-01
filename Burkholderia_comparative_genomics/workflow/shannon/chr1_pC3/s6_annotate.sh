#!/usr/bin/env bash
# Stage 6 - uniform Bakta annotation of the dereplicated set + MF6 itself.
#
# --keep-contig-headers is REQUIRED: Stage 7 subsets each GFF3 to the c3 contigs
# by their original NCBI accession, and Bakta renames contigs without it.
# --skip-plot saves several minutes per genome producing circular figures we
# never use.
# Re-runnable: a genome with a non-empty .gff3 is skipped.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
DB=/mnt/LargeStorageNoBackup/Datasets/Moshea/Databases/bakta_db/db
BAKTA=$W/envs/bakta/bin/bakta
# Bakta resolves tRNAscan-SE/aragorn/cmscan/pilercr/diamond from PATH, not from
# its own prefix -- calling the binary by absolute path alone makes it fail with
# "tRNAscan-SE not found". The env bin MUST be on PATH.
export PATH=$W/envs/bakta/bin:$PATH
JOBS=${JOBS:-14}          # concurrent genomes
THREADS=${THREADS:-8}     # threads each -> 112 cores of 128
mkdir -p $W/annot $W/logs/bakta

# ---- build the job list: accession <TAB> genus <TAB> species ----
python3 - <<'PY'
import csv, os
W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
keep = [l.strip() for l in open(f"{W}/results/derep_accessions.txt") if l.strip()]
asm = {r["accession"]: r for r in
       csv.DictReader(open(f"{W}/metadata/assemblies.tsv"), delimiter="\t")}
with open(f"{W}/annot_jobs.tsv", "w") as fh:
    for a in keep:
        name = (asm.get(a, {}).get("organism_name") or "Burkholderia sp.").split()
        genus = name[0] if name else "Burkholderia"
        sp = name[1] if len(name) > 1 and not name[1].startswith("sp.") else "sp."
        fh.write(f"{a}\t{W}/genomes/{a}.fna\t{genus}\t{sp}\n")
    # MF6 is not in NCBI; annotate it identically so it can join the pangenome
    fh.write(f"MF6\t{W}/MF6.fna\tBurkholderia\tsola\n")
print("jobs:", len(keep) + 1)
PY

run_one() {
  IFS=$'\t' read -r acc fna genus sp <<<"$1"
  out=$W/annot/$acc
  if [ -s "$out/$acc.gff3" ]; then echo "skip $acc"; return 0; fi
  mkdir -p "$out"
  $BAKTA --db $DB --threads $THREADS --genus "$genus" --species "$sp" \
         --prefix "$acc" --output "$out" --complete --keep-contig-headers \
         --skip-plot --force --tmp-dir $TMPDIR "$fna" \
         > $W/logs/bakta/$acc.log 2>&1
  if [ -s "$out/$acc.gff3" ]; then echo "ok $acc"; else echo "FAIL $acc"; fi
}
export -f run_one; export W BAKTA DB THREADS TMPDIR PATH

echo "=== annotating $(wc -l < $W/annot_jobs.tsv) genomes, $JOBS concurrent x $THREADS threads ==="
date
cat $W/annot_jobs.tsv | xargs -d '\n' -P $JOBS -I{} bash -c 'run_one "$@"' _ {}
date
echo "=== summary ==="
echo "gff3 produced: $(ls $W/annot/*/*.gff3 2>/dev/null | wc -l)"
echo "ANNOT_DONE"
