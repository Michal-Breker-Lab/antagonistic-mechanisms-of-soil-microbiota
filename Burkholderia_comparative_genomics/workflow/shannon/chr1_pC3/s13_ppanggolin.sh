#!/usr/bin/env bash
# PPanGGOLiN re-clustering of the SAME Bakta GFF3 subsets Panaroo consumed.
# Two identities: 0.60 (requested by Moshe) and 0.80 (PPanGGOLiN default = control arm).
# Panaroo outputs are NOT touched -- this is an added method, not a replacement.
#
# Uses the PRE-EXISTING ppanggolin 2.2.6 env from the within_genus_hgt project.
# Do not build a new one: a fresh solve on 2026-08-09 resolved to 0.3.88.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
PPGENV=/mnt/scratch/within_genus_hgt/env_ppanggolin
PPG=$PPGENV/bin/ppanggolin
# PPanGGOLiN resolves mmseqs from PATH, not from its own prefix -- calling the
# binary by absolute path alone dies with "Command 'mmseqs' not found". Same
# class of bug as the recorded Bakta/tRNAscan-SE gotcha. The env bundles its own
# matching mmseqs; use that rather than envs/burk's 18.8cc5c.
export PATH=$PPGENV/bin:$PATH
CPU=${CPU:-32}
command -v mmseqs >/dev/null || { echo "FATAL: mmseqs still not on PATH"; exit 1; }
echo "mmseqs: $(command -v mmseqs) $(mmseqs version 2>&1 | head -1)"
mkdir -p "$W/ppanggolin" "$W/logs"

# ---- build the --anno list files: genome_name <TAB> path/to/gff3 ----
# Bakta GFF3s carry their FASTA after ##FASTA, which is what PPanGGOLiN needs.
for SET in c3 chr1; do
  : > "$W/ppanggolin/list_$SET.tsv"
  for f in "$W"/pangenome/gff_$SET/*.gff3; do
    n=$(basename "$f" .gff3)
    printf '%s\t%s\n' "$n" "$f" >> "$W/ppanggolin/list_$SET.tsv"
  done
  echo "$SET: $(wc -l < "$W/ppanggolin/list_$SET.tsv") genomes"
done

run_one () {                       # $1 = set (c3|chr1)   $2 = identity
  local SET=$1 ID=$2
  local O=$W/ppanggolin/${SET}_id${ID/./}
  local H=$O/pangenome.h5
  mkdir -p "$O"
  echo "=== START $SET identity=$ID -> $O ==="
  $PPG annotate --anno "$W/ppanggolin/list_$SET.tsv" -o "$O" -f -c "$CPU"       || return 1
  $PPG cluster   -p "$H" --identity "$ID" --coverage 0.8 --mode 1 -c "$CPU"     || return 1
  $PPG graph     -p "$H"                                                        || return 1  # no --cpu flag on 2.2.6
  $PPG partition -p "$H" -c "$CPU"                                              || return 1
  $PPG write_pangenome -p "$H" -o "$O/out" -f --Rtab --csv --stats --partitions || return 1
  echo "=== DONE $SET id=$ID ==="
}

for SET in c3 chr1; do
  for ID in 0.60 0.80; do
    run_one "$SET" "$ID" 2>&1 | tee -a "$W/logs/s13_${SET}_${ID/./}.log"
  done
done
echo "S13_ALL_DONE"
exit
