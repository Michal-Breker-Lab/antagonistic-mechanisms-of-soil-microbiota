#!/bin/bash
# Exercise the exact per-task parsing the array job uses.
W=/sci/backup/ofinkel/moshea/burkholderia_c3
for SLURM_ARRAY_TASK_ID in 1 2 400 771 772 773; do
  read -r acc genus sp < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$W/genus_species.tsv" | tr '\t' ' ')
  case "$acc" in
    MF6) fna="$W/MF6.fna"; cflag="--complete" ;;
    MF7) fna="$W/MF7.fna"; cflag="" ;;
    *)   fna="$W/genomes/$acc.fna"; cflag="--complete" ;;
  esac
  printf "task %-4s acc=%-18s genus=%-16s sp=%-12s complete=%-11s fasta=%s\n" \
     "$SLURM_ARRAY_TASK_ID" "$acc" "$genus" "$sp" "${cflag:-NO}" "$([ -s "$fna" ] && echo present || echo 'not yet)')"
done
echo
echo "blank/short lines in genus_species.tsv: $(awk -F'\t' 'NF!=3' "$W/genus_species.tsv" | wc -l)"
echo "duplicate accessions: $(cut -f1 "$W/genus_species.tsv" | sort | uniq -d | wc -l)"
