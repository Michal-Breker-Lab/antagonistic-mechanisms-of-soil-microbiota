#!/bin/bash
W=/sci/backup/ofinkel/moshea/burkholderia_c3
echo "gff3 present: $(ls $W/annot/*/*.gff3 2>/dev/null | wc -l) / 773"
echo "still missing:"
n=0
while read -r acc; do
  [ -z "$acc" ] && continue
  if [ ! -s "$W/annot/$acc/$acc.gff3" ]; then echo "  $acc"; n=$((n+1)); fi
done < "$W/genome_list_full.txt"
echo "  total: $n"
echo "retry queue: running=$(squeue -j $(cat $W/logs/retry_jobid.txt) -h -t RUNNING 2>/dev/null | wc -l) pending=$(squeue -j $(cat $W/logs/retry_jobid.txt) -h -t PENDING 2>/dev/null | wc -l)"
