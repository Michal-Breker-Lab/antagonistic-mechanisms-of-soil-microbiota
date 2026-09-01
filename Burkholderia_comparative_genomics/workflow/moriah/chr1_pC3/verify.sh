#!/bin/bash
W=/sci/backup/ofinkel/moshea/burkholderia_c3
low=0; high=0; empty=0; n=0
: > $W/results/annotation_cds_counts.tsv
while read -r acc; do
  [ -z "$acc" ] && continue
  f="$W/annot/$acc/$acc.gff3"; n=$((n+1))
  if [ ! -s "$f" ]; then echo "EMPTY $acc"; empty=$((empty+1)); continue; fi
  c=$(grep -cP "\tCDS\t" "$f"); t=$(grep -cP "\ttRNA\t" "$f"); s=$(grep -c "^##sequence-region" "$f")
  printf "%s\t%s\t%s\t%s\n" "$acc" "$c" "$t" "$s" >> $W/results/annotation_cds_counts.tsv
  [ "$c" -lt 2000 ] && { echo "LOW  $acc CDS=$c"; low=$((low+1)); }
  [ "$c" -gt 12000 ] && { echo "HIGH $acc CDS=$c"; high=$((high+1)); }
done < "$W/genome_list_full.txt"
echo "checked=$n empty=$empty low=$low high=$high"
echo
echo "CDS distribution:"
awk -F'\t' '{print $2}' $W/results/annotation_cds_counts.tsv | sort -n | awk '
  {a[NR]=$1; s+=$1}
  END{printf "  min=%d  q1=%d  median=%d  q3=%d  max=%d  mean=%.0f\n",
      a[1], a[int(NR*0.25)], a[int(NR*0.5)], a[int(NR*0.75)], a[NR], s/NR}'
echo
echo "lab isolates:"
grep -E "^(MF6|MF7)\b" $W/results/annotation_cds_counts.tsv | awk -F'\t' '{printf "  %-5s CDS=%-6s tRNA=%-4s contigs=%s\n",$1,$2,$3,$4}'
