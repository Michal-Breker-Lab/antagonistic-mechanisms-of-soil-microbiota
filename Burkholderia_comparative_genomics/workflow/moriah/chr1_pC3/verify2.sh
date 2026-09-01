#!/bin/bash
# Single-pass verification. The earlier loop spawned 3 greps per genome across a
# networked filesystem; one awk pass over each file is far cheaper, and writing
# to a fresh file avoids the interleaving that corrupted the first attempt.
W=/sci/backup/ofinkel/moshea/burkholderia_c3
OUT=$W/results/annotation_qc.tsv
printf "accession\tCDS\ttRNA\trRNA\tcontigs\n" > $OUT
while read -r acc; do
  [ -z "$acc" ] && continue
  f="$W/annot/$acc/$acc.gff3"
  [ -s "$f" ] || { printf "%s\tEMPTY\t0\t0\t0\n" "$acc" >> $OUT; continue; }
  awk -F'\t' -v a="$acc" '
    /^##sequence-region/ {s++; next}
    /^#/ {next}
    $3=="CDS"  {c++}
    $3=="tRNA" {t++}
    $3=="rRNA" {r++}
    END {printf "%s\t%d\t%d\t%d\t%d\n", a, c, t, r, s}' "$f" >> $OUT
done < "$W/genome_list_full.txt"
echo "rows written: $(( $(wc -l < $OUT) - 1 )) (expect 773)"
awk -F'\t' 'NR>1{
  if($2=="EMPTY"){e++; next}
  n++; s+=$2; v[n]=$2
  if($2<2000) low++; if($2>12000) high++
} END{
  asort(v)
  printf "empty=%d  low(<2000)=%d  high(>12000)=%d\n", e, low, high
  printf "CDS: min=%d q1=%d median=%d q3=%d max=%d mean=%.0f\n",
    v[1], v[int(n*0.25)], v[int(n*0.5)], v[int(n*0.75)], v[n], s/n
}' $OUT
echo "lab isolates:"
grep -E "^(MF6|MF7)	" $OUT | awk -F'\t' '{printf "  %-5s CDS=%-6s tRNA=%-4s rRNA=%-3s contigs=%s\n",$1,$2,$3,$4,$5}'
