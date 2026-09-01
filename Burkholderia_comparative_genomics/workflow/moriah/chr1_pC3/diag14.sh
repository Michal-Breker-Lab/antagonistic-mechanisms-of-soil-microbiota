#!/bin/bash
W=/sci/backup/ofinkel/moshea/burkholderia_c3
R=$(cat $W/logs/retry_jobid.txt)
echo "retry job $R"
for acc in GCF_000755825.1 GCF_000959345.1 GCF_009578005.1; do
  idx=$(grep -n "^${acc}	" $W/genus_species.tsv | cut -d: -f1)
  echo "=== $acc (task $idx) ==="
  echo "  sacct: $(sacct -j ${R}_${idx} --format=State,ExitCode,Elapsed -n -P 2>/dev/null | head -1)"
  echo "  stdout: $(tail -4 $W/logs/bakta_${R}_${idx}.out 2>/dev/null | tr '\n' ' | ')"
  echo "  stderr: $(tail -3 $W/logs/bakta_${R}_${idx}.err 2>/dev/null | tr '\n' ' ' | cut -c1-220)"
  echo "  fasta:  $(ls -la $W/genomes/${acc}.fna 2>/dev/null | awk '{print $5" bytes"}' || echo MISSING)"
done
