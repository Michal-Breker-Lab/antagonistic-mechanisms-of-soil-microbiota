#!/usr/bin/env bash
set -euo pipefail

OUT=$1; shift
IFS=',' read -r -a QS <<< "$1"; shift

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

{
  printf 'clone\tplatform\tcontig\tlength_bp'
  for q in "${QS[@]}"; do
      printf '\tnumreads_q%s\tcovbases_q%s\tbreadth_pct_q%s\tmeandepth_q%s\tmeanmapq_q%s' \
             "$q" "$q" "$q" "$q" "$q"
  done
  printf '\n'

  for spec in "$@"; do
      plat=${spec%%:*}
      bam=${spec#*:}
      clone=$(basename "$bam" .bam)

      for q in "${QS[@]}"; do
          samtools coverage -q "$q" "$bam" | awk 'NR>1' > "$tmp/q$q"
      done

      cut -f1,3 "$tmp/q${QS[0]}" | while IFS=$'\t' read -r contig len; do
          printf '%s\t%s\t%s\t%s' "$clone" "$plat" "$contig" "$len"
          for q in "${QS[@]}"; do
              awk -F'\t' -v c="$contig" -v q="$q" '
                  $1 == c { printf "\t%s\t%s\t%s\t%s\t%s", $4, $5, $6, $7, $9; found = 1 }
                  END { if (!found) { print "missing in q" q ": " c > "/dev/stderr"; exit 1 } }
              ' "$tmp/q$q"
          done
          printf '\n'
      done
  done
} > "$OUT"

echo "wrote $OUT ($(($(wc -l < "$OUT") - 1)) rows)" >&2
