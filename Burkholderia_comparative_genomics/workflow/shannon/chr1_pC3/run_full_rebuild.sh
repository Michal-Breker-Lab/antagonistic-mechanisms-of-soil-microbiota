#!/usr/bin/env bash
# Driver for the non-dereplicated pC3 rebuild on Shannon.
#
# Stages are gated by .done markers under $W/state, so a re-run resumes rather
# than repeats. Run one stage at a time with --stage N, or --all to chain them.
# Every stage verifies its own output and exits non-zero on failure; --all stops
# at the first failure rather than carrying a broken input forward.
#
# Stage 5 (pangenomes) and 6 (trees) are the multi-day steps. Nothing here
# submits to a scheduler -- this is Shannon, not Moriah.
set -uo pipefail

W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
export TMPDIR="$W/tmp"            # Shannon's / is ~100% full; never use /tmp
mkdir -p "$TMPDIR" "$W/logs" "$W/state" "$W/results"

THREADS=${THREADS:-64}
PANAROO_CEILING_H=${PANAROO_CEILING_H:-72}
STAGE=""; ALL=0
while [ $# -gt 0 ]; do
  case "$1" in
    --stage) STAGE="$2"; shift 2;;
    --all) ALL=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

log()  { echo "[$(date +%F' '%T)] $*" | tee -a "$W/logs/rebuild.log"; }
done_marker() { echo "$W/state/stage$1.done"; }
is_done() { [ -f "$(done_marker "$1")" ]; }
mark()    { touch "$(done_marker "$1")"; }

# ------------------------------------------------------------------ stage 0
stage0() {
  log "STAGE 0 preflight"
  local fail=0
  local avail
  avail=$(df -BG --output=avail /mnt/LargeStorageNoBackup | tail -1 | tr -dc '0-9')
  log "  disk free on LargeStorageNoBackup: ${avail} GB"
  # Bakta output for ~464 genomes plus Panaroo intermediates at 771 genomes.
  if [ "$avail" -lt 250 ]; then
    log "  FAIL: need >=250 GB free, have ${avail} GB"; fail=1
  fi
  local la
  la=$(awk '{print int($1)}' /proc/loadavg)
  log "  load average: ${la} (cores: $(nproc))"
  if [ "$la" -gt 200 ]; then
    log "  FAIL: Shannon is saturated - do not start"; fail=1
  fi
  log "  TMPDIR=$TMPDIR"
  [ -d "$W/genomes" ] || { log "  FAIL: $W/genomes missing"; fail=1; }
  [ -s "$W/MF6.fna" ] || { log "  FAIL: MF6.fna missing"; fail=1; }
  [ -s "$W/MF7.fna" ] || { log "  FAIL: MF7.fna missing - stage it first"; fail=1; }
  [ "$fail" -eq 0 ] || return 1
  mark 0; log "STAGE 0 OK"
}

# ------------------------------------------------------------------ stage 1
stage1() {
  log "STAGE 1 refresh genome set, build genome_list_full.txt"
  # s1_metadata.py overwrites metadata/ in place. Keep the 2026-08-03 snapshot so
  # "did the set actually change" stays answerable after the refresh.
  if [ ! -d "$W/metadata_20260803" ]; then
    cp -a "$W/metadata" "$W/metadata_20260803"
    log "  backed up original metadata/ -> metadata_20260803/"
  fi
  python3 "$W/s1_metadata.py" 2>&1 | tee -a "$W/logs/s1_full.log"
  # s1 writes accessions.txt (the GCA/GCF-deduplicated list), assemblies.tsv,
  # replicons.tsv and dataset_reports.json -- there is no genome_set.tsv.
  local src="$W/metadata/accessions.txt"
  [ -s "$src" ] || { log "  FAIL: $src missing"; return 1; }
  sort -u "$src" > "$W/genome_list_full.txt"
  printf "MF6\nMF7\n" >> "$W/genome_list_full.txt"
  local n; n=$(wc -l < "$W/genome_list_full.txt")
  log "  full set size: $n"
  if [ "$n" -lt 760 ] || [ "$n" -gt 820 ]; then
    log "  FAIL: unexpected genome count $n - investigate before continuing"; return 1
  fi
  if [ "$(sort "$W/genome_list_full.txt" | uniq -d | wc -l)" -ne 0 ]; then
    log "  FAIL: duplicate accessions in the list"; return 1
  fi
  # what actually changed vs the 2026-08-03 snapshot
  if [ -s "$W/metadata_20260803/accessions.txt" ]; then
    log "  new since 2026-08-03: $(comm -13 <(sort -u "$W/metadata_20260803/accessions.txt") \
                                            <(sort -u "$W/metadata/accessions.txt") | tr '\n' ' ')"
    log "  gone since 2026-08-03: $(comm -23 <(sort -u "$W/metadata_20260803/accessions.txt") \
                                            <(sort -u "$W/metadata/accessions.txt") | tr '\n' ' ')"
  fi
  mark 1; log "STAGE 1 OK"
}

# ------------------------------------------------------------------ stage 2
stage2() {
  log "STAGE 2 annotate the remainder (Bakta)"
  bash "$W/s6b_annotate_rest.sh" --confirm 2>&1 | tee -a "$W/logs/s6b.log"
  local rc=${PIPESTATUS[0]}
  [ "$rc" -eq 0 ] || { log "  FAIL: annotation stage returned $rc"; return 1; }
  mark 2; log "STAGE 2 OK"
}

# ------------------------------------------------------------------ stage 3
stage3() {
  log "STAGE 3 replicon census + pC3 calls on the full set"
  python3 "$W/s3_census.py"   2>&1 | tee -a "$W/logs/s3_full.log"
  python3 "$W/s4a_typing.py"  2>&1 | tee -a "$W/logs/s4a_full.log"
  python3 "$W/s12_extend_c3.py" 2>&1 | tee -a "$W/logs/s12_full.log"
  # MF7 is a draft: place its contigs against MF6 rather than by size rank.
  # blastn, not minimap2 -- minimap2 is not installed in any project env, while
  # blastn is (envs/panaroo). The -outfmt below is chosen so its columns line up
  # with PAF's first eleven, so s30 parses either without a code path of its own.
  if [ -s "$W/MF7.fna" ]; then
    BL=$W/envs/panaroo/bin
    [ -s "$W/results/mf6_replicon_map.tsv" ] || {
      log "  FAIL: results/mf6_replicon_map.tsv missing"; return 1; }
    "$BL/makeblastdb" -in "$W/MF6.fna" -dbtype nucl \
       -out "$W/tmp/mf6db" > "$W/logs/mf6_makeblastdb.log" 2>&1 || {
       log "  FAIL: makeblastdb"; return 1; }
    "$BL/blastn" -query "$W/MF7.fna" -db "$W/tmp/mf6db" \
       -outfmt "6 qseqid qlen qstart qend sstrand sseqid slen sstart send nident length" \
       -evalue 1e-20 -perc_identity 90 -max_target_seqs 20 \
       -num_threads "$THREADS" > "$W/results/mf7_vs_mf6.tsv" \
       2> "$W/logs/mf7_blastn.log" || { log "  FAIL: blastn"; return 1; }
    log "  blastn HSPs: $(wc -l < "$W/results/mf7_vs_mf6.tsv")"
    python3 "$W/s30_mf7_replicon_assign.py" \
      --paf "$W/results/mf7_vs_mf6.tsv" \
      --replicon-map "$W/results/mf6_replicon_map.tsv" \
      --out "$W/results/mf7_replicon_assignment.tsv" \
      2>&1 | tee -a "$W/logs/s30.log"
  fi
  mark 3; log "STAGE 3 OK"
}

# ------------------------------------------------------------------ stage 4
stage4() {
  log "STAGE 4 clone-cluster labels"
  # ANI of the post-dereplication genomes (incl. MF6/MF7) to the labelled set,
  # so they join an existing 99% cluster instead of becoming false singletons.
  if [ ! -s "$W/results/ani_newcomers.tsv" ]; then
    log "  computing skani dist for newcomers"
    "$W/envs/burk/bin/skani" dist -t "$THREADS" \
        -q "$W/MF6.fna" "$W/MF7.fna" -r "$W"/genomes/*.fna \
        > "$W/results/ani_newcomers.tsv" 2> "$W/logs/skani_newcomers.log"
  fi
  python3 "$W/s28_clone_labels.py" \
    --membership "$W/results/derep_cluster_membership.tsv" \
    --genomes    "$W/genome_list_full.txt" \
    --ani        "$W/results/ani_newcomers.tsv" \
    --out        "$W/results/clone_cluster.tsv" 2>&1 | tee -a "$W/logs/s28.log"
  [ -s "$W/results/clone_cluster.tsv" ] || { log "  FAIL: no clone_cluster.tsv"; return 1; }
  mark 4; log "STAGE 4 OK"
}

# ------------------------------------------------------------------ stage 5
stage5() {
  log "STAGE 5 pangenomes at full size (ceiling ${PANAROO_CEILING_H} h)"
  # Expected direction: core RISES and cloud fraction FALLS relative to the
  # 306-genome run, because most added genomes are near-clones. A falling core
  # means something is wrong -- stop rather than build trees on it.
  timeout "${PANAROO_CEILING_H}h" bash "$W/s78_pangenomes.sh" full \
    2>&1 | tee -a "$W/logs/panaroo_full.log"
  local rc=${PIPESTATUS[0]}
  if [ "$rc" -eq 124 ]; then
    log "  Panaroo hit the ${PANAROO_CEILING_H} h ceiling."
    log "  FALLBACK REQUIRED: PPanGGOLiN (s13_ppanggolin.sh), and the 306-genome"
    log "  set must be re-run through PPanGGOLiN too so old-vs-new is tool-matched."
    log "  This is a scientific decision point - stop and report."
    return 1
  fi
  [ "$rc" -eq 0 ] || { log "  FAIL: panaroo returned $rc"; return 1; }
  mark 5; log "STAGE 5 OK"
}

# ------------------------------------------------------------------ stage 6
stage6() {
  log "STAGE 6 trees (model fixed to GTR+F+R10, no ModelFinder)"
  # The 304-tip run selected GTR+F+R10 at 836 CPU-h / 18 h wall. Fixing the
  # model keeps the new tree comparable and removes ModelFinder's share of a
  # job that is already 2,500-4,200 CPU-h at 771 taxa.
  bash "$W/s89_trees.sh" full 2>&1 | tee -a "$W/logs/trees_full.log"
  local rc=${PIPESTATUS[0]}
  [ "$rc" -eq 0 ] || { log "  FAIL: tree stage returned $rc"; return 1; }
  for t in chr1_core_full c3_core_full; do
    [ -s "$W/trees/$t.treefile" ] || { log "  FAIL: $t.treefile missing"; return 1; }
    local n; n=$(tr ',' '\n' < "$W/trees/$t.treefile" | grep -c ":")
    log "  $t: ~$n tips"
  done
  grep -q "MF7" "$W/trees/chr1_core_full.treefile" || log "  WARN: MF7 not on chr1 tree"
  grep -q "MF6" "$W/trees/chr1_core_full.treefile" || log "  WARN: MF6 not on chr1 tree"
  mark 6; log "STAGE 6 OK"
}

# ------------------------------------------------------------------ stage 7
stage7() {
  log "STAGE 7 searches on the full set (toxins, warhead, effectors)"
  python3 "$W/s23_pyrodigal.py"     2>&1 | tee -a "$W/logs/s23_full.log"
  python3 "$W/s24_search.py"        2>&1 | tee -a "$W/logs/s24_full.log"
  python3 "$W/s27_search_toxins.py" 2>&1 | tee -a "$W/logs/s27_full.log"
  python3 "$W/s20_screen.py"        2>&1 | tee -a "$W/logs/s20_full.log"
  mark 7; log "STAGE 7 OK"
}

# ------------------------------------------------------------------ stage 8
stage8() {
  log "STAGE 8 downstream stats with dual nulls"
  python3 "$W/s10_hosts.py" 2>&1 | tee -a "$W/logs/s10_full.log"
  Rscript  "$W/s11_phylo_stats.R" 2>&1 | tee -a "$W/logs/s11_full.log"
  # Both nulls for every randomisation test. The stratified one is the
  # interpretable test; the unrestricted one is kept for comparability.
  rm -f "$W/results/null_dual.tsv"
  for set in warhead toxin12; do
    [ -s "$W/results/carriers_${set}.txt" ] || { log "  WARN: carriers_${set}.txt missing"; continue; }
    python3 "$W/s29_null_stratified.py" \
      --tree "$W/trees/chr1_core_full.treefile" \
      --clusters "$W/results/clone_cluster.tsv" \
      --carriers "$W/results/carriers_${set}.txt" \
      --label "${set}_alltips" \
      --out "$W/results/null_dual.tsv" 2>&1 | tee -a "$W/logs/s29.log"
    if [ -s "$W/results/pc3_tips.txt" ]; then
      python3 "$W/s29_null_stratified.py" \
        --tree "$W/trees/chr1_core_full.treefile" \
        --clusters "$W/results/clone_cluster.tsv" \
        --carriers "$W/results/carriers_${set}.txt" \
        --pool "$W/results/pc3_tips.txt" \
        --label "${set}_pc3only" \
        --out "$W/results/null_dual.tsv" 2>&1 | tee -a "$W/logs/s29.log"
    fi
  done
  mark 8; log "STAGE 8 OK"
}

run_stage() {
  local s="$1"
  if is_done "$s"; then log "stage $s already done (rm $(done_marker "$s") to force)"; return 0; fi
  "stage$s"
}

if [ -n "$STAGE" ]; then
  run_stage "$STAGE"; exit $?
elif [ "$ALL" -eq 1 ]; then
  for s in 0 1 2 3 4 5 6 7 8; do
    run_stage "$s" || { log "STOPPED at stage $s"; exit 1; }
  done
  log "ALL STAGES COMPLETE"
  exit 0
else
  echo "usage: $0 --stage N | --all"
  echo "stages: 0 preflight 1 genome-set 2 annotate 3 census 4 clone-labels"
  echo "        5 pangenomes 6 trees 7 searches 8 stats"
  exit 2
fi
