#!/usr/bin/env bash
# Stage 7 + 8 - two pangenomes from the same annotations:
#   chr1  -> strict   -> core alignment -> species tree   (no plasmids expected;
#                                                          phylogenetics is the goal)
#   c3    -> moderate -> c3 pangenome + core alignment    (strict's 5% floor would
#                                                          delete real c3 accessory
#                                                          families; see TOOLS.md)
#   c3    -> sensitive -> reported sensitivity check
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
PANAROO=$W/envs/panaroo/bin/panaroo
export PATH=$W/envs/panaroo/bin:$W/envs/burk/bin:$PATH
mkdir -p $W/pangenome/{gff_chr1,gff_c3} $W/logs

# ---------- build per-replicon GFF3 subsets ----------
echo "=== subsetting GFF3s ==="
python3 - <<'PY'
import csv, os, subprocess, collections
W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
types = list(csv.DictReader(open(f"{W}/results/replicon_types.tsv"), delimiter="\t"))
c3 = {}
for r in csv.DictReader(open(f"{W}/results/secondary_replicon_clusters.tsv"), delimiter="\t"):
    if r.get("is_c3") == "True":
        c3.setdefault(r["accession"], []).append(r["contig"])
chr1 = collections.defaultdict(list)
for r in types:
    if r["replicon_type"].startswith("chromosome1"):
        chr1[r["accession"]].append(r["contig"])

sub = f"{W}/subset_gff.py"
for label, mapping, outdir in (("chr1", chr1, f"{W}/pangenome/gff_chr1"),
                               ("c3", c3, f"{W}/pangenome/gff_c3")):
    n = 0
    for acc, ctgs in mapping.items():
        src = f"{W}/annot/{acc}/{acc}.gff3"
        if not os.path.exists(src) or not ctgs:
            continue
        subprocess.run(["python3", sub, src, f"{outdir}/{acc}.gff3", *ctgs],
                       stdout=subprocess.DEVNULL)
        n += 1
    print(f"{label}: wrote {n} GFF3 subsets")
PY

ls $W/pangenome/gff_chr1/*.gff3 > $W/pangenome/list_chr1.txt 2>/dev/null
ls $W/pangenome/gff_c3/*.gff3   > $W/pangenome/list_c3.txt   2>/dev/null
echo "chr1 GFFs: $(wc -l < $W/pangenome/list_chr1.txt)"
echo "c3   GFFs: $(wc -l < $W/pangenome/list_c3.txt)"

run_panaroo () {  # $1=list $2=outdir $3=mode $4=extra
  local list=$1 out=$2 mode=$3; shift 3
  echo "=== panaroo $mode -> $out ==="
  date
  $PANAROO -i "$list" -o "$out" --clean-mode "$mode" -t 48 \
           --remove-invalid-genes "$@" > $W/logs/panaroo_$(basename $out).log 2>&1
  local rc=$?
  echo "exit=$rc"
  [ -f "$out/summary_statistics.txt" ] && cat "$out/summary_statistics.txt"
  date
}

# Stage 8: chromosome-1 core, aligned, for the species tree
run_panaroo $W/pangenome/list_chr1.txt $W/pangenome/chr1_strict strict \
            -a core --aligner mafft --core_threshold 0.95

# Stage 7: c3 pangenome (primary) + sensitivity
run_panaroo $W/pangenome/list_c3.txt $W/pangenome/c3_moderate moderate \
            -a core --aligner mafft --core_threshold 0.95
run_panaroo $W/pangenome/list_c3.txt $W/pangenome/c3_sensitive sensitive

echo "PANGENOME_DONE"
