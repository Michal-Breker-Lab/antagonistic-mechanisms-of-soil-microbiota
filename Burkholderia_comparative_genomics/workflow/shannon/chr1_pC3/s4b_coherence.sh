#!/usr/bin/env bash
# Stage 4b - build ortholog groups across ALL large secondary replicons, so the
# question "is c3 one evolutionary entity or just a size class?" can be answered
# by shared gene content.
#
# Gene content rather than ANI, because skani is documented reliable only above
# ~82% ANI with >=15% aligned fraction, and cross-genus c3 comparisons fall far
# below that floor (see TOOLS.md).
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
MMSEQS=$W/envs/burk/bin/mmseqs
mkdir -p $W/replicons $W/results

echo "=== collecting proteins from large secondary replicons ==="
python3 - <<'PY'
import collections, os
W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
want = collections.defaultdict(dict)          # accession -> locus_tag -> contig
with open(f"{W}/results/secondary_replicon_loci.tsv") as fh:
    next(fh)
    for line in fh:
        acc, ctg, lt = line.rstrip("\n").split("\t")
        want[acc][lt] = ctg

out = open(f"{W}/replicons/secondary_proteins.faa", "w")
n_written = n_missing = 0
for acc, lmap in want.items():
    faa = f"{W}/annot/{acc}/{acc}.faa"
    if not os.path.exists(faa):
        n_missing += 1
        continue
    keep, hdr = False, None
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                lt = line[1:].split()[0]
                ctg = lmap.get(lt)
                keep = ctg is not None
                if keep:
                    # header encodes replicon identity: accession|contig|locus_tag
                    out.write(f">{acc}|{ctg}|{lt}\n")
                    n_written += 1
            elif keep:
                out.write(line)
out.close()
print(f"proteins written: {n_written}")
print(f"genomes with no .faa: {n_missing}")
PY

echo "=== mmseqs clustering (50% identity, 80% coverage) ==="
rm -rf $W/replicons/mmseqs_tmp $W/replicons/clu*
$MMSEQS easy-cluster $W/replicons/secondary_proteins.faa \
    $W/replicons/clu $W/replicons/mmseqs_tmp \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 48 -v 1

echo "clusters: $(cut -f1 $W/replicons/clu_cluster.tsv | sort -u | wc -l)"
echo "members : $(wc -l < $W/replicons/clu_cluster.tsv)"
echo "S4B_CLUSTER_DONE"
