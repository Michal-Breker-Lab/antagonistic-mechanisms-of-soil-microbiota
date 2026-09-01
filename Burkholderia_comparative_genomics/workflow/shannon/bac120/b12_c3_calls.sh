#!/bin/bash
# pC3 presence/absence for the 142 bac120 genomes the original pipeline never
# typed, so the Figure 12b pC3 column has a measured value for every tip.
#
# METHOD IS s15b's, NOT A NEW ONE. s15b_assign_replicons.py faced exactly this
# situation for 3 genomes and its reasoning applies unchanged to 142:
#     chr1 := the longest contig
#     c3   := a secondary contig >= 300 kb, not chr1, whose best skani hit against
#             the KNOWN c3 contigs clears ANI >= 95% and aligned fraction >= 60%
#     anything failing that bar is UNASSIGNED and reported, never forced
# Identifying c3 by similarity to already-typed c3 replicons leaves the frozen
# classifier untouched. RE-TRAINING the c3 classifier on a larger set would
# perturb every published c3 call in the report, which is not acceptable for a
# figure fix.
#
# ONE DELIBERATE DEPARTURE, stated because it changes the reference set: s15b
# used 6 hand-picked c3 contigs, chosen as relatives of its 3 targets. 142
# genomes span far more of the genus than 3 near-clones do, so the reference is
# ALL 140 contigs marked is_c3=True in secondary_replicon_clusters.tsv. A
# 6-contig reference would under-call c3 in lineages distant from those 6, and
# an under-call would read on the figure as a confident absence.
set -euo pipefail
C=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
SK=$C/envs/burk/bin/skani
WORK=$W/c3_assign
mkdir -p "$WORK" "$W/logs"
cd "$W"
CPUS=${CPUS:-32}

echo "=== preflight ==="
"$SK" --version 2>&1 | head -1
[ -s need_prot.txt ] || { echo "FAIL: no need_prot.txt (run b10 first)" >&2; exit 1; }
echo "  genomes to type: $(wc -l < need_prot.txt)"

python3 - <<'PY'
import csv, os, subprocess, sys
from pathlib import Path
C = Path("/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
W = Path("/mnt/LargeStorageNoBackup/Moshea/burk_bac120")
WORK = W / "c3_assign"
SK = str(C / "envs" / "burk" / "bin" / "skani")
MIN_SECONDARY, MIN_ANI, MIN_AF = 300_000, 95.0, 60.0


def read_fasta(p):
    recs, name, buf = [], None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    recs.append((name, buf))
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        recs.append((name, buf))
    return recs


# ---- reference: every contig already typed as c3 ----------------------------
refdir = WORK / "ref"; refdir.mkdir(parents=True, exist_ok=True)
reflist = WORK / "refs.txt"
if not reflist.exists() or reflist.stat().st_size == 0:
    want = {}
    with open(C / "results" / "secondary_replicon_clusters.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["is_c3"] == "True":
                want.setdefault(r["accession"], set()).add(r["contig"])
    n = 0
    with open(reflist, "w") as out:
        for acc, contigs in want.items():
            fna = W / "genomes" / f"{acc}.fna"
            if not fna.exists():
                continue
            for name, buf in read_fasta(fna):
                if name in contigs:
                    p = refdir / f"ref__{acc}__{name}.fna"
                    if not p.exists():
                        p.write_text(f">{name}\n" + "\n".join(buf) + "\n")
                    out.write(str(p) + "\n"); n += 1
    print(f"  reference c3 contigs written: {n}", flush=True)
print(f"  reference list: {sum(1 for _ in open(reflist))} contigs", flush=True)

# ---- split each target's secondary contigs ----------------------------------
targets = [l.strip() for l in open(W / "need_prot.txt") if l.strip()]
qdir = WORK / "q"; qdir.mkdir(exist_ok=True)
info = {}                       # query path -> (acc, contig, length)
qlist = WORK / "queries.txt"
with open(qlist, "w") as out:
    for acc in targets:
        recs = read_fasta(W / "genomes" / f"{acc}.fna")
        lens = sorted(((n, sum(len(x) for x in b)) for n, b in recs),
                      key=lambda t: -t[1])
        if not lens:
            continue
        chr1 = lens[0][0]
        for name, L in lens[1:]:
            if L < MIN_SECONDARY:
                continue
            p = qdir / f"{acc}__{name}.fna"
            if not p.exists():
                seq = dict((n, b) for n, b in recs)[name]
                p.write_text(f">{name}\n" + "\n".join(seq) + "\n")
            info[str(p)] = (acc, name, L)
            out.write(str(p) + "\n")
print(f"  secondary contigs >={MIN_SECONDARY//1000} kb to test: {len(info)}", flush=True)

# ---- skani: every candidate secondary contig vs every known c3 --------------
out_tsv = WORK / "dist.tsv"
if not out_tsv.exists() or out_tsv.stat().st_size == 0:
    subprocess.run([SK, "dist", "--ql", str(qlist), "--rl", str(reflist),
                    "-t", os.environ.get("CPUS", "32"), "-o", str(out_tsv)],
                   check=True)
best = {}
with open(out_tsv) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        q = r["Query_file"]
        if q not in info:
            continue
        ani, af = float(r["ANI"]), float(r["Align_fraction_query"])
        cur = best.get(q)
        if cur is None or ani > cur[0]:
            best[q] = (ani, af, os.path.basename(r["Ref_file"]))

# ---- call ------------------------------------------------------------------
rows, npos = [], 0
for acc in targets:
    cands = [(p, v) for p, v in info.items() if v[0] == acc]
    hit = None
    for p, (_, contig, L) in cands:
        b = best.get(p)
        if b and b[0] >= MIN_ANI and b[1] >= MIN_AF:
            if hit is None or b[0] > hit[3]:
                hit = (contig, L, b[2], b[0], b[1])
    if hit:
        npos += 1
        rows.append([acc, "True", "measured_skani_vs_known_c3",
                     hit[0], hit[1], f"{hit[3]:.2f}", f"{hit[4]:.1f}", hit[2]])
    else:
        rows.append([acc, "False", "measured_skani_vs_known_c3",
                     "", "", "", "", ""])
out = W / "c3_calls_new.tsv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "c3_present", "evidence", "c3_contig", "c3_len",
                "c3_ani_to_ref", "c3_aligned_fraction", "c3_best_ref"])
    w.writerows(rows)
print(f"  typed {len(rows)} genomes: {npos} pC3-positive, {len(rows)-npos} negative")
print(f"  wrote {out}")
PY
echo "=== C3 CALLS OK ==="
