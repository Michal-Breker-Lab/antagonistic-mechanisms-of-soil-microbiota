#!/usr/bin/env python3
"""Stage 4a - type every replicon from its Bakta annotation.

Chromosome 1 is identified by CONTENT, not by size (decision D2): it carries
dnaA and the ribosomal-protein superoperon. This matters for two reasons:

  * a genome could in principle have a chromosome 2 larger than chromosome 1;
  * fused assemblies (chr1+chr2 collapsed into one 6-7 Mb molecule) are only
    detectable this way, and Stage 3 found 40 candidates. A fused replicon
    carries chromosome-1 markers AND is far larger than any real chr1, so it
    is flagged rather than silently counted as a c3 loss.

Outputs one row per replicon with its markers, plus a protein FASTA per
secondary replicon for the Stage 4b gene-content clustering.
"""
import collections
import csv
import glob
import os
import re
import sys

W = os.environ.get("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
RES, ANN = f"{W}/results", f"{W}/annot"
MIN_LARGE = 300_000
FUSED_MAX_CHR1 = 4_500_000     # p95 of largest-replicon size across the genus

RIBO = re.compile(r"^rp[sml][A-Z]\d*$", re.I)
CORE_MARKERS = {"dnaA", "gyrB", "rpoB", "rpoC", "recA", "ftsZ", "dnaN"}


def parse_gff(path):
    """-> contig -> {length, cds, ribo, markers:set, loci:[...]}"""
    info = collections.defaultdict(
        lambda: {"length": 0, "cds": 0, "ribo": 0, "markers": set(), "loci": []})
    with open(path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    info[p[1]]["length"] = int(p[3])
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            ctg = f[0]
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            d = info[ctg]
            d["cds"] += 1
            lt = attrs.get("locus_tag", "")
            if lt:
                d["loci"].append(lt)
            gene = attrs.get("gene", "")
            if gene:
                if RIBO.match(gene):
                    d["ribo"] += 1
                if gene in CORE_MARKERS:
                    d["markers"].add(gene)
    return info


def main():
    gffs = sorted(glob.glob(f"{ANN}/*/*.gff3"))
    print(f"annotations found: {len(gffs)}")
    if not gffs:
        sys.exit("no GFF3 files - run Stage 6 first")

    rows, loci_map = [], {}
    for i, g in enumerate(gffs, 1):
        acc = os.path.basename(os.path.dirname(g))
        info = parse_gff(g)
        if not info:
            print(f"  WARN empty annotation: {acc}")
            continue
        # chromosome 1 = most ribosomal proteins; ties broken by dnaA then size
        def chr1_score(kv):
            c, d = kv
            return (d["ribo"], "dnaA" in d["markers"], d["length"])
        chr1 = max(info.items(), key=chr1_score)[0]

        order = sorted(info.items(), key=lambda kv: -kv[1]["length"])
        rank = {c: n + 1 for n, (c, _) in enumerate(order)}

        for ctg, d in order:
            is_chr1 = (ctg == chr1)
            fused = is_chr1 and d["length"] > FUSED_MAX_CHR1
            if is_chr1:
                rtype = "chromosome1_fused" if fused else "chromosome1"
            elif d["length"] >= MIN_LARGE:
                rtype = "secondary_large"
            else:
                rtype = "small_plasmid"
            rows.append({
                "accession": acc, "contig": ctg, "size_rank": rank[ctg],
                "length": d["length"], "cds": d["cds"], "ribosomal_proteins": d["ribo"],
                "core_markers": ";".join(sorted(d["markers"])),
                "n_core_markers": len(d["markers"]),
                "replicon_type": rtype,
            })
            if rtype == "secondary_large":
                loci_map[(acc, ctg)] = d["loci"]
        if i % 50 == 0:
            print(f"  parsed {i}/{len(gffs)}", flush=True)

    cols = list(rows[0].keys())
    os.makedirs(RES, exist_ok=True)
    with open(f"{RES}/replicon_types.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {RES}/replicon_types.tsv  ({len(rows)} replicons)")

    c = collections.Counter(r["replicon_type"] for r in rows)
    print("\n=== replicon types ===")
    for k, v in c.most_common():
        print(f"  {v:>5}  {k}")

    # sanity: does chr1 always look like chr1?
    chr1s = [r for r in rows if r["replicon_type"].startswith("chromosome1")]
    rib = sorted(r["ribosomal_proteins"] for r in chr1s)
    print(f"\nribosomal proteins on called chromosome 1: "
          f"min={rib[0]} median={rib[len(rib)//2]} max={rib[-1]}")
    bad = [r for r in chr1s if r["ribosomal_proteins"] < 20]
    print(f"chromosome-1 calls with <20 ribosomal proteins (suspect): {len(bad)}")
    for r in bad[:10]:
        print(f"   {r['accession']} {r['contig']} len={r['length']:,} ribo={r['ribosomal_proteins']}")

    # is chr1 ever NOT the largest replicon?
    notbig = [r for r in chr1s if r["size_rank"] != 1]
    print(f"\ngenomes where chromosome 1 is NOT the largest replicon: {len(notbig)}")
    for r in notbig[:10]:
        print(f"   {r['accession']} {r['contig']} rank={r['size_rank']} len={r['length']:,}")

    fused = [r for r in chr1s if r["replicon_type"] == "chromosome1_fused"]
    print(f"\nfused chromosome-1 candidates (>{FUSED_MAX_CHR1:,} bp): {len(fused)}")
    for r in sorted(fused, key=lambda x: -x["length"])[:12]:
        print(f"   {r['accession']} len={r['length']:,} ribo={r['ribosomal_proteins']} cds={r['cds']}")

    # per-genome count of large secondary replicons -> the annotation-based census
    per = collections.Counter()
    for r in rows:
        if r["replicon_type"] == "secondary_large":
            per[r["accession"]] += 1
    allacc = {r["accession"] for r in rows}
    dist = collections.Counter(per.get(a, 0) for a in allacc)
    print("\n=== large secondary replicons per genome (annotation-based) ===")
    for k in sorted(dist):
        print(f"  {k} secondary large: {dist[k]} genomes")

    with open(f"{RES}/secondary_replicon_loci.tsv", "w") as fh:
        fh.write("accession\tcontig\tlocus_tag\n")
        for (acc, ctg), loci in loci_map.items():
            for lt in loci:
                fh.write(f"{acc}\t{ctg}\t{lt}\n")
    print(f"\nsecondary-replicon CDS: {sum(len(v) for v in loci_map.values())} "
          f"across {len(loci_map)} replicons")
    print("STAGE4A_DONE")


if __name__ == "__main__":
    sys.exit(main())
