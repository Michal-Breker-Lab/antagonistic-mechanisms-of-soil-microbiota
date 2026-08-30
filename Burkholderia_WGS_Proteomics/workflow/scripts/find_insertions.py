"""Call Tn5 insertion sites per clone from junction and spanning-read evidence."""
import csv
import os
import sys
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tn_annotate import annotate, parse_gff  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")

CFG = snakemake.params.detect
PLATFORM = snakemake.params.platform
MIN_INS = CFG["min_ins"]
CLUSTER = CFG["cluster"]
MIN_AF = CFG["min_af"]
CASSETTE_LEN = CFG["cassette_len"]
CASSETTE = CFG["cassette_name"]
ABSENT = snakemake.params.absent  # noqa: F821
FLANK = CFG["presets"][PLATFORM]["FLANK"]
ADJACENT = CFG["presets"][PLATFORM]["ADJACENT"]
OVERLAP = CFG["presets"][PLATFORM]["OVERLAP"]
MIN_READS = CFG["presets"][PLATFORM]["MIN_READS"]




def read_key(aln):
    """Mate index keeps R1 and R2, which share a QNAME, as separate molecules."""
    return f"{aln.query_name}/{2 if aln.is_read2 else 1}"


def fwd_query(aln):
    """Query interval in original read orientation, comparable across strands."""
    ql = aln.infer_read_length()
    s, e = aln.query_alignment_start, aln.query_alignment_end
    return (ql - e, ql - s) if aln.is_reverse else (s, e)


def detect(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ev = []
    orient = {}

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_duplicate or aln.reference_name == CASSETTE:
            continue
        ref = aln.reference_start
        qpos = 0
        for op, n in aln.cigartuples or []:
            if op == 1:
                if n >= MIN_INS:
                    seq = aln.query_sequence
                    ev.append({
                        "contig": aln.reference_name, "pos": ref + 1,
                        "read": read_key(aln), "how": "spanning_I",
                        "ins_len": n, "side": "both",
                        "seq": seq[qpos:qpos + n] if seq else "",
                        "cassette_rev": None,
                    })
                qpos += n
            elif op in (0, 7, 8):
                ref += n; qpos += n
            elif op in (2, 3):
                ref += n
            elif op == 4:
                qpos += n

    cas = defaultdict(list)
    for aln in bam.fetch(CASSETTE):
        if aln.is_unmapped or aln.is_duplicate:
            continue
        q0, q1 = fwd_query(aln)
        cas[read_key(aln)].append((q0, q1, aln.is_reverse))
    if cas:
        bam2 = pysam.AlignmentFile(bam_path, "rb")
        for aln in bam2.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_duplicate or aln.reference_name == CASSETTE:
                continue
            hits = cas.get(read_key(aln))
            if not hits:
                continue
            g0, g1 = fwd_query(aln)
            for t0, t1, trev in hits:
                rev = trev != aln.is_reverse
                if -OVERLAP <= t0 - g1 <= ADJACENT:
                    at_end = not aln.is_reverse
                elif -OVERLAP <= g0 - t1 <= ADJACENT:
                    at_end = aln.is_reverse
                elif g0 <= t0 and t1 <= g1:
                    orient[read_key(aln)] = rev
                    continue
                else:
                    continue
                orient[read_key(aln)] = rev
                pos = aln.reference_end if at_end else aln.reference_start + 1
                ev.append({"contig": aln.reference_name, "pos": pos,
                           "read": read_key(aln), "how": "junction",
                           "ins_len": None, "side": "left" if at_end else "right",
                           "seq": "", "cassette_rev": rev})
        bam2.close()

    depth = {}
    cas_core = []
    for ref, ln in zip(bam.references, bam.lengths):
        tot = 0
        for col in bam.pileup(ref, min_mapping_quality=10, truncate=True):
            tot += col.nsegments
            if ref == CASSETTE and 50 <= col.reference_pos + 1 <= CASSETTE_LEN - 50:
                cas_core.append(col.nsegments)
        depth[ref] = (tot / ln, ln)
    bam.close()
    cas_core.sort()
    cas_med = cas_core[len(cas_core) // 2] if cas_core else 0
    return ev, depth, orient, cas_med


def cluster(ev):
    out = []
    by_contig = defaultdict(list)
    for e in ev:
        by_contig[e["contig"]].append(e)
    for contig, items in by_contig.items():
        items.sort(key=lambda d: d["pos"])
        cur = [items[0]]
        for e in items[1:]:
            if e["pos"] - cur[-1]["pos"] <= CLUSTER:
                cur.append(e)
            else:
                out.append((contig, cur)); cur = [e]
        out.append((contig, cur))
    return out


def ref_spanning(bam, contig, pos, support_reads):
    n = 0
    for aln in bam.fetch(contig, max(0, pos - FLANK - 1), pos + FLANK):
        if (aln.is_unmapped or aln.is_supplementary or aln.is_secondary
                or aln.is_duplicate or read_key(aln) in support_reads):
            continue
        if aln.reference_start > pos - FLANK or aln.reference_end < pos + FLANK:
            continue
        if any(op == 1 and ln >= MIN_INS for op, ln in aln.cigartuples or []):
            continue
        n += 1
    return n




def main():
    print(f"platform {PLATFORM}: FLANK={FLANK} ADJACENT={ADJACENT} "
          f"OVERLAP={OVERLAP} MIN_READS={MIN_READS}", file=sys.stderr)
    feats = parse_gff(snakemake.input.gff)
    rows, allrows, cov, summary = [], [], [], []
    fa = open(snakemake.output.inserted_seqs, "w")

    for path in sorted(snakemake.input.bams):
        clone = os.path.basename(path)[:-4]
        ev, depth, orient, cas_med = detect(path)
        for ref, (mean, ln) in depth.items():
            cov.append({"clone": clone, "reference": ref, "length_bp": ln,
                        "mean_depth": round(mean, 2)})
        bam = pysam.AlignmentFile(path, "rb")
        genome_depth = [m for r, (m, _) in depth.items() if r != CASSETTE]
        med_depth = sorted(genome_depth)[len(genome_depth) // 2] if genome_depth else 0

        n_kept = 0
        for contig, items in cluster(ev):
            reads = {e["read"] for e in items}
            lens = sorted(e["ins_len"] for e in items if e["ins_len"])
            lefts = sorted(e["pos"] for e in items if e["side"] == "left")
            rights = sorted(e["pos"] for e in items if e["side"] == "right")
            hows = sorted({e["how"] for e in items})
            lm = lefts[len(lefts) // 2] if lefts else None
            rm = rights[len(rights) // 2] if rights else None

            tsd_j = (lm - rm + 1) if (lm is not None and rm is not None) else ""
            med_ins = lens[len(lens) // 2] if lens else ""
            tsd_i = (med_ins - CASSETTE_LEN) if lens else ""

            if lm is not None and rm is not None:
                tsd_start, tsd_end = min(lm, rm), max(lm, rm)
                pos = tsd_start
            else:
                tsd_start = tsd_end = ""
                pos = int(sorted(e["pos"] for e in items)[len(items) // 2])
            revs = [orient[r] for r in reads if r in orient]
            refn = ref_spanning(bam, contig, pos, reads)
            af = len(reads) / (len(reads) + refn) if (len(reads) + refn) else 0.0
            ann = annotate(feats, contig, pos)
            row = {
                "clone": clone, "platform": PLATFORM,
                "contig": contig, "insertion_pos": pos,
                "tsd_start": tsd_start, "tsd_end": tsd_end,
                "tsd_len_junction": tsd_j, "tsd_len_insert": tsd_i,
                "left_junction": lm if lm is not None else "",
                "right_junction": rm if rm is not None else "",
                "insert_len_median": med_ins,
                "support_reads": len(reads), "ref_spanning_reads": refn,
                "allele_fraction": round(af, 3),
                "support_vs_depth": round(len(reads) / med_depth, 2) if med_depth else "",
                "cassette_orientation": ("reverse" if revs and sum(revs) > len(revs) / 2
                                         else "forward" if revs else ""),
                "evidence": "+".join(hows),
                **ann,
            }
            keep = len(reads) >= MIN_READS and af >= MIN_AF
            row["verdict"] = "accepted" if keep else (
                f"rejected: only {len(reads)} read(s)" if len(reads) < MIN_READS
                else f"rejected: allele fraction {af:.2f} < {MIN_AF}")
            allrows.append(row)
            if keep:
                n_kept += 1
                rows.append(row)
                for i, e in enumerate(items):
                    if e["seq"]:
                        fa.write(f">{clone}_{contig}_{pos}_read{i} len={len(e['seq'])}\n"
                                 f"{e['seq']}\n")
        bam.close()
        c2 = next((m for r, (m, _) in depth.items() if r == ABSENT), None)
        copies = (cas_med / med_depth) if med_depth else 0
        summary.append({
            "clone": clone, "platform": PLATFORM,
            "genome_median_depth": round(med_depth, 1),
            "cassette_core_depth": cas_med,
            "est_cassette_copies": round(copies, 2),
            "n_sites_accepted": n_kept,
            "absent_replicon_mean_depth": round(c2, 2) if c2 is not None else "",
            "absent_replicon_present": ("no" if c2 is not None and c2 < 1
                                 else "yes" if c2 is not None else ""),
        })
        print(f"{clone}: {len(ev)} raw evidence -> {n_kept} accepted site(s); "
              f"genome depth ~{med_depth:.0f}x, cassette core {cas_med}x "
              f"-> ~{copies:.2f} copies"
              + (f"; {ABSENT} {c2:.1f}x" if c2 is not None else ""), file=sys.stderr)
    fa.close()

    order = ["clone", "platform", "contig", "insertion_pos", "tsd_start", "tsd_end",
             "tsd_len_junction", "tsd_len_insert", "left_junction", "right_junction",
             "insert_len_median", "support_reads", "ref_spanning_reads",
             "allele_fraction", "support_vs_depth", "cassette_orientation",
             "evidence", "locus_tag", "gene", "product", "cds_start", "cds_end",
             "cds_strand", "cds_len_bp", "bp_into_cds", "pct_into_cds", "codon",
             "aa_pos_of", "effect", "context", "verdict"]
    for name, data in ((snakemake.output.insertions, rows),
                       (snakemake.output.insertions_all, allrows)):
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, order, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(sorted(data, key=lambda r: (r["clone"], r["contig"],
                                                    r["insertion_pos"])))
    with open(snakemake.output.clone_summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(summary)
    with open(snakemake.output.coverage, "w", newline="") as fh:
        w = csv.DictWriter(fh, ["clone", "reference", "length_bp", "mean_depth"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(cov)
    print(f"\n{len(rows)} accepted sites across {len(snakemake.input.bams)} clones",
          file=sys.stderr)


main()
