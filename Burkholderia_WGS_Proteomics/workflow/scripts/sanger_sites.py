"""Call the insertion site for the Sanger-only clones from their aligned flanks."""
import csv
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tn_annotate import annotate, parse_gff  # noqa: E402

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

MIN_MAPQ = int(snakemake.params.min_mapq)  # noqa: F821
CASSETTE = snakemake.params.cassette  # noqa: F821

feats = parse_gff(snakemake.input.gff)  # noqa: F821

wgs = {}
with open(snakemake.input.wgs) as fh:  # noqa: F821
    for r in csv.DictReader(fh, delimiter="\t"):
        wgs.setdefault((r["contig"], int(r["insertion_pos"])), []).append(r)

placed = {}
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # noqa: F821
for aln in bam:
    if aln.is_unmapped or aln.is_supplementary or aln.is_secondary:
        continue
    clone, direction, _how = aln.query_name.split("|")
    contig = bam.get_reference_name(aln.reference_id)
    if contig == CASSETTE:
        print(f"{aln.query_name}: flank aligned to {CASSETTE}, not genomic - skipped",
              file=sys.stderr)
        continue
    if aln.mapping_quality < MIN_MAPQ:
        print(f"{aln.query_name}: MAPQ {aln.mapping_quality} < {MIN_MAPQ} - skipped",
              file=sys.stderr)
        continue
    pos = (aln.reference_end if aln.is_reverse else aln.reference_start + 1)
    placed.setdefault((clone, contig, pos), []).append(
        {"direction": direction, "mapq": aln.mapping_quality,
         "strand": "-" if aln.is_reverse else "+",
         "aligned_bp": aln.query_alignment_length})
    print(f"{aln.query_name}: {contig}:{pos} ({'-' if aln.is_reverse else '+'}, "
          f"MAPQ {aln.mapping_quality}, {aln.query_alignment_length} bp aligned)",
          file=sys.stderr)
bam.close()

best = {}
for (clone, contig, pos), hits in placed.items():
    score = (len(hits), sum(h["aligned_bp"] for h in hits))
    if clone not in best or score > best[clone][0]:
        best[clone] = (score, contig, pos, hits)

rows = []
for clone in sorted(best):
    (n_reads, _bp), contig, pos, hits = best[clone]
    same = wgs.get((contig, pos), [])
    tsd = same[0]["tsd_len_junction"] if same else ""
    row = {
        "clone": clone, "platform": "Sanger",
        "contig": contig, "insertion_pos": pos,
        "tsd_start": "", "tsd_end": "",
        "tsd_len_junction": tsd, "tsd_len_insert": "",
        "left_junction": "", "right_junction": "", "insert_len_median": "",
        "support_reads": n_reads, "ref_spanning_reads": "",
        "allele_fraction": "", "support_vs_depth": "",
        "cassette_orientation": "",
        "evidence": "sanger_flank:" + ",".join(sorted(h["direction"] for h in hits)),
        "tsd_source": (f"WGS call in {','.join(sorted(r['clone'] for r in same))}"
                       if same else ""),
        **annotate(feats, contig, pos),
        "verdict": "accepted",
    }
    rows.append(row)
    print(f"\n{clone}: {contig}:{pos} from {n_reads} flank(s)"
          + (f"; TSD {tsd} carried from {row['tsd_source']}" if tsd else ""),
          file=sys.stderr)

if not rows:
    sys.exit("no Sanger clone could be placed")

cols = list(rows[0])
with open(snakemake.output.sites, "w", newline="") as fh:  # noqa: F821
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"\n{len(rows)} Sanger site(s)", file=sys.stderr)
