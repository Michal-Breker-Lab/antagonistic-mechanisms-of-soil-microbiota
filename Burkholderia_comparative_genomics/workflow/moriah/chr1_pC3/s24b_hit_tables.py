#!/usr/bin/env python3
"""Build the per-genome / carrier tables the toxin and RHS figures consume.

Takes the tiered hit table from s24_search.py and joins it to the replicon
census (organism, architecture) and the pC3 calls, then emits one row per
(query, genome) -- including genomes with NO hit, which the figures rely on to
compute prevalence denominators correctly.

Column layouts reproduce the retained tables exactly:
  <prefix>_search_per_genome.tsv   20 columns
  <prefix>_carriers.tsv             8 columns
  <prefix>_carrier_replicon.tsv     4 columns   (tier-1 carriers only)

Note on MF7 (D10): MF7 stays IN these tables. Contig fragmentation does not
distort a gene-content measurement, and MF7's toxin profile is a real result
(identical to MF6, its 100%-ANI clone). Only size and architecture panels
exclude it -- see rebuild_rules.py. Its `architecture` value is however the
meaningless size-rank output, so it is overwritten with the alignment-derived
call where a --mf7-architecture override is supplied.
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

PER_GENOME_COLS = ["query", "accession", "organism", "architecture",
                   "pC3_present", "n_hits_any", "n_hits_tier1", "best_tier",
                   "best_pident", "best_qcovs", "best_qstart", "best_qend",
                   "best_bitscore", "best_evalue", "best_contig",
                   "best_contig_len", "best_contig_rank", "best_subject_len",
                   "best_start", "best_end"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, type=Path)
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--calls", type=Path, default=None,
                    help="c3_calls_all_genomes.tsv; without it pC3_present is blank")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--mf7-architecture", default=None,
                    help="override MF7's size-rank architecture, e.g. '3_large'")
    a = ap.parse_args(argv)
    a.outdir.mkdir(parents=True, exist_ok=True)

    census = {r["accession"]: r for r in
              csv.DictReader(open(a.census), delimiter="\t")}
    if a.mf7_architecture and "MF7" in census:
        census["MF7"]["architecture"] = a.mf7_architecture

    calls = {}
    if a.calls and a.calls.exists():
        for r in csv.DictReader(open(a.calls), delimiter="\t"):
            calls[r["accession"]] = r.get("c3_present", "")

    hits = collections.defaultdict(list)
    queries = set()
    for r in csv.DictReader(open(a.hits), delimiter="\t"):
        hits[(r["query"], r["acc"])].append(r)
        queries.add(r["query"])
    print(f"queries {len(queries)}  genomes in census {len(census)}")

    def best_of(rows):
        return max(rows, key=lambda r: float(r["bits"]))

    per_rows, carrier_rows, carrep_rows = [], [], []
    for q in sorted(queries):
        for acc in sorted(census):
            rs = hits.get((q, acc), [])
            c = census[acc]
            pc3 = calls.get(acc, "")
            n_any = len(rs)
            n_t1 = sum(1 for r in rs if r["tier"] == "1")
            if rs:
                b = best_of(rs)
                per_rows.append([q, acc, c.get("organism_name", ""),
                                 c.get("architecture", ""), pc3, n_any, n_t1,
                                 b["tier"], b["pident"], b["qcovs"],
                                 b.get("qstart", ""), b.get("qend", ""),
                                 b["bits"], b["evalue"], b["contig"],
                                 b["contig_len"], b["contig_rank"], b["slen"],
                                 b["start"], b["end"]])
                carrier_rows.append([q, acc, c.get("organism_name", ""), pc3,
                                     str(n_t1 > 0), b["pident"], b["qcovs"], n_any])
                if n_t1 > 0:
                    bt1 = best_of([r for r in rs if r["tier"] == "1"])
                    carrep_rows.append([q, acc, bt1["contig_rank"],
                                        bt1["contig_len"]])
            else:
                per_rows.append([q, acc, c.get("organism_name", ""),
                                 c.get("architecture", ""), pc3, 0, 0] + [""] * 13)
                carrier_rows.append([q, acc, c.get("organism_name", ""), pc3,
                                     "False", "0.0", "0.0", 0])

    def dump(name, header, rows):
        p = a.outdir / name
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header); w.writerows(rows)
        print(f"wrote {p}  ({len(rows)} rows)")

    dump(f"{a.prefix}_search_per_genome.tsv", PER_GENOME_COLS, per_rows)
    dump(f"{a.prefix}_carriers.tsv",
         ["locus", "accession", "organism", "pC3_present", "tier1_carrier",
          "best_pident", "best_qcov", "n_hits"], carrier_rows)
    dump(f"{a.prefix}_carrier_replicon.tsv",
         ["locus", "accession", "best_contig_rank", "best_contig_len"],
         carrep_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
