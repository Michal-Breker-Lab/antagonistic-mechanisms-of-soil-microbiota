#!/usr/bin/env python3
"""blastp a query set against the all-genome pyrodigal protein DB, and tier the hits.

Reproduces the original search recorded in TOOLS.md:

    blastp -query queries.faa -db all_proteins -outfmt "6 qseqid sseqid pident
        length nident mismatch gapopen qstart qend sstart send evalue bitscore
        qcovhsp qcovs slen qlen" -evalue 1e-5 -max_target_seqs 1000000
        -comp_based_stats 2

Hit positions come from the pyrodigal coords table written in the same pass as
the proteins (s23_pyrodigal.py), so the gene index in the subject id and the
coordinates cannot drift apart.

Tier rule, recovered from the retained tables and validated at 0 mismatches on
both toxin12_search_all_hits.tsv (10,754 rows) and rhs_search_all_hits.tsv
(3,792 rows):

    tier 1  pident >= 60 and qcovs >= 70     strong, near-full-length
    tier 2  pident >= 60 and qcovs >= 30     strong but partial
    tier 3  pident >= 40 and qcovs >= 70     remote homolog, good coverage
    tier 4  otherwise

CAVEAT: tier 3's coverage cut is NOT identifiable from the retained data -- 70,
80 and 85 all reproduce it exactly, because no retained hit falls between 70 and
85 at pident 40-60. 70 is used here as the parsimonious reading (symmetric with
tier 1). Any NEW hit landing in that window is reported by --flag-ambiguous so
it can be inspected rather than silently tiered.
"""
from __future__ import annotations

import argparse
import collections
import csv
import subprocess
import sys
from pathlib import Path

OUTFMT = ("6 qseqid sseqid pident length nident mismatch gapopen qstart qend "
          "sstart send evalue bitscore qcovhsp qcovs slen qlen")
FIELDS = OUTFMT.split()[1:]


def tier_of(pident: float, qcovs: float) -> str:
    if pident >= 60 and qcovs >= 70:
        return "1"
    if pident >= 60 and qcovs >= 30:
        return "2"
    if pident >= 40 and qcovs >= 70:
        return "3"
    return "4"


def ambiguous(pident: float, qcovs: float) -> bool:
    """Hit whose tier depends on the unidentifiable tier-3 coverage cut."""
    return 40 <= pident < 60 and 70 <= qcovs < 85


def load_coords(path: Path, wanted: set | None = None) -> dict:
    """(acc, contig, gene) -> (start, end, strand, contig_len).

    `wanted` restricts the table to the keys the hits actually reference. The
    full table is ~5M rows / 3 GB resident; a search returns ~10k hits, so
    filtering keeps this comfortably small and lets the join run anywhere.
    """
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["acc"], r["contig"], r["gene"])
            if wanted is not None and k not in wanted:
                continue
            out[k] = (r["start"], r["end"], r["strand"], r["contig_len"])
    return out


def load_ranks(path: Path) -> dict:
    """(acc, contig) -> size rank, from the validated contig_ranks table."""
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        # the retained table calls this column "rank"; the typing table calls the
        # same quantity "size_rank". Accept either rather than depending on which
        # file the caller passed.
        col = next((c for c in ("size_rank", "rank") if c in (rd.fieldnames or [])), None)
        if col is None:
            raise SystemExit(f"no rank column in {path}: {rd.fieldnames}")
        for r in rd:
            out[(r["accession"], r["contig"])] = r[col]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queries", required=True, type=Path)
    ap.add_argument("--db", required=True)
    ap.add_argument("--coords", required=True, type=Path)
    ap.add_argument("--ranks", required=True, type=Path)
    ap.add_argument("--prefix", required=True, help="output basename, e.g. toxin12")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--blastp", default="blastp")
    ap.add_argument("--raw", type=Path, default=None,
                    help="reuse an existing raw blast table instead of running")
    a = ap.parse_args(argv)

    a.outdir.mkdir(parents=True, exist_ok=True)
    raw = a.raw or (a.outdir / f"{a.prefix}_blastp_raw.tsv")

    if not a.raw:
        cmd = [a.blastp, "-query", str(a.queries), "-db", a.db,
               "-outfmt", OUTFMT, "-evalue", a.evalue,
               "-max_target_seqs", "1000000",
               "-num_threads", str(a.threads),
               "-comp_based_stats", "2", "-out", str(raw)]
        print("running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

    wanted = set()
    with open(raw) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            sid = ln.split("\t")[1]
            parts = sid.split("|")
            if len(parts) >= 3:
                wanted.add((parts[0], parts[1], parts[2]))
    print(f"distinct subject genes referenced by hits: {len(wanted):,d}", flush=True)

    coords = load_coords(a.coords, wanted)
    ranks = load_ranks(a.ranks)
    print(f"coords {len(coords):,d}  ranks {len(ranks):,d}", flush=True)

    out_path = a.outdir / f"{a.prefix}_search_all_hits.tsv"
    per_query = collections.defaultdict(lambda: collections.defaultdict(set))
    n = miss = amb = 0
    with open(raw) as fh, open(out_path, "w", newline="") as oh:
        w = csv.writer(oh, delimiter="\t")
        w.writerow(["query", "acc", "contig", "gene", "start", "end", "strand",
                    "pident", "qcovs", "qcovhsp", "qstart", "qend", "bits",
                    "evalue", "slen", "qlen", "tier", "contig_len",
                    "contig_rank", "tier_ambiguous"])
        for ln in fh:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            # FIELDS already begins with qseqid, so zip against the WHOLE row.
            # Zipping against p[1:] silently shifts every name by one and puts
            # pident into rec["sseqid"] -- which then fails to split on "|" and
            # discards every hit.
            rec = dict(zip(FIELDS, p))
            q = rec["qseqid"]
            sid = rec["sseqid"]
            try:
                acc, contig, gene = sid.split("|")
            except ValueError:
                miss += 1
                continue
            key = (acc, contig, gene)
            if key not in coords:
                miss += 1
                continue
            start, end, strand, clen = coords[key]
            pid, qcs = float(rec["pident"]), float(rec["qcovs"])
            t = tier_of(pid, qcs)
            am = ambiguous(pid, qcs)
            amb += am
            w.writerow([q, acc, contig, gene, start, end, strand,
                        rec["pident"], rec["qcovs"], rec["qcovhsp"],
                        rec["qstart"], rec["qend"], rec["bitscore"],
                        rec["evalue"], rec["slen"], rec["qlen"], t, clen,
                        ranks.get((acc, contig), ""), "True" if am else "False"])
            per_query[q][t].add(acc)
            n += 1

    print(f"hits written {n:,d}  (subject ids unresolved: {miss})")
    if amb:
        print(f"WARNING {amb} hits fall in the unidentifiable tier-3 window "
              f"(pident 40-60, qcovs 70-85); flagged tier_ambiguous=True")

    genomes = {a_ for a_ in {k[0] for k in coords}}
    summ = a.outdir / f"{a.prefix}_search_summary.tsv"
    with open(summ, "w", newline="") as oh:
        w = csv.writer(oh, delimiter="\t")
        w.writerow(["query", "n_genomes_any_hit", "n_tier1", "n_tier2_only",
                    "n_tier3_only", "n_genomes_total"])
        for q in sorted(per_query):
            d = per_query[q]
            t1, t2, t3 = d.get("1", set()), d.get("2", set()), d.get("3", set())
            any_hit = set().union(*d.values()) if d else set()
            w.writerow([q, len(any_hit), len(t1), len(t2 - t1),
                        len(t3 - t1 - t2), len(genomes)])
    print(f"wrote {out_path}\nwrote {summ}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
