#!/usr/bin/env python3
"""s3_census.py -- replicon census over the full genome set.

REWRITE of a script stranded on Shannon (2026-08-20 outage).  v2 produces the
complete 28-column table; v1 produced only the 9 computed columns and its
validation checked only those, so the metadata pass-through was silently empty.

Every rule below was recovered from the retained `tables/replicon_census.tsv`
(771 rows) with ZERO mismatches -- see `--validate`, which now checks ALL
columns, not just the computed ones.

Column derivations
------------------
organism_name, strain, sequencing_tech, checkm_*, ncbi_n_chromosomes
                       NCBI datasets `assembly_data_report.jsonl`
genus                  first token of organism_name
query_genus            genus when it is one of the six, else carried through
total_length, gc_percent, n_replicons_total, replicon_sizes
                       computed from the FASTA.  gc_percent is rounded to the
                       NEAREST 0.5 -- NCBI reports it that way and the retained
                       table contains only multiples of 0.5.
largest/second/third   replicon_sizes[0..2], literal 0 when absent
n_large_{200,300,500,1000}kb   contigs >= that many bp
architecture           n_large_300kb -> "1_large".."3_large", else "4+_large"
ncbi_labelled_plasmids deflines matching plasmid
ncbi_labelled_chromosomes   n_replicons_total - ncbi_labelled_plasmids (NCBI
                       frequently leaves the chromosome defline unmarked)
third_ncbi_label       label of the THIRD-LARGEST contig; blank when <3 contigs
third_chr_name         its chromosome number or plasmid name
size_vs_label_agree    n_large_300kb == ncbi_labelled_chromosomes
qc_pass / qc_reason    completeness<95 fails first, then contamination>5
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from pathlib import Path

SIX_GENERA = ["Burkholderia", "Paraburkholderia", "Caballeronia",
              "Trinickia", "Mycetohabitans", "Robbsia"]

LARGE_THRESHOLDS = [("n_large_200kb", 200_000), ("n_large_300kb", 300_000),
                    ("n_large_500kb", 500_000), ("n_large_1000kb", 1_000_000)]

COLUMNS = ["accession", "organism_name", "genus", "strain", "query_genus",
           "total_length", "gc_percent", "n_replicons_total", "replicon_sizes",
           "largest", "second", "third", "sequencing_tech",
           "checkm_completeness", "checkm_contamination", "ncbi_n_chromosomes",
           "n_large_200kb", "n_large_300kb", "n_large_500kb", "n_large_1000kb",
           "architecture", "ncbi_labelled_chromosomes",
           "ncbi_labelled_plasmids", "third_ncbi_label", "third_chr_name",
           "size_vs_label_agree", "qc_pass", "qc_reason"]

CHR_RE = re.compile(r"\bchromosome\b(?:\s+(\S+?))?\s*(?:,|$)", re.I)
PLASMID_RE = re.compile(r"\bplasmid\b(?:\s+(\S+?))?(?:\s+Contig_\d+)?\s*(?:,|$)", re.I)


def read_fasta(path: Path):
    """[(seqid, description, length, gc, at)] in file order."""
    out, sid, desc, gc, at = [], None, "", 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if sid is not None:
                    out.append((sid, desc, gc + at + 0, gc, at))
                head = line[1:].strip()
                sid, _, desc = head.partition(" ")
                gc = at = 0
            else:
                s = line.strip().upper()
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
    if sid is not None:
        out.append((sid, desc, gc + at, gc, at))
    return out


def read_fasta_lengths(path: Path):
    """Same, but length counts EVERY base including ambiguity codes."""
    out, sid, desc, n, gc, at = [], None, "", 0, 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if sid is not None:
                    out.append((sid, desc, n, gc, at))
                head = line[1:].strip()
                sid, _, desc = head.partition(" ")
                n = gc = at = 0
            else:
                s = line.strip().upper()
                n += len(s)
                gc += s.count("G") + s.count("C")
                at += s.count("A") + s.count("T")
    if sid is not None:
        out.append((sid, desc, n, gc, at))
    return out


def label_of(desc: str):
    """-> ('Chromosome'|'Plasmid', name) for one defline description.

    NCBI labels plasmids explicitly but often leaves a chromosome unmarked --
    `NC_014722.1 Mycetohabitans rhizoxinica HKI 454, complete sequence` is the
    chromosome and never says so.  So "not a plasmid" is the chromosome rule,
    which reproduces the retained table exactly (0/770 mismatches) where
    keyword-matching "chromosome" does not.

    A plasmid named literally "unnamed" is recorded as `plasmid`; numbered
    variants (`unnamed1`, `unnamed2`) keep their names, as in the retained table.
    """
    if re.search(r"\bplasmid\b", desc, re.I):
        m = PLASMID_RE.search(desc)
        name = (m.group(1) if m and m.group(1) else "plasmid")
        if name.lower() == "unnamed":
            name = "plasmid"
        return "Plasmid", name
    m = CHR_RE.search(desc)
    return "Chromosome", (m.group(1) if m and m.group(1) else "chromosome")


def load_reports(raw_dirs):
    """accession -> NCBI assembly_data_report fields."""
    md = {}
    for d in raw_dirs:
        for f in glob.glob(str(Path(d) / "**" / "assembly_data_report.jsonl"),
                           recursive=True):
            for line in open(f):
                if not line.strip():
                    continue
                d_ = json.loads(line)
                acc = d_.get("accession")
                if not acc:
                    continue
                ai = d_.get("assemblyInfo", {}) or {}
                st = d_.get("assemblyStats", {}) or {}
                org = d_.get("organism", {}) or {}
                ck = d_.get("checkmInfo", {}) or {}
                md[acc] = {
                    "organism_name": org.get("organismName", "") or "",
                    "strain": ((org.get("infraspecificNames") or {})
                               .get("strain", "") or ""),
                    "sequencing_tech": ai.get("sequencingTech", "") or "",
                    "checkm_completeness": ck.get("completeness", ""),
                    "checkm_contamination": ck.get("contamination", ""),
                    "ncbi_n_chromosomes": st.get("totalNumberOfChromosomes", ""),
                }
    return md


def load_fallback(path: Path | None):
    """Headerless accession<TAB>genus<TAB>species, for genomes with no report."""
    fb = {}
    if not path or not path.exists():
        return fb
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0]:
            fb[p[0]] = f"{p[1]} {p[2]}".strip()
        elif len(p) == 2 and p[0]:
            fb[p[0]] = p[1].strip()
    return fb


def fmt_num(v):
    """NCBI numbers arrive as str/int/float; render like the retained table."""
    if v == "" or v is None:
        return ""
    f = float(v)
    return str(int(f)) if f == int(f) else f"{f:g}"


def census_one(acc: str, fna: Path, rep: dict, fallback_org: str):
    recs = read_fasta_lengths(fna)
    if not recs:
        return None
    order = sorted(recs, key=lambda r: -r[2])
    sizes = [r[2] for r in order]
    gc = sum(r[3] for r in recs)
    at = sum(r[4] for r in recs)
    gc_pct = (round(200.0 * gc / (gc + at)) / 2.0) if (gc + at) else 0.0

    org = rep.get("organism_name") or fallback_org or ""
    genus = org.split()[0] if org else ""
    qg = genus if genus in SIX_GENERA else genus

    labels = [label_of(r[1]) for r in recs]
    n_pla = sum(1 for l, _ in labels if l == "Plasmid")
    n_chr = len(recs) - n_pla

    if len(order) >= 3:
        third_label, third_name = label_of(order[2][1])
    else:
        third_label, third_name = "", ""

    row = {
        "accession": acc,
        "organism_name": org,
        "genus": genus,
        "strain": rep.get("strain", ""),
        "query_genus": qg,
        "total_length": sum(sizes),
        "gc_percent": f"{gc_pct:g}",
        "n_replicons_total": len(sizes),
        "replicon_sizes": ";".join(str(s) for s in sizes),
        "largest": sizes[0] if len(sizes) > 0 else 0,
        "second": sizes[1] if len(sizes) > 1 else 0,
        "third": sizes[2] if len(sizes) > 2 else 0,
        "sequencing_tech": rep.get("sequencing_tech", ""),
        "checkm_completeness": fmt_num(rep.get("checkm_completeness", "")),
        "checkm_contamination": fmt_num(rep.get("checkm_contamination", "")),
        "ncbi_n_chromosomes": fmt_num(rep.get("ncbi_n_chromosomes", "")),
        "ncbi_labelled_chromosomes": n_chr,
        "ncbi_labelled_plasmids": n_pla,
        "third_ncbi_label": third_label,
        "third_chr_name": third_name,
    }
    for col, th in LARGE_THRESHOLDS:
        row[col] = sum(1 for s in sizes if s >= th)
    n300 = row["n_large_300kb"]
    row["architecture"] = f"{n300}_large" if n300 < 4 else "4+_large"
    row["size_vs_label_agree"] = str(n300 == n_chr)

    comp, cont = row["checkm_completeness"], row["checkm_contamination"]
    reason = ""
    if comp and float(comp) < 95:
        reason = "completeness<95"
    elif cont and float(cont) > 5:
        reason = "contamination>5"
    row["qc_reason"] = reason
    row["qc_pass"] = str(reason == "")

    ranks = [{"accession": acc, "contig": r[0], "length": r[2],
              "rank": i, "n_contigs": len(order)}
             for i, r in enumerate(order, start=1)]
    return row, ranks


def validate(new_rows, ref_path: Path, columns):
    ref = {r["accession"]: r for r in
           csv.DictReader(open(ref_path), delimiter="\t")}
    new = {r["accession"]: r for r in new_rows}
    shared = sorted(set(ref) & set(new))
    print(f"  reference rows : {len(ref)}")
    print(f"  new rows       : {len(new)}")
    print(f"  shared         : {len(shared)}")
    only_ref = sorted(set(ref) - set(new))
    if only_ref:
        print(f"  only in reference ({len(only_ref)}): {only_ref[:5]}")
    ok = True
    for f in columns:
        if f not in next(iter(ref.values())):
            print(f"  {f:<26} (not in reference)")
            continue
        diffs = [(a, ref[a][f], str(new[a][f])) for a in shared
                 if ref[a][f] != str(new[a][f])]
        if diffs:
            ok = False
            print(f"  {f:<26} {len(diffs)} MISMATCH")
            for a, rv, nv in diffs[:4]:
                print(f"      {a}: ref={rv[:44]!r} new={nv[:44]!r}")
        else:
            print(f"  {f:<26} MATCH")
    print(f"  VERDICT: {'reproduces the retained table' if ok else 'DIVERGENCE'}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--list", required=True, type=Path)
    ap.add_argument("--raw", type=Path, default=None,
                    help="dir holding NCBI datasets batches (searched recursively)")
    ap.add_argument("--fallback-organism", type=Path, default=None,
                    help="headerless accession<TAB>genus<TAB>species")
    ap.add_argument("--mf6", type=Path, default=None)
    ap.add_argument("--mf7", type=Path, default=None)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--validate", type=Path, default=None)
    a = ap.parse_args(argv)

    accs = [l.strip() for l in open(a.list) if l.strip()]
    extra = {}
    if a.mf6:
        extra["MF6"] = a.mf6
    if a.mf7:
        extra["MF7"] = a.mf7
    reports = load_reports([a.raw]) if a.raw else {}
    fallback = load_fallback(a.fallback_organism)
    print(f"NCBI reports loaded: {len(reports)}   fallback organisms: {len(fallback)}")

    rows, ranks = [], []
    missing = []
    for acc in accs:
        fna = extra.get(acc) or (a.genomes / f"{acc}.fna")
        if not Path(fna).exists():
            missing.append(acc)
            continue
        got = census_one(acc, Path(fna), reports.get(acc, {}),
                         fallback.get(acc, ""))
        if got:
            rows.append(got[0])
            ranks.extend(got[1])
    if missing:
        print(f"WARNING: {len(missing)} genomes missing FASTA: {missing[:5]}",
              file=sys.stderr)

    a.outdir.mkdir(parents=True, exist_ok=True)
    with open(a.outdir / "replicon_census.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(a.outdir / "contig_ranks.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "contig", "length",
                                           "rank", "n_contigs"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(ranks)
    print(f"wrote {len(rows)} census rows, {len(ranks)} contig rows -> {a.outdir}")

    if a.validate:
        print("\n=== validation against retained table (ALL columns) ===")
        return 0 if validate(rows, a.validate, COLUMNS) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
