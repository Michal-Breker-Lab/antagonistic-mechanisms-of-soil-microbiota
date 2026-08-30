#!/usr/bin/env python3
"""Fetch each KEGG map's class from the BRITE hierarchy, so non-bacterial and
overview maps can be dropped before testing.
"""

import argparse
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path

BRITE = "https://rest.kegg.jp/get/br:br08901"
NON_BACTERIAL = ["Organismal Systems", "Human Diseases", "Drug Development"]
OVERVIEW_SUBCLASS = "Global and overview maps"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--url", default=BRITE)
    p.add_argument("--timeout", type=int, default=120)
    return p.parse_args()


def fetch(url, timeout):
    r = subprocess.run(["curl", "-sS", "--max-time", str(timeout), url],
                       capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        sys.exit(f"failed to fetch {url}: {r.stderr.strip()[:300]}")
    return r.stdout


def parse(text):
    """-> [(map_id, A class, B subclass, name)] in file order."""
    rows, a, b = [], None, None
    for line in text.splitlines():
        if line.startswith("A"):
            a = line[1:].strip()
        elif line.startswith("B"):
            s = line[1:].strip()
            if s:
                b = s
        elif line.startswith("C"):
            m = re.match(r"C\s+(\d{5})\s+(.*)", line)
            if m:
                rows.append(("ko" + m.group(1), a, b, m.group(2).strip()))
    return rows


def main():
    args = parse_args()
    rows = parse(fetch(args.url, args.timeout))
    if len(rows) < 300:
        sys.exit(f"only {len(rows)} maps parsed from br08901 - refusing to write a "
                 f"truncated class table")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_nb = sum(1 for r in rows if r[1] in NON_BACTERIAL)
    n_ov = sum(1 for r in rows if r[2] == OVERVIEW_SUBCLASS)
    with args.out.open("w") as fh:
        fh.write(f"# KEGG pathway classes from BRITE br08901, fetched "
                 f"{dt.date.today().isoformat()} from {args.url}\n"
                 f"# class / subclass are KEGG's own; nothing here is curated by hand\n"
                 f"# is_non_bacterial = class in {NON_BACTERIAL}\n"
                 f"#   these maps cannot describe a bacterium; eggNOG reaches them\n"
                 f"#   through KOs shared with eukaryote-specific pathways\n"
                 f"# is_overview     = subclass is '{OVERVIEW_SUBCLASS}'\n"
                 f"#   supersets, not pathways - ko01100 alone covers ~46% of the\n"
                 f"#   KEGG-annotated proteome\n"
                 f"# {len(rows)} maps, {n_nb} non-bacterial, {n_ov} overview\n")
        fh.write("pathway_id\tclass\tsubclass\tname\tis_non_bacterial\tis_overview\n")
        for pid, a, b, name in rows:
            fh.write(f"{pid}\t{a}\t{b}\t{name}\t"
                     f"{'yes' if a in NON_BACTERIAL else 'no'}\t"
                     f"{'yes' if b == OVERVIEW_SUBCLASS else 'no'}\n")
    print(f"{len(rows)} maps -> {args.out}")
    print(f"  non-bacterial (Organismal Systems / Human Diseases / Drug Development): {n_nb}")
    print(f"  global and overview maps: {n_ov}")


if __name__ == "__main__":
    main()
