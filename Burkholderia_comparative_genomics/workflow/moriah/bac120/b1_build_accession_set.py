#!/usr/bin/env python3
"""b1 -- build the bac120 tree's genome set.

Burkholderiaceae complete+chromosome -> sensu lato genera -> GCA/GCF paired ->
isolation-source/host metadata present.  Writes final_accessions.json.

THREE things here are load-bearing and each was a bug waiting to happen:

1. QUERY THE FAMILY, FILTER BY GENUS NAME. Querying genus taxids directly
   silently returns ZERO for Trinickia -- its genus node is not resolvable
   through the datasets taxon endpoint (both 2571227 and 3716038 give 0, while
   the species taxids 28094 and 2571746 return records). Trinickia, Robbsia and
   Pararobbsia are the DEEPEST outgroup members, so a genus-taxid query would
   have quietly dropped the root of the tree while looking like it worked. Same
   class of trap as GTDB Patescibacteriota vs NCBI Patescibacteria.

2. ASSERT THE PAGE COUNTS. This API truncates paged streams at round numbers and
   still exits clean. Every level asserts got == total_count.

3. PAIR GCA/GCF. 1,816 sensu lato records collapse to 919 genomes; without this
   the tree would be almost exactly 2x too big -- the error this project already
   made once (see the rebuild plan's Task 1).

Metadata filter: isolation_source OR env_medium OR host (biosample attribute or
top-level field), rejecting the nullish vocabulary NCBI actually contains
("missing", "not applicable", "-", ...). Including host matches how
host_categories.tsv is built -- plant_compartment.py reads raw_isolation_source
AND raw_host -- so this is the project's existing convention, not a new one.

Observed 2026-08-25: 3,159 family records -> 1,816 sensu lato -> 919 paired ->
790 with metadata. Plus MF6 and MF7 = 792 genomes.
"""
from __future__ import annotations

import collections
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/119060/dataset_report"
GENERA = {"Burkholderia", "Paraburkholderia", "Caballeronia", "Trinickia",
          "Mycetohabitans", "Robbsia", "Pararobbsia"}
LEVELS = ("complete_genome", "chromosome")
NULLISH = re.compile(
    r"^\s*(|-|--|n/?a|na|none|null|unknown|unspecified|missing|not\s*"
    r"(applicable|collected|provided|determined|available)|nan|\?)\s*$", re.I)
OUT = Path(__file__).resolve().parent.parent


def fetch(level, token=None):
    u = (f"{BASE}?filters.assembly_version=current&filters.assembly_level={level}"
         f"&page_size=1000&returned_content=COMPLETE")
    if token:
        u += f"&page_token={token}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(u, timeout=180) as r:
                return json.load(r)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))


def genus(r):
    return ((r.get("organism", {}).get("organism_name", "?") or "?").split() or ["?"])[0]


def attrs(r):
    return {x.get("name"): (x.get("value") or "") for x in
            r.get("assembly_info", {}).get("biosample", {}).get("attributes", []) or []}


def good(v):
    return bool(v) and not NULLISH.match(v)


def has_meta(r):
    a = attrs(r)
    bs = r.get("assembly_info", {}).get("biosample", {})
    return any(good(v) for v in (a.get("isolation_source", ""), a.get("env_medium", ""),
                                 a.get("host", ""), bs.get("host") or ""))


def main():
    rows = []
    for lvl in LEVELS:
        tok, got, tot = None, 0, None
        while True:
            d = fetch(lvl, tok)
            if tot is None:
                tot = d.get("total_count", 0)
            batch = d.get("reports", [])
            rows += batch
            got += len(batch)
            tok = d.get("next_page_token")
            if not tok:
                break
        assert got == tot, f"{lvl}: got {got}, expected {tot} - TRUNCATED STREAM"
        print(f"  Burkholderiaceae {lvl:<16} {got}", file=sys.stderr)

    sl = [r for r in rows if genus(r) in GENERA]
    print(f"  -> sensu lato genera: {len(sl)} of {len(rows)} family records", file=sys.stderr)

    by = {r["accession"]: r for r in sl}
    keep = {}
    for r in sl:
        acc = r["accession"]
        pair = r.get("paired_accession") or ""
        # keep RefSeq always; keep a GenBank record only when it has no RefSeq twin here
        if acc.startswith("GCF_") or not (pair and pair in by):
            keep[acc] = r
    print(f"  -> after GCA/GCF pairing: {len(keep)}", file=sys.stderr)

    final = {k: v for k, v in keep.items() if has_meta(v)}
    print(f"  -> with isolation_source/host metadata: {len(final)}", file=sys.stderr)
    for k, v in collections.Counter(genus(r) for r in final.values()).most_common():
        print(f"       {k:<18} {v}", file=sys.stderr)

    assert final, "empty genome set - refusing to write"
    (OUT / "final_accessions.json").write_text(json.dumps(sorted(final), indent=0))
    meta = {a: {"organism": r.get("organism", {}).get("organism_name"),
                "level": r.get("assembly_info", {}).get("assembly_level"),
                "isolation_source": attrs(r).get("isolation_source", ""),
                "env_medium": attrs(r).get("env_medium", ""),
                "host": attrs(r).get("host", "")
                        or (r.get("assembly_info", {}).get("biosample", {}).get("host") or "")}
            for a, r in final.items()}
    (OUT / "genome_metadata.json").write_text(json.dumps(meta, indent=1, sort_keys=True))
    print(f"\nFINAL: {len(final)} NCBI genomes + MF6 + MF7 = {len(final) + 2}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
