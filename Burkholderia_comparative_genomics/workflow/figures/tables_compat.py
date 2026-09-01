"""Small compatibility shims for tables whose schema moved in the rebuild."""
import csv
from pathlib import Path


def organism_names(tab: Path) -> dict:
    """accession -> organism_name.

    The original c3_calls_all_genomes.tsv carried organism_name; the rebuilt one
    carries only accession/c3_present/evidence/n_secondary_large, and the names
    live in host_categories.tsv. Try each in turn so a script runs against either
    layout instead of dying on a KeyError.
    """
    out = {}
    for src in ("host_categories.tsv", "c3_calls_all_genomes.tsv"):
        f = Path(tab) / src
        if not f.exists():
            continue
        with open(f, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if "organism_name" not in (rd.fieldnames or []):
                continue
            for r in rd:
                out.setdefault(r["accession"], r["organism_name"])
        if out:
            return out
    raise SystemExit(f"FAIL: no organism_name column found under {tab}")
