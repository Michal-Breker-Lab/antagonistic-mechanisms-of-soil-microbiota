#!/usr/bin/env python3
"""Shared parsing of the FragPipe combined_protein table."""
from __future__ import annotations

import re
from pathlib import Path


MEASURES = [
    "MaxLFQ Intensity",
    "Unique Spectral Count",
    "Total Spectral Count",
    "Spectral Count",
    "Intensity",
]

RUN_RE = re.compile(r"^(?P<cond>(?P<strain>[^_]+)_(?P<cult>alone|withC)_(?P<day>d\d+))_(?P<rep>\d+)$")


def classify(protein: str, organism: str = "", locus_prefix: str = "MF6_") -> str:
    """Bucket a FragPipe protein entry into bacterial / host / contaminant."""
    if protein.startswith("contam_"):
        return "contaminant"
    if protein.startswith(locus_prefix):
        return "MF6"
    if "CHLRE" in protein or "Chlamydomonas" in organism:
        return "Chlamydomonas"
    return "other"


def read_combined_protein(path: Path):
    """Return (header, rows, measure_columns).

    measure_columns maps measure -> [(run_name, column_index), ...] in file
    order, so downstream code never has to re-parse column names.
    """
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    measure_columns: dict[str, list[tuple[str, int]]] = {m: [] for m in MEASURES}
    for i, col in enumerate(header):
        for m in MEASURES:
            suffix = " " + m
            if col.endswith(suffix):
                run = col[: -len(suffix)]
                if RUN_RE.match(run):
                    measure_columns[m].append((run, i))
                break

    n = len(header)
    for r in rows:
        if len(r) < n:
            r.extend([""] * (n - len(r)))
    return header, rows, measure_columns


def to_float(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def median(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")
