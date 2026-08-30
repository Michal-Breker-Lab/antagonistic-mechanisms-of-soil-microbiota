#!/usr/bin/env python3
"""Extract the per-run MaxLFQ intensity matrix for MF6 from the FragPipe table."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fp_common import (  # noqa: E402
                       classify, read_combined_protein, to_float, write_tsv)

MEASURE = "MaxLFQ Intensity"
KLASS = "MF6"
ANNOT_COLS = ["Protein", "Gene", "Protein Length", "Combined Total Peptides",
              "Combined Unique Spectral Count", "Protein Probability",
              "Description"]


def main() -> int:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    src = Path(sm.input.combined)
    out = Path(sm.output.maxlfq)
    locus_prefix = sm.params.locus_prefix

    header, rows, measure_columns = read_combined_protein(src)
    idx = {c: i for i, c in enumerate(header)}
    for c in ANNOT_COLS:
        if c not in idx:
            sys.exit(f"expected column missing from {src}: {c!r}")
    org_i = idx.get("Organism", -1)

    kept = [r for r in rows
            if classify(r[idx["Protein"]], r[org_i] if org_i >= 0 else "",
                        locus_prefix) == KLASS]
    print(f"{src}")
    print(f"  {len(rows)} proteins -> {len(kept)} after keeping class {KLASS!r}")

    cols = measure_columns.get(MEASURE, [])
    if not cols:
        sys.exit(f"no per-run columns found for {MEASURE!r} in {src}")

    out_rows = [[r[idx[c]] for c in ANNOT_COLS] + [r[i] for _, i in cols]
                for r in kept]
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(out, ANNOT_COLS + [run for run, _ in cols], out_rows)

    n_runs = len(cols)
    detected = [sum(1 for v in row[len(ANNOT_COLS):] if to_float(v) > 0)
                for row in out_rows]
    print(f"  {out.name:34s} {len(out_rows):5d} proteins x {n_runs} runs"
          f"   in all runs: {sum(1 for d in detected if d == n_runs)},"
          f" in none: {sum(1 for d in detected if d == 0)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
