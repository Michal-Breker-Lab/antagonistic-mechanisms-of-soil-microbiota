"""Collapse each genome's MacSyFinder best solution into one row per system."""
import csv
import sys

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

COLS = ["genome", "system", "system_id", "replicon", "model_fqn", "wholeness",
        "score", "n_loci", "n_genes", "components"]


def read_hits(path):
    """macsyfinder's best_solution.tsv, comment and blank lines dropped."""
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def collapse(hits, genome):
    """One row per system, components in the order the hits were reported."""
    systems, order = {}, []
    for h in hits:
        sid = h["sys_id"]
        if sid not in systems:
            order.append(sid)
            systems[sid] = {
                "genome": genome,
                "system": h["model_fqn"].rsplit("/", 1)[-1],
                "system_id": sid,
                "replicon": h["replicon"],
                "model_fqn": h["model_fqn"],
                "wholeness": h["sys_wholeness"],
                "score": h["sys_score"],
                "n_loci": h["sys_loci"],
                "_components": [],
            }
        systems[sid]["_components"].append(h["gene_name"])
    out = []
    for sid in order:
        s = systems[sid]
        comps = s.pop("_components")
        s["n_genes"] = len(comps)
        s["components"] = ",".join(comps)
        out.append(s)
    return out


rows = []
for genome, path in zip(snakemake.params.genomes,      # noqa: F821
                        snakemake.input.relatives):    # noqa: F821
    hits = read_hits(path)
    got = collapse(hits, genome)
    rows.extend(got)
    print(f"{genome}: {len(hits)} hits -> {len(got)} system(s)", file=sys.stderr)

mf6 = collapse(read_hits(snakemake.input.mf6), "MF6")  # noqa: F821
rows.extend(mf6)
print(f"MF6: {len(mf6)} system(s) across its replicons", file=sys.stderr)

with open(snakemake.output.systems, "w", newline="") as fh:  # noqa: F821
    w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"\n{len(rows)} system(s) across {len(snakemake.params.genomes) + 1} genomes",  # noqa: F821
      file=sys.stderr)
