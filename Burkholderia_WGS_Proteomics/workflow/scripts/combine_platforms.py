"""Concatenate the per-platform insertion calls into one site table and one clone
table.
"""
import csv
import sys

sys.stderr = open(snakemake.log[0], "w")


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


sites, clones = [], []
for platform, ins, summ in zip(snakemake.params.platforms,
                               snakemake.input.insertions,
                               snakemake.input.summaries):
    for r in rd(ins):
        r.setdefault("platform", platform)
        sites.append(r)
    for r in rd(summ):
        r.setdefault("platform", platform)
        clones.append(r)

seen = {}
for r in sites:
    seen.setdefault(r["clone"], set()).add(r["platform"])
for r in sites + clones:
    r["both_platforms"] = "yes" if len(seen.get(r["clone"], ())) > 1 else "no"

at = {}
for r in sites:
    at.setdefault((r["contig"], r["insertion_pos"]), []).append(r)
for (contig, pos), group in at.items():
    names = sorted(g["clone"] for g in group)
    plats = sorted({g["platform"] for g in group})
    for g in group:
        g["site_id"] = f"{contig}:{pos}"
        g["clones_at_site"] = ",".join(names)
        g["n_clones_at_site"] = len(group)
        g["cross_platform_confirmed"] = "yes" if len(plats) > 1 else "no"

for r in clones:
    n = int(r["n_sites_accepted"])
    cp = float(r["est_cassette_copies"])
    r["copy_call"] = (
        "no insertion detected" if n == 0 else
        f"{n} site(s), depth consistent" if abs(cp - n) <= 0.5 else
        f"{n} site(s) but ~{cp:.1f}x cassette depth - possible tandem/multimer"
        if cp > n else
        f"{n} site(s) but only ~{cp:.1f}x cassette depth - check subclonality")

for path, data, cols in ((snakemake.output.sites, sites, list(sites[0]) if sites else []),
                         (snakemake.output.clones, clones, list(clones[0]) if clones else [])):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(data, key=lambda r: (r["platform"], r["clone"])))

print(f"{len(sites)} sites across {len(clones)} clones", file=sys.stderr)
