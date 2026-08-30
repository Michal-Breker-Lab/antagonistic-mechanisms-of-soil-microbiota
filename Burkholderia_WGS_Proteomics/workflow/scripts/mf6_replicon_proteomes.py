"""Split MF6's proteome by replicon for the per-replicon TXSScan searches."""
import re
import sys
from collections import OrderedDict

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821


def read_fasta(path):
    seqs, name = OrderedDict(), None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.rstrip("\n"))
    return seqs


prot = read_fasta(snakemake.input.faa)  # noqa: F821
print(f"proteome: {len(prot)} sequences", file=sys.stderr)

cds, seen = [], set()
dups = []
for line in open(snakemake.input.gff):  # noqa: F821
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS":
        continue
    m = re.search(r"locus_tag=([^;]+)", f[8])
    if not m:
        continue
    tag = m.group(1)
    if tag in seen:
        dups.append(tag)
        continue
    seen.add(tag)
    if tag in prot:
        cds.append((f[0], int(f[3]), tag))
if dups:
    print(f"locus tags appearing more than once, kept at first occurrence: "
          f"{', '.join(sorted(set(dups)))}", file=sys.stderr)

by_replicon = {}
for contig, start, tag in sorted(cds, key=lambda c: (c[0], c[1])):
    by_replicon.setdefault(contig, []).append(tag)

for out in snakemake.output:  # noqa: F821
    replicon = out.rsplit("/", 1)[-1][: -len(".faa")]
    tags = by_replicon.get(replicon)
    if not tags:
        sys.exit(f"no CDS on {replicon!r}; annotation has: "
                 f"{', '.join(sorted(by_replicon))}")
    with open(out, "w") as fh:
        for t in tags:
            fh.write(f">{t}\n" + "\n".join(prot[t]) + "\n")
    print(f"{replicon}: {len(tags)} proteins", file=sys.stderr)
