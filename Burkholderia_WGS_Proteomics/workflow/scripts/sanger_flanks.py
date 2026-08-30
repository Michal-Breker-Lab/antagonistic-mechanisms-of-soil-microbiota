"""Cut the genomic flank out of each Sanger read at the Tn5 mosaic end."""
import csv
import re
import sys

sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

ME = "AGATGTGTATAAGAGACAG"
ME_RC = "CTGTCTCTTATACACATCT"
MIN_FLANK = int(snakemake.params.min_flank)  # noqa: F821


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq).upper()


def flank_of(seq):
    """Genomic sequence 3' of the mosaic end, in genomic orientation.

    Returns (flank, orientation) or (None, reason).  A read may present the
    cassette in either direction; taking the reverse complement first when the
    mosaic end appears as ME_RC puts the flank the same way round either way.
    """
    i = seq.find(ME)
    if i >= 0:
        return seq[i + len(ME):], "as-read"
    j = seq.find(ME_RC)
    if j >= 0:
        return revcomp(seq[:j]), "revcomp"
    return None, "no mosaic end in read"


with open(snakemake.input.sheet) as fh:  # noqa: F821
    rows = [r for r in csv.DictReader(
        (l for l in fh if not l.startswith("#")), delimiter="\t")]

n_written = 0
with open(snakemake.output.fasta, "w") as out:  # noqa: F821
    for r in rows:
        for direction, path in (("F", r["fwd"]), ("R", r["rev"])):
            _, seq = read_fasta(path)
            flank, how = flank_of(seq)
            if flank is None:
                print(f"{r['clone']}_{direction}: skipped - {how}", file=sys.stderr)
                continue
            m = re.search(r"N{10,}", flank)
            if m:
                flank = flank[:m.start()]
            if len(flank) < MIN_FLANK:
                print(f"{r['clone']}_{direction}: skipped - flank {len(flank)} bp "
                      f"< min_flank {MIN_FLANK}", file=sys.stderr)
                continue
            out.write(f">{r['clone']}|{direction}|{how}\n{flank}\n")
            n_written += 1
            print(f"{r['clone']}_{direction}: {len(flank)} bp flank ({how})",
                  file=sys.stderr)

print(f"\n{n_written} flank(s) from {len(rows)} clone(s)", file=sys.stderr)
if not n_written:
    sys.exit("no Sanger flank could be cut - is the cassette mosaic end correct?")
