#!/usr/bin/env python3
"""Subset a Bakta GFF3 to a chosen set of contigs, keeping the file valid.

Bakta emits a combined GFF3: directives, feature lines, then '##FASTA' followed
by the genome sequence. Panaroo needs both halves to agree -- dropping feature
lines while leaving the full FASTA (or vice versa) produces errors that are
confusing to trace back. This subsets directives, features and sequences
together.

usage: subset_gff.py in.gff3 out.gff3 contig [contig ...]
"""
import sys


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    src, dst, keep = sys.argv[1], sys.argv[2], set(sys.argv[3:])

    n_feat = n_seq = 0
    with open(src) as fh, open(dst, "w") as out:
        in_fasta, writing = False, False
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                out.write(line)
                continue
            if not in_fasta:
                if line.startswith("##sequence-region"):
                    p = line.split()
                    if len(p) >= 2 and p[1] in keep:
                        out.write(line)
                    continue
                if line.startswith("#"):
                    out.write(line)
                    continue
                f = line.split("\t", 1)
                if f[0] in keep:
                    out.write(line)
                    n_feat += 1
            else:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    writing = name in keep
                    if writing:
                        n_seq += 1
                if writing:
                    out.write(line)
    if n_seq != len(keep):
        print(f"WARN {src}: wanted {len(keep)} contigs, wrote {n_seq} sequences",
              file=sys.stderr)
    if n_feat == 0:
        print(f"WARN {src}: no features retained", file=sys.stderr)
    print(f"{dst}\t{n_feat}\t{n_seq}")


if __name__ == "__main__":
    main()
