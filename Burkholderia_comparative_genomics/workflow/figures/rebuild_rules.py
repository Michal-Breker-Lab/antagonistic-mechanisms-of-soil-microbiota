"""Shared exclusion rules for the rebuild figures.

DRAFT_GENOMES -- genomes whose assembly is a fragmented draft rather than a
closed genome. Their CONTIG sizes are not REPLICON sizes, so they must never
enter a size distribution, a size-by-rank plot, a replicon-count, or an
architecture class. Decision D10 (approved 2026-08-21) makes this binding.

MF7 is the only such genome in the set: the other 772 are complete. Its census
row shows exactly why -- size-rank typing of its 29-contig draft yields
`architecture=4+_large` and `n_large_300kb=8`, when the genome actually has
three large replicons (chr1, chr2, pC3), and `total_length=8,695,339` is
inflated by roughly 940 kb of real tandem duplication on chromosome 2.

MF7's true architecture (3 large, pC3-positive) was established by alignment to
MF6 -- a different method from the size-rank typing applied to every other
genome -- so it is reported as a footnote, never as a counted data point.

MF7 is NOT dropped from the study: it stays in the 773 for gene-content work
(toxin/RHS searches, pangenome membership), where contig fragmentation does not
distort the measurement. This module governs size and architecture only.
"""

DRAFT_GENOMES = {"MF7"}

FOOTNOTE = (
    "MF7 is excluded from replicon-size and architecture panels: it is a "
    "29-contig draft, so its contig sizes are not replicon sizes. Its "
    "architecture (3 large replicons, pC3-positive) was assigned by alignment "
    "to MF6, of which it is a 100%-ANI clone."
)


def drop_drafts(rows, key="accession"):
    """Filter census-like rows, returning (kept, n_dropped)."""
    kept = [r for r in rows if r.get(key) not in DRAFT_GENOMES]
    return kept, len(rows) - len(kept)
