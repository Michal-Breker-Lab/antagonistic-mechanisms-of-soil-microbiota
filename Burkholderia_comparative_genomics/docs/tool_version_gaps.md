# Tool version gaps — rebuild (Moriah) vs original run (Shannon)

The original pipeline ran on Shannon in `envs/burk`, `envs/bakta`, `envs/macsy`.
Shannon went down mid-project; the rebuild runs on Moriah against the
environments available there. Where a version differs, it is recorded here and
must be stated in the Methods rather than silently inheriting the original
version numbers from the old report.

| Tool | Original (Shannon) | Rebuild (Moriah) | Direction | Notes |
|---|---|---|---|---|
| Panaroo | 1.8.0 | **1.7.0** | downgrade | Minor version. Clustering/graph-correction defaults can differ; pangenome sizes are not guaranteed identical to the original tables. |
| IQ-TREE | 3.1.3 | **3.1.1** | downgrade | Patch level only; model selection and inference unchanged in this range. |
| pyrodigal | 3.7.1 | **3.5.1** | downgrade | Gene calls validated indirectly: MF6 gives 6,902 CDS vs Bakta's 6,899. |
| skani | (unrecorded) | 0.3.2 | — | Clone clustering reproduces the original 256-cluster partition exactly, so any version difference is immaterial at ANI>=99 / AF>=50. |
| BLAST+ | 2.17.0+ | **2.15.0+** | downgrade | blastp defaults unchanged for the recorded command; `-comp_based_stats 2` present in both. |
| MAFFT | (unrecorded) | 7.526 | — | |
| Bakta | 1.12.0 / DB 6.0 | 1.12.0 / DB 6.0 | same | |

## Where a version gap can actually change a number

- **Panaroo 1.8.0 -> 1.7.0 is the only one that can move a reported result.**
  Pangenome gene counts (core / soft-core / shell / cloud) may differ from the
  original tables even on identical input. The rebuild's counts are the ones to
  report; do not mix them with numbers carried over from the old report.
- The BLAST downgrade does not affect the recorded search command, and the tier
  rule reproduces both retained hit tables at zero mismatches.
- skani, IQ-TREE and pyrodigal differences are below the level at which any
  reported figure changes.

## Additional environment note

Several tools in the relocated Moriah environments carry **stale
`<old-home>/...` shebangs** from when `miniforge3` was moved to
`<lab-filesystem>/`. Confirmed for tRNAscan-SE (broke bakta, 574 failed tasks) and
the pyrodigal launcher. Preflight in every sbatch therefore RUNS each tool
rather than testing `-x`, which is what caught the pyrodigal case before it
could waste an array.

## InterProScan 5.69-101.0 (original, Shannon) -> 5.76-107.0 (rebuild, Moriah)

Shannon is down, so fig7's domain annotation moves to the Moriah install at
`<tools>/interproscan-5.76-107.0/` (full 52 GB `data/`
present: antifam cdd funfam gene3d hamap ncbifam panther pfam phobius pirsf
pirsr prints prosite sfld smart superfamily tmhmm).

**This gap CAN move numbers**, unlike the BLAST and pyrodigal gaps. A minor
InterProScan release bumps the bundled member databases (Pfam, PANTHER, CDD,
...), so a protein can gain or lose a domain call between 5.69 and 5.76
independently of anything this project changed. The "did InterProScan annotate
this protein at all" counts that fig7 rests on are therefore not strictly
comparable to the original's.

Mitigation: report the version actually used, and treat any change in the
"no-domain" fraction as confounded by the database bump rather than as a
result of the 140 -> 253 genome expansion. If the difference matters to a
conclusion, re-run the ORIGINAL setB protein set through 5.76 as a control --
that isolates the database effect from the gene-set effect at negligible cost.
