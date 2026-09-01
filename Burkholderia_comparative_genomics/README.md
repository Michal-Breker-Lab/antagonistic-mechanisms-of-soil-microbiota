# Burkholderia comparative genomics

Genus-wide comparative genomics of *Burkholderia* sensu lato, framing the lab isolate
*Burkholderia sola* MF6 against every complete public genome of the group: where its
third replicon (chromosome 3 / pC3) sits in the genus, what that replicon carries, and
how the candidate toxin and effector loci are distributed across relatives.

Unlike the two Snakemake workflows in this repository, this analysis ran as a staged
pipeline of shell, `sbatch` and Python steps across two compute hosts, because several
stages are multi-day jobs on 64–128 cores. The scripts are published **exactly as they
were executed**, including their hard-coded cluster paths, environment names and SLURM
directives: they are the record of what actually produced the trees and tables, not a
turnkey installer. Re-running them elsewhere means re-pointing the paths at the top of
each script.

## The two genome sets, and which figure uses which

The analysis produced **two independent phylogenies over overlapping genome sets**.
They are not interchangeable, and the counts differ accordingly:

| | Chromosome-1 core tree | bac120 marker tree |
| --- | --- | --- |
| Markers | Core genes of chromosome 1 (1,337,868 aligned nucleotide columns) | 120 GTDB single-copy markers (5,010 masked amino-acid columns) |
| Genomes in set | 773 (771 NCBI complete-level + MF6 + MF7) | 791 |
| Tips on the tree | 763 | 791 built, **754 drawn** |
| Chromosome 3 present | 252 / 763 (33.0%) | 282 / 754 (37.4%) |
| Inference | IQ-TREE, GTR+F+I+R5, 1000 UFBoot | IQ-TREE, LG+F+R4, 1000 UFBoot |
| Built by | `workflow/{shannon,moriah}/chr1_pC3/` | `workflow/{shannon,moriah}/bac120/` |
| Trees | `data/trees/chr1_core_763.*` | `data/trees/bac120_791.*` |

The 37-tip difference between 791 and 754 is a set of isogenic transposon mutants of a
single *B. sola* strain, R-12632 (BioProject PRJEB40633), each of which inherits the
parent's isolate metadata. They are excluded from the drawn phylogeny so that one strain
does not contribute 38 near-identical tips to the isolation-source and chromosome-3
distributions; the parental genome is retained in their place. No other dereplication
was performed on either tree — see `docs/analysis_decisions.md` (D1) for why clone
structure is handled as a covariate rather than a filter.

A third tree, over the pC3 core genes of the 253 carriers, supports the replicon-history
comparisons (`data/trees/pC3_core_253.*`).

## Layout

| Path | Contents |
| --- | --- |
| `config/` | The genome sets: the 773 accessions of the chromosome-1 analysis, their genus/species labels, and the bac120 accession list with its assembly metadata. |
| `resources/queries/` | Protein queries for the homology searches: the RHS toxin, the twelve candidate toxin/effector peptides, and the MF6 candidate set. |
| `workflow/shannon/chr1_pC3/` | Stages 1–30 of the chromosome-1 / pC3 pipeline as run on the workstation: metadata acquisition, download, replicon census and typing, dereplication, Bakta annotation, pangenomes, trees, host-vocabulary mapping, functional screens, homology searches. |
| `workflow/moriah/chr1_pC3/` | The same pipeline's SLURM arm — `m1`–`m17` job scripts plus `s23`–`s46` — used for the rebuild after the workstation became unavailable. `s46_section7_comparative.R` holds the phylogenetic comparative tests (`phylo.d`, `make.simmap`, `phyloglm`). |
| `workflow/shannon/bac120/` | The bac120 marker-tree pipeline: GTDB-Tk environment build, genome fetch, `identify`/`align`, FastTree, IQ-TREE, ANI, and the gene calls / toxin / pC3 screens for the genomes the chromosome-1 pipeline never annotated. |
| `workflow/moriah/bac120/` | Accession-set construction and the SLURM download / identify / align jobs. |
| `workflow/figures/` | Figure scripts and the statistics behind them, run locally against `data/`. Also the derived-table builders (`s33`–`s46`) and the workbook generator. |
| `workflow/drivers/` | Orchestration run from the workstation: staging to a remote host, pulling results back, job watchers, and the script-mirroring helper. |
| `data/trees/` | The three phylogenies, in IQ-TREE `.treefile` (ML) and `.contree` (majority-rule consensus with UFBoot support) form, plus the FastTree approximations used for the congruence check. |
| `data/tables/` | The small result tables the figure scripts read — pC3 calls, replicon census, host categories, ANI to MF6, clone clusters, homology-search hits, functional and screen summaries. |
| `docs/` | Analysis decisions D1–D20 with their reasoning, per-tool provenance (versions, flags, papers, caveats), and the version gaps between the original run and the rebuild. |

## Reproducing the figures

The figure scripts read from `data/trees/` and `data/tables/` and need only a
scientific Python stack (matplotlib, numpy, pandas, ete3/Bio.Phylo). `figstyle.py`
enforces the project's vector-text rule — PDF/SVG output keeps live, editable text
rather than outlined glyphs — and should be imported before any figure is drawn.

```bash
cd Burkholderia_comparative_genomics
python workflow/figures/fig13_genus_tree_bac120.py    # the genus-wide marker tree
python workflow/figures/fig1_species_tree.py          # chromosome-1 tree, pC3 + host
```

Scripts resolve their inputs relative to the analysis directory they ran in; adjust the
path constants at the top of each if you run them from elsewhere.

## What is not here

Genome FASTAs, Bakta and pyrodigal annotations, reference databases (GTDB r232, Bakta
DB 6.0, InterProScan, antiSMASH), the pangenome presence/absence matrices and the
concatenated alignments. All of these are either public and re-downloadable from the
accessions in `config/`, or regenerate from the scripts here; the largest of them run to
hundreds of gigabytes. `docs/tool_provenance.md` records every version and flag needed
to reproduce them.

Two earlier states of this analysis are also not published: a 304-tip dereplicated
version of the chromosome-1 tree, superseded when dereplication was replaced by
clone-cluster covariates, and the pre-MF7 771-genome set. Where a decision reversed an
earlier one, `docs/analysis_decisions.md` records both.
