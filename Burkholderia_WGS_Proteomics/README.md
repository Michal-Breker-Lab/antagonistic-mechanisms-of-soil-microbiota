# Proteomics analysis for MF6

A reproducible [Snakemake](https://snakemake.github.io/) workflow for the label-free proteomics presented in the MF6 paper. The pipeline takes FragPipe protein quantifications for *Burkholderia sola* MF6 and two transposon mutants, grown alone and in co-culture with *Chlamydomonas reinhardtii*, through normalisation, differential abundance and presence/absence calling, clustering and functional enrichment, and adds secretion-system inventories and Tn5 insertion verification, producing the supplementary tables and publication figures end-to-end.

## Experimental design

Three strains (MF6 and the transposon mutants 27D6 and 34F7) were grown alone or with *C. reinhardtii*, sampled on days 2 and 3, in four biological replicates — 48 runs in total (`config/samples.tsv`, two of which fail QC and are excluded by the `qc_pass` column). Models are fitted on day 2, the only day with mutants.

Proteins are identified against the closed three-replicon MF6 genome (GenBank [GCA_056856995.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_056856995.2/), `AC1V0C_*` locus tags on `chr1`/`chr2`/`chr3`); both mutants have lost `chr3`, which is censored in their columns.

Two side analyses share the same annotation: secretion systems are called in MF6 and seven IMG co-isolate genomes (`config/config.yaml`), and Tn5 insertions are verified from whole-genome nanopore and Illumina runs fetched from SRA for the clones listed in `config/clones.tsv`, plus clones with Sanger evidence only (`config/sanger.tsv`).

## Workflow

```
FragPipe combined_protein.tsv
  → MaxLFQ matrix, core-set median centring
  → limma          (differential abundance, day 2, seven contrasts)
  → ON/OFF calling (present in one group, absent in the other; percentile-gated)
  → hierarchical clustering of the core set
  → COG / KEGG over-representation (BH-corrected)

MF6 + 7 relative genomes
  → MacSyFinder / TXSScan  (secretion systems, pili, flagella)

SRA runs + Sanger reads
  → minimap2 / BWA-MEM → insertion calling → coverage breadth

  → supplementary tables (xlsx) + volcano, PCA, heatmap and bubble figures
```

## Requirements

- [Snakemake](https://snakemake.github.io/) (≥ 9.12; developed on 9.24)
- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) — per-rule software environments are defined under `workflow/envs/` and resolved automatically.

Network access is needed on the first run: to build the conda environments, to install the TXSScan model set, to fetch the Tn5 sequencing runs from SRA, and to fetch the KEGG BRITE hierarchy.

## Usage

All inputs are already staged under `resources/` and declared in `config/config.yaml`. To run everything:

```bash
snakemake
```

`workflow/profiles/default/` is picked up automatically and supplies `--sdm conda` and `--cores 12`; the aligners are pinned to 12 threads there because their thread count changes the alignments. Override with `--cores <N>` if you accept that.

All analysis parameters (DE thresholds, the ON/OFF percentile gate, cluster count, KEGG filtering, insertion-calling presets, published table names) are exposed in `config/config.yaml`.

## Repository structure

```
config/      sample sheets (proteomics runs, Tn5 clones, Sanger reads) and analysis configuration
workflow/    Snakefile, rules (*.smk), scripts, conda environments, schemas, and the run profile
resources/   genome, annotation, FragPipe quantifications, and the other staged inputs
results/     generated matrices, statistics, supplementary tables, and figures
```

Finished workbooks are written straight to `results/published/` as `Table_XX…Table_XX`.

## Citation

If you use this workflow, please cite the MF6 paper (citation to be added upon publication).
