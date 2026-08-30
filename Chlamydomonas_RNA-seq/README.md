# RNA-seq analysis for MF6

A reproducible [Snakemake](https://snakemake.github.io/) workflow for the time-course RNA-seq analysis presented in the MF6 paper. The pipeline processes raw single-end reads from *Chlamydomonas reinhardtii* through quality control, quantification, differential and temporal expression analysis, soft clustering, and functional enrichment, producing the publication figures end-to-end.

## Experimental design

Cells were sampled across five time points (T0–T4) in three biological replicates (A, B, C), for a total of 15 libraries (`config/samples.csv`). Reads were aligned to the *C. reinhardtii* CC-4532 v6.1 reference genome ([Phytozome, genome ID 707](https://phytozome-next.jgi.doe.gov/)).

## Workflow

```
raw reads
  → fastp (adapter/quality trimming) → FastQC
  → SortMeRNA (rRNA depletion)
  → STAR (alignment) + RSEM (transcript quantification)
  → DESeq2        (differential expression, rlog, PCA)
  → ImpulseDE2    (temporal differential expression)
  → Mfuzz         (soft clustering of expression trajectories)
  → topGO         (GO term enrichment per cluster)
  → figures + MultiQC report
```

## Requirements

- [Snakemake](https://snakemake.github.io/) (≥ 8)
- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) — per-rule software environments are defined under `workflow/envs/` and resolved automatically.

## Usage

Place the input FASTQ files in `resources/fastq/` and the reference genome/annotation in `resources/ref/` as specified in `config/config.yaml`, then run:

```bash
snakemake --sdm conda --cores <N>
```

All analysis parameters (reference paths, DESeq2 design, ImpulseDE2 thresholds, number of Mfuzz clusters, enrichment settings) are exposed in `config/config.yaml`.

## Repository structure

```
config/      sample sheet and analysis configuration
workflow/    Snakefile, rules (*.smk), scripts, and conda environments
resources/   reference genome, annotation, and input reads
results/     generated quantifications, statistics, and figures
```

## Citation

If you use this workflow, please cite the MF6 paper (citation to be added upon publication).
