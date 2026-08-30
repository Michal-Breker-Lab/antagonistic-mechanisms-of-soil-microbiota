# Antagonistic mechanisms of soil microbiota

Analysis code for the paper on antagonistic mechanisms of soil microbiota, covering the interaction between the soil bacterium *Burkholderia sola* MF6 and the green alga *Chlamydomonas reinhardtii*. Both analyses are [Snakemake](https://snakemake.github.io/) workflows that run end-to-end from the staged inputs to the publication figures and supplementary tables.

## Contents

| Directory | Analysis |
| --- | --- |
| [`Burkholderia_WGS_Proteomics/`](Burkholderia_WGS_Proteomics/) | Label-free proteomics of MF6 and two transposon mutants, alone and in co-culture: differential abundance, ON/OFF calling, clustering and COG/KEGG enrichment, plus secretion-system inventories (TXSScan) and Tn5 insertion verification from WGS and Sanger reads. |
| [`Chlamydomonas_RNA-seq/`](Chlamydomonas_RNA-seq/) | Time-course RNA-seq of *C. reinhardtii* (T0–T4, three replicates): QC, rRNA depletion, STAR/RSEM quantification, DESeq2, ImpulseDE2, Mfuzz clustering and topGO enrichment. |

Each directory has its own README with the experimental design, workflow diagram, and usage instructions.

## Requirements

- [Snakemake](https://snakemake.github.io/) — ≥ 9.12 for the proteomics workflow, ≥ 8 for the RNA-seq workflow
- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) — per-rule software environments are defined under each workflow's `workflow/envs/` and resolved automatically

Run either workflow from its own directory:

```bash
cd Burkholderia_WGS_Proteomics && snakemake        # profile supplies --sdm conda --cores 12
cd Chlamydomonas_RNA-seq   && snakemake --sdm conda --cores <N>
```

## Data

Code, configuration and small staged inputs are tracked here. Generated outputs (`results/`, `logs/`) and large reference and raw-read data are not; the RNA-seq workflow expects its FASTQ files and reference genome to be placed under `resources/` as described in its README, and the proteomics workflow fetches the Tn5 sequencing runs from SRA on first run.

## Citation

If you use this code, please cite the MF6 paper (citation to be added upon publication).
