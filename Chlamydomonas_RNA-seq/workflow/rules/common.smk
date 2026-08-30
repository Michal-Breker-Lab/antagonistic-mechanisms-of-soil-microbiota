import pandas as pd
from itertools import product

samples = pd.read_csv(config["samples"], dtype=str, comment="#").set_index(
    "sample_name", drop=False
)

impulsede2_directions = ["up", "down"]


wildcard_constraints:
    sample="|".join(samples["sample_name"].tolist()),
    subset="|".join(impulsede2_directions),


def get_fastq(wildcards):
    return (samples.loc[wildcards.sample, "fastq_path"],)


def get_mfuzz_clusters(wildcards):
    if wildcards.subset.startswith("transient"):
        return config["Mfuzz"]["transient_nclusters"]
    if wildcards.subset.startswith("transition"):
        return config["Mfuzz"]["transition_nclusters"]


def multiqc_input(wildcards):
    output = []

    output.extend(
        expand("results/qc/fastqc_raw/{sample}_fastqc.zip", sample=samples.index)
    )
    output.extend(expand("results/qc/fastp/{sample}.json", sample=samples.index))
    output.extend(
        expand("results/qc/fastqc_trimmed/{sample}_fastqc.zip", sample=samples.index)
    )
    output.extend(expand("results/qc/sortmerna/{sample}.log", sample=samples.index))
    output.extend(
        expand("results/qc/star/{sample}/Log.final.out", sample=samples.index)
    )

    return output


def get_final_output(wildcards):
    output = []
    output.extend(["results/figures/pca-plot.svg"])
    output.extend(["results/qc/deseq2/dispersion-plot.svg"])

    output.extend(
        expand(
            "results/figures/{direction}_expression_enrich.svg",
            direction=impulsede2_directions,
        )
    )

    output.extend(
        [
            "results/figures/heatmap_metal_transporters.svg",
            "results/figures/heatmap_manual_list.svg",
            "results/figures/heatmap_PCD.svg",
            "results/figures/heatmap_autophagy.svg",
            "results/figures/heatmap_metal_transporters_horizontal.pdf",
        ]
    )

    output.extend(
        [
            "results/figures/enrichment_plot_all.svg",
            "results/figures/enrichment_plot_up1_down6.svg",
        ]
    )

    output.extend(["results/qc/multiqc.html"])
    output.extend(["results/deseq2/lfc_annotated.tsv"])
    output.extend(["results/impulseDE2/impulseDE2_results.tsv"])

    return output
