rule deseq2_init:
    input:
        counts="results/rsem/gene_counts.tsv",
    output:
        rds="results/deseq2/all.rds",
    params:
        condition=config["deseq2"]["condition"],
        baselevel=config["deseq2"]["baselevel"],
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/deseq2_rds.log",
    script:
        "../scripts/deseq2_init.R"


rule dispersion_plot:
    input:
        rds=rules.deseq2_init.output.rds,
    output:
        dispersion_plot="results/qc/deseq2/dispersion-plot.svg",
        dispersion_stats="results/qc/deseq2/dispersion_stats.tsv",
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/dispersion-plot.log",
    script:
        "../scripts/deseq2_dispersion_plot.R"


rule deseq2_zscores:
    input:
        rds=rules.deseq2_init.output.rds,
    output:
        zscores="results/deseq2/z-scores.tsv",
    params:
        condition=config["deseq2"]["condition"],
        transformation=config["deseq2"]["transformation"],
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/z-scores.log",
    script:
        "../scripts/deseq2_z_scores.R"


rule deseq2_de:
    input:
        rds="results/deseq2/all.rds",
    output:
        long="results/deseq2/Ts_vs_T0.csv",
    params:
        condition=config["deseq2"]["condition"],
        baselevel=config["deseq2"]["baselevel"],
    log:
        "logs/deseq2/Ts_vs_T0.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_lfc.R"


rule annotate_merged_clusters:
    input:
        table="results/deseq2/lfc_all_vs_T0_merged_clusters.csv",
        all_annotations=rules.all_annotations.output.all_annotations,
    output:
        annotated="results/deseq2/lfc_annotated.tsv",
    threads: 1
    log:
        "logs/deseq2/lfc_annotated.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_annotations_to_table.py"
