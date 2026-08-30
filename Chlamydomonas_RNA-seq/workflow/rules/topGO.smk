rule topgo_clusters:
    input:
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        all_annotations="resources/ref/all_annotations.tsv",
        deseq2_rds=rules.deseq2_init.output.rds,
    output:
        enrichment_table="results/topGO/enrichment_by_clusters.tsv",
    params:
        minMembership=config["Enrichment"]["minMembership"],
        minNodeSize=config["Enrichment"]["minNodeSize"],
        minGeneCount=config["Enrichment"]["minGeneCount"],
    threads: 1
    log:
        "logs/topGO/enrichment_by_clusters.log",
    conda:
        "../envs/topGO.yaml"
    script:
        "../scripts/topGO.R"


rule merge_up1_down6:
    input:
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
    output:
        merged_table="results/Mfuzz/up1_down6.tsv",
    threads: 1
    log:
        "logs/topGO/merge_up1_down6.log",
    params:
        up_clusters=[1],
        down_clusters=[6],
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/merge_selected_clusters.py"


rule topgo_up1_down6:
    input:
        single_file="results/Mfuzz/up1_down6.tsv",
        all_annotations="resources/ref/all_annotations.tsv",
        deseq2_rds=rules.deseq2_init.output.rds,
    output:
        enrichment_table="results/topGO/up1_down6.tsv",
    params:
        minMembership=0,
        minNodeSize=config["Enrichment"]["minNodeSize"],
        minGeneCount=config["Enrichment"]["minGeneCount"],
    threads: 1
    log:
        "logs/topGO/up1_down6.log",
    conda:
        "../envs/topGO.yaml"
    script:
        "../scripts/topGO.R"


# rule clusters2direction:
#     input:
#         clusters="results/Mfuzz/{direction}_annotated.tsv",
#     output:
#         clusters="results/topGO/just_{direction}.tsv"
#     params:
#         use_clusters=lambda wc: config["Enrichment"][f"use_{wc.direction}_clusters"],
#     threads: 1
#     log:
#         "logs/topGO/clusters2direction_{direction}.log"
#     conda:
#         "../envs/pandas.yaml"
#     script:
#         "../scripts/clusters2direction.py"
# rule topgo_direction:
#     input:
#         up_clusters="results/topGO/just_up.tsv",
#         down_clusters="results/topGO/just_down.tsv",
#         all_annotations="resources/ref/all_annotations.tsv",
#         deseq2_rds=rules.deseq2_init.output.rds,
#     output:
#         enrichment_table="results/topGO/enrichment_by_direction.tsv"
#     params:
#         minMembership=0,
#         minNodeSize=config["Enrichment"]["minNodeSize"],
#         minGeneCount=config["Enrichment"]["minGeneCount"],
#     threads: 1
#     log:
#         "logs/topGO/enrichment_by_direction.log"
#     conda:
#         "../envs/topGO.yaml"
#     script:
#         "../scripts/topGO.R"
