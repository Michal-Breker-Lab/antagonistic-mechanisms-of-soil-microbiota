rule mfuzz:
    input:
        zscores="results/impulseDE2/zscores_{direction}.tsv",
    output:
        centers="results/Mfuzz/{direction}_centers.tsv",
        clusters="results/Mfuzz/{direction}_clusters.tsv",
        membership="results/Mfuzz/{direction}_membership.tsv",
    params:
        seed=config["seed"],
        Nclusters=config["Mfuzz"]["Nclusters"],
    threads: 1
    log:
        "logs/Mfuzz/{direction}.log",
    conda:
        "../envs/Mfuzz.yaml"
    script:
        "../scripts/Mfuzz.R"


rule add_info:
    input:
        clusters="results/Mfuzz/{direction}_clusters.tsv",
        membership="results/Mfuzz/{direction}_membership.tsv",
        all_annotations=rules.all_annotations.output.all_annotations,
    output:
        clusters_annotated="results/Mfuzz/{direction}_annotated.tsv",
    params:
        rename_clusters=config["Mfuzz"]["rename_clusters"],
    threads: 1
    log:
        "logs/Mfuzz/{direction}_add_info.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/add_info_to_clusters.py"
