rule impulseDE2_init:
    input:
        deseq2_rds=rules.deseq2_init.output.rds,
        counts=rules.genecounts_matrix.output[0],
    output:
        rds="results/impulseDE2/impulseDE2.rds",
    params:
        condition=config["deseq2"]["condition"],
        scaQThres=config["impulseDE2"]["scaQThres"],
    threads: 11
    log:
        "logs/impulseDE2/impulseDE2_init.log",
    conda:
        "../envs/impulseDE2.yaml"
    script:
        "../scripts/impulseDE2_init.R"


rule impulseDE2_subset_zscores:
    input:
        impulsede2_rds=rules.impulseDE2_init.output.rds,
        zscores=rules.deseq2_zscores.output.zscores,
    output:
        zscores_direction="results/impulseDE2/zscores_{direction}.tsv",
    params:
        direction=lambda wc: wc.direction,
        scaQThres=config["impulseDE2"]["scaQThres"],
    threads: 1
    log:
        "logs/impulseDE2/impulseDE2_subset_zscores_{direction}.log",
    conda:
        "../envs/impulseDE2.yaml"
    script:
        "../scripts/impulseDE2_zscores.R"


rule impulseDE2_results_table:
    input:
        impulsede2_rds=rules.impulseDE2_init.output.rds,
        all_annotations=rules.all_annotations.output.all_annotations,
    output:
        table="results/impulseDE2/impulseDE2_results.tsv",
    params:
        scaQThres=config["impulseDE2"]["scaQThres"],
    threads: 1
    log:
        "logs/impulseDE2/impulseDE2_results_table.log",
    conda:
        "../envs/impulseDE2.yaml"
    script:
        "../scripts/impulseDE2_results_table.R"
