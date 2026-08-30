rule split_matrices:
    input:
        combined=res("fragpipe"),
        code=deps("split_matrices.py"),
    output:
        maxlfq="results/Proteomics/MF6_maxlfq_intensity.tsv",
    log:
        "logs/split_matrices.log",
    conda:
        "../envs/python.yaml"
    params:
        locus_prefix=config["locus_prefix"],
    script:
        "../scripts/split_matrices.py"
