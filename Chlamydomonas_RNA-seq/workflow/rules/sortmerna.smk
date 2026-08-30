rule sortmerna_build_index:
    input:
        ref=config["sortmernaDB"],
    output:
        idx=directory("resources/sortmerna_idx"),
    threads: 1
    log:
        "logs/sortmerna/build_idx.log",
    conda:
        "../envs/sortmerna.yaml"
    script:
        "../scripts/sortmerna_build_index.py"


rule pigz_decompres:
    input:
        fq=rules.fastp.output.trimmed,
    output:
        fq=temp("results/pigz_decompres/{sample}.fastq"),
    threads: 6
    log:
        "logs/pigz_decompres/{sample}.log",
    conda:
        "../envs/pigz.yaml"
    shell:
        "pigz -p {threads} -c -d {input.fq} > {output.fq} 2> {log}"


rule sortmerna:
    input:
        ref=config["sortmernaDB"],
        reads=rules.pigz_decompres.output.fq,
        idx=rules.sortmerna_build_index.output.idx,
    output:
        aligned=temp("results/sortmerna/{sample}/aligned_{sample}.fastq"),
        other=temp("results/sortmerna/{sample}/unpaired_{sample}.fastq"),
        stats="results/qc/sortmerna/{sample}.log",
    params:
        extra="--zip-out",
    threads: 12
    log:
        "logs/sortmerna/{sample}.log",
    conda:
        "../envs/sortmerna.yaml"
    script:
        "../scripts/sortmerna.py"


rule pigz_compress:
    input:
        fq=rules.sortmerna.output.other,
    output:
        fq="results/sortmerna/{sample}.fastq.gz",
    threads: 6
    log:
        "logs/pigz_compress/{sample}.log",
    conda:
        "../envs/pigz.yaml"
    shell:
        "pigz -p {threads} -c {input.fq} > {output.fq} 2> {log}"
