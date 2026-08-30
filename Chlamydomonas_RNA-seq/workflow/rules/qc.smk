rule fastqc_raw:
    input:
        get_fastq,
    output:
        html="results/qc/fastqc_raw/{sample}.html",
        zip="results/qc/fastqc_raw/{sample}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_raw/{sample}.log",
    threads: 1
    resources:
        mem_mb=1024,
    wrapper:
        "v9.0.0/bio/fastqc"


rule fastqc_trimmed:
    input:
        rules.fastp.output.trimmed,
    output:
        html="results/qc/fastqc_trimmed/{sample}.html",
        zip="results/qc/fastqc_trimmed/{sample}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_trimmed/{sample}.log",
    threads: 1
    resources:
        mem_mb=1024,
    wrapper:
        "v9.0.0/bio/fastqc"


rule multiqc:
    input:
        multiqc_input,
        config=config["multiqc_config"],
    output:
        "results/qc/multiqc.html",
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    wrapper:
        "v9.0.0/bio/multiqc"
