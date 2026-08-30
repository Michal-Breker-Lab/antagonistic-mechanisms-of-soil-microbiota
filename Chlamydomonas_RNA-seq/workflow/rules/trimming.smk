rule fastp:
    input:
        sample=get_fastq,
    output:
        trimmed=temp("results/trimming/{sample}/{sample}.fastq.gz"),
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        adapters="",
        extra=config["extra"]["fastp"],
    threads: 6
    wrapper:
        "v9.0.0/bio/fastp"
