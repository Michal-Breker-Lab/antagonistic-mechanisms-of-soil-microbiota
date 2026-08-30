rule star:
    input:
        fq1=rules.pigz_compress.output.fq,
        idx=rules.star_index.output[0],
    output:
        aln="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="results/qc/star/{sample}/Log.out",
        log_final="results/qc/star/{sample}/Log.final.out",
        transcriptome_bam=temp("results/star/{sample}/Aligned.toTranscriptome.out.bam"),
    log:
        "logs/star/{sample}.log",
    params:
        extra=config["extra"]["star"],
    threads: 6
    conda:
        "../envs/star.yaml"
    script:
        "../scripts/star.py"


rule index_aln_star:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/index_aligned_star/{sample}.log",
    params:
        extra="",
    threads: 6
    wrapper:
        "v9.0.0/bio/samtools/index"


rule calculate_expression:
    input:
        bam=rules.star.output.transcriptome_bam,
        reference=multiext(
            "resources/ref/rsem_index/reference",
            ".grp",
            ".ti",
            ".transcripts.fa",
            ".seq",
            ".idx.fa",
            ".n2g.idx.fa",
        ),
    output:
        genes_results="results/rsem/{sample}/{sample}.genes.results",
        isoforms_results="results/rsem/{sample}/{sample}.isoforms.results",
    params:
        extra=config["extra"]["rsem_calculate_expression"],
    log:
        "logs/rsem/calculate_expression/{sample}.log",
    threads: 6
    wrapper:
        "v9.0.0/bio/rsem/calculate-expression"


rule genecounts_matrix:
    input:
        expand("results/rsem/{sample}/{sample}.genes.results", sample=samples.index),
    output:
        "results/rsem/gene_counts.tsv",
    params:
        column_name="expected_count",
    threads: 1
    log:
        "logs/rsem/gene_counts.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/concat_rsem.py"


rule tpm_matrix:
    input:
        expand("results/rsem/{sample}/{sample}.genes.results", sample=samples.index),
    output:
        "results/rsem/gene_TPM.tsv",
    params:
        column_name="TPM",
    threads: 1
    log:
        "logs/rsem/gene_TPM.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/concat_rsem.py"
