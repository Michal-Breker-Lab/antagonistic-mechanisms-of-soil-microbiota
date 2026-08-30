rule gff2gtf:
    input:
        gff=config["ref_gff"],
    output:
        gtf="resources/ref/reference.gtf",
    log:
        "logs/ref/gff2gtf.log",
    conda:
        "../envs/gffread.yaml"
    shell:
        r"""
        gffread -T {input.gff:q} > {output.gtf:q} 2> {log}
        """


rule prepare_reference:
    input:
        reference_genome=config["ref_fasta"],
        gtf=rules.gff2gtf.output.gtf,
    output:
        seq="resources/ref/rsem_index/reference.seq",
        grp="resources/ref/rsem_index/reference.grp",
        ti="resources/ref/rsem_index/reference.ti",
        transcripts_fa="resources/ref/rsem_index/reference.transcripts.fa",
        idx_fa="resources/ref/rsem_index/reference.idx.fa",
        n2g_idx_fa="resources/ref/rsem_index/reference.n2g.idx.fa",
    params:
        extra=f"--gtf {rules.gff2gtf.output.gtf}",
    log:
        "logs/rsem/prepare-reference.log",
    wrapper:
        "v9.0.0/bio/rsem/prepare-reference"


rule star_index:
    input:
        fasta=config["ref_fasta"],
        gtf=rules.gff2gtf.output.gtf,
    output:
        directory("resources/ref/star_index"),
    threads: 6
    params:
        extra="--genomeSAindexNbases 12",
    log:
        "logs/star/star_index.log",
    wrapper:
        "v9.0.0/bio/star/index"


rule all_annotations:
    input:
        master_annotation=config["master_annotation"],
        annotation_info=config["annotation_info"],
    output:
        all_annotations="resources/ref/all_annotations.tsv",
    threads: 1
    log:
        "logs/ref/all_annotations.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/all_annotations.py"
