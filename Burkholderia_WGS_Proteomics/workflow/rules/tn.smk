# Tn5 insertion verification: Tables S12 and S13.
#
# Ported from snakemake_insertion_verification_chr3, minus its assembly branch
# (flye / dnadiff / assembly_check) - neither table depends on it.
#
# Reads come from SRA by accession, so nothing here points at a path outside the
# workflow.  The annotation is res("gff"), the SAME staged file the proteomics
# rules read, so locus tags and replicon names cannot diverge between the two
# halves of this pipeline.


rule fetch_reads_nanopore:
    """One nanopore SRA run -> a single gzipped FASTQ."""
    output:
        fq=temp("results/Tn/reads/{clone}.fastq.gz"),
    log:
        "logs/tn/fetch_reads/{clone}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 6
    params:
        sra=sra_of,
        outdir=reads_dir,
    shell:
        r"""
        (
            set -e
            fasterq-dump --threads {threads} --split-3 \
                --outdir {params.outdir:q} {params.sra}
            pigz -p {threads} -f {params.outdir:q}/{params.sra}.fastq
            mv {params.outdir:q}/{params.sra}.fastq.gz {output.fq:q}
        ) >{log} 2>&1
        """


rule fetch_reads_illumina:
    """One Illumina SRA run -> the two mate files."""
    output:
        r1=temp("results/Tn/reads/{clone}_1.fastq.gz"),
        r2=temp("results/Tn/reads/{clone}_2.fastq.gz"),
    log:
        "logs/tn/fetch_reads/{clone}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 6
    params:
        sra=sra_of,
        outdir=reads_dir,
    shell:
        r"""
        (
            set -e
            fasterq-dump --threads {threads} --split-3 \
                --outdir {params.outdir:q} {params.sra}
            pigz -p {threads} -f {params.outdir:q}/{params.sra}_1.fastq \
                {params.outdir:q}/{params.sra}_2.fastq
            mv {params.outdir:q}/{params.sra}_1.fastq.gz {output.r1:q}
            mv {params.outdir:q}/{params.sra}_2.fastq.gz {output.r2:q}
        ) >{log} 2>&1
        """


rule build_reference:
    """Genome plus the cassette as a fourth contig.

    Reads that cross a junction align to `transposon` instead of being clipped
    away, which is what makes the junction callable at all.
    """
    input:
        genome=res("genome"),
        tn5=res("tn5"),
    output:
        fasta="results/Tn/ref/MF6_plus_Tn5.fa",
    log:
        "logs/tn/build_reference.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        # awk 1, not cat: the genome FASTA has no trailing newline and cat would
        # fuse the Tn5 header onto the last contig's final sequence line.
        awk 1 {input.genome:q} {input.tn5:q} >{output.fasta:q} 2>{log}
        """


rule faidx_reference:
    input:
        rules.build_reference.output.fasta,
    output:
        "results/Tn/ref/MF6_plus_Tn5.fa.fai",
    log:
        "logs/tn/faidx.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input:q} 2> {log}"


rule bwa_index:
    input:
        rules.build_reference.output.fasta,
    output:
        multiext(
            "results/Tn/ref/MF6_plus_Tn5.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/tn/bwa_index.log",
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input:q} 2> {log}"


rule minimap2_nanopore:
    """Nanopore alignment.

    minimap2 2.30 and samtools 1.19.2 pin incompatible libzlib and cannot share
    an environment, so the SAM is written here and sorted by the next rule.
    """
    input:
        fq=nanopore_fastq,
        ref=rules.build_reference.output.fasta,
    output:
        sam=temp("results/Tn/bam/{clone}.sam"),
    log:
        "logs/tn/minimap2/{clone}.log",
    conda:
        "../envs/minimap2.yaml"
    threads: 12
    shell:
        r"""
        minimap2 -ax map-ont --MD -Y -L --secondary=no -t {threads} \
            {input.ref:q} {input.fq:q} >{output.sam:q} 2>{log}
        """


rule sort_nanopore:
    input:
        sam=rules.minimap2_nanopore.output.sam,
    output:
        bam="results/Tn/bam/{clone}.bam",
    log:
        "logs/tn/sort_nanopore/{clone}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools sort -@ {threads} -m 1G -o {output.bam:q} {input.sam:q} 2> {log}"


rule bwa_mem_illumina:
    """Illumina alignment, duplicate-marked.

    bwa mem batches 10M x nthreads bases and estimates the insert size per batch,
    so the thread count changes the alignments.  It is pinned in
    workflow/profiles/default rather than left to whatever --cores happens to be.
    """
    input:
        fq=illumina_fastq,
        ref=rules.build_reference.output.fasta,
        idx=rules.bwa_index.output,
    output:
        bam="results/Tn/bam_illumina/{clone}.bam",
    log:
        "logs/tn/bwa/{clone}.log",
    conda:
        "../envs/bwa.yaml"
    threads: 12
    shell:
        r"""
        bwa mem -t {threads} -Y {input.ref:q} {input.fq:q} 2>{log} \
            | samtools fixmate -u -m - - \
            | samtools sort -u -@ 4 -m 1G - \
            | samtools markdup -@ 4 - {output.bam:q} 2>>{log}
        """


rule samtools_index:
    input:
        "results/Tn/{bamdir}/{clone}.bam",
    output:
        "results/Tn/{bamdir}/{clone}.bam.bai",
    log:
        "logs/tn/samtools_index/{bamdir}/{clone}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input:q} 2> {log}"


rule find_insertions:
    """Call insertions per clone from junction and spanning-read evidence."""
    input:
        bams=platform_bams,
        bais=platform_bais,
        gff=res("gff"),
        code=deps("find_insertions.py"),
    output:
        insertions="results/Tn/tables_{platform}/insertions.tsv",
        insertions_all="results/Tn/tables_{platform}/insertions_all.tsv",
        inserted_seqs="results/Tn/tables_{platform}/inserted_seqs.fasta",
        clone_summary="results/Tn/tables_{platform}/clone_summary.tsv",
        coverage="results/Tn/tables_{platform}/coverage_per_contig.tsv",
    log:
        "logs/tn/find_insertions/{platform}.log",
    conda:
        "../envs/pysam.yaml"
    params:
        platform=platform_of,
        detect=tn_detect,
        absent=config["absent_replicon"],
    script:
        "../scripts/find_insertions.py"


rule combine_platforms:
    input:
        insertions=expand(
            "results/Tn/tables_{platform}/insertions.tsv", platform=TN_PLATFORMS
        ),
        summaries=expand(
            "results/Tn/tables_{platform}/clone_summary.tsv", platform=TN_PLATFORMS
        ),
        code=deps("combine_platforms.py"),
    output:
        sites=temp("results/Tn/tables/insertions_combined.tsv"),
        clones="results/Tn/tables/clones_complete.tsv",
    log:
        "logs/tn/combine_platforms.log",
    conda:
        "../envs/python.yaml"
    params:
        platforms=TN_PLATFORMS,
    script:
        "../scripts/combine_platforms.py"


rule cds_notation:
    """HGVS-style CDS notation and nearest-gene columns."""
    input:
        sites=rules.combine_platforms.output.sites,
        gff=res("gff"),
        code=deps("cds_notation.py"),
    output:
        sites="results/Tn/tables/insertions_complete.tsv",
    log:
        "logs/tn/cds_notation.log",
    conda:
        "../envs/python.yaml"
    params:
        cassette_len=CASSETTE_LEN,
    script:
        "../scripts/cds_notation.py"


# --------------------------------------------------------------- Sanger clones


rule sanger_flanks:
    input:
        sheet=config["sanger"],
        reads=SANGER_READS,
        code=deps("sanger_flanks.py"),
    output:
        fasta="results/Tn/sanger/flanks.fasta",
    log:
        "logs/tn/sanger_flanks.log",
    conda:
        "../envs/python.yaml"
    params:
        min_flank=config["tn"]["sanger"]["min_flank"],
    script:
        "../scripts/sanger_flanks.py"


rule sanger_map:
    input:
        fasta=rules.sanger_flanks.output.fasta,
        ref=rules.build_reference.output.fasta,
    output:
        sam=temp("results/Tn/sanger/flanks.sam"),
    log:
        "logs/tn/sanger_map.log",
    conda:
        "../envs/minimap2.yaml"
    threads: 4
    shell:
        r"""
        minimap2 -ax map-ont --MD -Y -L --secondary=no -t {threads} \
            {input.ref:q} {input.fasta:q} >{output.sam:q} 2>{log}
        """


rule sanger_sort:
    input:
        sam=rules.sanger_map.output.sam,
    output:
        bam="results/Tn/sanger/flanks.bam",
        bai="results/Tn/sanger/flanks.bam.bai",
    log:
        "logs/tn/sanger_sort.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        r"""
        samtools sort -@ {threads} -o {output.bam:q} {input.sam:q} 2>{log}
        samtools index -@ {threads} {output.bam:q} 2>>{log}
        """


rule sanger_sites:
    input:
        bam=rules.sanger_sort.output.bam,
        bai=rules.sanger_sort.output.bai,
        gff=res("gff"),
        wgs=rules.cds_notation.output.sites,
        code=deps("sanger_sites.py"),
    output:
        sites="results/Tn/tables/insertions_sanger.tsv",
    log:
        "logs/tn/sanger_sites.log",
    conda:
        "../envs/pysam.yaml"
    params:
        min_mapq=config["tn"]["sanger"]["min_mapq"],
        cassette=CASSETTE_NAME,
    script:
        "../scripts/sanger_sites.py"


# -------------------------------------------------------------- coverage depth


rule coverage_breadth_raw:
    input:
        bams=all_tn_bams,
        bais=all_tn_bais,
        script="workflow/scripts/coverage_breadth.sh",
    output:
        table="results/Tn/tables/coverage_breadth_raw.tsv",
    log:
        "logs/tn/coverage_breadth_raw.log",
    conda:
        "../envs/samtools.yaml"
    params:
        mapq=tn_mapq_csv,
        specs=tn_bam_specs,
    shell:
        "bash {input.script:q} {output.table:q} {params.mapq} {params.specs} 2> {log}"


rule coverage_breadth:
    input:
        table=rules.coverage_breadth_raw.output.table,
        code=deps("coverage_breadth_chr.py"),
    output:
        table="results/Tn/tables/coverage_breadth.tsv",
    log:
        "logs/tn/coverage_breadth.log",
    conda:
        "../envs/python.yaml"
    params:
        cassette=CASSETTE_NAME,
    script:
        "../scripts/coverage_breadth_chr.py"


# ------------------------------------------------------------ published tables


rule table_tn_insertion_sites:
    input:
        wgs=rules.cds_notation.output.sites,
        sanger=rules.sanger_sites.output.sites,
        gff=res("gff"),
        code=deps("table_tn_insertion_sites.py"),
    output:
        tsv="results/supp_tables/Table_06_Tn_insertion_sites.tsv",
        xlsx=published("tn_insertion_sites"),
    log:
        "logs/tn/table_tn_insertion_sites.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/table_tn_insertion_sites.py"


rule table_coverage_breadth:
    input:
        table=rules.coverage_breadth.output.table,
        code=deps("table_coverage_breadth.py"),
    output:
        tsv="results/supp_tables/Table_07_coverage_breadth.tsv",
        xlsx=published("coverage_breadth"),
    log:
        "logs/tn/table_coverage_breadth.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/table_coverage_breadth.py"
