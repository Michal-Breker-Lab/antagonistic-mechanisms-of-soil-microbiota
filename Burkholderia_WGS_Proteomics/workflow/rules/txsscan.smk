# TXSScan: secretion systems, pili and flagella, for MF6 and its relatives.
#
# Two run modes, both `--models TXSScan all` so every model is searched in one
# pass, and both reproducing the invocation recorded in the existing output:
#
#   relatives  one run per genome over its whole proteome, --db-type gembase,
#              --replicon-topology linear (draft assemblies, many contigs)
#   MF6        one run PER REPLICON, --db-type ordered_replicon,
#              --replicon-topology circular (a closed genome, three circular
#              chromosomes), merged afterwards by macsy_merge_results
#
# MF6's proteomes are staged already relabelled into the submission id space
# (AC1V0C_*, chr1/chr2/chr3), so the systems this produces are named the way every
# other table in this workflow names them, with no post-hoc rewrite.


rule mf6_replicon_proteomes:
    """MF6's proteome split by replicon, from the staged genome annotation."""
    input:
        gff=res("gff"),
        faa=res("proteome"),
        code=deps("mf6_replicon_proteomes.py"),
    output:
        expand("results/txsscan/MF6/proteomes/{replicon}.faa", replicon=MF6_REPLICONS),
    log:
        "logs/txsscan/mf6_replicon_proteomes.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mf6_replicon_proteomes.py"


rule txsscan_models:
    """Fetch the TXSScan model set at the pinned version.

    macsydata is MacSyFinder's own package manager, so what lands here is exactly
    what `--models-dir` expects - the same tree this workflow used to store, which
    was checked file-for-file against this download.

    The version is asserted after installing rather than assumed: macsydata writes
    the installed version into metadata.yml, and a mismatch means the tag moved.
    """
    output:
        models=directory("results/txsscan/models/TXSScan"),
    log:
        "logs/txsscan/models.log",
    conda:
        "../envs/macsyfinder.yaml"
    params:
        org=config["secretion_systems"]["models"]["org"],
        package=config["secretion_systems"]["models"]["package"],
        version=config["secretion_systems"]["models"]["version"],
        target=txsscan_models_target,
    shell:
        r"""
        (
            set -e
            rm -rf {params.target:q}
            mkdir -p {params.target:q}
            macsydata install --org {params.org:q} --target {params.target:q} \
                "{params.package}=={params.version}"
            got=$(sed -n 's/^vers: *//p' {output.models:q}/metadata.yml)
            if [ "$got" != "{params.version}" ]; then
                echo "installed {params.package} $got, expected {params.version}" >&2
                exit 1
            fi
            echo "{params.package} $got installed"
        ) >{log} 2>&1
        """


rule relative_gembase:
    """A relative's proteome in gembase form, and its id map, from its annotation."""
    input:
        gff="resources/secretion_systems/genomes/{genome}/{genome}.gff",
        faa="resources/secretion_systems/genomes/{genome}/{genome}.genes.faa",
        code=deps("gembase_from_gff.py"),
    output:
        faa="results/txsscan/gembase/{genome}.faa",
        map="results/txsscan/gembase/{genome}.map.tsv",
    log:
        "logs/txsscan/gembase/{genome}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gembase_from_gff.py"


rule txsscan_relative:
    """One relative genome, whole proteome, gembase ids."""
    input:
        faa=rules.relative_gembase.output.faa,
        models=rules.txsscan_models.output.models,
    output:
        best="results/txsscan/relatives/{genome}/best_solution.tsv",
    log:
        "logs/txsscan/relatives/{genome}.log",
    conda:
        "../envs/macsyfinder.yaml"
    threads: 12
    params:
        outdir=txsscan_outdir,
        models_dir=txsscan_models_dir,
    shell:
        r"""
        # macsyfinder refuses to write into a non-empty directory, so a rerun
        # has to clear the previous attempt first.
        rm -rf {params.outdir:q}
        macsyfinder \
            --sequence-db {input.faa:q} \
            --db-type gembase \
            --replicon-topology linear \
            --models-dir {params.models_dir:q} \
            --models TXSScan all \
            --out-dir {params.outdir:q} \
            --index-dir {params.outdir:q} \
            --worker {threads} \
            --mute >{log} 2>&1
        """


rule txsscan_mf6_replicon:
    """One MF6 replicon: closed and circular, so gene order is meaningful."""
    input:
        faa="results/txsscan/MF6/proteomes/{replicon}.faa",
        models=rules.txsscan_models.output.models,
    output:
        best="results/txsscan/MF6/{replicon}/best_solution.tsv",
    log:
        "logs/txsscan/MF6/{replicon}.log",
    conda:
        "../envs/macsyfinder.yaml"
    threads: 12
    params:
        outdir=txsscan_outdir,
        models_dir=txsscan_models_dir,
    shell:
        r"""
        # macsyfinder refuses to write into a non-empty directory, so a rerun
        # has to clear the previous attempt first.
        rm -rf {params.outdir:q}
        macsyfinder \
            --sequence-db {input.faa:q} \
            --db-type ordered_replicon \
            --replicon-topology circular \
            --models-dir {params.models_dir:q} \
            --models TXSScan all \
            --out-dir {params.outdir:q} \
            --index-dir {params.outdir:q} \
            --worker {threads} \
            --mute >{log} 2>&1
        """


rule txsscan_mf6_merge:
    """Merge the three per-replicon MF6 runs, as macsyfinder itself provides."""
    input:
        best=expand(
            "results/txsscan/MF6/{replicon}/best_solution.tsv",
            replicon=MF6_REPLICONS,
        ),
    output:
        best="results/txsscan/MF6/merged/merged_best_solution.tsv",
    log:
        "logs/txsscan/MF6_merge.log",
    conda:
        "../envs/macsyfinder.yaml"
    params:
        indirs=txsscan_mf6_dirs,
        outdir=txsscan_outdir,
    shell:
        r"""
        msf_merge -o {params.outdir:q} {params.indirs} >{log} 2>&1
        """


rule txsscan_all_systems:
    """One all_systems.tsv over every genome, the shape Table 01 joins against."""
    input:
        relatives=expand(
            "results/txsscan/relatives/{genome}/best_solution.tsv",
            genome=SS_ISOLATES,
        ),
        mf6=rules.txsscan_mf6_merge.output.best,
        code=deps("txsscan_collect.py"),
    output:
        systems="results/txsscan/all_systems.tsv",
    log:
        "logs/txsscan/all_systems.log",
    conda:
        "../envs/python.yaml"
    params:
        genomes=SS_ISOLATES,
    script:
        "../scripts/txsscan_collect.py"


rule mf6_txsscan_genes:
    """MF6's per-protein TXSScan hits, from its own run plus the annotation."""
    input:
        best=rules.txsscan_mf6_merge.output.best,
        gff=res("gff"),
        code=deps("mf6_txsscan_genes.py"),
    output:
        genes="results/txsscan/MF6/txsscan_genes.tsv",
    log:
        "logs/txsscan/mf6_txsscan_genes.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mf6_txsscan_genes.py"
