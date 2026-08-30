rule supp_tables:
    input:
        matrix=rules.split_matrices.output.maxlfq,
        de=DE_MF6,
        mutants=de_mutant_files(),
        enrichment=ENRICHMENT,
        clusters=rules.core_clusters.output,
        heatmap=rules.candidate_heatmap.output.table,
        gff=res("gff"),
        biolib=res("deeplocpro"),
        candidates=rules.candidate_list.output.table,
        pred=res("t6ss"),
        txss=rules.mf6_txsscan_genes.output.genes,
        txsscan_tree=res("txsscan_tree"),
        txsscan_systems="results/txsscan/all_systems.tsv",
        txsscan_files=txsscan_files,
        code=deps("build_supp_tables.py"),
    output:
        SUPP,
        lfq_raw=published("lfq_raw"),
        lfq_normalized=published("lfq_normalized"),
        potential_toxins=published("potential_toxins"),
        secretion_systems=published("secretion_systems"),
    log:
        "logs/supp_tables.log",
    conda:
        "../envs/figures.yaml"
    params:
        root=subpath(input.matrix, ancestor=2),  # results/Proteomics/x.tsv -> results
        variant=ONOFF_VARIANT,
        alpha=config["figures"]["bubble"]["alpha"],
        lfq_decimals=config["supp_tables"]["lfq_decimals"],
        render=render_targets,
        isolates=SS_ISOLATES,
        txsscan_runs="results/txsscan/relatives",
    script:
        "../scripts/build_supp_tables.py"


rule table_de_proteomics:
    input:
        s2="results/supp_tables/DE_all_contrasts.tsv",
        onoff="results/DE/MF6/on_off_withC_d2_vs_alone_d2.tsv",
        genes=res("antismash_genes"),
        code=deps("table_de_proteomics.py"),
    output:
        xlsx=published("de_proteomics"),
        tsv="results/supp_tables/Table_04_DE_proteomics.tsv",
    log:
        "logs/table_de_proteomics.log",
    conda:
        "../envs/pandas.yaml"
    params:
        gate=config["onoff_percentile_gate"],
        locus_prefix=config["locus_prefix"],
    script:
        "../scripts/table_de_proteomics.py"


rule table_enrichment:
    input:
        enrichment=ENRICHMENT,
        s4="results/supp_tables/functional_enrichment.tsv",
        code=deps("table_enrichment.py"),
    output:
        xlsx=published("enrichment"),
        tsv="results/supp_tables/Table_05_enrichment.tsv",
    log:
        "logs/table_enrichment.log",
    conda:
        "../envs/pandas.yaml"
    params:
        enrich_dir=enrich_dir_of,
        variant=ONOFF_VARIANT,
    script:
        "../scripts/table_enrichment.py"
