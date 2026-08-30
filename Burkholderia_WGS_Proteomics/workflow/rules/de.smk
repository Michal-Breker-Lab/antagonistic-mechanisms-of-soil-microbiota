rule de_mf6:
    input:
        matrix=rules.split_matrices.output.maxlfq,
        gff=res("gff"),
        code=deps("de_mf6.R"),
    output:
        DE_MF6,
    log:
        "logs/de_mf6.log",
    conda:
        "../envs/limma.yaml"
    params:
        absent=config["absent_replicon"],
        drop=QC_EXCLUDE,
        day=config["de"]["day"],
        min_valid=config["de"]["min_valid"],
        min_on=config["de"]["min_on"],
        fdr=config["de"]["fdr"],
        lfc=config["de"]["lfc"],
        min_core=config["de"]["min_core"],
        n_runs=N_DE_RUNS,
        n_mf6=N_MF6_RUNS,
    script:
        "../scripts/de_mf6.R"


rule de_mutants:
    input:
        matrix=rules.split_matrices.output.maxlfq,
        gff=res("gff"),
        code=deps("de_mutants.R"),
    output:
        de_mutant_files(),
    log:
        "logs/de_mutants.log",
    conda:
        "../envs/limma.yaml"
    params:
        absent=config["absent_replicon"],
        drop=QC_EXCLUDE,
        day=config["de"]["day"],
        min_on=config["de"]["min_on"],
        fdr=config["de"]["fdr"],
        lfc=config["de"]["lfc"],
        min_core=config["de"]["min_core"],
        n_runs=N_DE_RUNS,
        mutants=MUTANTS,
    script:
        "../scripts/de_mutants.R"


rule enrichment_mf6:
    input:
        de=DE_MF6,
        eggnog=res("eggnog"),
        kegg=res("kegg_names"),
        kegg_classes=res("kegg_classes"),
        code=deps("enrichment_mf6.R"),
    output:
        ENRICHMENT,
    log:
        "logs/enrichment_mf6.log",
    conda:
        "../envs/limma.yaml"
    params:
        de_dir=de_dir_of,
        drop_nb=config["kegg"]["drop_non_bacterial"],
        drop_ov=config["kegg"]["drop_overview"],
        min_term=config["kegg"]["min_term_size"],
        min_hits=config["kegg"]["min_hits"],
        drop_ret=config["kegg"]["drop_retired"],
        pct_gate=config["onoff_percentile_gate"],
    script:
        "../scripts/enrichment_mf6.R"
