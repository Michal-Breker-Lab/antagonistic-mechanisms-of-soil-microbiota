rule candidate_list:
    input:
        de=DE_MF6,
        mutants=de_mutant_files(),
        matrix=rules.split_matrices.output.maxlfq,
        gff=res("gff"),
        topology=res("deeplocpro"),
        t6ss=res("t6ss"),
        pfams=res("all_pfams"),
        interpro=res("interpro"),
        eggnog=res("eggnog"),
        code=deps("candidate_list.py"),
    output:
        table="results/toxins/candidate_list.tsv",
    log:
        "logs/candidate_list.log",
    conda:
        "../envs/figures.yaml"
    params:
        root=subpath(output.table, ancestor=2),
        absent=config["absent_replicon"],
        locus_prefix=config["locus_prefix"],
        min_reps=config["candidates"]["min_reps"],
        pct_gate=config["onoff_percentile_gate"],
    script:
        "../scripts/candidate_list.py"


rule candidate_heatmap:
    input:
        candidates=rules.candidate_list.output.table,
        centred="results/DE/day2_log2_maxlfq_core_centred.tsv",
        gff=res("gff"),
        de=DE_MF6,
        mutants=de_mutant_files(),
        code=deps("candidate_heatmap.py"),
    output:
        table=f"results/toxins/plots/{PREFIX}_abundance_heatmap.tsv",
        pdf=f"results/toxins/plots/{PREFIX}_abundance_heatmap.pdf",
        svg=f"results/toxins/plots/{PREFIX}_abundance_heatmap.svg",
    log:
        "logs/candidate_heatmap.log",
    conda:
        "../envs/figures.yaml"
    params:
        root=subpath(input.candidates, ancestor=2),
        outdir=subpath(output.table, parent=True),
        prefix=PREFIX,
        absent=config["absent_replicon"],
        gff=heatmap_gff,
    script:
        "../scripts/candidate_heatmap.py"
