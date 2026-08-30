# Volcano figures.  Everything else in this workflow emits tables; these are the
# two figures the manuscript uses, ported from scripts/proteomics/06*.py.
#
# Only what the figures actually read is declared.  Both scripts used to also open
# the standalone SignalP / DeepLocPro tables, but their values were consumed only
# by the `extracellular` highlight mode, which nothing selects: the DeepLoc and
# SignalP columns in these figures come from candidate_list.tsv, which takes them
# from biolib_topology_summary.tsv.  Same for the broad toxin list, read only by
# the `toxin` mode.  Three inputs no rule could ever open, so all three are gone.
#
# The ON/OFF percentile gate is the script default (50) and is passed explicitly
# anyway, so the figure and config/config.yaml cannot drift apart.

VOLCANO_INPUTS = dict(
    candidates=rules.candidate_list.output.table,
    final=res("final_list"),
    # MF6's TXSScan hits, computed by the txsscan rules - not a staged table
    txss=rules.mf6_txsscan_genes.output.genes,
)

MF6_VOLCANO_TAG = "candidates_plus_txsscan"
MUT_VOLCANO_TAG = "candidates_plus_txsscan_final"


rule volcano_mf6:
    """MF6 co-culture volcano: curated candidates coloured by evidence, TXSScan
    components by system, with the ON/OFF band above it."""
    input:
        **VOLCANO_INPUTS,
        de=DE_MF6,
        code=deps("volcano_mf6.py"),
    output:
        pdf=f"results/DE/MF6/plots/volcano_{MF6_VOLCANO_TAG}_coculture_d2.pdf",
        svg=f"results/DE/MF6/plots/volcano_{MF6_VOLCANO_TAG}_coculture_d2.svg",
        table=f"results/DE/MF6/plots/volcano_{MF6_VOLCANO_TAG}_coculture_d2_highlighted.tsv",
        onoff=f"results/DE/MF6/plots/volcano_{MF6_VOLCANO_TAG}_coculture_d2_onoff_highlighted.tsv",
    log:
        "logs/volcano_mf6.log",
    conda:
        "../envs/figures.yaml"
    params:
        de_dir=de_dir_of,
        outdir=subpath(output.pdf, parent=True),
        pct=config["onoff_percentile_gate"],
        fdr=config["figures"]["volcano"]["fdr"],
        lfc=config["figures"]["volcano"]["lfc"],
        highlight="candidates_txss",
    script:
        "../scripts/volcano_mf6.py"


rule volcano_mutant:
    """The three contrasts of one mutant, curated-shortlist colouring."""
    input:
        **VOLCANO_INPUTS,
        de=de_mutant_files(),
        code=deps("volcano_mutants.py"),
    output:
        pdf=expand(
            f"results/DE/{{{{strain}}}}/plots/volcano_{MUT_VOLCANO_TAG}_{{{{strain}}}}_{{c}}.pdf",
            c=["withC_vs_alone", "withC_vs_MF6_withC", "alone_vs_MF6_alone"],
        ),
        svg=expand(
            f"results/DE/{{{{strain}}}}/plots/volcano_{MUT_VOLCANO_TAG}_{{{{strain}}}}_{{c}}.svg",
            c=["withC_vs_alone", "withC_vs_MF6_withC", "alone_vs_MF6_alone"],
        ),
        table=expand(
            f"results/DE/{{{{strain}}}}/plots/volcano_{MUT_VOLCANO_TAG}_{{{{strain}}}}_{{c}}_highlighted.tsv",
            c=["withC_vs_alone", "withC_vs_MF6_withC", "alone_vs_MF6_alone"],
        ),
    log:
        "logs/volcano_{strain}.log",
    conda:
        "../envs/figures.yaml"
    params:
        de_root=mutant_de_root,
        outdir=pdf_dir,
        # volcano_mutants.py loads volcano_mf6.py by path, so this rule must
        # supply the thresholds that module reads at import time too.
        fdr=config["figures"]["volcano"]["fdr"],
        lfc=config["figures"]["volcano"]["lfc"],
        highlight="candidates_txss_final",
    script:
        "../scripts/volcano_mutants.py"


rule enrichment_bubbles:
    """COG + KEGG bubble plot of the MF6 co-culture enrichment.

    Drawn from the de_plus_onoff_p50 variant, not the script's own default
    (de_plus_onoff): every other output in this workflow gates ON/OFF calls at the
    50th percentile, and Table S9 already takes its DE block from the p50 tables.
    A figure on one foreground beside a table on another is the kind of mismatch
    nobody catches in review.

    The variant is read from config, so the figure cannot drift from the gate the
    rest of the pipeline uses.
    """
    input:
        enrichment=ENRICHMENT,
        code=deps("enrichment_bubbles.py"),
    output:
        pdf="results/DE/MF6/enrichment/plots/bubble_coculture_d2_onoff_p50.pdf",
        svg="results/DE/MF6/enrichment/plots/bubble_coculture_d2_onoff_p50.svg",
    log:
        "logs/enrichment_bubbles.log",
    conda:
        "../envs/figures.yaml"
    params:
        indir=enrich_dir_of,
        outdir=subpath(output.pdf, parent=True),
        variant=ONOFF_VARIANT,
        alpha=config["figures"]["bubble"]["alpha"],
        fdr_cap=config["figures"]["bubble"]["fdr_cap"],
    script:
        "../scripts/enrichment_bubbles.py"


rule cluster_heatmap:
    """Day-2 core clustermap: KEGG enrichment | z-scored heatmap | COG enrichment.

    Consumes the three tables core_clusters writes, so it picks up the KEGG class
    filtering automatically - the panels are drawn from core_pathway_enrichment.tsv,
    which no longer carries the Human Diseases / Organismal Systems / overview maps.

    min_count is deliberately NOT the script's default of 2.  The enrichment now
    gates at k >= config[kegg][min_hits], and a figure that draws terms the table
    would not test is the mismatch this whole change set exists to remove.
    """
    input:
        clusters=rules.core_clusters.output.clusters,
        cog=rules.core_clusters.output.cog,
        pathway=rules.core_clusters.output.pathway,
        code=deps("cluster_heatmap.py"),
    output:
        pdf="results/Clustermap/core_heatmap_enrichment.pdf",
        svg="results/Clustermap/core_heatmap_enrichment.svg",
    log:
        "logs/cluster_heatmap.log",
    conda:
        "../envs/figures.yaml"
    params:
        indir=subpath(input.clusters, parent=True),
        prefix="core",
        fdr=config["clusters"]["fdr"],
        min_count=config["kegg"]["min_hits"],
    script:
        "../scripts/cluster_heatmap.py"


rule pca:
    input:
        matrix=rules.split_matrices.output.maxlfq,
        gff=res("gff"),
        code=deps("pca.py"),
    output:
        pdf=f"results/PCA/{{variant}}.pdf",
        svg=f"results/PCA/{{variant}}.svg",
        scores=f"results/PCA/{{variant}}_scores.tsv",
        variance=f"results/PCA/{{variant}}_variance.tsv",
    log:
        "logs/{variant}.log",
    conda:
        "../envs/figures.yaml"
    params:
        root=pca_root,
        outdir=pca_outdir,
        data=pca_data,
        gff=pca_gff,
        days=pca_days,
        cond=pca_cond,
        # Both QC-failed runs, straight from the sheet.  This used to add
        # 34F7_alone_d3_1 by hand because config carried only the other one;
        # qc_pass marks both, so nothing is appended here any more.
        exclude=QC_EXCLUDE,
        qc_failed=QC_EXCLUDE,
        absent=config["absent_replicon"],
    script:
        "../scripts/pca.py"
