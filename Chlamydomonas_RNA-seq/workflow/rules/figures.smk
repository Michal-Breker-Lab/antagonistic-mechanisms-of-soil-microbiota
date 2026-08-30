rule plot_pca:
    input:
        rds=rules.deseq2_init.output.rds,
    output:
        pca_plot="results/figures/pca-plot.svg",
    params:
        condition=config["deseq2"]["condition"],
        ntop=config["deseq2"]["ntop"],
        transformation=config["deseq2"]["transformation"],
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/pca-plot.log",
    script:
        "../scripts/deseq2_plot_pca.R"


rule enrichment_plot_all:
    input:
        enrichment_table=rules.topgo_clusters.output.enrichment_table,
    output:
        enrichment_plot="results/figures/enrichment_plot_all.svg",
    params:
        minGeneCount=config["Enrichment"]["minGeneCount"],
        minPval=config["Enrichment"]["minPval"],
    threads: 1
    log:
        "logs/figures/enrichment_plot_all.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/enrichment_plot.py"


rule enrichment_plot_up1_down6:
    input:
        enrichment_table=rules.topgo_up1_down6.output.enrichment_table,
    output:
        enrichment_plot="results/figures/enrichment_plot_up1_down6.svg",
    params:
        minGeneCount=config["Enrichment"]["minGeneCount"],
        minPval=config["Enrichment"]["minPval"],
    threads: 1
    log:
        "logs/figures/enrichment_plot_up1_down6.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/enrichment_plot_up1_down6.py"


rule expr_enrich_plot:
    input:
        cluster_table="results/Mfuzz/{direction}_annotated.tsv",
        zscores="results/impulseDE2/zscores_{direction}.tsv",
        enrichment_table=rules.topgo_clusters.output.enrichment_table,
    output:
        expression_plot="results/figures/{direction}_expression_enrich.svg",
    params:
        minGeneCount=config["Enrichment"]["minGeneCount"],
        minPval=config["Enrichment"]["minPval"],
        topNterms=config["Enrichment"]["topNterms"],
    threads: 1
    log:
        "logs/figures/{direction}_expression_enrich.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/enrich_expr_plot.py"


rule heatmap_metal_transporters:
    input:
        deseq2_long=rules.deseq2_de.output.long,
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        all_annotations=rules.all_annotations.output.all_annotations,
        gene_list=config["heatmap_de"]["metal_transporters"],
    output:
        heatmap="results/figures/heatmap_metal_transporters.svg",
        table="results/figures/heatmap_metal_transporters.tsv",
    params:
        heatmap_arguments={
            "col_as_name": "Alias",
            "info_cols": ["Alias", "Cluster"],
            "sort_by": "Cluster",
            "annotate": "padj",
            "inch_y": 8,
            "inch_x": 6,
            "info_cols_y": -0.42,
            "spacing": 0.5,
            "vmax": 2,
            "vmin": -2,
            "left_ticks": False,
        },
    threads: 1
    log:
        "logs/figures/heatmap_metal_transporters.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/heatmap_de.py"


rule heatmap_metal_transporters_horizontal:
    input:
        deseq2_long=rules.deseq2_de.output.long,
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        gene_list=config["heatmap_de"]["metal_transporters"],
    output:
        heatmap="results/figures/heatmap_metal_transporters_horizontal.pdf",
    params:
        heatmap_arguments={
            "col_as_name": "Alias",
            "sort_by": "Cluster",
            "annotate": "padj",
            "vmax": 2,
            "vmin": -2,
        },
    threads: 1
    log:
        "logs/figures/heatmap_metal_transporters_horizontal.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/heatmap_de_horizontal.py"


rule heatmap_manual:
    input:
        deseq2_long=rules.deseq2_de.output.long,
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        all_annotations=rules.all_annotations.output.all_annotations,
        gene_list=config["heatmap_de"]["manual_list"],
    output:
        heatmap="results/figures/heatmap_manual_list.svg",
        table="results/figures/heatmap_manual_list.tsv",
    params:
        heatmap_arguments={
            "inch_y": 20,
            "spacing": 0.6,
            "inch_x": 6,
            "annotate": "padj",
            "col_as_name": "Alias",
            "info_cols": ["Cluster"],
            "sort_by": "Cluster",
            "info_cols_y": -0.54,
            "vmax": 2,
            "vmin": -2,
        },
    threads: 1
    log:
        "logs/figures/heatmap_manual_list.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/heatmap_de.py"


rule heatmap_autophagy:
    input:
        deseq2_long=rules.deseq2_de.output.long,
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        all_annotations=rules.all_annotations.output.all_annotations,
        gene_list=config["heatmap_de"]["autophagy_list"],
    output:
        heatmap="results/figures/heatmap_autophagy.svg",
        table="results/figures/heatmap_autophagy.tsv",
    params:
        heatmap_arguments={
            "inch_y": 10,
            "spacing": 0.6,
            "inch_x": 6,
            "annotate": "padj",
            "col_as_name": "Alias",
            "info_cols": ["Cluster"],
            "sort_by": "Cluster",
            "info_cols_y": -0.5,
            "vmax": 2,
            "vmin": -2,
        },
    threads: 1
    log:
        "logs/figures/heatmap_autophagy.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/heatmap_de.py"


rule heatmap_PCD:
    input:
        deseq2_long=rules.deseq2_de.output.long,
        up_clusters="results/Mfuzz/up_annotated.tsv",
        down_clusters="results/Mfuzz/down_annotated.tsv",
        all_annotations=rules.all_annotations.output.all_annotations,
        gene_list=config["heatmap_de"]["PCD_list"],
    output:
        heatmap="results/figures/heatmap_PCD.svg",
        table="results/figures/heatmap_PCD.tsv",
    params:
        heatmap_arguments={
            "inch_y": 10,
            "spacing": 0.6,
            "inch_x": 6,
            "annotate": "padj",
            "col_as_name": "Alias",
            "info_cols": ["Cluster"],
            "info_cols_y": -0.4,
            "sort_by": "Cluster",
            "vmax": 2,
            "vmin": -2,
        },
    threads: 1
    log:
        "logs/figures/heatmap_PCD.log",
    conda:
        "../envs/python_figures.yaml"
    script:
        "../scripts/heatmap_de.py"


# rule enrichment_plot_by_direction:
#     input:
#         enrichment_table=rules.topgo_direction.output.enrichment_table,
#     output:
#         enrichment_plot="results/figures/enrichment_plot_by_direction.svg"
#     params:
#         minGeneCount=config["Enrichment"]["minGeneCount"],
#         minPval=config["Enrichment"]["minPval"],
#     threads: 1
#     log:
#         "logs/figures/enrichment_plot_by_direction.log"
#     conda:
#         "../envs/python_figures.yaml"
#     script:
#         "../scripts/enrichment_plot.py"
# rule expression_plot:
#     input:
#         cluster_table="results/Mfuzz/{direction}_annotated.tsv",
#         zscores="results/impulseDE2/zscores_{direction}.tsv",
#     output:
#         expression_plot="results/figures/{direction}_expression.svg"
#     threads: 1
#     log:
#         "logs/figures/{direction}_expression.log"
#     conda:
#         "../envs/python_figures.yaml"
#     script:
#         "../scripts/expression_plot.py"
# rule enrichment_plot:
#     input:
#         enrichment_table=rules.topgo_clusters.output.enrichment_table,
#     output:
#         enrichment_plot="results/figures/enrichment_plot.svg"
#     params:
#         minGeneCount=config["Enrichment"]["minGeneCount"],
#         minPval=config["Enrichment"]["minPval"],
#         topNterms=config["Enrichment"]["topNterms"]
#     threads: 1
#     log:
#         "logs/figures/enrichment_plot.log"
#     conda:
#         "../envs/python_figures.yaml"
#     script:
#         "../scripts/enrichment_plot.py"
