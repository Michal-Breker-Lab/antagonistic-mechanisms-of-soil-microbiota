import os

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

# Outputs are written under results/; inputs come from res().


def res(name):
    """Path of a declared resource, by name.  See `resources:` in config.yaml."""
    try:
        return config["resources"][name]
    except KeyError:
        raise WorkflowError(
            f"unknown resource {name!r}; config/config.yaml declares: "
            + ", ".join(sorted(config["resources"]))
        ) from None


def up(path, n=1):
    """`n` directories up from a file in the results tree."""
    for _ in range(n):
        path = os.path.dirname(path)
    return path


SRC = "workflow/scripts"

_IMPORTS = {
    "split_matrices.py": ["fp_common.py"],
    "enrichment_bubbles.py": ["svg_lib.py"],
    "volcano_mf6.py": ["svg_lib.py"],
    "volcano_mutants.py": ["svg_lib.py", "volcano_mf6.py"],
    "table_de_proteomics.py": ["supp_xlsx_style.py"],
    "table_enrichment.py": ["supp_xlsx_style.py"],
    "build_supp_tables.py": ["supp_xlsx_style.py"],
    "table_tn_insertion_sites.py": ["supp_xlsx_style.py", "tn_annotate.py"],
    "table_coverage_breadth.py": ["supp_xlsx_style.py"],
    "find_insertions.py": ["tn_annotate.py"],
    "sanger_sites.py": ["tn_annotate.py"],
}


def deps(name):
    """The local modules a script imports, as rule inputs.

    Rules run their scripts through `script:`, so snakemake already hashes the
    script itself - but not the siblings that script imports.  Declaring those
    keeps a change to svg_lib.py or supp_xlsx_style.py from leaving figures and
    workbooks stale-but-considered-current.
    """
    return [f"{SRC}/{m}" for m in _IMPORTS.get(name, [])]


# Three PCA views of the same matrix, selected by the `variant` wildcard:
#   all_samples  every run, both days - justifies the day-2-only DE design
#   day2         the runs the DE uses
#   withC_d2_d3  co-cultures only, both days
PCA_VARIANTS = {
    "pca_all_samples": {"days": ["2", "3"], "cond": ["alone", "withC"]},
    "pca_day2": {"days": ["2"], "cond": ["alone", "withC"]},
    "pca_withC_d2_d3": {"days": ["2", "3"], "cond": ["withC"]},
}


# Named input functions: rules carry no lambdas.
def pdf_dir(wildcards, output):
    return os.path.dirname(output.pdf[0])


def mutant_de_root(wildcards, output):
    """results/DE/<strain>/plots/x.pdf -> results/DE"""
    return up(output.pdf[0], 3)


def de_dir_of(wildcards, input):
    return os.path.dirname(input.de[0])


def enrich_dir_of(wildcards, input):
    return os.path.dirname(input.enrichment[0])


# Paths relative to the results root: the scripts join them against --root.
def pca_root(wildcards, input):
    return up(input.matrix, 2)


def pca_data(wildcards, input):
    return os.path.relpath(input.matrix, up(input.matrix, 2))


def pca_gff(wildcards, input):
    return os.path.relpath(input.gff, up(input.matrix, 2))


def pca_outdir(wildcards, output):
    return os.path.relpath(os.path.dirname(output.pdf), up(output.pdf, 2))


def pca_days(wildcards):
    return PCA_VARIANTS[wildcards.variant]["days"]


def pca_cond(wildcards):
    return PCA_VARIANTS[wildcards.variant]["cond"]


def heatmap_gff(wildcards, input):
    return os.path.relpath(input.gff, up(input.candidates, 2))


ONOFF_VARIANT = f"de_plus_onoff_p{config['onoff_percentile_gate']}"

PREFIX = config["heatmap"]["prefix"]

samples = pd.read_csv(config["samples"], sep="\t", dtype=str)
validate(samples, "../schemas/samples.schema.yaml")

# Runs dropped from every analysis, from the sample sheet's qc_pass column.
QC_EXCLUDE = sorted(samples.loc[samples["qc_pass"] == "no", "run"])


# Mutants and the run counts the DE scripts guard on, derived from the sheet.
def _day_runs(day):
    ok = samples[samples["qc_pass"] == "yes"]
    return ok[ok["day"] == f"d{day}"]


DE_DAY = config["de"]["day"]
_d = _day_runs(DE_DAY)
MUTANTS = sorted(set(_d["strain"]) - {"MF6"})
N_DE_RUNS = len(_d)  # all strains, the shared core
N_MF6_RUNS = int((_d["strain"] == "MF6").sum())  # the MF6-only fit

DE_MF6 = expand(
    f"results/DE/MF6/{{f}}",
    f=[
        "DE_coculture_d2.tsv",
        "DE_summary.tsv",
        "on_off_withC_d2_vs_alone_d2.tsv",
        "MF6_log2_maxlfq_core_centred.tsv",
    ],
)

MUT_CONTRASTS = [
    "{m}_withC_vs_alone",
    "{m}_alone_vs_MF6_alone",
    "{m}_withC_vs_MF6_withC",
]


def de_mutant_files():
    out = expand(
        f"results/DE/{{f}}",
        f=[
            "day2_core_proteins.tsv",
            "day2_median_offsets.tsv",
            "day2_log2_maxlfq_core_centred.tsv",
            "DE_summary_mutants_day2.tsv",
        ],
    )
    for m in MUTANTS:
        out.append(f"results/DE/{m}/DE_summary.tsv")
        out.append(f"results/DE/{m}/{m}_log2_maxlfq_core_centred.tsv")
        for c in MUT_CONTRASTS:
            out.append(f"results/DE/{m}/DE_{c.format(m=m)}.tsv")
            out.append(f"results/DE/{m}/on_off_{c.format(m=m)}.tsv")
    return out


ENRICHMENT = expand(
    f"results/DE/MF6/enrichment/{{f}}",
    f=[
        "COG_enrichment.tsv",
        "KEGG_enrichment.tsv",
        "COG_enrichment_onoff_p50.tsv",
        "KEGG_enrichment_onoff_p50.tsv",
    ],
)

SUPP = expand(
    f"results/supp_tables/{{f}}",
    f=[
        "Table_02_LFQ_raw_proteomics.tsv",
        "Table_03_LFQ_normalized_proteomics.tsv",
        "DE_all_contrasts.tsv",
        "Table_08_potential_toxins.tsv",
        "functional_enrichment.tsv",
        "Table_01_secretion_systems.tsv",
        "MF6_supplementary_tables.xlsx",
    ],
)

FINAL = expand(
    f"results/supp_tables/{{f}}",
    f=[
        "Table_04_DE_proteomics.tsv",
        "Table_05_enrichment.tsv",
    ],
)


MF6_VOLCANO = [
    f"results/DE/MF6/plots/volcano_candidates_plus_txsscan_coculture_d2.{ext}"
    for ext in ("pdf", "svg")
] + [
    f"results/DE/MF6/plots/volcano_candidates_plus_txsscan_coculture_d2_{t}.tsv"
    for t in ("highlighted", "onoff_highlighted")
]

MUT_VOLCANO = [
    f"results/DE/{m}/plots/volcano_candidates_plus_txsscan_final_{m}_{c}.{ext}"
    for m in MUTANTS
    for c in ("withC_vs_alone", "withC_vs_MF6_withC", "alone_vs_MF6_alone")
    for ext in ("pdf", "svg")
] + [
    f"results/DE/{m}/plots/volcano_candidates_plus_txsscan_final_{m}_{c}_highlighted.tsv"
    for m in MUTANTS
    for c in ("withC_vs_alone", "withC_vs_MF6_withC", "alone_vs_MF6_alone")
]

BUBBLE = [
    f"results/DE/MF6/enrichment/plots/bubble_coculture_d2_onoff_p50.{ext}"
    for ext in ("pdf", "svg")
]

CLUSTERMAP = [
    f"results/Clustermap/core_heatmap_enrichment.{ext}" for ext in ("pdf", "svg")
]

PCA = [
    f"results/PCA/{v}{sfx}"
    for v in ("pca_all_samples", "pca_day2", "pca_withC_d2_d3")
    for sfx in (".pdf", ".svg", "_scores.tsv", "_variance.tsv")
]

# The abundance heatmap figure behind the supplementary toxin table.
HEATMAP = [
    f"results/toxins/plots/{PREFIX}_abundance_heatmap.{ext}" for ext in ("pdf", "svg")
]

FIGURES = MF6_VOLCANO + MUT_VOLCANO + BUBBLE + CLUSTERMAP + PCA + HEATMAP


# --------------------------------------------------------------- Tn5 insertions
# Reads are fetched from SRA by accession; the annotation is res("gff").


def _read_single_fasta(path):
    """(name, length) of the one record in a FASTA. Errors if there is not exactly one."""
    lengths, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    if len(lengths) != 1:
        raise WorkflowError(
            f"{path}: expected exactly one record, found {len(lengths)}: "
            + ", ".join(lengths)
        )
    return next(iter(lengths.items()))


# Cassette name and length are read off the FASTA, not declared in config.
CASSETTE_NAME, CASSETTE_LEN = _read_single_fasta(res("tn5"))

TN_PLATFORMS = ["illumina", "nanopore"]

tn_clones = pd.read_csv(config["clones"], sep="\t", dtype=str, comment="#")
validate(tn_clones, "../schemas/clones.schema.yaml")
tn_clones = tn_clones.set_index("clone", drop=False)

sanger_clones = pd.read_csv(config["sanger"], sep="\t", dtype=str, comment="#")
validate(sanger_clones, "../schemas/sanger.schema.yaml")

NANOPORE = sorted(tn_clones.loc[tn_clones["platform"] == "nanopore", "clone"])
ILLUMINA = sorted(tn_clones.loc[tn_clones["platform"] == "illumina", "clone"])
TN_CLONES = NANOPORE + ILLUMINA

# Both Sanger reads of every Sanger-only clone, as rule inputs.
SANGER_READS = sorted(sanger_clones["fwd"]) + sorted(sanger_clones["rev"])

TN_BAMDIR = {"nanopore": "results/Tn/bam", "illumina": "results/Tn/bam_illumina"}


def sra_of(wildcards):
    """The SRA run accession backing a clone."""
    return tn_clones.loc[wildcards.clone, "sra"]


def reads_dir(wildcards, output):
    """Scratch directory fasterq-dump writes into, taken from the rule's output."""
    return os.path.dirname(output[0])


def nanopore_fastq(wildcards):
    return f"results/Tn/reads/{wildcards.clone}.fastq.gz"


def illumina_fastq(wildcards):
    return expand(
        "results/Tn/reads/{clone}_{mate}.fastq.gz", clone=wildcards.clone, mate=[1, 2]
    )


def tn_clones_of(platform):
    return NANOPORE if platform == "nanopore" else ILLUMINA


def platform_of(wildcards):
    return wildcards.platform


def platform_bams(wildcards):
    d = TN_BAMDIR[wildcards.platform]
    return [f"{d}/{c}.bam" for c in tn_clones_of(wildcards.platform)]


def platform_bais(wildcards):
    d = TN_BAMDIR[wildcards.platform]
    return [f"{d}/{c}.bam.bai" for c in tn_clones_of(wildcards.platform)]


def all_tn_bams(wildcards):
    return [f"{TN_BAMDIR[p]}/{c}.bam" for p in TN_PLATFORMS for c in tn_clones_of(p)]


def all_tn_bais(wildcards):
    return [
        f"{TN_BAMDIR[p]}/{c}.bam.bai" for p in TN_PLATFORMS for c in tn_clones_of(p)
    ]


def tn_bam_specs(wildcards, input):
    """`platform:bam` per clone, in the order coverage_breadth.sh expects.

    Built from the rule's own resolved input list rather than re-deriving the
    paths, so the specs cannot name a BAM the rule did not depend on.
    """
    order = [(p, c) for p in TN_PLATFORMS for c in tn_clones_of(p)]
    return [f"{p}:{bam}" for (p, _c), bam in zip(order, input.bams)]


def tn_mapq_csv(wildcards):
    return ",".join(str(q) for q in config["tn"]["coverage_mapq"])


def tn_detect(wildcards):
    """The detect block, with the cassette identity folded in.

    find_insertions.py reads cassette_name/cassette_len out of the same dict as
    the thresholds; both are read off the cassette FASTA, so the values the
    caller uses cannot disagree with the sequence it aligned against.
    """
    return dict(
        config["tn"]["detect"],
        cassette_name=CASSETTE_NAME,
        cassette_len=CASSETTE_LEN,
    )


TN_TABLES = [
    "results/supp_tables/Table_06_Tn_insertion_sites.tsv",
    "results/supp_tables/Table_07_coverage_breadth.tsv",
]


# Table 01 declares the individual staged files it reads, not the directory,
# so a change to any of them re-runs the table.
SS_ISOLATES = config["secretion_systems"]["isolates"]


def txsscan_files(wildcards):
    """What Table 01 reads.

    best_solution.tsv and all_systems.tsv are COMPUTED by the txsscan rules -
    macsyfinder 2.1.6 over the staged TXSScan 1.1.4 model set - and only the id
    table is staged: the isolate -> IMG Genome ID mapping, which is an external
    accession and cannot come from a genome.  Everything else - proteomes, id
    maps, hit tables - is computed.
    """
    root = res("txsscan_tree")
    out = [
        f"{root}/IMG_info.csv",
        "results/txsscan/all_systems.tsv",
    ]
    for g in SS_ISOLATES:
        out.append(f"results/txsscan/relatives/{g}/best_solution.tsv")
        out.append(f"results/txsscan/gembase/{g}.map.tsv")
    return out


# --------------------------------------------------------------------- TXSScan
# MF6 is closed, so it is searched one circular replicon at a time; the relatives
# are drafts searched whole, linear.

MF6_REPLICONS = ["chr1", "chr2", "chr3"]


def txsscan_outdir(wildcards, output):
    """macsyfinder writes a directory of files; it is named by the rule's output."""
    return os.path.dirname(output[0])


def txsscan_models_dir(wildcards, input):
    """--models-dir is the PARENT of the model package, not the package itself."""
    return os.path.dirname(input.models)


def txsscan_models_target(wildcards, output):
    """macsydata --target is the parent too: it creates <target>/<package>."""
    return os.path.dirname(output.models)


def txsscan_mf6_dirs(wildcards, input):
    """The per-replicon result directories msf_merge combines."""
    return [os.path.dirname(f) for f in input.best]


# ------------------------------------------------------- published table names
# Rules write straight to results/published/ under the name config gives them;
# published() is the only place a rule asks for that path.

PUB_GENERATED = config["published_tables"]["generated"]


def published(key):
    """Where the workbook `key` is published, per config."""
    try:
        return f"results/published/{PUB_GENERATED[key]}"
    except KeyError:
        raise WorkflowError(
            f"unknown published table {key!r}; config declares: "
            + ", ".join(sorted(PUB_GENERATED))
        ) from None


def render_targets(wildcards, output):
    """TSV basename -> published workbook path, for the three build_supp_tables
    renders.  Taken from the rule's own outputs so the script cannot write a file
    the rule did not declare."""
    return {
        "Table_02_LFQ_raw_proteomics.tsv": output.lfq_raw,
        "Table_03_LFQ_normalized_proteomics.tsv": output.lfq_normalized,
        "Table_08_potential_toxins.tsv": output.potential_toxins,
        "Table_01_secretion_systems.tsv": output.secretion_systems,
    }


PUBLISHED_TABLES = [f"results/published/{fn}" for fn in PUB_GENERATED.values()]


def get_final_output(wildcards):
    return (
        DE_MF6
        + de_mutant_files()
        + ENRICHMENT
        + SUPP
        + FINAL
        + FIGURES
        + TN_TABLES
        + PUBLISHED_TABLES
    )
