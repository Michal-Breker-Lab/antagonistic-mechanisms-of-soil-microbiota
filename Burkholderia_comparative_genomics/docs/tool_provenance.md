# Tool provenance

One section per external tool used in this analysis: repository, publication, the version
actually installed, the flags this project relies on, and the caveats found while using it.
Documentation and the primary paper were read before each tool was used, and cross-checked
against the on-disk `--version`. Unless stated otherwise, versions are those in the project
environments on the workstation (`<workdir>/envs/`); where the rebuild ran on the cluster
instead, and the version differs, `tool_version_gaps.md` records the gap.

First compiled 2026-08-03; extended as stages were added.

---

## NCBI Datasets CLI — genome + metadata acquisition

- Docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- API: `https://api.ncbi.nlm.nih.gov/datasets/v2alpha`
- Paper: O'Leary NA *et al.* (2024) Exploring and retrieving sequence and metadata for species
  across the tree of life with NCBI Datasets. *Sci Data* 11:732. DOI 10.1038/s41597-024-03571-y
- On-disk version: **18.34.0** (`envs/burk/bin/datasets`)

**Flags relied on**

| Flag | Purpose |
|---|---|
| `--inputfile` | batch accession list (we batch at 100) |
| `--include genome` | FASTA only; PGAP annotation deliberately *not* used (see Panaroo caveat) |
| `--no-progressbar` | required for clean logs under `nohup` |
| `filters.assembly_level=complete_genome` | REST filter for the census set |
| `returned_content=COMPLETE` | **required** — without it the BioSample block is omitted |

**Caveats discovered (these cost real errors):**

1. **GCA/GCF double counting.** A taxon query returns the GenBank *and* RefSeq record for the
   same physical genome as two rows. Raw counts across the six genera came to **1,523**; after
   dropping each GCA whose `assembly_info.paired_assembly.accession` is a retained GCF, the
   distinct-genome count is **771** (752 GCF + 19 GenBank-only). Any count not deduplicated this
   way is inflated roughly 2×.
2. **Replicon data is not in `dataset_report`.** Per-sequence molecule type and length come from
   the separate `/genome/accession/{acc}/sequence_reports` endpoint, one call per assembly.
3. **Rate limit** 3 req/s without an API key, 10 req/s with. We use a 0.40 s delay and
   exponential-backoff retry.
4. **Transient HTTP/2 failures.** `Download error: stream error: ...; INTERNAL_ERROR` occurs
   sporadically and leaves a truncated zip. The download script is therefore batched with a
   per-batch `.done` marker so re-running retries only failed batches.

Useful fields present in the report that we exploit: `assembly_stats.total_number_of_chromosomes`
(NCBI's own chromosome count — used as a *cross-check* against our size-based census, not as the
metric), `checkm_info.{completeness,contamination}`, `assembly_info.sequencing_tech` (needed to
down-weight c3 absence in short-read-only assemblies), and
`average_nucleotide_identity.best_ani_match` (precomputed type-strain ANI).

---

## skani — ANI and replicon-level identity

- Repo: https://github.com/bluenote-1577/skani
- Paper: Shaw J & Yu YW (2023) Fast and robust metagenomic sequence comparison through sparse
  chaining with skani. *Nature Methods* 20:1661–1665. DOI 10.1038/s41592-023-02018-3
- On-disk version: **0.3.2** (`envs/burk/bin/skani`)

**Flags relied on:** `dist` (pairwise), `triangle` (all-vs-all, `-E` for edge list),
`--qi`/`--ri` (**treat each FASTA record as its own genome** — this is how we compare individual
replicons rather than whole assemblies), `-t` threads.

**Caveat that changed the analysis design.** skani is documented as reliable only above
**~82% ANI** and where aligned fraction is **≥15%**; below that it declines to report. Cross-genus
c3 comparisons (*Burkholderia* vs *Paraburkholderia*) will fall well under this floor. Therefore
**ANI cannot be the primary signal for Stage 4 c3-coherence clustering across the genus** —
gene-content orthology (MMseqs2) is primary, and skani ANI is used only for within-species
resolution and for the Stage 5 dereplication at 99% ANI, both of which sit comfortably in range.

---

## Bakta — uniform annotation

- Repo: https://github.com/oschwengers/bakta
- Paper: Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A (2021) Bakta: rapid
  and standardized annotation of bacterial genomes via alignment-free sequence identification.
  *Microbial Genomics* 7(11):000685. DOI 10.1099/mgen.0.000685
- On-disk version: **1.12.0** (`envs/bakta/bin/bakta`)
- Database: full DB, **84 GB**, `<db-root>/bakta_db/db`

**Flags relied on:** `--db`, `--threads`, `--genus`/`--species` (from NCBI metadata),
`--prefix`, `--output`, `--complete` (all our inputs are complete replicons),
**`--keep-contig-headers`** — essential, because Stage 7 subsets the GFF3 to the c3 contigs by
their original NCBI accession; Bakta renames contigs by default and would break that join.

**Runtime:** documented at 10 ± 5 min per genome. Budget accordingly and cap concurrency to
respect Shannon's RAM headroom.

---

## Panaroo — pangenome

- Repo: https://github.com/gtonkinhill/panaroo · Docs: https://gthlab.au/panaroo/
- Paper: Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S,
  Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J (2020) Producing polished
  prokaryotic pangenomes with the Panaroo pipeline. *Genome Biology* 21(1):180.
  DOI 10.1186/s13059-020-02090-4 · PMID 32698896 · PMCID PMC7376924
- On-disk version: **1.8.0** (`envs/panaroo/bin/panaroo`)
- Verified defaults from `--help`: `--family_threshold 0.7`, `--len_dif_percent 0.98`,
  `--threads 1` (must be raised), `--codon-table 11`, aligner choices
  `muscle|muscle-super5|famsa|prank|clustal|mafft|none`, `-a {core,pan}`.

**Why Panaroo over Roary — and why we re-annotate with Bakta (decision D5).** The paper reports a
*M. tuberculosis* outbreak set where Panaroo called **~130** accessory genes while other tools
called **2,584–3,670** — close to a tenfold inflation, of which 59% traced to fragmented genes and
10% to identical sequences carrying inconsistent annotations. This is the direct evidence that
mixing NCBI PGAP annotations of different vintages would corrupt the pangenome, and it justifies
uniform re-annotation.

**Clean-mode guidance — this contradicts the originally approved D9.** The authors recommend
`strict` "when investigating genomes where rare plasmids are not expected or when phylogenetic
parameters such as gene gain and loss rates are of interest," and recommend `sensitive`
(no cluster removal) "if a researcher is interested in rare plasmids which may be hard to
distinguish from contamination." **Our Stage 7 target is a megaplasmid**, so strict is the mode
its authors warn against for this exact use.

The exact thresholds, from `panaroo --help` on the installed 1.8.0, make the risk quantitative
rather than rhetorical:

| Mode | Evidence needed to keep a suspected-contaminant gene | Refound genes |
|---|---|---|
| `strict` | present in **≥5%** of genomes | **removed** if refound more often than originally called |
| `moderate` | present in **≥1%** of genomes | kept |
| `sensitive` | deletes nothing; merge + refind only | kept |

With a dereplicated set of a few hundred genomes, a genuine c3 accessory family carried by ~10
genomes sits near **3%** — **below strict's 5% floor and above moderate's 1% floor**. Strict would
therefore delete exactly the class of variable, lineage-restricted c3 genes this project exists to
characterise. Revised plan:

| Run | Mode | Reason |
|---|---|---|
| Stage 8, chromosome-1 core (species tree) | `strict` | no rare plasmids expected on chr1; phylogenetic inference is the goal — the authors' stated case for strict |
| Stage 7, c3 pangenome | `moderate` (default) | plasmid content must not be culled as contamination |
| Stage 7 sensitivity check | `sensitive` | report how core/accessory counts shift; cheap honesty check |

**Input requirement:** GFF3 in Prokka/Bakta style. Bakta emits a combined GFF3 with the FASTA
appended after `##FASTA`, which is what Panaroo expects; when subsetting to c3 contigs the FASTA
block must be subset in step with the feature lines or Panaroo will fail on missing sequence.

---

## IQ-TREE — phylogenetic inference

- Repo: https://github.com/iqtree/iqtree3
- Paper: Wong TKF, Ly-Trong N, Ren H, Demotte P, Banos H, Roger AJ, Susko E, Bielow C, De Maio N,
  Goldman N, Hahn MW, dos Reis M, Vinh LS, Huttley G, Lanfear R, Minh BQ (2026) IQ-TREE 3:
  Phylogenomic Inference Software using Complex Evolutionary Models. *Mol Biol Evol* msag117.
  DOI 10.1093/molbev/msag117
- On-disk version: **3.1.3** (`envs/burk/bin/iqtree3`; `iqtree` also present)
- **Version actually used for the rebuilt trees: 3.1.1** (built 2026-04-08). `m10_iqtree.sbatch`
  and `m16_iqtree_finalize.sbatch` invoke Moriah's `envs/anvio-9/bin/iqtree3`, not `envs/burk`.
  Both `tree_chr1.iqtree` and `tree_pc3.iqtree` record 3.1.1 in their headers, so Methods must
  cite 3.1.1 — the 3.1.3 above is a different binary on a different host.
- Supporting methods to cite: ModelFinder (Kalyaanamoorthy *et al.* 2017, *Nat Methods*
  14:587–589, DOI 10.1038/nmeth.4285) and UFBoot2 (Hoang *et al.* 2018, *Mol Biol Evol*
  35:518–522, DOI 10.1093/molbev/msx281).

**Note:** Shannon also carries a stale `iqtree 1.6.12` at `<archive-mount>/` earlier on
`PATH`. Always invoke the env binary by absolute path.

---

## MAFFT — alignment

- Paper: Katoh K & Standley DM (2013) MAFFT multiple sequence alignment software version 7.
  *Mol Biol Evol* 30:772–780. DOI 10.1093/molbev/mst010
- On-disk version: **v7.526 (2024/Apr/26)** (`envs/burk/bin/mafft`)
- Invoked via Panaroo's `--aligner mafft`.

## MMseqs2 — ortholog clustering for replicon gene content

- Paper: Steinegger M & Söding J (2017) MMseqs2 enables sensitive protein sequence searching for
  the analysis of massive data sets. *Nat Biotechnol* 35:1026–1028. DOI 10.1038/nbt.3988
- On-disk version: **18.8cc5c** (`envs/burk/bin/mmseqs`)
- Used in Stage 4 to build ortholog groups across all secondary replicons, from which the
  gene-content Jaccard distance is computed (the primary c3-coherence signal, since ANI is out of
  range cross-genus).

## HMMER — replicon typing

- Paper: Eddy SR (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.
  DOI 10.1371/journal.pcbi.1002195
- On-disk version: **3.4 (Aug 2023)** (`envs/burk/bin/hmmsearch`)
- Used to identify chromosome 1 by *dnaA* + the ribosomal-protein superoperon (decision D2),
  rather than assuming the largest replicon is chromosome 1.

## seqkit — FASTA manipulation

- Paper: Shen W, Le S, Li Y, Hu F (2016) SeqKit: a cross-platform and ultrafast toolkit for
  FASTA/Q file manipulation. *PLoS ONE* 11:e0163962. DOI 10.1371/journal.pone.0163962
- On-disk version: **v2.13.0** (`envs/burk/bin/seqkit`)

---

## PPanGGOLiN — second pangenome method (clustering-identity sensitivity)

- Repo: https://github.com/labgem/PPanGGOLiN · Docs: https://ppanggolin.readthedocs.io/
- Paper: Gautreau G, Bazin A, Gachet M, Planel R, Burlot L, Dubois M, Perrin A, Médigue C,
  Calteau A, Cruveiller S, Matias C, Ambroise C, Rocha EPC, Vallenet D (2020) PPanGGOLiN:
  Depicting microbial diversity via a partitioned pangenome graph. *PLoS Comput Biol*
  16(3):e1007732. DOI 10.1371/journal.pcbi.1007732 · Correction: DOI 10.1371/journal.pcbi.1009687
- On-disk version: **2.2.6** — `<env>/bin/ppanggolin`
  (pre-existing env from the within-genus-HGT project; **reuse it**, see caveat 3)
- Verified defaults from `ppanggolin cluster --help` on the installed 2.2.6:
  `--identity 0.8`, `--coverage 0.8`, `--mode 1` (MMseqs2 single-linkage / connected component;
  0 = set-cover, 2 = CD-HIT-like, 3 = CD-HIT-like low-mem). Defrag strategy is ON unless
  `--no_defrag`.
- Flags we set: `--identity {0.60,0.80} --coverage 0.8 --mode 1`, then
  `msa --partition softcore --soft_core 0.95 --phylo --source dna`.

**Why PPanGGOLiN and not Panaroo for the 60%-identity question.** On the installed Panaroo 1.8.0,
`-c/--threshold` (sequence identity) defaults to **0.98** and `-f/--family_threshold` to **0.70**.
Setting `-c 0.60` would place the *initial* clustering below the family-merge threshold that runs
after it — internally incoherent. Panaroo is designed to recover diverged and fragmented genes
through the graph neighbourhood, not by relaxing identity. PPanGGOLiN exposes `--identity` as a
documented first-class MMseqs2 parameter, and additionally partitions persistent/shell/cloud with
a **statistical model rather than a fixed prevalence cutoff**, which sidesteps the arbitrary-95%
question entirely.

**Partition vocabulary (differs from Panaroo/Roary — do not conflate):** `core` = present in
**100%** of genomes; `softcore` = the `--soft_core` threshold, default **0.95**, i.e. the
equivalent of what this project calls "core"; `persistent` = the statistically inferred partition.
The 140-genome pC3 run gives exact_core 30, softcore 58, persistent 66 — quoting "core = 30"
beside a 58-family alignment looks inconsistent and is not.

**Caveats discovered (these cost real errors):**

1. **PPanGGOLiN resolves `mmseqs` from `PATH`, not from its own prefix.** Invoking the binary by
   absolute path alone dies with `FileNotFoundError: Command 'mmseqs' not found`. The env bin
   **must** be on `PATH` — exactly the same failure mode as the Bakta/tRNAscan-SE caveat above.
   The env bundles its own matching mmseqs (15.6f452); use that rather than `envs/burk`'s 18.8cc5c.
2. **`graph` is the one subcommand with no `--cpu`/`-c` flag.** Passing `-c 32` to it aborts the
   run with `unrecognized arguments`. `annotate`, `cluster`, `partition`, `write_pangenome` and
   `msa` all accept it.
3. **Do not build a fresh env for this.** A `mamba create -c conda-forge -c bioconda ppanggolin`
   on 2026-08-09 solved to **0.3.88** — an ancient build — while the pre-existing
   `within_genus_hgt` env carries 2.2.6. Always search the filesystem for an existing env before
   installing.
4. **Output filenames are not Roary's.** The Roary-compatible matrix is written as `matrix.csv`,
   not `gene_presence_absence.csv`; `write_pangenome` also emits `genomes_statistics.tsv` and
   `mean_persistent_duplication.tsv`. The latter carries a per-family `duplication_ratio` and an
   `is_single_copy_marker` flag — use it directly as the paralogue-collapse diagnostic rather than
   reimplementing one.
5. **Statistical partitioning is unreliable at small n.** For the 5- and 11-genome MF6-neighbour
   sets, `partition` was allowed to fail and the prevalence-threshold core was read from the Rtab
   instead.

**The risk this tool review was run to quantify:** at 60% identity with single-linkage clustering,
paralogues can fuse into their orthologue family, trivially inflating the core and making a core
alignment concatenate non-orthologues. Measured on the pC3 set, going 0.80 → 0.60 moved the core
single-copy fraction 0.720 → 0.603 and the persistent single-copy-marker fraction 0.983 → 0.924 —
both within the pre-declared 0.15 tolerance, so the 0.60 core alignment was accepted for
phylogenetic use. The diagnostic is in `tables/collapse_diagnostic.tsv`.

---

## Functional annotation — Bakta DbXrefs, and KEGG module completeness

### Annotation source: no new tool

The functional characterisation in §6.4 adds **no annotation tool and no database download**.
Bakta v1.12.0 / DB v6.0 (already documented above, and already run over every genome in this
project) writes **COG identifiers, COG functional-category letters, KEGG KO, EC, GO and Pfam**
into the `DbXrefs` column of its `<genome>.tsv`. All five *B. sola* genomes carry the same Bakta
version and database, so the comparison rests on one consistent annotation source. Two practical
notes: `COG:COG0592` (an orthologous group) and `COG:L` (a functional category) both appear under
the `COG:` prefix and must be told apart by pattern; and the Roary-style `matrix.csv` from
PPanGGOLiN gives the family → per-genome locus-tag mapping that joins the pangenome to these TSVs.

**Caveat that matters for interpretation:** COG coverage is *not uniform across replicons*. Only
28.7% of pC3 clade-core families receive a COG category against 67.1% on chromosome 1. Any
category-level comparison is therefore conditioned on an unevenly annotated subset, which is why
every contrast in §6.4 is computed twice — once with "no COG" as an explicit category and once
restricted to COG-assigned families — and why the coverage rate is reported as a result.

### KEGG MODULE definitions — data source

- **Source on disk**: `<db-root>/anvioKEGG/modules/`
  on Shannon — the anvi'o KEGG snapshot (dated 2025-11-08), **557 module flat files**, each with
  KEGG's own `ENTRY` / `NAME` / `DEFINITION` / `CLASS` records. Already present; nothing was
  downloaded.
- **KEGG citation**: Kanehisa M, Goto S (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes.
  *Nucleic Acids Res* 28(1):27–30. DOI 10.1093/nar/28.1.27.
- **Completeness definition**: anvi'o `anvi-estimate-metabolism`
  (https://anvio.org/help/main/programs/anvi-estimate-metabolism/), described in
  Veseli I, Chen YT, Schechter MS, et al. (2023) *eLife* 12:RP89862.
  DOI 10.7554/eLife.89862.

### The completeness algorithm as implemented

**Stepwise completeness** was used (not pathwise): the `DEFINITION` string is split on
**top-level spaces** into steps, and completeness is the fraction of steps that evaluate true.
Within a step the operators are `,` = OR, ` ` and `+` = AND, `-K#####` / `-(...)` = non-essential
and therefore stripped before evaluation, parentheses group. Steps written `--` (no KO) always
score false. Nested module references (`M#####`) are resolved iteratively, a nested module
counting as true once it clears the threshold. **Threshold 0.75**, anvi'o's default.

Implemented directly against the flat files in `s17_functional.py` on Shannon rather than by
running `anvi-estimate-metabolism`, because the latter would require building contigs databases
and re-running `anvi-run-kegg-kofams` — hours of compute to regenerate KO calls the Bakta run
already provides.

**Five caveats, all of which constrain how these numbers may be quoted:**

1. **KOs come from Bakta (UniRef-based inference), not from KOfam HMMs with per-family bit-score
   thresholds.** Completeness values are therefore **not comparable to published
   `anvi-estimate-metabolism` figures**, which use KOfam. Bakta's KO coverage is also sparser —
   157 KOs on MF6's pC3, 1,131 genome-wide.
2. **Many KEGG "modules" are single-step complexes, not pathways.** M00155 (cytochrome *c*
   oxidase) and M00153 (cytochrome *bd*) each have exactly one top-level step, so their
   completeness is binary and "complete" means "all essential subunits present" — it must not be
   reported as "a complete pathway resides on this replicon".
3. **Absence of a KO is not absence of the function.** Unannotated genes are frequent on pC3
   (see above), so low completeness is a statement about annotation as much as about biology.
4. **Replicon comparisons must use the whole genome as the complement.** "Complete on pC3 and
   absent from chromosome 1" does *not* establish uniqueness, because MF6's second replicon holds
   2,636 of its 6,895 CDS. The uniqueness results in §6.4 are computed against all 6,895.
5. **Module completeness is computed on a gene-content union**, so it says nothing about whether
   the genes are expressed, co-regulated, or functional.

---

## InterProScan — domain annotation of the pC3 unknowns

- **Version on disk**: 5.69-101.0 at `<tools>/interproscan-5.69-101.0/` on
  Shannon. Already installed; nothing was downloaded. Java 17 present (5.69 needs ≥11).
- **Paper**: Jones P, Binns D, Chang H-Y, et al. (2014) InterProScan 5: genome-scale protein
  function classification. *Bioinformatics* 30(9):1236–1240. DOI 10.1093/bioinformatics/btu031.
- **Docs**: https://interproscan-docs.readthedocs.io/ ·
  https://github.com/ebi-pf-team/interproscan-docs

**Command used** (verified against the installed `interproscan.sh --help`):

```
interproscan.sh -i pC3_clade_pangenome.faa -b pC3_ips -f TSV,GFF3 \
                -goterms -pa -iprlookup -dp -cpu 32 -T <TMPDIR under <large-storage>>
```

- `-appl` omitted, so **all 17 member databases run**: AntiFam, CDD, Coils, FunFam, Gene3D,
  Hamap, MobiDBLite, NCBIfam, PANTHER, Pfam, PIRSF, PIRSR, PRINTS, PROSITE (patterns+profiles),
  SFLD, SMART, SUPERFAMILY.
- `-dp` **disables the precalculated match lookup service**. This makes the run slower (every
  match is computed locally) but removes a dependency on an external EBI service and makes the
  result reproducible from the on-disk data alone. Chosen deliberately.
- `-iprlookup -goterms -pa` add InterPro entry, GO and pathway mappings.
- `-T` is mandatory here: Shannon's `/` is ~97% full and InterProScan writes large temporaries.

**Caveats:**

1. **InterProScan takes the FASTA identifier as the first whitespace-delimited token.** Headers
   in this project pack `family|locus|set|COG=x|product` into that token, so the product string
   is truncated at its first space in the output. Fields 1–4 are reliable; the product must be
   re-joined from `setB_functional_families.tsv`.
2. **Some member databases are structurally uninformative for function** — MobiDBLite (disorder),
   Coils, Phobius/TMHMM (topology). These are excluded when computing "did InterProScan
   characterise this unknown?", otherwise nearly every protein counts as annotated.
3. No SignalP or TMHMM licence-restricted binaries were used.

## antiSMASH — biosynthetic gene clusters on pC3

- **Version**: 7.1.0, built into `envs/antismash7` with mamba. The pre-existing base-environment
  install was **broken** — antiSMASH 7.1 calls `Bio.Seq.UnknownSeq`, removed in Biopython ≥1.80,
  and base carries 1.83. The bioconda package also does not pull its external binaries, so
  `hmmer hmmer2 fasttree prodigal glimmerhmm diamond muscle blast` were installed alongside.
- **Database**: the existing snapshot at
  `<db-root>/antismash_db` (clusterblast,
  knownclusterblast, nrps_pks, pfam, resfam, tigrfam). antiSMASH 7 was pinned specifically to
  match this database rather than installing v8 and re-downloading multiple GB.
- **Paper**: Blin K, Shaw S, Augustijn HE, et al. (2023) antiSMASH 7.0. *Nucleic Acids Research*
  51(W1):W46–W50. DOI 10.1093/nar/gkad344.

**Command used:**

```
antismash --taxon bacteria --genefinding-tool none --cpus 16 --databases <db> \
          --cb-general --cb-knownclusters --asf --pfam2go \
          --output-dir <out> <genome>_pC3.gbk
```

- **`--genefinding-tool none` is required** because the input is a Bakta GenBank that already
  carries CDS features; letting antiSMASH re-predict genes would discard the annotation the rest
  of this report is built on and break the locus-tag join.
- Input is a **GenBank subset of the pC3 contig only**, built by filtering the Bakta `.gbff` —
  antiSMASH needs nucleotide sequence with genomic context and cannot be run on a protein set.

**Caveats:**

1. `cassis` fails a prerequisite check (MEME 5.5.9 installed, 4.11.2 expected). It is a
   fungal-promoter module, off by default for bacteria, and unused here.
2. **Region counts are not comparable across replicons of different length** without
   normalising per Mb.
3. antiSMASH detects *clusters that resemble known biosynthetic logic*. A "region" is a
   prediction, not a demonstrated product, and `knownclusterblast` similarity to a
   characterised cluster is a similarity score, not an identification.

## MacSyFinder / TXSScan — secretion systems, at system level

- **Version**: MacSyFinder **2.1.6**, built into `envs/macsy` with mamba. The pre-existing
  MacSyFinder 2.1.4 at `<conda-root>/` **cannot be
  used**: TXSScan 1.1.4 models declare grammar `vers="2.1"` and 2.1.4 only accepts `2.0`,
  aborting on the first model parsed.
- **Models**: TXSScan 1.1.4, installed with `macsydata install --target <dir> TXSScan` into
  `macsy_models/`. Small HMM/XML model set, not a reference database.
- **Papers**: Néron B, Denise R, Coluzzi C, Touchon M, Rocha EPC, Abby SS (2023) MacSyFinder v2.
  *Peer Community Journal* 3:e28. DOI 10.24072/pcjournal.250 · Abby SS, Cury J, Guglielmini J,
  Néron B, Touchon M, Rocha EPC (2016) Identification of protein secretion systems in bacterial
  genomes. *Sci Rep* 6:23080. DOI 10.1038/srep23080 · Denise R, Abby SS, Rocha EPC (2019)
  *PLoS Biol* 17(7):e3000390.

**Command used:**

```
macsyfinder --db-type ordered_replicon --sequence-db <genome>_pC3.faa \
            --models TXSScan/bacteria/diderm all --models-dir macsy_models -o <out> -w 8
```

- `ordered_replicon` is correct here because each pC3 is a **single contig** and the input FASTA
  is written in genomic coordinate order, which enables the co-localisation criterion.
- `hmmsearch` must be on `PATH` — MacSyFinder aborts otherwise, the same failure mode already
  recorded for PPanGGOLiN/mmseqs and Bakta/tRNAscan-SE.

**Caveats — the important one is methodological:**

1. **Scope makes "no complete system" the expected answer for anything the chromosome hosts.**
   The screen was run on pC3 alone, so any system whose mandatory genes are chromosomal
   *cannot* reach quorum. Reporting only `best_solution.tsv` would convert that scoping choice
   into a false negative. `rejected_candidates.tsv` is therefore reported alongside, and is where
   partial/orphan modules appear.
2. A quorum-satisfying call means the model's gene content and co-localisation rules are met —
   it is not proof the system is expressed or functional.
3. The Flagellum and T3SS models share the `sct` export core, so a call of one should be checked
   against the system-specific genes (`flgB`, `flgC`, `fliE` for the flagellum) before being
   believed.

## Pre-declared screen for effectors, toxins and secretion systems

`effector_toxin_signatures.tsv` (on Shannon, mirrored in `tables/`) fixes the accession and
description patterns for every reported category **before any InterProScan output existed**.
This is the safeguard against a fishing expedition: screening ~1,400 largely uncharacterised
proteins for "toxins" will always return something.

Rules, all set in advance:

- Every declared category is reported **including those scoring zero**, so the denominator of
  the screen is visible.
- **Accession-based matches** (exact Pfam/InterPro/NCBIfam accession) are treated as confident;
  **free-text description matches** are reported separately as low-specificity, because
  substring matching produces false positives.
- **No thresholds of my own**: matches are taken at each member database's curated cutoff as
  applied by InterProScan.
- Anything notable found *outside* the declared list goes to a separate table and must be
  labelled **post hoc** wherever it is reported.
- One amendment was made before results existed: the pattern `NAD` was replaced by `NADase` and
  `NAD+ phosphorylase`, because bare `NAD` matches every NAD-binding domain in the proteome.

---

## Pyrodigal — uniform gene calls across all 771 assemblies

- **Version on disk**: 3.7.1 (present in `envs/bakta` and `envs/panaroo` on Shannon; nothing was
  installed). Wraps **Prodigal v2.6.3+31b300a** and is tested to reproduce Prodigal's output
  exactly.
- **Papers**: Larralde M (2022) Pyrodigal: Python bindings and interface to Prodigal.
  *J Open Source Softw* 7(72):4296. DOI 10.21105/joss.04296 · Hyatt D, Chen G-L, LoCascio PF,
  Land ML, Larimer FW, Hauser LJ (2010) Prodigal: prokaryotic gene recognition and translation
  initiation site identification. *BMC Bioinformatics* 11:119. DOI 10.1186/1471-2105-11-119.
- **Docs**: https://github.com/althonos/pyrodigal · https://pyrodigal.readthedocs.io/

**Why it was used.** Only **309 of the 771** genomes carry a Bakta proteome (the dereplicated set
built for the trees). Searching only those would answer a genome-content question for 40% of the
data *and* would confound "gene absent" with "genome never annotated". Pyrodigal re-calls genes on
all 771 with one caller and one parameter set, so presence/absence is not an artefact of which
pipeline happened to touch a genome.

**Settings used** (`s23_pyrodigal.py`):

```python
gf = pyrodigal.GeneFinder(meta=False)      # SINGLE mode, self-training
gf.train(*contigs, translation_table=11)
```

- **Single (self-training) mode, not `meta=True`.** Prodigal's guidance is that single mode is
  more accurate whenever ≥20 kb is available to train on; these are 7–9 Mb genomes. Meta mode is
  for short/fragmented metagenomic contigs and would lose accuracy here.
- **`train()` accepts all contigs at once** and merges them internally with `TTAATTAATTAA`
  linkers — passing the whole assembly is the documented, correct usage, not a hack.
- Defaults kept: `min_gene=90`, `min_edge_gene=60`, `max_overlap=60`, `closed=False`,
  `translation_table=11`.
- The script falls back to `meta=True` **only** if training raises, and **logs every fallback**
  rather than silently degrading.

**Caveats:**

1. Pyrodigal coordinates are its own — they do **not** correspond to Bakta locus tags. Hits are
   mapped back to a replicon by *contig* and *coordinate*, not by locus tag.
2. Gene calls on draft assemblies are truncated at contig edges; a toxin split across a contig
   break appears as two short partial proteins. This is why the search is also reported by
   coverage tier rather than presence/absence alone.
3. Re-calling genes means the protein set differs slightly from Bakta's for the 309 overlapping
   genomes. That is deliberate (uniformity), but it means counts here are not directly
   interchangeable with the Bakta-based pangenome tables.

## BLAST+ — homology search for the RHS toxin

- **Version on disk**: 2.17.0+ (`envs/bakta/bin`, build Aug 2025).
- **Paper**: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009)
  BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. DOI 10.1186/1471-2105-10-421.
- **Docs**: https://www.ncbi.nlm.nih.gov/books/NBK279684/ (command-line appendix).

**Command** (`s24_search.py`):

```
blastp -query queries.faa -db all_pyrodigal.faa -outfmt "6 qseqid sseqid pident length nident
       mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs slen qlen"
       -evalue 1e-5 -max_target_seqs 1000000 -num_threads 48 -comp_based_stats 2
```

**Caveats that shaped the settings:**

1. **`-max_target_seqs` does not return the best N hits.** It keeps the first N hits above
   threshold, in database order — so a small value silently biases results by database ordering
   (Shah N, Nute MG, Warnow T, Pop M (2019) Misunderstood parameter of NCBI BLAST impacts the
   correctness of bioinformatics workflows. *Bioinformatics* 35(9):1613–1614.
   DOI 10.1093/bioinformatics/bty833). NCBI patched the worst of it in 2.8.1 and we run 2.17.0,
   but the ceiling is set to 10^6 so the question does not arise.
2. **`pident` is computed over the aligned segment only**, not the query. For this project's
   query that is the central hazard: the first 180 aa of the RHS C-terminal domain are **45%
   Gly+Ala** with `AGGGAG` repeats, so an unrelated Gly-rich membrane protein can clear 60%
   identity over a ~50-residue window. **Identity is therefore always paired with `qcovs`**, and
   the primary call requires ≥60% identity over ≥70% of the query.
3. `-comp_based_stats 2` (blastp default) corrects E-values for composition bias and is kept
   explicitly rather than by default, because it is doing real work on this query. SEG filtering
   was **not** enabled — masking the Gly tract would also discard it as evidence; the coverage
   floor handles it instead, and both decisions are recorded rather than assumed.
4. A one-way best hit identifies **homologs, not orthologs**. Results are reported as "best hit
   per genome above threshold"; no reciprocal-best-hit or gene-tree test was run, so the word
   "ortholog" is avoided in the output tables.

### Extensions for the ten-locus search (§6.7, `s27_search_toxins.py`)

5. **Two databases are searched at once**, using BLAST's multi-database syntax:

   ```
   -db "search_rhs/all_pyrodigal.faa pyrodigal_mf6/MF6.faa"
   ```

   MF6 is not among the 771 downloaded genomes, so without the second database the source of every
   query could not be scored on the same footing as the genomes it is compared against. BLAST sums
   the effective database sizes across listed databases, so E-values remain comparable; appending
   MF6 to the 1.9 GB FASTA and rebuilding would have been equivalent but wasteful.

6. **`tblastn` is used for one locus only.** `MF6_003686` is a 120-aa ORF that Bakta does not call
   and that pyrodigal truncates to 78 aa by selecting an internal start codon — MF6 scores 65%
   coverage against its own gene, below the tier-1 floor, so a protein-level search cannot answer
   this locus at all. It is therefore searched against a nucleotide database of the assemblies:

   ```
   tblastn -query q_003686.faa -db all_genomes.fna -db_gencode 11 -evalue 1e-5
           -max_target_seqs 1000000 -num_threads 48
   ```

   This is a **different unit of search** — six-frame translation of DNA rather than predicted
   proteins — and is flagged as such in every table it touches. `qcovhsp` replaces `qcovs` for
   tiering, because tblastn HSPs are not whole genes.

7. **The query alignment span (`qstart`, `qend`) is exported for every hit.** Four of the ten
   queries begin with a signal peptide or lipobox, which aligns to the signal peptide of any
   secreted protein, and one 104-aa query is proline-rich. The coverage floor protects the long
   queries; the span is exported so a hit carried entirely by such a region is visible rather than
   silently counted. Audit result: **zero** tier-1 hits confined to residues 1–40, for any locus.

## PGAP + eggNOG annotation of MF6 — the canonical MF6 annotation

> **As of 2026-08-12, the collaborator-supplied PGAP annotation is THE MF6 annotation for this
> project.** Use it for MF6 locus tags, products, coordinates and functional calls, and label
> figures with `MF6_######` tags. The Bakta run remains the annotation for the other 771 genomes,
> where uniformity across the genus is what matters.
>
> The two describe the **same assembly**: PGAP's `contig_1` / `contig_2` / `contig_4` are
> byte-identical (matching MD5) to Bakta's `cluster_001` / `cluster_003` / `cluster_002`
> — chr1 / pC3 / chr2. PGAP omits two small contigs Bakta kept (33,003 bp and 9,111 bp), giving
> 6,841 PGAP proteins against 6,902 from the pyrodigal call over the larger contig set. No locus
> discussed in this report is affected.

Not run in this project. A collaborator in the Finkel group annotated the MF6 assembly independently with
**NCBI PGAP** and functionally annotated the result with **eggNOG-mapper 2.1.12** (DIAMOND,
`--sensmode sensitive`, `--tax_scope auto`). Outputs live under
the shared Breker–Finkel collaboration folder.

- **PGAP**: Tatusova T *et al.* (2016) NCBI prokaryotic genome annotation pipeline. *Nucleic Acids
  Res* 44(14):6614–6624. DOI 10.1093/nar/gkw569.
- **eggNOG-mapper**: Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J (2021)
  eggNOG-mapper v2. *Mol Biol Evol* 38(12):5825–5829. DOI 10.1093/molbev/msab293.

**Why it matters here:** the query file's `MF6_######` tags are PGAP tags, and they do **not**
align with this project's Bakta `CFFIHE_` tags — the offset varies by locus (`MF6_002734` =
`CFFIHE_02739`, `MF6_004947` = `CFFIHE_04987`). Cross-referencing by number produced one wrong
identification before this was noticed. **Always map between the two by sequence**;
`data/tables/toxin12_query_map.tsv` records the correspondence.

**Where the two annotations disagree**, both calls are recorded rather than silently reconciled:
Bakta calls `MF6_004284` a Peptidase C39 while PGAP and eggNOG call it an RHS / COG3209 Rhs family
protein with an `RHS_repeat` domain, and `MF6_003686` is called by PGAP but not by Bakta. The
CheckM summary shipped with the annotation reports 94.11% completeness and 4.83% contamination.

## Biological references framing the analysis

- **Agnoli K, Schwager S, Uehlinger S, Vergunst A, Viteri DF, Nguyen DT, Sokol PA, Carlier A,
  Eberl L (2012)** Exposing the third chromosome of *Burkholderia cepacia* complex strains as a
  virulence plasmid. *Mol Microbiol* 83(2):362–378. DOI 10.1111/j.1365-2958.2011.07937.x ·
  PMID 22171913.
  The premise of this project. Established c3 as **not an essential chromosomal element but a
  large plasmid**, renamed **pC3**. Cured from strains using a constructed c3 mini-replicon
  exploiting **plasmid incompatibility**; **nine c3-null strains across seven Bcc species** were
  obtained, i.e. loss is tolerated genus-wide within the Bcc. c3-null mutants showed reduced
  virulence in multiple hosts and lost D-xylose, fatty-acid and pyrimidine utilisation,
  exopolysaccharide production, and proteolytic activity in some strains.
  *Relevance:* if MF6's 1.18 Mb Replicon 3 is a pC3, its loss in four independent Tn5 clones is
  expected instability rather than an assembly or mapping artifact.
- **Agnoli K *et al.* (2014)** The third replicon of members of the *Burkholderia cepacia*
  complex, plasmid pC3, plays a role in stress tolerance. *Appl Environ Microbiol* 80(4):1340–1348.
  PMID 24334662. Follow-up extending pC3 function to stress tolerance.

---

## Pending review

R packages for Stage 11 (`caper::phylo.d`, `phytools::make.simmap`, `phylolm::phyloglm`,
`ggtree`) — to be reviewed and versioned when that stage is set up.

## GTDB-Tk — bac120 single-copy marker extraction and alignment

- Repo: https://github.com/Ecogenomics/GTDBTk · Docs: https://ecogenomics.github.io/GTDBTk/
- Paper: Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH (2022) GTDB-Tk v2: memory friendly
  classification with the genome taxonomy database. *Bioinformatics* 38:5315–5316.
  DOI 10.1093/bioinformatics/btac672
- GTDB itself: Parks DH *et al.* (2022) GTDB: an ongoing census of bacterial and archaeal
  diversity through a phylogenetically consistent, rank normalized and complete genome-based
  taxonomy. *Nucleic Acids Res* 50:D785–D794. DOI 10.1093/nar/gkab776
- **IN USE (Shannon): 2.7.2** at `<workdir>/envs/gtdbtk-2.7.2`, built by
  `workflow/shannon/bac120/b2_setup_env.sh`. Reference data: **GTDB r232** (94 GB) at
  `<db-root>/gtdbtk/release232` — `metadata/metadata.txt` reads `VERSION_DATA=r232`.
  A **matched** toolkit/reference pair.
- Superseded: Moriah had 2.4.0 + r226, a toolkit one patch *below* the minimum supported for its
  own reference data. The pipeline moved to Shannon (plan addendum, decisions S1–S5).
- **The mask is release-specific, so releases are NOT interchangeable.** `align` applies
  `masks/gtdb_r<ver>_bac120.mask`; a different release retains a different column set. This is
  acceptable here only because the tree is standalone and never compared column-wise against a
  GTDB reference tree. Report the alignment length; do not claim it matches another run.
- **A version mismatch does NOT crash on a missing mask** — a claim worth correcting because it
  is tempting and wrong. Verified in `gtdbtk/config/common.py:284`: the mask filename is
  *derived* from `VERSION_DATA` in the reference package's own `metadata.txt`
  (`f"gtdb_{self.VERSION_DATA}_bac120.mask"`), so an older toolkit resolves it too. Matching
  versions buys a supported pairing, not protection from a loud failure.
- **`identify` refuses inputs whose id matches a GTDB reference id.** For a Burkholderia set this
  is not an edge case: **68 of our 790 accessions are GTDB representative genomes**, because GTDB
  picks representatives from the same assemblies. The reference id set is the raw accession
  (`gtdbtk/tools.py:27`), so any prefix clears it. Fix used: symlink **all** inputs into
  `genomes_u/` as `u_<id>.fna`, strip `u_` from the MSA headers afterwards. Prefix all, never
  only the clashing ones — a two-class naming scheme is what silently mislabels tips later.
- Bundled deps in that env, all confirmed: **Prodigal V2.6.3**, **HMMER 3.4**,
  **FastTree 2.2.0 (double precision)**, `pplacer`, `skani`, `mash`.
- Subcommands (`gtdbtk` with no args prints them): workflows `classify_wf`, `de_novo_wf`;
  methods `identify`, `align`, `classify`, `infer`, `root`, `decorate`.
- **We use `identify` + `align` only, NOT `de_novo_wf`.** `de_novo_wf` chains
  identify→align→infer→root→decorate, and the `root`/`decorate` steps read the r226 taxonomy and
  RED files — exactly the parts a version mismatch corrupts, and exactly the parts we do not
  need. Marker HMMs are stable across releases. Rooting is done with this project's own
  ingroup-stem convention (`workflow/figures/tree_display.py`), not `--outgroup_taxon`.
- **`align --skip_gtdb_refs`** keeps the MSA to our genomes only. Without it, thousands of GTDB
  reference genomes join the alignment and the tree stops being about our genome set.
- Gene calling is GTDB-Tk's own Prodigal run, so **no Bakta annotation is needed** for this
  pipeline — a large saving versus the chromosome-1 route.
- bac120 = 120 bacterial single-copy marker proteins, concatenated amino acids, masked to
  ~5,000 columns. Single-copy is the property being bought: it is what stops a duplicated gene
  being mistaken for an orthologue (compare D17's single-copy-core convention).

## FastTree — approximate-ML tree at 800+ tips

- Repo/docs: http://www.microbesonline.org/fasttree/
- Paper: Price MN, Dehal PS, Arkin AP (2010) FastTree 2 — approximately maximum-likelihood trees
  for large alignments. *PLoS ONE* 5:e9490. DOI 10.1371/journal.pone.0009490
- On-disk version: **2.2.0 Double precision, OpenMP (128 threads)** at
  `<conda-root>/envs/anvio-9/bin/FastTreeMP` (Shannon); also in the gtdbtk env.
- **Double precision matters here and is not the default build.** Single-precision FastTree
  cannot resolve very short branches, and this genome set is full of them — the
  B. pseudomallei/mallei block spans 0.0055 subs/site across 232 genomes. Confirm the banner
  reads "Double precision" at run time; do not assume.
- Support values are **SH-like local supports, not bootstraps** (on by default; `-nosupport`
  disables). They are not comparable to the chromosome-1 tree's UFBoot and must not be reported
  as though they were.
- `FastTreeMP` parallelises with `OMP_NUM_THREADS`; speed-up is sublinear, so do not request a
  whole 128-core node for it.

## Settled — GTDB-Tk / FastTree / IQ-TREE on Shannon

- `check_install` was **deliberately not run**: the reference package is 190 GB, the check is slow,
  and the package was independently confirmed good. Verified instead, cheaply and specifically, that
  the toolkit resolves the files it will actually read — `VERSION_DATA=r232`,
  `gtdb_r232_bac120.mask` present, Pfam dir and TIGRFAM HMM present.
- **IQ-TREE: use 3.1.1**, installed alongside gtdbtk in the same env
  (`envs/gtdbtk-2.7.2/bin/iqtree`). Shannon also ships **3.0.1** at
  `envs/anvio-9/bin/iqtree3` — do NOT use it here: 3.1.1 is the build that made this project's
  chromosome-1 and pC3 trees, and those trees get compared to this one.

## Still open — GTDB-Tk / FastTree

- Confirm the bac120 marker count actually recovered for MF7 (a 63-contig draft) — the whole
  reason this tree can include MF7 when the chromosome-1 tree cannot.
