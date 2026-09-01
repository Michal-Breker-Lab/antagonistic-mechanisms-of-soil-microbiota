# Analysis decisions — pC3 rebuild (no dereplication, MF7 added)

The running record of every analysis choice made during the rebuild of the
*Burkholderia* sensu lato comparative-genomics pipeline on all 773 complete
genomes, with the reasoning and the evidence behind each. Decisions are numbered
D1–D20 in the order they were taken (2026-08-19 to 2026-08-25); each was reviewed
and signed off before the affected stage ran. Where a decision was later found to
be wrong, the correction is recorded in place rather than by rewriting history.

## The fact that drives all of it

99% ANI dereplication collapsed 771 genomes into 256 clusters. One cluster holds
**233 genomes (30.2% of the set)** — the *B. pseudomallei* / *B. mallei* complex.
472 genomes (61%) sit in clusters of >=5. 194 clusters are singletons.
That skew is sequencing effort on a biothreat pathogen, not biology.

## D1 — Keep every genome; use clone clusters as a covariate, not a filter

Nothing is removed. The retained `derep_cluster_membership.tsv` stops being a filter
and becomes a label, used to produce a companion clone-stratified estimate wherever a
statistic is sensitive to how many times a strain happened to be sequenced.
Affects: Fig 4 (Heaps' gamma), Fig 10/11 (randomisation nulls), host association,
SIMMAP rate estimates.

## D2 — Randomisation nulls computed twice

Unrestricted (as originally) and clone-stratified (draws match the observed
carrier set's cluster composition). Both p-values reported side by side; the
verdict stated explicitly when they disagree.

**Direction corrected 2026-08-19 after implementation.** The rationale first
given at the approval gate — that clone structure makes results look *more*
significant — is wrong. p is P(a random set is at least as tightly clustered as
the carriers), so only the null's LEFT tail matters, and a large clone cluster
fills that tail with all-clone draws of ~zero mean distance. That puts a FLOOR
under p. Measured across clone fractions from 10% to 80% of the tree, the
unrestricted p was >= the stratified p at every point, with the gap widening as
the clone fraction grew (`s29_null_stratified.py --selftest`, test 2), and
p_unrestricted tracked P(all tips drawn from the clone cluster) (test 3).

So the unrestricted null is **conservative**: the hazard is a real clustering
result being erased by an artefactual floor, not a false positive. The dual-null
design is unchanged and is if anything more clearly needed, since the direction
is not predictable from the tree alone.

**Magnitude for the actual data**, computed from the retained cluster table:

| Pool | n | Largest clone cluster | Floor on p (8 carriers) |
|---|---|---|---|
| All tips | 771 | 233 | 6.4e-05 |
| pC3-positive only | 254 | 28 | 8.1e-09 |

The 233-genome *B. pseudomallei* cluster contains **zero** pC3-positive genomes,
so it is absent from the pC3-positive pool entirely — and that pool is the
sharper test the report already designates as the one that isolates the warhead
from its replicon. Published warhead p = 0.0004 on the dereplicated tree; the
all-tips floor of 6.4e-05 is ~16% of that, so the result is expected to survive
with p drifting slightly upward, and the pC3-only test is effectively unaffected.

## D3 — Collapse clone clusters in tree *rendering* only

Every tip stays in the alignment, the tree and every statistic. Main-panel tree
figures draw large clone clusters as triangles annotated with member counts; the
fully expanded tree ships as a supplementary panel.

## D4 — Fix the substitution model to GTR+F+R10; skip ModelFinder

The 304-tip run selected GTR+F+R10 and cost 836 CPU-h / 18 h wall. Fixing it buys
comparability with the published tree and removes a large fraction of the cost at
771 taxa.

## D5 — Panaroo first, PPanGGOLiN as fallback, and if we fall back we re-run the old set too

Switching clustering software changes every pangenome number for reasons unrelated
to the genome set. If PPanGGOLiN is used at full scale, the 306-genome set is
re-run through it so old-vs-new is tool-matched.

## D6 — Same Bakta version/DB for all ~464 new annotations

Mixing Bakta and pyrodigal gene calls inside one Panaroo run is not permitted.

## D7 — MF7 eligibility rule

If MF7 is a draft (>~10 contigs), it enters the tree, pangenome and all
presence/absence searches, and is excluded with an explicit footnote from the
replicon-architecture panels, where replicon assignment is undeterminable.

## D8 — Refresh the NCBI snapshot to the rebuild date

Live check 2026-08-19: 1,525 GCA+GCF records vs 1,523 on 2026-08-03 — about one
new distinct genome. Trivial cost; makes "all available" true as stated.

---

## D9 — Re-deriving the pC3 diagnostic gene families

**Status: APPROVED 2026-08-21** (D9a clone-collapsed training, D9b self-consistency
check, D9c thresholds held at original absolute values + 95% sensitivity).

**Why this came up.** Reconstructing the pC3 calling path showed that
`tables/c3_diagnostic_orthogroups.txt` names its 1,013 families by Panaroo
cluster-representative ID (`GCA_024298805.1|CP101281.1|CFKEBF_07249`). Those IDs
are artefacts of one specific Panaroo run. A new run mints new ones, so the
retained lists cannot be carried over without a sequence-level remapping
exercise. Re-derivation is forced, not chosen.

**The original two-step logic, to be reproduced:**
1. Training set = genomes with exactly two large secondary replicons, where
   larger = chromosome 2 and smaller = pC3 is unambiguous by position (151
   genomes originally).
2. Families with >=90% of occurrences on the smaller replicon (min 10 genomes)
   -> pC3-diagnostic; >=90% on the larger -> chromosome-2-diagnostic.
3. Classify every large secondary replicon in every genome by which diagnostic
   set it carries; a genome is pC3-positive if it has >=1 pC3-classified replicon.

**D9a (proposed).** Derive the diagnostics from a **clone-collapsed training
set** -- one genome per 99% ANI cluster from `clone_cluster.tsv` -- then apply
the resulting classifier to all 773 genomes unweighted. Without derep the
training set becomes clone-dominated; a family carried by 200 near-identical Bcc
strains and nothing else would score as "diagnostic" on what is effectively one
observation. This is the same clone-inflation problem as D1/D2, reaching a step
the plan did not cover. Keeps "no dereplication" for the results while keeping
it out of the training.

**D9b (proposed).** Add a **self-consistency check the original did not report**:
how many training genomes does the classifier reclassify against their own
positional assignment? The bootstrap is only sound if that number is small.

**D9c.** Keep the >=90% and >=10-genome thresholds at their original absolute
values for comparability; report sensitivity at 95%.

**Consequence, independent of the above.** With all 773 genomes annotated, every
pC3 call is direct. The "466 of 771 calls are inherited" caveat in Section 5.3 --
and limitation #1 in the report's limitations list -- resolve rather than being
renumbered. The direct-only subset the report says to "quote if challenged"
becomes the whole table.

---

## D10 — MF7 handling (APPROVED 2026-08-21, option A)

**Finding.** MF7 is the same strain as MF6 (skani ANI = 100.00, clone cluster 22),
not an independent isolate. It additionally carries a real ~940 kb tandem
duplication on chr2 that MF6 lacks (two copies at 99.9878% identity, 115
mismatches; MF6-lineage read depth flat at 21.0x vs 22.5x, so MF6 did not
collapse an array). Full evidence:
a separate finding memo (not distributed with this repository).

**Decision (A).** MF7 remains a member of the full 773-genome set. The
pseudo-replication it introduces is absorbed by the 99% ANI clone-cluster
covariate (D9), consistent with the no-dereplication design applied to the other
232 clones in cluster 0. MF6/MF7 are NOT collapsed.

**Two binding constraints that follow:**

1. **MF7's assembly size (8,681,028 bp) never enters a size distribution,
   size-vs-architecture comparison, or replicon-size statistic.** It overstates
   the strain's genome by ~940 kb of duplicated sequence and is not a genome
   size. Where a size is required for this strain, MF6's is authoritative.
2. **MF7's architecture call (3 large, pC3-positive) is footnoted as derived by
   alignment to MF6**, not by the size-rank typing applied to the other 772
   genomes. Its 63-contig draft types as a meaningless `4+_large` under size-rank
   rules and must not be counted that way.

**Rejected alternative (B):** collapsing MF6/MF7 to one genome. Cleaner for
frequencies, but inconsistent with retaining every other clone in the set, and it
would discard the only within-strain structural comparison available.

**Open, not blocking:** a mixed MF7 culture would give the same signature as a
tandem array. MF7's reads would settle it (2x vs 1x depth over the interval) but
are not present on Moriah or the Drive. Reported as a stated limitation.

---

## D11 — Stage 4 tool correction and the pC3 classification rule (APPROVED 2026-08-22, option 3)

**Tool correction.** The pC3 diagnostics were rebuilt on **MMseqs2**, not Panaroo.
The project's own methods draft states it plainly: "Proteins from all large secondary
replicons were clustered into orthologous groups with MMseqs2 v18.8cc5c
(easy-cluster, 50% identity, 80% coverage)." Panaroo is the *pangenome* tool
(Stage 7, id 0.98) and PPanGGOLiN the *clustering-identity sensitivity* arms
(§6.2, id 0.80 / 0.60) -- three tools, three jobs. A 3h25m Panaroo run
was spent against the wrong stage before this was caught, because the tool was
inferred from output-file shapes rather than read from the project's own methods.

MMseqs2 on Moriah is **18.8cc5c** -- the exact version of the original run.

**Two implementation details recovered from the retained tables, not documented:**

1. `c3_content` / `c2_content` are the fraction of the **diagnostic set** a
   replicon carries, counted over **distinct families** -- not the fraction of
   the replicon's own families that are diagnostic. The two differ ~10x because
   the c2-diagnostic set is ~5x the c3 set; the wrong reading collapses
   `other_megaplasmid` into `chromosome2` (208 of 289 in testing).
2. The retained table's operational rule is
   `c2>=0.30 -> chromosome2; c3>=0.10 & c2<=0.07 -> c3; else other_megaplasmid`,
   which reproduces it at **0 mismatches on all 592 replicons**. Those absolute
   cuts do NOT transfer to the rebuild, because the original's 592 replicons were
   dereplicated and the rebuild's 1,219 are not -- the 233 pseudomallei clones
   alone shift every distribution the cuts were tuned against.

**Decision (option 3, hybrid).**

- The c3-vs-chromosome2 boundary is the **midpoint of the two reference means**
  of `score = c3_content - c2_content`, exactly as the methods draft specifies.
  Being a midpoint of the observed references it is self-calibrating, so it is
  immune to the content-scale difference above. Boundary -0.0816; accuracy on the
  reference sets **98.5% (pC3) / 97.1% (chromosome 2)**.
- The documented size-rank tie-breaker applies within +/-0.05 of the boundary.
- `other_megaplasmid` (carries neither signature) is retained via **one declared
  fitted parameter**: below the **20th percentile** of both reference signature
  distributions. 20 was chosen to maximise agreement with the retained
  classification -- **565/586 = 96.4%** on shared replicons. This is the only
  fitted parameter and must be stated as such in the Methods.

**Results.** 253 c3 replicons, 262 chromosome2, 704 other_megaplasmid of 1,219;
**253 of 773 genomes pC3-positive**, every call direct.

**D9b self-consistency, reported as two numbers.** A *hard inversion* (the
classifier swapping pC3 and chromosome 2) is the failure the check exists to
catch: **1 of 272, 0.37%** -- `GCF_003812585.1` NZ_CP033747.1, whose smaller
secondary replicon is genuinely chromosome-2-like (c2=0.4037 over 1,844
families). *Declining to call* (-> other_megaplasmid) is a legitimate outcome,
not a contradiction, and accounts for the remaining 53 of 272 (19.49%).

**D9c sensitivity passes outright:** the pC3-positive count is **253 at both
frac>=0.90 and frac>=0.95**. The headline number does not depend on that cut.

---

## D12 — Substitution model for both trees (APPROVED 2026-08-22)

**Decision: GTR+F+I+R5 for both the pC3 and chromosome-1 trees, with 1,000
ultrafast bootstrap replicates. The model is INHERITED from the original run's
ModelFinder output, not re-searched.**

**Why not re-run ModelFinder.** At 771 taxa x 1,337,868 sites the model search is
a large fraction of total runtime. The original's ModelFinder results survived
the Shannon outage on the Drive (`logs/c3_core.iqtree`, `logs/chr1_core.iqtree`)
and can simply be adopted -- faithful *and* free.

| Tree | Original data | ModelFinder best-fit (BIC) |
|---|---|---|
| pC3 core | 140 seq x 52,834 sites | **GTR+F+I+R5** |
| chr1 core | 304 seq x 490,720 sites | **GTR+F+R10** |

Both selected **GTR+F** -- general time-reversible matrix with empirical base
frequencies. `+F` matters here: these genomes average **66.6% GC** (range
58.5-69.0), so assuming equal base frequencies would be wrong by ~2x.

**Rate heterogeneity is strongly multimodal**, which is why the category count
matters. chr1's fitted R10 is:

    (0.4586, 0.0124)  <- 46% of sites effectively invariant
    (0.1207, 0.295) (0.0236, 0.566) (0.0393, 0.823) (0.0589, 0.880)
    (0.1119, 1.578) (0.1106, 2.690) (0.0598, 4.344)
    (0.0154, 7.025) (0.0014, 13.77) <- thin, very fast tail

A 4-category model cannot represent that shape. An earlier plan to use
`GTR+F+R4` was therefore **wrong and was corrected before the chr1 run** -- it
was under-parameterised, and the pC3 selection (`+I+R5`) also shows FreeRate did
NOT absorb the invariant class on this data, contrary to the usual expectation.

**Why I+R5 for chr1 rather than the original's R10.** R10 costs roughly 2.5x R4
per site and ~2x I+R5; on 771 x 1.34 Mb that is material. `+I` captures the 46%
near-invariant class explicitly with one parameter instead of spending free
categories on it. Stated in Methods as a deliberate, cost-motivated departure
from the chr1 selection; the pC3 tree uses the original's exact model.

**Also inherited:** 1,000 UFBoot replicates, as the methods draft specifies.
UFBoot >= 95 is the support threshold used in the figure captions (roughly as
stringent as classical bootstrap >= 70).

---

## D13 — Eight taxa dropped from the chromosome-1 tree (APPROVED 2026-08-23)

**Decision.** The chromosome-1 species tree is built on **763 taxa**, not 771:
eight sequences carrying **>90% gaps** in the core alignment are excluded.

| Genome | Gaps in chr1 core alignment |
|---|---|
| GCA_040954445.1 | 99.18% |
| GCF_050430075.1 | 95.97% |
| GCF_034047095.1 | 95.97% |
| GCF_039852015.1 | 95.62% |
| GCF_004842085.1 | 93.97% |
| **MF7** | **93.96%** |
| GCF_003812585.1 | 93.86% |
| GCF_003568605.1 | 91.06% |

A tip inferred from <10% of the alignment cannot be placed reliably and invites
long-branch artefacts. Figures 1, 3 and 12 all rest on this topology.

**Why MF7 in particular, and why it costs nothing.** MF7 is a 100% ANI clone of
MF6 (D10), and MF6 is present at a normal 24.4% gaps, so the strain is fully
represented. Dropping MF7 removes pseudo-replication, not information.

**Root cause, worth recording.** Replicon typing assigns exactly ONE contig per
genome as `chromosome1` (argmax ribosomal proteins). For MF7's 63-contig draft
that rule picked `BMF7_NODE_12` (216,555 bp, 32 ribosomal proteins) -- correct in
kind, since alignment confirms it is chr1 at 0.9992 -- but MF7's true chr1 is
**3,573,788 bp across 8 contigs**. The chr1 pangenome therefore received 216 kb
of 3.57 Mb, i.e. 6%, producing the 93.96% gaps. This is the draft-assembly
failure mode D10 anticipated; the D10 exclusion had been applied to figures via
`rebuild_rules.py` but NOT to pangenome set construction. Fixed here.

**Stage 4 is NOT affected.** MF7's misfiled chr1 fragments entered the secondary
-replicon pool, but the diagnostic classifier handled them correctly: its true
pC3 contigs score strongly (NODE_3 c3=0.393 vs c2=0.0005; NODE_8 0.060 vs
0.0005), its chr2 contig scores chromosome2 (c2=0.310), and the chr1 fragments
score weakly on both and fall to `other_megaplasmid`. **The 253 pC3-positive
count stands.**

**Note for the report.** `GCF_003812585.1` appears here AND as the single D9b
hard inversion -- two independent signals that this genome's architecture is
genuinely unusual, not a processing artefact. Worth a sentence.

---

## D14 — MF6 and MF7 count as RHS warhead carriers (APPROVED 2026-08-23)

The original excluded MF6 from the RHS search entirely: the warhead query
`CT354_CFFIHE_03684` *is* an MF6 locus, so MF6 scores 100%/100% against itself.
Its 771-genome denominator was NCBI-only. The rebuild searched all 773, so MF6
and MF7 (MF6's 100%-ANI clone) both appear as carriers and the count goes 7 -> 9.

**Decision: include them.** The association therefore reads 9/253 pC3-positive
vs 0/520 pC3-negative, OR inf, Fisher p = 3.91e-05.

**Methods must state the caveat.** Two of the nine carriers are the query strain
and its clone, so the ten-fold improvement over the original p (3.98e-04) is
partly circular and is NOT independent evidence. The biological claim is
unchanged and does not rest on it: zero pC3-negative carriers either way, and
excluding both still gives 7/251 vs 0/520, p ~ 4e-04. Panel C shows MF6 and MF7
at 100.0% identity, which makes the self-hit visible rather than hidden.

**Two constants independently re-derived.** fig9 hard-coded QLEN = 2867 and
LAST_SHARED = 2645 (last residue any CT-negative genome aligns to), plus "247
CT-negative genomes". The rebuild recomputes all three from a different genome
set and gets exactly 2867, 2645 and 247. They are now read from the coverage
profile rather than hard-coded.

**Large-replicon threshold recovered: 300,000 bp.** Undocumented in the methods.
It reproduces the original `architecture` label for 770 of 771 genomes; 295 kb
and 305 kb each give 3 errors. The single disagreement at 300 kb is
GCF_054166145.1 (original `0_large`, contigs 3.97 Mb / 2.96 Mb / 454 kb /
331 kb) -- the known-bad label already flagged for correction.

---

## D15 — MF7 joins Set B, the *B. sola* clade (APPROVED 2026-08-23)

Set B is MF6's ANI >=95% clade, the basis of the setB_* functional tables and
the screens behind figs 7-8. MF7 sits at ANI 100.00, so the rule admits it:
Set B goes 5 -> 6 genomes (MF6, MF7, GCF_016899425.1, GCF_905400185.1,
GCF_053038975.1, GCF_053209605.1).

**The initial objection, later overruled.** The objection was that Set B's
core is an INTERSECTION, so adding a 63-contig draft of the same strain can only
hold the core flat or push it down, and a family lost to a contig break would be
recorded as "not conserved in the clade" when it is really "not assembled".

**The reason for including it: MF7 was isolated separately.** That is a
distinct biological sampling event, not a re-sequencing of MF6, so it is an
independent observation of what a *B. sola* pC3 carries -- which is exactly the
question Set B asks. Clone-ness in ANI does not make it a duplicate record.

**The assembly-fragmentation hazard is real but must be handled by reporting,
not by exclusion.** Any family present in the other five and absent only in MF7
must be checked against MF7's contig boundaries before being called a loss.
Report the Set B core both ways -- with and without MF7 -- so the reader can see
how much of any change is assembly rather than biology. Contrast with D13, where
MF7 was dropped from the chr1 TREE: there, fragmentation corrupted the alignment
itself; here it only risks a false absence, which is measurable.

---

## D16 — MF7's pC3 and chr1 are COMPLETE; s36/s37 carved one contig each (2026-08-23)

**Supersedes the fragmentation premise of D15.**

`s36_setb_definition.py` wrote `setB_members.tsv` with a single `pc3_contig`
column and `s37_setb_inputs.py` carved that one contig. For the five closed Set B
genomes one contig *is* the replicon, so nothing was lost. MF7 is a draft
assembly, and its replicons are split:

| MF7 replicon | contigs | total bp | CDS | vs MF6 |
|---|---|---|---|---|
| pC3 | 2 (NODE_3 + NODE_8) | 1,176,232 | 1,011 | 99.7% bp, 99.9% CDS |
| chr1 | 8 | 3,573,788 | 3,247 | 99.5% bp, 99.9% CDS |

Everything downstream inherited the omission. The "contiguous 253-CDS block
missing from MF7" measured on 2026-08-23 **is NODE_8**: MF6 was being diffed
against 73% of MF7's pC3.

**Consequences.**
1. D15's stated hazard — MF7's fragmentation shrinking the clade core — was
   false. With the complete pC3, MF7 costs the core **zero** families (833 both
   ways); chr1 core moves 2,833 → 2,832.
2. MF7's fig8 screen results from job 45911143 are incomplete by construction
   (5 antiSMASH regions vs MF6's 7; 21 MacSyFinder rows vs 37; 759 InterProScan
   proteins instead of 1,011). That spread was previously reported as MF7's
   "assembly-fragmentation signature across three independent measures". It is
   one missing file, three times. `m15_mf7_screens.sbatch` re-runs MF7 alone.
3. D13 dropped MF7 from the chr1 tree at 93.96% gaps. That is an artefact of the
   core-genome **alignment**, not of missing sequence — MF7 carries 99.5% of
   chr1. D13 deserves re-examination on those grounds.
4. Replicon membership is now `setB_contigs.tsv` (`s40_setb_contigs.py`), and
   inputs are carved by `s41_setb_inputs.py`, which takes a contig *list*.

**Decision (2026-08-23):** keep the five closed genomes as the primary
core definition with MF7 reported separately, and additionally produce a
6-genome version for reference. Both are built; they agree to within one family.

## D17 — Set B clade-core method, recovered by reproduction (2026-08-23)

The generator behind `tables/setB_*` was lost and its settings were never
documented. Recovered by sweeping MMseqs2 `easy-cluster` over the five closed
genomes until the retained counts reproduced:

- **`--min-seq-id 0.8 -c 0.8 --cov-mode 1`** (PPanGGOLiN's defaults) — the unique
  cell giving 1,382 pC3 families.
- **Core is SINGLE-COPY core**: present in every genome with exactly one gene
  each. This is what makes the split 833 + 549 rather than 843 + 539 — the 10
  all-genome families that are multi-copy in at least one genome were counted as
  accessory. Confirmed exactly.
- **Family annotation by STRICT majority** (> half the members) for COG and KO,
  which are controlled vocabularies; **most-common** for the free-text product,
  whose wording varies between genomes. Strict-majority throughout gave 30.9%
  COG coverage against the retained 28.7%; the split rule gives 28.3%.
- Thread count does not affect the result (1/4/16 threads all give 2,833).

**Two residuals not reproduced**, both documented rather than forced:

- **chr1 core 2,833 vs the retained 2,847** (0.5%). Not thread nondeterminism.
  Immaterial to panels A–B, which are percentages of ~2,800.
- **genus-derived 62 vs the retained 50.** The bridge to the genus-wide pC3 core
  is through MF6's locus tags, but the original's genus reference set is
  unrecoverable: PPanGGOLiN id 0.60 soft-core gives 75, id 0.80 gives 66, and
  the Panaroo ≥95%-of-253-carriers soft-core (55 families) gives 62. The last is
  used, being the same Panaroo matrix the rest of the report uses. The report
  sentence "94% of the pC3 core is clade-specific" becomes 93%.

**Panel D is carried over unchanged.** It is computed on MF6 alone (pC3 vs the
rest of MF6's genome), so neither MF7 nor the re-clustering can touch it.

**One caption correction.** Panel A's legend read "pC3 clade core (n = 833) /
chromosome 1 clade core (n = 2,847)", but the panel plots percentages of
COG-ANNOTATED families, whose totals are 236 and 1,896. The legend is now
data-driven and reports those.

## D18 — MF7's screens re-run; it is indistinguishable from MF6 (2026-08-23)

Job 45921383 (`m15_mf7_screens.sbatch`, 16 cores, COMPLETED 00:42:51)
re-screened MF7 against its complete two-contig pC3. Every measure that had
looked like MF7 "fragmentation" now matches MF6:

| | MF7 before | MF7 after | MF6 |
|---|---|---|---|
| antiSMASH regions | 5 | **7** | 7 |
| secretion systems | (21 rows) | **7** | 7 |
| InterProScan proteins with signatures | 759 in | **923** | 922 |

Panel B of fig8 now shows MF7's column identical to MF6's across all nine region
types, which is what two assemblies of one strain at ANI 100.00 must give. The
earlier three-measure "fragmentation signature" was one missing contig, counted
three times.

Screen tables are rebuilt by `s44_screen_tables.py`, which reads MF7 from a
second directory and merges its per-contig MacSyFinder runs. `s45_ips_coverage.py`
(Drive-side) joins the per-protein screen to the gene families for fig8 panel A,
using the same strict-majority rule as s43; it gives 88.7% annotated, reproducing
the report's 89%.

---

## D19 — The chromosome-1 typer mis-calls 5 genomes; documented, not fixed (2026-08-25)

**The defect.** Replicon typing picks chromosome 1 as `argmax(ribosomal_proteins)`
over a genome's contigs. In five genomes the **ribosomal superoperon and *dnaA*
sit on different replicons**, and the rule follows the superoperon — demoting the
real chromosome to `secondary_large`.

| genome | typed as chr1 | CDS | r-prot | the *dnaA* contig | CDS | r-prot | on tree? |
|---|---|---|---|---|---|---|---|
| GCF_003812585.1 | 1,030,740 | 919 | 33 | 3,239,301 | 2,926 | 20 | dropped |
| GCF_003568605.1 | 2,735,759 | 2,261 | 34 | 3,986,340 | 3,394 | 19 | dropped |
| GCF_004842085.1 | 3,363,129 | 2,645 | 33 | 3,800,466 | 3,226 | 20 | dropped |
| MF7             |   216,555 |  203 | 32 |   263,971 |   242 |  1 | dropped |
| GCF_033870355.1 | 2,668,693 | 2,299 | 44 | 3,152,963 | 2,675 |  9 | **ON TREE** |

The affected genomes are not defective — `GCF_003812585.1` is a closed,
complete 3-replicon assembly. The rule is.

**Two-step consequence.** Wrong replicon -> wrong gene set into the chr1
pangenome -> the genome is absent from nearly every core family -> its row in
the core alignment is >90% gaps -> D13 drops it as unalignable. The cause was
invisible at the point the symptom was caught, which is why D13 recorded these
as divergent taxa.

**So half of D13's eight exclusions are this bug, not divergence.** The other
four are genuine: three *Robbsia* (a different genus, formerly *B. andropogonis*)
and `GCA_040954445.1`, a 2.5 Mb assembly where a *Burkholderia* should be 6–9 Mb.

**`GCF_033870355.1` is the one that got through.** Its wrong replicon still
shared 623/1,018 alignment families (61.2%), clearing D13's >90%-gap threshold,
so nothing flagged it. Its tip on the 763-taxon tree is inferred from the wrong
replicon. Retained genomes hold a median 1,018/1,018.

**Decision (2026-08-25): document and move on — option 1 of three.**
All four affected genomes sit at **84.7–89.0% ANI to MF6** (*B. mallei*,
*B. multivorans*, *B. thailandensis*, *B. pseudomallei*). None is in Set B, none
is inside the >=94% clade of Figure 12, and none is near the pC3, warhead or
toxin analyses. What is affected is the genus-wide **backbone** of Figures 1 and
3: one tip in 763 built from the wrong replicon, and four legitimate
*Burkholderia* missing from a genus survey.

Rejected: (2) pruning `GCF_033870355.1` from the existing tree, and (3) fixing
the rule and rebuilding the chr1 pangenome, alignment and tree (~2-3 days of compute,
and it would invalidate Figures 1, 3, 10, 11, 12 meanwhile).

**The fix, when it is worth doing.** Give *dnaA* precedence over the ribosomal
superoperon: chromosome 1 is the contig carrying *dnaA*, using r-protein count
only to break ties or when no contig carries *dnaA*. CLAUDE.md's D2 already
describes chr1 as "*dnaA* + the ribosomal-protein superoperon"; the code treats
them as interchangeable, and they are not.

## D20 — Two genomes never entered the chr1 alignment (2026-08-25)

`GCA_059696275.1` and `GCF_050955445.1` hold **zero** of the 1,018 alignment
families and appear on no tree. The taxon count reconciles as
773 − 2 (these) − 8 (D13) = **763**, so nothing is wrong — but D13 documents
only the eight, leaving the report unable to explain 773 → 763. Recorded here so
the Methods can account for all ten.
