#!/usr/bin/env python3
"""
s17_functional.py - functional composition of the B. sola (Set B) clade core.

Annotation source: Bakta v1.12.0 / DB v6.0 DbXrefs (COG id + COG category, KEGG KO,
EC, Pfam). No re-annotation: every genome in the project carries the same Bakta run.

Contrasts (all approved 2026-08-10):
  A  pC3 clade core        vs  chr1 clade core        - what pC3 is for
  B  pC3 clade core        vs  pC3 clade accessory    - conserved vs variable on pC3
  C  pC3 clade core, genus-core-derived vs clade-specific - what the clade added

Statistics are DESCRIPTIVE: the two sets are two replicons of the same 5 genomes,
not independent draws, so log2 odds ratios with CIs lead and BH-corrected Fisher
is supporting detail only.
"""
import csv, os, re, sys, json, subprocess
from collections import defaultdict, Counter
from math import log, sqrt, log2

import numpy as np
from scipy.stats import fisher_exact, norm

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
KEGGDIR = "/mnt/LargeStorageNoBackup/Datasets/Moshea/Databases/anvioKEGG/modules"
OUT = f"{W}/results"
os.makedirs(OUT, exist_ok=True)

SETB = ["GCF_016899425.1", "GCF_905400185.1", "GCF_053038975.1", "GCF_053209605.1", "MF6"]
MODULE_THRESHOLD = 0.75          # anvi'o default
GENUS_CORE_FRAC = 0.95           # matches the project's "core" definition

csv.field_size_limit(10 ** 8)

# --------------------------------------------------------------------------
# 1. Bakta annotation: locus_tag -> functional labels
# --------------------------------------------------------------------------
COG_ID_RE = re.compile(r"^COG\d+$")

def load_bakta(genome):
    """locus_tag -> dict(product, gene, cog_id, cog_cats, kos, ec, pfam)"""
    path = f"{W}/annot/{genome}/{genome}.tsv"
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[1] != "cds":
                continue
            locus, gene, product = f[5], f[6], f[7]
            dbx = f[8] if len(f) > 8 else ""
            cog_id, cog_cats, kos, ec, pfam = None, set(), set(), set(), set()
            for tok in dbx.split(","):
                tok = tok.strip()
                if tok.startswith("COG:"):
                    v = tok[4:]
                    if COG_ID_RE.match(v):
                        cog_id = v
                    else:                       # functional category letter(s)
                        cog_cats.update(ch for ch in v if ch.isalpha())
                elif tok.startswith("KEGG:"):
                    kos.add(tok[5:])
                elif tok.startswith("EC:"):
                    ec.add(tok[3:])
                elif tok.startswith("PFAM:"):
                    pfam.add(tok[5:])
            out[locus] = dict(product=product, gene=gene, cog_id=cog_id,
                              cog_cats=cog_cats, kos=kos, ec=ec, pfam=pfam)
    return out

print("[1] loading Bakta annotation for 5 Set B genomes ...", flush=True)
ANN = {}
for g in SETB:
    a = load_bakta(g)
    ANN.update(a)
    print(f"    {g}: {len(a)} CDS", flush=True)

# --------------------------------------------------------------------------
# 2. Pangenome matrices: family -> {genome: [locus_tags]}
# --------------------------------------------------------------------------
def load_matrix(path):
    fams = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        genomes = hdr[14:]
        for row in rd:
            fam = row[0]
            cells = {}
            for gname, cell in zip(genomes, row[14:]):
                toks = [t for t in re.split(r"[\s;,]+", cell.strip()) if t]
                if toks:
                    cells[gname] = toks
            fams[fam] = cells
    return fams, genomes

print("[2] loading pangenome matrices ...", flush=True)
c3_fams, c3_genomes = load_matrix(f"{W}/ppanggolin/neighbours/setB_c3_id080/out/matrix.csv")
chr1_fams, chr1_genomes = load_matrix(f"{W}/ppanggolin/neighbours/setB_chr1_id080/out/matrix.csv")
print(f"    pC3  families {len(c3_fams)} over {len(c3_genomes)} genomes", flush=True)
print(f"    chr1 families {len(chr1_fams)} over {len(chr1_genomes)} genomes", flush=True)

def split_core(fams, genomes):
    core, acc = [], []
    for fam, cells in fams.items():
        (core if len(cells) == len(genomes) else acc).append(fam)
    return core, acc

c3_core, c3_acc = split_core(c3_fams, c3_genomes)
chr1_core, chr1_acc = split_core(chr1_fams, chr1_genomes)
print(f"    pC3  core {len(c3_core)}  accessory {len(c3_acc)}", flush=True)
print(f"    chr1 core {len(chr1_core)} accessory {len(chr1_acc)}", flush=True)

# --------------------------------------------------------------------------
# 3. Genus bridge: which clade-core families derive from the genus-wide core?
#    MF6 is present in BOTH the 140-genome run and the 5-genome run, so its
#    locus tags are the bridge between two independent clusterings.
# --------------------------------------------------------------------------
print("[3] bridging to the 140-genome pC3 pangenome via MF6 locus tags ...", flush=True)
genus_core_mf6_loci = set()
n_genus_core = 0
with open(f"{W}/ppanggolin/c3_id080/out/matrix.csv", newline="") as fh:
    rd = csv.reader(fh)
    hdr = next(rd)
    gnames = hdr[14:]
    mf6_i = gnames.index("MF6")
    ngen = len(gnames)
    need = GENUS_CORE_FRAC * ngen
    for row in rd:
        cells = row[14:]
        present = sum(1 for c in cells if c.strip())
        if present >= need:
            n_genus_core += 1
            for t in re.split(r"[\s;,]+", cells[mf6_i].strip()):
                if t:
                    genus_core_mf6_loci.add(t)
print(f"    genus soft-core families (>={GENUS_CORE_FRAC:.0%} of {ngen}): {n_genus_core}", flush=True)
print(f"    MF6 loci in genus core: {len(genus_core_mf6_loci)}", flush=True)

def is_genus_derived(fam):
    for t in c3_fams[fam].get("MF6", []):
        if t in genus_core_mf6_loci:
            return True
    return False

c3_core_genus = [f for f in c3_core if is_genus_derived(f)]
c3_core_clade = [f for f in c3_core if not is_genus_derived(f)]
print(f"    clade core: {len(c3_core_genus)} genus-derived, {len(c3_core_clade)} clade-specific",
      flush=True)

# --------------------------------------------------------------------------
# 4. Family-level functional labels (majority vote across members)
# --------------------------------------------------------------------------
def family_label(cells):
    cats, kos, cogids, prods = Counter(), Counter(), Counter(), Counter()
    n_members = 0
    n_with_cog = 0
    for g, loci in cells.items():
        for lt in loci:
            a = ANN.get(lt)
            if a is None:
                continue
            n_members += 1
            if a["cog_cats"]:
                n_with_cog += 1
            for c in a["cog_cats"]:
                cats[c] += 1
            for k in a["kos"]:
                kos[k] += 1
            if a["cog_id"]:
                cogids[a["cog_id"]] += 1
            prods[a["product"]] += 1
    half = n_members / 2.0
    # a category / KO is assigned to the family if carried by a majority of members
    fam_cats = {c for c, n in cats.items() if n > half} or set()
    fam_kos = {k for k, n in kos.items() if n > half}
    return dict(
        cats=fam_cats,
        kos=fam_kos,
        cog_id=cogids.most_common(1)[0][0] if cogids else "",
        product=prods.most_common(1)[0][0] if prods else "",
        n_members=n_members,
    )

print("[4] labelling families ...", flush=True)
LAB = {}
for fam in c3_fams:
    LAB[("c3", fam)] = family_label(c3_fams[fam])
for fam in chr1_fams:
    LAB[("chr1", fam)] = family_label(chr1_fams[fam])

SETS = {
    "pC3_clade_core":       [("c3", f) for f in c3_core],
    "pC3_clade_accessory":  [("c3", f) for f in c3_acc],
    "chr1_clade_core":      [("chr1", f) for f in chr1_core],
    "pC3_core_genus_derived": [("c3", f) for f in c3_core_genus],
    "pC3_core_clade_specific": [("c3", f) for f in c3_core_clade],
}

# --------------------------------------------------------------------------
# 5. Annotation coverage - reported as a RESULT, not silently dropped
# --------------------------------------------------------------------------
COG_NAMES = {
 "J":"Translation, ribosomal structure and biogenesis","A":"RNA processing and modification",
 "K":"Transcription","L":"Replication, recombination and repair",
 "B":"Chromatin structure and dynamics","D":"Cell cycle control, cell division",
 "Y":"Nuclear structure","V":"Defense mechanisms","T":"Signal transduction mechanisms",
 "M":"Cell wall/membrane/envelope biogenesis","N":"Cell motility","Z":"Cytoskeleton",
 "W":"Extracellular structures","U":"Intracellular trafficking and secretion",
 "O":"Posttranslational modification, protein turnover, chaperones",
 "X":"Mobilome: prophages, transposons",
 "C":"Energy production and conversion","G":"Carbohydrate transport and metabolism",
 "E":"Amino acid transport and metabolism","F":"Nucleotide transport and metabolism",
 "H":"Coenzyme transport and metabolism","I":"Lipid transport and metabolism",
 "P":"Inorganic ion transport and metabolism",
 "Q":"Secondary metabolites biosynthesis, transport and catabolism",
 "R":"General function prediction only","S":"Function unknown",
}

cov_rows = []
for sname, keys in SETS.items():
    n = len(keys)
    n_cog = sum(1 for k in keys if LAB[k]["cats"])
    n_ko = sum(1 for k in keys if LAB[k]["kos"])
    n_hyp = sum(1 for k in keys if "hypothetical" in LAB[k]["product"].lower())
    cov_rows.append(dict(gene_set=sname, n_families=n,
                         n_with_COG=n_cog, pct_with_COG=round(100*n_cog/n, 1) if n else 0,
                         n_with_KO=n_ko, pct_with_KO=round(100*n_ko/n, 1) if n else 0,
                         n_hypothetical=n_hyp,
                         pct_hypothetical=round(100*n_hyp/n, 1) if n else 0))

with open(f"{OUT}/setB_annotation_coverage.tsv", "w") as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(cov_rows[0].keys()), delimiter="\t")
    wtr.writeheader()
    wtr.writerows(cov_rows)
print("[5] annotation coverage:", flush=True)
for r in cov_rows:
    print(f"    {r['gene_set']:26s} n={r['n_families']:5d}  COG {r['pct_with_COG']:5.1f}%  "
          f"KO {r['pct_with_KO']:5.1f}%  hypothetical {r['pct_hypothetical']:5.1f}%", flush=True)

# --------------------------------------------------------------------------
# 6. COG category profile + descriptive effect sizes
# --------------------------------------------------------------------------
def profile(keys, restrict_annotated):
    """category -> count. Families may carry >1 category (standard COG practice)."""
    ks = [k for k in keys if LAB[k]["cats"]] if restrict_annotated else keys
    cnt = Counter()
    for k in ks:
        cats = LAB[k]["cats"]
        if cats:
            for c in cats:
                cnt[c] += 1
        else:
            cnt["-"] += 1          # explicit 'no COG' category
    return cnt, len(ks)

def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, p[i] * n / r)
        q[i] = prev
    return q

CONTRASTS = [
    ("A_pC3core_vs_chr1core",   "pC3_clade_core",         "chr1_clade_core"),
    ("B_pC3core_vs_pC3acc",     "pC3_clade_core",         "pC3_clade_accessory"),
    ("C_genusderived_vs_cladespecific", "pC3_core_genus_derived", "pC3_core_clade_specific"),
]

enrich_rows = []
for mode in ("all_families", "COG_annotated_only"):
    restrict = (mode == "COG_annotated_only")
    for cname, sa, sb in CONTRASTS:
        pa, na = profile(SETS[sa], restrict)
        pb, nb = profile(SETS[sb], restrict)
        cats = sorted(set(pa) | set(pb))
        rows, pvals = [], []
        for c in cats:
            a = pa.get(c, 0); b = na - a
            cc = pb.get(c, 0); d = nb - cc
            # Haldane-Anscombe correction keeps the OR finite at zero cells
            lor = log((a + .5) * (d + .5) / ((b + .5) * (cc + .5)))
            se = sqrt(1/(a+.5) + 1/(b+.5) + 1/(cc+.5) + 1/(d+.5))
            lo, hi = lor - 1.96*se, lor + 1.96*se
            _, p = fisher_exact([[a, b], [cc, d]])
            rows.append(dict(mode=mode, contrast=cname, cog_category=c,
                             cog_name=COG_NAMES.get(c, "no COG assigned" if c == "-" else c),
                             set_a=sa, n_a=a, total_a=na, pct_a=round(100*a/na, 2) if na else 0,
                             set_b=sb, n_b=cc, total_b=nb, pct_b=round(100*cc/nb, 2) if nb else 0,
                             log2_OR=round(lor/log(2), 3),
                             ci_low=round(lo/log(2), 3), ci_high=round(hi/log(2), 3),
                             p_fisher=p))
            pvals.append(p)
        for r, q in zip(rows, bh(pvals)):
            r["q_BH"] = q
        enrich_rows.extend(rows)

with open(f"{OUT}/setB_cog_enrichment.tsv", "w") as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(enrich_rows[0].keys()), delimiter="\t")
    wtr.writeheader()
    wtr.writerows(enrich_rows)

# raw profile table for plotting
prof_rows = []
for sname, keys in SETS.items():
    cnt, n = profile(keys, False)
    for c, v in sorted(cnt.items()):
        prof_rows.append(dict(gene_set=sname, cog_category=c,
                              cog_name=COG_NAMES.get(c, "no COG assigned" if c == "-" else c),
                              n_families=v, total=n, pct=round(100*v/n, 3) if n else 0))
with open(f"{OUT}/setB_cog_profile.tsv", "w") as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(prof_rows[0].keys()), delimiter="\t")
    wtr.writeheader()
    wtr.writerows(prof_rows)
print("[6] COG profile + enrichment written", flush=True)

# --------------------------------------------------------------------------
# 7. Per-family table (the auditable record behind every number above)
# --------------------------------------------------------------------------
with open(f"{OUT}/setB_functional_families.tsv", "w") as fh:
    fh.write("gene_set\tfamily\tcog_id\tcog_categories\tkos\tn_members\tproduct\n")
    for sname in ("pC3_clade_core", "pC3_clade_accessory", "chr1_clade_core"):
        for k in SETS[sname]:
            L = LAB[k]
            tag = ""
            if sname == "pC3_clade_core":
                tag = "genus_derived" if k in set(SETS["pC3_core_genus_derived"]) else "clade_specific"
            fh.write(f"{sname}{('/' + tag) if tag else ''}\t{k[1]}\t{L['cog_id']}\t"
                     f"{''.join(sorted(L['cats']))}\t{','.join(sorted(L['kos']))}\t"
                     f"{L['n_members']}\t{L['product']}\n")

# --------------------------------------------------------------------------
# 8. KEGG module stepwise completeness (anvi'o definition, threshold 0.75)
#    KOs come from Bakta, NOT from KOfam HMMs - see TOOLS.md caveat.
# --------------------------------------------------------------------------
print("[8] KEGG module completeness ...", flush=True)

def parse_module(path):
    rec, key = {}, None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(" "):
            if key:
                rec[key] += " " + line.strip()
        else:
            parts = line.rstrip("\n").split(None, 1)
            if not parts:
                continue
            key = parts[0]
            rec[key] = parts[1].strip() if len(parts) > 1 else ""
    return rec

def strip_optional(s):
    """remove '-K#####', '-M#####' and '-(...)' non-essential components"""
    out, i = [], 0
    while i < len(s):
        if s[i] == "-":
            j = i + 1
            if j < len(s) and s[j] == "(":
                depth = 0
                while j < len(s):
                    if s[j] == "(":
                        depth += 1
                    elif s[j] == ")":
                        depth -= 1
                        if depth == 0:
                            j += 1
                            break
                    j += 1
            else:
                while j < len(s) and (s[j].isalnum() or s[j] == "_"):
                    j += 1
            i = j
        else:
            out.append(s[i]); i += 1
    return "".join(out)

def split_top(s, seps):
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if depth == 0 and ch in seps:
            parts.append("".join(cur)); cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip() != ""]

def eval_expr(s, kos, mstatus):
    s = s.strip()
    if not s or s == "--":
        return False
    parts = split_top(s, ",")
    if len(parts) > 1:
        return any(eval_expr(p, kos, mstatus) for p in parts)
    parts = split_top(s, " +")
    if len(parts) > 1:
        return all(eval_expr(p, kos, mstatus) for p in parts)
    if s.startswith("(") and s.endswith(")") and len(split_top(s[1:-1], "")) >= 0:
        inner = s[1:-1]
        depth = 0; ok = True
        for i, ch in enumerate(inner):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth < 0:
                    ok = False; break
        if ok and depth == 0:
            return eval_expr(inner, kos, mstatus)
    if re.match(r"^M\d{5}$", s):
        return mstatus.get(s, False)
    return s in kos

MODULES = {}
for fn in sorted(os.listdir(KEGGDIR)):
    if not re.match(r"^M\d{5}$", fn):
        continue
    rec = parse_module(os.path.join(KEGGDIR, fn))
    d = rec.get("DEFINITION", "")
    if not d:
        continue
    MODULES[fn] = dict(name=rec.get("NAME", ""), klass=rec.get("CLASS", ""),
                       definition=d, steps=split_top(strip_optional(d), " "))

print(f"    parsed {len(MODULES)} KEGG modules", flush=True)

def stepwise_completeness(kos):
    mstatus = {}
    for _ in range(3):                      # resolve nested M-references
        newstat = {}
        for m, rec in MODULES.items():
            steps = rec["steps"]
            if not steps:
                newstat[m] = False; continue
            done = sum(1 for st in steps if eval_expr(st, kos, mstatus))
            newstat[m] = (done / len(steps)) >= MODULE_THRESHOLD
        if newstat == mstatus:
            break
        mstatus = newstat
    res = {}
    for m, rec in MODULES.items():
        steps = rec["steps"]
        if not steps:
            continue
        done = sum(1 for st in steps if eval_expr(st, kos, mstatus))
        res[m] = (done, len(steps), done / len(steps))
    return res

KOSETS = {}
for sname in ("pC3_clade_core", "chr1_clade_core", "pC3_clade_accessory"):
    ks = set()
    for k in SETS[sname]:
        ks |= LAB[k]["kos"]
    KOSETS[sname] = ks
    print(f"    {sname}: {len(ks)} distinct KOs", flush=True)

COMP = {s: stepwise_completeness(k) for s, k in KOSETS.items()}

with open(f"{OUT}/setB_kegg_module_completeness.tsv", "w") as fh:
    fh.write("module\tname\tclass\tn_steps\t"
             "pC3_core_steps\tpC3_core_completeness\t"
             "chr1_core_steps\tchr1_core_completeness\tpC3_complete\tchr1_complete\n")
    for m in sorted(MODULES):
        if m not in COMP["pC3_clade_core"]:
            continue
        a = COMP["pC3_clade_core"][m]
        b = COMP["chr1_clade_core"][m]
        fh.write(f"{m}\t{MODULES[m]['name']}\t{MODULES[m]['klass']}\t{a[1]}\t"
                 f"{a[0]}\t{a[2]:.4f}\t{b[0]}\t{b[2]:.4f}\t"
                 f"{int(a[2] >= MODULE_THRESHOLD)}\t{int(b[2] >= MODULE_THRESHOLD)}\n")

n_pc3_complete = sum(1 for m in COMP["pC3_clade_core"]
                     if COMP["pC3_clade_core"][m][2] >= MODULE_THRESHOLD)
n_chr1_complete = sum(1 for m in COMP["chr1_clade_core"]
                      if COMP["chr1_clade_core"][m][2] >= MODULE_THRESHOLD)
print(f"    modules >={MODULE_THRESHOLD:.0%} complete: pC3 core {n_pc3_complete}, "
      f"chr1 core {n_chr1_complete}", flush=True)

# --------------------------------------------------------------------------
# 9. Redundancy: does the pC3 core duplicate chromosome-1 functions?
# --------------------------------------------------------------------------
print("[9] pC3-core vs chr1 homology (DIAMOND, MF6) ...", flush=True)
DIAMOND = f"{W}/envs/bakta/bin/diamond"
tmp = f"{W}/tmp/func"
os.makedirs(tmp, exist_ok=True)

def gff_loci(path):
    s = set()
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS":
            m = re.search(r"locus_tag=([^;]+)", f[8])
            if m:
                s.add(m.group(1))
    return s

mf6_c3_loci = gff_loci(f"{W}/pangenome/gff_c3/MF6.gff3")
mf6_chr1_loci = gff_loci(f"{W}/pangenome/gff_chr1/MF6.gff3")
print(f"    MF6 CDS: pC3 {len(mf6_c3_loci)}, chr1 {len(mf6_chr1_loci)}", flush=True)

core_mf6_loci = set()
for f in c3_core:
    core_mf6_loci.update(c3_fams[f].get("MF6", []))

def write_faa(loci, path):
    keep, out = False, open(path, "w")
    n = 0
    for line in open(f"{W}/annot/MF6/MF6.faa"):
        if line.startswith(">"):
            lt = line[1:].split()[0]
            keep = lt in loci
            if keep:
                n += 1
        if keep:
            out.write(line)
    out.close()
    return n

nq = write_faa(core_mf6_loci, f"{tmp}/pc3_core.faa")
ns = write_faa(mf6_chr1_loci, f"{tmp}/chr1.faa")
print(f"    query {nq} pC3-core proteins vs {ns} chr1 proteins", flush=True)

subprocess.run([DIAMOND, "makedb", "--in", f"{tmp}/chr1.faa", "-d", f"{tmp}/chr1",
                "--quiet"], check=True)
subprocess.run([DIAMOND, "blastp", "-q", f"{tmp}/pc3_core.faa", "-d", f"{tmp}/chr1",
                "-o", f"{tmp}/pc3_vs_chr1.tsv", "--id", "30", "--query-cover", "70",
                "--subject-cover", "70", "-e", "1e-5", "--max-target-seqs", "1",
                "--threads", "16", "--quiet"], check=True)

hits = set()
for line in open(f"{tmp}/pc3_vs_chr1.tsv"):
    hits.add(line.split("\t")[0])
n_dup = len(hits)
print(f"    pC3-core proteins with a chr1 homologue (>=30% id, >=70% cov): "
      f"{n_dup}/{nq} ({100*n_dup/nq:.1f}%)", flush=True)

# --------------------------------------------------------------------------
# 10. Summary
# --------------------------------------------------------------------------
summary = dict(
    setB_genomes=SETB,
    pC3_clade_core=len(c3_core), pC3_clade_accessory=len(c3_acc),
    chr1_clade_core=len(chr1_core),
    pC3_core_genus_derived=len(c3_core_genus),
    pC3_core_clade_specific=len(c3_core_clade),
    genus_softcore_families=n_genus_core,
    MF6_pC3_CDS=len(mf6_c3_loci), MF6_chr1_CDS=len(mf6_chr1_loci),
    pC3_core_as_frac_of_MF6_pC3=round(len(c3_core) / len(mf6_c3_loci), 4),
    kegg_modules_parsed=len(MODULES),
    modules_complete_pC3_core=n_pc3_complete,
    modules_complete_chr1_core=n_chr1_complete,
    pC3_core_proteins_with_chr1_homologue=n_dup,
    pC3_core_proteins_tested=nq,
    coverage=cov_rows,
)
with open(f"{OUT}/setB_functional_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print("[10] done. Summary:", flush=True)
print(json.dumps(summary, indent=2), flush=True)
