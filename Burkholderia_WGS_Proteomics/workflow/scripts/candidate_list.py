#!/usr/bin/env python3
"""Genome-wide toxin candidates, scored by extracellular localisation, tox-labelled
Pfam domains and T6SS effector prediction.
"""
import csv, re
from collections import defaultdict
from pathlib import Path

import sys

sm = snakemake  # noqa: F821
sys.stderr = sys.stdout = open(sm.log[0], "w")

ROOT = Path(sm.params.root)
ABSENT_IN_MUTANTS = sm.params.absent
LOCUS_PREFIX = sm.params.locus_prefix
LOCUS_RE = re.compile(re.escape(LOCUS_PREFIX) + r'\d+')
OUT = Path(sm.output.table).parent
MIN_REPS = int(sm.params.min_reps)
ONOFF_PCT_GATE = float(sm.params.pct_gate)


def rd(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith('#')), delimiter='\t'))


def pfacc(v):
    m = re.search(r'(\d{4,6})', v or '')
    return f'PF{m.group(1)}' if m else None


coords, product = {}, {}
for line in open(sm.input.gff):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'CDS':
        continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    if 'locus_tag' in a:
        coords[a['locus_tag']] = (f[0], int(f[3]), int(f[4]), f[6])
        product[a['locus_tag']] = a.get('product', '').replace('%2C', ',')

loc = {}
for r in rd(sm.input.topology):
    p = re.sub(r'^gnl\|extdb\|', '', r['protein_id'])
    loc[p] = r['deeplocpro_localization']

pred = {}
for r in rd(sm.input.t6ss):
    if r['protein_id']:
        pred[r['protein_id']] = (r['gene_id'], r['IMG_locus_tag'], float(r['xg_probas']))

cls, pname = {}, {}
for r in csv.DictReader(open(sm.input.pfams)):
    k = pfacc(r['pfamID'])
    if k:
        cls[k] = (r['pfam_toxity'] or '').strip().lower() or 'nan'
        pname[k] = r['pfam_real_name']

pfam = defaultdict(dict)
for line in open(sm.input.interpro):
    f = line.rstrip('\n').split('\t')
    if len(f) < 6 or f[3] != 'Pfam':
        continue
    m = LOCUS_RE.search(f[0])
    if m:
        pfam[m.group(0)][f[4].split('.')[0]] = f[5]

def de(p):
    out = {}
    for r in rd(p):
        try:
            out[r['Protein']] = (float(r['logFC']), float(r['adj.P.Val']))
        except (ValueError, KeyError):
            pass
    return out

d2 = de(ROOT / 'DE/MF6/DE_coculture_d2.tsv')
mut   = {m: de(ROOT / f'DE/{m}/DE_{m}_withC_vs_MF6_withC.tsv') for m in ('27D6', '34F7')}
alone = {m: de(ROOT / f'DE/{m}/DE_{m}_alone_vs_MF6_alone.tsv') for m in ('27D6', '34F7')}

def onoff(path, ga, gb):
    out = {}
    if not Path(path).exists():
        return out
    for r in rd(path):
        pres = r['present_in']
        full = pres if '_' in pres and not pres.startswith(('withC', 'alone')) else f'MF6_{pres}'
        state = 'ON' if full == ga else 'OFF' if full == gb else None
        if state and float(r['pct_in_opposite_group']) > ONOFF_PCT_GATE:
            out[r['Protein']] = (f"{state}, {r['evidence'].split()[0]} "
                                 f"(log2 LFQ {float(r['mean_log2_LFQ_present']):.2f})")
    return out

OO = {
    ('MF6_withC_d2', 'MF6_alone_d2'):
        onoff(ROOT / 'DE/MF6/on_off_withC_d2_vs_alone_d2.tsv',
              'MF6_withC_d2', 'MF6_alone_d2'),
}
for m in ('27D6', '34F7'):
    OO[(f'{m}_withC_d2', 'MF6_withC_d2')] = onoff(
        ROOT / f'DE/{m}/on_off_{m}_withC_vs_MF6_withC.tsv',
        f'{m}_withC_d2', 'MF6_withC_d2')
    OO[(f'{m}_alone_d2', 'MF6_alone_d2')] = onoff(
        ROOT / f'DE/{m}/on_off_{m}_alone_vs_MF6_alone.tsv',
        f'{m}_alone_d2', 'MF6_alone_d2')
cog = {}
for line in open(sm.input.eggnog):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    m = LOCUS_RE.search(f[0]) if f else None
    if m and len(f) > 6 and f[6] not in ('', '-'):
        cog[m.group(0)] = f[6]
_lfq = {r['Protein']: r for r in rd(ROOT / 'Proteomics/MF6_maxlfq_intensity.tsv')}
QC_DROP = {'34F7_alone_d2_1'}
_cols = [c for c in next(iter(_lfq.values()))
         if re.search(r'_d2_[1-4]$', c) and c not in QC_DROP]
GROUP = {}
for c in _cols:
    GROUP.setdefault(c.rsplit('_', 1)[0], []).append(c)

def nvalid(p, g):
    r = _lfq.get(p)
    cols = GROUP[g]
    if r is None:
        return 0, len(cols)
    return sum(1 for c in cols if float(r[c] or 0) > 0), len(cols)

def expressed(p):
    """'+' when the protein is quantified in >= MIN_REPS replicates of at least
    one day-2 group.  Any group counts - a protein seen only in a mutant, or
    only in monoculture, is still expressed somewhere.  34F7_alone_d2 has 3
    runs after QC, so there '>= 3' means all three."""
    return '+' if any(nvalid(p, g)[0] >= MIN_REPS for g in GROUP) else '-'

def star(p, tbl, ga, gb, mutant):
    """One contrast cell.  `NA (contig 2)` belongs ONLY to mutant comparisons -
    both mutants have lost that replicon, so a mutant number there would be a
    deletion reported as regulation.  contig_2 is present in MF6, so the
    wild-type co-culture contrast is measurable and an absent value there is an
    ordinary ND.

    The contig_2 guard has to come FIRST, before the DE lookup.  MaxLFQ's
    match-between-runs transfers intensities from the wild-type runs onto genes
    the mutants do not carry, so a deleted gene can still appear in a mutant DE
    table with a fold change and a p-value - e.g. MF6_003684 at -0.174 ns,
    reading as "unchanged in the mutant" for a gene that is not in its genome.
    Those numbers are artefacts and must never be reported."""
    if mutant and coords[p][0] == ABSENT_IN_MUTANTS:
        return f'NA ({ABSENT_IN_MUTANTS})'
    if p in tbl:
        v, q = tbl[p]
        return f'{v:.3f}' + ('*' if q < 0.05 else '')
    cell = OO.get((ga, gb), {}).get(p)
    if cell:
        return cell
    a, na = nvalid(p, ga)
    b, nb = nvalid(p, gb)
    return f'ND ({a}/{na} vs {b}/{nb})'

HDR = ['MF6_ID', 'Description', 'chr', 'start', 'stop', 'strand',
       'n_evidence', 'evidence_extracellular', 'evidence_tox_pfam',
       'tox_pfams', 'tox_pfam_names',
       'evidence_t6ss_prediction',
       'DeepLoc', 'expressed',
       'T6SS_effector_probability',
       'MF6C_vs_MF6_d2_LFC',
       '27D6C_vs_MF6C_LFC', '27D6_vs_MF6_LFC',
       '34F7C_vs_MF6C_LFC', '34F7_vs_MF6_LFC',
       ]
rows = []
for p in sorted(coords):
    hit = {k: v for k, v in pfam.get(p, {}).items() if k in cls}
    ks = sorted(hit)
    is_extra = loc.get(p) == 'Extracellular'
    tox_ks = [k for k in ks if cls[k] == 'tox']
    is_tox = bool(tox_ks)
    is_imm = any(cls[k] == 'immunity' for k in ks)
    is_pred = p in pred
    n_ev = sum((is_extra, is_tox, is_pred))
    if not n_ev:
        continue
    proba = pred[p][2] if p in pred else None
    c, s, e, st = coords[p]
    rows.append(dict(zip(HDR, [
        p, product.get(p, ''), c, s, e, st,
        n_ev, 'Yes' if is_extra else 'No', 'Yes' if is_tox else 'No',
        ','.join(tox_ks), ','.join(pname.get(k, '') for k in tox_ks),
        'Yes' if is_pred else 'No',
        loc.get(p, ''), expressed(p),
        f'{proba:.3f}' if proba is not None else '',
        star(p, d2, 'MF6_withC_d2', 'MF6_alone_d2', mutant=False),
        star(p, mut['27D6'], '27D6_withC_d2', 'MF6_withC_d2', mutant=True),
        star(p, alone['27D6'], '27D6_alone_d2', 'MF6_alone_d2', mutant=True),
        star(p, mut['34F7'], '34F7_withC_d2', 'MF6_withC_d2', mutant=True),
        star(p, alone['34F7'], '34F7_alone_d2', 'MF6_alone_d2', mutant=True),
    ])))
def _lfc(v):
    try:
        return float(v.rstrip('*'))
    except ValueError:
        return -99.0
rows.sort(key=lambda r: (-int(r['n_evidence']),
                         -float(r['T6SS_effector_probability'] or 0),
                         -_lfc(r['MF6C_vs_MF6_d2_LFC'])))

with open(OUT / 'candidate_list.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, HDR, delimiter='\t'); w.writeheader(); w.writerows(rows)
print(f'extracellular {sum(1 for v in loc.values() if v=="Extracellular")} | '
      f'T6SS predictions mapped {len(pred)} | proteins with an all_pfams domain '
      f'{sum(1 for p in pfam if set(pfam[p]) & set(cls))}')
from collections import Counter
print(f'-> toxins/candidate_list.tsv: {len(rows)} genes (extracellular OR tox-Pfam '
      f'OR T6SS-predicted, whole genome)')
print('   by evidence count:', dict(sorted(Counter(r['n_evidence'] for r in rows).items())))
for k, lbl in [('evidence_extracellular', 'extracellular'),
               ('evidence_tox_pfam', 'tox-labelled Pfam'),
               ('evidence_t6ss_prediction', 'T6SS prediction')]:
    print(f'   {lbl:20} {sum(1 for r in rows if r[k] == "Yes")}')
print(f'   {sum(1 for r in rows if r["evidence_tox_pfam"] == "Yes")} carry a '
      f'tox-labelled Pfam (kept as evidence_tox_pfam)')
