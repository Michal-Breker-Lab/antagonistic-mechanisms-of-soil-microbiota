#!/usr/bin/env python3
"""Volcano plot of the MF6 co-culture contrast, with the ON/OFF band above it."""
import math
import sys
from types import SimpleNamespace
import zlib
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import svg_lib as S  # noqa: E402

FDR = float(snakemake.params.fdr)  # noqa: F821
LFC = float(snakemake.params.lfc)  # noqa: F821

COMPARISONS = [
    ("DE_coculture_d2", "co-culture vs monoculture, day 2",
     "on_off_withC_d2_vs_alone_d2", "withC_d2", "alone_d2", None),
]

PCT_DEFAULT = 50.0

TXSS_SYSTEM = {
    "T2SS": ("T2SS", "s3", "diamond"),
    "T3SS": ("T3SS", "s5", "triangle"),
    "T5aSS": ("T5SS", "s6", "plus"),
    "T5bSS": ("T5SS", "s6", "plus"),
    "T5cSS": ("T5SS", "s6", "plus"),
    "T6SSi": ("T6SS", "s1", "circle"),
    "T6SSii": ("T6SS", "s1", "circle"),
    "T6SSiii": ("T6SS", "s1", "circle"),
    "T4aP": ("T4aP", "s2", "triangle_down"),
    "T4bP": ("T4bP", "s4", "cross"),
    "Tad": ("Tad", "s8", "square"),
}
TXSS_ORDER = ["T6SSi", "T6SSii", "T6SSiii", "T2SS", "T3SS", "T5aSS", "T5bSS",
              "T5cSS", "T4aP", "T4bP", "Tad"]
TXSS_LEGEND_ORDER = ["T2SS", "T3SS", "T5SS", "T6SS", "T4aP", "T4bP", "Tad"]

CANDIDATE_GROUPS = [
    ("evidence_t6ss_prediction", "t6ss effector", "s1", "cross"),
    ("evidence_tox_pfam", "toxicity_domain", "s2", "cross"),
    ("evidence_extracellular", "extracellular", "s3", "cross"),
]

CAND_TXSS_SYSTEM = {
    "T2SS":    ("T2SS", "s6", "circle"),
    "T3SS":    ("T3SS", "s4", "circle"),
    "T5aSS":   ("T5SS", "s7", "circle"),
    "T5bSS":   ("T5SS", "s7", "circle"),
    "T5cSS":   ("T5SS", "s7", "circle"),
    "T6SSi":   ("T6SS", "s8", "circle"),
    "T6SSii":  ("T6SS", "s8", "circle"),
    "T6SSiii": ("T6SS", "s8", "circle"),
    "T4aP":    ("T4aP", "s2", "circle"),
    "T4bP":    ("T4bP", "s4", "circle"),
    "Tad":     ("Tad",  "s5", "circle"),
}
CAND_TXSS_LEGEND_ORDER = ["T2SS", "T3SS", "T5SS", "T6SS", "T4aP", "T4bP", "Tad"]

FINAL_SHAPE, FINAL_R = "star", 6.6

W, HB, HV = 940, 104, 452
PL, PR, PT = 66, 250, 16


def read_tsv(path: Path, comment="#") -> list[dict]:
    rows, hdr = [], None
    with path.open() as fh:
        for line in fh:
            if line.startswith(comment):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))
    return rows


def disp(g: str) -> str:
    """Group label as it should read on a figure: alone -> Chlamy-, withC -> Chlamy+."""
    return (g.replace("withC_", "Chlamy+ ").replace("alone_", "Chlamy- "))


def num(s, default=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return default



DROP_COLS = ("Protein Probability", "Combined Total Peptides")


def table_rows(rows_in, extra_cols, extra_fn):
    """(header, rows) for the highlighted proteins, extras after Description."""
    hl = [r for r in rows_in if r["_hl"]]
    cols = [c for c in rows_in[0] if not c.startswith("_") and c not in DROP_COLS] \
        if rows_in else []
    i = cols.index("Description") + 1 if "Description" in cols else len(cols)
    hdr = cols[:i] + list(extra_cols) + cols[i:]
    out = []
    for r in hl:
        e = extra_fn(r)
        out.append([str(r.get(c, "")) if c in r else str(e.get(c, "")) for c in hdr])
    return hdr, out


def write_table(path: Path, hdr, rows, header_note: str) -> None:
    with path.open("w") as fh:
        fh.write(header_note)
        fh.write("\t".join(hdr) + "\n")
        for row in rows:
            fh.write("\t".join(c.replace("\t", " ") for c in row) + "\n")


ROUND2 = ("logFC", "AveExpr", "t", "B", "mean_A", "mean_B",
          "mean_log2_LFQ_present", "pct_in_own_group", "pct_in_opposite_group")
SCI2 = ("P.Value", "adj.P.Val")



def build(de, onoff, title, sub, note, labA, labB, pct_cut,
          legend_title, legend_sub, legend_items, band_shapes=False):
    x0, x1 = PL, W - PR
    band_top, band_bot = PT + 26, PT + 26 + HB
    vol_top, vol_bot = band_bot + 34, band_bot + 34 + HV
    h = vol_bot + 62

    lfcs = [r["_lfc"] for r in de]
    m = max(abs(min(lfcs)), abs(max(lfcs))) * 1.06
    lo, hi = -m, m
    ymax = max(r["_nlp"] for r in de) * 1.08

    def sx(v):
        return x0 + (v - lo) / (hi - lo) * (x1 - x0)

    def sy(v):
        return vol_bot - v / ymax * (vol_bot - vol_top)

    p = S.head(W, h, title)

    p.append(S.rect(x0, band_top, x1 - x0, HB, cls="band", rx=3))
    p.append(S.text(x0, band_top - 8,
                    f"ON/OFF  —  x = percentile of the group it is missing "
                    f"from ({pct_cut:g} at centre → 100 at each edge)",
                    cls="axb"))
    xc = sx(0.0)
    p.append(S.line(xc, band_top, xc, band_bot, cls="base"))
    span = 100 - pct_cut
    for v in range(int(pct_cut), 101, 10):
        frac = (v - pct_cut) / span if span else 0
        for edge in ((x0, x1) if v > pct_cut else (x1,)):
            xx = xc + (edge - xc) * frac
            if v > pct_cut:
                p.append(S.line(xx, band_top, xx, band_bot, cls="grid"))
            p.append(S.text(xx, band_bot + 13, str(v), cls="ax", anchor="middle"))
    p.append(S.text(x0 + 6, band_top + 14, f"present in {disp(labB)}", cls="ax"))
    p.append(S.text(x1 - 6, band_top + 14, f"present in {disp(labA)}", cls="ax",
                    anchor="end"))

    yb = (band_top + band_bot) / 2 + 6
    shown = 0
    for r in sorted(onoff, key=lambda z: z["_pct"]):
        frac = (r["_pct"] - pct_cut) / span if span else 0
        frac = min(max(frac, 0), 1)
        right = r.get("present_in") == labA
        xx = xc + ((x1 if right else x0) - xc) * frac
        yy = yb + ((zlib.crc32(r["Protein"].encode()) % 5) - 2) * 7.0
        cls = r["_cls"] if r["_hl"] else "pt0"
        hover = (f"{r['Protein']}  present in {disp(r.get('present_in',''))}, "
                 f"absent from {disp(r.get('absent_from',''))}  "
                 f"pct {r['_pct']:.0f}  {r['_sig']}")
        if band_shapes:
            p.append(S.mark(r.get("_shape", "circle"), xx, yy,
                            5.2 if r["_hl"] else 2.6, cls=cls,
                            stroke="#fcfcfb" if r["_hl"] else None,
                            stroke_width=1.6 if r["_hl"] else None, title=hover))
        else:
            p.append(S.triangle(xx, yy, 5.2 if r["_hl"] else 3.4, cls=cls,
                                up=right, stroke="#fcfcfb" if r["_hl"] else None,
                                stroke_width=1.6 if r["_hl"] else None, title=hover))
        shown += 1
    nA = sum(1 for r in onoff if r.get("present_in") == labA)
    p.append(S.text(x1 - 6, band_bot - 8, f"{nA} proteins", cls="ax", anchor="end"))
    p.append(S.text(x0 + 6, band_bot - 8, f"{shown - nA} proteins", cls="ax"))
    if shown == 0:
        p.append(S.text((x0 + x1) / 2, yb, f"no ON/OFF proteins above "
                                           f"percentile {pct_cut:g}",
                        cls="ax", anchor="middle"))

    for t in S.nice_ticks(0, ymax, 5):
        p.append(S.line(x0, sy(t), x1, sy(t)))
        p.append(S.text(x0 - 8, sy(t) + 3.5, S.fmt(t), cls="ax", anchor="end"))
    for t in S.nice_ticks(lo, hi, 8):
        if lo <= t <= hi:
            p.append(S.text(sx(t), vol_bot + 16, S.fmt(t), cls="ax",
                            anchor="middle"))
    p.append(S.line(x0, vol_bot, x1, vol_bot, cls="base"))
    p.append(S.text((x0 + x1) / 2, vol_bot + 36, "log2 fold change", cls="axb",
                    anchor="middle"))
    p.append(S.text(16, (vol_top + vol_bot) / 2, "-log10 adjusted P", cls="axb",
                    anchor="middle", rotate=-90))
    p.append(S.line(x0, sy(-math.log10(FDR)), x1, sy(-math.log10(FDR)), cls="thr"))
    for v in (-LFC, LFC):
        if lo <= v <= hi:
            p.append(S.line(sx(v), vol_top, sx(v), vol_bot, cls="thr"))

    for r in de:
        if not r["_hl"]:
            p.append(S.circle(sx(r["_lfc"]), sy(r["_nlp"]), 1.7, cls="pt0",
                              opacity=0.5))
    labelled = []
    for r in de:
        if r["_hl"]:
            ring = r.get("_ring", 1.6)
            p.append(S.mark(r.get("_shape", "circle"),
                            sx(r["_lfc"]), sy(r["_nlp"]), r.get("_r", 4.6),
                            cls=r["_cls"],
                            stroke="#fcfcfb" if ring else None,
                            stroke_width=ring or None,
                            title=f"{r['Protein']}  logFC {r['_lfc']:+.2f}  "
                                  f"adjP {r['adj.P.Val']}  {r['_sig']}  "
                                  + (r["_note"] + "  " if r.get("_note") else "")
                                  + f"{r['Description'][:70]}"))
            if r.get("_label"):
                labelled.append(r)

    CH, LINE = 5.6, 15.0
    placed = []
    for r in sorted(labelled, key=lambda z: sy(z["_nlp"])):
        xx, yy = sx(r["_lfc"]), sy(r["_nlp"])
        w = len(r["_label"]) * CH
        right = xx + 12 + w <= x1
        ty = yy
        for q in placed:
            if abs(ty - q) < LINE:
                ty = q + LINE
        placed.append(ty)
        tx = xx + 12 if right else xx - 12
        p.append(S.text(tx, ty + 3.5, r["_label"], cls="lgd",
                        anchor="start" if right else "end"))
        if abs(ty - yy) > 4:
            p.append(S.line(xx + (7 if right else -7), yy, tx - (2 if right else -2),
                            ty, cls="base"))

    lx, ly = x1 + 26, band_top + 6
    p.append(S.text(lx, ly, legend_title, cls="axb"))
    p.append(S.text(lx, ly + 15, legend_sub, cls="ax"))
    for i, (lab, cls, shape) in enumerate(legend_items):
        yy = ly + 38 + i * 21
        p.append(S.mark(shape, lx + 6, yy - 4, 4.6, cls=cls, stroke="#fcfcfb",
                        stroke_width=1.6))
        p.append(S.text(lx + 18, yy, lab, cls="lgd"))
    yy = ly + 38 + len(legend_items) * 21 + 16
    if band_shapes:
        p.append(S.text(lx, yy, "ON/OFF band", cls="axb"))
        p.append(S.text(lx, yy + 15, f"right = ON in {disp(labA)}", cls="lgd"))
        p.append(S.text(lx, yy + 31, f"left = ON in {disp(labB)}", cls="lgd"))
        yy += 20
    else:
        p.append(S.triangle(lx + 6, yy - 4, 5.2, cls="pt4", up=True))
        p.append(S.text(lx + 18, yy, f"ON in {disp(labA)}", cls="lgd"))
        yy += 20
        p.append(S.triangle(lx + 6, yy - 4, 5.2, cls="pt4", up=False))
        p.append(S.text(lx + 18, yy, f"ON in {disp(labB)}", cls="lgd"))

    n_hl = sum(1 for r in de if r["_hl"])
    n_oo = len(onoff)
    n_oo_hl = sum(1 for r in onoff if r["_hl"])
    yy += 30
    p.append(S.text(lx, yy, f"{len(de):,} tested", cls="ax"))
    p.append(S.text(lx, yy + 15, f"{n_hl} highlighted", cls="ax"))
    p.append(S.text(lx, yy + 32, f"{n_oo} ON/OFF ({n_oo_hl} highlighted)", cls="ax"))
    if note:
        p.append(S.text(18, h - 14, "Note: " + note, cls="ax"))
    p.append("</svg>")
    return "\n".join(p)


def main() -> int:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        de_dir=Path(sm.params.de_dir),
        out=Path(sm.params.outdir),
        pct=float(sm.params.pct),
        no_pdf=False,
        highlight=sm.params.highlight,
        candidate_list=Path(sm.input.candidates),
        final_list=Path(sm.input.final),
        txss_table=Path(sm.input.txss),
    )
    args.out.mkdir(parents=True, exist_ok=True)


    if args.highlight.startswith("candidates"):
        cand = {r["MF6_ID"]: r
                for r in read_tsv(args.candidate_list) if r.get("MF6_ID")}
        if not cand:
            return print(f"no rows in {args.candidate_list}") or 2
        with_final = args.highlight == "candidates_txss_final"
        with_txss = args.highlight in ("candidates_txss", "candidates_txss_final")
        final = {}
        if with_final:
            rows = [l.rstrip("\n") for l in open(args.final_list) if l.strip()]
            hdr = rows[0].split("\t")
            for line in rows[1:]:
                r = dict(zip(hdr, line.split("\t")))
                final[r["MF6_ID"]] = (r.get("Description", "").strip(),
                                      (r.get("note") or "").strip())
            print(f"shortlist: {len(final)} proteins in {args.final_list.name}")
        txss: dict[str, list[dict]] = {}
        if with_txss:
            for t in read_tsv(args.txss_table):
                if t.get("MF6_ID") and t.get("system") in CAND_TXSS_SYSTEM:
                    txss.setdefault(t["MF6_ID"], []).append(t)
        legend_title = "Candidate list" + (" + TXSScan" if with_txss else "")
        legend_sub = ("strongest evidence, then secretion system" if with_txss
                      else "coloured by strongest evidence")
        legend_items = []
        tag = ("candidates_plus_txsscan_final" if with_final else
               "candidates_plus_txsscan" if with_txss else "candidates")
        extra_cols = ["evidence_group", "n_evidence", "txss_system", "manual_worth",
                      "final_note", "DeepLoc", "SignalP", "T6SS_effector_probability"]
        keep_cols = ["Protein", "Gene", "evidence_group", "manual_worth", "txss_system",
                     "logFC", "adj.P.Val", "pct_in_opposite_group", "final_note",
                     "Description"]

        def group_of(c):
            """First matching group in CANDIDATE_GROUPS order = strongest claim."""
            for col, lab, cls, shape in CANDIDATE_GROUPS:
                if (c.get(col) or "").strip().lower() in ("yes", "y", "true", "1"):
                    return lab, cls, shape
            return None

        def sys_rank(h):
            s = h.get("system", "")
            return TXSS_ORDER.index(s) if s in TXSS_ORDER else len(TXSS_ORDER)

        def decorate(r):
            c = cand.get(r["Protein"])
            g = group_of(c) if c else None
            hits = sorted(txss.get(r["Protein"], []), key=sys_rank)
            if g is None and hits:
                g = CAND_TXSS_SYSTEM[hits[0]["system"]]
            r["_hl"] = g is not None
            r["_sys"] = " / ".join(dict.fromkeys(h["system"] for h in hits))
            if not g:
                r["_cls"], r["_shape"], r["_sig"], r["_grp"] = "pt0", "circle", "", None
                return r
            r["_grp"], r["_cls"], r["_shape"] = g
            all_g = [lab for col, lab, _c, _s in CANDIDATE_GROUPS
                     if c and (c.get(col) or "").strip().lower() == "yes"]
            if hits:
                all_g.append(f"TXSScan {r['_sys']}")
            r["_sig"] = " + ".join(all_g)
            if r["Protein"] in final:
                desc, note = final[r["Protein"]]
                r["_shape"], r["_r"] = FINAL_SHAPE, FINAL_R
                r["_label"], r["_note"] = desc or r["Description"][:40], note
            return r

        def extra_fn(r):
            c = cand.get(r["Protein"], {})
            return {"evidence_group": r.get("_grp") or "", "n_evidence": r["_sig"],
                    "txss_system": r.get("_sys", ""),
                    "manual_worth": "Yes" if r["Protein"] in final else "",
                    "final_note": r.get("_note", ""),
                    "DeepLoc": c.get("DeepLoc", ""), "SignalP": c.get("SignalP", ""),
                    "T6SS_effector_probability": c.get("T6SS_effector_probability", "")}

        def legend_from(rows):
            """Candidate groups first, then only the systems actually present, with
            their counts.  Colour and mark come from the fixed maps, never from position
            in this list, so a system dropping out does not repaint the survivors."""
            n = Counter(r["_grp"] for r in rows if r["_hl"])
            out = [(f"{lab}   {n[lab]}", cls, shape)
                   for _col, lab, cls, shape in CANDIDATE_GROUPS if n.get(lab)]
            seen = {lab: (cls, shape) for lab, cls, shape in CAND_TXSS_SYSTEM.values()}
            out += [(f"{lab}   {n[lab]}", *seen[lab])
                    for lab in CAND_TXSS_LEGEND_ORDER if n.get(lab)]
            n_star = sum(1 for r in rows if r.get("_label"))
            if n_star:
                out.append((f"shortlist, named   {n_star}", "pt4", FINAL_SHAPE))
            return out
    elif args.highlight == "txss":
        txss: dict[str, list[dict]] = {}
        for t in read_tsv(args.txss_table):
            if t.get("MF6_ID"):
                txss.setdefault(t["MF6_ID"], []).append(t)
        legend_title = "Secretion systems (TXSScan)"
        legend_sub = "one colour + one mark per system"
        legend_items = []
        tag = "txss"
        extra_cols = ["system", "components"]
        keep_cols = ["Protein", "Gene", "system", "components", "logFC",
                     "adj.P.Val", "pct_in_opposite_group", "evidence",
                     "Description"]

        def extra_fn(r):
            return {"system": r.get("_grp") or "", "components": r["_sig"]}

        def rank(h):
            s = h.get("system", "")
            return TXSS_ORDER.index(s) if s in TXSS_ORDER else len(TXSS_ORDER)

        def decorate(r):
            hits = sorted((h for h in txss.get(r["Protein"], [])
                           if h["system"] in TXSS_SYSTEM), key=rank)
            r["_hl"] = bool(hits)
            if not hits:
                r["_cls"], r["_shape"], r["_sig"], r["_grp"] = \
                    "pt0", "circle", "", None
                return r
            r["_grp"], r["_cls"], r["_shape"] = TXSS_SYSTEM[hits[0]["system"]]
            r["_sig"] = " / ".join(
                f"{h['system']} {h['component'].split('_', 1)[-1]} "
                f"({h['hit_status']})" for h in hits)
            return r

        def legend_from(rows):
            """Only the systems actually present, with their counts.

            Colour and mark come from TXSS_SYSTEM, never from position in this
            list, so a system that drops out does not repaint the survivors.
            """
            n = Counter(r["_grp"] for r in rows if r["_hl"])
            seen = {}
            for s, (lab, cls, shape) in TXSS_SYSTEM.items():
                seen[lab] = (cls, shape)
            return [(f"{lab}   {n[lab]}", *seen[lab])
                    for lab in TXSS_LEGEND_ORDER if n.get(lab)]
    else:
        return print(f"unknown --highlight: {args.highlight}") or 2

    made = []
    for de_name, sub, oo_name, labA, labB, note in COMPARISONS:
        de = read_tsv(args.de_dir / f"{de_name}.tsv")
        for r in de:
            r["_lfc"] = num(r["logFC"], 0.0)
            p = num(r["adj.P.Val"])
            r["_nlp"] = -math.log10(max(p, 1e-300)) if p and p > 0 else 0.0
            decorate(r)
        oo_path = args.de_dir / f"{oo_name}.tsv"
        oo = read_tsv(oo_path) if oo_path.exists() else []
        oo = [decorate(r) for r in oo
              if (num(r.get("pct_in_opposite_group")) or -1) > args.pct]
        for r in oo:
            r["_pct"] = num(r["pct_in_opposite_group"], 0.0)

        svg = build(de, oo, f"MF6 — {sub}",
                    f"volcano: {de_name}   •   ON/OFF: {oo_name} "
                    f"(percentile > {args.pct:g})", note, labA, labB, args.pct,
                    legend_title, legend_sub,
                    legend_from(de + oo) if legend_from else legend_items,
                    band_shapes=args.highlight in ("candidates", "candidates_txss"))
        stem = f"volcano_{tag}_{de_name.replace('DE_','')}"
        if args.pct != PCT_DEFAULT:
            stem += f"_pct{args.pct:g}"
        f = args.out / f"{stem}.svg"
        f.write_text(svg)

        note = (f"# highlighted proteins in {stem}.pdf\n"
                f"# {legend_title} — {legend_sub}\n")
        hdr_v, rows_v = table_rows(de, extra_cols, extra_fn)
        write_table(args.out / f"{stem}_highlighted.tsv", hdr_v, rows_v,
                    note + f"# panel: volcano ({de_name}.tsv), "
                           f"logFC = {disp(labA)} − {disp(labB)}\n"
                           f"# {len(rows_v)} of {len(de)} tested proteins, "
                           f"sorted by adj.P.Val\n")
        hdr_o, rows_o = (table_rows(oo, extra_cols, extra_fn) if oo else ([], []))
        if rows_o:
            write_table(args.out / f"{stem}_onoff_highlighted.tsv", hdr_o, rows_o,
                        note + f"# panel: ON/OFF band ({oo_name}.tsv), "
                               f"percentile > {args.pct:g}\n"
                               f"# {len(rows_o)} of {len(oo)} shown proteins, "
                               f"sorted by pct_in_opposite_group\n")
        made.append((f, hdr_v, rows_v, hdr_o, rows_o))
        print(f"  {f.name:46s} {len(de):5d} tested, "
              f"{len(rows_v):3d} highlighted, "
              f"{len(oo):3d} on/off shown ({len(rows_o)} highlighted)")


    if not args.no_pdf:
        write_pdfs([f for f, *_ in made], args.out)
    return 0


def write_pdfs(svgs, out: Path) -> None:
    """Vector PDF per SVG, with selectable, editable text - for Illustrator.

    Rendering lives in svg_lib.svg_to_pdf, which uses librsvg from this rule's
    conda env.  It replaced Chrome print-to-pdf: same real-text output, but a
    declared dependency instead of whatever browser happened to be on PATH.
    """
    for svg_path in svgs:
        pdf = svg_path.with_suffix(".pdf")
        S.svg_to_pdf(svg_path, pdf)
        print(f"  {pdf.name}")


if __name__ == "__main__":
    sys.exit(main())
