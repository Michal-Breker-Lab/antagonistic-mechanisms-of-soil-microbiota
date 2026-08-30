#!/usr/bin/env python3
"""Volcano plots of the three contrasts of one mutant."""
import importlib.util
import math
import sys
from types import SimpleNamespace
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import svg_lib as S  # noqa: E402

FDR = float(snakemake.params.fdr)  # noqa: F821  config figures.volcano
LFC = float(snakemake.params.lfc)  # noqa: F821


def comparisons(strain: str):
    """(DE table stem, subtitle, group A label, group B label).

    logFC is A - B in every table, so A is always the right-hand side.
    """
    return [
        (f"DE_{strain}_withC_vs_alone",
         f"{strain} co-culture vs monoculture, day 2",
         f"{strain} Chlamy+", f"{strain} Chlamy−"),
        (f"DE_{strain}_withC_vs_MF6_withC",
         f"{strain} vs MF6, both in co-culture, day 2",
         f"{strain} Chlamy+", "MF6 Chlamy+"),
        (f"DE_{strain}_alone_vs_MF6_alone",
         f"{strain} vs MF6, both in monoculture, day 2",
         f"{strain} Chlamy−", "MF6 Chlamy−"),
    ]


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

_spec = importlib.util.spec_from_file_location(
    "volcano_onoff", str(HERE / "volcano_mf6.py"))
_V6 = importlib.util.module_from_spec(_spec)
_V6.snakemake = snakemake  # noqa: F821
_spec.loader.exec_module(_V6)
CANDIDATE_GROUPS = _V6.CANDIDATE_GROUPS
CAND_TXSS_SYSTEM = _V6.CAND_TXSS_SYSTEM
CAND_TXSS_LEGEND_ORDER = _V6.CAND_TXSS_LEGEND_ORDER
FINAL_SHAPE, FINAL_R = _V6.FINAL_SHAPE, _V6.FINAL_R

W, HV = 940, 520
PL, PR, PT = 66, 250, 34


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


def num(s, default=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return default



DROP_COLS = ("Protein Probability", "Combined Total Peptides")


def table_rows(de, extra_cols, extra_fn):
    """(header, rows) for the highlighted proteins, extras after Description."""
    hl = [r for r in de if r["_hl"]]
    de_cols = [c for c in de[0] if not c.startswith("_") and c not in DROP_COLS] \
        if de else []
    i = de_cols.index("Description") + 1 if "Description" in de_cols else len(de_cols)
    hdr = de_cols[:i] + list(extra_cols) + de_cols[i:]
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



def build(de, title, labA, labB, legend_title, legend_sub, legend_items):
    x0, x1 = PL, W - PR
    vol_top, vol_bot = PT + 10, PT + 10 + HV
    h = vol_bot + 74

    lfcs = [r["_lfc"] for r in de]
    m = max(abs(min(lfcs)), abs(max(lfcs))) * 1.06
    lo, hi = -m, m
    ymax = max(r["_nlp"] for r in de) * 1.08

    def sx(v):
        return x0 + (v - lo) / (hi - lo) * (x1 - x0)

    def sy(v):
        return vol_bot - v / ymax * (vol_bot - vol_top)

    p = S.head(W, h, title)

    for t in S.nice_ticks(0, ymax, 6):
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
    p.append(S.text(x0, vol_bot + 56, f"← higher in {labB}", cls="ax"))
    p.append(S.text(x1, vol_bot + 56, f"higher in {labA} →", cls="ax",
                    anchor="end"))
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
        right = xx + 12 + len(r["_label"]) * CH <= x1
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

    lx, ly = x1 + 26, vol_top + 12
    p.append(S.text(lx, ly, legend_title, cls="axb"))
    p.append(S.text(lx, ly + 15, legend_sub, cls="ax"))
    for i, (lab, cls, shape) in enumerate(legend_items):
        yy = ly + 38 + i * 21
        p.append(S.mark(shape, lx + 6, yy - 4, 4.6, cls=cls, stroke="#fcfcfb",
                        stroke_width=1.6))
        p.append(S.text(lx + 18, yy, lab, cls="lgd"))

    n_hl = sum(1 for r in de if r["_hl"])
    n_sig = sum(1 for r in de if r["_nlp"] > -math.log10(FDR)
                and abs(r["_lfc"]) >= LFC)
    yy = ly + 38 + len(legend_items) * 21 + 24
    p.append(S.text(lx, yy, f"{len(de):,} tested", cls="ax"))
    p.append(S.text(lx, yy + 15, f"{n_hl} highlighted", cls="ax"))
    p.append(S.text(lx, yy + 30, f"{n_sig} FDR<{FDR:g} & |log2FC|≥1", cls="ax"))
    p.append("</svg>")
    return "\n".join(p)


def main() -> int:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        strain=[sm.wildcards.strain],
        de_root=Path(sm.params.de_root),
        out=Path(sm.params.outdir),
        no_pdf=False,
        highlight=sm.params.highlight,
        candidate_list=Path(sm.input.candidates),
        final_list=Path(sm.input.final),
        txss_table=Path(sm.input.txss),
    )


    if args.highlight == "txss":
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
                     "adj.P.Val", "Description"]

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
            """Only the systems actually present, with their counts.  Colour and
            mark come from TXSS_SYSTEM, never from position in this list."""
            n = Counter(r["_grp"] for r in rows if r["_hl"])
            seen = {}
            for s, (lab, cls, shape) in TXSS_SYSTEM.items():
                seen[lab] = (cls, shape)
            return [(f"{lab}   {n[lab]}", *seen[lab])
                    for lab in TXSS_LEGEND_ORDER if n.get(lab)]
    elif args.highlight.startswith("candidates"):
        cand = {r["MF6_ID"]: r
                for r in read_tsv(args.candidate_list) if r.get("MF6_ID")}
        if not cand:
            return print(f"no rows in {args.candidate_list}") or 2
        with_final = args.highlight == "candidates_txss_final"
        with_txss = args.highlight in ("candidates_txss", "candidates_txss_final")
        final = {}
        if with_final:
            rws = [l.rstrip("\n") for l in open(args.final_list) if l.strip()]
            hdr = rws[0].split("\t")
            for line in rws[1:]:
                q = dict(zip(hdr, line.split("\t")))
                final[q["MF6_ID"]] = (q.get("Description", "").strip(),
                                      (q.get("note") or "").strip())
            print(f"shortlist: {len(final)} proteins")
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
                     "logFC", "adj.P.Val", "final_note", "Description"]

        def group_of(c):
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
    else:
        return print(f"unknown --highlight: {args.highlight}") or 2

    for strain in args.strain:
        de_dir = args.de_root / strain
        out = args.out or (de_dir / "plots")
        out.mkdir(parents=True, exist_ok=True)
        print(f"{strain}  ->  {out}")

        made = []
        for de_name, sub, labA, labB in comparisons(strain):
            de = read_tsv(de_dir / f"{de_name}.tsv")
            for r in de:
                r["_lfc"] = num(r["logFC"], 0.0)
                pv = num(r["adj.P.Val"])
                r["_nlp"] = -math.log10(max(pv, 1e-300)) if pv and pv > 0 else 0.0
                decorate(r)

            svg = build(de, f"{strain} — {sub}", labA, labB,
                        legend_title, legend_sub,
                        legend_from(de) if legend_from else legend_items)
            stem = f"volcano_{tag}_{de_name.replace('DE_', '')}"
            f = out / f"{stem}.svg"
            f.write_text(svg)

            hdr, rows = table_rows(de, extra_cols, extra_fn)
            tbl = out / f"{stem}_highlighted.tsv"
            write_table(tbl, hdr, rows,
                        f"# highlighted proteins in {stem}.pdf\n"
                        f"# {legend_title} — {legend_sub}\n"
                        f"# contrast: {labA} vs {labB}; logFC = {labA} − {labB}\n"
                        f"# {len(rows)} of {len(de)} tested proteins, "
                        f"sorted by adj.P.Val as in {de_name}.tsv\n")
            made.append((f, tbl, hdr, rows))
            print(f"  {f.name:52s} {len(de):5d} tested, "
                  f"{len(rows):3d} highlighted -> {tbl.name}")


        if not args.no_pdf:
            write_pdfs([f for f, *_ in made], out)
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
