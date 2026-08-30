#!/usr/bin/env python3
"""Bubble plot of the COG and KEGG enrichment for the MF6 co-culture contrast."""
from types import SimpleNamespace
import csv
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svg_lib as S                                            # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_IN = ROOT / "DE/MF6/enrichment"

VARIANT = "de_plus_onoff"
VARIANT_SUFFIX = {"de_only": "_no_onoff",
                  "de_plus_onoff_all": "_onoff_all",
                  "de_plus_onoff_p50": "_onoff_p50"}
ALPHA = float(snakemake.params.alpha)  # noqa: F821
FDR_CAP = float(snakemake.params.fdr_cap)  # noqa: F821  -log10 ramp saturation

W, PAD_L, PAD_R, PAD_T = 740, 322, 26, 48
ROW_H, SUB_HDR, BLK_HDR, PAD_B = 21, 24, 17, 100
R_MAX, R_MIN = 11.0, 3.2
LABEL_CHARS = 50

RAMP = S.SEQ[3:]
BAR_W, BAR_H = 176, 10
BAR_TICKS = [1.5, 2, 2.5, 3, 3.5, 4]

CONTRASTS = [
    ("coculture_d2", "Up in co-culture", "Down in co-culture"),
]

STYLE = """
  .surface { fill: #ffffff; }
  text { font-family: %(font)s; fill: #000000; }
  .sub  { font-size: 11px; }
  .ax   { font-size: 10px; }
  .axb  { font-size: 11px; }
  .hdr  { font-size: 11.5px; font-weight: 600; }
  .val  { font-size: 10.5px; }
  .grid { stroke: #e1e0d9; stroke-width: 1; }
  .base { stroke: #c3c2b7; stroke-width: 1; }
  .ring { stroke: #ffffff; stroke-width: 1.5; }
  .thr  { stroke: #898781; stroke-width: 1; stroke-dasharray: 4 3; }
""" % {"font": S.FONT}


def read_rows(indir: Path, variant: str = VARIANT):
    rows = []
    for ont in ("COG", "KEGG"):
        f = indir / f"{ont}_enrichment.tsv"
        if not f.is_file():
            sys.exit(f"not found: {f}  (run 06_enrichment_MF6.R first)")
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["variant"] != variant or float(r["p_adj"]) >= ALPHA:
                    continue
                r["p_adj"] = float(r["p_adj"])
                r["mlog"] = -math.log10(r["p_adj"])
                r["fold"] = float(r["fold_enrichment"])
                r["k"], r["K"] = int(r["k_fg"]), int(r["K_bg"])
                r["n"], r["N"] = int(r["n_fg"]), int(r["N_bg"])
                rows.append(r)
    return rows


def trunc(s: str, n: int) -> str:
    return s if len(s) <= n else s[: n - 1].rstrip(" ,") + "…"


def is_junk_map(term: str) -> bool:
    """KEGG global/overview and human-disease maps carry generic bacterial genes
    (LPS, secretion, central metabolism), so they light up on any bacterial hit
    list.  Flagged rather than dropped, so the call stays visible."""
    return term.startswith(("map05", "map011", "map012"))


def fmt_p(p: float) -> str:
    if p < 1e-4:
        e = int(math.floor(math.log10(p)))
        return f"{p / 10 ** e:.1f}e{e}"
    return f"{p:.3f}".lstrip("0")


def ramp_t(mlog: float) -> float:
    lo = -math.log10(ALPHA)
    return max(0.0, min(1.0, (min(mlog, FDR_CAP) - lo) / (FDR_CAP - lo)))


def ramp_color(mlog: float) -> str:
    """Continuous colour, interpolated between the ramp stops in sRGB - the same
    space and the same stops the <linearGradient> in the legend uses, so a bubble
    and the point of the colourbar it sits under are the identical colour."""
    t = ramp_t(mlog) * (len(RAMP) - 1)
    i = min(int(t), len(RAMP) - 2)
    f = t - i
    a = [int(RAMP[i][j:j + 2], 16) for j in (1, 3, 5)]
    b = [int(RAMP[i + 1][j:j + 2], 16) for j in (1, 3, 5)]
    return "#" + "".join(f"{round(a[j] + (b[j] - a[j]) * f):02x}" for j in range(3))


def head(w: int, h: int, label: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
        f'width="{w}" height="{h}" role="img" aria-label="{S.esc(label)}">',
        f"<style>{STYLE}</style>",
        f'<rect class="surface" x="0" y="0" width="{w}" height="{h}"/>',
    ]


def bubble_plot(rows, up_lab, down_lab, fold_max, k_max, uid="fdr"):
    """One figure: an `up` panel above a `down` panel, sharing one x axis."""
    panels = []
    for direction, lab in (("up", up_lab), ("down", down_lab)):
        blocks = []
        for ont, hdr in (("COG", "COG category"), ("KEGG", "KEGG pathway")):
            sel = sorted([r for r in rows
                          if r["direction"] == direction and r["ontology"] == ont],
                         key=lambda r: -r["fold"])
            if sel:
                blocks.append((hdr, sel))
        if blocks:
            panels.append((lab, blocks))

    body = sum(SUB_HDR + sum(BLK_HDR + len(s) * ROW_H + 6 for _, s in bl)
               for _, bl in panels)
    foot = any(is_junk_map(r["term"]) for r in rows)
    pad_b = PAD_B + (16 if foot else 0)
    h = int(PAD_T + body + 10 + pad_b)
    x0, x1 = PAD_L, W - PAD_R

    plot_w = x1 - x0
    margin = R_MAX + 8
    lo = 1.0
    hi = lo + max(fold_max - lo, 0.2) * plot_w / max(plot_w - margin, 1.0)

    def xs(v):
        return x0 + (v - lo) / (hi - lo) * plot_w

    def radius(k):
        return round(R_MIN + (R_MAX - R_MIN) * math.sqrt(k / k_max), 2)

    p = head(W, h, "COG and KEGG enrichment")
    p.append(S.text(14, 16, f"Fisher exact, BH < {ALPHA:g}. Bubble area = proteins "
                            "in term; fill = −log10 FDR.", cls="sub"))

    plot_bot = h - pad_b + 4
    ticks = [t for t in S.nice_ticks(lo, hi, 7) if t > 1.4 and t <= hi + 1e-9]
    ticks = [1.0] + ticks
    for t in ticks:
        if t > 1:
            p.append(S.line(xs(t), PAD_T - 10, xs(t), plot_bot))
        p.append(S.text(xs(t), PAD_T - 16, f"{t:g}×", cls="ax", anchor="middle"))
    p.append(S.line(xs(1), PAD_T - 10, xs(1), plot_bot, cls="thr"))

    y = PAD_T
    for lab, blocks in panels:
        p.append(S.text(14, y + 6, lab, cls="hdr"))
        y += SUB_HDR
        for hdr, sel in blocks:
            p.append(S.text(26, y + 4, hdr, cls="axb"))
            y += BLK_HDR
            for r in sel:
                cy = y + ROW_H / 2
                name = r["term_name"] or r["term"]
                label = f"[{r['term']}] {name}"
                label = trunc(label, LABEL_CHARS) + (" †" if is_junk_map(r["term"]) else "")
                p.append(S.text(x0 - 12, cy + 3.6, label, cls="val", anchor="end"))
                p.append(S.line(x0 - 7, cy, xs(r["fold"]) - radius(r["k"]) - 3, cy,
                                cls="grid"))
                tip = (f"{r['term']}  {r['term_name']}\n"
                       f"{r['k']} of {r['n']} {lab} proteins carry this term\n"
                       f"{r['K']} of {r['N']} in the tested universe\n"
                       f"fold enrichment {r['fold']:g}×   FDR {fmt_p(r['p_adj'])}")
                p.append(S.circle(xs(r["fold"]), cy, radius(r["k"]),
                                  fill=ramp_color(r["mlog"]), cls=None,
                                  stroke="#ffffff", stroke_width=1.5, title=tip))
                y += ROW_H
            y += 6

    p.append(S.line(x0, plot_bot, x1, plot_bot, cls="base"))
    p.append(S.text((x0 + x1) / 2, plot_bot + 24, "fold enrichment",
                    cls="axb", anchor="middle"))

    size_w = 46 + 3 * (2 * R_MAX + 20) - 20
    total = BAR_W + 54 + size_w
    lx = (W - total) / 2
    ly = plot_bot + 48

    p.append(f'<defs><linearGradient id="{uid}" x1="0" y1="0" x2="1" y2="0">'
             + "".join(f'<stop offset="{i / (len(RAMP) - 1):.4f}" stop-color="{c}"/>'
                       for i, c in enumerate(RAMP))
             + "</linearGradient></defs>")
    p.append(S.text(lx + BAR_W / 2, ly - 9, "−log10 FDR", cls="ax", anchor="middle"))
    p.append(S.rect(lx, ly, BAR_W, BAR_H, fill=f"url(#{uid})"))
    p.append(S.rect(lx, ly, BAR_W, BAR_H, fill="none", stroke="#c3c2b7",
                    stroke_width=1))
    for t in BAR_TICKS:
        tx = lx + ramp_t(t) * BAR_W
        p.append(S.line(tx, ly + BAR_H, tx, ly + BAR_H + 3, cls="base"))
        p.append(S.text(tx, ly + BAR_H + 14,
                        (f"≥{t:g}" if t >= FDR_CAP else f"{t:g}"),
                        cls="ax", anchor="middle"))

    lx += BAR_W + 54
    p.append(S.text(lx, ly + 9, "proteins", cls="ax"))
    lx += 46
    for k in (5, 20, k_max):
        p.append(S.circle(lx + R_MAX, ly + 5, radius(k), fill="#c3c2b7",
                          stroke="#ffffff", stroke_width=1.5))
        p.append(S.text(lx + R_MAX, ly + R_MAX + 19, str(k), cls="ax",
                        anchor="middle"))
        lx += 2 * R_MAX + 20
    if foot:
        p.append(S.text(14, h - 7, "†  KEGG global/overview or human-disease map — "
                                   "enriches on generic bacterial signal", cls="ax"))
    p.append("</svg>")
    return "\n".join(p), h


def to_pdf(svg: str, w: int, h: int, out: Path) -> bool:
    """Vector PDF from the SVG markup, text kept as text.

    Rendering lives in svg_lib.svg_to_pdf, which uses librsvg from this rule's
    conda env.  It replaced Chrome print-to-pdf, which was invoked through a
    hardcoded "google-chrome" - an undeclared dependency that no conda env
    could satisfy.  w and h are unused now: librsvg takes the page box from the
    SVG's own width/height, which is what the @page rule was reproducing.
    """
    try:
        S.svg_to_pdf(svg, out)
        return True
    except (RuntimeError, subprocess.TimeoutExpired) as e:
        print(f"  ! PDF export failed for {out.name}: {e}")
        return False


def main() -> int:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        input=Path(sm.params.indir),
        out=Path(sm.params.outdir),
        variant=sm.params.variant,
        no_pdf=False,
        shared_x=False,
    )
    outdir = args.out or (args.input / "plots")
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(args.input, args.variant)
    sfx = VARIANT_SUFFIX.get(args.variant, "")
    if not rows:
        sys.exit("no terms pass the threshold - nothing to plot")
    fold_max = max(r["fold"] for r in rows)
    k_max = max(r["k"] for r in rows)

    plots = []
    for cname, up_lab, down_lab in CONTRASTS:
        sel = [r for r in rows if r["contrast"] == cname]
        if not sel:
            print(f"  {cname}: no significant terms, skipped")
            continue
        fmax = fold_max if args.shared_x else max(r["fold"] for r in sel)
        svg, h = bubble_plot(sel, up_lab, down_lab, fmax, k_max,
                             uid=f"fdr_{cname}")
        (outdir / f"bubble_{cname}{sfx}.svg").write_text(svg)
        note = ""
        if not args.no_pdf:
            ok = to_pdf(svg, W, h, outdir / f"bubble_{cname}{sfx}.pdf")
            note = "  + pdf" if ok else "  (pdf FAILED)"
        plots.append((cname, svg))
        print(f"  bubble_{cname}{sfx}.svg   {len(sel)} terms   x axis 1x-{fmax:.1f}x{note}")


    print(f"\nwrote {outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
