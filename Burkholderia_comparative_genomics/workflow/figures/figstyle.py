"""Shared figure style enforcing the Illustrator-editable-text rule.

Per CLAUDE.md: matplotlib's defaults silently break Illustrator editing in BOTH
vector formats -- PDF gets Type 3 fonts (imported as uneditable art) and SVG gets
svg.fonttype='path' (every glyph converted to a <path>, zero <text> elements).
These rcParams must be set BEFORE any figure is created, not just before savefig.

Arial is also registered explicitly: with svg.fonttype='none' the SVG only NAMES
the font, so matplotlib's default DejaVu Sans -- not a system font on Windows or
macOS -- would make Illustrator substitute and reflow every label.
"""
from pathlib import Path

import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42      # TrueType, not Type 3
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"  # keep <text>; do not outline glyphs

from matplotlib import font_manager as fm  # noqa: E402

for _stem in ("arial", "arialbd", "ariali", "arialbi"):
    _p = Path("/mnt/c/Windows/Fonts") / f"{_stem}.ttf"
    if _p.exists():
        try:
            fm.fontManager.addfont(_p)
        except Exception:  # noqa: BLE001
            pass

if "Arial" in {f.name for f in fm.fontManager.ttflist}:
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    ARIAL = True
else:
    ARIAL = False

matplotlib.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
})

# colourblind-safe palette for host categories
HOST_COLORS = {
    "human_clinical": "#D55E00",
    "animal": "#E69F00",
    "insect": "#F0E442",
    "plant": "#009E73",
    "rhizosphere": "#56B4E9",
    "soil": "#8C6D31",
    "water": "#0072B2",
    "fungus": "#CC79A7",
    "amoeba_protist": "#7B3294",
    "industrial": "#999999",
    "environmental_other": "#BBBBBB",
    "unknown": "#EEEEEE",
}
C3_COLORS = {True: "#1B3A6B", False: "#D9D9D9"}


def save(fig, path_noext, formats=("pdf", "svg", "png")):
    """Save and verify the vector outputs really contain editable text."""
    import re
    out = []
    for ext in formats:
        p = f"{path_noext}.{ext}"
        fig.savefig(p, format=ext)
        out.append(p)
    ok = True
    svg = f"{path_noext}.svg"
    if Path(svg).exists():
        n = open(svg, encoding="utf-8", errors="ignore").read().count("<text")
        print(f"  {Path(svg).name}: {n} <text> elements", end="")
        if n == 0:
            print("  <-- FAIL: glyphs were outlined"); ok = False
        else:
            print("  OK")
    pdf = f"{path_noext}.pdf"
    if Path(pdf).exists():
        d = open(pdf, "rb").read()
        subs = sorted({s.decode() for s in re.findall(rb"/Subtype\s*/(\w+)", d)})
        fonty = [s for s in subs if "Type" in s or "Font" in s]
        print(f"  {Path(pdf).name}: font subtypes {fonty}", end="")
        if any(s == "Type3" for s in subs):
            print("  <-- FAIL: Type 3 fonts"); ok = False
        else:
            print("  OK")
    if not ARIAL:
        print("  NOTE: Arial not registered - labels will reflow in Illustrator")
    return ok
