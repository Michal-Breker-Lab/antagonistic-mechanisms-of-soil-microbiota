#!/usr/bin/env python3
"""Minimal dependency-free SVG builder, plus PDF export via librsvg."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
       "#4a3aa7", "#e34948"]
CAT_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300",
            "#9085e9", "#e66767"]
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
       "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
       "#0d366b"]
STATUS = {"good": "#0ca30c", "warning": "#fab219",
          "serious": "#ec835a", "critical": "#d03b3b"}

FONT = 'Arial, Helvetica, "Liberation Sans", sans-serif'

STYLE = """
  .surface { fill: #fcfcfb; }
  text { font-family: %(font)s; fill: #0b0b0b; }
  .ttl  { font-size: 15px; font-weight: 600; }
  .sub  { font-size: 12px; fill: #52514e; }
  .ax   { font-size: 10.5px; fill: #898781; }
  .axb  { font-size: 11px; fill: #52514e; }
  .lgd  { font-size: 11.5px; fill: #52514e; }
  .val  { font-size: 10.5px; fill: #0b0b0b; }
  .grid { stroke: #e1e0d9; stroke-width: 1; }
  .base { stroke: #c3c2b7; stroke-width: 1; }
  .gap  { stroke: #fcfcfb; }
  .ring { stroke: #fcfcfb; stroke-width: 2; }
  .s1 { fill: #2a78d6; } .s2 { fill: #eb6834; } .s3 { fill: #1baf7a; }
  .s4 { fill: #eda100; } .s5 { fill: #e87ba4; } .s6 { fill: #008300; }
  .s7 { fill: #4a3aa7; } .s8 { fill: #e34948; }
  .st1 { stroke: #2a78d6; }
  .pt0 { fill: #b9b8b2; }   /* unhighlighted background points */
  .pt4 { fill: #52514e; }   /* highlighted but no categorical hue */
  /* Ordinal ramp for ranked tiers (A best -> D weakest): ONE hue, stepped by
     rank, never categorical hues - the categories are ordered, so hue would
     throw that ordering away.  On the light surface the palest step is 250,
     the documented floor for an ordinal ramp (2.06:1 against the surface). */
  .t1 { fill: #0d366b; } .t2 { fill: #1c5cab; }
  .t3 { fill: #3987e5; } .t4 { fill: #86b6ef; }
  .thr { stroke: #898781; stroke-width: 1; stroke-dasharray: 4 3; }
  .band { fill: #f0efec; }
  @media (prefers-color-scheme: dark) {
    .surface { fill: #1a1a19; }
    text { fill: #ffffff; }
    .sub, .axb, .lgd { fill: #c3c2b7; }
    .ax { fill: #898781; }
    .val { fill: #ffffff; }
    .grid { stroke: #2c2c2a; }
    .base { stroke: #383835; }
    .gap, .ring { stroke: #1a1a19; }
    .s1 { fill: #3987e5; } .s2 { fill: #d95926; } .s3 { fill: #199e70; }
    .s4 { fill: #c98500; } .s5 { fill: #d55181; } .s6 { fill: #008300; }
    .s7 { fill: #9085e9; } .s8 { fill: #e66767; }
    .st1 { stroke: #3987e5; }
    .pt0 { fill: #4a4a47; }
    .pt4 { fill: #c3c2b7; }
    .band { fill: #232321; }
    /* mirrored for the dark surface: the ramp runs the other way, and the
       darkest step is 600 - the documented floor there (2.15:1) */
    .t1 { fill: #cde2fb; } .t2 { fill: #9ec5f4; }
    .t3 { fill: #5598e7; } .t4 { fill: #184f95; }
  }
""" % {"font": FONT}


def esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def head(w: int, h: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
        f'width="{w}" height="{h}" role="img" aria-label="{esc(title)}">',
        f"<style>{STYLE}</style>",
        f'<rect class="surface" x="0" y="0" width="{w}" height="{h}"/>',
    ]


def text(x, y, s, cls="ax", anchor="start", rotate=None, extra="",
         fill=None) -> str:
    """fill is emitted as an inline style, NOT a presentation attribute.

    A CSS class rule always beats a `fill=` attribute, so `fill="red"` on an
    element that also carries `class="val"` silently loses to `.val{fill:...}`.
    Inline style wins over both, which is what a per-element override needs.
    """
    tr = f' transform="rotate({rotate},{x:.1f},{y:.1f})"' if rotate else ""
    st = f' style="fill:{fill}"' if fill else ""
    return (f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" '
            f'text-anchor="{anchor}"{tr}{st}{extra}>{esc(s)}</text>')


def _close(tag: str, body: str, title: str | None) -> str:
    """Self-close, or wrap a <title> child so hover works with no JS."""
    if title:
        return f"<{tag}{body}><title>{esc(title)}</title></{tag}>"
    return f"<{tag}{body}/>"


def rect(x, y, w, h, fill=None, cls=None, rx=0, extra="", title=None,
         stroke=None, stroke_width=None) -> str:
    b = (f' x="{x:.2f}" y="{y:.2f}" width="{max(w,0):.2f}"'
         f' height="{max(h,0):.2f}"')
    if cls:
        b = f' class="{cls}"' + b
    if fill:
        b += f' fill="{fill}"'
    if stroke:
        b += f' stroke="{stroke}"'
    if stroke_width:
        b += f' stroke-width="{stroke_width}"'
    if rx:
        b += f' rx="{rx}"'
    return _close("rect", b + extra, title)


def line(x1, y1, x2, y2, cls="grid", extra="", stroke=None, stroke_width=None,
         dash=None, title=None) -> str:
    b = (f' class="{cls}" x1="{x1:.1f}" y1="{y1:.1f}"'
         f' x2="{x2:.1f}" y2="{y2:.1f}"')
    if stroke:
        b += f' stroke="{stroke}"'
    if stroke_width:
        b += f' stroke-width="{stroke_width}"'
    if dash:
        b += f' stroke-dasharray="{dash}"'
    return _close("line", b + extra, title)


def circle(cx, cy, r, fill="none", extra="", title=None, stroke=None,
           stroke_width=None, cls=None, opacity=None) -> str:
    b = f' cx="{cx:.2f}" cy="{cy:.2f}" r="{r}"'
    b += f' fill="{fill}"' if fill and not cls else ""
    if cls:
        b = f' class="{cls}"' + b
    if opacity is not None:
        b += f' opacity="{opacity}"'
    if stroke:
        b += f' stroke="{stroke}"'
    if stroke_width:
        b += f' stroke-width="{stroke_width}"'
    return _close("circle", b + extra, title)


def triangle(cx, cy, r, fill=None, up=True, cls=None, stroke=None,
             stroke_width=None, title=None, extra="") -> str:
    """Up/down triangle - used where a mark must read as a different kind of
    thing from the volcano's circles, not merely a different colour."""
    d = 1 if up else -1
    pts = (f"{cx:.2f},{cy-d*r*1.1:.2f} {cx-r:.2f},{cy+d*r*0.75:.2f} "
           f"{cx+r:.2f},{cy+d*r*0.75:.2f}")
    b = f' points="{pts}"'
    if cls:
        b = f' class="{cls}"' + b
    elif fill:
        b += f' fill="{fill}"'
    if stroke:
        b += f' stroke="{stroke}"'
    if stroke_width:
        b += f' stroke-width="{stroke_width}"'
    return _close("polygon", b + extra, title)


def square(cx, cy, r, fill=None, extra="", title=None, stroke=None,
           stroke_width=None, cls=None) -> str:
    s = r * 0.9
    return rect(cx - s, cy - s, 2 * s, 2 * s, fill=fill, cls=cls, rx=1.5,
                extra=extra, title=title, stroke=stroke,
                stroke_width=stroke_width)


def _poly(pts, cls, fill, stroke, stroke_width, title, extra) -> str:
    b = ' points="' + " ".join(f"{x:.2f},{y:.2f}" for x, y in pts) + '"'
    if cls:
        b = f' class="{cls}"' + b
    elif fill:
        b += f' fill="{fill}"'
    if stroke:
        b += f' stroke="{stroke}"'
    if stroke_width:
        b += f' stroke-width="{stroke_width}"'
    return _close("polygon", b + extra, title)


def diamond(cx, cy, r, fill=None, cls=None, stroke=None, stroke_width=None,
            title=None, extra="") -> str:
    d = r * 1.28
    return _poly([(cx, cy - d), (cx + d, cy), (cx, cy + d), (cx - d, cy)],
                 cls, fill, stroke, stroke_width, title, extra)


def _crossing(cx, cy, r, arm, rot) -> list[tuple[float, float]]:
    """The 12 corners of a + (rot=0) or x (rot=45), as a single polygon."""
    import math as _m
    a, b, out = r * 1.3, r * arm, []
    c, s = _m.cos(_m.radians(rot)), _m.sin(_m.radians(rot))
    for px, py in [(b, a), (-b, a), (-b, b), (-a, b), (-a, -b), (-b, -b),
                   (-b, -a), (b, -a), (b, -b), (a, -b), (a, b), (b, b)]:
        out.append((cx + px * c - py * s, cy + px * s + py * c))
    return out


def plus(cx, cy, r, fill=None, cls=None, stroke=None, stroke_width=None,
         title=None, extra="") -> str:
    return _poly(_crossing(cx, cy, r, 0.36, 0), cls, fill, stroke,
                 stroke_width, title, extra)


def cross(cx, cy, r, fill=None, cls=None, stroke=None, stroke_width=None,
          title=None, extra="") -> str:
    return _poly(_crossing(cx, cy, r, 0.36, 45), cls, fill, stroke,
                 stroke_width, title, extra)


def star(cx, cy, r, fill=None, cls=None, stroke=None, stroke_width=None,
         title=None, extra="") -> str:
    """Five-pointed star - the most shape-distinct mark available here.

    Reach for it when a series must be unmistakable against every other mark on
    the figure, including the up/down triangles that carry direction.
    """
    import math as _m
    R, ri, pts = r * 1.5, r * 1.5 * 0.42, []
    for i in range(10):
        a = _m.radians(-90 + i * 36)
        rad = R if i % 2 == 0 else ri
        pts.append((cx + rad * _m.cos(a), cy + rad * _m.sin(a)))
    return _poly(pts, cls, fill, stroke, stroke_width, title, extra)


MARKS = {"circle": circle, "square": square, "diamond": diamond,
         "plus": plus, "cross": cross, "star": star, "triangle": triangle,
         "triangle_down": lambda cx, cy, r, **kw: triangle(cx, cy, r,
                                                           up=False, **kw)}


def mark(shape: str, cx, cy, r, **kw) -> str:
    """Dispatch to one of MARKS; unknown shapes fall back to a circle."""
    return MARKS.get(shape, circle)(cx, cy, r, **kw)


def nice_ticks(lo: float, hi: float, n: int = 5) -> list[float]:
    """Human-readable tick positions covering [lo, hi]."""
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / n
    mag = 10 ** (len(str(int(abs(raw)))) - 1) if abs(raw) >= 1 else \
        10 ** -(len(str(abs(raw)).split(".")[1]) - len(str(abs(raw)).split(".")[1].lstrip("0")) + 1)
    for m in (1, 2, 2.5, 5, 10):
        if raw <= mag * m:
            step = mag * m
            break
    else:
        step = mag * 10
    t, out = (int(lo / step)) * step, []
    while t <= hi + step * 0.001:
        if t >= lo - step * 0.001:
            out.append(round(t, 10))
        t += step
    return out


def fmt(v: float) -> str:
    a = abs(v)
    if a >= 1_000_000:
        return f"{v/1_000_000:g}M"
    if a >= 1_000:
        return f"{v/1_000:g}k"
    if a == int(a):
        return str(int(v))
    return f"{v:g}"


def legend(x, y, items: list[tuple[str, str]], gap=118) -> list[str]:
    """items = [(label, fill), ...] - always present for >= 2 series."""
    out = []
    for i, (lab, fill) in enumerate(items):
        cx = x + i * gap
        out.append(rect(cx, y - 8, 10, 10, fill=fill, rx=2))
        out.append(text(cx + 15, y + 1, lab, cls="lgd"))
    return out



def svg_to_pdf(svg, out) -> None:
    """Vector PDF with selectable, editable text - for Illustrator.

    rsvg-convert (librsvg) keeps every <text> as a real text object with an
    embedded font subset, which is the property these figures exist for;
    rasterisers and path-outlining exporters do not.  It replaced Chrome's
    print-to-pdf, which produced the same thing but only if a browser happened
    to be installed on the machine - an undeclared dependency no conda env
    could satisfy, since no channel packages Chrome or Chromium.  librsvg is in
    figures.yaml, so the rule now carries its own renderer.

    Verified against the Chrome output it replaces, over all eight figures:
    identical extractable text, same embedded font subset, mean pixel
    difference 0.04-0.11%.  The page box is 0.02% smaller (Chrome rounded the
    px page up to whole pixels) - immaterial at figure scale.

    The SVG's prefers-color-scheme block needs no stripping: librsvg has no
    media-query state and renders the light-mode declarations, which is exactly
    what the Chrome path forced by hand.

    `svg` is a Path to an SVG file, or the SVG markup as a str.
    Raises if librsvg is missing - these rules declare .pdf outputs, so a
    silent skip would surface later as a confusing missing-output error.
    """
    exe = shutil.which("rsvg-convert")
    if exe is None:
        raise RuntimeError(
            "rsvg-convert not found: PDF export needs librsvg, which is "
            "declared in workflow/envs/figures.yaml. Run snakemake with "
            "--sdm conda so the rule gets its environment."
        )
    out = Path(out)
    if isinstance(svg, Path):
        markup = None
    elif isinstance(svg, str) and svg.lstrip()[:1] == "<":
        markup = svg
    else:
        markup = None
    cmd = [exe, "-f", "pdf", "-o", str(out)]
    if markup is None:
        cmd.append(str(svg))
        r = subprocess.run(cmd, capture_output=True, timeout=300)
    else:
        r = subprocess.run(cmd, input=markup.encode(), capture_output=True, timeout=300)
    if r.returncode != 0 or not out.is_file() or out.stat().st_size == 0:
        raise RuntimeError(
            f"rsvg-convert failed for {out.name}: "
            f"{r.stderr.decode(errors='replace').strip()[:400]}"
        )
