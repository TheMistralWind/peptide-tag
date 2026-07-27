"""Drawing the peptide.

One geometry routine feeds three surfaces: the inline hero on the page, a
specimen card you can download or print, and a 1200x630 social image.

The chain wraps like text, left to right on every row. An earlier prototype
wrapped serpentine, which is prettier as a chain but means the second row of
ROOSAMARIA reads AIRAM, and the entire point of the toy is reading your own
name.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path

import chemistry as C

FONT_DIR = Path(__file__).resolve().parent / "static" / "fonts"
MONO_TTF = FONT_DIR / "JetBrainsMono-Regular.ttf"
SERIF_TTF = FONT_DIR / "Newsreader.ttf"

# Measured for contrast: on paper, ink is 17.5:1, ink_70 is 9.3:1, ink_50 is
# 4.9:1 and accent is 4.6:1. Dark mode is 16.8 / 5.6 / 5.9:1. Nothing on the
# page falls below WCAG AA.
LIGHT = {
    "paper": "#F7F5F0", "ink": "#12100E", "ink70": "#45413B",
    "ink50": "#6E6A62", "ink20": "#DEDAD1", "accent": "#C2452D",
}
DARK = {
    "paper": "#0F0E0C", "ink": "#F2EFE9", "ink70": "#C4BFB6",
    "ink50": "#8F8A80", "ink20": "#35322B", "accent": "#E4694D",
}

SERIF = ("ui-serif, Charter, 'Bitstream Charter', 'Iowan Old Style', "
         "'Palatino Linotype', Georgia, serif")
MONO = ("ui-monospace, 'SF Mono', SFMono-Regular, Menlo, Consolas, "
        "'Roboto Mono', 'DejaVu Sans Mono', monospace")

# Tile treatment per residue class. Three redundant channels carry meaning
# (fill tone, the +/- glyph, and the letter itself), so the artwork survives
# greyscale printing and every form of colour blindness.
SOLID = "solid"      # knocked-out letter on ink
TINT = "tint"        # letter on a light fill
OUTLINE = "outline"  # letter on paper, ruled edge

TILE_STYLE = {
    C.HYDROPHOBIC: SOLID, C.AROMATIC: SOLID, C.SPECIAL: SOLID,
    C.POLAR: TINT, C.ACIDIC: OUTLINE, C.BASIC: OUTLINE, C.RARE: OUTLINE,
}


@dataclass
class Layout:
    cols: int
    rows: int
    tile: int
    gap: int
    row_gap: int
    width: int
    height: int


def plan(n: int, max_width: int, max_cols: int = 10,
         min_tile: int = 32, max_tile: int = 96) -> Layout:
    """Balance the chain into rows that fit max_width without orphan tiles."""
    n = max(n, 1)
    rows = max(1, -(-n // max_cols))
    cols = -(-n // rows)  # balances 10 into 5+5 rather than 8+2
    tile = int(max_width / (1.2 * cols - 0.2)) if cols else max_tile
    tile = max(min_tile, min(max_tile, tile))
    gap = round(0.20 * tile)
    row_gap = round(0.34 * tile)
    width = cols * (tile + gap) - gap
    height = rows * (tile + row_gap) - row_gap
    return Layout(cols, rows, tile, gap, row_gap, width, height)


def _esc(s: str) -> str:
    return html.escape(s, quote=True)


def chain_svg(sequence: str, max_width: int = 760, max_cols: int = 10,
              theme: dict | None = None, inline: bool = True,
              pad: int = 0) -> str:
    """The residue chain. Inline mode inherits colour from CSS variables."""
    if not sequence:
        return ""
    lay = plan(len(sequence), max_width, max_cols)
    t = theme or LIGHT

    def col(name: str) -> str:
        # Inline on the page, the page's own custom properties drive colour, so
        # light and dark mode need no second rendering.
        return f"var(--{name})" if inline else t[name]

    # Room for the H2N and OH labels on either end.
    label = max(11, round(0.24 * lay.tile))
    lead = label * 3
    w = lay.width + lead * 2 + pad * 2
    h = lay.height + pad * 2
    ox, oy = lead + pad, pad

    parts: list[str] = []
    letters = ", ".join(C.RESIDUES[a].name.lower() for a in sequence)
    title = f"Peptide {sequence}: {letters}."

    for i, a in enumerate(sequence):
        r, c = divmod(i, lay.cols)
        x = ox + c * (lay.tile + lay.gap)
        y = oy + r * (lay.tile + lay.row_gap)
        res = C.RESIDUES[a]
        style = TILE_STYLE[res.klass]

        if style == SOLID:
            fill, stroke, text_col = col("ink"), "none", col("paper")
        elif style == TINT:
            fill, stroke, text_col = col("ink20"), "none", col("ink")
        else:
            edge = col("accent") if res.klass == C.RARE else col("ink")
            fill, stroke, text_col = col("paper"), edge, col("ink")

        if stroke == "none":
            parts.append(
                f'<rect x="{x}" y="{y}" width="{lay.tile}" height="{lay.tile}" '
                f'rx="3" fill="{fill}"/>'
            )
        else:
            inset = 0.75
            parts.append(
                f'<rect x="{x + inset}" y="{y + inset}" '
                f'width="{lay.tile - 2 * inset}" height="{lay.tile - 2 * inset}" '
                f'rx="3" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>'
            )

        parts.append(
            f'<text x="{x + lay.tile / 2}" y="{y + lay.tile / 2 + 0.185 * lay.tile}" '
            f'font-family="{MONO}" font-size="{0.50 * lay.tile:.1f}" '
            f'text-anchor="middle" fill="{text_col}">{a}</text>'
        )

        charge = C.CHARGE.get(a, 0)
        if charge:
            parts.append(
                f'<text x="{x + 0.84 * lay.tile}" y="{y + 0.28 * lay.tile}" '
                f'font-family="{MONO}" font-size="{0.26 * lay.tile:.1f}" '
                f'text-anchor="middle" fill="{col("accent")}">'
                f'{"+" if charge > 0 else "−"}</text>'
            )

        # Peptide bond to the next residue, within a row only.
        if i + 1 < len(sequence) and c + 1 < lay.cols:
            ym = y + lay.tile / 2
            parts.append(
                f'<line x1="{x + lay.tile}" y1="{ym}" '
                f'x2="{x + lay.tile + lay.gap}" y2="{ym}" '
                f'stroke="{col("ink")}" stroke-width="2"/>'
            )

    ym = oy + lay.tile / 2
    parts.append(
        f'<text x="{ox - label}" y="{ym + label * 0.35}" font-family="{MONO}" '
        f'font-size="{label}" text-anchor="end" fill="{col("ink50")}">H₂N</text>'
    )
    last_r, last_c = divmod(len(sequence) - 1, lay.cols)
    lx = ox + last_c * (lay.tile + lay.gap) + lay.tile + label
    ly = oy + last_r * (lay.tile + lay.row_gap) + lay.tile / 2 + label * 0.35
    parts.append(
        f'<text x="{lx}" y="{ly}" font-family="{MONO}" font-size="{label}" '
        f'text-anchor="start" fill="{col("ink50")}">OH</text>'
    )

    body = "".join(parts)
    # A viewBox is what makes this scale. Its absence is why the old tag broke.
    # width/height come from CSS, not attributes: height="auto" is not valid SVG
    # and makes the drawing stretch to fill its parent.
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
        f'style="display:block;width:100%;height:auto;max-width:{w}px" '
        f'role="img" aria-label="{_esc(title)}">'
        f'<title>{_esc(title)}</title>{body}</svg>'
    )


def specimen_svg(name: str, sequence: str, props: dict,
                 theme: dict | None = None, width: int = 1000) -> str:
    """The downloadable card: name, letter chain, the helix, and the numbers.

    Everything is stacked from an explicit running cursor. The previous version
    hard-coded a 150px header and then translated the chain to exactly 150,
    which put the tiles straight through the "THE PEPTIDE OF THAT NAME" line and
    left a large void underneath.
    """
    t = theme or LIGHT
    margin = 80
    inner = width - margin * 2

    lay = plan(len(sequence), inner - 120, max_cols=10)
    chain = chain_svg(sequence, max_width=inner - 120, theme=t, inline=False)
    chain_w = lay.width + max(11, round(0.24 * lay.tile)) * 6

    helix_w = inner
    helix_h = 300 if lay.rows <= 2 else 260
    helix = helix_svg(sequence, width=helix_w, height=helix_h, theme=t, inline=False)

    head_h = 178                       # name plus eyebrow, with air beneath
    chain_h = lay.height + 56
    foot_h = 152                       # label, sequence, then the numbers
    height = head_h + chain_h + helix_h + foot_h

    # The sequence gets its own line and shrinks to fit. Sharing a line with the
    # numbers meant a 24 residue name ran straight through the molecular weight.
    letter_spacing = 2.5 if len(sequence) <= 28 else 1.2
    seq_size = 21.0
    if sequence:
        fits = (inner / len(sequence) - letter_spacing) / 0.6
        seq_size = max(10.0, min(21.0, fits))

    mw = props.get("molecular_weight", 0.0)
    pi = props.get("isoelectric_point", 0.0)
    formula = C.format_formula_unicode(props.get("formula", ""))

    def inner_of(svg: str) -> str:
        return svg.split(">", 1)[1].rsplit("</svg>", 1)[0] if svg else ""

    chain_y = head_h
    helix_y = head_h + chain_h

    # SVG collapses runs of whitespace, so the metrics are three positioned
    # fields rather than one string padded with spaces.
    third = inner / 3.0
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">
<rect width="{width}" height="{height}" fill="{t['paper']}"/>
<text x="{margin}" y="{margin + 34}" font-family="{SERIF}" font-size="46" fill="{t['ink']}">{_esc(name)}</text>
<text x="{margin}" y="{margin + 68}" font-family="{MONO}" font-size="14" letter-spacing="1.6" fill="{t['ink50']}">THE PEPTIDE OF THAT NAME</text>
<g transform="translate({(width - chain_w) / 2:.1f}, {chain_y})">{inner_of(chain)}</g>
<g transform="translate({margin}, {helix_y})">{inner_of(helix)}</g>
<text x="{margin}" y="{height - 100}" font-family="{MONO}" font-size="13" letter-spacing="1.6" fill="{t['ink50']}">IDEALISED ALPHA HELIX, A STYLISATION</text>
<text x="{margin}" y="{height - 62}" font-family="{MONO}" font-size="{seq_size:.1f}" letter-spacing="{letter_spacing}" fill="{t['ink']}">{_esc(sequence)}</text>
<text x="{margin}" y="{height - 30}" font-family="{MONO}" font-size="15" fill="{t['ink50']}">{mw:.1f} Da</text>
<text x="{margin + third:.0f}" y="{height - 30}" font-family="{MONO}" font-size="15" fill="{t['ink50']}">pI {pi:.2f}</text>
<text x="{width - margin}" y="{height - 30}" font-family="{MONO}" font-size="15" text-anchor="end" fill="{t['ink50']}">{_esc(formula)}</text>
</svg>"""


# ---------------------------------------------------------------------------
# Social card
#
# This is drawn with Pillow rather than rasterised from the SVG above. The
# obvious route, cairosvg, needs libcairo present at runtime, which means an
# apt layer in a container that is deliberately kept lean. Pillow is already
# in the tree as an RDKit dependency and needs nothing but the two bundled
# font files. Slack, iMessage, WhatsApp and X will not render an SVG in
# og:image, so a PNG is not optional here.
# ---------------------------------------------------------------------------

OG_W, OG_H = 1200, 630


def _font(path, size: int):
    from PIL import ImageFont
    return ImageFont.truetype(str(path), size)


def social_png(name: str, sequence: str, props: dict) -> bytes:
    from PIL import Image, ImageDraw

    t = LIGHT
    img = Image.new("RGB", (OG_W, OG_H), t["paper"])
    d = ImageDraw.Draw(img)

    margin = 80
    f_name = _font(SERIF_TTF, 68)
    f_label = _font(MONO_TTF, 17)
    f_seq = _font(MONO_TTF, 26)
    f_meta = _font(MONO_TTF, 20)

    display = name.strip() or sequence
    while d.textlength(display, font=f_name) > OG_W - margin * 2 and len(display) > 4:
        display = display[:-2]
        if len(display) < len(name.strip()):
            display = display.rstrip() + "…"

    d.text((margin, margin - 12), display, font=f_name, fill=t["ink"])
    d.text((margin, margin + 76), "THE PEPTIDE OF THAT NAME", font=f_label,
           fill=t["ink50"])

    # Chain, using the same geometry as the SVG so the two always agree.
    avail = OG_W - margin * 2
    lay = plan(len(sequence), avail, max_cols=10, min_tile=34, max_tile=88)
    # Centre the chain in the band between the heading and the footer rather
    # than pinning it, so a one-row chain does not sit high with a void beneath.
    band_top, band_bottom = 200, 470
    top = band_top + max(0, (band_bottom - band_top - lay.height) / 2)
    left = margin + (avail - lay.width) / 2

    for i, a in enumerate(sequence):
        r, c = divmod(i, lay.cols)
        x = left + c * (lay.tile + lay.gap)
        y = top + r * (lay.tile + lay.row_gap)
        res = C.RESIDUES[a]
        style = TILE_STYLE[res.klass]
        box = (x, y, x + lay.tile, y + lay.tile)

        if style == SOLID:
            d.rounded_rectangle(box, radius=3, fill=t["ink"])
            text_col = t["paper"]
        elif style == TINT:
            d.rounded_rectangle(box, radius=3, fill=t["ink20"])
            text_col = t["ink"]
        else:
            edge = t["accent"] if res.klass == C.RARE else t["ink"]
            d.rounded_rectangle(box, radius=3, fill=t["paper"], outline=edge, width=2)
            text_col = t["ink"]

        f_tile = _font(MONO_TTF, int(lay.tile * 0.52))
        d.text((x + lay.tile / 2, y + lay.tile / 2), a, font=f_tile,
               fill=text_col, anchor="mm")

        charge = C.CHARGE.get(a, 0)
        if charge:
            f_ch = _font(MONO_TTF, max(11, int(lay.tile * 0.28)))
            d.text((x + lay.tile * 0.84, y + lay.tile * 0.26),
                   "+" if charge > 0 else "−", font=f_ch,
                   fill=t["accent"], anchor="mm")

        if i + 1 < len(sequence) and c + 1 < lay.cols:
            ym = y + lay.tile / 2
            d.line((x + lay.tile, ym, x + lay.tile + lay.gap, ym),
                   fill=t["ink"], width=2)

    seq_text = sequence if len(sequence) <= 40 else sequence[:39] + "…"
    d.text((margin, OG_H - 128), seq_text, font=f_seq, fill=t["ink"])
    meta = (f"{props.get('molecular_weight', 0):.1f} Da     "
            f"pI {props.get('isoelectric_point', 0):.2f}     "
            f"{C.format_formula_unicode(props.get('formula', ''))}")
    d.text((margin, OG_H - 86), meta, font=f_meta, fill=t["ink50"])
    d.text((OG_W - margin, OG_H - 86), "peptide tag", font=f_meta,
           fill=t["ink50"], anchor="ra")

    buf = BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# The helix
#
# An orthographic projection of the ideal alpha helix from structure.py, drawn
# ball and stick in the same two inks as everything else. Atoms and bonds get a
# paper-coloured halo so that overlapping geometry still reads, which is the
# trick that makes a flat monochrome projection legible without colour.
# ---------------------------------------------------------------------------

# Relative ball sizes. Carbon is the reference.
_ATOM_R = {"C": 1.0, "N": 0.95, "O": 0.92, "S": 1.25, "SE": 1.35, "H": 0.6}
_BOND_MAX = 1.85


def _parse_pdb(pdb: str):
    import numpy as np
    elements, coords = [], []
    for line in pdb.splitlines():
        if line.startswith("ATOM"):
            elements.append((line[76:78].strip() or line[12:16].strip()[:1]).upper())
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return elements, np.array(coords) if coords else None


def _orient(coords):
    """Rotate so the long axis of the molecule runs left to right."""
    import numpy as np
    centred = coords - coords.mean(axis=0)
    # Principal axes. The first is the helix axis, which we want horizontal.
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    rotated = centred @ vt.T
    return rotated


def helix_svg(sequence: str, width: int = 760, height: int = 260,
              theme: dict | None = None, inline: bool = True) -> str:
    """Ball and stick projection of the idealised alpha helix."""
    import numpy as np

    import structure

    pdb = structure.model_pdb(sequence, "helix")
    if not pdb:
        return ""
    elements, coords = _parse_pdb(pdb)
    if coords is None or len(coords) < 2:
        return ""

    t = theme or LIGHT

    def col(name: str) -> str:
        return f"var(--{name})" if inline else t[name]

    pts = _orient(coords)
    xy = pts[:, :2]
    depth = pts[:, 2]

    pad = 26
    span = xy.max(axis=0) - xy.min(axis=0)
    span[span == 0] = 1.0
    scale = min((width - pad * 2) / span[0], (height - pad * 2) / span[1])
    screen = (xy - xy.min(axis=0)) * scale
    screen[:, 0] += (width - span[0] * scale) / 2
    # SVG y grows downward; flip so the molecule is not upside down.
    screen[:, 1] = (height - span[1] * scale) / 2 + (span[1] * scale - screen[:, 1])

    ball = max(4.5, scale * 0.42)
    stick = max(1.8, scale * 0.13)
    halo = max(1.4, stick * 0.85)

    bonds = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if float(np.linalg.norm(pts[i] - pts[j])) < _BOND_MAX:
                bonds.append((i, j))

    # Painter's algorithm: everything sorted back to front, so near geometry
    # overwrites far geometry and the halos do the occlusion work.
    items = [("bond", i, j, (depth[i] + depth[j]) / 2) for i, j in bonds]
    items += [("atom", i, i, depth[i]) for i in range(len(pts))]
    items.sort(key=lambda it: it[3])

    lo, hi = float(depth.min()), float(depth.max())
    rng = (hi - lo) or 1.0

    parts: list[str] = []
    for kind, i, j, d in items:
        # Far geometry is drawn slightly lighter, which reads as depth without
        # introducing a second colour.
        fade = 0.55 + 0.45 * ((d - lo) / rng)
        if kind == "bond":
            x1, y1 = screen[i]
            x2, y2 = screen[j]
            # No halo on the sticks. Giving them one carved paper channels
            # everywhere and the whole thing read as outlined noodles.
            parts.append(
                f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                f'stroke="{col("ink")}" stroke-width="{stick:.1f}" '
                f'stroke-linecap="round" opacity="{fade:.2f}"/>'
            )
        else:
            x, y = screen[i]
            r = ball * _ATOM_R.get(elements[i], 1.0)
            # Every atom is the same ink. Element is carried by size, not by a
            # second colour, and the paper ring is what separates one ball from
            # the one behind it.
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r + halo:.1f}" '
                f'fill="{col("paper")}"/>'
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                f'fill="{col("ink")}" opacity="{fade:.2f}"/>'
            )

    body = "".join(parts)
    label = f"Idealised alpha helix of the peptide {sequence}."
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'style="display:block;width:100%;height:auto;max-width:{width}px" '
        f'role="img" aria-label="{_esc(label)}">'
        f'<title>{_esc(label)}</title>{body}</svg>'
    )
