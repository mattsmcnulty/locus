#!/usr/bin/env python3
"""Generate the Locus macOS app icon: a DNA double-helix on a blue→purple squircle.

Renders supersampled with Pillow, writes a full .iconset (all required sizes) plus a
preview PNG. Run via: uv run --with pillow python scripts/make_app_icon.py
Then build the .icns with: iconutil -c icns <iconset> -o Locus.icns
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/locus-icon")
S = 2048  # supersample master; downscaled per icon size
BG_TOP, BG_BOT = (20, 32, 78), (52, 26, 92)   # deep blue -> indigo background
COL_A, COL_B = (90, 178, 255), (181, 138, 255)  # blue + purple strands


def lerp(a: tuple, b: tuple, f: float) -> tuple:
    return tuple(int(a[i] + (b[i] - a[i]) * f) for i in range(3))


def render_master() -> Image.Image:
    # Background: vertical gradient masked to a rounded-rect (squircle-ish), like a macOS icon.
    grad = Image.new("RGBA", (S, S))
    gd = ImageDraw.Draw(grad)
    for y in range(S):
        gd.line([(0, y), (S, y)], fill=(*lerp(BG_TOP, BG_BOT, y / (S - 1)), 255))
    mask = Image.new("L", (S, S), 0)
    pad, rad = int(S * 0.055), int(S * 0.225)
    ImageDraw.Draw(mask).rounded_rectangle([pad, pad, S - pad, S - pad], radius=rad, fill=255)
    base = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    base.paste(grad, (0, 0), mask)

    cx = S / 2
    amp = S * 0.17           # helix half-width
    turns = 2.0              # number of crossings
    y0, y1 = S * 0.20, S * 0.80
    n = 700

    def strand(phase: float) -> list[tuple[float, float]]:
        return [(cx + amp * math.sin(turns * 2 * math.pi * (i / n) + phase), y0 + (y1 - y0) * (i / n))
                for i in range(n + 1)]

    a, b = strand(0.0), strand(math.pi)

    # Soft colored glow behind the strands (draw wide colored lines, then blur + composite).
    glow = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow)
    gw = int(S * 0.05)
    gdraw.line(a, fill=(*COL_A, 255), width=gw, joint="curve")
    gdraw.line(b, fill=(*COL_B, 255), width=gw, joint="curve")
    glow = glow.filter(ImageFilter.GaussianBlur(int(S * 0.022)))
    base = Image.alpha_composite(base, glow)

    d = ImageDraw.Draw(base)

    # Thin base-pair rungs; length shrinks near crossings (reads as 3D twist).
    rungs = 9
    rw = int(S * 0.010)
    for i in range(rungs):
        f = (i + 0.5) / rungs
        y = y0 + (y1 - y0) * f
        t = turns * 2 * math.pi * f
        ax, bx = cx + amp * math.sin(t), cx + amp * math.sin(t + math.pi)
        if abs(ax - bx) < amp * 0.25:   # skip near-crossings (would be dots)
            continue
        d.line([(ax, y), (bx, y)], fill=(231, 236, 255, 150), width=rw)

    # Two crisp two-tone strands with rounded ends.
    lw = int(S * 0.028)
    cap = lw / 2
    for pts, col in ((b, COL_B), (a, COL_A)):
        d.line(pts, fill=(*col, 255), width=lw, joint="curve")
        for x, y in (pts[0], pts[-1]):
            d.ellipse([x - cap, y - cap, x + cap, y + cap], fill=(*col, 255))
    return base


def main() -> None:
    master = render_master()
    iconset = OUT / "Locus.iconset"
    iconset.mkdir(parents=True, exist_ok=True)
    sizes = [
        (16, "16x16"), (32, "16x16@2x"), (32, "32x32"), (64, "32x32@2x"),
        (128, "128x128"), (256, "128x128@2x"), (256, "256x256"), (512, "256x256@2x"),
        (512, "512x512"), (1024, "512x512@2x"),
    ]
    for px, name in sizes:
        master.resize((px, px), Image.LANCZOS).save(iconset / f"icon_{name}.png")
    master.resize((512, 512), Image.LANCZOS).save(OUT / "preview.png")
    print(f"iconset: {iconset}")
    print(f"preview: {OUT / 'preview.png'}")


if __name__ == "__main__":
    main()
