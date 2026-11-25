"""Calibrated reader for Fig. 5A bottom (p. 5): 1D smooth epicuticle surface height profile h(x) at 35% dev."""
from __future__ import annotations

import sys
from pathlib import Path
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_box, find_ticks_mask, neutral_dark_mask, trace_black_curve_at


def read(image, at):
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    # Inward x-ticks along bottom border: 0, 2, 4, 6 um
    xt = find_ticks_mask(dark, r1, c0 + 15, c1 - 15, "x", -1, offset_lo=10, offset_hi=22)
    if len(xt) != 4:
        raise RuntimeError(f"Fig. 5A bottom expected 4 x-ticks (0,2,4,6), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [0.0, 2.0, 4.0, 6.0])))

    # Inward y-ticks along left border: 150, 100, 50, 0 nm
    yt = find_ticks_mask(dark, c0, r0 + 15, r1 - 15, "y", +1, offset_lo=10, offset_hi=22)
    if len(yt) != 4:
        raise RuntimeError(f"Fig. 5A bottom expected 4 y-ticks (150,100,50,0), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [150.0, 100.0, 50.0, 0.0])))

    xs = [float(v) for v in at] if at else [round(0.1 * i, 4) for i in range(1, 71)]
    return trace_black_curve_at(image, x_axis, y_axis, (c0, r0, c1, r1), xs, exclude_red=False, tick_pad=28)
