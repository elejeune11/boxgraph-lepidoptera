"""Calibrated reader for Fig. 2J (p. 3): 1D unwrapped ridge surface height profile h(x) at 83% dev."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_box, find_ticks_mask, neutral_dark_mask, trace_black_curve_at


def read(image, at):
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    # Bottom x-ticks extend outward below r1: 0, 1, ..., 9 um
    xt = find_ticks_mask(dark, r1, c0 + 10, c1 - 10, "x", +1, offset_lo=10, offset_hi=24)
    if len(xt) != 10:
        raise RuntimeError(f"Fig. 2J expected 10 x-ticks (0..9), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [float(i) for i in range(10)])))

    # Left inward y-ticks: 200, 100, 0 nm (avoids top-border stroke merge at 300 nm)
    yt = find_ticks_mask(dark, c0, r0 + 15, r1 - 15, "y", +1, offset_lo=10, offset_hi=24)
    if len(yt) != 3:
        raise RuntimeError(f"Fig. 2J expected 3 inward y-ticks (200,100,0), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [200.0, 100.0, 0.0])))

    # Mask out the panel letter 'J' in the upper-left interior so it cannot confuse tracing near x < 0.7
    clean = image.copy()
    clean[int(r0) + 10:int(r0) + 330, int(c0) + 10:int(c0) + 240] = 255

    xs = [float(v) for v in at] if at else [round(0.1 * i, 4) for i in range(1, 91)]
    return trace_black_curve_at(clean, x_axis, y_axis, (c0, r0, c1, r1), xs, exclude_red=False, tick_pad=26)
