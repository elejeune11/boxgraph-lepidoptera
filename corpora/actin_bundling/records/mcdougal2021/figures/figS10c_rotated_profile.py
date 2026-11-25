"""Calibrated reader for SI Fig. S10C (p. 23): slope-rotated 1D surface height profile h(x) at 41% dev."""
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

    xt = find_ticks_mask(dark, r1, c0 + 15, c1 - 15, "x", -1, offset_lo=10, offset_hi=24)
    if len(xt) != 8:
        raise RuntimeError(f"Fig. S10C expected 8 x-ticks (0..7), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [float(i) for i in range(8)])))

    yt = find_ticks_mask(dark, c0, r0 + 15, r1 - 15, "y", +1, offset_lo=10, offset_hi=24)
    if len(yt) != 3:
        raise RuntimeError(f"Fig. S10C expected 3 y-ticks (200,100,0), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [200.0, 100.0, 0.0])))

    # 79 uniform samples across [0.00, 6.905] um (dx = 88.5 nm, matching S10D Nyquist 5.65 um^-1);
    # exclude_red=True strips the red 'peak height' arrow & dashed baseline cleanly.
    xs = [float(v) for v in at] if at else [float(v) for v in np.linspace(0.0, 6.905, 79)]
    return trace_black_curve_at(image, x_axis, y_axis, (c0, r0, c1, r1), xs, exclude_red=True, tick_pad=30)
