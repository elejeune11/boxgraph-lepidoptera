"""Calibrated reader for SI Fig. S10D (p. 23): 1D FFT spatial frequency spectrum magnitude (nm) vs frequency (1/um)."""
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

    # Top inward x-ticks (1, 2, 3, 4, 5) + left border (0) + right border (6)
    xt = find_ticks_mask(dark, r0, c0 + 15, c1 - 15, "x", +1, offset_lo=10, offset_hi=22)
    if len(xt) != 5:
        raise RuntimeError(f"Fig. S10D expected 5 interior top x-ticks (1..5), found {len(xt)}: {xt}")
    x_axis = Axis("linear", [(c0, 0.0)] + list(zip(xt, [1.0, 2.0, 3.0, 4.0, 5.0])) + [(c1, 6.0)])

    # Right inward y-ticks (4000, 3000, 2000, 1000 nm) + bottom border (0 nm); right border avoids DC curve
    yt = find_ticks_mask(dark, c1, r0 + 15, r1 - 15, "y", -1, offset_lo=10, offset_hi=22)
    if len(yt) != 4:
        raise RuntimeError(f"Fig. S10D expected 4 right y-ticks (4,3,2,1 x10^3 nm), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [4000.0, 3000.0, 2000.0, 1000.0])) + [(r1, 0.0)])

    # Mask out the printed annotation text ('freq = s^-1 = 0.684...', 's = 1.4618um') and the
    # vertical dashed pointer line above the dominant peak (f > 0.35 and mag > 1350 nm).
    clean = image.copy()
    c_cut = int(round(x_axis.to_pixel(0.35)))
    r_cut = int(round(y_axis.to_pixel(1350.0)))
    clean[int(r0) + 10:r_cut, c_cut:int(c1) - 10] = 255

    # Allow curve right down to bottom border (zero troughs touch r1; use small tick_pad near bottom)
    # Wait: bottom ticks at 1..5 extend ~22px up from r1, so exclude bottom-tick cols only when at integer f!
    xs = [float(v) for v in at] if at else [
        0.02, 0.05, 0.10, 0.15, 0.20, 0.24, 0.30, 0.39, 0.48, 0.55, 0.62,
        0.684, 0.74, 0.80, 0.85, 0.94, 1.04, 1.13, 1.24, 1.34, 1.44, 1.52,
        1.64, 1.75, 1.92, 2.05, 2.13, 2.24, 2.32, 2.41, 2.52, 2.69, 2.85,
        3.05, 3.35, 3.65, 3.95, 4.25, 4.65, 5.05, 5.45, 5.60,
    ]
    return trace_black_curve_at(clean, x_axis, y_axis, (c0, r0, c1, r1), xs,
                                exclude_red=True, y_min_row=int(r0) + 6, y_max_row=int(r1) - 4, tick_pad=8)
