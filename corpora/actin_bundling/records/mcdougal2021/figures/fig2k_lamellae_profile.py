"""Calibrated reader for Fig. 2K (p. 3): 1D unwrapped ridge-lamellae surface height profile h(x) at 83% dev."""
from __future__ import annotations

import sys
from pathlib import Path
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_box, find_ticks_mask, neutral_dark_mask, trace_black_curve_at


def read(image, at):
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    xt = find_ticks_mask(dark, r1, c0 + 10, c1 - 10, "x", +1, offset_lo=10, offset_hi=24)
    if len(xt) != 10:
        raise RuntimeError(f"Fig. 2K expected 10 x-ticks (0..9), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [float(i) for i in range(10)])))

    # Left outward y-ticks excluding bottom border r1 (no tick is drawn at 0 in Fig. 2K): 300, 200, 100 nm
    yt = find_ticks_mask(dark, c0, r0 - 15, r1 - 25, "y", -1, offset_lo=10, offset_hi=24)
    if len(yt) != 3:
        raise RuntimeError(f"Fig. 2K expected 3 y-ticks (300,200,100), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [300.0, 200.0, 100.0])))

    # Mask out the panel letter 'K' in the upper-left interior
    clean = image.copy()
    clean[int(r0) + 10:int(r0) + 360, int(c0) + 10:int(c0) + 270] = 255

    xs = [float(v) for v in at] if at else [round(0.1 * i, 4) for i in range(1, 91)]
    return trace_black_curve_at(clean, x_axis, y_axis, (c0, r0, c1, r1), xs, exclude_red=False, tick_pad=26)
