"""Calibrated reader for Main Fig. 2E (histogram of myosin-to-droplet alignment angle theta_MD)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    find_box_spines,
    find_ticks,
    match_ticks_regular,
    runs_1d,
    verify_axis,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # Top spine has 9 ticks at 0, 10, 20, 30, 40, 50, 60, 70, 80 deg; right spine is 90 deg.
    xt_raw = find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5)
    xt = match_ticks_regular(xt_raw, 9, tol=5.0)
    x_marks = [(xt[i], float(10 * i)) for i in range(9)] + [(right[0], 90.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Right spine has 4 ticks at 0.4, 0.3, 0.2, 0.1; top = 0.5, bot = 0.0.
    yt_raw = find_ticks(dark, right, (top[2] + 10, bot[1] - 10), "y", side=-1, band=(5, 25), frac=0.5)
    yt = match_ticks_regular(yt_raw, 4, tol=5.0)
    y_marks = [(top[0], 0.5)] + [(yt[i], 0.4 - 0.1 * i) for i in range(4)] + [(bot[0], 0.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    # Detect solid black histogram bars rising from the bottom spine (width >= 40 px excludes 9-px ticks).
    band = dark[bot[1] - 25 : bot[1] - 5, left[2] + 2 : right[1] - 2].mean(axis=0) > 0.85
    bar_runs = [(c0 + left[2] + 2, c1 + left[2] + 2) for c0, c1 in runs_1d(band) if (c1 - c0 + 1) >= 40]
    if len(bar_runs) != 4:
        raise RuntimeError(f"Fig. 2E: expected 4 histogram bars, found {len(bar_runs)}: {bar_runs}")

    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
    out = []
    for c0, c1 in bar_runs:
        tops = []
        for c in range(c0 + 8, c1 - 7):
            r = bot[1] - 2
            while r > top[2] + 2 and dark[r - 1, c]:
                r -= 1
            tops.append(float(r))
        r_top = float(np.median(tops))
        c_mid = 0.5 * (c0 + c1) + 0.5
        half_px = 2.0 + y_axis.residual_px
        y_val = float(y_axis.to_value(r_top))
        y_a = float(y_axis.to_value(r_top - half_px))
        y_b = float(y_axis.to_value(r_top + half_px))
        out.append(
            {
                "x": float(x_axis.to_value(c_mid)),
                "y": y_val,
                "lo": min(y_a, y_b),
                "hi": max(y_a, y_b),
                "px": [c_mid, r_top],
                "px_lo": r_top - half_px,
                "px_hi": r_top + half_px,
                "note": f"bar w={c1 - c0 + 1}px, {cal_note}",
            }
        )
    return out
