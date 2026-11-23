"""Calibrated reader for SI Fig. S14 (p. 22): steady-state numerical ring wavelength probability density histogram."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def read(image, at):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 22), margin=15, max_width=14)

    # X ticks: left = 0.5, center tick = 1.0, right = 1.5
    if not xt_top:
        raise RuntimeError("Fig. S14: could not locate center x-tick at 1.0")
    x_mid = min(xt_top, key=lambda t: abs(t - 0.5 * (v_left[0] + v_right[0])))
    x_axis = verify_axis("linear", [(v_left[0], 0.5), (x_mid, 1.0), (v_right[0], 1.5)], max_residual_px=3.0)

    # Y ticks: top = 14, 6 interior ticks = 12, 10, 8, 6, 4, 2 (exclude bar-thickened bottom border)
    yt_int = match_ticks_regular(yt_left, 6, tol=3.0)
    y_axis = verify_axis(
        "linear",
        [(h_top[0], 14.0)] + list(zip(yt_int, [12.0, 10.0, 8.0, 6.0, 4.0, 2.0])),
        max_residual_px=1.5,
    )

    rgb = image.astype(float)
    # Dark navy blue MATLAB histogram bar fill: B > 90, R < 40, G < 40
    bar_fill = (rgb[:, :, 2] > 90.0) & (rgb[:, :, 0] < 40.0) & (rgb[:, :, 1] < 40.0)
    r_lo, r_hi = int(h_top[2]) + 8, int(h_bot[1]) - 2
    c_lo, c_hi = int(v_left[2]) + 4, int(v_right[1]) - 4

    col_bars = clusters(np.where(bar_fill[r_lo:r_hi, c_lo:c_hi].any(axis=0))[0] + c_lo, gap=2)
    half_px = 3.0 + y_axis.residual_px

    out = []
    for c_mid, c0, c1 in col_bars:
        if (c1 - c0 + 1) < 15:
            continue
        # Top edge of the navy blue bar interior (excluding thin black outline border, so top border center is ~1.5 px above)
        interior_cols = range(c0 + 4, c1 - 3)
        tops = []
        for cc in interior_cols:
            rs = np.where(bar_fill[r_lo:r_hi, cc])[0] + r_lo
            if len(rs) > 0:
                tops.append(float(rs.min()))
        if not tops:
            continue
        r_top_bar = float(np.median(tops)) - 1.5
        x_val = float(x_axis.to_value(c_mid))
        y_val, lo, hi = value_interval(y_axis, r_top_bar, half_px)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_mid, r_top_bar],
            "px_lo": float(r_top_bar - half_px),
            "px_hi": float(r_top_bar + half_px),
            "note": f"histogram bar width={c1-c0+1}px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
