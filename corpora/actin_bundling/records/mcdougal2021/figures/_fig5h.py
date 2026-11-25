"""Shared calibration and SD error-bar reader for Fig. 5H (p. 5; identical to SI Fig. S8D p. 20)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

from _calib import find_box, find_ticks_mask, interval_from_px


def calibrate_fig5h(image: np.ndarray):
    """Calibrate Fig. 5H x-axis (34..42% dev), purple height y-axis (0..150 nm), and teal spacing y-axis (0..2+ um)."""
    c0, r0, c1, r1 = find_box(image, dark_thr=150.0, neutral_only=False)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    gray = (r + g + b) / 3.0

    # Bottom x-ticks (black/dark): 36, 38, 40 plus left/right frame borders 34, 42
    dark = gray < 150.0
    xt = find_ticks_mask(dark, r1, c0 + 25, c1 - 25, "x", -1, offset_lo=10, offset_hi=22)
    if len(xt) != 3:
        raise RuntimeError(f"Fig. 5H expected 3 interior x-ticks (36,38,40), found {len(xt)}: {xt}")
    x_axis = Axis("linear", [(c0, 34.0)] + list(zip(xt, [36.0, 38.0, 40.0])) + [(c1, 42.0)])

    # Purple left y-ticks: 150, 100, 50, 0 nm
    purple = (b - g > 25.0) & (r - g > 15.0) & (g < 150.0)
    pt = find_ticks_mask(purple, c0, r0 + 10, r1 - 10, "y", +1, offset_lo=10, offset_hi=22)
    if len(pt) != 4:
        raise RuntimeError(f"Fig. 5H expected 4 purple y-ticks (150,100,50,0), found {len(pt)}: {pt}")
    y_height = Axis("linear", list(zip(pt, [150.0, 100.0, 50.0, 0.0])))

    # Teal left y-tick at 2 um + bottom border at 0 um (ignore the dashed separator near purple 0)
    teal = (g - r > 25.0) & (b - r > 25.0) & (r < 120.0)
    r_sep = float(pt[-1])
    tt = find_ticks_mask(teal, c0, r_sep + 30.0, r1 - 15.0, "y", +1, offset_lo=10, offset_hi=22)
    if len(tt) != 1:
        raise RuntimeError(f"Fig. 5H expected 1 teal interior y-tick (2 um), found {len(tt)}: {tt}")
    y_spacing = Axis("linear", [(tt[0], 2.0), (r1, 0.0)])

    return x_axis, y_height, y_spacing, purple, teal, (c0, r0, c1, r1), r_sep


def find_purple_stems(purple: np.ndarray, box: tuple[float, float, float, float], r_sep: float) -> list[float]:
    """Locate the 10 error-bar vertical stem columns from the purple height panel."""
    c0, r0, c1, _r1 = box
    roi = np.zeros_like(purple)
    roi[int(r0) + 60:int(r_sep) - 5, int(c0) + 30:int(c1) - 30] = (
        purple[int(r0) + 60:int(r_sep) - 5, int(c0) + 30:int(c1) - 30]
    )
    stems = ndimage.binary_opening(roi, structure=np.ones((45, 1), bool))
    idx = np.nonzero(stems.any(axis=0))[0]
    if idx.size == 0:
        raise RuntimeError("Fig. 5H: no purple error-bar stems found")
    sp = np.nonzero(np.diff(idx) > 10)[0]
    cols = [float(np.mean(idx[s:e + 1]) + 0.5) for s, e in zip(np.r_[0, sp + 1], np.r_[sp, idx.size - 1])]
    if len(cols) != 10:
        raise RuntimeError(f"Fig. 5H expected 10 stems, found {len(cols)}: {cols}")
    return cols


def read_bars_from_stems(mask: np.ndarray, x_axis: Axis, y_axis: Axis,
                         stem_cols: list[float], r_min: int, r_max: int) -> list[dict]:
    """Read midpoint and exact top/bottom cap centers of each vertical SD bar at its stem column."""
    cal_err = max(y_axis.residual_px, 0.5)
    cap_half_err = 2.0 + cal_err  # conservative 95% reading precision on cap center (20-px thick vector cap)
    mid_half_err = 2.0 + cal_err

    out = []
    for sc in stem_cols:
        ci = int(round(sc))
        rows = np.nonzero(mask[r_min:r_max, ci - 2:ci + 3].any(axis=1))[0] + r_min
        if rows.size == 0:
            out.append({"x": float(x_axis.to_value(sc)), "y": None, "note": "error bar not found at stem"})
            continue
        r_top, r_bot = int(rows.min()), int(rows.max())
        # Measure exact center row of top horizontal cap and bottom horizontal cap (~20 px thick, ~121 px wide).
        # Restrict top cap to upper half (at most 21 px from r_top) and bottom cap to lower half (at most 21 px from r_bot)
        # so very short error bars (where caps are only a few px apart) never mix pixels between top and bottom caps.
        r_half = (r_top + r_bot) // 2
        top_end = min(r_top + 21, r_half + 1)
        bot_start = max(r_bot - 20, r_half)
        w_top = mask[r_top:top_end, ci - 25:ci + 26].sum(axis=1)
        w_bot = mask[bot_start:r_bot + 1, ci - 25:ci + 26].sum(axis=1)
        tr = np.nonzero(w_top > 25)[0] + r_top
        br = np.nonzero(w_bot > 25)[0] + bot_start
        if tr.size == 0 or br.size == 0:
            out.append({"x": float(x_axis.to_value(sc)), "y": None, "note": "error bar cap not found"})
            continue
        r_tc = float(tr.mean() + 0.5)
        r_bc = float(br.mean() + 0.5)
        r_mid = 0.5 * (r_tc + r_bc)

        xv = float(x_axis.to_value(sc))
        yv, lo, hi = interval_from_px(y_axis, r_mid, mid_half_err)
        bar_hi = float(y_axis.to_value(r_tc))
        bar_lo = float(y_axis.to_value(r_bc))
        bar_hi_outer = float(y_axis.to_value(r_tc - cap_half_err))
        bar_lo_outer = float(y_axis.to_value(r_bc + cap_half_err))
        out.append({
            "x": xv,
            "y": yv,
            "lo": lo,
            "hi": hi,
            "bar_lo": min(bar_lo, bar_hi),
            "bar_hi": max(bar_lo, bar_hi),
            "bar_lo_outer": min(bar_lo_outer, bar_hi_outer),
            "bar_hi_outer": max(bar_lo_outer, bar_hi_outer),
            "px": [sc, r_mid],
            "px_lo": r_tc,
            "px_hi": r_bc,
            "note": f"SD={(bar_hi - bar_lo) / 2.0:.3g}; caps ({tr.size},{br.size}) px; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px",
        })
    return out
