"""Calibrated reader for Main Fig. 1A (p. 2): relative cortical actin intensity (%) vs time (s) with error bars."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, verify_axis, value_interval


def calibrate_fig1a(image: np.ndarray):
    gray = image.astype(float).mean(axis=2)
    H, W = gray.shape

    # Left gray vertical y-axis spine (gray < 210 spanning > 45% of crop height)
    v_cols = clusters(np.where((gray[:, : W // 2] < 210).mean(axis=0) > 0.45)[0])
    if not v_cols:
        raise RuntimeError("Fig. 1A: could not locate y-axis spine")
    c_left = v_cols[0][0]

    # Bottom gray x-axis line and 4 clean horizontal gray gridlines measured on the right (no data there)
    c_r0 = int(round(c_left)) + int(0.75 * (W - c_left))
    c_r1 = min(W - 40, int(round(c_left)) + int(0.95 * (W - c_left)))
    row_min_r = gray[:, c_r0:c_r1].mean(axis=1)
    g_rows = clusters(np.where(row_min_r[: int(0.92 * H)] < 248.0)[0], gap=4)
    # Expect top (170%), 136%, 68%, 34%, bottom (0%) — the 102% gridline is fainter (~250)
    # Match the 5 evenly spaced intervals from top (95.5) to bottom (1287.0)
    if len(g_rows) < 5:
        raise RuntimeError(f"Fig. 1A: expected >=5 horizontal gridlines, got {g_rows}")
    r_top = g_rows[0][0]
    r_bot = g_rows[-1][0]
    step = (r_bot - r_top) / 5.0
    marks_y = [(r_top, 170.0), (r_bot, 0.0)]
    for k in (1, 2, 3, 4):
        target = r_top + k * step
        near = min(g_rows, key=lambda g: abs(g[0] - target))
        if abs(near[0] - target) < 8.0:
            marks_y.append((near[0], 170.0 * (5 - k) / 5.0))
    y_axis = verify_axis("linear", marks_y, max_residual_px=4.0)

    # Right extent of x-axis line at r_bot -> t = 5000 s (verified against 72.6 min micrograph endpoint)
    rb = int(round(r_bot))
    x_line_cols = np.where(gray[rb - 4 : rb + 5, int(c_left) : W - 30].min(axis=0) < 235.0)[0] + int(c_left)
    c_right = float(x_line_cols.max()) + 0.5
    x_axis = verify_axis("linear", [(c_left, 0.0), (c_right, 5000.0)], max_residual_px=2.0)
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def read(image, at):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig1a(image)
    gray = image.astype(float).mean(axis=2)
    black = gray < 95.0
    bar_ink = (gray >= 95.0) & (gray < 222.0)

    out = []
    for x_val in at:
        c_f = float(x_axis.to_pixel(x_val))
        c = int(round(c_f))
        # Local black curve mean row in +/- 3 px window
        rr_b = np.where(black[int(r_top) - 30 : int(r_bot) - 30, max(int(c_left), c - 3) : c + 4])[0] + (int(r_top) - 30)
        if len(rr_b) == 0:
            out.append({"x": float(x_val), "y": None, "found": False, "note": f"no curve ink near t={x_val}"})
            continue
        r_mid = float(np.median(rr_b)) + 0.5
        thick = float(np.percentile(rr_b, 90) - np.percentile(rr_b, 10) + 1.0)
        half_px = max(0.5 * thick, 2.0) + 1.5 + y_axis.residual_px
        y_val, lo, hi = value_interval(y_axis, r_mid, half_px)

        # Find nearest vertical gray error bar within +/- 18 px
        best_bar = None
        for dc in range(0, 20):
            for sgn in ((0,) if dc == 0 else (-1, 1)):
                cc = c + sgn * dc
                if cc <= int(c_left) + 2 or cc >= int(c_right) - 2:
                    continue
                # Check vertical stem above and below r_mid
                r_lo_s = max(int(r_top) - 40, int(r_mid) - 180)
                r_hi_s = min(int(r_bot) - 20, int(r_mid) + 180)
                col_bar = bar_ink[r_lo_s:r_hi_s, cc - 1 : cc + 2].any(axis=1)
                idx = np.where(col_bar)[0] + r_lo_s
                if len(idx) >= 25 and idx.min() < r_mid - 15 and idx.max() > r_mid + 15:
                    # Verify it is contiguous-ish vertical bar (not stray gridline)
                    best_bar = (float(idx.min()) + 0.5, float(idx.max()) + 0.5, cc)
                    break
            if best_bar is not None:
                break

        rec = {
            "x": float(x_val),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_f, r_mid],
            "px_lo": float(r_mid - half_px),
            "px_hi": float(r_mid + half_px),
            "note": f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        }
        if best_bar is not None:
            r_up, r_dn, cc_bar = best_bar
            err_px = 3.0 + y_axis.residual_px
            b_hi = float(y_axis.to_value(r_up))
            b_lo = float(y_axis.to_value(r_dn))
            b_hi_o = float(y_axis.to_value(r_up - err_px))
            b_lo_o = float(y_axis.to_value(r_dn + err_px))
            rec.update({
                "bar_lo": min(b_lo, b_hi),
                "bar_hi": max(b_lo, b_hi),
                "bar_lo_outer": min(b_lo_o, b_hi_o),
                "bar_hi_outer": max(b_lo_o, b_hi_o),
                "px_lo": r_up,
                "px_hi": r_dn,
                "note": f"error bar at col {cc_bar} (rows {r_up:.0f}..{r_dn:.0f}); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            })
        out.append(rec)
    return out
