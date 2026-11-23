"""Calibrated reader for Main Fig. 3B left (p. 4): Wild-type ring spacing probability distribution."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_fig3b_dist(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    # Exclude bottom figure panel border if present
    h_box = [h for h in hl if h[0] < image.shape[0] - 80]
    h_top, h_bot = h_box[0], h_box[-1]
    v_left, v_right = vl[0], vl[-1]

    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 14), margin=12, max_width=10)
    # 7 regular interior x-ticks at 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4 (left border = 0.0, right border = 1.5)
    xt_int = match_ticks_regular(xt_top, 7, tol=3.0)
    x_marks = [(v_left[0], 0.0)] + list(zip(xt_int, [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])) + [(v_right[0], 1.5)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=1.5)

    # 8 interior y-ticks at 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, bottom border = 0.00
    yt_int = match_ticks_regular([t for t in yt_left if t < h_bot[0] - 25], 8, tol=2.5)
    y_marks = list(zip(yt_int, [0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05])) + [(h_bot[0], 0.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=1.5)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read_circle_series(mask: np.ndarray, x_axis, y_axis, box, at):
    cl, rt, cr, rb = box
    # Exclude legend in top right
    m = mask.copy()
    m[: int(rt) + 320, int(cl) + int(0.55 * (cr - cl)) :] = False
    m[: int(rt) - 6, :] = False
    m[int(rb) + 28 :, :] = False

    half_px = 4.0 + y_axis.residual_px
    out = []
    for x_val in at:
        c_f = float(x_axis.to_pixel(x_val))
        c0 = max(int(cl) + 2, int(round(c_f)) - 22)
        c1 = min(int(cr) - 2, int(round(c_f)) + 22)
        # Look for circular marker pixels (height ~ 38 px) centered near c_f
        sub = m[:, c0 : c1 + 1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 40:
            out.append({"x": float(x_val), "y": None, "found": False, "note": f"no circle ink near x={x_val}"})
            continue
        # At steep connecting lines, restrict to narrow column band +/- 6 px around c_f to find marker center
        narrow = m[:, max(0, int(round(c_f)) - 4) : int(round(c_f)) + 5]
        nrr = np.where(narrow.any(axis=1))[0]
        if len(nrr) == 0:
            r_mid = float(np.median(rr)) + 0.5
        else:
            # Circle top & bottom edges span ~36-44 px; use midpoint of outer vertical span of circle near c_f
            # Check vertical runs in narrow band; if multiple (circle + line or baseline), pick cluster around circle
            # At y=0 baseline, circle is centered right on rb
            if abs(nrr.max() - rb) < 25 and (nrr.max() - nrr.min()) < 55:
                r_mid = 0.5 * (float(nrr.min()) + float(nrr.max())) + 0.5
            else:
                # Find the dense 38-px circle ring (vertical span where marker ink spans width >= 28 px)
                row_w = np.array([
                    (np.nonzero(sub[r])[0].max() - np.nonzero(sub[r])[0].min() + 1) if sub[r].any() else 0
                    for r in range(int(nrr.min()), int(nrr.max()) + 1)
                ])
                wide_rows = np.where(row_w >= 26)[0] + int(nrr.min())
                if len(wide_rows) >= 5:
                    r_mid = 0.5 * (float(wide_rows.min()) + float(wide_rows.max())) + 0.5
                else:
                    r_mid = 0.5 * (float(nrr.min()) + float(nrr.max())) + 0.5
        y_val, lo, hi = value_interval(y_axis, r_mid, half_px)
        # Clamp tiny negative quantization jitter at baseline to 0.0 within interval
        if abs(y_val) < 0.002 and lo <= 0.0 <= hi:
            y_val = 0.0
        out.append({
            "x": float(x_val),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_f, r_mid],
            "px_lo": float(r_mid - half_px),
            "px_hi": float(r_mid + half_px),
            "note": f"circle center at x={x_val}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return out


def read(image, at):
    x_axis, y_axis, box = calibrate_fig3b_dist(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 80) & (rgb[:, :, 2] < 80)
    return read_circle_series(red, x_axis, y_axis, box, at)
