"""Calibrated reader for Fig. 4E (bundle contractile strain box-and-whisker plot vs filamin %).

Note: The published y-axis misprints the 0.3 tick label as '0.6' (reading 0.2, 0.6, 0.4, 0.5, 0.6
from bottom to top) and shifts that tick ~8 px upward. We calibrate from the four unambiguous
ticks [0.6, 0.5, 0.4, 0.2] (residual < 1.0 px) and extract each box's white median line, Q1-Q3
gray box edges, and black whisker endpoints.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    if len(xt) != 3:
        raise RuntimeError(f"Fig4E: expected 3 categorical x ticks (2.5, 5.0, 10.0%), got {xt}")
    if len(yt) != 5:
        raise RuntimeError(f"Fig4E: expected 5 y ticks (0.6, 0.5, 0.4, [misprinted 0.3], 0.2), got {yt}")
    # Exclude index 3 (the misprinted '0.6' tick at ~0.3, shifted ~8 px high) from linear axis fit
    y_marks = [(yt[0], 0.6), (yt[1], 0.5), (yt[2], 0.4), (yt[4], 0.2)]
    ya = C.verify_axis("linear", y_marks, max_residual_px=1.5)
    return gray, xt, ya, int(ht[2]) + 25, int(hb[1]) - 25


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray, xt, ya, r0, r1 = calibrate(image)
    fil_vals = [2.5, 5.0, 10.0]
    out = []
    for c_tick, fil in zip(xt, fil_vals):
        ci = int(round(c_tick))
        # Gray box region spans ~+-80 px horizontally around x-tick; sample core +-50 px
        box_strip = (gray[r0:r1, ci - 50 : ci + 51] > 140.0) & (gray[r0:r1, ci - 50 : ci + 51] < 210.0)
        row_frac_gray = box_strip.mean(axis=1)
        gray_rows = np.where(row_frac_gray > 0.7)[0] + r0
        if len(gray_rows) < 20:
            out.append({"x": fil, "y": None, "found": False, "note": f"no gray box at fil={fil}%"})
            continue
        r_q3 = int(gray_rows[0])
        r_q1 = int(gray_rows[-1])

        # White horizontal median line inside the gray box
        white_strip = gray[r_q3 + 6 : r_q1 - 5, ci - 50 : ci + 51] > 230.0
        white_rows = np.where(white_strip.mean(axis=1) > 0.7)[0] + (r_q3 + 6)
        if len(white_rows) == 0:
            out.append({"x": fil, "y": None, "found": False, "note": f"no white median line at fil={fil}%"})
            continue
        r_med = float(white_rows.mean()) + 0.5

        # Black vertical whiskers above r_q3 and below r_q1
        blk = gray < 70.0
        bar_lo, bar_hi, bar_lo_out, bar_hi_out, row_hi, row_lo = C.trace_dashed_error_bar(
            blk, c_tick, r_q3 - 2, r_q1 + 2, r0, r1, ya, stem_half_w=8, gap=10
        )

        q1_val = float(ya.to_value(r_q1 + 0.5))
        q3_val = float(ya.to_value(r_q3 + 0.5))
        med_val = float(ya.to_value(r_med))
        half_px = max(0.5 * (white_rows[-1] - white_rows[0] + 1), 2.0) + ya.residual_px
        ya_, yb_ = float(ya.to_value(r_med - half_px)), float(ya.to_value(r_med + half_px))
        out.append(
            {
                "x": fil,
                "y": med_val,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "bar_lo": bar_lo,
                "bar_hi": bar_hi,
                "bar_lo_outer": bar_lo_out,
                "bar_hi_outer": bar_hi_out,
                "px": [float(c_tick), r_med],
                "px_lo": float(row_hi),
                "px_hi": float(row_lo),
                "note": (
                    f"Fig4E box median={med_val:.3f}, box=[{q1_val:.3f},{q3_val:.3f}], "
                    f"whiskers=[{bar_lo:.3f},{bar_hi:.3f}], cal_res_y={ya.residual_px:.2f}px"
                ),
            }
        )
    return out
