"""Calibrated reader for Fig. 4I (p. 7): barbed-end dwell time (open circles, left axis) with SEM error bars."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya_l, _, fr = C.calibrate_fig4i(image)
    g = C.gray(image)
    rT, rB = int(round(fr["rT"])), int(round(fr["rB"]))
    cal_err = max(ya_l.residual_px, 0.5)
    concs = [0.5, 0.75, 1.0, 1.5, 2.0]
    out = []
    for xv in concs:
        c = xa.to_pixel(xv)
        ci = int(round(c))
        # Open circles all lie in dwell-time band [1.10, 1.72] s
        r0, r1 = int(ya_l.to_pixel(1.70)), int(ya_l.to_pixel(1.10))
        sub = g[r0:r1, ci - 13 : ci + 14] < 175
        runs = C.runs_1d(sub.any(axis=1))
        if not runs:
            raise ValueError(f"Fig. 4I dwell time: no marker found at [Actin]={xv} uM")
        # At x = 1.0, the bottom error-bar cap of the filled circle may intrude at top of band (y ~ 1.70); take lowest run!
        s, e = runs[-1]
        # Within [s, e], marker body has widest horizontal extent (width >= 13 px) while error-bar stems are thinner
        row_w = sub[s : e + 1].sum(axis=1)
        body_rows = np.nonzero(row_w >= max(8, 0.65 * row_w.max()))[0]
        r_mid = r0 + s + 0.5 * (body_rows[0] + body_rows[-1]) + 0.5
        r_top_err = r0 + float(s) + 0.5
        r_bot_err = r0 + float(e) + 0.5
        y, lo, hi = C.val_interval(ya_l, r_mid, 1.5 + cal_err)
        err_px = 1.5 + cal_err
        bar_hi = float(ya_l.to_value(r_top_err))
        bar_lo = float(ya_l.to_value(r_bot_err))
        bar_hi_out = float(ya_l.to_value(r_top_err - err_px))
        bar_lo_out = float(ya_l.to_value(r_bot_err + err_px))
        note = f"open circle [{r_top_err:.1f}, {r_bot_err:.1f}]px"
        if abs(xv - 0.75) < 0.01:
            note += " (overlaps filled monomer marker at 0.75 uM)"
        out.append({
            "x": xv, "y": y, "lo": lo, "hi": hi,
            "bar_lo": min(bar_lo, bar_hi), "bar_hi": max(bar_lo, bar_hi),
            "bar_lo_outer": min(bar_lo_out, bar_hi_out), "bar_hi_outer": max(bar_lo_out, bar_hi_out),
            "px": [c, r_mid], "px_lo": r_top_err, "px_hi": r_bot_err,
            "note": note,
        })
    return out
