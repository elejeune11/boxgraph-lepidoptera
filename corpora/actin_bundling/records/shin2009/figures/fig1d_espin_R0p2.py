"""Calibrated reader for Fig. 1(d) (p. 2): espin:actin SAXS intensity I(q) at R = 0.2 (black, upper panel)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import black_mask, calibrate_fig1d, runs_1d, track_black_curve, value_interval


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB), lr_diff_px = calibrate_fig1d(image, "upper")
    blk = black_mask(image, gray_max=85.0, sat_max=30.0)

    # Locate the two vertical dashed guide lines near q = 0.120 and q = 0.135
    excl = []
    for q_lo, q_hi in ((0.118, 0.122), (0.133, 0.137)):
        c0 = int(round(x_axis.to_pixel(q_lo)))
        c1 = int(round(x_axis.to_pixel(q_hi)))
        proj = blk[rT + 150 : rB - 50, c0:c1].sum(axis=0)
        c_pk = c0 + int(np.argmax(proj))
        excl.append((c_pk - 13, c_pk + 13))

    # Find curve seed right inside left border (I in [105, 145])
    r_seed_lo = int(round(y_axis.to_pixel(145.0)))
    r_seed_hi = int(round(y_axis.to_pixel(105.0)))
    c_start = cL + 15
    seed_rows = np.where(blk[r_seed_lo:r_seed_hi, c_start])[0] + r_seed_lo
    if len(seed_rows) == 0:
        raise RuntimeError("Fig. 1(d) R=0.2: could not locate black curve seed near left border")
    r_seed = int(round(runs_1d(seed_rows)[0][2]))

    tracked = track_black_curve(
        blk, c_start=c_start, c_end=cR - 15, r_seed=r_seed,
        r_min=rT + 10, r_max=int(round(y_axis.to_pixel(75.0))),
        exclude_cols=excl, max_step=30, max_thick=55,
    )

    cal_y = max(y_axis.residual_px, 0.5) + 0.5 * lr_diff_px
    out = []
    for xv in at:
        c_f = float(x_axis.to_pixel(xv))
        c0 = int(round(c_f))
        if c0 in tracked:
            r_c, thick = tracked[c0]
            h_l = tracked.get(c0 - 3)
            h_r = tracked.get(c0 + 3)
            slope_px = abs(h_r[0] - h_l[0]) / 6.0 if (h_l and h_r) else 0.0
            half_px = 0.5 * thick + 1.5 + cal_y + 0.5 * slope_px
            extra = ""
        else:
            cl = next((c0 - d for d in range(1, 45) if (c0 - d) in tracked), None)
            cr = next((c0 + d for d in range(1, 45) if (c0 + d) in tracked), None)
            if cl is None or cr is None:
                out.append({"x": float(xv), "y": None, "note": f"fig1d R=0.2 black: no tracked pixel near col {c0}"})
                continue
            rl, tl = tracked[cl]
            rr, tr = tracked[cr]
            t = (c_f - cl) / (cr - cl)
            r_c = rl + t * (rr - rl)
            thick = 0.5 * (tl + tr)
            half_px = 0.5 * thick + 2.5 + cal_y + 0.25 * abs(rr - rl)
            extra = f"; bridged [{cl},{cr}]"
        y_val, y_lo, y_hi = value_interval(y_axis, r_c, half_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "px": [c_f, r_c],
            "px_lo": r_c - half_px,
            "px_hi": r_c + half_px,
            "note": f"fig1d espin R=0.2 black (upper); thick={thick:.0f}px, lr_tick_diff={lr_diff_px:.1f}px{extra}",
        })
    return out
