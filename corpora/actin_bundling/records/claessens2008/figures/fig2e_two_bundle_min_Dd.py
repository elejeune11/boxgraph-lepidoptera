"""Calibrated reader for Fig. 2E (p. 2): lower droplet-diameter onset of the two-bundle phase (leftmost green circle per c_a)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2e import calibrate_fig2e, circle_radius_px
from _calib import value_interval


def read(image, at):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig2e(image)
    rgb = image.astype(float)
    # Exact green circle RGB in Fig. 2E is ~(120, 200, 50); excludes darkened yellow/black overlaps
    dist_green = np.sqrt(((rgb - np.array([120.0, 200.0, 50.0])) ** 2).sum(axis=2))
    m_green = dist_green <= 40.0
    rad = circle_radius_px(image, y_axis, c_right)

    ca_levels = [0.02, 0.03, 0.04, 0.05, 0.10]
    readings = []
    for ca_nom in ca_levels:
        r_nom = int(round(y_axis.to_pixel(ca_nom)))
        sl_r = slice(r_nom - 40, r_nom + 40)
        sl_c = slice(int(c_left) + 15, int(c_right) - 15)
        sub = m_green[sl_r, sl_c]
        cols = np.nonzero(sub.sum(axis=0) >= 4)[0]
        if len(cols) == 0:
            continue
        c_leftmost_edge = int(cols.min()) + int(c_left) + 15
        c_center = float(c_leftmost_edge + rad + 0.5)
        arc_sub = m_green[sl_r, c_leftmost_edge : min(image.shape[1], c_leftmost_edge + 14)]
        rr_arc = np.nonzero(arc_sub)[0]
        r_center = float(rr_arc.mean() + 0.5 + (r_nom - 40)) if len(rr_arc) else float(r_nom)

        dd_val = float(x_axis.to_value(c_center))
        half_px = 4.0 + y_axis.residual_px
        ca_val, ca_lo, ca_hi = value_interval(y_axis, r_center, half_px)
        readings.append({
            "x": dd_val,
            "y": ca_val,
            "lo": ca_lo,
            "hi": ca_hi,
            "px": [c_center, r_center],
            "px_lo": r_center - half_px,
            "px_hi": r_center + half_px,
            "note": f"leftmost green two-bundle circle at nominal c_a={ca_nom:g} mg/ml (rad={rad:.1f} px, c_a*D_d^2={ca_val*dd_val*dd_val:.2f}); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })

    readings.sort(key=lambda r: r["x"])
    return readings
