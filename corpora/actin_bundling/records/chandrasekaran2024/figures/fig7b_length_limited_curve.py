"""Calibrated reader for the solid black theoretical curve ('Filament length-limited aspect ratio') in Fig. 7B (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import math
import numpy as np
from scipy import ndimage
import _plot_common as pc

FRACS = [2 / 3, 5 / 6, 11 / 12, 1.0, 13 / 12, 7 / 6, 5 / 4, 4 / 3, 5 / 3, 2.0]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [f * math.pi for f in FRACS])
    y_axis = pc.fit_axis(info["y_ticks"], [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    c0, r0, c1, r1 = info["box"]

    # Pure black pixels inside the plot box
    dark = (image.max(axis=2) < 60) & ((image.max(axis=2) - image.min(axis=2)) < 25)
    mask = np.zeros_like(dark)
    mask[r0 + 4:r1 - 4, c0 + 4:c1 - 4] = dark[r0 + 4:r1 - 4, c0 + 4:c1 - 4]

    # Isolate only long horizontal-span connected components (the theoretical curve spans >1500 px total,
    # broken only where the red error bar at 5/4 crosses it; dots and text letters are < 70 px wide).
    lab, n = ndimage.label(mask)
    curve_mask = np.zeros_like(mask)
    for k in range(1, n + 1):
        cols = np.flatnonzero((lab == k).any(axis=0))
        if cols.size and (cols[-1] - cols[0]) > 250:
            curve_mask |= (lab == k)

    # For any column where a replicate dot touches the curve vertically, trim vertical runs to the
    # expected ~22-25 px stroke thickness via local median across +-12 columns.
    readings = []
    for xv in at:
        pc_x = x_axis.to_pixel(xv)
        col = int(round(pc_x))
        # Find local clean curve columns (single vertical run of thickness 14..34 px)
        neigh_c, neigh_h = [], []
        for dc in range(-45, 46, 3):
            c2 = col + dc
            if c0 < c2 < c1:
                rs2 = np.flatnonzero(curve_mask[:, c2])
                if 14 <= rs2.size <= 34 and (rs2[-1] - rs2[0] + 1) <= 36:
                    neigh_c.append((c2, float(np.median(rs2)) + 0.5))
                    neigh_h.append(0.5 * (float(rs2.max()) - float(rs2.min()) + 1.0))
        rs = np.flatnonzero(curve_mask[:, col])
        if rs.size and 14 <= rs.size <= 34 and (rs[-1] - rs[0] + 1) <= 36:
            center = 0.5 * (float(rs.min()) + float(rs.max())) + 0.5
            half = 0.5 * (float(rs.max()) - float(rs.min()) + 1.0)
            amb = False
        elif len(neigh_c) >= 2:
            xs_n = np.array([p[0] for p in neigh_c], dtype=float)
            ys_n = np.array([p[1] for p in neigh_c], dtype=float)
            m_sl, b_ic = np.polyfit(xs_n, ys_n, 1)
            pred_r = float(m_sl * pc_x + b_ic)
            # Check if there is a clean separate vertical run at `col` close to pred_r
            matched = None
            if rs.size:
                sp = np.flatnonzero(np.diff(rs) > 2)
                st = np.concatenate(([0], sp + 1))
                en = np.concatenate((sp, [rs.size - 1]))
                for a2, b2 in zip(st, en):
                    rc = 0.5 * (float(rs[a2]) + float(rs[b2])) + 0.5
                    rh = 0.5 * (float(rs[b2]) - float(rs[a2]) + 1.0)
                    if 8.0 <= rh <= 18.0 and abs(rc - pred_r) < 12.0:
                        matched = (rc, rh)
                        break
            if matched is not None:
                center, half = matched
                amb = False
            else:
                center = pred_r
                half = float(np.median(neigh_h)) + 2.0
                amb = True
        else:
            readings.append({"x": float(xv), "y": None, "note": f"no curve pixels near x={xv:g}"})
            continue
        slope = 0.0
        if len(neigh_c) >= 2:
            xs_n = np.array([p[0] for p in neigh_c], dtype=float)
            ys_n = np.array([p[1] for p in neigh_c], dtype=float)
            slope = float(np.polyfit(xs_n, ys_n, 1)[0])
        eff_half = half + abs(slope) * abs(pc_x - (col + 0.5))
        note = f"Fig. 7B length-limited curve (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        if amb:
            note += "; interpolated across overlapping dot/error bar"
        readings.append(pc.make_reading(xv, center, y_axis, pc_x, eff_half, x_axis=x_axis, slope_px=slope, note=note))
    return readings
