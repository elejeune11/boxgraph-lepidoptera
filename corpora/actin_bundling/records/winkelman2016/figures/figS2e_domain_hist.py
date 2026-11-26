"""Calibrated reader for Fig. S2E (histogram of TMR-fascin domain sizes at 125 nM fascin + 400 nM alpha-actinin)."""

from __future__ import annotations

import numpy as np

import sys
from pathlib import Path

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, 2400), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 15, 10, 5, 0 counts
    yt_raw = C.find_ticks(dark, y_spine, (40, int(x_spine[1]) + 15), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 4, tol=8.0, min_step=600.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [15.0, 10.0, 5.0, 0.0])), max_residual_px=2.5)

    # Inward x ticks occur every 0.5 um (spacing ~67.5 px); left side is partly occluded by bars,
    # whereas y_spine is x=0 and ticks at 9.5, 10.0, 10.5, 11.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 16.0 are unobstructed!
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[0]) - 5, image.shape[1] - 20), "x", side=-1, band=(4, 20), frac=0.35, max_width=25)
    x0 = float(y_spine[0])
    x16 = max(xt_raw)
    step = (x16 - x0) / 32.0  # 32 half-um steps from 0 to 16 um
    x_marks = [(x0, 0.0), (x16, 16.0)]
    for k in range(1, 32):
        pred = x0 + k * step
        near = [t for t in xt_raw if abs(t - pred) <= 6.0]
        if near:
            x_marks.append((min(near, key=lambda t: abs(t - pred)), 0.5 * k))
    x_marks.sort(key=lambda t: t[1])
    x_axis = C.verify_axis("linear", x_marks, max_residual_px=4.5)

    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    cyan = (B > 170) & (G > 160) & (R < 140)
    r_base = int(x_spine[1]) - 8
    c0 = int(y_spine[2]) + 2
    cols_cyan = cyan[r_base - 90 : r_base, c0 : image.shape[1] - 20].mean(axis=0) > 0.5
    runs = [r for r in C.runs_1d(cols_cyan, gap=3) if (r[1] - r[0] + 1) >= 40]

    out = []
    for a, b in runs:
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        # Nothing cyan exists above any bar: take the topmost row with solid cyan fill across interior
        cyan_rows = np.where(cyan[60:r_base, ca + 10 : cb - 10].mean(axis=1) > 0.35)[0] + 60
        if len(cyan_rows) == 0:
            continue
        top_r = int(cyan_rows.min())
        blk_strip = dark[max(0, top_r - 20) : top_r + 2, ca + 10 : cb - 10].mean(axis=1) > 0.4
        blk_rows = np.where(blk_strip)[0] + max(0, top_r - 20)
        bar_top = float(blk_rows.mean()) + 0.5 if len(blk_rows) > 0 else float(top_r)
        xv = float(x_axis.to_value(c_mid))
        out.append(
            C.make_point_reading(
                xv,
                c_mid,
                bar_top,
                y_axis,
                err_px=3.0 + y_axis.residual_px,
                note=f"Fig. S2E histogram bin [{x_axis.to_value(ca):.2f}, {x_axis.to_value(cb):.2f}] um; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            )
        )
    return out
