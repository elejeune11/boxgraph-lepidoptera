"""Shared calibration and curve readers for Fig. 3A (p. 3) and its magnified Inset."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks, sample_curve_at


def calibrate_fig3a(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]
    r_in_bot = h[1][2]  # horizontal bottom border of the inset

    # Main axes: X linear (0.4, 0.8, 1.2, 1.6, 2.0), Y log10 (0.05..0.005)
    xt_bot = find_ticks(image, r_bot, int(c_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=25)
    if len(xt_bot) != 4:
        raise RuntimeError(f"Fig. 3A main expected 4 interior x-ticks (0.4..1.6), found {len(xt_bot)}")
    x_main = Axis("linear", [(xt_bot[0][0], 0.4), (xt_bot[1][0], 0.8), (xt_bot[2][0], 1.2), (xt_bot[3][0], 1.6), (c_right, 2.0)])

    yt_l = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=25)
    y_vals_main = [0.05, 0.04, 0.03, 0.02, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005]
    if len(yt_l) != len(y_vals_main):
        raise RuntimeError(f"Fig. 3A main expected {len(y_vals_main)} y-ticks, found {len(yt_l)}")
    y_main = Axis("log10", list(zip([p for p, _, _ in yt_l], y_vals_main)))

    # Inset vertical left border
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    v_in = np.nonzero(dark[int(r_top) : int(r_in_bot), :].mean(axis=0) > 0.85)[0]
    splits = np.where(np.diff(v_in) > 1)[0] + 1
    v_in_centers = [float(g.mean() + 0.5) for g in np.split(v_in, splits)]
    # Second vertical line is inset left border (first is main left border)
    c_in_left = [c for c in v_in_centers if c_left + 100 < c < c_right - 100][0]

    xt_in = find_ticks(image, r_in_bot, int(c_in_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=25)
    x_vals_in = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    if len(xt_in) != 6:
        raise RuntimeError(f"Fig. 3A inset expected 6 interior x-ticks (0.8..1.8), found {len(xt_in)}")
    x_inset = Axis("linear", list(zip([p for p, _, _ in xt_in] + [c_right], x_vals_in)))

    yt_in = find_ticks(image, c_right, int(r_top) + 8, int(r_in_bot) - 8, "y", -1, offset_lo=8, offset_hi=25)
    y_vals_in = [0.63, 0.62, 0.61, 0.60, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.50]
    if len(yt_in) != 13:
        raise RuntimeError(f"Fig. 3A inset expected 13 interior y-ticks (0.51..0.63), found {len(yt_in)}")
    y_inset = Axis("log10", list(zip([p for p, _, _ in yt_in] + [r_in_bot], y_vals_in)))

    boxes = {
        "main": (c_left, r_top, c_right, r_bot),
        "inset": (c_in_left, r_top, c_right, r_in_bot),
    }
    return x_main, y_main, x_inset, y_inset, boxes


def read_fig3a_series(image: np.ndarray, at: list[float], which: str):
    x_main, y_main, x_inset, y_inset, boxes = calibrate_fig3a(image)
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    c_l, r_t, c_r, r_b = boxes["main"]
    ci_l, ri_t, ci_r, ri_b = boxes["inset"]

    if which == "inset_R0p5":
        mask = np.zeros(image.shape[:2], bool)
        r0, r1 = int(ri_t) + 6, int(ri_b) - 6
        c0, c1 = int(ci_l) + 6, int(ci_r) - 6
        mask[r0:r1, c0:c1] = (R[r0:r1, c0:c1] > 160.0) & (G[r0:r1, c0:c1] < 80.0) & (B[r0:r1, c0:c1] < 80.0)
        return [
            sample_curve_at(mask, x_inset, y_inset, q, r0, r1, window=2, extra_px=1.0, note_prefix="inset R*=0.5; ")
            for q in at
        ]

    # Main plot region: exclude inset box (with small margin for inset tick labels below r_in_bot? No, inset labels are black!)
    # Only red curve appears inside inset; excluding inset rectangle prevents mixing inset red with main red.
    r0, r1 = int(r_t) + 6, int(r_b) - 6
    c0, c1 = int(c_l) + 6, int(c_r) - 6
    roi = np.zeros(image.shape[:2], bool)
    roi[r0:r1, c0:c1] = True
    roi[int(ri_t) : int(ri_b) + 12, int(ci_l) - 12 : int(ci_r) + 5] = False

    if which == "R0p5":
        mask = roi & (R > 160.0) & (G < 80.0) & (B < 80.0)
        prefix = "main R*=0.5 (red); "
    elif which == "R0p05":
        mask = roi & (B > 150.0) & (R < 80.0) & (G < 80.0)
        prefix = "main R*=0.05 (blue); "
    elif which == "R0p005":
        mask = roi & (G > 85.0) & (R < 60.0) & (B < 60.0)
        prefix = "main R*=0.005 (green); "
    else:
        raise ValueError(f"unknown Fig. 3A series {which!r}")

    return [sample_curve_at(mask, x_main, y_main, q, r0, r1, window=2, extra_px=1.0, note_prefix=prefix) for q in at]
