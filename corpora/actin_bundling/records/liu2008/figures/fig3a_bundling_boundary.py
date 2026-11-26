"""Calibrated reader for Fig. 3a (p. 3 / printed p. 791): membrane-induced bundling phase boundary."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.6, min_v_frac=0.45)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Inward x ticks: 40, 60, 80, 100, 120, 140, 160, 180 nm
    xt_all = C.scan_ticks(dark, "x", bot[0] - 26, bot[0] - 6, left[1] + 5, right[0] - 5, min_frac=0.55, max_len=20)
    yt_all = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 6, left[1] + 26, min_frac=0.6, max_len=20)

    xt_reg = C.select_regular_ticks(xt_all, 8, tol=2.5)
    xa = C.make_axis("linear", list(zip(xt_reg, [40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0])), max_res_px=1.5)

    yt_reg = C.select_regular_ticks(yt_all, 5, tol=2.5)
    ya = C.make_axis("linear", list(zip(yt_reg, [100.0, 80.0, 60.0, 40.0, 20.0])), max_res_px=1.5)

    # Dark-blue phase boundary stroke has crisp RGB ~ (29, 92, 154)
    rgb = image.astype(np.float32)
    dist_blue = np.abs(rgb[:, :, 0] - 29.0) + np.abs(rgb[:, :, 1] - 92.0) + np.abs(rgb[:, :, 2] - 154.0)
    r0, r1 = top[1] + 6, bot[0] - 2
    mask = np.zeros(image.shape[:2], dtype=bool)
    mask[r0:r1, left[1] + 6 : right[0] - 6] = dist_blue[r0:r1, left[1] + 6 : right[0] - 6] < 75.0
    # Exclude the two inset schematic boxes
    mask[: int(round(ya.to_pixel(45.0))), : int(round(xa.to_pixel(100.0)))] = False
    mask[:, int(round(xa.to_pixel(138.0))) :] = False

    # Build unobstructed column profile: skip columns clipped by inward black ticks (at multiples of 20 nm)
    # or the vertical dashed black line (at L_prot^B = 100 nm and L_prot^F = 130.4 nm), then interpolate cleanly.
    c_cols = []
    c_rows = []
    c_half = []
    c_end = min(right[0] - 6, int(round(xa.to_pixel(136.0))))
    for c in range(left[1] + 6, c_end):
        runs = C.find_runs(mask[r0:r1, c], min_len=36)  # unobstructed blue stroke is >= 45 px thick vertically
        if runs:
            best = max(runs, key=lambda t: t[1] - t[0])
            c_cols.append(float(c) + 0.5)
            c_rows.append(r0 + best[2])
            c_half.append(0.5 * (best[1] - best[0]))
    cols_arr = np.array(c_cols)
    rows_arr = np.array(c_rows)
    half_arr = np.array(c_half)

    out = []
    for xv in at:
        cf = xa.to_pixel(xv)
        idx = int(np.searchsorted(cols_arr, cf))
        if idx == 0 or idx >= len(cols_arr):
            k = 0 if idx == 0 else len(cols_arr) - 1
            if abs(cols_arr[k] - cf) > 25.0:
                out.append({"x": float(xv), "y": None, "note": f"Fig. 3a: x={xv:g} outside curve"})
                continue
            row_ctr = float(rows_arr[k])
            half_thick = float(half_arr[k])
            gap = abs(cols_arr[k] - cf)
            slope_px = 0.0
        else:
            cL, cR = cols_arr[idx - 1], cols_arr[idx]
            rL, rR = rows_arr[idx - 1], rows_arr[idx]
            frac = (cf - cL) / max(cR - cL, 1e-9)
            row_ctr = float(rL + frac * (rR - rL))
            half_thick = float(max(half_arr[idx - 1], half_arr[idx]))
            gap = min(cf - cL, cR - cf)
            slope_px = abs(rR - rL) / max(cR - cL, 1.0)
        col_err = max(xa.residual_px, 0.5) + 0.25 * gap
        half_px = math.hypot(half_thick, slope_px * col_err)
        note = (
            f"Fig. 3a bundling boundary (sigma=0.005): width={2*half_thick:.1f}px, "
            f"cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
            + (f", interp_gap={gap:.1f}px" if gap > 2.0 else "")
        )
        out.append(C.reading_point(x_val=xv, col=cf, row=row_ctr, y_axis=ya, half_px=half_px, note=note, clamp_min=0.0))
    return out
