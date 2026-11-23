"""Calibration and curve isolation for Fig. 3(a): single-species dispersion relation for Pe = 25, 10, 7, 3."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _calib import (
    find_frame,
    find_inward_ticks,
    match_ticks_regular,
    runs_1d,
    sample_solid_curve_y,
    verify_axis,
)


def calibrate_fig3a(image: np.ndarray):
    # Threshold at 95 excludes the vertical gray lines (gray ~ 149) while keeping frame + black curves
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=95.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    # x ticks: left border (0), 2 interior ticks (1, 2), right border (3)
    x_inner = match_ticks_regular(xt_bot, 2, tol=2.0)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], 0.0), (x_inner[0], 1.0), (x_inner[1], 2.0), (v_right[0], 3.0)],
        max_residual_px=1.0,
    )
    # y ticks: 5 interior ticks on right border (2, 1, 0, -1, -2)
    y_inner = match_ticks_regular(yt_right, 5, tol=2.0)
    y_axis = verify_axis(
        "linear",
        list(zip(y_inner, [2.0, 1.0, 0.0, -1.0, -2.0])),
        max_residual_px=1.0,
    )
    # Horizontal zero axis line is the middle h_line (~row 751.5)
    h_zero = h_lines[1] if len(h_lines) >= 3 else None
    box = (v_left[2] + 15, h_top[2] + 15, v_right[1] - 15, h_bot[1] - 15)
    return dark, x_axis, y_axis, box, h_zero


def isolate_thick_curves(dark: np.ndarray, box, h_zero) -> np.ndarray:
    """Keep only thick dispersion curves (~18-30 px) and remove thin zero-line, arrow, and text (~6-8 px)."""
    c0, r0, c1, r1 = box
    interior = np.zeros_like(dark)
    interior[r0 : r1 + 1, c0 : c1 + 1] = dark[r0 : r1 + 1, c0 : c1 + 1]

    # Morphological opening with a 11x11 disk removes thin 6-8 px strokes (zero line, arrow, 'increasing Pe' text)
    # while leaving the 18-22 px dispersion curves completely intact, then restore original stroke edges inside dilation.
    yy, xx = np.ogrid[-5:6, -5:6]
    disk = (xx * xx + yy * yy) <= 25
    core = ndimage.binary_opening(interior, structure=disk)
    # Keep only large curve components (area > 800 px) to drop any thick serif letter blobs
    lab, n = ndimage.label(core)
    counts = np.bincount(lab.ravel())
    clean_core = np.zeros_like(core)
    for idx in range(1, n + 1):
        if counts[idx] >= 800:
            clean_core |= lab == idx
    # Dilate clean_core slightly and intersect with interior (excluding pure thin horizontal zero line)
    # to recover exact un-eroded curve stroke edges without re-attaching arrow/text.
    # First remove horizontal zero line from interior where it's not near clean_core
    dil = ndimage.binary_dilation(clean_core, structure=np.ones((15, 15), bool))
    mask = interior & dil
    if h_zero is not None:
        _, hz0, hz1 = h_zero
        # Remove horizontal zero line pixels that are only thin (no core within 12 px vertically)
        for col in range(c0, c1 + 1):
            if not np.any(clean_core[max(r0, hz0 - 6) : min(r1, hz1 + 6), col]):
                mask[hz0 - 1 : hz1 + 2, col] = False
    return mask


def read_curve(image: np.ndarray, at: list[float], rank_from_top: int, label: str) -> list[dict]:
    """Read one of the 4 strictly ordered dispersion curves in Fig. 3(a):
    rank_from_top = 0 (Pe=25), 1 (Pe=10), 2 (Pe=7), 3 (Pe=3).
    """
    dark, x_axis, y_axis, box, h_zero = calibrate_fig3a(image)
    mask = isolate_thick_curves(dark, box, h_zero)
    c0, r0, c1, r1 = box

    # Trace each curve continuously across columns starting from k*l = 0.8 (where all 4 are well-separated)
    c_seed = int(round(x_axis.to_pixel(0.8)))
    seed_runs = [(a + r0, b + r0) for a, b in runs_1d(mask[r0 : r1 + 1, c_seed])]
    if len(seed_runs) != 4:
        raise RuntimeError(f"Expected 4 dispersion curves at k*l=0.8 (col {c_seed}), found {len(seed_runs)}: {seed_runs}")

    # Trace target curve left and right from c_seed to record its center row at every column
    H, W = mask.shape
    trace = np.full(W, np.nan)
    r_cur = 0.5 * (seed_runs[rank_from_top][0] + seed_runs[rank_from_top][1]) + 0.5
    trace[c_seed] = r_cur

    for direction in (1, -1):
        r_prev = r_cur
        col = c_seed + direction
        while c0 <= col <= c1:
            rs = [(a + r0, b + r0) for a, b in runs_1d(mask[r0 : r1 + 1, col])]
            if not rs:
                break
            mids = [0.5 * (a + b) + 0.5 for a, b in rs]
            best_i = int(np.argmin([abs(m - r_prev) for m in mids]))
            if abs(mids[best_i] - r_prev) > 35.0:
                break
            r_prev = mids[best_i]
            trace[col] = r_prev
            col += direction

    out = []
    for x in at:
        c_f = float(x_axis.to_pixel(float(x)))
        c_int = int(round(c_f))
        r_tr = trace[c_int]
        if np.isnan(r_tr):
            out.append({"x": float(x), "y": None, "found": False, "note": f"{label} trace lost at k*l={x}"})
            continue
        res = sample_solid_curve_y(
            mask,
            x_axis,
            y_axis,
            float(x),
            r0,
            r1,
            pick_fn=lambda rs, col, target=r_tr: min(rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - target)),
            note_prefix=f"{label}, ",
        )
        out.append(res)
    return out
