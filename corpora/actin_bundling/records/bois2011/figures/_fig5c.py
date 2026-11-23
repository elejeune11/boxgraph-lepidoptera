"""Shared calibration and curve isolation for Fig. 5(c): ASDM dispersion relations."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _calib import (
    find_frame,
    find_inward_ticks,
    match_ticks_regular,
    runs_1d,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def calibrate_fig5c(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=110.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    # x ticks: left border (0), 2 interior ticks (2, 4), right border (6)
    x_inner = match_ticks_regular(xt_bot, 2, tol=2.0)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], 0.0), (x_inner[0], 2.0), (x_inner[1], 4.0), (v_right[0], 6.0)],
        max_residual_px=1.0,
    )
    # y ticks: top border (10), 2 interior ticks (5, 0), bottom border (-5)
    y_inner = match_ticks_regular(yt_right, 2, tol=2.0)
    y_axis = verify_axis(
        "linear",
        [(h_top[0], 10.0), (y_inner[0], 5.0), (y_inner[1], 0.0), (h_bot[0], -5.0)],
        max_residual_px=1.0,
    )
    h_zero = h_lines[1] if len(h_lines) >= 3 else None
    box = (v_left[2] + 15, h_top[2] + 15, v_right[1] - 15, h_bot[1] - 15)
    return dark, x_axis, y_axis, box, h_zero


def isolate_curves(dark: np.ndarray, box, h_zero):
    """Separate thick solid curves and thick dashed curve from thin zero line and text labels."""
    c0, r0, c1, r1 = box
    interior = np.zeros_like(dark)
    interior[r0 : r1 + 1, c0 : c1 + 1] = dark[r0 : r1 + 1, c0 : c1 + 1]

    # Erase thin horizontal zero line (~6 px) except where thick curves cross it
    # Morphological opening with 11x11 disk preserves ~16-20 px curve strokes (both solid & dashes)
    # and removes 6-8 px zero line and thin serif label strokes ("Pe=50", "Pe=25", "Pe=0").
    yy, xx = np.ogrid[-5:6, -5:6]
    disk = (xx * xx + yy * yy) <= 25
    core = ndimage.binary_opening(interior, structure=disk)

    lab, n = ndimage.label(core)
    counts = np.bincount(lab.ravel())

    # Solid curves have core area > 3500 px; dashes of the Pe=0 (+) curve have core area ~80..1200 px
    # (thin serif letters have area 0 after opening with 11x11 disk, or tiny blobs < 70 px)
    solid_core = np.zeros_like(core)
    dash_core = np.zeros_like(core)
    for idx in range(1, n + 1):
        area = int(counts[idx])
        if area >= 3500:
            solid_core |= lab == idx
        elif area >= 80:
            dash_core |= lab == idx

    # Reconstruct full un-eroded strokes by intersecting interior with dilated core,
    # while masking thin zero-line pixels away from crossings.
    dil_s = ndimage.binary_dilation(solid_core, structure=np.ones((15, 15), bool))
    solid_mask = interior & dil_s
    dil_d = ndimage.binary_dilation(dash_core, structure=np.ones((13, 13), bool))
    dash_mask = interior & dil_d & (~solid_mask)

    if h_zero is not None:
        _, hz0, hz1 = h_zero
        for col in range(c0, c1 + 1):
            if not np.any(solid_core[max(r0, hz0 - 5) : min(r1, hz1 + 5), col]):
                solid_mask[hz0 - 1 : hz1 + 2, col] = False
            if not np.any(dash_core[max(r0, hz0 - 5) : min(r1, hz1 + 5), col]):
                dash_mask[hz0 - 1 : hz1 + 2, col] = False

    return solid_mask, dash_mask


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read one of 'Pe50', 'Pe25', 'Pe0_cross' (solid), or 'Pe0_plus' (dashed) from Fig. 5(c)."""
    dark, x_axis, y_axis, box, h_zero = calibrate_fig5c(image)
    solid_mask, dash_mask = isolate_curves(dark, box, h_zero)
    c0, r0, c1, r1 = box

    if curve == "Pe0_plus":
        # Exclude any stray text-label blobs (from "Pe=25" above y=1.5 and "Pe=0" below y=-2.5)
        # by keeping only dash components whose centroid lies within the vertical span of the dashed curve
        # (starting near y=0 at k*l=2 and declining monotonically to y~-1.8 at k*l=6).
        r_hi_band = int(round(y_axis.to_pixel(1.0)))
        r_lo_band = int(round(y_axis.to_pixel(-2.3)))
        clean_dash = np.zeros_like(dash_mask)
        clean_dash[r_hi_band : r_lo_band + 1, c0 : c1 + 1] = dash_mask[r_hi_band : r_lo_band + 1, c0 : c1 + 1]
        # Also remove any stray text blob of "Pe=0" below the dashed curve in k*l in [2.2, 3.6]
        r_cut = int(round(y_axis.to_pixel(-0.5)))
        c_p0_a = int(round(x_axis.to_pixel(2.2)))
        c_p0_b = int(round(x_axis.to_pixel(3.6)))
        clean_dash[r_cut : r_lo_band + 1, c_p0_a : c_p0_b + 1] = False
        return [
            sample_dashed_curve_y(
                clean_dash, x_axis, y_axis, float(x), r0, r1, max_search_dx=65, note_prefix="Pe=0 (+) dashed, "
            )
            for x in at
        ]

    # Trace the 3 solid curves starting from k*l = 2.0 (where Pe=50, Pe=25, Pe=0 are cleanly ordered top-to-bottom)
    rank_map = {"Pe50": 0, "Pe25": 1, "Pe0_cross": 2}
    if curve not in rank_map:
        raise ValueError(f"Unknown Fig. 5(c) curve {curve!r}")
    rank = rank_map[curve]

    c_seed = int(round(x_axis.to_pixel(2.0)))
    seed_runs = [(a + r0, b + r0) for a, b in runs_1d(solid_mask[r0 : r1 + 1, c_seed])]
    if len(seed_runs) != 3:
        raise RuntimeError(f"Expected 3 solid curves at k*l=2.0 (col {c_seed}), found {len(seed_runs)}: {seed_runs}")

    H, W = solid_mask.shape
    trace = np.full(W, np.nan)
    r_cur = 0.5 * (seed_runs[rank][0] + seed_runs[rank][1]) + 0.5
    trace[c_seed] = r_cur
    for direction in (1, -1):
        r_prev = r_cur
        col = c_seed + direction
        while c0 <= col <= c1:
            rs = [(a + r0, b + r0) for a, b in runs_1d(solid_mask[r0 : r1 + 1, col])]
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
        c_int = int(round(x_axis.to_pixel(float(x))))
        r_tr = trace[c_int]
        if np.isnan(r_tr):
            out.append({"x": float(x), "y": None, "found": False, "note": f"{curve} trace lost at k*l={x}"})
            continue
        res = sample_solid_curve_y(
            solid_mask,
            x_axis,
            y_axis,
            float(x),
            r0,
            r1,
            pick_fn=lambda rs, col, target=r_tr: min(rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - target)),
            note_prefix=f"{curve} solid, ",
        )
        out.append(res)
    return out
