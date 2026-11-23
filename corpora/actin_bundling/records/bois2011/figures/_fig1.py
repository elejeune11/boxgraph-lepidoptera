"""Shared calibration and curve isolation for Fig. 1(a) and Fig. 1(b) (two stacked subpanels)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _calib import (
    find_frame,
    find_inward_ticks,
    match_ticks_regular,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def calibrate_fig1(image: np.ndarray, subpanel: str):
    """Calibrate top ('top') or bottom ('bottom') subpanel of Fig. 1(a) or 1(b)."""
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=110.0)
    if len(h_lines) < 3:
        raise RuntimeError(f"Expected 3 horizontal frame lines in Fig. 1, found {len(h_lines)}: {h_lines}")
    h_top, h_mid, h_bot = h_lines[0], h_lines[1], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    if subpanel == "top":
        ht, hb = h_top, h_mid
    elif subpanel == "bottom":
        ht, hb = h_mid, h_bot
    else:
        raise ValueError(f"Unknown subpanel {subpanel!r}")

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, ht, hb, v_left, v_right)
    # x ticks: left border (0), 6 interior ticks (1, 2, 3, 4, 5, 6)
    x_inner = match_ticks_regular(xt_bot, 6, tol=2.0)
    x_marks = [(v_left[0], 0.0)] + list(zip(x_inner, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
    x_axis = verify_axis("linear", x_marks, max_residual_px=1.0)

    if subpanel == "top":
        # Left y-axis (c): 4 interior ticks at 8, 6, 4, 2 (tick at 0 is slightly above bottom border and hit by c)
        y_left_inner = match_ticks_regular(yt_left, 4, tol=2.0)
        y_left = verify_axis("linear", list(zip(y_left_inner, [8.0, 6.0, 4.0, 2.0])), max_residual_px=1.0)
        # Right y-axis (v/U): 3 regular interior ticks at +0.2, 0.0, -0.2
        y_right_inner = match_ticks_regular(yt_right, 3, tol=2.0)
        y_right = verify_axis("linear", list(zip(y_right_inner, [0.2, 0.0, -0.2])), max_residual_px=1.0)
    else:
        # Bottom subpanel (stress): 5 interior ticks on right border at 0.8, 0.6, 0.4, 0.2, 0.0
        y_stress_inner = match_ticks_regular(yt_right, 5, tol=2.0)
        y_left = verify_axis("linear", list(zip(y_stress_inner, [0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=1.0)
        y_right = y_left

    # Note: keep r0 close to top border (ht[2]+4) because peak of c in Fig. 1(a) approaches the top border
    box = (v_left[2] + 15, ht[2] + 4, v_right[1] - 15, hb[1] - 12)
    return dark, x_axis, y_left, y_right, box


def clear_legend(mask: np.ndarray, x_axis, y_axis, box, subpanel: str) -> np.ndarray:
    """Detect and clear the rectangular legend box in the upper region of the subpanel."""
    c0, r0, c1, r1 = box
    out = mask.copy()
    # Legend box sits in upper 42% of the subpanel, away from left/right borders
    search_r0 = r0 + 10
    search_r1 = r0 + int(0.42 * (r1 - r0))
    search_c0 = c0 + int(0.15 * (c1 - c0))
    search_c1 = c1 - 15

    # Find horizontal border lines of the legend box (long thin horizontal segments > 220 px)
    sub = out[search_r0:search_r1, search_c0:search_c1]
    h_rows = []
    for r in range(sub.shape[0]):
        row = sub[r]
        # longest dark run in this row
        idx = np.nonzero(row)[0]
        if len(idx) >= 220:
            splits = np.where(np.diff(idx) > 2)[0] + 1
            for g in np.split(idx, splits):
                if len(g) >= 220:
                    h_rows.append((r + search_r0, int(g[0]) + search_c0, int(g[-1]) + search_c0))
    if len(h_rows) >= 2:
        r_min = min(t[0] for t in h_rows)
        r_max = max(t[0] for t in h_rows)
        c_min = min(t[1] for t in h_rows)
        c_max = max(t[2] for t in h_rows)
        if r_max - r_min >= 100 and c_max - c_min >= 220:
            out[max(r0, r_min - 12) : min(r1 + 1, r_max + 12), max(c0, c_min - 12) : min(c1 + 1, c_max + 12)] = False
    return out


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read 'c' (top solid), 'v' (top dashed), 'sigma' (bottom solid), or 'zdm' (bottom dashed)."""
    subpanel = "top" if curve in ("c", "v") else "bottom"
    dark, x_axis, y_left, y_right, box = calibrate_fig1(image, subpanel)
    c0, r0, c1, r1 = box

    interior = np.zeros_like(dark)
    interior[r0 : r1 + 1, c0 : c1 + 1] = dark[r0 : r1 + 1, c0 : c1 + 1]
    interior = clear_legend(interior, x_axis, y_left, box, subpanel)

    # Separate solid curve (single huge connected component spanning width, area > 10000 px)
    # from dashed curve (discrete dashes of area ~150..2000 px)
    lab, n = ndimage.label(interior)
    counts = np.bincount(lab.ravel())
    solid_mask = np.zeros_like(interior)
    dash_mask = np.zeros_like(interior)
    for idx in range(1, n + 1):
        area = int(counts[idx])
        if area >= 6000:
            solid_mask |= lab == idx
        elif area >= 80:
            dash_mask |= lab == idx

    is_solid = curve in ("c", "sigma")
    y_axis = y_right if curve == "v" else y_left

    out = []
    for x in at:
        if is_solid:
            res = sample_solid_curve_y(
                solid_mask,
                x_axis,
                y_axis,
                float(x),
                r0,
                r1,
                note_prefix=f"Fig. 1 {curve} solid, ",
            )
        else:
            res = sample_dashed_curve_y(
                dash_mask,
                x_axis,
                y_axis,
                float(x),
                r0,
                r1,
                max_search_dx=60,
                note_prefix=f"Fig. 1 {curve} dashed, ",
            )
        out.append(res)
    return out
