"""Shared calibrated reader for Figs. 5 and 6 (pages 9-10, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_main_box,
    find_outward_ticks,
    mask_interior,
    match_ticks_regular,
    rgb_mask,
    runs_1d,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)

RGB_BLUE = (31.0, 119.0, 180.0)
RGB_SKY = (135.0, 206.0, 250.0)
RGB_GREEN = (44.0, 160.0, 44.0)
RGB_ORANGE = (255.0, 127.0, 14.0)
RGB_RED = (214.0, 39.0, 40.0)
RGB_PURPLE = (148.0, 103.0, 189.0)


def calibrate_fig56(image: np.ndarray, panel: str):
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    if panel == "5a":
        x_vals = [-0.5, 0.0, 0.5, 1.0]
        y_vals = [1.50, 1.25, 1.00, 0.75, 0.50, 0.25, 0.00]
    elif panel == "5b":
        x_vals = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5]
        y_vals = [2.0, 1.5, 1.0, 0.5, 0.0]
    elif panel == "5c":
        x_vals = [-0.4, -0.2, 0.0, 0.2, 0.4]
        y_vals = [1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0]
    elif panel == "6a":
        x_vals = [-0.5, 0.0, 0.5, 1.0, 1.5]
        y_vals = [2.0, 1.5, 1.0, 0.5, 0.0]
    elif panel == "6b":
        x_vals = [1.30, 1.35, 1.40, 1.45, 1.50]
        y_vals = [0.6, 0.4, 0.2, 0.0]
    elif panel == "6c":
        x_vals = [-0.4, -0.2, 0.0, 0.2, 0.4]
        y_vals = [0.75, 0.50, 0.25, 0.00, -0.25, -0.50, -0.75]
    else:
        raise ValueError(f"Unknown panel {panel!r}")

    xt_m = match_ticks_regular(xt, len(x_vals), tol=3.5)
    yt_m = match_ticks_regular(yt, len(y_vals), tol=3.5)
    x_axis = verify_axis("linear", list(zip(xt_m, x_vals)), max_residual_px=1.0)
    y_axis = verify_axis("linear", list(zip(yt_m, y_vals)), max_residual_px=1.0)
    return ht, hb, vl, vr, x_axis, y_axis


def _sample_ith_dashed_branch(
    mask: np.ndarray,
    x_axis,
    y_axis,
    x_val: float,
    row_min: int,
    row_max: int,
    branch_idx: int,
    half_win: int = 35,
    note_prefix: str = "",
) -> dict:
    """For multiple dashed curves of the same color (Fig. 5a light-blue n=2,3,...),
    locate distinct curve bands in a vertical window around col(x_val) and sample the
    branch_idx-th curve from the top."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    # Gather row projections in a small horizontal band to cluster separate branches
    c0 = max(0, c - half_win)
    c1 = min(W - 1, c + half_win)
    # Use a narrow strip +/- 18 px to separate nearby curves vertically before slope tilts them
    cs0, cs1 = max(0, c - 18), min(W - 1, c + 18)
    row_proj = mask[row_min : row_max + 1, cs0 : cs1 + 1].any(axis=1)
    bands = [(a + row_min, b + row_min) for a, b in runs_1d(row_proj, gap=14)]
    if len(bands) <= branch_idx:
        return {"x": float(x_val), "y": None, "found": False, "note": f"only {len(bands)} bands near col {c}"}
    b_lo, b_hi = bands[branch_idx]
    # Expand window vertically by 25 px (clamped halfway to neighboring bands)
    r_lo = max(row_min, b_lo - 25)
    r_hi = min(row_max, b_hi + 25)
    if branch_idx > 0:
        prev_mid = 0.5 * (bands[branch_idx - 1][1] + b_lo)
        r_lo = max(r_lo, int(round(prev_mid)) + 2)
    if branch_idx + 1 < len(bands):
        next_mid = 0.5 * (b_hi + bands[branch_idx + 1][0])
        r_hi = min(r_hi, int(round(next_mid)) - 2)
    return sample_dashed_curve_y(
        mask, x_axis, y_axis, x_val, r_lo, r_hi, max_search_dx=half_win, note_prefix=note_prefix
    )


def read_fig56_curve(
    image: np.ndarray,
    at: list[float],
    panel: str,
    color_name: str,
    style: str = "solid",
    branch_idx: int = 0,
    dashed_above_x: float | None = None,
) -> list[dict]:
    ht, hb, vl, vr, x_axis, y_axis = calibrate_fig56(image, panel)
    r_min, r_max = ht[2] + 5, hb[1] - 5

    rgb_map = {
        "blue": RGB_BLUE,
        "sky": RGB_SKY,
        "green": RGB_GREEN,
        "orange": RGB_ORANGE,
        "red": RGB_RED,
        "purple": RGB_PURPLE,
    }
    from scipy import ndimage
    tol = 32.0 if color_name == "sky" else (50.0 if color_name == "blue" else 60.0)
    mask = mask_interior(rgb_mask(image, rgb_map[color_name], tol=tol, min_chroma=35.0), ht, hb, vl, vr, inset=4)
    if color_name == "sky":
        dark_blue = mask_interior(rgb_mask(image, RGB_BLUE, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
        mask &= ~ndimage.binary_dilation(dark_blue, iterations=7)

    # Mask legend in Fig. 5(c) (bottom-left: x < -0.25, y < -1.0)
    if panel == "5c":
        c_leg = int(round(x_axis.to_pixel(-0.25)))
        r_leg = int(round(y_axis.to_pixel(-1.0)))
        mask[r_leg:, :c_leg] = False

    # Avoid y = 0 horizontal axis clutter in bifurcation diagrams
    if panel in ("5a", "5b", "6a", "6b"):
        r_max = min(r_max, int(round(y_axis.to_pixel(0.03))))

    out = []
    for xv in at:
        xf = float(xv)
        if panel == "5a" and color_name == "sky":
            out.append(
                _sample_ith_dashed_branch(
                    mask, x_axis, y_axis, xf, r_min, r_max, branch_idx,
                    half_win=38, note_prefix=f"Fig. {panel} n={branch_idx+2}, "
                )
            )
            continue

        use_dashed = (style == "dashed") or (dashed_above_x is not None and xf >= dashed_above_x)
        if use_dashed:
            out.append(
                sample_dashed_curve_y(
                    mask, x_axis, y_axis, xf, r_min, r_max,
                    max_search_dx=40, note_prefix=f"Fig. {panel} {color_name}, "
                )
            )
        else:
            out.append(
                sample_solid_curve_y(
                    mask, x_axis, y_axis, xf, r_min, r_max,
                    pick_fn=lambda rs, col: rs[0],
                    note_prefix=f"Fig. {panel} {color_name}, "
                )
            )
    return out
