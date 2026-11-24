"""Shared calibrated reader for Figs. 2, 3, 4 (pages 6-7, 1500 dpi)."""

from __future__ import annotations

import math
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
RGB_GREEN = (44.0, 160.0, 44.0)
RGB_RED = (214.0, 39.0, 40.0)

PANEL_TICKS = {
    "2a": ([0.0, 2 * math.pi, 4 * math.pi], [0.2, 0.1, 0.0, -0.1, -0.2]),
    "2b": ([0.0, 2 * math.pi, 4 * math.pi], [1.0, 0.5, 0.0, -0.5, -1.0]),
    "3a": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [0.02, 0.01, 0.0, -0.01]),
    "3b": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2]),
    "4a": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [0.01, 0.0, -0.01]),
    "4b": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [0.4, 0.2, 0.0, -0.2, -0.4, -0.6, -0.8, -1.0]),
    "4c": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [0.6, 0.4, 0.2, 0.0, -0.2]),
    "4d": ([0.0, 2 * math.pi, 4 * math.pi, 6 * math.pi], [1.0, 0.5, 0.0, -0.5, -1.0, -1.5]),
}


def calibrate_panel(image: np.ndarray, panel: str):
    x_vals, y_vals = PANEL_TICKS[panel]
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    xt_m = match_ticks_regular(xt, len(x_vals), tol=3.5)
    yt_m = match_ticks_regular(yt, len(y_vals), tol=3.5)
    x_axis = verify_axis("linear", list(zip(xt_m, x_vals)), max_residual_px=1.0)
    y_axis = verify_axis("linear", list(zip(yt_m, y_vals)), max_residual_px=1.0)
    return ht, hb, vl, vr, x_axis, y_axis


def _mask_legend_and_boxes(mask: np.ndarray, panel: str, x_axis, y_axis):
    """Mask out legend entries inside the upper corners so they never interfere with curves."""
    if panel == "2b":
        # Legend in top-right (k > 8, a > 0.78)
        mask[: int(round(y_axis.to_pixel(0.78))), int(round(x_axis.to_pixel(8.0))) :] = False
    elif panel == "3b":
        # Legend in top-left (k < 10, a > 1.84)
        mask[: int(round(y_axis.to_pixel(1.84))), : int(round(x_axis.to_pixel(10.0)))] = False
    elif panel == "4a":
        # Legend in top-right above curve peak (y > 0.0080, k > 6.5)
        mask[: int(round(y_axis.to_pixel(0.0080))), int(round(x_axis.to_pixel(6.5))) :] = False
    elif panel == "4b":
        # Legend in top-left (k < 10, a > 0.15)
        mask[: int(round(y_axis.to_pixel(0.15))), : int(round(x_axis.to_pixel(10.0)))] = False
    elif panel == "4c":
        # Legend in top-right above y = 0.585 (peak is at y ~ 0.545, k ~ 12.0)
        mask[: int(round(y_axis.to_pixel(0.585))), int(round(x_axis.to_pixel(6.5))) :] = False
    elif panel == "4d":
        # Legend in top-right (k > 8.0, a > 0.45)
        mask[: int(round(y_axis.to_pixel(0.45))), int(round(x_axis.to_pixel(8.0))) :] = False


def read_panel_curve(
    image: np.ndarray,
    at: list[float],
    panel: str,
    color_name: str,
    style: str = "solid",
    branch: str = "top",
    dashed_x_range: tuple[float, float] | None = None,
) -> list[dict]:
    ht, hb, vl, vr, x_axis, y_axis = calibrate_panel(image, panel)
    r_min, r_max = ht[2] + 5, hb[1] - 5

    rgb_map = {"blue": RGB_BLUE, "green": RGB_GREEN, "red": RGB_RED}
    mask = mask_interior(rgb_mask(image, rgb_map[color_name], tol=60.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
    _mask_legend_and_boxes(mask, panel, x_axis, y_axis)

    pick_fn = None
    if branch == "top":
        pick_fn = lambda rs, col: rs[0]
    elif branch == "bottom":
        pick_fn = lambda rs, col: rs[-1]

    out = []
    for xv in at:
        xf = float(xv)
        use_dashed = (style == "dashed") or (
            dashed_x_range is not None and dashed_x_range[0] <= xf <= dashed_x_range[1]
        )
        if use_dashed:
            r_lo, r_hi = r_min, r_max
            if panel == "4c" and branch == "top":
                r_hi = min(r_max, int(round(y_axis.to_pixel(-0.15))))
            out.append(
                sample_dashed_curve_y(
                    mask,
                    x_axis,
                    y_axis,
                    xf,
                    r_lo,
                    r_hi,
                    max_search_dx=60,
                    note_prefix=f"Fig. {panel} {color_name}, ",
                )
            )
        else:
            r_lo, r_hi = r_min, r_max
            if panel == "2a" and branch == "bottom":
                # Lower branch lam_-(k) is strictly negative for k >= 4.5
                r_lo = max(r_min, int(round(y_axis.to_pixel(-0.005))))
            out.append(
                sample_solid_curve_y(
                    mask,
                    x_axis,
                    y_axis,
                    xf,
                    r_lo,
                    r_hi,
                    pick_fn=pick_fn,
                    note_prefix=f"Fig. {panel} {color_name}, ",
                )
            )
    return out
