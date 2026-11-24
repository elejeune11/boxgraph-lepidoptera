"""Shared calibrated reader for Figs. 14(a)-(d), 15(a)-(b), and 16(a) (pages 23, 24, 26)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_main_box,
    find_outward_ticks,
    mask_interior,
    match_ticks_regular,
    rgb_mask,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)

RGB_BLUE = (31.0, 119.0, 180.0)
RGB_ORANGE = (255.0, 127.0, 14.0)
RGB_GREEN = (44.0, 160.0, 44.0)


def read_fig14(image: np.ndarray, at: list[float], panel: str) -> list[dict]:
    """Read heavy solid I-IV coexistence binodal (14a, 14c) or phase boundary (14b, 14d)."""
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt_maj, yt_maj = find_outward_ticks(dark, ht, hb, vl, vr, band=(10, 15), max_width=14)
    if panel in ("14a", "14c"):
        x_vals = [-2.0, -1.0, 0.0, 1.0, 2.0]
        y_vals = [2.0, 1.0, 0.0, -1.0, -2.0]
    elif panel == "14b":
        x_vals = [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
        y_vals = [4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0]
    elif panel == "14d":
        x_vals = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
        y_vals = [3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0]
    else:
        raise ValueError(f"Unknown panel {panel}")

    xt_m = match_ticks_regular(xt_maj, len(x_vals), tol=4.0)
    yt_m = match_ticks_regular(yt_maj, len(y_vals), tol=4.0)
    x_axis = verify_axis("linear", list(zip(xt_m, x_vals)), max_residual_px=1.5)
    y_axis = verify_axis("linear", list(zip(yt_m, y_vals)), max_residual_px=1.5)

    if panel in ("14a", "14b"):
        # Black heavy solid curve
        rgb = image.astype(np.float32)
        chroma = rgb.max(axis=2) - rgb.min(axis=2)
        mask = mask_interior((rgb.max(axis=2) < 90.0) & (chroma < 25.0), ht, hb, vl, vr, inset=6)
    else:
        # Blue heavy solid curve in 14c and 14d
        mask = mask_interior(rgb_mask(image, RGB_BLUE, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=6)

    if panel == "14a":
        r_min = int(round(y_axis.to_pixel(1.60)))
        r_max = int(round(y_axis.to_pixel(1.20)))
        pick = lambda rs, col: rs[0]
    elif panel == "14c":
        r_min = int(round(y_axis.to_pixel(1.70)))
        r_max = int(round(y_axis.to_pixel(1.38)))
        pick = lambda rs, col: rs[0]
    elif panel == "14b":
        r_min = int(round(y_axis.to_pixel(-1.20)))
        r_max = int(round(y_axis.to_pixel(-2.15)))
        pick = lambda rs, col: max(rs, key=lambda r: r[1] - r[0])
    else:  # 14d
        r_min = int(round(y_axis.to_pixel(-1.30)))
        r_max = int(round(y_axis.to_pixel(-2.05)))
        pick = lambda rs, col: max(rs, key=lambda r: r[1] - r[0])

    out = []
    for xv in at:
        out.append(
            sample_solid_curve_y(
                mask,
                x_axis,
                y_axis,
                float(xv),
                r_min,
                r_max,
                pick_fn=pick,
                note_prefix=f"Fig. {panel} I-IV, ",
            )
        )
    return out


def read_fig15a(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    """Read black solid bifurcation curve ||delta phi||(phi_bar_2) in Fig. 15(a)."""
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    xt_m = match_ticks_regular(xt, 5, tol=4.0)
    yt_m = match_ticks_regular(yt, 5, tol=4.0)
    x_axis = verify_axis("linear", list(zip(xt_m, [-2.0, -1.5, -1.0, -0.5, 0.0])), max_residual_px=1.5)
    y_axis = verify_axis("linear", list(zip(yt_m, [2.0, 1.8, 1.6, 1.4, 1.2])), max_residual_px=1.5)

    rgb = image.astype(np.float32)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    # Keep solid black line only (gray '+' markers have rgb ~ 128)
    mask = mask_interior((rgb.max(axis=2) < 90.0) & (chroma < 25.0), ht, hb, vl, vr, inset=6)
    # Mask panel label '(a)' in upper left
    mask[ht[2] : int(y_axis.to_pixel(1.90)), vl[2] : int(x_axis.to_pixel(-1.60))] = False

    if branch == "left":
        r_min = int(round(y_axis.to_pixel(1.40)))
        r_max = hb[1] - 6
        pick = lambda rs, col: rs[-1]
    else:
        r_min = ht[2] + 6
        r_max = int(round(y_axis.to_pixel(1.40)))
        pick = lambda rs, col: rs[0]

    out = []
    for xv in at:
        out.append(
            sample_solid_curve_y(
                mask,
                x_axis,
                y_axis,
                float(xv),
                r_min,
                r_max,
                pick_fn=pick,
                note_prefix=f"Fig. 15a {branch}, ",
            )
        )
    return out


def read_fig15b_profile(image: np.ndarray, at: list[float], field: str) -> list[dict]:
    """Read phi_1(x) (blue dashed) or phi_2(x) (orange solid) in Fig. 15(b) at phi_bar_2 = 0."""
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    xt_m = match_ticks_regular(xt, 3, tol=4.0)
    yt_m = match_ticks_regular(yt, 3, tol=4.0)
    x_axis = verify_axis("linear", list(zip(xt_m, [-0.5, 0.0, 0.5])), max_residual_px=1.5)
    y_axis = verify_axis("linear", list(zip(yt_m, [1.0, 0.0, -1.0])), max_residual_px=1.5)

    r_min, r_max = ht[2] + 5, hb[1] - 5
    if field == "phi1":
        mask = mask_interior(rgb_mask(image, RGB_BLUE, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
        out = [
            sample_dashed_curve_y(
                mask, x_axis, y_axis, float(xv), r_min, r_max,
                max_search_dx=35, note_prefix="Fig. 15b phi1, ",
            )
            for xv in at
        ]
    else:
        mask = mask_interior(rgb_mask(image, RGB_ORANGE, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
        out = [
            sample_solid_curve_y(
                mask, x_axis, y_axis, float(xv), r_min, r_max,
                pick_fn=lambda rs, col: rs[0],
                note_prefix="Fig. 15b phi2, ",
            )
            for xv in at
        ]
    return out


def read_fig16a(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    """Read blue in-phase or green antiphase stationary n = 1 branch in Fig. 16(a)."""
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    xt_m = match_ticks_regular(xt, 3, tol=4.0)
    yt_m = match_ticks_regular(yt, 8, tol=4.0)
    x_axis = verify_axis("linear", list(zip(xt_m, [-1.5, -1.0, -0.5])), max_residual_px=1.5)
    y_axis = verify_axis("linear", list(zip(yt_m, [1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=1.5)

    rgb = RGB_BLUE if branch == "inphase" else RGB_GREEN
    mask = mask_interior(rgb_mask(image, rgb, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=5)
    r_min = ht[2] + 5
    r_max = int(round(y_axis.to_pixel(0.05)))
    out = []
    for xv in at:
        out.append(
            sample_solid_curve_y(
                mask,
                x_axis,
                y_axis,
                float(xv),
                r_min,
                r_max,
                pick_fn=lambda rs, col: rs[0],
                note_prefix=f"Fig. 16a {branch}, ",
            )
        )
    return out
