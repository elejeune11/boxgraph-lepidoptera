"""Shared calibrated reader for Fig. 1(a)-(c) (p. 5, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_main_box,
    find_outward_ticks,
    mask_interior,
    match_ticks_regular,
    rgb_mask,
    sample_solid_curve_y,
    verify_axis,
)

# Matplotlib tab10 RGBs observed in Fig. 1
RGB_BLUE = (31.0, 119.0, 180.0)
RGB_ORANGE = (255.0, 127.0, 14.0)
RGB_GREEN = (44.0, 160.0, 44.0)


def calibrate_fig1(image: np.ndarray):
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    xt5 = match_ticks_regular(xt, 5, tol=3.0)
    yt7 = match_ticks_regular(yt, 7, tol=3.0)
    x_vals = [-2.0, -1.0, 0.0, 1.0, 2.0]
    y_vals = [3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0]
    x_axis = verify_axis("linear", list(zip(xt5, x_vals)), max_residual_px=1.0)
    y_axis = verify_axis("linear", list(zip(yt7, y_vals)), max_residual_px=1.0)
    return ht, hb, vl, vr, x_axis, y_axis


def read_fig1_curve(image: np.ndarray, at: list[float], panel: str, curve: str) -> list[dict]:
    ht, hb, vl, vr, x_axis, y_axis = calibrate_fig1(image)
    r_min, r_max = ht[2] + 5, hb[1] - 5

    if curve == "hopf":
        mask = mask_interior(rgb_mask(image, RGB_BLUE, tol=60.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
    elif curve == "turing":
        mask = mask_interior(rgb_mask(image, RGB_ORANGE, tol=65.0, min_chroma=40.0), ht, hb, vl, vr, inset=4)
    elif curve == "ch":
        mask = mask_interior(rgb_mask(image, RGB_GREEN, tol=60.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    # Mask out legend box in Fig. 1(a) (bottom right: f2'' > 0.65, f1'' < -0.95)
    if panel == "1a":
        c_leg = int(round(x_axis.to_pixel(0.65)))
        r_leg = int(round(y_axis.to_pixel(-0.95)))
        mask[r_leg:, c_leg:] = False

    pick_fn = None
    if panel == "1a" and curve == "turing":
        # Upper/solid Turing line for kappa < 1/Q (f2'' < 0); avoid dotted line if any
        pick_fn = lambda rs, col: rs[0]
    elif panel == "1b" and curve == "turing":
        # Lower/solid Turing line for kappa > 1/Q (lies at f1'' < 0, below dotted line at f1'' ~ 0)
        r_zero = y_axis.to_pixel(0.0)
        def pick_lower_orange(rs, col):
            below = [seg for seg in rs if 0.5 * (seg[0] + seg[1]) > r_zero + 20.0]
            return below[-1] if below else rs[-1]
        pick_fn = pick_lower_orange

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
                pick_fn=pick_fn,
                note_prefix=f"Fig. {panel} {curve}, ",
            )
        )
    return out
