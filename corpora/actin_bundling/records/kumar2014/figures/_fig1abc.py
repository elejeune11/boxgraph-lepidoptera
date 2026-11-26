"""Shared calibrated reader for Fig. 1(a)-(c): 1D active advection-diffusion dispersion curves."""

from __future__ import annotations

import numpy as np

from _calib import (
    find_frame,
    find_inward_ticks,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def read_dispersion(image: np.ndarray, at: list[float], series: str, panel: str) -> list[dict]:
    """Read Re[tau*lambda_+] (blue solid) or Im[tau*lambda_+] (red-orange dashed) vs k*ell.

    Args:
        image: Rendered RGB crop at 2400 dpi of panel (a), (b), or (c).
        at: List of k*ell positions to sample.
        series: 're' (solid blue) or 'im' (dashed red-orange).
        panel: 'a', 'b', or 'c' (for diagnostic notes).
    """
    dark, _, h_lines, v_lines = find_frame(image, row_frac=0.45, col_frac=0.35, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    _, xt_bot, yt_left, _ = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(8, 22), margin=25, max_width=20, frac=0.55
    )

    if len(yt_left) != 3:
        raise RuntimeError(f"Fig. 1({panel}): expected 3 left y-ticks (+1, 0, -1), found {yt_left}")
    y_axis = verify_axis(
        "linear",
        [(yt_left[0], 1.0), (yt_left[1], 0.0), (yt_left[2], -1.0)],
        max_residual_px=1.0,
    )

    if len(xt_bot) == 1:
        x_marks = [(v_left[0], 0.0), (xt_bot[0], 1.0), (v_right[0], 2.0)]
    else:
        x_marks = [(v_left[0], 0.0), (v_right[0], 2.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=1.0)

    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    rmin = h_top[2] + 4
    rmax = h_bot[1] - 4

    if series == "re":
        mask = (B - R > 60) & (B > 180)
        return [
            sample_solid_curve_y(
                mask,
                x_axis,
                y_axis,
                float(xv),
                rmin,
                rmax,
                note_prefix=f"Fig. 1({panel}) Re, ",
            )
            for xv in at
        ]
    elif series == "im":
        mask = (R - B > 100) & (R > 180)
        return [
            sample_dashed_curve_y(
                mask,
                x_axis,
                y_axis,
                float(xv),
                rmin,
                rmax,
                max_search_dx=38,
                note_prefix=f"Fig. 1({panel}) Im, ",
            )
            for xv in at
        ]
    raise ValueError(f"Unknown series {series!r}")
