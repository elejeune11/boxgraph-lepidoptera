"""Shared calibrated reader for Fig. 2(a)-(f): 1D PDE spatial snapshots at t, t+T/4, t+T/2."""

from __future__ import annotations

import math
import numpy as np

from _calib import (
    find_frame,
    find_inward_ticks,
    sample_solid_curve_y,
    verify_axis,
)


def read_fig2_profile(image: np.ndarray, at: list[float], field: str, panel: str) -> list[dict]:
    """Read one spatial profile curve from a single-panel 1800-dpi crop of Fig. 2.

    Args:
        image: Rendered RGB crop at 1800 dpi of panel (a)-(f).
        at: List of x/ell positions in [0, 2*pi].
        field: 'A' (green), 'I' (red-orange), 'v' (blue), or 'f' (black).
        panel: 'a', 'b', 'c', 'd', 'e', or 'f'.
    """
    dark, _, h_lines, v_lines = find_frame(image, row_frac=0.60, col_frac=0.50, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    _, xt_bot, yt_left, _ = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(6, 20), margin=20, max_width=20, frac=0.55
    )

    if len(xt_bot) != 1:
        raise RuntimeError(f"Fig. 2({panel}): expected 1 bottom x-tick (pi), found {xt_bot}")
    x_axis = verify_axis(
        "linear",
        [(v_left[0], 0.0), (xt_bot[0], math.pi), (v_right[0], 2.0 * math.pi)],
        max_residual_px=1.0,
    )

    if panel in ("a", "c", "e"):
        if len(yt_left) != 2:
            raise RuntimeError(f"Fig. 2({panel}): expected 2 interior left y-ticks (4, 2), found {yt_left}")
        y_axis = verify_axis(
            "linear",
            [(yt_left[0], 4.0), (yt_left[1], 2.0), (h_bot[0], 0.0)],
            max_residual_px=1.0,
        )
    elif panel in ("b", "d", "f"):
        if len(yt_left) != 3:
            raise RuntimeError(f"Fig. 2({panel}): expected 3 left y-ticks (+1, 0, -1), found {yt_left}")
        y_axis = verify_axis(
            "linear",
            [(yt_left[0], 1.0), (yt_left[1], 0.0), (yt_left[2], -1.0)],
            max_residual_px=1.0,
        )
    else:
        raise ValueError(f"Unknown panel {panel!r}")

    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    if field == "A":
        mask = (G > 130) & (G - B > 80) & (G > R + 20)
        rmin, rmax = h_top[2] + 3, h_bot[1] - 2
    elif field == "I":
        mask = (R > 190) & (R - B > 120) & (R > G + 60)
        # Allow reading down to h_bot[2] where I/I0 approaches 0 along the bottom axis
        rmin, rmax = h_top[1], h_bot[2]
    elif field == "v":
        mask = (B > 170) & (B - R > 60)
        rmin, rmax = h_top[2] + 3, h_bot[1] - 3
    elif field == "f":
        mask = (R < 50) & (G < 50) & (B < 50)
        # Active stress curve stays strictly in [0.15, 0.90]: exclude top ticks, panel label, and dashed zero-line
        rmin = h_top[2] + 35
        rmax = int(round(y_axis.to_pixel(0.06)))
        # Mask out panel label "(b)" / "(d)" / "(f)" in upper right corner
        mask[: int(round(y_axis.to_pixel(0.75))), int(round(x_axis.to_pixel(4.7))) :] = False
    else:
        raise ValueError(f"Unknown field {field!r}")

    res = [
        sample_solid_curve_y(
            mask,
            x_axis,
            y_axis,
            float(xv),
            rmin,
            rmax,
            note_prefix=f"Fig. 2({panel}) {field}, ",
        )
        for xv in at
    ]
    if field == "I":
        # Where I/I_0 -> 0 along the bottom axis, the black border (h_bot) occludes the lower half of the red stroke;
        # honestly extend the reading interval down to y = 0 (h_bot[0]) whenever the red stroke touches h_bot[1].
        for row in res:
            if row.get("y") is not None and row.get("px_hi") is not None and row["px_hi"] >= h_bot[1] - 3:
                row["px_hi"] = max(float(row["px_hi"]), float(h_bot[0]))
                row["lo"] = min(float(row["lo"]), 0.0)
                row["note"] += ", clipped by y=0 border"
    return res
