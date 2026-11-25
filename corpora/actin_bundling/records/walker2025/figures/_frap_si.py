"""Shared calibrated reader for Fig. S3B, S3D, S3F, S3H 2-channel FRAP recovery curves
(protein channel in teal-blue + actin channel in red, with crisp black +/- SD error bars, 900 dpi).
"""

from __future__ import annotations

import numpy as np
import _common as C
from fig1f_frap_minilpd import _read_frap_sd_at_col


def read_frap_si(
    image: np.ndarray,
    channel: str,
    x_tick_vals: list[float],
    at: list[float],
    label: str = "",
) -> list[dict]:
    ax = C.find_l_axes(image, dark_thresh=115.0, min_h_frac=0.45, min_v_frac=0.45, max_tick_len=40)
    xt = C.match_regular_ticks(ax["x_ticks"], len(x_tick_vals), tol=6.0)
    yt = C.match_regular_ticks(ax["y_ticks"], 6, tol=6.0)
    xa = C.make_axis("linear", list(zip(xt, x_tick_vals)), max_res_px=4.0)
    ya = C.make_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=4.0)

    rgb = image.astype(np.float64)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)

    if channel == "protein":
        # Teal-blue Lpd / Eps15 curve: B and G both exceed R clearly
        curve = (B - R > 28.0) & (G - R > 15.0) & (B > 85.0) & (R < 145.0)
    elif channel == "actin":
        # Bright red Actin curve: R strongly exceeds G and B
        curve = (R - G > 75.0) & (R - B > 75.0) & (R > 150.0)
    else:
        raise ValueError(f"Unknown channel: {channel}")

    # Mask upper-right legend region
    x_max = max(x_tick_vals)
    c_leg = int(round(xa.to_pixel(0.65 * x_max)))
    r_leg = int(round(ya.to_pixel(0.74)))
    curve[:r_leg, c_leg:] = False

    # Pure black error bars (exclude axis spines by working strictly right of left spine)
    black = (gray < 115.0) & (sat < 30.0)
    black[:, : ax["left"]["hi"] + 3] = False
    black[:r_leg, c_leg:] = False

    r_top = max(0, int(round(ya.to_pixel(1.06))))
    r_bot = min(image.shape[0] - 1, ax["bot"]["lo"] + 12)

    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        rec = _read_frap_sd_at_col(
            image,
            pc,
            curve,
            black,
            r_top,
            r_bot,
            ya,
            float(xv),
            xa,
            note_prefix=f"{label} {channel} FRAP; ",
        )
        # Clamp normalized intensity values and bars to [0, 1.05] (slightly >1 allowed for pre-bleach SD)
        if rec.get("y") is not None:
            rec["y"] = float(max(0.0, rec["y"]))
            rec["lo"] = float(max(0.0, rec["lo"]))
            if rec.get("bar_lo") is not None:
                rec["bar_lo"] = float(max(0.0, rec["bar_lo"]))
                rec["bar_lo_outer"] = float(max(0.0, rec["bar_lo_outer"]))
        out.append(rec)
    return out
