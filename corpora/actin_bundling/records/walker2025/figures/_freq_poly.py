"""Shared calibrated reader for droplet diameter frequency polygon plots (Figs. 1I and 2D)."""

from __future__ import annotations

import numpy as np
import _common as C


def read_freq_curve(
    img: np.ndarray,
    x_tick_vals: list[float],
    y_tick_vals: list[float],
    color_kind: str,
    at: list[float],
    label: str = "",
) -> list[dict]:
    """Read one colored frequency polygon curve at requested x values (diameters in um)."""
    ax = C.find_l_axes(img, dark_thresh=115.0, min_h_frac=0.35, min_v_frac=0.5, max_tick_len=65)
    xt = C.match_regular_ticks(ax["x_ticks"], len(x_tick_vals), tol=4.5)
    yt = C.match_regular_ticks(ax["y_ticks"], len(y_tick_vals), tol=4.5)
    xa = C.make_axis("linear", list(zip(xt, x_tick_vals)), max_res_px=4.0)
    ya = C.make_axis("linear", list(zip(yt, y_tick_vals)), max_res_px=4.0)

    rgb = img.astype(np.float64)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)

    if color_kind == "blue":
        mask = (B - R > 80.0) & (B - G > 80.0) & (B > 130.0)
    elif color_kind == "red":
        mask = (R - G > 80.0) & (R - B > 80.0) & (R > 130.0)
    elif color_kind == "grey":
        mask = (gray > 65.0) & (gray < 155.0) & (sat < 20.0)
    elif color_kind == "black":
        mask = (gray < 55.0) & (sat < 20.0)
    else:
        raise ValueError(f"Unknown color_kind: {color_kind}")

    # Mask out legend region in upper right and axes spines
    H, W = mask.shape
    mask = mask.copy()
    r_zero = float(ya.to_pixel(0.0))
    c_left = ax["left"]["hi"] + 2
    mask[:, :c_left] = False
    mask[ax["bot"]["lo"] - 1 :, :] = False
    # Exclude upper-right legend region cleanly without touching curves:
    # In Fig. 1I (y_max=500): legend is x >= 4.4, y >= 220
    # In Fig. 2D (y_max=160): legend is x >= 3.8, y >= 138 (all curves right of x=3.8 stay below 135)
    if max(y_tick_vals) <= 200.0:
        c_leg = int(round(xa.to_pixel(3.8)))
        r_leg = int(round(ya.to_pixel(137.0)))
    else:
        c_leg = int(round(xa.to_pixel(4.4)))
        r_leg = int(round(ya.to_pixel(220.0)))
    mask[: max(0, r_leg), max(0, c_leg) :] = False

    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        col = int(round(pc))
        # Check narrow column window [-2, +2] for the curve
        ys = np.nonzero(mask[:, max(c_left, col - 1) : min(W, col + 2)].any(axis=1))[0]
        if len(ys) == 0:
            # If no pixels above the baseline, curve lies flush on the y=0 baseline!
            # Verify there is curve ink on or near the baseline (or report 0.0 explicitly)
            out.append(
                C.reading_point(
                    x_val=float(xv),
                    col=pc,
                    row=r_zero,
                    y_axis=ya,
                    half_px=3.0,
                    note=f"{label}: flush with baseline y=0, res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                    clamp_min=0.0,
                )
            )
            continue
        cs = C.clusters_1d(ys, gap=4)
        # Prefer run at exact column if available
        ys_ex = np.nonzero(mask[:, col])[0]
        if len(ys_ex) > 0:
            cs_ex = C.clusters_1d(ys_ex, gap=4)
            best = max(cs_ex, key=lambda t: t[2] - t[1])
        else:
            best = max(cs, key=lambda t: t[2] - t[1])
        y_row = best[0]
        half_px = 0.5 * (best[2] - best[1] + 1)

        # Estimate local slope over +/- 4 px
        slope = 0.0
        yl = np.nonzero(mask[:, max(c_left, col - 4)])[0]
        yr = np.nonzero(mask[:, min(W - 1, col + 4)])[0]
        if len(yl) > 0 and len(yr) > 0:
            c_l = min(C.clusters_1d(yl, gap=4), key=lambda t: abs(t[0] - y_row))[0]
            c_r = min(C.clusters_1d(yr, gap=4), key=lambda t: abs(t[0] - y_row))[0]
            slope = (c_r - c_l) / 8.0

        h_px = half_px + 0.5 + abs(slope) * max(xa.residual_px, 0.5)
        out.append(
            C.reading_point(
                x_val=float(xv),
                col=pc,
                row=y_row,
                y_axis=ya,
                half_px=h_px,
                note=f"{label}: stroke={2*half_px:.1f}px, res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                clamp_min=0.0,
            )
        )
    return out
