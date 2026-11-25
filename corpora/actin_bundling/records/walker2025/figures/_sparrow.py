"""Shared calibrated reader for Fig. S1B SPARROW Disorder & Predicted pLDDT profiles
(Lpd aa850-1250 and FL Lpd, 900 dpi).

Calibrates both axes directly from the internal light-grey gridlines inside the plot box
(sub-pixel linearity < 1.5 px across 3,200 px width) and tracks the blue Disorder and
red Predicted pLDDT curves at requested residue indices.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import median_filter
import _common as C


def read_sparrow_panel(
    img: np.ndarray,
    channel: str,
    x_tick_vals: list[float],
    at: list[float],
    label: str = "",
) -> list[dict]:
    rgb = img.astype(np.float64)
    g = rgb.mean(axis=2)
    sat = rgb.max(axis=2) - rgb.min(axis=2)
    H, W = g.shape

    # Locate bottom spine (y = 0.0) and left spine
    row_prof = ((g < 210.0) & (sat < 15.0)).sum(axis=1)
    r_bot = int(np.argmax(row_prof[H // 2 :]) + H // 2)
    col_prof = ((g[150 : r_bot + 1, : W // 4] < 210.0) & (sat[150 : r_bot + 1, : W // 4] < 15.0)).sum(axis=0)
    c_left = int(np.argmax(col_prof))

    sub_bg = (g[200 : r_bot - 10, c_left + 5 : int(0.85 * W)] < 245.0) & (
        sat[200 : r_bot - 10, c_left + 5 : int(0.85 * W)] < 15.0
    )
    c_right = c_left + 5 + int(np.flatnonzero(sub_bg.mean(axis=0) > 0.25)[-1])

    # Horizontal gridlines at y = [1.0, 0.8, 0.6, 0.4, 0.2] + bottom spine y = 0.0
    row_med = np.median(g[150:r_bot, c_left + 20 : c_right - 20], axis=1)
    rdiff = median_filter(row_med, size=25) - row_med
    h_grid = [c[0] for c in C.clusters_1d(np.flatnonzero(rdiff > 2.0) + 150, gap=4)]
    # Match 5 evenly spaced gridlines (1.0, 0.8, 0.6, 0.4, 0.2) + bottom spine
    if len(h_grid) >= 6:
        yt = C.match_regular_ticks(h_grid, 6, tol=8.0)
        ya = C.make_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=4.5)
    else:
        yt5 = C.match_regular_ticks(h_grid, 5, tol=8.0)
        ya = C.make_axis("linear", list(zip(yt5 + [float(r_bot)], [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=5.5)

    # Vertical gridlines at x_tick_vals using neutral-pixel column median + x-tick label cross-check
    sub_g = g[220 : r_bot - 15, c_left + 10 : c_right - 10].copy()
    sub_sat = sat[220 : r_bot - 15, c_left + 10 : c_right - 10]
    sub_g[sub_sat > 15.0] = np.nan
    col_med = np.nanmedian(sub_g[int(0.50 * sub_g.shape[0]) :, :], axis=0)
    col_med = np.where(np.isnan(col_med), np.nanmedian(sub_g, axis=0), col_med)
    cdiff = median_filter(col_med, size=31) - col_med
    vg = [c[0] for c in C.clusters_1d(np.flatnonzero(cdiff > 1.5) + c_left + 10, gap=5)]

    # Use x-tick labels below r_bot as fallback / validation if needed
    n_x = len(x_tick_vals)
    try:
        xt = C.match_regular_ticks(vg, n_x, tol=10.0)
    except ValueError:
        lbl_dark = g[r_bot + 12 : min(H, r_bot + 110), c_left : c_right + 30] < 150.0
        lbl_cls = [c[0] for c in C.clusters_1d(np.flatnonzero(lbl_dark.sum(axis=0) >= 4) + c_left, gap=35)]
        xt = C.match_regular_ticks(lbl_cls, n_x, tol=14.0)
    xa = C.make_axis("linear", list(zip(xt, x_tick_vals)), max_res_px=6.0)

    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    if channel == "disorder":
        # Steel blue curve: B significantly > R, not dark text
        mask = (B - R > 26.0) & (B > 105.0) & (R < 175.0)
    elif channel == "plddt":
        # Bright red curve: R strongly > G and B
        mask = (R - B > 75.0) & (R - G > 65.0) & (R > 140.0)
    else:
        raise ValueError(f"Unknown channel: {channel}")

    # Restrict strictly to plot box interior
    r_top = max(0, int(round(ya.to_pixel(1.04))))
    r_lo = min(H - 1, int(round(ya.to_pixel(-0.02))))
    box_mask = np.zeros_like(mask)
    box_mask[r_top : r_lo + 1, c_left + 2 : c_right - 1] = mask[r_top : r_lo + 1, c_left + 2 : c_right - 1]

    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        c_i = int(round(pc))
        ys = np.nonzero(box_mask[:, max(c_left + 2, c_i - 2) : min(c_right, c_i + 3)].any(axis=1))[0]
        if ys.size == 0:
            ys = np.nonzero(box_mask[:, max(c_left + 2, c_i - 6) : min(c_right, c_i + 7)].any(axis=1))[0]
        if ys.size == 0:
            out.append({"x": float(xv), "y": None, "note": f"{label}: no {channel} pixels near res={xv}"})
            continue
        cls = C.clusters_1d(ys, gap=4)
        best = max(cls, key=lambda t: t[2] - t[1] + 1)
        row_c = best[0]
        half_px = max(2.0, 0.5 * (best[2] - best[1] + 1))

        # Local slope across +/- 6 px
        slope = 0.0
        yl = np.nonzero(box_mask[:, max(c_left + 2, c_i - 6)])[0]
        yr = np.nonzero(box_mask[:, min(c_right - 2, c_i + 6)])[0]
        if yl.size > 0 and yr.size > 0:
            cl_l = min(C.clusters_1d(yl, gap=4), key=lambda t: abs(t[0] - row_c))[0]
            cl_r = min(C.clusters_1d(yr, gap=4), key=lambda t: abs(t[0] - row_c))[0]
            slope = (cl_r - cl_l) / 12.0

        h_px = half_px + 1.0 + abs(slope) * max(xa.residual_px, 1.0)
        out.append(
            C.reading_point(
                x_val=float(xv),
                col=pc,
                row=row_c,
                y_axis=ya,
                half_px=h_px,
                note=f"{label} {channel}: stroke={2*half_px:.1f}px, res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                clamp_min=0.0,
                clamp_max=1.0,
            )
        )
    return out
