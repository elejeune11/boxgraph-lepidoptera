"""Calibrated reader for Fig. S6B (Cy5-fascin intensity vs [TMR-Fimbrin] competition curve, IC50 = 47.4 nM)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 40), col_range=(120, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 0.8, 0.6, 0.4, 0.2, 0.0 (1.0 lies near row ~312, touching the x=0 open circle)
    yt_raw = C.find_ticks(dark, y_spine, (550, int(x_spine[1]) - 80), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 5, tol=8.0, min_step=280.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=2.5)

    # Inward x ticks: 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 uM
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 15, image.shape[1] - 50), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 7, tol=8.0, min_step=220.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _find_open_circles(dark: np.ndarray, x_spine, y_spine) -> list[tuple[float, float, float]]:
    """Fast detection of all 8 open circular markers (including those bisected by the black fit curve!)."""
    # Downsample 2x for super-fast hole filling (< 0.15 s!)
    sub = dark[::2, ::2].copy()
    sub[:100, :] = False
    sub[int(x_spine[1] // 2) - 5 :, :] = False
    sub[:, : int(y_spine[2] // 2) + 2] = False
    # Exclude 'IC50 = 47.4 nM TMR-Fimbrin' text region
    sub[375:750, 525:] = False

    closed = ndimage.binary_dilation(sub, iterations=1)
    filled = ndimage.binary_fill_holes(closed)
    holes = filled & (~closed)
    lbl, n = ndimage.label(holes)
    frags = []
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        area = len(rr)
        # At 2x downsampling, a full circle interior has area ~600..1600 px^2; a bisected half-moon has ~100..600 px^2!
        if 70 <= area <= 2000:
            frags.append((float(cc.mean()) * 2.0, float(rr.mean()) * 2.0, area, rr, cc))

    # Merge pairs of half-moon fragments belonging to the same bisected circle (centers within 75 px!)
    used = [False] * len(frags)
    circles = []
    for i in range(len(frags)):
        if used[i]:
            continue
        group_rr = [frags[i][3]]
        group_cc = [frags[i][4]]
        tot_area = frags[i][2]
        used[i] = True
        for j in range(i + 1, len(frags)):
            if used[j]:
                continue
            if np.hypot(frags[i][0] - frags[j][0], frags[i][1] - frags[j][1]) <= 75.0:
                group_rr.append(frags[j][3])
                group_cc.append(frags[j][4])
                tot_area += frags[j][2]
                used[j] = True
        if tot_area < 280:
            continue
        all_r = np.concatenate(group_rr) * 2.0
        all_c = np.concatenate(group_cc) * 2.0
        h = all_r.max() - all_r.min() + 2.0
        w = all_c.max() - all_c.min() + 2.0
        if 0.65 <= (h / max(w, 1.0)) <= 1.55 and 35.0 <= h <= 115.0 and 35.0 <= w <= 115.0:
            c_mid = 0.5 * (all_c.min() + all_c.max()) + 0.5
            r_mid = 0.5 * (all_r.min() + all_r.max()) + 0.5
            rad = 0.25 * (h + w) + 12.0
            circles.append((float(c_mid), float(r_mid), float(rad)))
    circles.sort(key=lambda t: t[0])
    return circles


def read_series(image: np.ndarray, kind: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    dark = image.mean(axis=2) < 125
    circles = _find_open_circles(dark, x_spine, y_spine)

    if kind == "data":
        out = []
        for c_mid, r_mid, _ in circles:
            xv = float(x_axis.to_value(c_mid))
            if abs(xv) < 0.008:
                xv = 0.0
            out.append(
                C.make_point_reading(
                    xv,
                    c_mid,
                    r_mid,
                    y_axis,
                    err_px=2.5 + y_axis.residual_px,
                    note=f"Fig. S6B open circle marker; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                )
            )
        return out

    if kind == "fit":
        H, W = dark.shape
        curve = dark.copy()
        row_min = 120
        # Stop strictly ABOVE inward x-ticks (which protrude upward from x_spine at y ~ -0.10 to -0.05!)
        row_max = int(round(y_axis.to_pixel(-0.04)))
        col_min = int(y_spine[2]) + 15
        curve[:row_min, :] = False
        curve[row_max + 1 :, :] = False
        curve[:, :col_min] = False
        curve[750:1500, 1050:] = False

        # Mask out the 8 open circle rings so the curve sampler never snaps onto a circle rim
        yy_g, xx_g = np.ogrid[:H, :W]
        for c_m, r_m, rad in circles:
            r_out = rad + 6.0
            y0, y1 = max(0, int(r_m - r_out) - 2), min(H, int(r_m + r_out) + 3)
            x0, x1 = max(0, int(c_m - r_out) - 2), min(W, int(c_m + r_out) + 3)
            d2 = (yy_g[y0:y1, :] - r_m) ** 2 + (xx_g[:, x0:x1] - c_m) ** 2
            curve[y0:y1, x0:x1] &= ~(d2 <= r_out**2)

        return [
            C.sample_solid_curve_y(
                curve,
                x_axis,
                y_axis,
                float(xv),
                row_min,
                row_max,
                search_dx=20,
                note_prefix="Fig. S6B IC50 fit curve; ",
            )
            for xv in at
        ]
    raise ValueError(f"Unknown kind {kind!r}")
