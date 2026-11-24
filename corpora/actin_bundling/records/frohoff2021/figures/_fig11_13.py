"""Shared calibrated reader for Figs. 11(a)-(f) and 13(a)-(c) (pages 16, 18, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from _calib import (
    clusters,
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


def calibrate_fig11_13(image: np.ndarray, panel: str):
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    if panel == "11f":
        # Inward ticks along bottom and left borders
        in_xt = [
            c[0]
            for c in clusters(
                np.where(dark[hb[1] - 14 : hb[1] - 3, vl[2] + 10 : vr[1] - 10].mean(axis=0) > 0.6)[0]
                + vl[2]
                + 10
            )
            if (c[2] - c[1] + 1) <= 14
        ]
        in_yt = [
            c[0]
            for c in clusters(
                np.where(dark[ht[2] + 10 : hb[1] - 10, vl[2] + 3 : vl[2] + 14].mean(axis=1) > 0.6)[0]
                + ht[2]
                + 10
            )
            if (c[2] - c[1] + 1) <= 14
        ]
        x_vals = [-0.70, -0.65, -0.60, -0.55, -0.50, -0.45]
        y_vals = [1.3125, 1.3100, 1.3075, 1.3050, 1.3025, 1.3000]
        xt_m = match_ticks_regular(in_xt, len(x_vals), tol=3.5)
        yt_m = match_ticks_regular(in_yt, len(y_vals), tol=3.5)
    else:
        xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
        x_vals = [-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2]
        y_vals = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
        xt_m = match_ticks_regular(xt, len(x_vals), tol=3.5)
        yt_m = match_ticks_regular(yt, len(y_vals), tol=3.5)

    x_axis = verify_axis("linear", list(zip(xt_m, x_vals)), max_residual_px=1.0)
    y_axis = verify_axis("linear", list(zip(yt_m, y_vals)), max_residual_px=1.0)
    return dark, ht, hb, vl, vr, x_axis, y_axis


def read_fig11_13_curve(
    image: np.ndarray,
    at: list[float],
    panel: str,
    color_name: str,
    style: str = "solid",
    y_min_cut: float = 0.04,
) -> list[dict]:
    dark, ht, hb, vl, vr, x_axis, y_axis = calibrate_fig11_13(image, panel)
    r_min = ht[2] + 5
    r_max = min(hb[1] - 5, int(round(y_axis.to_pixel(y_min_cut))))

    rgb = RGB_BLUE if color_name == "blue" else RGB_GREEN
    mask = mask_interior(rgb_mask(image, rgb, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)

    # Mask insets in 11b, 11c, 11d
    if panel == "11b":
        mask[int(y_axis.to_pixel(0.98)) : int(y_axis.to_pixel(0.49)), int(x_axis.to_pixel(-0.50)) : int(x_axis.to_pixel(-0.33))] = False
    elif panel == "11c":
        mask[int(y_axis.to_pixel(0.98)) : int(y_axis.to_pixel(0.45)), int(x_axis.to_pixel(-0.55)) : int(x_axis.to_pixel(-0.24))] = False
    elif panel == "11d":
        mask[int(y_axis.to_pixel(0.98)) : int(y_axis.to_pixel(0.48)), int(x_axis.to_pixel(-0.55)) : int(x_axis.to_pixel(-0.24))] = False

    out = []
    for xv in at:
        xf = float(xv)
        if style == "dashed":
            out.append(
                sample_dashed_curve_y(
                    mask, x_axis, y_axis, xf, r_min, r_max,
                    max_search_dx=35, note_prefix=f"Fig. {panel} {color_name}, "
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


def read_fig11f_locus(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read horizontal a(alpha) locus in Fig. 11(f) for curve in ('sn', 'dp_left', 'dp_right').

    `at` holds alpha values (vertical axis); returned `y` is the horizontal a coordinate.
    Uses local quadratic regression c(r) on the isolated dotted curve + marker pixels.
    """
    dark, ht, hb, vl, vr, x_axis, y_axis = calibrate_fig11_13(image, "11f")
    rgb = RGB_BLUE if curve == "dp_right" else RGB_GREEN
    mask = mask_interior(rgb_mask(image, rgb, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=6)

    # Split left vs right curve of the same color cleanly using known topological separation:
    # Green: left = dp_left (triangles), right = sn (crosses)
    # Blue:  left = hopf (diamonds), right = dp_right (triangles)
    c_min, c_max = vl[2] + 6, vr[1] - 6
    H, W = mask.shape
    out = []
    for alpha_v in at:
        af = float(alpha_v)
        r_f = float(y_axis.to_pixel(af))
        r0 = int(round(r_f))
        r_lo = max(ht[2] + 6, r0 - 30)
        r_hi = min(hb[1] - 6, r0 + 30)

        # In a thin horizontal band +/- 18 px around r0, find horizontal clusters to isolate left vs right branch
        rb0, rb1 = max(ht[2] + 6, r0 - 18), min(hb[1] - 6, r0 + 18)
        col_proj = mask[rb0 : rb1 + 1, c_min : c_max + 1].any(axis=0)
        bands = [(a + c_min, b + c_min) for a, b in runs_1d(col_proj, gap=24)]
        if not bands:
            out.append({"x": af, "y": None, "found": False, "note": f"no pixels near alpha={af}"})
            continue

        if curve == "dp_left":
            b0, b1 = bands[0]
            if len(bands) >= 2:
                next_mid = int(round(0.5 * (b1 + bands[1][0])))
                w_lo, w_hi = max(c_min, b0 - 45), min(next_mid - 2, b1 + 45)
            else:
                w_lo, w_hi = max(c_min, b0 - 45), min(c_max, b1 + 45)
        else:
            # 'sn' (rightmost green) or 'dp_right' (rightmost blue)
            b0, b1 = bands[-1]
            if len(bands) >= 2:
                prev_mid = int(round(0.5 * (bands[-2][1] + b0)))
                w_lo, w_hi = max(prev_mid + 2, b0 - 45), min(c_max, b1 + 45)
            else:
                w_lo, w_hi = max(c_min, b0 - 45), min(c_max, b1 + 45)

        sub = mask[r_lo : r_hi + 1, w_lo : w_hi + 1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 10:
            out.append({"x": af, "y": None, "found": False, "note": f"too few pixels for {curve} at alpha={af}"})
            continue
        ys_r = rr.astype(float) + r_lo + 0.5
        xs_c = cc.astype(float) + w_lo + 0.5
        w = np.exp(-0.5 * ((ys_r - r_f) / 18.0) ** 2)
        p1 = np.poly1d(np.polyfit(ys_r, xs_c, 1, w=w))
        resid1 = np.abs(xs_c - p1(ys_r))
        keep = resid1 <= max(10.0, 2.0 * float(np.median(resid1)) + 4.0)
        ys_r, xs_c, w = ys_r[keep], xs_c[keep], w[keep]
        p_lin = np.poly1d(np.polyfit(ys_r, xs_c, 1, w=w))
        p_quad = np.poly1d(np.polyfit(ys_r, xs_c, 2, w=w))
        c_est = float(p_quad(r_f))
        slope = float(np.polyder(p_quad)(r_f))
        resid = np.abs(xs_c - p_quad(ys_r))
        half_thick = max(float(np.percentile(resid, 95)), 4.0)
        curv_diff = 0.25 * abs(c_est - float(p_lin(r_f)))
        half_px = half_thick + 1.0 + curv_diff + x_axis.residual_px + abs(slope) * (0.5 + y_axis.residual_px)
        a_val = float(x_axis.to_value(c_est))
        a_lo, a_hi = float(x_axis.to_value(c_est - half_px)), float(x_axis.to_value(c_est + half_px))
        out.append({
            "x": af,
            "y": a_val,
            "lo": min(a_lo, a_hi),
            "hi": max(a_lo, a_hi),
            "px": [c_est, r_f],
            "note": f"Fig. 11f {curve}, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, half={half_px:.1f}px",
        })
    return out
