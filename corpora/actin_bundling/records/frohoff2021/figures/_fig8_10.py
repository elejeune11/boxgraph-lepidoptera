"""Shared calibrated reader for Figs. 8, 9, and 10(a) (pages 13-15, 1200 dpi)."""

from __future__ import annotations

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
RGB_ORANGE = (255.0, 127.0, 14.0)
RGB_RED = (214.0, 39.0, 40.0)
RGB_PURPLE = (148.0, 103.0, 189.0)
RGB_OLIVE = (188.0, 189.0, 34.0)
RGB_CYAN = (23.0, 190.0, 207.0)


def calibrate_fig8_10(image: np.ndarray, panel: str):
    dark, gray, ht, hb, vl, vr = find_main_box(image, row_frac=0.35, col_frac=0.35)
    xt, yt = find_outward_ticks(dark, ht, hb, vl, vr, band=(4, 12), max_width=14)
    specs = {
        "8a": ([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], [1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]),
        "8b": ([0.4, 0.6, 0.8, 1.0, 1.2, 1.4], [1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]),
        "8c": ([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8], [1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]),
        "8d": ([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0], [0.8, 0.6, 0.4, 0.2, 0.0]),
        "9":  ([0.0, 2.0, 4.0, 6.0], [1.0, 0.75, 0.5, 0.25, 0.0]),
        "10a": ([-20.0, -15.0, -10.0, -5.0, 0.0], [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]),
    }
    x_vals, y_vals = specs[panel]
    xt_m = match_ticks_regular(xt, len(x_vals), tol=3.5)
    yt_m = match_ticks_regular(yt, len(y_vals), tol=3.5)
    x_axis = verify_axis("linear", list(zip(xt_m, x_vals)), max_residual_px=1.0)
    y_axis = verify_axis("linear", list(zip(yt_m, y_vals)), max_residual_px=1.0)
    return dark, ht, hb, vl, vr, x_axis, y_axis


def read_fig8_10_curve(
    image: np.ndarray,
    at: list[float],
    panel: str,
    color_name: str,
    style: str = "solid",
    dashed_above_x: float | None = None,
    y_min_cut: float = 0.03,
) -> list[dict]:
    dark, ht, hb, vl, vr, x_axis, y_axis = calibrate_fig8_10(image, panel)
    r_min, r_max = ht[2] + 5, min(hb[1] - 5, int(round(y_axis.to_pixel(y_min_cut))))

    rgb_map = {
        "blue": RGB_BLUE,
        "green": RGB_GREEN,
        "orange": RGB_ORANGE,
        "red": RGB_RED,
        "purple": RGB_PURPLE,
        "olive": RGB_OLIVE,
        "cyan": RGB_CYAN,
    }
    mask = mask_interior(rgb_mask(image, rgb_map[color_name], tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)

    # Mask insets so inset curves can never interfere with main panel curves
    if panel == "8a":
        c0, c1 = int(x_axis.to_pixel(0.24)), int(x_axis.to_pixel(0.92))
        r0, r1 = int(y_axis.to_pixel(0.82)), int(y_axis.to_pixel(0.03))
        mask[r0:r1, c0:c1] = False
    elif panel == "8b":
        c0, c1 = int(x_axis.to_pixel(0.32)), int(x_axis.to_pixel(0.96))
        r0, r1 = int(y_axis.to_pixel(0.82)), int(y_axis.to_pixel(0.03))
        mask[r0:r1, c0:c1] = False
    elif panel == "8d":
        c0, c1 = int(x_axis.to_pixel(0.28)), int(x_axis.to_pixel(0.92))
        r0, r1 = int(y_axis.to_pixel(0.75)), int(y_axis.to_pixel(0.46))
        mask[r0:r1, c0:c1] = False
    elif panel == "10a":
        c0, c1 = int(x_axis.to_pixel(-22.5)), int(x_axis.to_pixel(-11.8))
        r0, r1 = int(y_axis.to_pixel(3.3)), int(y_axis.to_pixel(0.4))
        mask[r0:r1, c0:c1] = False

    out = []
    for xv in at:
        xf = float(xv)
        r_lo, r_hi = r_min, r_max
        if panel == "8a" and color_name == "blue":
            # Stay on upper stable branch above the saddle-node fold (norm > 0.35)
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.38))))
        elif panel == "8b" and color_name == "blue":
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.38))))
        elif panel == "8b" and color_name == "green":
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.45))))
        elif panel == "8c" and color_name == "blue":
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.35))))
        elif panel == "8c" and color_name == "green":
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.45))))
        elif panel == "8c" and color_name == "orange":
            r_hi = min(r_max, int(round(y_axis.to_pixel(0.35))))

        use_dashed = (style == "dashed") or (dashed_above_x is not None and xf >= dashed_above_x)
        if use_dashed:
            out.append(
                sample_dashed_curve_y(
                    mask, x_axis, y_axis, xf, r_lo, r_hi,
                    max_search_dx=42, note_prefix=f"Fig. {panel} {color_name}, "
                )
            )
        else:
            out.append(
                sample_solid_curve_y(
                    mask, x_axis, y_axis, xf, r_lo, r_hi,
                    pick_fn=lambda rs, col: rs[0],
                    note_prefix=f"Fig. {panel} {color_name}, "
                )
            )
    return out


def read_fig9_subcrit_plus(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read lower boundary Delta_+(Sigma) of the blue n+ subcriticality region in Fig. 9."""
    dark, ht, hb, vl, vr, x_axis, y_axis = calibrate_fig8_10(image, "9")
    blue = mask_interior(rgb_mask(image, RGB_BLUE, tol=55.0, min_chroma=35.0), ht, hb, vl, vr, inset=4)
    # Restrict above insets (Delta > 0.63) and below Hopf line (Delta < 0.90)
    r_top_band = int(round(y_axis.to_pixel(0.90)))
    r_bot_band = int(round(y_axis.to_pixel(0.63)))
    out = []
    for xv in at:
        xf = float(xv)
        c_f = float(x_axis.to_pixel(xf))
        c = int(round(c_f))
        # Find bottom-most blue pixel and the black boundary stroke immediately below it
        b_rows = np.nonzero(blue[r_top_band:r_bot_band, c])[0]
        if len(b_rows) == 0:
            out.append({"x": xf, "y": None, "found": False, "note": f"no blue fill at col {c}"})
            continue
        r_blue_bot = r_top_band + int(b_rows.max())
        # Look for black boundary line right at/below r_blue_bot
        d_sub = dark[r_blue_bot - 2 : r_blue_bot + 16, c]
        rs = runs_1d(d_sub)
        if rs:
            r0, r1 = rs[0][0] + r_blue_bot - 2, rs[0][1] + r_blue_bot - 2
            r_mid = 0.5 * (r0 + r1) + 0.5
            thick = float(r1 - r0 + 1)
        else:
            r_mid = float(r_blue_bot) + 1.0
            thick = 6.0
        half_px = max(0.5 * thick, 2.0) + 0.5 + y_axis.residual_px + 0.5 * x_axis.residual_px
        y_val = float(y_axis.to_value(r_mid))
        ya_v, yb_v = float(y_axis.to_value(r_mid - half_px)), float(y_axis.to_value(r_mid + half_px))
        out.append({
            "x": xf, "y": y_val, "lo": min(ya_v, yb_v), "hi": max(ya_v, yb_v),
            "px": [c_f, r_mid], "px_lo": r_mid - half_px, "px_hi": r_mid + half_px,
            "note": f"Fig. 9 blue boundary, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, thick={thick:.0f}px",
        })
    return out


def read_fig9_path_rho03(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    """Read dashed black rho = 0.3 coupling path Delta(Sigma) in Fig. 9 [Eq. (23)].

    branch='pos': alpha > -rho (right curve, Sigma in [0.4, 1.75], lies BELOW diagonal Delta = Sigma).
    branch='neg': alpha < -rho (left curve, Sigma in [0.05, 0.50], lies ABOVE diagonal Delta = Sigma).
    """
    from scipy import ndimage
    dark, ht, hb, vl, vr, x_axis, y_axis = calibrate_fig8_10(image, "9")
    blk = mask_interior(dark, ht, hb, vl, vr, inset=6)
    r_hopf = int(round(y_axis.to_pixel(0.87)))
    r_zero = int(round(y_axis.to_pixel(0.03)))
    c_right = int(round(x_axis.to_pixel(1.80)))
    blk[:r_hopf, :] = False
    blk[r_zero:, :] = False
    blk[:, c_right:] = False
    gray = image.astype(float).mean(axis=2)
    chroma = image.max(axis=2).astype(float) - image.min(axis=2).astype(float)
    blk &= (gray < 95.0) & (chroma < 30.0)
    if branch == "neg":
        # Remove the solid black border of the orange shaded region
        orange = mask_interior(rgb_mask(image, RGB_ORANGE, tol=65.0, min_chroma=40.0), ht, hb, vl, vr, inset=2)
        blk &= ~ndimage.binary_dilation(orange, iterations=18)
    H, W = blk.shape
    cc_grid, rr_grid = np.meshgrid(np.arange(W), np.arange(H))
    sig_grid = x_axis.to_value(cc_grid + 0.5)
    del_grid = y_axis.to_value(rr_grid + 0.5)
    if branch == "pos":
        blk &= (del_grid < sig_grid - 0.12)
    else:
        blk &= (del_grid > sig_grid + 0.10)

    out = []
    for xv in at:
        out.append(
            sample_dashed_curve_y(
                blk, x_axis, y_axis, float(xv), r_hopf, r_zero,
                max_search_dx=35, note_prefix=f"Fig. 9 rho=0.3 {branch}, "
            )
        )
    return out
