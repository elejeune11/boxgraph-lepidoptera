"""Shared calibrated reader for Fig. 4(d) main panels (p. 8): h(x), m(x), n(x) at t = 0, 1, 2, 3 s."""

from __future__ import annotations

import numpy as np
from _calib import (
    clusters,
    find_inward_ticks,
    find_neutral_frame,
    match_ticks_regular,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_fig4d(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_neutral_frame(
        image, row_frac=0.55, col_frac=0.55, gray_thresh=140.0, sat_thresh=40
    )
    if len(v_lines) < 2 or len(h_lines) < 4:
        raise RuntimeError(f"Fig. 4(d): unexpected frame lines h={h_lines}, v={v_lines}")

    v_left, v_right = v_lines[0], v_lines[-1]
    h0, h1, h2, h3 = h_lines[0], h_lines[1], h_lines[2], h_lines[3]
    c0, c1 = int(round(v_left[0])), int(round(v_right[0]))

    # --- Shared x-axis: x = 0, 0.5, 1.0, 1.5, 2.0, 2.5 um ---
    # Use interior band (9, 20) strictly inside the 12-px-tall rasterized border
    dark_t = (gray < 175.0) & ((image.max(axis=2).astype(int) - image.min(axis=2).astype(int)) < 45)
    xt_t, _, _, _ = find_inward_ticks(
        dark_t, int(round(h0[0])), int(round(h1[0])), c0, c1, band=(9, 20), margin=15, max_width=24, frac=0.40
    )
    x_ticks = match_ticks_regular([v_left[0]] + xt_t + [v_right[0]], 6, tol=8.0)
    x_ax = verify_axis("linear", list(zip(x_ticks, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5])), max_residual_px=4.0)

    # Helper: locate the horizontal light-blue t=0 zero line in a panel [ra..rb]
    r_ch = image[:, :, 0].astype(np.int16)
    b_ch = image[:, :, 2].astype(np.int16)

    def find_zero_line(ra: int, rb: int) -> float:
        blue_diff = (b_ch[ra + 15 : rb - 15, c0 + 30 : c1 - 30] - r_ch[ra + 15 : rb - 15, c0 + 30 : c1 - 30]) > 35
        cs = clusters(np.where(blue_diff.mean(axis=1) > 0.45)[0] + ra + 15)
        if not cs:
            raise RuntimeError(f"Fig. 4(d): could not locate zero line in rows [{ra}, {rb}]")
        return cs[0][0]

    # --- Top panel h(x) left y-axis: 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6 um ---
    _, _, yt_h_l, _ = find_inward_ticks(
        gray < 175.0, int(round(h0[0])), int(round(h1[0])), c0, c1, band=(9, 20), margin=10, max_width=24, frac=0.40
    )
    h_ticks_l = match_ticks_regular(yt_h_l, 7, tol=8.0)
    y_h_left = verify_axis(
        "linear", list(zip(h_ticks_l, [0.0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6])), max_residual_px=5.0
    )

    # --- Top panel h(x) right y-axis (for t = 0): 0, -0.5, -1 um ---
    # 0 is anchored by the left/right aligned zero tick; -0.5 and -1 are read from right-side label centroids
    dark_r = gray[int(h0[0]) + 10 : int(h1[0]) + 40, int(v_right[2]) + 20 :] < 150.0
    r_lbls = clusters(np.where(dark_r.any(axis=1))[0] + int(h0[0]) + 10, gap=20)
    if len(r_lbls) < 3:
        raise RuntimeError(f"Fig. 4(d): could not detect 3 right-axis labels (0, -0.5, -1): {r_lbls}")
    y_h_right = verify_axis(
        "linear", [(h_ticks_l[0], 0.0), (r_lbls[1][0], -0.5), (r_lbls[2][0], -1.0)], max_residual_px=8.0
    )

    # --- Detect inset bounding box inside top panel h(x) so inset curves are masked ---
    # Look for horizontal gray frame lines spanning > 40% of the right half of h(x)
    sub_dark = dark[int(h0[2]) + 100 : int(h1[1]) - 10, int(c0 + 0.35 * (c1 - c0)) : c1 - 10]
    ins_h = clusters(np.where(sub_dark.mean(axis=1) > 0.50)[0] + int(h0[2]) + 100)
    sub_v_dark = dark[int(h0[2]) + 400 : int(h1[1]) - 30, int(c0 + 0.30 * (c1 - c0)) : c1 - 10]
    ins_v = clusters(np.where(sub_v_dark.mean(axis=0) > 0.50)[0] + int(c0 + 0.30 * (c1 - c0)))
    if not ins_h or not ins_v:
        raise RuntimeError(f"Fig. 4(d): could not detect inset frame box: ins_h={ins_h}, ins_v={ins_v}")
    inset_box = (int(ins_h[0][1]) - 6, int(h1[1]), int(ins_v[0][1]) - 6, int(v_right[1]))

    # --- Middle panel m(x) y-axis: +0.5 (top border), 0 (zero line), -0.5 (left tick), -1.0 (bottom border) ---
    m_zero = find_zero_line(int(h1[0]), int(h2[0]))
    _, _, yt_m_l, _ = find_inward_ticks(
        dark_t, int(round(h1[0])), int(round(h2[0])), c0, c1, band=(9, 20), margin=15, max_width=24, frac=0.40
    )
    # Pick the tick closest to expected -0.5 position (m_zero + (m_zero - h1[0]))
    exp_m05 = m_zero + (m_zero - h1[0])
    m_tick_05 = min(yt_m_l, key=lambda r: abs(r - exp_m05))
    y_m = verify_axis(
        "linear", [(h1[0], 0.5), (m_zero, 0.0), (m_tick_05, -0.5), (h2[0], -1.0)], max_residual_px=4.0
    )

    # --- Bottom panel n(x) y-axis: +1.0 (top border), +0.5, 0, -0.5, -1.0 (ticks), -1.5 (bottom border) ---
    n_zero = find_zero_line(int(h2[0]), int(h3[0]))
    _, _, yt_n_l, _ = find_inward_ticks(
        dark_t, int(round(h2[0])), int(round(h3[0])), c0, c1, band=(9, 20), margin=15, max_width=24, frac=0.40
    )
    n_cands = [h2[0], n_zero, h3[0]] + yt_n_l
    n_ticks = match_ticks_regular(n_cands, 6, tol=8.0)
    y_n = verify_axis(
        "linear", list(zip(n_ticks, [1.0, 0.5, 0.0, -0.5, -1.0, -1.5])), max_residual_px=4.0
    )

    boxes = {
        "h": (int(h0[2]) + 2, int(h1[1]) - 2, int(v_left[2]) + 2, int(v_right[1]) - 2),
        "m": (int(h1[2]) + 2, int(h2[1]) - 2, int(v_left[2]) + 2, int(v_right[1]) - 2),
        "n": (int(h2[2]) + 2, int(h3[1]) - 2, int(v_left[2]) + 2, int(v_right[1]) - 2),
    }
    return x_ax, y_h_left, y_h_right, y_m, y_n, boxes, inset_box


def _color_masks(image: np.ndarray) -> dict[str, np.ndarray]:
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    # Distinct RGB prototypes measured on crop_fig4d:
    #   cyan  (t=0): ~[99, 186, 201]  (g > 155, b > 160, g - r > 45, b - g < 22)
    #   blue  (t=1): ~[68, 145, 198] .. [153, 197, 232] (b > 140, b - r > 45, b - g >= 24)
    #   green (t=2): ~[138, 157, 148] .. [178, 198, 189] (g > 110, g - r >= 10, g - b >= 4, r > 75, b < 215)
    #   red   (t=3): ~[242, 129, 145] (r > 175, r - g > 45, r - b > 40)
    return {
        "t0": (b > 160) & (g > 155) & (g - r > 45) & (b - g < 22),
        "t1": (b > 140) & (b - r > 45) & (b - g >= 24),
        "t2": (g > 110) & (g - r >= 10) & (g - b >= 4) & (r > 75) & (b < 215),
        "t3": (r > 175) & (r - g > 45) & (r - b > 40),
    }


def read_fig4d(image: np.ndarray, at: list[float], panel: str, time_key: str) -> list[dict]:
    x_ax, y_h_left, y_h_right, y_m, y_n, boxes, inset_box = _calibrate_fig4d(image)
    masks = _color_masks(image)
    mask = masks[time_key].copy()

    rmin, rmax, cmin, cmax = boxes[panel]
    if panel == "h":
        ir0, ir1, ic0, ic1 = inset_box
        mask[ir0 : ir1 + 1, ic0 : ic1 + 1] = False
        y_ax = y_h_right if time_key == "t0" else y_h_left
    elif panel == "m":
        y_ax = y_m
    elif panel == "n":
        y_ax = y_n
    else:
        raise ValueError(f"Unknown Fig. 4(d) panel: {panel}")

    pick_fn = None
    if time_key == "t1":
        # The horizontal zero reference line shares blue ink with t = 1 s (especially in h(x));
        # whenever the wave profile lies off zero, prefer the run outside the zero-line band.
        z_row = float(y_ax.to_pixel(0.0))

        def _pick_off_zero(rs, c):
            off = [seg for seg in rs if abs(0.5 * (seg[0] + seg[1]) + 0.5 - z_row) > 14.0]
            if off:
                return max(off, key=lambda t: t[1] - t[0])
            return rs[0]

        pick_fn = _pick_off_zero

    return [
        sample_solid_curve_y(
            mask, x_ax, y_ax, x, rmin, rmax,
            pick_fn=pick_fn,
            search_dx=8,
            note_prefix=f"fig4d {panel}({time_key}); ",
        )
        for x in at
    ]
