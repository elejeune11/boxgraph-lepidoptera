"""Shared calibrated reader for SI Fig. S1 (myosin cluster width vs length at 0, 8, 20 min)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_blobs,
    find_box_spines,
    find_ticks,
    make_point_reading,
    match_ticks_regular,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # 6 interior x-ticks at 0.8, 1.0, 1.2, 1.4, 1.6, 1.8; left = 0.6, right = 2.0.
    xt_raw = find_ticks(dark, bot, (left[2] + 10, right[1] - 10), "x", side=-1, band=(5, 25), frac=0.5)
    xt = match_ticks_regular(xt_raw, 6, tol=5.0)
    x_marks = [(left[0], 0.6)] + [(xt[i], 0.8 + 0.2 * i) for i in range(6)] + [(right[0], 2.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # 5 interior y-ticks at 1.2, 1.0, 0.8, 0.6, 0.4; top = 1.4, bot = 0.2.
    yt_raw = find_ticks(dark, left, (top[2] + 10, bot[1] - 10), "y", side=+1, band=(5, 25), frac=0.5)
    yt = match_ticks_regular(yt_raw, 5, tol=5.0)
    y_marks = [(top[0], 1.4)] + [(yt[i], 1.2 - 0.2 * i) for i in range(5)] + [(bot[0], 0.2)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_series(image: np.ndarray, series: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    H, W = image.shape[:2]
    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 15 : bot[1] - 15, left[2] + 15 : right[1] - 15] = True
    # Exclude legend box (top-left)
    roi[: top[2] + 470, : left[2] + 740] = False

    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"

    if series == "0min":
        black = (R < 50) & (G < 50) & (B < 50) & roi
        # Exclude text labels "puncta" and "cluster"
        black[int(y_axis.to_pixel(0.56)) : int(y_axis.to_pixel(0.41)), int(x_axis.to_pixel(0.61)) : int(x_axis.to_pixel(0.86))] = False
        black[int(y_axis.to_pixel(0.65)) : int(y_axis.to_pixel(0.53)), int(x_axis.to_pixel(1.24)) : int(x_axis.to_pixel(1.50))] = False
        blobs = [b for b in find_blobs(black, min_area=1500) if 70 <= b["w"] <= 95]
        pts = []
        for b in blobs:
            if b["h"] > 120:
                # Two vertically touching 82x82 black x's (total height 158 px)
                # Measure exact single-x half-height from isolated blobs (~41 px)
                pts.append((b["c"], b["r0"] + 41.0, f"doublet upper x ({b['w']}x{b['h']}px)"))
                pts.append((b["c"], b["r1"] - 40.0, f"doublet lower x ({b['w']}x{b['h']}px)"))
            else:
                pts.append((b["c"], b["r"], f"black x {b['w']}x{b['h']}px"))
        pts.sort(key=lambda t: (t[0], -t[1]))
        if len(pts) != 6:
            raise RuntimeError(f"Fig. S1 0min: expected 6 black x's, found {len(pts)}: {pts}")
        return [make_point_reading(c, r, x_axis, y_axis, half_px_y=3.5, note=f"{info}, {cal_note}") for c, r, info in pts]

    elif series == "8min":
        blue = (B > 120) & ((B - R) > 40) & ((B - G) > 40) & roi
        blobs = sorted(find_blobs(blue, min_area=1500), key=lambda b: b["c"])
        if len(blobs) != 3:
            raise RuntimeError(f"Fig. S1 8min: expected 3 blue circles, found {len(blobs)}")
        return [
            make_point_reading(b["c"], b["r"], x_axis, y_axis, half_px_y=3.0,
                               note=f"blue circle {b['w']}x{b['h']}px, {cal_note}")
            for b in blobs
        ]

    elif series == "20min":
        red = (R > 160) & (G < 80) & (B < 80) & roi
        blobs = sorted(find_blobs(red, min_area=1500), key=lambda b: b["c"])
        if len(blobs) != 1:
            raise RuntimeError(f"Fig. S1 20min: expected 1 red diamond, found {len(blobs)}")
        b = blobs[0]
        return [
            make_point_reading(b["c"], b["r"], x_axis, y_axis, half_px_y=3.0,
                               note=f"red diamond {b['w']}x{b['h']}px, {cal_note}")
        ]
    else:
        raise ValueError(f"Unknown series {series!r}")
