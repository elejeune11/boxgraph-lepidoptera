"""Shared calibrated reader for Main Fig. 2B (myosin spatial probability distribution)."""

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

    # Top spine has 4 unobstructed ticks at 0.2, 0.4, 0.6, 0.8; left = 0.0, right = 1.0.
    xt_raw = find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5)
    xt = match_ticks_regular(xt_raw, 4, tol=5.0)
    x_marks = [(left[0], 0.0)] + [(xt[i], 0.2 * (i + 1)) for i in range(4)] + [(right[0], 1.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Right spine has 6 unobstructed ticks at 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 (left tick 0.3 is occluded by blue x).
    yt_raw = find_ticks(dark, right, (top[2] + 10, bot[1] - 10), "y", side=-1, band=(5, 25), frac=0.5)
    yt = match_ticks_regular(yt_raw, 6, tol=5.0)
    y_marks = [(top[0], 0.7)] + [(yt[i], 0.6 - 0.1 * i) for i in range(6)] + [(bot[0], 0.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_series(image: np.ndarray, series: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    H, W = image.shape[:2]
    # Expected x bin positions (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    bin_xs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"

    if series == "t_lt_2min":
        # Blue 'x' markers
        blue = (B > 100) & ((B - R) > 40) & ((B - G) > 30)
        # Exclude text label "t < 2 min" (roughly c in [930, 1380], r in [760, 880])
        out = []
        for bx in bin_xs:
            c_exp = x_axis.to_pixel(bx)
            c0 = max(0, int(round(c_exp)) - 45)
            c1 = min(W, int(round(c_exp)) + 55)
            sub = np.zeros((H, W), dtype=bool)
            sub[top[1] - 35 : bot[2] + 40, c0:c1] = blue[top[1] - 35 : bot[2] + 40, c0:c1]
            # Mask out "t < 2 min" text region around x = 0.4
            sub[750:895, 920:1390] = False
            blobs = [b for b in find_blobs(sub, min_area=600) if 40 <= b["w"] <= 75 and 40 <= b["h"] <= 75]
            if len(blobs) != 1:
                raise RuntimeError(f"Fig. 2B t<2min at bin {bx}: expected 1 blue x, found {blobs}")
            b = blobs[0]
            out.append(make_point_reading(b["c"], b["r"], x_axis, y_axis, half_px_y=3.0,
                                          note=f"blue x {b['w']}x{b['h']}px, {cal_note}"))
        return out

    elif series == "t_gt_5min":
        # Red filled circles (diameter ~49 px). Note at x=0.8 and 1.0 the black bottom spine bisects the dot,
        # leaving its unobstructed lower crown below the bottom spine.
        red = (R > 180) & (G < 80) & (B < 80)
        # First measure exact radius from the 4 clean interior dots at bins 0.0, 0.2, 0.4, 0.6
        clean_dots = {}
        for bx in bin_xs[:4]:
            c_exp = x_axis.to_pixel(bx)
            c0 = max(0, int(round(c_exp)) - 45)
            c1 = min(W, int(round(c_exp)) + 55)
            sub = np.zeros((H, W), dtype=bool)
            sub[top[1] - 35 : bot[2] + 45, c0:c1] = red[top[1] - 35 : bot[2] + 45, c0:c1]
            # Mask out "t > 5 min" label
            sub[270:410, 480:940] = False
            blobs = [b for b in find_blobs(sub, min_area=1200) if 40 <= b["w"] <= 60 and 40 <= b["h"] <= 60]
            if len(blobs) != 1:
                raise RuntimeError(f"Fig. 2B t>5min at bin {bx}: expected 1 red dot, found {blobs}")
            clean_dots[bx] = blobs[0]

        dot_radius = 0.5 * float(np.median([b["h"] for b in clean_dots.values()]))
        out = []
        for bx in bin_xs:
            if bx in clean_dots:
                b = clean_dots[bx]
                out.append(make_point_reading(b["c"], b["r"], x_axis, y_axis, half_px_y=3.0,
                                              note=f"red dot {b['w']}x{b['h']}px, {cal_note}"))
            else:
                # Reconstruct center from unobstructed lower rim below bottom spine + measured radius
                c_exp = x_axis.to_pixel(bx)
                c0 = max(0, int(round(c_exp)) - 45)
                c1 = min(W, int(round(c_exp)) + 55)
                sub = np.zeros((H, W), dtype=bool)
                sub[bot[2] - 5 : bot[2] + 45, c0:c1] = red[bot[2] - 5 : bot[2] + 45, c0:c1]
                blobs = [b for b in find_blobs(sub, min_area=400) if b["w"] <= 60]
                if len(blobs) != 1:
                    raise RuntimeError(f"Fig. 2B t>5min at bin {bx}: expected 1 lower red arc, found {blobs}")
                b = blobs[0]
                r_center = float(b["r1"]) + 0.5 - dot_radius
                out.append(make_point_reading(b["c"], r_center, x_axis, y_axis, half_px_y=4.0,
                                              note=f"red dot on spine (lower rim r1={b['r1']}, R={dot_radius:.1f}px), {cal_note}"))
        return out

    else:
        raise ValueError(f"Unknown series {series!r}")
