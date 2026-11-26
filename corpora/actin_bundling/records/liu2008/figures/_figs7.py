"""Shared calibration, circle detection, and dotted-fit tracing for Fig. s7 (p. 15 / SI p. 10)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import ndimage
import _common as C


def calibrate_figs7(image: np.ndarray):
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.55, min_v_frac=0.65)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # X ticks: 0, 50, 100, 150, 200, 250, 300, 350 s (t=0 clean on top; t=350 clean on bottom)
    xt_top = C.scan_ticks(dark, "x", top[1] + 10, top[1] + 40, left[1] + 5, right[0] - 5, min_frac=0.5, min_len=6, max_len=25)
    xt_bot = C.scan_ticks(dark, "x", bot[0] - 40, bot[0] - 10, left[1] + 5, right[0] - 5, min_frac=0.5, min_len=6, max_len=25)
    xt_comb = sorted(xt_top + [x for x in xt_bot if all(abs(x - t) > 20.0 for t in xt_top)])
    xt_reg = C.select_regular_ticks(xt_comb, 8, tol=6.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0])), max_res_px=2.5)

    # Left inward y ticks: 0.52, 0.51, 0.50, 0.49, 0.48
    yt_l = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 10, left[1] + 40, min_frac=0.5, min_len=6, max_len=25)
    yt_reg = C.select_regular_ticks(yt_l, 5, tol=6.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [0.52, 0.51, 0.50, 0.49, 0.48])), max_res_px=3.0)

    return box, xa, ya


def read_exp(image: np.ndarray) -> list[dict]:
    box, xa, ya = calibrate_figs7(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    blk = C.pure_black_mask(image, max_ch=90.0)

    # Interior excluding ticks along all 4 spines, EXCEPT keeping circles near bottom-left (t~0) and top-right (t~367)
    mask = np.zeros_like(blk)
    mask[top[1] + 6 : bot[0] - 6, left[1] + 6 : right[0] - 6] = blk[top[1] + 6 : bot[0] - 6, left[1] + 6 : right[0] - 6]

    # Erase tick rectangles precisely around calibrated tick positions (width 22 px, length 52 px)
    for xv in [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0]:
        cx = int(round(xa.to_pixel(xv)))
        mask[top[1] : top[1] + 56, max(0, cx - 12) : cx + 13] = False
        if xv > 10.0:  # preserve the circle near (0 s, 0.477) which sits just above the bottom spine at x=0
            mask[bot[0] - 58 : bot[0], max(0, cx - 12) : cx + 13] = False
        else:
            mask[bot[0] - 12 : bot[0], max(0, cx - 12) : cx + 13] = False
    for yv in [0.52, 0.51, 0.50, 0.49, 0.48]:
        ry = int(round(ya.to_pixel(yv)))
        mask[max(0, ry - 14) : ry + 15, left[1] : left[1] + 52] = False
        mask[max(0, ry - 14) : ry + 15, right[0] - 52 : right[0]] = False

    # Close across the red dotted fit line (dots ~16 px wide) so split halves of each black open circle fuse
    yy, xx = np.ogrid[-14:15, -14:15]
    closed = ndimage.binary_closing(mask, structure=(xx * xx + yy * yy) <= 14 * 14)
    lab, n = ndimage.label(closed)
    pts = []
    for i in range(1, n + 1):
        # Measure bounding box using original unclosed black pixels belonging to this component
        comp_orig = mask & (lab == i)
        rr, cc = np.nonzero(comp_orig)
        if len(rr) < 400:
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        if h < 35 or w < 35 or h > 95 or w > 95:
            continue
        c0 = 0.5 * (cc.min() + cc.max()) + 0.5
        r0 = 0.5 * (rr.min() + rr.max()) + 0.5
        pts.append((c0, r0, len(rr), w, h))

    pts.sort(key=lambda t: t[0])
    if len(pts) != 5:
        raise RuntimeError(f"Fig. s7 exp: expected 5 black open circles, found {len(pts)}: {pts}")

    out = []
    for c0, r0, area, w, h in pts:
        note = f"Fig. s7 FRAP circle (box={w}x{h}px, area={area}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px)"
        out.append(C.reading_point(x_val=xa.to_value(c0), col=c0, row=r0, y_axis=ya, half_px=3.0, note=note))
    return out


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    box, xa, ya = calibrate_figs7(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    r0, r1 = top[1] + 6, bot[0] - 6
    red = (R > 140.0) & ((R - G) > 75.0) & ((R - B) > 75.0)

    # Column-by-column red stroke center (handles both separated dots at t < 180 s and merged dots at t >= 180 s)
    c_cols = []
    c_rows = []
    c_half = []
    for c in range(left[1] + 6, right[0] - 6):
        runs = C.find_runs(red[r0:r1, c], min_len=4)
        if runs:
            best = max(runs, key=lambda t: t[1] - t[0])
            c_cols.append(float(c) + 0.5)
            c_rows.append(r0 + best[2])
            c_half.append(0.5 * (best[1] - best[0]))
    if len(c_cols) < 50:
        raise RuntimeError(f"Fig. s7 fit: expected red fit columns, found only {len(c_cols)}")

    cols = np.array(c_cols)
    rows = np.array(c_rows)
    halfs = np.array(c_half)

    out = []
    for xv in at:
        cf = xa.to_pixel(xv)
        idx = int(np.searchsorted(cols, cf))
        if idx == 0 or idx >= len(cols):
            k = 0 if idx == 0 else len(cols) - 1
            if abs(cols[k] - cf) > 25.0:
                out.append({"x": float(xv), "y": None, "note": f"Fig. s7 fit: x={xv:g} outside dotted curve"})
                continue
            row_ctr = float(rows[k])
            half_thick = float(halfs[k])
            gap = abs(cols[k] - cf)
            slope_px = 0.0
        else:
            cL, cR = cols[idx - 1], cols[idx]
            rL, rR = rows[idx - 1], rows[idx]
            frac = (cf - cL) / max(cR - cL, 1e-9)
            row_ctr = float(rL + frac * (rR - rL))
            half_thick = float(max(halfs[idx - 1], halfs[idx]))
            gap = min(cf - cL, cR - cf)
            slope_px = abs(rR - rL) / max(cR - cL, 1.0)
        col_err = max(xa.residual_px, 0.5) + 0.25 * gap
        half_px = math.hypot(half_thick, slope_px * col_err)
        note = f"Fig. s7 dotted red fit: width={2*half_thick:.1f}px, gap={gap:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
        out.append(C.reading_point(x_val=xv, col=cf, row=row_ctr, y_axis=ya, half_px=half_px, note=note))
    return out
