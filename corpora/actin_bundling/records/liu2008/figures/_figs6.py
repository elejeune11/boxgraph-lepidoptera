"""Shared calibration, marker detection, and curve tracing for Fig. s6 (p. 14 / SI p. 9)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import ndimage
import _common as C


def calibrate_figs6(image: np.ndarray):
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.55, min_v_frac=0.65)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Bottom inward x ticks: 0, 500, 1000, 1500, 2000 (the 4 unobstructed ticks 500..2000 are crisp;
    # t=0 tick intersects y=0 tick near bottom-left corner, so use regular ticks including 500..2000)
    xt_bot = C.scan_ticks(dark, "x", bot[0] - 42, bot[0] - 12, left[1] + 5, right[0] - 5, min_frac=0.5, min_len=6, max_len=25)
    xt_reg = C.select_regular_ticks(xt_bot, 5, tol=6.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [0.0, 500.0, 1000.0, 1500.0, 2000.0])), max_res_px=3.0)

    # Left inward y ticks: 12, 10, 8, 6, 4, 2, 0
    yt_l = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 10, left[1] + 45, min_frac=0.5, min_len=6, max_len=30)
    yt_reg = C.select_regular_ticks(yt_l, 7, tol=8.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0])), max_res_px=3.5)

    return box, xa, ya


def read_markers(image: np.ndarray, kind: str) -> list[dict]:
    box, xa, ya = calibrate_figs6(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    blk = C.pure_black_mask(image, max_ch=90.0)

    # Interior excluding top/bottom/left ticks, while keeping the rightmost square at t ~ 2010 s
    # (right y ticks have height <= 20 px and area < 750 px, whereas markers have height >= 36 px and area >= 1000 px)
    mask = np.zeros_like(blk)
    mask[top[1] + 50 : bot[0] - 50, left[1] + 50 : right[0] - 6] = blk[top[1] + 50 : bot[0] - 50, left[1] + 50 : right[0] - 6]
    lab, n = ndimage.label(mask)

    squares = []
    crosses = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        area = len(rr)
        if area < 950:
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        if h < 32 or w < 30 or h > 75 or w > 75:
            continue
        fill = area / float(h * w)
        c0 = float(cc.mean()) + 0.5
        r0 = float(rr.mean()) + 0.5
        if area >= 1650 and fill >= 0.74:
            squares.append((c0, r0, area, fill))
        elif area < 1650 and fill < 0.75:
            crosses.append((c0, r0, area, fill))

    pts = sorted(squares if kind == "squares" else crosses, key=lambda t: t[0])
    expected = 11 if kind == "squares" else 6
    if len(pts) != expected:
        raise RuntimeError(f"Fig. s6 {kind}: expected {expected} markers, found {len(pts)}")

    lbl = "Fig. s6 protrusion 1 (black square)" if kind == "squares" else "Fig. s6 protrusion 2 (black cross)"
    out = []
    for c0, r0, area, fill in pts:
        note = f"{lbl}: area={area}px, fill={fill:.2f}, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
        out.append(C.reading_point(x_val=xa.to_value(c0), col=c0, row=r0, y_axis=ya, half_px=3.0, note=note, clamp_min=0.0))
    return out


def read_curve(image: np.ndarray, at: list[float], color: str) -> list[dict]:
    box, xa, ya = calibrate_figs6(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    interior = np.zeros(image.shape[:2], dtype=bool)
    interior[top[1] + 6 : bot[0] - 6, left[1] + 6 : right[0] - 6] = True

    if color == "red":
        mask = interior & (R > 140.0) & ((R - G) > 65.0) & ((R - B) > 65.0)
        lbl = "Fig. s6 kinetic model prediction (red)"
    elif color == "blue":
        mask = interior & (B > 100.0) & ((B - R) > 45.0) & ((B - G) > 25.0)
        lbl = "Fig. s6 kinetic model prediction (blue)"
    else:
        raise ValueError(f"Unknown color {color!r}")

    # Collect all unobstructed column centers along the curve so any column occluded by a black square/cross
    # is interpolated cleanly from immediately flanking colored pixels on both sides.
    r0, r1 = top[1] + 6, bot[0] - 6
    c_cols = []
    c_rows = []
    c_half = []
    for c in range(left[1] + 8, right[0] - 8):
        runs = C.find_runs(mask[r0:r1, c], min_len=8)
        if runs:
            best = max(runs, key=lambda t: t[1] - t[0])
            c_cols.append(float(c) + 0.5)
            c_rows.append(r0 + best[2])
            c_half.append(0.5 * (best[1] - best[0]))
    cols_arr = np.array(c_cols)
    rows_arr = np.array(c_rows)
    half_arr = np.array(c_half)

    out = []
    for xv in at:
        cf = xa.to_pixel(xv)
        idx = int(np.searchsorted(cols_arr, cf))
        if idx == 0 or idx >= len(cols_arr):
            if abs(cols_arr[min(idx, len(cols_arr) - 1)] - cf) <= 10.0:
                k = min(idx, len(cols_arr) - 1)
                row_ctr = float(rows_arr[k])
                half_thick = float(half_arr[k])
                gap = abs(cols_arr[k] - cf)
                slope_px = 0.0
            else:
                out.append({"x": float(xv), "y": None, "note": f"{lbl}: x={xv:g} outside curve range"})
                continue
        else:
            cL, cR = cols_arr[idx - 1], cols_arr[idx]
            rL, rR = rows_arr[idx - 1], rows_arr[idx]
            frac = (cf - cL) / max(cR - cL, 1e-9)
            row_ctr = float(rL + frac * (rR - rL))
            half_thick = float(max(half_arr[idx - 1], half_arr[idx]))
            gap = min(cf - cL, cR - cf)
            slope_px = abs(rR - rL) / max(cR - cL, 1.0)
        # Slope-propagated x-calibration error + interpolation allowance across markers
        col_err = max(xa.residual_px, 0.5) + 0.25 * gap
        half_px = math.hypot(half_thick, slope_px * col_err)
        note = (
            f"{lbl}: width={2*half_thick:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
            + (f", interp_gap={gap:.1f}px" if gap > 2.0 else "")
        )
        out.append(C.reading_point(x_val=xv, col=cf, row=row_ctr, y_axis=ya, half_px=half_px, note=note, clamp_min=0.0))
    return out
