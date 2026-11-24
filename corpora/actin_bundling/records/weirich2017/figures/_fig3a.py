"""Shared calibrated reader for Fig. 3A (tactoid coarsening length and aspect ratio vs time, log-log)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
import _calib as C


def calibrate(image: np.ndarray):
    dark200, gray, hl, vl = C.find_box_spines(image, thresh=200.0, row_frac=0.45, col_frac=0.45)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    blk = gray < 80.0

    xt = C.find_ticks(blk, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(4, 25))
    x_vals = [9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    if len(xt) != len(x_vals):
        raise RuntimeError(f"Fig3A: expected {len(x_vals)} x ticks, got {len(xt)}: {xt}")
    xa = C.verify_axis("log10", list(zip(xt, x_vals)), max_residual_px=2.0)

    yt_l = C.find_ticks(blk, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(4, 25))
    yl_vals = [9, 8, 7, 6, 5, 4, 3, 2, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
    if len(yt_l) != len(yl_vals):
        raise RuntimeError(f"Fig3A: expected {len(yl_vals)} left y ticks, got {len(yt_l)}")
    ya_l = C.verify_axis("log10", [(ht[0], 10.0)] + list(zip(yt_l, yl_vals)) + [(hb[0], 0.1)], max_residual_px=2.0)

    yt_r = C.find_ticks(dark200, vr0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=-1, band=(4, 25))
    yr_vals = [9, 8, 7, 6, 5, 4, 3, 2]
    if len(yt_r) != len(yr_vals):
        raise RuntimeError(f"Fig3A: expected {len(yr_vals)} right y ticks, got {len(yt_r)}")
    ya_r = C.verify_axis("log10", [(ht[0], 10.0)] + list(zip(yt_r, yr_vals)) + [(hb[0], 1.0)], max_residual_px=2.0)

    return gray, xa, ya_l, ya_r, int(ht[2]) + 45, int(hb[1]) - 45, int(vl0[2]) + 45, int(vr0[1]) - 45


def _black_circle_centers(gray: np.ndarray, ya_l, r0: int, c0: int, c1: int) -> list[dict]:
    blk = gray < 75.0
    m = np.zeros_like(blk)
    r_cut = int(ya_l.to_pixel(0.6))
    m[r0:r_cut, c0:c1] = blk[r0:r_cut, c0:c1]
    # Distance transform >= 32 cleanly splits the touching circles at t=80 and t=85 min (radius ~45 px)
    cores = ndimage.distance_transform_edt(m) >= 31.0
    blobs = C.find_blobs(cores, min_area=120)
    blobs.sort(key=lambda b: b["c"])
    return blobs


def read_length_data(image: np.ndarray) -> list[dict]:
    gray, xa, ya_l, _, r0, r1, c0, c1 = calibrate(image)
    blobs = _black_circle_centers(gray, ya_l, r0, c0, c1)
    # Upper error bars are dark gray/black (gray < 140), distinct from lower gray aspect-ratio series
    stem = np.zeros_like(gray, dtype=bool)
    r_bot_lim = int(ya_l.to_pixel(0.55))
    stem[r0:r_bot_lim, c0:c1] = gray[r0:r_bot_lim, c0:c1] < 145.0

    out = []
    for b in blobs:
        c_m, r_m = b["c"], b["r"]
        # Circle radius is ~45 px; start tracing error bar 46 px above/below center
        bar_lo, bar_hi, bar_lo_out, bar_hi_out, row_hi, row_lo = C.trace_dashed_error_bar(
            stem, c_m, int(r_m) - 46, int(r_m) + 46, r0, r_bot_lim, ya_l, stem_half_w=5, gap=8
        )
        half_px = 3.0 + ya_l.residual_px
        yv = float(ya_l.to_value(r_m))
        ya_, yb_ = float(ya_l.to_value(r_m - half_px)), float(ya_l.to_value(r_m + half_px))
        out.append(
            {
                "x": float(xa.to_value(c_m)),
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "bar_lo": bar_lo,
                "bar_hi": bar_hi,
                "bar_lo_outer": bar_lo_out,
                "bar_hi_outer": bar_hi_out,
                "px": [c_m, r_m],
                "px_lo": float(row_hi),
                "px_hi": float(row_lo),
                "note": f"Fig3A length circle, ±1 SD=[{bar_lo:.3f},{bar_hi:.3f}], cal_res=({xa.residual_px:.2f},{ya_l.residual_px:.2f})px",
            }
        )
    return out


def read_length_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    gray, xa, ya_l, _, r0, r1, c0, c1 = calibrate(image)
    blobs = _black_circle_centers(gray, ya_l, r0, c0, c1)
    # Black dashed fit line is pure black (gray < 55), darker than the medium-gray error bar stems (~80-120)!
    # Inset well past top inward ticks (r0_fit = r0 + 50, i.e. ~95 px below top spine) and erase circles + stems.
    dash = np.zeros_like(gray, dtype=bool)
    r_top_fit = r0 + 50
    r_cut = int(ya_l.to_pixel(0.6))
    dash[r_top_fit:r_cut, c0:c1] = gray[r_top_fit:r_cut, c0:c1] < 55.0
    H, W = dash.shape
    for b in blobs:
        cm, rm = b["c"] - 0.5, b["r"] - 0.5
        r_lo, r_hi = max(0, int(rm) - 62), min(H, int(rm) + 63)
        c_lo, c_hi = max(0, int(cm) - 62), min(W, int(cm) + 63)
        yy, xx = np.ogrid[r_lo:r_hi, c_lo:c_hi]
        dash[r_lo:r_hi, c_lo:c_hi] &= ((xx - cm) ** 2 + (yy - rm) ** 2) > 58.0**2
        # Erase vertical error-bar stem strip around each dot
        dash[r_top_fit:r_cut, max(0, int(cm) - 24) : min(W, int(cm) + 25)] = False

    return [
        C.sample_dashed_curve_y(
            dash, xa, ya_l, float(x), r_top_fit, r_cut, max_search_dx=350, sigma_dx=160.0, deg=1, note_prefix="Fig3A power-law fit: "
        )
        for x in at
    ]


def read_aspect_ratio_data(image: np.ndarray) -> list[dict]:
    gray, xa, ya_l, ya_r, r0, r1, c0, c1 = calibrate(image)
    # Open squares sit at the same 16 timestamps as the length circles, in the lower half (L/r in [1.3, 3.5])
    blobs_t = _black_circle_centers(gray, ya_l, r0, c0, c1)
    gmask = np.zeros_like(gray, dtype=bool)
    r_top_lim = int(ya_r.to_pixel(3.8))
    r_bot_lim = r1
    gmask[r_top_lim:r_bot_lim, c0:c1] = (gray[r_top_lim:r_bot_lim, c0:c1] > 105.0) & (
        gray[r_top_lim:r_bot_lim, c0:c1] < 210.0
    )

    # Locate square vertical center at each timestamp using matched-filter correlation with open-square frame
    # Square frame is ~80x80 px with ~8 px stroke
    yy, xx = np.ogrid[-45:46, -45:46]
    sq_kern = (
        ((np.abs(np.abs(xx) - 38) <= 5) & (np.abs(yy) <= 42))
        | ((np.abs(np.abs(yy) - 38) <= 5) & (np.abs(xx) <= 42))
    ).astype(np.float32)
    # Subtract interior to avoid matching horizontal line alone
    sq_kern[(np.abs(xx) < 28) & (np.abs(yy) < 28)] = -0.35

    sub_r0 = int(ya_r.to_pixel(2.8))
    sub_r1 = int(ya_r.to_pixel(1.6))
    patch = gmask[sub_r0:sub_r1, :].astype(np.float32)
    score = signal.fftconvolve(patch, sq_kern[::-1, ::-1], mode="same")

    out = []
    for b in blobs_t:
        c_est = int(round(b["c"]))
        c_lo = max(0, c_est - 25)
        c_hi = min(score.shape[1], c_est + 26)
        local = score[:, c_lo:c_hi]
        idx = int(np.argmax(local))
        dr, dc = divmod(idx, local.shape[1])
        r_m = float(sub_r0 + dr) + 0.5
        c_m = float(c_lo + dc) + 0.5

        bar_lo, bar_hi, bar_lo_out, bar_hi_out, row_hi, row_lo = C.trace_dashed_error_bar(
            gmask, c_m, int(r_m) - 44, int(r_m) + 44, r_top_lim, r_bot_lim, ya_r, stem_half_w=6, gap=8
        )
        half_px = 4.0 + ya_r.residual_px
        yv = float(ya_r.to_value(r_m))
        ya_, yb_ = float(ya_r.to_value(r_m - half_px)), float(ya_r.to_value(r_m + half_px))
        out.append(
            {
                "x": float(xa.to_value(c_m)),
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "bar_lo": bar_lo,
                "bar_hi": bar_hi,
                "bar_lo_outer": bar_lo_out,
                "bar_hi_outer": bar_hi_out,
                "px": [c_m, r_m],
                "px_lo": float(row_hi),
                "px_hi": float(row_lo),
                "note": f"Fig3A aspect-ratio square, ±1 SD=[{bar_lo:.3f},{bar_hi:.3f}], cal_res=({xa.residual_px:.2f},{ya_r.residual_px:.2f})px",
            }
        )
    return out


def read_aspect_ratio_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    gray, xa, ya_l, ya_r, r0, r1, c0, c1 = calibrate(image)
    blobs_t = _black_circle_centers(gray, ya_l, r0, c0, c1)
    gmask = np.zeros_like(gray, dtype=bool)
    r_top = int(ya_r.to_pixel(2.35))
    r_bot = int(ya_r.to_pixel(1.85))
    gmask[r_top:r_bot, c0:c1] = (gray[r_top:r_bot, c0:c1] > 105.0) & (gray[r_top:r_bot, c0:c1] < 205.0)
    # Erase columns around all 16 open squares (square half-width ~45 px) so only the clean horizontal gray line remains
    W = gmask.shape[1]
    for b in blobs_t:
        cm = int(round(b["c"]))
        gmask[:, max(0, cm - 50) : min(W, cm + 51)] = False

    return [
        C.sample_dashed_curve_y(
            gmask, xa, ya_r, float(x), r_top, r_bot, max_search_dx=450, sigma_dx=220.0, deg=1, note_prefix="Fig3A AR fit line: "
        )
        for x in at
    ]
