"""Shared calibration and boundary extraction for Fig. s2 (p. 9 / SI p. 4)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def calibrate_figs2(image: np.ndarray):
    box = C.find_plot_box(image, dark_thr=100.0, sat_thr=35.0, min_h_frac=0.65, min_v_frac=0.55)
    dark = C.dark_neutral_mask(image, dark_thr=100.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Top inward x ticks: 20, 40, 60, 80, 100, 120 (dashed lines at ~69, 100, 130 are also present,
    # so select_regular_ticks cleanly picks the 6 evenly spaced major ticks 20..120; x=140 merges with right spine).
    xt_all = C.scan_ticks(dark, "x", top[1] + 5, top[1] + 35, left[1] + 10, right[0] - 10, min_frac=0.55, min_len=10, max_len=55)
    xt_reg = C.select_regular_ticks(xt_all, 6, tol=8.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [20.0, 40.0, 60.0, 80.0, 100.0, 120.0])), max_res_px=4.0)

    # Right inward y ticks: interior ticks 70, 60, 50, 40, 30, 20 (80 and 10 touch top/bottom spines)
    yt_all = C.scan_ticks(dark, "y", top[1] + 40, bot[0] - 40, right[0] - 55, right[0] - 8, min_frac=0.55, min_len=10, max_len=55)
    yt_reg = C.select_regular_ticks(yt_all, 6, tol=8.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [70.0, 60.0, 50.0, 40.0, 30.0, 20.0])), max_res_px=3.0)

    return box, xa, ya


def read_boundary(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    box, xa, ya = calibrate_figs2(image)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]
    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    r0, r1 = top[1] + 4, bot[0] - 2

    if curve == "sigma_002":
        # Red boundary stroke has RGB ~ (205, 5, 12) with B < 32, whereas crimson shaded fill has B > 40.
        stroke_m = (R > 165.0) & (G < 55.0) & (B < 35.0)
        fill_m = (R > 165.0) & ((R - G) > 90.0) & ((R - B) > 80.0)
        lbl = "Fig. s2 red boundary (sigma=0.02)"
    elif curve == "sigma_0005":
        # Dark-blue boundary stroke has RGB ~ (20, 41, 130); light-blue shaded fill (>100 nm) has R > 120, G > 115.
        # Below 100 nm, dark-blue fill extends down to the bottom spine, so we track the stroke top edge +
        # calibrated perpendicular stroke half-thickness when the dark-blue run reaches the bottom spine.
        stroke_m = (B > 95.0) & (R < 75.0) & (G < 85.0) & ((B - R) > 40.0)
        fill_m = stroke_m
        lbl = "Fig. s2 blue boundary (sigma=0.005)"
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    # Estimate typical stroke thickness on the unobstructed region (blue: x in [108, 118]; red: x in [45, 58])
    ref_xs = [108.0, 112.0, 115.0, 118.0] if curve == "sigma_0005" else [45.0, 50.0, 55.0]
    ref_thicks = []
    for rx in ref_xs:
        cc = int(round(xa.to_pixel(rx)))
        runs = C.find_runs(stroke_m[r0:r1, cc - 2 : cc + 3].mean(axis=1) >= 0.4, min_len=8)
        if runs:
            best_r = runs[0]
            if r0 + best_r[1] < r1 - 15:  # clearly separated from bottom spine
                ref_thicks.append(best_r[1] - best_r[0])
    # Default vertical stroke width at moderate slope ~ 50 px -> perpendicular half-width ~ 22 px
    base_half_w = 0.5 * float(np.median(ref_thicks)) * 0.78 if ref_thicks else 22.0

    # First pass: locate raw top edge & stroke center at each requested x
    raw = []
    for xv in at:
        cf = xa.to_pixel(xv)
        c0 = int(round(cf))
        # Avoid dashed black vertical lines (L_prot^F,sigma_H ~69, L_prot^B = 100, L_prot^F,sigma_L ~130)
        hit = None
        for off in [0, -18, 18, -30, 30, -45, 45]:
            cc = c0 + off
            if cc < left[1] + 5 or cc >= right[0] - 5:
                continue
            col_s = stroke_m[r0:r1, cc - 2 : cc + 3].mean(axis=1) >= 0.4
            runs = C.find_runs(col_s, min_len=6)
            if not runs and curve == "sigma_002":
                col_f = fill_m[r0:r1, cc - 2 : cc + 3].mean(axis=1) >= 0.4
                runs = C.find_runs(col_f, min_len=6)
            if runs:
                # Topmost substantial run is the boundary stroke
                a, b, ctr = runs[0]
                hit = (r0 + a, r0 + b, r0 + ctr, float(cc) + 0.5)
                break
        raw.append((xv, cf, hit))

    out = []
    for i, (xv, cf, hit) in enumerate(raw):
        if hit is None:
            out.append({"x": float(xv), "y": None, "note": f"{lbl}: no pixels found at x={xv:g}"})
            continue
        ra, rb, rctr, used_c = hit
        # Estimate local slope from neighboring top edges
        slope_px = 0.0
        for nb in (i - 1, i + 1):
            if 0 <= nb < len(raw) and raw[nb][2] is not None:
                dc = abs(raw[nb][1] - cf)
                if dc > 1.0:
                    slope_px = max(slope_px, abs(raw[nb][2][0] - ra) / dc)

        touches_bot = rb >= r1 - 12
        if curve == "sigma_0005" and touches_bot:
            # Dark-blue fill merges with stroke below L_prot <= 100 nm: use top edge + calibrated stroke half-width
            # (clamped inside the visible dark-blue band)
            est_half = min(base_half_w * math.hypot(1.0, slope_px), 0.5 * (rb - ra))
            row_ctr = ra + est_half
            half_thick = max(est_half, 12.0)
            mode_note = "top-edge+half-stroke"
        else:
            row_ctr = rctr
            half_thick = 0.5 * (rb - ra)
            mode_note = "stroke-center"

        col_err = abs(used_c - cf) + max(xa.residual_px, 0.5)
        half_px = math.hypot(half_thick, slope_px * col_err)
        note = (
            f"{lbl} ({mode_note}, width={rb-ra}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
            + (f", col_shift={used_c - cf:+.1f}px" if abs(used_c - cf) > 1.5 else "")
            + ")"
        )
        out.append(
            C.reading_point(
                x_val=xv,
                col=cf,
                row=row_ctr,
                y_axis=ya,
                half_px=half_px,
                note=note,
                clamp_min=0.0,
            )
        )
    return out
