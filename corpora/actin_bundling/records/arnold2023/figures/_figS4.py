"""Shared reader for Fig. S4 (PDF p. 18 / SI-12): near-critical passive bilayer a/a_0 vs. t [s]."""

from __future__ import annotations

import numpy as np

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_figS4(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    # 7 major interior log10 x-ticks (length >= 20 px: 20, 30, 40, 50, 60, 70, 80)
    xt_all = match_ticks_both_sides(
        scan_ticks(image, box, "bottom_in", offset_min=6, offset_max=35, min_len=6),
        scan_ticks(image, box, "top_in", offset_min=6, offset_max=35, min_len=6),
    )
    xt_maj = [(p, l) for p, l in xt_all if l >= 20]
    xa = build_axis("log10", [p for p, _ in xt_maj], [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0], max_residual_px=2.0)

    # 6 major interior log10 y-ticks (2.2, 2.0, 1.8, 1.6, 1.4, 1.2) plus bottom edge 1.0
    yt_l = scan_ticks(image, box, "left_in", offset_min=6, offset_max=35, min_len=20)
    y_pos = [p for p, _ in yt_l] + [r1]
    y_vals = [2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0]
    ya = build_axis("log10", y_pos, y_vals, max_residual_px=2.0)
    return box, xa, ya


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya = calibrate_figS4(image)
    c0, r0, c1, r1 = box
    rgb = image.astype(np.float64)

    if series == "early":
        target = np.array([112.0, 31.0, 117.0])
    elif series == "late":
        target = np.array([241.0, 67.0, 163.0])
    else:
        raise ValueError(f"unknown series {series!r}")

    dist = np.sqrt(((rgb - target) ** 2).sum(axis=2))
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    mask = (dist < 80.0) & (chroma > 40.0)

    # Restrict inside plot box and exclude lower-right legend (t > 52 s and a/a_0 < 1.25)
    valid = np.zeros_like(mask)
    valid[int(r0) + 5 : int(r1) - 5, int(c0) + 5 : int(c1) - 5] = True
    leg_c = int(xa.to_pixel(52.0))
    leg_r = int(ya.to_pixel(1.25))
    valid[leg_r : int(r1), leg_c : int(c1)] = False
    mask &= valid

    out = []
    half_win = 20
    for x in at:
        col = xa.to_pixel(x)
        ci = int(round(col))
        c_lo = max(int(c0) + 5, ci - half_win)
        c_hi = min(int(c1) - 5, ci + half_win + 1)
        strip = mask[int(r0) + 5 : int(r1) - 5, c_lo:c_hi]
        rr, _ = np.nonzero(strip)
        if len(rr) < 10:
            out.append({"x": float(x), "y": None, "found": False, "note": f"no {series} pixels near t={x}"})
            continue
        rows_abs = rr + int(r0) + 5 + 0.5
        row_med = float(np.median(rows_abs))
        p10, p90 = float(np.percentile(rows_abs, 10)), float(np.percentile(rows_abs, 90))
        half_px = max(0.5 * (p90 - p10), 10.0)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_med,
                half_px,
                x_val=float(x),
                note=f"Fig. S4 {series} median of {len(rows_abs)} px (half-span {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out
