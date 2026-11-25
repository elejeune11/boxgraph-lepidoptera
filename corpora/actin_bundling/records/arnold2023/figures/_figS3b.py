"""Shared reader for Fig. S3(B) (PDF p. 17 / SI-11): a/a_0 vs. t [s] for 'Fast' and 'Slow' regions."""

from __future__ import annotations

import numpy as np

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_figS3b(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    # 11 interior log10 x-ticks (8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90) plus right box edge 100
    xt = match_ticks_both_sides(
        scan_ticks(image, box, "bottom_in", offset_min=6, offset_max=35, min_len=6),
        scan_ticks(image, box, "top_in", offset_min=6, offset_max=35, min_len=6),
    )
    x_pos = [p for p, _ in xt] + [c1]
    x_vals = [8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    xa = build_axis("log10", x_pos, x_vals, max_residual_px=2.0)

    # 5 major interior log10 y-ticks (length >= 20 px: 2.0, 1.8, 1.6, 1.4, 1.2) plus bottom edge 1.0
    yt_all = match_ticks_both_sides(
        scan_ticks(image, box, "left_in", offset_min=6, offset_max=35, min_len=6),
        scan_ticks(image, box, "right_in", offset_min=6, offset_max=35, min_len=6),
    )
    yt_maj = [(p, l) for p, l in yt_all if l >= 20]
    y_pos = [p for p, _ in yt_maj] + [r1]
    y_vals = [2.0, 1.8, 1.6, 1.4, 1.2, 1.0]
    ya = build_axis("log10", y_pos, y_vals, max_residual_px=2.0)
    return box, xa, ya


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya = calibrate_figS3b(image)
    c0, r0, c1, r1 = box
    rgb = image.astype(np.float64)

    if series == "fast":
        target = np.array([122.0, 1.0, 119.0])
    elif series == "slow":
        target = np.array([247.0, 104.0, 161.0])
    else:
        raise ValueError(f"unknown series {series!r}")

    dist = np.sqrt(((rgb - target) ** 2).sum(axis=2))
    mask = dist < 85.0

    # Restrict inside plot box and exclude lower-right legend (t > 32 s and a/a_0 < 1.32)
    valid = np.zeros_like(mask)
    valid[int(r0) + 5 : int(r1) - 5, int(c0) + 5 : int(c1) - 5] = True
    leg_c = int(xa.to_pixel(32.0))
    leg_r = int(ya.to_pixel(1.32))
    valid[leg_r : int(r1), leg_c : int(c1)] = False
    mask &= valid

    out = []
    half_win = 22
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
        half_px = max(0.5 * (p90 - p10), 12.0)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_med,
                half_px,
                x_val=float(x),
                note=f"Fig. S3(B) {series} median of {len(rows_abs)} px (half-span {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out
