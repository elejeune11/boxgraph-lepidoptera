"""Shared reader for Fig. S6 (PDF p. 19 / SI-13): phase-field mass <phi>(t) vs. t [s]."""

from __future__ import annotations

import numpy as np

from _common import build_axis, find_box, reading_from_px, scan_ticks


def read_panel(image: np.ndarray, at: list[float], x_tick_vals: list[float], label: str) -> list[dict]:
    box = find_box(image)
    c0, r0, c1, r1 = box

    xt = scan_ticks(image, box, "bottom_out", offset_min=8, offset_max=38, min_len=15)
    yt = scan_ticks(image, box, "left_out", offset_min=8, offset_max=38, min_len=15)
    xa = build_axis("linear", [p for p, _ in xt], x_tick_vals, max_residual_px=2.0)
    ya = build_axis("linear", [p for p, _ in yt], [0.40, 0.38, 0.36, 0.34, 0.32, 0.30], max_residual_px=2.0)

    rgb = image.astype(np.float64)
    dark = (rgb.mean(axis=2) < 90.0) & ((rgb.max(axis=2) - rgb.min(axis=2)) < 35.0)
    r_top, r_bot = int(r0) + 25, int(r1) - 25

    out = []
    for x in at:
        col = xa.to_pixel(x)
        ci = int(round(col))
        strip = dark[r_top:r_bot, max(int(c0) + 10, ci - 3) : min(int(c1) - 10, ci + 4)]
        rr, _ = np.nonzero(strip)
        if len(rr) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": f"no curve pixels at t={x}"})
            continue
        rows_abs = rr + r_top + 0.5
        row_mid = float(rows_abs.mean())
        half_px = max(0.5 * float(rows_abs.max() - rows_abs.min()), 2.0)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_mid,
                half_px,
                x_val=float(x),
                note=f"Fig. S6({label}) line half-width {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px",
            )
        )
    return out
