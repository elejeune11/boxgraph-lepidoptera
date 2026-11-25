"""Shared reader for Fig. S1 (PDF p. 15 / SI-9): growth exponent alpha vs. domain area fraction phi_A."""

from __future__ import annotations

import numpy as np
from _common import build_axis, connected_centroids, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_figS1(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    # 6 interior x-ticks at 0.15, 0.20, 0.25, 0.30, 0.35, 0.40 plus box edges 0.10 and 0.45
    xt = match_ticks_both_sides(
        scan_ticks(image, box, "bottom_in", offset_min=8, offset_max=35, min_len=12),
        scan_ticks(image, box, "top_in", offset_min=8, offset_max=35, min_len=12),
    )
    x_pos = [c0] + [p for p, _ in xt] + [c1]
    x_vals = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    xa = build_axis("linear", x_pos, x_vals, max_residual_px=2.0)

    # 5 interior y-ticks at 0.7, 0.6, 0.5, 0.4, 0.3 plus box edges 0.8 and 0.2
    yt = match_ticks_both_sides(
        scan_ticks(image, box, "left_in", offset_min=8, offset_max=35, min_len=12),
        scan_ticks(image, box, "right_in", offset_min=8, offset_max=35, min_len=12),
    )
    y_pos = [r0] + [p for p, _ in yt] + [r1]
    y_vals = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
    ya = build_axis("linear", y_pos, y_vals, max_residual_px=2.0)
    return box, xa, ya


def _find_legend_rect(image: np.ndarray, box: tuple[float, float, float, float]) -> tuple[int, int, int, int]:
    """Detect the framed legend box in the lower-right interior of Fig. S1."""
    c0, r0, c1, r1 = [int(round(v)) for v in box]
    # Search lower-right quadrant for neutral dark border lines of legend box
    r_mid = (r0 + r1) // 2
    c_mid = (c0 + c1) // 2
    sub = image[r_mid : r1 - 15, c_mid : c1 - 15].astype(np.float64)
    gray = sub.mean(axis=2)
    chroma = sub.max(axis=2) - sub.min(axis=2)
    dark = (gray < 150) & (chroma < 35)
    # Horizontal lines longer than 300 px (legend top & bottom borders)
    row_sums = dark.sum(axis=1)
    h_rows = np.where(row_sums > 300)[0]
    col_sums = dark.sum(axis=0)
    v_cols = np.where(col_sums > 200)[0]
    if len(h_rows) < 2 or len(v_cols) < 2:
        raise RuntimeError("Could not detect Fig. S1 legend box frame")
    return (
        r_mid + int(h_rows[0]) - 5,
        c_mid + int(v_cols[0]) - 5,
        r_mid + int(h_rows[-1]) + 5,
        c_mid + int(v_cols[-1]) + 5,
    )


def read_series(image: np.ndarray, series: str) -> list[dict]:
    box, xa, ya = calibrate_figS1(image)
    c0, r0, c1, r1 = box
    leg_r0, leg_c0, leg_r1, leg_c1 = _find_legend_rect(image, box)

    rgb = image.astype(np.float64)
    if series == "active":
        # Dark purple open circles (~RGB 122, 1, 119)
        target = np.array([122.0, 1.0, 119.0])
        dist = np.sqrt(((rgb - target) ** 2).sum(axis=2))
        mask = dist < 85.0
    elif series == "passive":
        # Pink open circles (~RGB 247, 104, 161)
        target = np.array([247.0, 104.0, 161.0])
        dist = np.sqrt(((rgb - target) ** 2).sum(axis=2))
        mask = dist < 85.0
    else:
        raise ValueError(f"unknown series {series!r}")

    # Restrict strictly inside plot box and exclude the legend box
    valid = np.zeros_like(mask)
    valid[int(r0) + 15 : int(r1) - 15, int(c0) + 15 : int(c1) - 15] = True
    valid[leg_r0:leg_r1, leg_c0:leg_c1] = False
    mask &= valid

    comps = connected_centroids(mask, min_area=500, max_area=15000)
    comps.sort(key=lambda t: t[0])
    if len(comps) != 5:
        raise RuntimeError(f"Fig. S1 ({series}): expected 5 circular markers, found {len(comps)}")

    out = []
    for col, row, area, (rmin, cmin, rmax, cmax) in comps:
        # Open circle centroid uncertainty: stroke thickness (~10 px) -> centroid accurate to ~3 px;
        # quote conservative half-width of 6 px (~0.0014 in alpha)
        half_px = 6.0
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row,
                half_px,
                note=f"{series} marker centroid (box {rmax-rmin}x{cmax-cmin} px, area {area} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out
