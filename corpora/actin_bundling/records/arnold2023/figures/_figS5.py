"""Shared reader for Fig. S5 (PDF p. 18 / SI-12): passive Cahn-Hilliard a/a_0 vs. t (Area/perimeter & Structure factor)."""

from __future__ import annotations

import numpy as np

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_figS5(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    # 16 interior log10 x-ticks (30..90, 100, 200..900) plus right box edge 1000
    xt = match_ticks_both_sides(
        scan_ticks(image, box, "bottom_in", offset_min=8, offset_max=35, min_len=6),
        scan_ticks(image, box, "top_in", offset_min=8, offset_max=35, min_len=6),
    )
    x_pos = [p for p, _ in xt] + [c1]
    x_vals = [30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    xa = build_axis("log10", x_pos, x_vals, max_residual_px=2.0)

    # 8 log10 y-ticks: 10, 9, 8, 7, 6, 5, 4, 3 (left_in has 10..4, right_in has 9,7,6,5,4,3; merge unique)
    yt_l = scan_ticks(image, box, "left_in", offset_min=8, offset_max=35, min_len=6)
    yt_r = scan_ticks(image, box, "right_in", offset_min=8, offset_max=35, min_len=6)
    merged = []
    for p, l in sorted(yt_l + yt_r, key=lambda t: t[0]):
        if not merged or abs(p - merged[-1][0]) > 15.0:
            merged.append([p])
        else:
            merged[-1].append(p)
    y_pos = [float(np.mean(g)) for g in merged]
    y_vals = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0]
    ya = build_axis("log10", y_pos, y_vals, max_residual_px=2.0)
    return box, xa, ya


def _red_mask(image: np.ndarray) -> np.ndarray:
    rgb = image.astype(np.float64)
    return (rgb[:, :, 0] > 180.0) & (rgb[:, :, 1] < 80.0) & (rgb[:, :, 2] < 80.0)


def _black_mask(image: np.ndarray) -> np.ndarray:
    rgb = image.astype(np.float64)
    gray = rgb.mean(axis=2)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    return (gray < 90.0) & (chroma < 35.0)


def _sample_curve(
    mask: np.ndarray,
    box: tuple[float, float, float, float],
    xa,
    ya,
    at: list[float],
    series: str,
    red_rows_by_col: dict[int, float] | None = None,
) -> list[dict]:
    c0, r0, c1, r1 = box
    out = []
    # Restrict below legend (r > r0 + 520) and strictly inside box
    r_top = int(r0) + 250
    r_bot = int(r1) - 35
    half_win = 10

    for x in at:
        col = xa.to_pixel(x)
        ci = int(round(col))
        c_lo = max(int(c0) + 20, ci - half_win)
        c_hi = min(int(c1) - 10, ci + half_win + 1)

        # For black structure-factor curve, search strictly below the red area/perimeter curve
        # so the black text labels ("alpha = 0.30", "alpha = 0.34") above the curves are never hit.
        if series == "sf" and red_rows_by_col is not None and ci in red_rows_by_col:
            r_search_top = max(r_top, int(red_rows_by_col[ci]) + 80)
        else:
            r_search_top = r_top

        strip = mask[r_search_top:r_bot, c_lo:c_hi]
        rr, cc = np.nonzero(strip)
        if len(rr) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": f"no {series} pixels near t={x}"})
            continue

        # In case of small stray text fragments above the black curve at large t, keep the lowest dense cluster
        rows_abs = rr + r_search_top + 0.5
        if series == "sf":
            # Split vertical gaps > 20 px and take the bottom-most cluster (the curve itself, excluding text above)
            u_rows = np.sort(np.unique(rr))
            gaps = np.where(np.diff(u_rows) > 20)[0]
            if len(gaps):
                cut = u_rows[gaps[-1] + 1]
                keep = rr >= cut
                rows_abs = rows_abs[keep]

        row_med = float(np.median(rows_abs))
        # Half-span of the plotted dots/line ribbon (95% range across the vertical cross-section)
        r_p05, r_p95 = float(np.percentile(rows_abs, 5)), float(np.percentile(rows_abs, 95))
        half_px = max(0.5 * (r_p95 - r_p05), 4.0)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_med,
                half_px,
                x_val=float(x),
                note=f"{series} median of {len(rows_abs)} px (ribbon half-width {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya = calibrate_figS5(image)
    red = _red_mask(image)
    c0, r0, c1, r1 = box
    # Exclude legend in upper-left (r < r0 + 560, c < c0 + 1650)
    red[int(r0) : int(r0) + 560, int(c0) : int(c0) + 1650] = False

    if series == "ap":
        return _sample_curve(red, box, xa, ya, at, "ap")
    elif series == "sf":
        # First locate red curve center row at each requested x to bound the black search strictly below it
        red_rows: dict[int, float] = {}
        for x in at:
            ci = int(round(xa.to_pixel(x)))
            strip = red[int(r0) + 250 : int(r1) - 35, max(int(c0) + 20, ci - 15) : min(int(c1) - 10, ci + 16)]
            rr, _ = np.nonzero(strip)
            if len(rr):
                red_rows[ci] = float(np.median(rr)) + int(r0) + 250
        blk = _black_mask(image)
        blk[int(r0) : int(r0) + 560, int(c0) : int(c0) + 1650] = False
        return _sample_curve(blk, box, xa, ya, at, "sf", red_rows_by_col=red_rows)
    else:
        raise ValueError(f"unknown series {series!r}")
