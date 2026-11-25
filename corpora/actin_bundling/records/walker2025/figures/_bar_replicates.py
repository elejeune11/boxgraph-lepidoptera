"""Shared calibrated reader for bar charts with black bar outlines and overlaid gray replicate circles
(Figs. 1H, 3F, 4D, 5J, 6K, 7D, 7G, S1C)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _common as C


def find_bars_by_vertical_walls(
    img: np.ndarray,
    bot: dict,
    left: dict,
    y_axis,
    y_top_val: float,
    n_bars: int,
    c_end: int | None = None,
    dark_thresh: float = 110.0,
) -> list[tuple[int, int]]:
    """Detect the n_bars vertical bar (bx0, bx1) column spans standing on the bottom spine."""
    rgb = img.astype(np.float64)
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    # Any non-white pixel (bar fill or black border) right above the baseline
    nonwhite = (gray < 210.0) | (sat > 30.0)

    c0 = left["hi"] + 15
    c1 = (c_end if c_end is not None else bot["end"]) - 5
    r_base = bot["lo"] - 1
    r_ceil = max(0, int(round(y_axis.to_pixel(y_top_val))))

    # Check a narrow band right above the baseline (rows r_base - 18 .. r_base - 3)
    # Even a very short bar (e.g., 0 uM actin at y ~ 0.008 -> ~20 px at 1200 dpi) sits right on the baseline!
    band_h = min(14, max(4, (r_base - r_ceil) // 40))
    base_slice = nonwhite[max(r_ceil, r_base - band_h) : r_base - 1, c0 : c1 + 1]
    col_on = base_slice.mean(axis=0) > 0.35
    segs = C.clusters_1d(np.flatnonzero(col_on) + c0, gap=8)
    # Filter to bar-like widths (>= 60 px at 1200 dpi)
    bars = [(s[1], s[2]) for s in segs if (s[2] - s[1] + 1) >= 60]
    if len(bars) != n_bars:
        # Fallback: scan taller band up to 40% height for bars with zero-height bars missing at baseline
        tall = nonwhite[int(round(y_axis.to_pixel(0.5 * y_top_val))) : r_base - 2, c0 : c1 + 1]
        col_on2 = tall.sum(axis=0) > 10
        segs2 = C.clusters_1d(np.flatnonzero(col_on2) + c0, gap=8)
        bars2 = [(s[1], s[2]) for s in segs2 if (s[2] - s[1] + 1) >= 60]
        if len(bars2) == n_bars:
            bars = bars2
    return bars


def read_bar_tops_and_replicates(
    img: np.ndarray,
    y_tick_vals: list[float],
    x_vals: list[float],
    n_bars: int | None = None,
    select_indices: list[int] | None = None,
    y_max_search: float | None = None,
    dark_thresh: float = 115.0,
    max_tick_len: int = 70,
    max_res_px: float = 5.5,
    y_ticks_bottom_n: int | None = None,
    label: str = "",
) -> list[dict]:
    """Read bar tops (and check overlaid replicate gray circles) for a standard Walker 2025 bar chart."""
    n_total = n_bars if n_bars is not None else len(x_vals)
    ax = C.find_l_axes(img, dark_thresh=dark_thresh, min_h_frac=0.42, min_v_frac=0.30, max_tick_len=max_tick_len)
    ycands = ax["y_ticks"][-y_ticks_bottom_n:] if y_ticks_bottom_n is not None else ax["y_ticks"]
    yt = C.match_regular_ticks(ycands, len(y_tick_vals), tol=9.0)
    ya = C.make_axis("linear", list(zip(yt, y_tick_vals)), max_res_px=max_res_px)

    bot = ax["bot"]
    left = ax["left"]
    y_top = y_max_search if y_max_search is not None else max(y_tick_vals)
    bars = find_bars_by_vertical_walls(img, bot, left, ya, y_top, n_total, dark_thresh=dark_thresh)
    if len(bars) != n_total:
        raise RuntimeError(f"{label}: expected {n_total} bars, found {len(bars)}: {bars}")

    rgb = img.astype(np.float64)
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    # Dark bar outline or bar fill (non-white)
    dark_border = (gray < 110.0) & (sat < 30.0)
    nonwhite = (gray < 215.0) | (sat > 30.0)

    r_zero = float(ya.to_pixel(0.0))
    r_ceil = max(0, int(round(ya.to_pixel(y_top))))

    indices = select_indices if select_indices is not None else list(range(n_total))
    out = []
    for xv, bi in zip(x_vals, indices):
        bx0, bx1 = bars[bi]
        bw = bx1 - bx0 + 1
        cx = 0.5 * (bx0 + bx1) + 0.5

        # Outer wings of bar (10%-30% from left and right walls) avoid centered replicate circles AND centered error bars
        wl0, wl1 = bx0 + max(8, int(0.08 * bw)), bx0 + max(20, int(0.28 * bw))
        wr0, wr1 = bx1 - max(20, int(0.28 * bw)), bx1 - max(8, int(0.08 * bw))

        # Walk upward from baseline along both wings while bar fill/border is present!
        # This is 100% immune to significance brackets sitting above the bar!
        r_start = min(bot["lo"] - 1, int(round(r_zero)) - 1)
        top_row = r_zero
        # Find highest row connected continuously to the baseline in BOTH wings
        r = r_start
        while r >= r_ceil:
            l_on = nonwhite[r, wl0 : wl1 + 1].mean() > 0.65
            r_on = nonwhite[r, wr0 : wr1 + 1].mean() > 0.65
            if l_on and r_on:
                r -= 1
            else:
                break
        # r + 1 is the topmost dark/filled pixel of the bar top border.
        # Measure thickness of the black top border stroke downward from r + 1 (typically 4-8 px)
        r_edge_outer = r + 1
        if r_edge_outer >= r_start:
            # Extremely short / near-zero bar: check if even a thin horizontal stroke exists right above baseline
            top_row = r_zero
            half_px = 1.5
        else:
            stroke_rows = 0
            for rr in range(r_edge_outer, min(r_start, r_edge_outer + 15)):
                if dark_border[rr, wl0 : wl1 + 1].mean() > 0.6 and dark_border[rr, wr0 : wr1 + 1].mean() > 0.6:
                    stroke_rows += 1
                else:
                    break
            # Center of the top horizontal border stroke (or outer edge + 0.5*stroke)
            if stroke_rows >= 2:
                top_row = r_edge_outer + 0.5 * stroke_rows
                half_px = max(1.0, 0.5 * stroke_rows)
            else:
                top_row = float(r_edge_outer) + 0.5
                half_px = 1.5

        # Also detect overlaid gray replicate circles near center column cx (for note verification)
        # Replicate circles have gray interior (~120-185, low sat) ringed by black border.
        rep_vals = []
        cc0 = max(bx0, int(round(cx - 0.16 * bw)))
        cc1 = min(bx1, int(round(cx + 0.16 * bw)))
        # Only search below any significance brackets (stop at first full-width black horizontal line above bar)
        r_rep_top = r_ceil
        for rr in range(max(r_ceil, int(round(top_row)) - 30), r_ceil, -1):
            # Significance bracket spans much wider than the bar
            if dark_border[rr, max(0, bx0 - 40) : min(img.shape[1], bx1 + 40)].mean() > 0.7:
                r_rep_top = rr + 15
                break
        gray_circle = (gray > 115.0) & (gray < 195.0) & (sat < 22.0)
        sub_gc = gray_circle[r_rep_top : bot["lo"] + 5, cc0 : cc1 + 1]
        lab, n_lab = ndimage.label(sub_gc)
        for li in range(1, n_lab + 1):
            rr_i, cc_i = np.nonzero(lab == li)
            if 120 <= len(rr_i) <= 4500:
                h_i = rr_i.max() - rr_i.min() + 1
                w_i = cc_i.max() - cc_i.min() + 1
                if 0.65 <= h_i / max(w_i, 1) <= 1.55 and abs((cc_i.mean() + cc0) - cx) <= 0.12 * bw:
                    rep_vals.append(ya.to_value(rr_i.mean() + r_rep_top + 0.5))
        rep_vals.sort()
        rep_str = ""
        if rep_vals:
            rep_str = "; reps=[" + ",".join(f"{v:.3g}" for v in rep_vals) + f"] (mean={np.mean(rep_vals):.3g})"

        note = f"{label} bar={bi+1} [{bx0},{bx1}], ya_res={ya.residual_px:.2f}px{rep_str}"
        out.append(
            C.reading_point(
                x_val=float(xv),
                col=cx,
                row=top_row,
                y_axis=ya,
                half_px=half_px,
                note=note,
                clamp_min=0.0,
            )
        )
    return out
