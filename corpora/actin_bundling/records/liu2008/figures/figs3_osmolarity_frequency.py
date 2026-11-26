"""Calibrated reader for Fig. s3 (p. 10 / SI p. 5): protrusion frequency bar chart with +/- s.e.m. error bars."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    H, W = image.shape[:2]
    blk = C.pure_black_mask(image, max_ch=90.0)
    rgb = image.astype(np.float32)
    red = (rgb[:, :, 0] > 150.0) & (rgb[:, :, 1] < 50.0) & (rgb[:, :, 2] < 50.0)

    # Find the 3 red bars first so we know the plot horizontal span unambiguously
    bar_band_r0, bar_band_r1 = int(0.55 * H), int(0.82 * H)
    bars = C.find_runs(red[bar_band_r0:bar_band_r1, :].mean(axis=0) > 0.5, min_len=200)
    if len(bars) != 3:
        raise RuntimeError(f"Fig. s3: expected 3 red bars, found {bars}")

    # Left y-axis spine is the strong pure-black vertical line immediately left of the first bar
    c_search_lo = max(0, bars[0][0] - 180)
    c_search_hi = bars[0][0] - 10
    v_runs = C.find_runs(blk[int(0.18 * H) : int(0.86 * H), c_search_lo:c_search_hi].mean(axis=0) > 0.70, min_len=4, max_len=40)
    if not v_runs:
        raise RuntimeError(f"Fig. s3: could not locate left y-axis spine in cols [{c_search_lo},{c_search_hi}]")
    c_left_lo = c_search_lo + v_runs[-1][0]

    # Outward y ticks to the left of the y-axis spine: 1.0, 0.8, 0.6, 0.4, 0.2, 0.0
    yt_all = C.scan_ticks(blk, "y", int(0.15 * H), int(0.90 * H), c_left_lo - 38, c_left_lo - 10, min_frac=0.55, min_len=6, max_len=35)
    yt_reg = C.select_regular_ticks(yt_all, 6, tol=6.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=3.0)

    labels = ["Hypoosmotic (N=3)", "Isoosmotic (N=4)", "Hyperosmotic (N=3)"]
    out: list[dict] = []
    r_ceil = int(round(ya.to_pixel(1.08)))
    r_base = int(round(ya.to_pixel(0.02)))

    for idx, (bx0, bx1, cx) in enumerate(bars):
        bw = bx1 - bx0
        c_int = int(round(cx))
        wl = slice(bx0 + int(0.10 * bw), bx0 + int(0.32 * bw))
        wr = slice(bx1 - int(0.32 * bw), bx1 - int(0.10 * bw))

        # Top edge of red fill on both wings away from the centered error bar
        red_rows = np.nonzero((red[r_ceil:r_base, wl].mean(axis=1) > 0.5) & (red[r_ceil:r_base, wr].mean(axis=1) > 0.5))[0]
        if red_rows.size == 0:
            raise RuntimeError(f"Fig. s3 bar {idx + 1}: could not find red fill")
        r_red_top = int(r_ceil + red_rows.min())

        # Walk upward across the bar's black top outline (outer top edge matches symmetric error-bar midpoint)
        # Check mean luminance < 165 on both wings so 1-2 px red/black anti-aliased transition rows are crossed cleanly
        r = r_red_top - 1
        while r > r_ceil and rgb[r, wl].mean() < 165.0 and rgb[r, wr].mean() < 165.0:
            r -= 1
        r_border_top = r + 1
        stroke = max(1, r_red_top - r_border_top)

        # Trace centered pure-black vertical error bar stem up (above r_border_top) and down (into red fill below r_red_top)
        stem = blk[:, c_int - 10 : c_int + 11].any(axis=1)

        ru = r_border_top - 1
        while ru > r_ceil and (stem[ru] or (ru - 1 > r_ceil and stem[ru - 1])):
            ru -= 1
        ru += 1
        cap_u = [rr for rr in range(ru, min(r_border_top, ru + 32)) if blk[rr, c_int - 28 : c_int + 29].sum() >= 26]
        bar_hi_row = 0.5 * (min(cap_u) + max(cap_u)) + 0.5 if cap_u else float(ru)

        rd = r_red_top + 1
        while rd < r_base and (stem[rd] or (rd + 1 < r_base and stem[rd + 1])):
            rd += 1
        rd -= 1
        cap_d = [rr for rr in range(max(r_red_top, rd - 32), rd + 1) if blk[rr, c_int - 28 : c_int + 29].sum() >= 26]
        bar_lo_row = 0.5 * (min(cap_d) + max(cap_d)) + 0.5 if cap_d else float(rd)

        # Outer top edge of bar aligns symmetrically with (bar_hi_row + bar_lo_row)/2;
        # half_px covers half the outline stroke honestly.
        top_row = float(r_border_top) + 0.5
        half_px = max(0.5 * stroke, 3.0)
        note = (
            f"Fig. s3 {labels[idx]}: bar_top={r_border_top} (stroke={stroke}px), "
            f"sem_caps=[{bar_hi_row:.1f},{bar_lo_row:.1f}], ya_res={ya.residual_px:.2f}px"
        )
        out.append(
            C.reading_point(
                x_val=float(idx + 1),
                col=cx,
                row=top_row,
                y_axis=ya,
                half_px=half_px,
                note=note,
                bar_hi_row=bar_hi_row,
                bar_lo_row=bar_lo_row,
                bar_err_px=2.5,
                clamp_min=0.0,
                clamp_max=1.0,
            )
        )
    return out
