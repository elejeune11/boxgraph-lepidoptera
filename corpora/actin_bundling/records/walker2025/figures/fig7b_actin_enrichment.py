"""Calibrated reader for Fig. 7B: Fold enrichment of G-actin in condensates (broken y-axis 0-30 // 160-190)."""

from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import ndimage
import _common as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(np.float64).mean(axis=2)
    dark = gray < 115.0
    H, W = dark.shape
    h_spines = C.group_spines(dark, min_len=0.5 * W)
    bot = h_spines[-1]
    # Broken vertical y-axis has two segments (upper and lower): detect leftmost spine column
    v_spines = C.group_spines(dark.T, min_len=0.25 * H)
    left = v_spines[0]

    # Collect outward y-ticks along the entire vertical span of the y-axis
    y_strip = dark[:, max(0, left["lo"] - 56) : left["lo"]][:, ::-1].T
    all_yt = C.find_ticks_on_strip(y_strip, 10, bot["hi"] + 6, min_len=5, max_len=50, max_width=28)
    # The broken axis has 4 upper ticks (190, 180, 170, 160), 2 break slashes (//), and 4 lower ticks (30, 20, 10, 0).
    # Split cleanly at the largest gap / middle break around row ~1100-1280
    # Upper ticks are the first 4 evenly spaced ticks; lower ticks are the last 4 evenly spaced ticks!
    upper_cands = [t for t in all_yt if t < 0.45 * H]
    lower_cands = [t for t in all_yt if t > 0.48 * H]
    yt_up = C.match_regular_ticks(upper_cands, 4, tol=5.0)
    yt_lo = C.match_regular_ticks(lower_cands, 4, tol=5.0)
    ya_up = C.make_axis("linear", list(zip(yt_up, [190.0, 180.0, 170.0, 160.0])), max_res_px=4.0)
    ya_lo = C.make_axis("linear", list(zip(yt_lo, [30.0, 20.0, 10.0, 0.0])), max_res_px=4.0)

    # Locate the 4 bars standing on the bottom baseline
    nonwhite = gray < 210.0
    dark_border = gray < 105.0
    c0, c1 = left["hi"] + 20, bot["end"] - 5
    base_slice = nonwhite[bot["lo"] - 25 : bot["lo"] - 2, c0 : c1 + 1]
    segs = C.clusters_1d(np.flatnonzero(base_slice.mean(axis=0) > 0.35) + c0, gap=8)
    bars = [(s[1], s[2]) for s in segs if (s[2] - s[1] + 1) >= 80]
    if len(bars) != 4:
        raise RuntimeError(f"Fig. 7B: expected 4 bars, found {bars}")

    # Break gap row range (between lower axis top and upper axis bottom)
    r_break_top = int(round(yt_up[-1] + 0.5 * (yt_up[1] - yt_up[0])))
    r_break_bot = int(round(yt_lo[0] - 0.1 * (yt_lo[1] - yt_lo[0])))

    out = []
    for k, (bx0, bx1) in enumerate(bars):
        bw = bx1 - bx0 + 1
        cx = 0.5 * (bx0 + bx1) + 0.5
        wl0, wl1 = bx0 + max(8, int(0.08 * bw)), bx0 + max(20, int(0.28 * bw))
        wr0, wr1 = bx1 - max(20, int(0.28 * bw)), bx1 - max(8, int(0.08 * bw))

        # Check if this bar continues into the upper broken segment (bar 0 = mini-Lpd)!
        has_upper = (
            nonwhite[int(round(yt_up[-1])) : r_break_top, wl0 : wl1 + 1].mean() > 0.5
            and nonwhite[int(round(yt_up[-1])) : r_break_top, wr0 : wr1 + 1].mean() > 0.5
        )
        if has_upper:
            ya = ya_up
            r_start = int(round(yt_up[-1]))
            r_ceil = int(round(yt_up[0])) - 10
        else:
            ya = ya_lo
            r_start = bot["lo"] - 1
            r_ceil = r_break_bot

        r = r_start
        while r >= r_ceil:
            if nonwhite[r, wl0 : wl1 + 1].mean() > 0.65 and nonwhite[r, wr0 : wr1 + 1].mean() > 0.65:
                r -= 1
            else:
                break
        r_edge = r + 1
        stroke = 0
        for rr in range(r_edge, min(r_start, r_edge + 20)):
            if dark_border[rr, wl0 : wl1 + 1].mean() > 0.6 and dark_border[rr, wr0 : wr1 + 1].mean() > 0.6:
                stroke += 1
            else:
                break
        top_row = r_edge + 0.5 * max(stroke, 1)
        half_px = max(1.0, 0.5 * stroke)

        # Detect replicate gray circles for verification note
        sat = np.maximum(np.maximum(image[:, :, 0], image[:, :, 1]), image[:, :, 2]).astype(float) - np.minimum(
            np.minimum(image[:, :, 0], image[:, :, 1]), image[:, :, 2]
        ).astype(float)
        gc = (gray > 125.0) & (gray < 200.0) & (sat < 22.0)
        cc0, cc1 = int(round(cx - 0.16 * bw)), int(round(cx + 0.16 * bw))
        sr0 = int(round(yt_up[0])) if has_upper else int(round(ya_lo.to_pixel(35.0)))
        sr1 = r_break_top if has_upper else bot["lo"]
        lab, n_lab = ndimage.label(gc[sr0:sr1, cc0 : cc1 + 1])
        reps = []
        for li in range(1, n_lab + 1):
            rr_i, cc_i = np.nonzero(lab == li)
            if 120 <= len(rr_i) <= 4000:
                reps.append(ya.to_value(rr_i.mean() + sr0 + 0.5))
        reps.sort()
        rep_str = "; reps=[" + ",".join(f"{v:.3g}" for v in reps) + f"] (mean={np.mean(reps):.3g})" if reps else ""

        seg_tag = "upper(160-190)" if has_upper else "lower(0-30)"
        note = f"Fig. 7B cat={k+1} {seg_tag} bar=[{bx0},{bx1}], ya_res={ya.residual_px:.2f}px{rep_str}"
        out.append(
            C.reading_point(
                x_val=float(k + 1),
                col=cx,
                row=top_row,
                y_axis=ya,
                half_px=half_px,
                note=note,
                clamp_min=0.0,
            )
        )
    return out
