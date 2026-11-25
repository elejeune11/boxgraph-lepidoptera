"""Shared calibrated reader for Fig. 1C (droplet diameter normalized frequency histograms, 4 conditions)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

SERIES_RGB = {
    "vasp": np.array([176.0, 16.0, 192.0]),
    "arp23": np.array([224.0, 160.0, 224.0]),
    "vca": np.array([152.0, 200.0, 136.0]),
    "vca_arp23": np.array([0.0, 144.0, 64.0]),
}


def read_series(image: np.ndarray, series: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=120.0, row_frac=0.45, col_frac=0.45)
    vl0 = vl[0]
    hb = hl[-1]

    # Outward y-ticks along left spine: 6 ticks at rows corresponding to 0.18, 0.14, 0.10, (0.06 misprinted 0.08), (0.02 misprinted 0.04), 0.00
    yt_all = C.find_outward_ticks(dark, vl0, (200, int(hb[2]) + 20), "y", side=-1, band=(5, 35), frac=0.45, max_width=28)
    if len(yt_all) != 6:
        raise RuntimeError(f"Expected 6 y-ticks in Fig. 1C, found {yt_all}")
    # Calibrate on the 4 mutually consistent ticks (0.18, 0.14, 0.10, 0.00); ticks 3 & 4 are half/1.5-step (0.06, 0.02) misprinted as 0.08, 0.04
    ya = C.verify_axis(
        "linear",
        [(yt_all[0], 0.18), (yt_all[1], 0.14), (yt_all[2], 0.10), (yt_all[5], 0.0)],
        max_residual_px=4.0,
    )

    # Outward x-ticks along bottom spine: 22 ticks from 0.0 to 8.4 um in steps of 0.4 um
    xt_all = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 20, image.shape[1] - 10), "x", side=+1, band=(5, 35), frac=0.45, max_width=28)
    xt = C.match_ticks_regular(xt_all, 22, tol=8.0)
    # Match the 5 labelled major ticks 0, 2, 4, 6, 8 um (indices 0, 5, 10, 15, 20)
    xa = C.verify_axis(
        "linear",
        [(xt[0], 0.0), (xt[5], 2.0), (xt[10], 4.0), (xt[15], 6.0), (xt[20], 8.0)],
        max_residual_px=6.0,
    )

    # Isolate bars of the requested series color inside plot interior
    target = SERIES_RGB[series]
    rgb = image.astype(float)
    r_top_lim = int(yt_all[0]) - 60
    r_base_edge = int(hb[1])  # top pixel row of x-axis spine
    r_zero = float(yt_all[5])  # center of x-axis spine (y = 0)
    c_left = int(vl0[2]) + 5
    c_right = int(xt[-1]) + 60

    sub = rgb[r_top_lim:r_base_edge, c_left:c_right]
    mask = np.linalg.norm(sub - target[None, None, :], axis=2) < 42.0

    # Find each bar in each 0.4 um bin [0.4*k, 0.4*(k+1)] for k = 0..20
    out = []
    for k in range(21):
        x_lo = 0.4 * k
        x_hi = 0.4 * (k + 1)
        x_center = round(0.4 * k + 0.2, 4)
        c0 = max(0, int(round(xa.to_pixel(x_lo))) - c_left)
        c1 = min(mask.shape[1], int(round(xa.to_pixel(x_hi))) - c_left)
        if c1 <= c0:
            continue
        bin_mask = mask[:, c0:c1]
        # Require a solid bar standing on the baseline (bottom 40 px has pixels)
        col_bot = bin_mask[-40:, :].any(axis=0)
        cols = np.where(col_bot)[0]
        if len(cols) < 15:
            continue
        # Use largest contiguous run of columns in this bin
        runs = C.runs_1d(col_bot, gap=3)
        seg = max(runs, key=lambda t: t[1] - t[0])
        if seg[1] - seg[0] + 1 < 15:
            continue
        # Inset by 3 px from bar side edges to avoid antialiasing
        s0, s1 = seg[0] + 3, seg[1] - 3
        if s1 <= s0:
            s0, s1 = seg[0], seg[1]
        bar_slice = bin_mask[:, s0 : s1 + 1]
        row_runs = C.runs_1d(bar_slice.mean(axis=1) > 0.5, gap=2)
        if not row_runs:
            continue
        # Pick the contiguous bar segment touching the baseline (ignoring floating legend swatches above)
        bot_seg = row_runs[-1]
        if (bin_mask.shape[0] - 1) - bot_seg[1] > 15:
            continue
        top_row = float(bot_seg[0] + r_top_lim)
        col_center = float(0.5 * (seg[0] + seg[1]) + 0.5 + c0 + c_left)
        # Skip tiny sub-pixel baseline slivers (< 3 px high)
        if r_zero - top_row < 3.0:
            continue
        out.append(
            C.make_bar_reading(
                x_val=x_center,
                col_center=col_center,
                top_row=top_row,
                y_axis=ya,
                err_px=2.0 + ya.residual_px,
                note=f"Fig1C {series}, bin=[{x_lo:.1f},{x_hi:.1f}]um, ya_res={ya.residual_px:.2f}px",
            )
        )
    if not out:
        raise RuntimeError(f"No bars detected for Fig. 1C series {series!r}")
    return out
