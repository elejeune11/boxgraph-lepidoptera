"""Calibrated reader for Fig. 5E (SNAP-647-espin-2B molecules/um F-actin bar chart)."""

from __future__ import annotations

import numpy as np

import sys
from pathlib import Path

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=85, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 80), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 7 (unlabeled top), 6, 5, 4, 3, 2, 1 (0 sits at x-spine center)
    yt_raw = C.find_ticks(dark, y_spine, (40, int(x_spine[1]) - 30), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 7, tol=8.0, min_step=160.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])), max_residual_px=2.0)

    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_above = int(x_spine[1]) - 15
    r_below = int(x_spine[2]) + 3

    # Detect bars: bars 1 (gray) and 3 (red) rise above x-spine; bar 2 (blue alpha-Actinin) extends slightly BELOW x-spine!
    band_up = gray[r_above - 60 : r_above, c0 : W - 20] < 235
    runs_up = [r for r in C.runs_1d(band_up.mean(axis=0) > 0.7, gap=3) if (r[1] - r[0] + 1) >= 250]
    R, G, B = image[:, :, 0].astype(np.float32), image[:, :, 1].astype(np.float32), image[:, :, 2].astype(np.float32)
    blue_below = (B[r_below : r_below + 45, c0 : W - 20] > 100) & (B[r_below : r_below + 45, c0 : W - 20] - R[r_below : r_below + 45, c0 : W - 20] > 30)
    runs_dn = [r for r in C.runs_1d(blue_below.any(axis=0), gap=40) if (r[1] - r[0] + 1) >= 250]
    if len(runs_up) != 2 or len(runs_dn) != 1:
        raise RuntimeError(f"Fig. 5E: expected 2 positive + 1 negative bar, got up={runs_up}, dn={runs_dn}")

    bars = sorted([(runs_up[0][0], runs_up[0][1], +1), (runs_dn[0][0], runs_dn[0][1], -1), (runs_up[1][0], runs_up[1][1], +1)], key=lambda t: t[0])
    labels = ["Single filaments (n=15)", "Parallel alpha-Actinin (n=15; slightly negative)", "Parallel Fascin (n=17)"]
    blk = (image[:, :, 0] < 55) & (image[:, :, 1] < 55) & (image[:, :, 2] < 55)

    out = []
    for idx, (a, b, sign) in enumerate(bars, start=1):
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        margin = int(0.12 * (cb - ca))
        if sign > 0:
            left_s = gray[:r_above, ca + margin : int(c_mid) - 35]
            right_s = gray[:r_above, int(c_mid) + 35 : cb - margin]
            is_bar = (left_s.mean(axis=1) < 235) & (right_s.mean(axis=1) < 235)
            top_r = r_above - 1
            while top_r > 0 and is_bar[top_r - 1]:
                top_r -= 1
            edge_row = float(top_r)
        else:
            # Negative bar extending downward from x-spine
            left_s = gray[int(x_spine[2]) + 1 : int(x_spine[2]) + 150, ca + margin : int(c_mid) - 35]
            right_s = gray[int(x_spine[2]) + 1 : int(x_spine[2]) + 150, int(c_mid) + 35 : cb - margin]
            is_bar = (left_s.mean(axis=1) < 235) & (right_s.mean(axis=1) < 235)
            bot_off = 0
            while bot_off + 1 < len(is_bar) and is_bar[bot_off + 1]:
                bot_off += 1
            edge_row = float(int(x_spine[2]) + 1 + bot_off)

        # Trace black vertical error bar around edge_row (masking out x-spine for negative bar)
        eb_mask = blk.copy()
        if sign < 0:
            eb_mask[int(x_spine[1]) - 1 : int(x_spine[2]) + 2, :] = False
        u_row, l_row = C.find_vertical_error_bar(
            eb_mask,
            c_mid,
            edge_row,
            row_min=40,
            row_max=H - 120,
            stem_half_w=12,
            max_gap=20 if sign < 0 else 8,
            cap_half_w=35,
            min_cap_width=18,
            start_up_offset=2,
            start_dn_offset=2,
        )
        out.append(
            C.make_point_reading(
                float(idx),
                c_mid,
                edge_row,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=True,
                err_px=2.0 + y_axis.residual_px,
                note=f"Fig. 5E bar {idx} ({labels[idx - 1]}); cal_res_y={y_axis.residual_px:.2f}px",
            )
        )
    return out
