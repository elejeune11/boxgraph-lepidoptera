"""Calibrated reader for Fig. 6D: WT Eps15 partitioning bar chart with +/- SD error bars (n = 3)."""

from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C
from _bar_replicates import find_bars_by_vertical_walls


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    ax = C.find_l_axes(image, dark_thresh=115.0, min_h_frac=0.45, min_v_frac=0.45, max_tick_len=55)
    y_vals = [float(2 * k) for k in range(10, -1, -1)]
    yt = C.match_regular_ticks(ax["y_ticks"], len(y_vals), tol=4.5)
    ya = C.make_axis("linear", list(zip(yt, y_vals)), max_res_px=4.0)

    bars = find_bars_by_vertical_walls(image, ax["bot"], ax["left"], ya, 20.0, 2)
    if len(bars) != 2:
        raise RuntimeError(f"Fig. 6D: expected 2 bars, found {bars}")

    rgb = image.astype(np.float64)
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    dark_ink = (gray < 55.0) & (sat < 20.0)
    nonwhite = gray < 210.0

    r_zero = float(ya.to_pixel(0.0))
    r_ceil = int(round(ya.to_pixel(17.0)))  # below the n.s. bracket at ~16.5? Let's check dynamically!
    out = []
    for k, (bx0, bx1) in enumerate(bars):
        bw = bx1 - bx0 + 1
        cx = 0.5 * (bx0 + bx1) + 0.5
        wl0, wl1 = bx0 + max(8, int(0.08 * bw)), bx0 + max(20, int(0.30 * bw))
        wr0, wr1 = bx1 - max(20, int(0.30 * bw)), bx1 - max(8, int(0.08 * bw))

        # Walk up from baseline on outer wings to find exact bar top
        r = ax["bot"]["lo"] - 1
        while r >= int(round(ya.to_pixel(20.0))):
            if nonwhite[r, wl0 : wl1 + 1].mean() > 0.65 and nonwhite[r, wr0 : wr1 + 1].mean() > 0.65:
                r -= 1
            else:
                break
        r_edge = r + 1
        stroke = 0
        for rr in range(r_edge, min(ax["bot"]["lo"], r_edge + 20)):
            if dark_ink[rr, wl0 : wl1 + 1].mean() > 0.6 and dark_ink[rr, wr0 : wr1 + 1].mean() > 0.6:
                stroke += 1
            else:
                break
        top_row = r_edge + 0.5 * max(stroke, 1)
        half_px = max(1.0, 0.5 * stroke)

        # Trace centered pure-black vertical SD error bar stem UP above r_edge and DOWN below r_edge+stroke
        c_int = int(round(cx))
        stem = dark_ink[:, c_int - 4 : c_int + 5].any(axis=1)

        # Walk UP from r_edge - 1
        up_last = r_edge
        miss = 0
        for rr in range(r_edge - 1, max(0, r_edge - 300), -1):
            if stem[rr]:
                up_last = rr
                miss = 0
            else:
                miss += 1
                if miss > 3:
                    break
        # Refine upper cap center (horizontal bar cap wider than stem)
        cap_rows_up = [rr for rr in range(up_last, min(r_edge, up_last + 18))
                       if dark_ink[rr, c_int - 18 : c_int + 19].sum() >= 18]
        bar_hi_row = 0.5 * (min(cap_rows_up) + max(cap_rows_up)) + 0.5 if cap_rows_up else float(up_last) + 2.0

        # Walk DOWN from r_edge + stroke + 1 (pure black stem stands out cleanly against medium grey bar fill ~85!)
        dn_start = r_edge + max(stroke, 4) + 1
        dn_last = dn_start
        miss = 0
        for rr in range(dn_start, min(ax["bot"]["lo"] - 5, dn_start + 300)):
            if stem[rr]:
                dn_last = rr
                miss = 0
            else:
                miss += 1
                if miss > 3:
                    break
        cap_rows_dn = [rr for rr in range(max(dn_start, dn_last - 18), dn_last + 1)
                       if dark_ink[rr, c_int - 18 : c_int + 19].sum() >= 18]
        bar_lo_row = 0.5 * (min(cap_rows_dn) + max(cap_rows_dn)) + 0.5 if cap_rows_dn else float(dn_last) - 2.0

        note = f"Fig. 6D cat={k+1} bar=[{bx0},{bx1}], ya_res={ya.residual_px:.2f}px"
        out.append(
            C.reading_point(
                x_val=float(k + 1),
                col=cx,
                row=top_row,
                y_axis=ya,
                half_px=half_px,
                note=note,
                bar_hi_row=bar_hi_row,
                bar_lo_row=bar_lo_row,
                bar_err_px=2.0,
                clamp_min=0.0,
            )
        )
    return out
