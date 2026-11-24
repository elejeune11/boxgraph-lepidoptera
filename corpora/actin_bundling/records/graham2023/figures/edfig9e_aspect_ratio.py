"""Calibrated reader for Extended Data Fig. 9e (droplet aspect ratio vs PEG %, mean + s.d., grey bars over teal lines)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

PEG_CONCS = [3.0, 5.0, 10.0, 15.0]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.35)
    hb, vl0 = hl[-1], vl[0]
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 15), "y", side=-1, band=(6, 30), frac=0.4, max_width=26)
    yt = C.match_ticks_regular(yt, 7, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0])), max_residual_px=5.0)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    # Pure black bar outline and error bars (completely excludes teal background data lines)
    black = (gray < 55.0) & (sat < 30.0)
    # Neutral grey bar fill (RGB ~180,180,180)
    bar_fill = (gray > 140.0) & (gray < 215.0) & (sat < 20.0)

    r_top_lim = int(ya.to_pixel(18.0))
    r_bot = int(hb[1]) - 3
    c_left = int(vl0[2]) + 25
    c_right = image.shape[1] - 20

    # Locate the 4 grey bars near the bottom (y in [0.1, 0.8])
    rb0, rb1 = int(ya.to_pixel(0.8)), int(ya.to_pixel(0.15))
    bar_cols = C.clusters(np.where(bar_fill[rb0:rb1, c_left:c_right].mean(axis=0) > 0.6)[0] + c_left, gap=20)
    bar_cols = [bc for bc in bar_cols if (bc[2] - bc[1] + 1) >= 250]
    if len(bar_cols) != 4:
        raise RuntimeError(f"Expected 4 PEG bars in ED Fig. 9e, found {bar_cols}")

    out = []
    for peg, (cx, bx0, bx1) in zip(PEG_CONCS, bar_cols):
        bw = bx1 - bx0 + 1
        # Bar top border: full-width black horizontal line spanning >= 80% of bar width
        top_row = None
        for r_i in range(r_top_lim, r_bot):
            if black[r_i, bx0 : bx1 + 1].mean() > 0.80:
                r_end = r_i
                while r_end + 1 < r_bot and black[r_end + 1, bx0 : bx1 + 1].mean() > 0.70:
                    r_end += 1
                top_row = 0.5 * (r_i + r_end) + 0.5
                break
        if top_row is None:
            raise RuntimeError(f"Could not find bar top border for {peg}% PEG in ED Fig. 9e")

        # Upper s.d. error cap strictly above top_row (centered, width 35%-75% of bw)
        cap_c0, cap_c1 = int(round(cx - 0.35 * bw)), int(round(cx + 0.35 * bw))
        cap_rows = []
        for r_i in range(r_top_lim, int(round(top_row)) - 8):
            rs_h = C.runs_1d(black[r_i, cap_c0 : cap_c1 + 1], gap=4)
            if rs_h and max(b_ - a_ + 1 for a_, b_ in rs_h) >= int(0.32 * bw):
                cap_rows.append(r_i)
        err_top = C.clusters(np.array(cap_rows, dtype=int), gap=4)[0][0] if cap_rows else None

        rec = C.make_bar_reading(
            x_val=peg,
            col_center=cx,
            top_row=top_row,
            y_axis=ya,
            upper_err_row=err_top,
            lower_err_row=None,
            mirror_lower=True,
            err_px=2.5 + ya.residual_px,
            note=f"ED9e {peg:g}% PEG (mean+s.d.), ya_res={ya.residual_px:.2f}px",
        )
        if "bar_lo" in rec and rec["bar_lo"] < 1.0:
            # Aspect ratio cannot be < 1 physically, keep raw symmetric mean-sd for s.d. fidelity
            pass
        out.append(rec)
    return out
