"""Calibrated reader for Supplemental Fig. S5b ([Arp2/3]=0 nM): median filament length L(600s) vs [VCA].

Reads the horizontal coordinate (Filament Length in um) of each of the 5 vertical CDF=0.5 steps
corresponding to [VCA] = 0, 1, 2, 5, 10 uM, directly testing SI Eqs. (2)-(3) VCA monomer competition.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

VCA_PALETTE = [
    (0.0, np.array([66.0, 151.0, 82.0])),
    (1.0, np.array([108.0, 165.0, 106.0])),
    (2.0, np.array([198.0, 161.0, 187.0])),
    (5.0, np.array([165.0, 90.0, 140.0])),
    (10.0, np.array([102.0, 34.0, 85.0])),
]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=100.0, row_frac=0.45, col_frac=0.45)
    # Leftmost and rightmost vertical borders of the plot frame
    vl_l, vl_r = vl[0], vl[-1]
    hb_t, hb_b = hl[0], hl[-1]

    # Calibrate y-axis (CDF: 1.00, 0.75, 0.50, 0.25, 0.00)
    yt_c = C.find_outward_ticks(
        dark, vl_l, (int(hb_t[1]) - 15, int(hb_b[2]) + 15),
        "y", side=-1, band=(4, 28), frac=0.45, max_width=22,
    )
    yt = C.match_ticks_regular(yt_c, 5, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.75, 0.50, 0.25, 0.0])), max_residual_px=4.0)

    # Calibrate x-axis (Filament Length: 0, 2, 4, 6 um)
    xt_c = C.find_outward_ticks(
        dark, hb_b, (int(vl_l[1]) - 15, int(vl_r[2]) + 15),
        "x", side=+1, band=(4, 28), frac=0.45, max_width=22,
    )
    xt = C.match_ticks_regular(xt_c, 4, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 2.0, 4.0, 6.0])), max_residual_px=4.0)

    # Read each VCA step's column position across CDF in [0.35, 0.65] (centered at CDF = 0.50)
    r_mid_lo = int(round(ya.to_pixel(0.65)))
    r_mid_hi = int(round(ya.to_pixel(0.35)))
    r_half = float(ya.to_pixel(0.50))
    c_left = int(vl_l[2]) + 5
    c_right = int(vl_r[1]) - 5

    sub = image[r_mid_lo:r_mid_hi, c_left:c_right].astype(np.float32)
    out = []
    for vca_uM, target in VCA_PALETTE:
        d2 = (sub[:, :, 0] - target[0]) ** 2 + (sub[:, :, 1] - target[1]) ** 2 + (sub[:, :, 2] - target[2]) ** 2
        mask = d2 < 18.0 ** 2
        # Keep only the main solid stroke column cluster (rejects 1-2 px antialiased fringes on darker curves)
        col_counts = mask.sum(axis=0)
        col_runs = C.runs_1d(col_counts >= max(10, int(0.25 * mask.shape[0])), gap=4)
        if not col_runs:
            raise RuntimeError(f"Fig. S5b (0 nM Arp2/3): VCA={vca_uM} uM curve not found")
        c0_s, c1_s = max(col_runs, key=lambda t: int(col_counts[t[0] : t[1] + 1].sum()))
        sub_mask = np.zeros_like(mask)
        sub_mask[:, c0_s : c1_s + 1] = mask[:, c0_s : c1_s + 1]
        rr, cc = np.nonzero(sub_mask)
        if len(cc) < 10:
            raise RuntimeError(f"Fig. S5b (0 nM Arp2/3): VCA={vca_uM} uM curve not found")
        cols = cc.astype(float) + c_left + 0.5
        rows = rr.astype(float) + r_mid_lo + 0.5
        # Linear fit col(row) across [0.35, 0.65] to evaluate exact crossing at CDF = 0.50
        slope, intercept = np.polyfit(rows, cols, 1)
        c_est = float(slope * r_half + intercept)
        # Half-width covers line stroke thickness + calibration residual + slight CDF slope width
        resid = np.abs(cols - (slope * rows + intercept))
        half_px = float(np.percentile(resid, 95)) + 1.5 + xa.residual_px
        l_val = float(xa.to_value(c_est))
        la, lb = float(xa.to_value(c_est - half_px)), float(xa.to_value(c_est + half_px))
        out.append({
            "x": vca_uM,
            "y": l_val,
            "lo": min(la, lb),
            "hi": max(la, lb),
            "px": [c_est, r_half],
            "note": f"FigS5b [Arp2/3]=0nM [VCA]={vca_uM:g}uM CDF=0.50 step, xa_res={xa.residual_px:.2f}px, n_px={len(cc)}",
        })
    return out
