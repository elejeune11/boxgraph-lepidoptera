"""Calibrated reader for Extended Data Fig. 4e (20 uM VASP droplet size histogram, n = 1656 droplets)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.35)
    hb, vl0 = hl[-1], vl[0]
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 15), "y", side=-1, band=(6, 30), frac=0.4, max_width=26)
    yt = C.match_ticks_regular(yt, 6, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [250.0, 200.0, 150.0, 100.0, 50.0, 0.0])), max_residual_px=5.0)

    # X-axis ticks sit between bins; 24 ticks from 0.6 to >7.3 (every other tick labelled: 0.6, 1.3, 1.9, ..., 6.9, >7.3)
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(6, 30), frac=0.4, max_width=26)
    # Every 2 ticks = 0.6333 um (bin width ~0.3167 um): 0.63, 0.95, 1.27, ..., 7.28
    # Fit linear diameter axis from the tick indices (tick 0 = left edge 0.6 um; tick 20 = 6.93 um -> step = 0.31667 um)
    if len(xt) < 22:
        raise RuntimeError(f"Expected >= 22 x-ticks in ED Fig. 4e, found {len(xt)}: {xt}")
    xt = sorted(xt)
    # Detect all pink histogram bars (r > 180, r - g > 50, r - b > 40) above hb
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    pink = (r > 180.0) & (r - g > 50.0) & (r - b > 40.0)
    r_top_lim = int(yt[0]) - 10
    r_bot = int(hb[1]) - 2
    c_left = int(vl0[2]) + 5
    pink[:r_top_lim, :] = False
    pink[r_bot:, :] = False
    pink[:, :c_left] = False

    # Find each pink bar column cluster
    bars = C.clusters(np.where(pink.any(axis=0))[0], gap=6)
    bars = [bc for bc in bars if (bc[2] - bc[1] + 1) >= 20]
    if len(bars) < 18:
        raise RuntimeError(f"Expected ~21 pink histogram bars in ED Fig. 4e, found {len(bars)}")

    # Calibrate x-axis from labelled alternating ticks: ticks[0]=0.6, ticks[2]=1.267, ..., ticks[20]=6.933 -> step 0.3167 um
    # Exact printed labels sit under bin centers at ticks[0], ticks[2], ... (0.6, 1.3, 1.9, 2.5, 3.2, 3.8, 4.4, 5.0, 5.7, 6.3, 6.9)
    # Match first 21 ticks (index 0..20) to 0.6333 + i * 0.31667 um (linear regression to printed labels 0.6..6.9)
    bin_step = (6.9 - 0.6) / 20.0  # 0.315 um per bin (matches 0.6, 1.23, 1.86, 2.49, 3.12, 3.75, 4.38, 5.01, 5.64, 6.27, 6.90)
    xa = C.verify_axis("linear", [(xt[0], 0.6), (xt[20], 6.9)], max_residual_px=2.0)

    out = []
    for cx, bx0, bx1 in bars:
        rows = np.where(pink[:, bx0 + 4 : bx1 - 3].mean(axis=1) > 0.6)[0]
        if len(rows) == 0:
            continue
        top_row = float(rows.min())
        x_um = round(float(xa.to_value(cx)), 2)
        out.append(
            C.make_bar_reading(
                x_val=x_um,
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                err_px=2.0 + ya.residual_px,
                note=f"ED4e bin D~{x_um:.2f}um, ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
