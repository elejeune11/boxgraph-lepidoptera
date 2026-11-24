"""Shared calibrated colorbar reader for Extended Data Fig. 6A (wetted fraction) & 6B (final energy)."""

from __future__ import annotations

import math
import numpy as np
import _calib as C

# Discrete 7-point logarithmic sweep grid on both axes (pN/nm)
GRID_VALS = [1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4]


def read_ed6_row(image: np.ndarray, at: list[float], panel: str, row_idx: int = 0) -> list[dict]:
    """Read row `row_idx` (0 = gamma_AVR = 1e-8 pN/nm) across gamma_AVD columns in ED Fig. 6A or 6B."""
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.45, col_frac=0.45)
    ht, hb = hl[0], hl[-1]
    vl0 = vl[0]
    r0, r1 = int(ht[2]) + 1, int(hb[1]) - 1
    c0 = int(vl0[2]) + 1

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    black = (gray < 65.0) & (sat < 25.0)

    # Find the right black spine of the 7x7 heatmap square (box is square, so right spine ~ c0 + box_h)
    box_h = r1 - r0 + 1
    c1_est = c0 + box_h
    v_prof = black[r0 : r1 + 1, c1_est - 80 : c1_est + 80].mean(axis=0)
    vc = C.clusters(np.where(v_prof > 0.75)[0] + c1_est - 80, gap=4)
    if not vc:
        raise RuntimeError(f"Could not locate right black border of ED Fig. 6{panel} heatmap near col {c1_est}")
    c1 = int(vc[0][1]) - 1
    c1_outer = int(vc[0][2])

    # Locate narrow vertical colorbar strictly to the right of the outer black border
    colored = sat > 8.0
    cb_cols = C.clusters(np.where(colored[r0 + 50 : r1 - 50, c1_outer + 10 :].mean(axis=0) > 0.5)[0] + c1_outer + 10, gap=10)
    if not cb_cols:
        raise RuntimeError(f"Could not find colorbar right of col {c1_outer} in ED Fig. 6{panel}")
    _, cb_c0, cb_c1 = cb_cols[0]
    cb_rows = np.where(colored[:, cb_c0 + 5 : cb_c1 - 5].mean(axis=1) > 0.5)[0]
    cb_r0, cb_r1 = int(cb_rows.min()), int(cb_rows.max())

    if panel == "A":
        # Linear colorbar: 0.90 (bottom) to 1.00 (top)
        lut_v, lut_rgb = C.decode_colorbar_rgb(image, (cb_c0 + 10, cb_c1 - 10), float(cb_r0 + 2), float(cb_r1 - 2), 1.00, 0.90, scale="linear")
    elif panel == "B":
        # Log10 colorbar: detect the 4 thin tick marks (10^6, 10^4, 10^2, 10^0) sticking right of cb_c1
        tick_dark = (gray < 175.0) & (sat < 30.0)
        cb_ticks = C.find_outward_ticks(tick_dark, (float(cb_c1), cb_c0, cb_c1), (cb_r0 - 20, cb_r1 + 20), "y", side=+1, band=(2, 16), frac=0.25, max_width=24)
        cb_ticks = C.match_ticks_regular(cb_ticks, 4, tol=10.0)
        cb_ax = C.verify_axis("linear", list(zip(cb_ticks, [6.0, 4.0, 2.0, 0.0])), max_residual_px=6.0)
        val_top = 10.0 ** float(cb_ax.to_value(cb_r0 + 2))
        val_bot = 10.0 ** float(cb_ax.to_value(cb_r1 - 2))
        lut_v, lut_rgb = C.decode_colorbar_rgb(image, (cb_c0 + 10, cb_c1 - 10), float(cb_r0 + 2), float(cb_r1 - 2), val_top, val_bot, scale="log10")
    else:
        raise ValueError(panel)

    cell_h = (r1 - r0 + 1) / 7.0
    cell_w = (c1 - c0 + 1) / 7.0
    ry = r0 + (row_idx + 0.5) * cell_h
    ri = int(round(ry))

    targets = at if at else GRID_VALS[row_idx + 1 :]
    out = []
    for xv in targets:
        # Match xv in log space to nearest grid column
        col_idx = min(range(7), key=lambda j: abs(math.log10(GRID_VALS[j]) - math.log10(xv)))
        xc = c0 + (col_idx + 0.5) * cell_w
        ci = int(round(xc))
        rgb_med = np.median(image[ri - 20 : ri + 21, ci - 20 : ci + 21, :].astype(float), axis=(0, 1))
        if np.all(rgb_med > 250.0):
            out.append({"x": float(xv), "y": None, "found": False, "note": f"white/uncomputed cell ({row_idx},{col_idx})"})
            continue
        if panel == "A":
            val, dist, hw = C.match_colorbar(rgb_med, lut_v, lut_rgb)
            hw_tot = max(hw, 0.002)
            lo = max(0.90, val - hw_tot)
            hi = min(1.00, val + hw_tot)
        else:
            # Log-scale colorbar uncertainty: all LUT values within +5 RGB units, at least +-25%
            dists = np.linalg.norm(lut_rgb - rgb_med[None, :], axis=1)
            idx = int(np.argmin(dists))
            dist = float(dists[idx])
            val = float(lut_v[idx])
            near = lut_v[dists <= dist + 5.0]
            lo = min(float(near.min()), val / 1.25)
            hi = max(float(near.max()), val * 1.25)
        out.append(
            {
                "x": float(xv),
                "y": float(val),
                "lo": float(lo),
                "hi": float(hi),
                "px": [float(xc), float(ry)],
                "px_lo": float(ry - 10.0),
                "px_hi": float(ry + 10.0),
                "note": f"ED6{panel} row={GRID_VALS[row_idx]:g}, col={GRID_VALS[col_idx]:g}: rgb_dist={dist:.1f}",
            }
        )
    return out
