"""Shared calibrated colorbar heatmap reader for Main Fig. 4g (gamma_I = 1.0 and 5.0 pN)."""

from __future__ import annotations

import numpy as np
import _calib as C


def _parse_fig4g(image: np.ndarray, block_idx: int):
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    colored = sat > 12.0
    dark = image.mean(axis=2) < 128.0

    col_blocks = C.clusters(np.where(colored[400:1600, :].mean(axis=0) > 0.03)[0], gap=40)
    if len(col_blocks) != 3:
        raise RuntimeError(f"Expected 2 heatmaps + 1 colorbar in Fig. 4g, found {len(col_blocks)} blocks")

    # Rightmost narrow block is the vertical colorbar (1.00 to 2.00)
    _, cb_c0, cb_c1 = col_blocks[2]
    cb_rows = np.where(colored[:, cb_c0 + 5 : cb_c1 - 5].mean(axis=1) > 0.5)[0]
    cb_r0, cb_r1 = int(cb_rows.min()), int(cb_rows.max())
    cb_ticks = C.find_outward_ticks(dark, (float(cb_c1), cb_c0, cb_c1), (cb_r0 - 20, cb_r1 + 20), "y", side=+1, band=(5, 35), frac=0.4)
    cb_ticks = C.match_ticks_regular(cb_ticks, 5, tol=4.0)
    lut_v, lut_rgb = C.decode_colorbar_rgb(image, (cb_c0 + 10, cb_c1 - 10), cb_ticks[0], cb_ticks[-1], 2.00, 1.00)

    _, c0, c1 = col_blocks[block_idx]
    hm_rows = np.where(colored[:, c0 : c1 + 1].any(axis=1))[0]
    r0, r1 = int(hm_rows.min()), int(hm_rows.max())
    yt = C.find_outward_ticks(dark, (float(c0), c0, c0), (r0 - 100, r1 + 30), "y", side=-1, band=(15, 75), frac=0.4)
    xt = C.find_outward_ticks(dark, (float(r1), r1, r1), (c0 - 10, c1 + 10), "x", side=+1, band=(15, 65), frac=0.4)
    yt = C.match_ticks_regular(yt, 5, tol=4.0)
    xt = C.match_ticks_regular(xt, 8, tol=4.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])), max_residual_px=2.5)
    ya = C.verify_axis("linear", list(zip(yt, [2.0, 1.5, 1.0, 0.5, 0.005])), max_residual_px=3.0)
    return xa, ya, lut_v, lut_rgb, colored


def read_slice(image: np.ndarray, at_nm: list[float], block_idx: int, d_um: float, label: str) -> list[dict]:
    xa, ya, lut_v, lut_rgb, colored = _parse_fig4g(image, block_idx)
    xc = float(xa.to_pixel(d_um))
    ci = int(round(xc))
    col_rows = np.where(colored[:, ci])[0]
    r_top_valid, r_bot_valid = int(col_rows.min()), int(col_rows.max())

    out = []
    for tr_nm in at_nm:
        ry_f = float(ya.to_pixel(tr_nm / 1000.0))
        ri = int(round(ry_f))
        if ri < r_top_valid + 2 or ri > r_bot_valid - 2:
            out.append({"x": float(tr_nm), "y": None, "found": False, "note": f"T_ring={tr_nm} nm outside valid colored bar [{r_top_valid},{r_bot_valid}]"})
            continue
        rgb_med = np.median(image[ri - 2 : ri + 3, ci - 15 : ci + 16, :].astype(float), axis=(0, 1))
        val, dist, hw = C.match_colorbar(rgb_med, lut_v, lut_rgb)
        if dist > 20.0:
            out.append({"x": float(tr_nm), "y": None, "found": False, "note": f"poor colorbar match dist={dist:.1f}"})
            continue
        hw_total = max(hw, 0.012)
        lo = max(1.0, val - hw_total)
        hi = min(2.0, val + hw_total)
        out.append(
            {
                "x": float(tr_nm),
                "y": float(val),
                "lo": float(lo),
                "hi": float(hi),
                "px": [xc, ry_f],
                "px_lo": float(ry_f - 4.0),
                "px_hi": float(ry_f + 4.0),
                "note": f"Fig4g {label} D={d_um}um: rgb_dist={dist:.1f}, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    return out
