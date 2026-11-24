"""Shared calibrated heatmap colorbar reader for Extended Data Fig. 7 panels a-g."""

from __future__ import annotations

import numpy as np
import _calib as C


def read_ed7_slice(image: np.ndarray, at_nm: list[float], d_um: float, label: str) -> list[dict]:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    colored = sat > 5.0
    dark = image.mean(axis=2) < 128.0

    # Find the heatmap block (left) and narrow colorbar (right)
    H, W = colored.shape
    col_blocks = C.clusters(np.where(colored[int(H * 0.12) : int(H * 0.85), :].mean(axis=0) > 0.03)[0], gap=30)
    # Filter to substantial blocks (width >= 40 px)
    col_blocks = [cb for cb in col_blocks if (cb[2] - cb[1] + 1) >= 40]
    if len(col_blocks) < 2:
        raise RuntimeError(f"Expected heatmap + colorbar in ED Fig. 7 ({label}), found {col_blocks}")

    _, hm_c0, hm_c1 = col_blocks[0]
    _, cb_c0, cb_c1 = col_blocks[-1]

    cb_rows = np.where(colored[:, cb_c0 + 5 : cb_c1 - 5].mean(axis=1) > 0.5)[0]
    cb_r0, cb_r1 = int(cb_rows.min()), int(cb_rows.max())
    # Top of colorbar = 2.00, bottom = 1.00 (verify ticks right of colorbar)
    cb_ticks = C.find_outward_ticks(dark, (float(cb_c1), cb_c0, cb_c1), (cb_r0 - 15, cb_r1 + 15), "y", side=+1, band=(2, 16), frac=0.35)
    if len(cb_ticks) >= 3:
        # Use colorbar colored extent aligned with ticks
        pass
    lut_v, lut_rgb = C.decode_colorbar_rgb(image, (cb_c0 + 8, cb_c1 - 8), float(cb_r0 + 2), float(cb_r1 - 2), 2.00, 1.00)

    hm_rows = np.where(colored[:, hm_c0 : hm_c1 + 1].any(axis=1))[0]
    r0, r1 = int(hm_rows.min()), int(hm_rows.max())
    # 15 equal-width diameter columns: D = 0.50, 0.75, ..., 4.00 um
    col_w = (hm_c1 - hm_c0 + 1) / 15.0
    col_idx = int(round((d_um - 0.5) / 0.25))
    xc = hm_c0 + (col_idx + 0.5) * col_w
    ci = int(round(xc))

    # Calibrate y-axis from left ticks (2.0, 1.5, 1.0, 0.5, 0.005); ticks are ~4-8 px wide immediately left of hm_c0
    yt = C.find_outward_ticks(dark, (float(hm_c0), hm_c0, hm_c0), (r0 - 80, r1 + 30), "y", side=-1, band=(1, 4), frac=0.4, max_width=20)
    yt = C.match_ticks_regular(yt, 5, tol=10.0)
    ya = C.verify_axis("linear", list(zip(yt, [2.0, 1.5, 1.0, 0.5, 0.005])), max_residual_px=6.5)

    col_rows = np.where(colored[:, ci])[0]
    r_top_valid, r_bot_valid = int(col_rows.min()), int(col_rows.max())

    out = []
    for tr_nm in at_nm:
        ry_f = float(ya.to_pixel(tr_nm / 1000.0))
        ri = int(round(ry_f))
        if ri < r_top_valid + 2 or ri > r_bot_valid - 2:
            out.append({"x": float(tr_nm), "y": None, "found": False, "note": f"T_ring={tr_nm} nm outside valid colored bar"})
            continue
        rgb_med = np.median(image[ri - 2 : ri + 3, ci - 12 : ci + 13, :].astype(float), axis=(0, 1))
        val, dist, hw = C.match_colorbar(rgb_med, lut_v, lut_rgb)
        hw_total = max(hw, 0.015)
        lo = max(1.0, val - hw_total)
        hi = min(2.0, val + hw_total)
        out.append(
            {
                "x": float(tr_nm),
                "y": float(val),
                "lo": float(lo),
                "hi": float(hi),
                "px": [float(xc), ry_f],
                "px_lo": float(ry_f - 4.0),
                "px_hi": float(ry_f + 4.0),
                "note": f"EDFig7 {label} D={d_um}um: rgb_dist={dist:.1f}, ya_res={ya.residual_px:.2f}px",
            }
        )
    return out
