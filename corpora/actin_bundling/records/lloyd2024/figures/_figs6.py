"""Shared readers for SI Fig. 6C regression + 95% CI and SI Fig. 6D 20 individual-by-region-by-type points (PDF p. 29, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D
from _common import run_centers, y_interval
from _fig2 import read_lm_and_ci


def read_figs6c_reg(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    c_cands = np.where((gray[600:2900, 200:700] < 140.0).mean(axis=0) > 0.85)[0] + 200
    c0 = float(c_cands.mean() + 0.5)
    r_cands = np.where((gray[2800:3150, 500:3100] < 140.0).mean(axis=1) > 0.85)[0] + 2800
    r1 = float(r_cands.mean() + 0.5)

    yt = run_centers(np.where(gray[500:int(r1) - 10, int(c_cands.min()) - 6] < 180.0)[0] + 500, gap=4)
    xt = run_centers(np.where(gray[int(round(r1)) + 8, int(round(c0)) + 10:3180] < 160.0)[0] + int(round(c0)) + 10, gap=4)
    if len(xt) != 4 or len(yt) != 5:
        raise RuntimeError(f"SI Fig. 6C ticks mismatch: xt={xt}, yt={yt}")
    xa = D.Axis("linear", list(zip(xt, [30.0, 40.0, 50.0, 60.0])))
    ya = D.Axis("linear", list(zip(yt, [70.0, 60.0, 50.0, 40.0, 30.0])))

    fit_rng = (26.0, 62.5) if series == "blue" else (25.5, 60.2)
    return read_lm_and_ci(image, xa, ya, c0, r1, r_top=480, c_right=3160,
                          series=series, at=at, x_fit_range=fit_rng)


def read_figs6d_markers(image: np.ndarray, region: str, scale_type: str) -> list[dict]:
    """Read the 5 isolated scatter points for a given (region: proximal|distal, scale_type: cover|ground) in SI Fig. 6D."""
    gray = image.astype(float).mean(axis=2)
    c0 = 300.0 + float(np.argmin(gray[600:2800, 300:600].mean(axis=0)))
    r1 = 2850.0 + float(np.argmin(gray[2850:3150, 600:3000].mean(axis=1)))
    yt = run_centers(np.where(gray[500:int(r1) - 10, int(round(c0)) - 8] < 160.0)[0] + 500, gap=4)
    xt = run_centers(np.where(gray[int(round(r1)) + 8, int(round(c0)) + 10:3180] < 160.0)[0] + int(round(c0)) + 10, gap=4)
    if len(xt) != 5 or len(yt) != 5:
        raise RuntimeError(f"SI Fig. 6D ticks mismatch: xt={xt}, yt={yt}")
    xa = D.Axis("linear", list(zip(xt, [0.85, 0.90, 0.95, 1.00, 1.05])))
    ya = D.Axis("linear", list(zip(yt, [70.0, 60.0, 50.0, 40.0, 30.0])))

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    # Strictly inside plot area (row >= 500 excludes top marginal density curve)
    in_box = np.zeros(image.shape[:2], dtype=bool)
    in_box[500:int(r1) - 15, int(c0) + 15:3180] = True
    if region == "proximal":
        mask = in_box & (B - R > 35.0) & (B > 170.0)
    else:
        mask = in_box & (gray < 190.0) & (np.abs(B - R) < 20.0)

    lab, n = ndimage.label(mask)
    out = []
    for k in range(1, n + 1):
        rr, cc = np.where(lab == k)
        if len(rr) < 1500:
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        fill = len(rr) / float(h * w)
        # Equilateral triangles have fill ~ 0.51-0.57; circles have fill ~ 0.76-0.80
        is_ground = fill < 0.65
        if (scale_type == "ground") != is_ground:
            continue
        # For an upward triangle, geometric center of ggplot2 shape=17 is ~0.58 of height from top (or centroid at 2/3;
        # ggplot2 places triangle center at its circumcenter / bounding-box ~0.62 from top, which closely equals centroid ~0.65).
        rc = float(rr.mean() + 0.5)
        cc_m = float(cc.mean() + 0.5)
        xv = float(xa.to_value(cc_m))
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, rc, 4.0)
        out.append({
            "x": xv, "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [cc_m, rc], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"H.e.dem {region} {scale_type} (fill={fill:.2f}, rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px)",
        })
    out.sort(key=lambda r: r["x"])
    return out
