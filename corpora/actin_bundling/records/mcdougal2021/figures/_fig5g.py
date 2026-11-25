"""Shared calibration and cross-detection helpers for Fig. 5G (p. 5)."""
from __future__ import annotations

import numpy as np
from boxgraph.digitize import Axis

from _calib import (
    detect_crosses,
    find_box,
    find_ticks_mask,
    interval_from_px,
    neutral_dark_mask,
)


def calibrate_fig5g(image: np.ndarray):
    """Calibrate Fig. 5G: x = 30..90% dev (7 ticks), y = 1.0..2.0 um (6 right-border ticks, avoiding broken left axis)."""
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    xt = find_ticks_mask(dark, r1, c0 + 15, c1 - 15, "x", -1, offset_lo=12, offset_hi=24)
    if len(xt) != 7:
        raise RuntimeError(f"Fig. 5G expected 7 x-ticks (30..90), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0])))

    # Right border y-ticks cleanly avoid the 'smooth' arrow and the bottom-left axis break
    yt = find_ticks_mask(dark, c1, r0 + 15, r1 - 15, "y", -1, offset_lo=12, offset_hi=24)
    if len(yt) != 6:
        raise RuntimeError(f"Fig. 5G expected 6 right y-ticks (2.0..1.0), found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, [2.0, 1.8, 1.6, 1.4, 1.2, 1.0])))
    return x_axis, y_axis, (c0, r0, c1, r1)


def read_black_smooth(image: np.ndarray) -> list[dict]:
    """Read all 8 black 'x' markers in the smooth early-actin stage (~38-44% dev)."""
    x_axis, y_axis, box = calibrate_fig5g(image)
    c0, r0, c1, r1 = box
    dark = neutral_dark_mask(image, dark_thr=80.0, max_chroma=30.0)
    # Restrict to the lower-left region below the 'smooth' arrow (s < 1.60 um, t < 46%)
    r_cut = int(round(y_axis.to_pixel(1.60)))
    c_cut = int(round(x_axis.to_pixel(46.0)))
    roi = np.zeros_like(dark)
    roi[r_cut:int(r1) - 18, int(c0) + 18:c_cut] = dark[r_cut:int(r1) - 18, int(c0) + 18:c_cut]

    pts = detect_crosses(roi, box, size=63, stroke=11.0, thr=0.65, min_dist=18,
                         margin=(18, 18, 18, 18), neg_radius=25.0, neg_weight=-0.35)
    half_px = 2.5 + max(y_axis.residual_px, 0.5)
    out = []
    for rr, cc, sc in pts:
        xv = float(x_axis.to_value(cc))
        yv, lo, hi = interval_from_px(y_axis, rr, half_px)
        out.append({
            "x": xv,
            "y": yv,
            "lo": lo,
            "hi": hi,
            "px": [cc, rr],
            "px_lo": float(rr - half_px),
            "px_hi": float(rr + half_px),
            "note": f"black cross score={sc:.2f}; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px",
        })
    return out


def read_teal_scale(image: np.ndarray, target_rgb: tuple[int, int, int], tol: float = 14.0,
                    recover_occluded_left: bool = False) -> list[dict]:
    """Read one individual scale's connected trajectory in Fig. 5G isolated by its exact RGB hue."""
    x_axis, y_axis, box = calibrate_fig5g(image)
    rgb = image.astype(float)
    dist = np.sqrt(((rgb - np.array(target_rgb, dtype=float)) ** 2).sum(axis=2))
    mask = dist <= tol

    pts = detect_crosses(mask, box, size=63, stroke=11.0, thr=0.60, min_dist=25,
                         margin=(20, 20, 20, 20), neg_radius=25.0, neg_weight=-0.35)
    half_px = 2.5 + max(y_axis.residual_px, 0.5)
    out = []
    for rr, cc, sc in pts:
        xv = float(x_axis.to_value(cc))
        yv, lo, hi = interval_from_px(y_axis, rr, half_px)
        out.append({
            "x": xv,
            "y": yv,
            "lo": lo,
            "hi": hi,
            "px": [cc, rr],
            "px_lo": float(rr - half_px),
            "px_hi": float(rr + half_px),
            "note": f"scale rgb={target_rgb} score={sc:.2f}; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px",
        })

    # For late scale 7 (0, 76, 76), the leftmost cross at ~61.3% sits directly behind scale 6 (0, 129, 129);
    # extrapolate scale 7's own visible straight line segment in [63.5, 74.0]% back to the common 61.30% x-column
    # of late scales 5 & 6 (found directly from late scale 6's cross at ~61.3%).
    if recover_occluded_left and len(out) == 2:
        d6 = np.sqrt(((rgb - np.array((0, 129, 129), dtype=float)) ** 2).sum(axis=2))
        pts6 = detect_crosses(d6 <= 14.0, box, size=63, stroke=11.0, thr=0.60, min_dist=25,
                              margin=(20, 20, 20, 20), neg_radius=25.0, neg_weight=-0.35)
        if pts6:
            c_left_anchor = pts6[0][1]
            c_a = int(round(x_axis.to_pixel(63.5)))
            c_b = int(round(x_axis.to_pixel(74.0)))
            r_lo = int(round(y_axis.to_pixel(1.88)))
            r_hi = int(round(y_axis.to_pixel(1.70)))
            seg_c, seg_r = [], []
            for col in range(c_a, c_b + 1):
                rows = np.nonzero(mask[r_lo:r_hi, col])[0] + r_lo
                if rows.size > 0:
                    seg_c.append(col + 0.5)
                    seg_r.append(float(rows.mean()) + 0.5)
            if len(seg_c) >= 20:
                sl, ic = np.polyfit(seg_c, seg_r, 1)
                r_est = float(sl * c_left_anchor + ic)
                xv = float(x_axis.to_value(c_left_anchor))
                half_ext = 5.0 + max(y_axis.residual_px, 0.5)
                yv, lo, hi = interval_from_px(y_axis, r_est, half_ext)
                out.insert(0, {
                    "x": xv,
                    "y": yv,
                    "lo": lo,
                    "hi": hi,
                    "px": [c_left_anchor, r_est],
                    "px_lo": float(r_est - half_ext),
                    "px_hi": float(r_est + half_ext),
                    "note": "leftmost vertex occluded behind scale 6; extrapolated from own rgb=(0,76,76) line segment",
                })
    return out
