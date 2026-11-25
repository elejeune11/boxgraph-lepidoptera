"""Shared calibration and cross/mean-curve readers for Fig. 5C-F (p. 5)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

from _calib import (
    detect_crosses,
    find_box,
    find_ticks_mask,
    interval_from_px,
    neutral_dark_mask,
)


def calibrate_fig5cf(image: np.ndarray, y_tick_vals: list[float]):
    """Calibrate any of Fig. 5C, 5D, 5E, 5F (common x-ticks 30, 40, 50%; y=0 at bottom border)."""
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    xt = find_ticks_mask(dark, r1, c0 + 15, c1 - 15, "x", -1, offset_lo=10, offset_hi=22)
    if len(xt) != 3:
        raise RuntimeError(f"Fig. 5C-F expected 3 x-ticks (30,40,50), found {len(xt)}: {xt}")
    x_axis = Axis("linear", list(zip(xt, [30.0, 40.0, 50.0])))

    yt = find_ticks_mask(dark, c0, r0 + 15, r1 - 15, "y", +1, offset_lo=10, offset_hi=22)
    if len(yt) != len(y_tick_vals):
        raise RuntimeError(f"Fig. 5C-F expected {len(y_tick_vals)} interior y-ticks, found {len(yt)}: {yt}")
    y_axis = Axis("linear", list(zip(yt, y_tick_vals)) + [(r1, 0.0)])
    return x_axis, y_axis, (c0, r0, c1, r1)


def color_mask_5cf(image: np.ndarray, series: str) -> np.ndarray:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    gray = (r + g + b) / 3.0
    if series == "black":
        chroma = np.maximum(r, np.maximum(g, b)) - np.minimum(r, np.minimum(g, b))
        return (gray < 95.0) & (chroma < 35.0)
    if series == "pink":
        return (r > 170.0) & (g < 165.0) & (b > 100.0) & ((r - g) > 50.0)
    raise ValueError(f"Unknown series {series!r}")


def fig5f_legend_box(box: tuple[float, float, float, float]) -> tuple[int, int, int, int]:
    """Locate the white-filled legend rectangle inside Fig. 5F dynamically relative to frame."""
    c0, r0, c1, r1 = box
    w = c1 - c0
    h = r1 - r0
    # Legend sits in the lower-right quadrant: ~[c0 + 0.36*w, c1 - 0.02*w] x [r0 + 0.58*h, r1 - 0.05*h]
    return (
        int(round(c0 + 0.365 * w)),
        int(round(r0 + 0.58 * h)),
        int(round(c1 - 0.02 * w)),
        int(round(r1 - 0.05 * h)),
    )


def read_scales(image: np.ndarray, y_tick_vals: list[float], series: str = "black",
                has_legend: bool = False) -> list[dict]:
    """Detect all individual scale 'x' markers in Fig. 5C-F."""
    x_axis, y_axis, box = calibrate_fig5cf(image, y_tick_vals)
    mask = color_mask_5cf(image, series)
    leg = fig5f_legend_box(box) if has_legend else None
    if series == "pink":
        pts = detect_crosses(mask, box, size=85, stroke=14.0, thr=0.54, min_dist=24,
                             margin=(14, 14, 14, 7), legend_box=leg, neg_radius=34.0, neg_weight=-0.5)
    else:
        pts = detect_crosses(mask, box, size=85, stroke=14.0, thr=0.72, min_dist=20,
                             margin=(18, 18, 18, 14), legend_box=leg, neg_radius=34.0, neg_weight=-0.6)

    half_px = 3.0 + max(y_axis.residual_px, 0.5)
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
            "note": f"cross score={sc:.2f}; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px",
        })
    return out


def read_means(image: np.ndarray, y_tick_vals: list[float], at: list[float],
               series: str = "black", has_legend: bool = False) -> list[dict]:
    """Read the 6 timepoint mean vertices of Fig. 5C-F at the requested timepoints `at`."""
    x_axis, y_axis, box = calibrate_fig5cf(image, y_tick_vals)
    mask = color_mask_5cf(image, series)
    leg = fig5f_legend_box(box) if has_legend else None
    if series == "pink":
        pts = detect_crosses(mask, box, size=85, stroke=14.0, thr=0.54, min_dist=24,
                             margin=(14, 14, 14, 7), legend_box=leg, neg_radius=34.0, neg_weight=-0.5)
    else:
        pts = detect_crosses(mask, box, size=85, stroke=14.0, thr=0.72, min_dist=20,
                             margin=(18, 18, 18, 14), legend_box=leg, neg_radius=34.0, neg_weight=-0.6)

    targets = [float(v) for v in at] if at else [29.70, 34.77, 37.13, 41.02, 47.01, 53.68]
    groups: list[list[tuple[float, float]]] = [[] for _ in targets]
    for rr, cc, _sc in pts:
        xv = float(x_axis.to_value(cc))
        yv = float(y_axis.to_value(rr))
        gi = int(np.argmin([abs(xv - tc) for tc in targets]))
        groups[gi].append((xv, yv))

    # For pink at the earliest timepoint (~29.7%), one pink cross is occluded behind a black cross;
    # read the exact dotted pink mean-curve start by linear extrapolation of the 4 clean pink dots in [30.8, 34.0]%.
    pink_t0_dot_val = None
    if series == "pink":
        c0, r0, c1, r1 = box
        m = np.zeros_like(mask)
        m[int(r0) + 12:int(r1) - 7, int(c0) + 14:int(c1) - 12] = mask[int(r0) + 12:int(r1) - 7, int(c0) + 14:int(c1) - 12]
        if leg is not None:
            m[leg[1]:leg[3], leg[0]:leg[2]] = False
        lab, n = ndimage.label(m)
        gap0_dots = []
        for i in range(1, n + 1):
            rr_d, cc_d = np.nonzero(lab == i)
            if 30 <= rr_d.size <= 380:
                xd = float(x_axis.to_value(float(cc_d.mean()) + 0.5))
                yd = float(y_axis.to_value(float(rr_d.mean()) + 0.5))
                if 30.8 <= xd <= 34.0:
                    gap0_dots.append((xd, yd))
        if len(gap0_dots) >= 2:
            sl, ic = np.polyfit([p[0] for p in gap0_dots], [p[1] for p in gap0_dots], 1)
            pink_t0_dot_val = float(sl * targets[0] + ic)

    out = []
    base_half_px = 3.5 + max(y_axis.residual_px, 0.5)
    for gi, tk in enumerate(targets):
        grp = groups[gi]
        if series == "pink" and gi == 0 and pink_t0_dot_val is not None:
            ym = pink_t0_dot_val
            rm = float(y_axis.to_pixel(ym))
            half_px = 7.0 + max(y_axis.residual_px, 0.5)
            note = f"pink dotted-line vertex (n=2 visible + 1 occluded behind black cross); cal res y={y_axis.residual_px:.2f} px"
        elif grp:
            ym = float(np.mean([p[1] for p in grp]))
            rm = float(y_axis.to_pixel(ym))
            half_px = base_half_px
            note = f"mean of n={len(grp)} crosses; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px"
        else:
            out.append({"x": tk, "y": None, "note": "no crosses found at timepoint"})
            continue
        cm = float(x_axis.to_pixel(tk))
        yv, lo, hi = interval_from_px(y_axis, rm, half_px)
        out.append({
            "x": tk,
            "y": yv,
            "lo": lo,
            "hi": hi,
            "px": [cm, rm],
            "px_lo": float(rm - half_px),
            "px_hi": float(rm + half_px),
            "note": note,
        })
    return out
