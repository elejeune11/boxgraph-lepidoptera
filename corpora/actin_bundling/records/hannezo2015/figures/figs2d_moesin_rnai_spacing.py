"""Calibrated reader for SI Fig. S2D (p. 7): taenidial ring spacing bar chart for Wild-type vs Moesin RNAi with error bars."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, verify_axis, value_interval


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 120.0
    H, W = dark.shape

    h_lines = clusters(np.where(dark.mean(axis=1) > 0.50)[0])
    v_lines = clusters(np.where(dark[:, : int(0.35 * W)].mean(axis=0) > 0.50)[0])
    if not h_lines or not v_lines:
        raise RuntimeError("Fig. S2D: could not locate axes")
    r_base, rb0, rb1 = h_lines[-1]
    c_left, cl0, cl1 = v_lines[0]

    # Leftward y-ticks at 2.5, 2.0, 1.5, 1.0, 0.5, 0.0
    yt = [c[0] for c in clusters(np.where(dark[: rb1 + 5, cl0 - 28 : cl0 - 8].mean(axis=1) > 0.5)[0]) if (c[2] - c[1] + 1) <= 20]
    if len(yt) != 6:
        raise RuntimeError(f"Fig. S2D: expected 6 y-ticks (2.5..0.0), found {len(yt)}: {yt}")
    y_axis = verify_axis("linear", list(zip(yt, [2.5, 2.0, 1.5, 1.0, 0.5, 0.0])), max_residual_px=1.5)

    # Two downward x-ticks for Wild-type (1) and Moesin RNAi (2)
    xt = [c[0] for c in clusters(np.where(dark[rb1 + 8 : rb1 + 32, cl1 + 20 :].mean(axis=0) > 0.5)[0] + cl1 + 20) if (c[2] - c[1] + 1) <= 20]
    if len(xt) != 2:
        raise RuntimeError(f"Fig. S2D: expected 2 x-ticks, found {len(xt)}: {xt}")

    specs = [
        (1, xt[0], "Wild-type", 0.0, 100.0),
        (2, xt[1], "Moesin RNAi", 120.0, 200.0),
    ]
    half_px = 2.0 + y_axis.residual_px
    out = []
    for idx, xc, label, g_lo, g_hi in specs:
        xc_i = int(round(xc))
        # Bar top measured outside error-bar cap width (xc-250..xc-160)
        bar_mask = (gray[int(yt[0]) - 20 : rb0 - 10, xc_i - 250 : xc_i - 160] >= g_lo) & (
            gray[int(yt[0]) - 20 : rb0 - 10, xc_i - 250 : xc_i - 160] <= g_hi
        )
        r_bar = float(np.where(bar_mask.mean(axis=1) > 0.7)[0].min() + (int(yt[0]) - 20))

        # Horizontal error-bar cap above bar top
        cap_mask = (gray[int(yt[0]) - 30 : int(r_bar) - 15, xc_i - 80 : xc_i + 80] >= g_lo) & (
            gray[int(yt[0]) - 30 : int(r_bar) - 15, xc_i - 80 : xc_i + 80] <= g_hi
        )
        cap_rows = np.where(cap_mask.mean(axis=1) > 0.5)[0] + (int(yt[0]) - 30)
        if len(cap_rows) == 0:
            raise RuntimeError(f"Fig. S2D: could not locate error-bar cap for {label}")
        r_cap = float(cap_rows.mean()) + 0.5

        y_val, lo, hi = value_interval(y_axis, r_bar, half_px)
        bar_hi = float(y_axis.to_value(r_cap))
        bar_hi_out = float(y_axis.to_value(r_cap - half_px))
        # Lower error-bar half hidden inside filled bar -> mirror of upper half
        bar_lo = 2.0 * y_val - bar_hi
        bar_lo_out = 2.0 * y_val - bar_hi_out
        out.append({
            "x": float(idx),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo_out,
            "bar_hi_outer": bar_hi_out,
            "px": [float(xc), r_bar],
            "px_lo": r_cap,
            "px_hi": r_bar,
            "note": f"{label}; lower error bar hidden inside bar (mirrored); cal_res={y_axis.residual_px:.2f}px",
        })
    return out
