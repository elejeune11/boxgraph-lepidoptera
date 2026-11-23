"""Calibrated reader for SI Fig. S16 top (p. 24): two-cell pattern correlation C vs intercellular coupling strength beta (log10 x-axis)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, match_ticks_regular, verify_axis, value_interval


def read(image, at):
    rgb = image.astype(float)
    gray = rgb.mean(axis=2)
    dark = gray < 100.0
    H, W = dark.shape

    # Locate heavy black arrow axes
    h_lines = clusters(np.where(dark.mean(axis=1) > 0.50)[0])
    v_lines = clusters(np.where(dark.mean(axis=0) > 0.50)[0])
    if not h_lines or not v_lines:
        raise RuntimeError("Fig. S16: could not locate arrow axes")
    _, rb0, _ = h_lines[-1]
    c_left, _, cl1 = v_lines[0]

    # Faint gray Mathematica major/minor ticks right beside the black axes
    # Y-axis: 24 regular ticks from C = 1.20 down to 0.05 in steps of 0.05 (every 4th is major: 1.2, 1.0, 0.8, 0.6, 0.4, 0.2)
    yt_raw = [c[0] for c in clusters(np.where((gray[280 : rb0 - 20, cl1 + 8 : cl1 + 28] < 230.0).mean(axis=1) > 0.40)[0] + 280) if (c[2] - c[1] + 1) <= 10]
    yt_24 = match_ticks_regular(yt_raw, 24, tol=2.5)
    y_vals = [1.20 - 0.05 * k for k in range(24)]
    y_axis = verify_axis("linear", list(zip(yt_24, y_vals)), max_residual_px=1.5)

    # X-axis: 4 visible decade ticks at 10^-4, 10^-3, 10^-2, 10^-1 (10^-5 is covered by the first blue dot; 10^-6 is v_axis!)
    xt_raw = [c[0] for c in clusters(np.where((gray[rb0 - 28 : rb0 - 8, cl1 + 10 : W - 150] < 230.0).mean(axis=0) > 0.40)[0] + cl1 + 10) if (c[2] - c[1] + 1) <= 10]
    xt_4 = match_ticks_regular(xt_raw, 4, tol=3.0)
    x_axis = verify_axis(
        "log10",
        [(c_left, 1e-6), (xt_4[0], 1e-4), (xt_4[1], 1e-3), (xt_4[2], 1e-2), (xt_4[3], 1e-1)],
        max_residual_px=2.0,
    )

    # Detect all 13 blue circular scatter markers
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    blue = (B > 140.0) & (B - R > 35.0) & (R > 50.0) & (R < 150.0)
    m = np.zeros_like(blue)
    m[250 : rb0 + 10, cl1 + 10 : W - 100] = blue[250 : rb0 + 10, cl1 + 10 : W - 100]
    lab, n = ndimage.label(m)

    half_px = 3.0 + y_axis.residual_px
    out = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 400:
            continue
        r_c = float(rr.mean()) + 0.5
        c_c = float(cc.mean()) + 0.5
        x_val = float(x_axis.to_value(c_c))
        y_val, lo, hi = value_interval(y_axis, r_c, half_px)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_c, r_c],
            "px_lo": float(r_c - half_px),
            "px_hi": float(r_c + half_px),
            "note": f"blue marker ({len(rr)} px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
