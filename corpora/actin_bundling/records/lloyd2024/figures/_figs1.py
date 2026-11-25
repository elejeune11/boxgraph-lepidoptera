"""Shared reader for SI Fig. 1E-H 1D electron-density correlation curves (PDF p. 23, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D
from _common import run_centers


def read_corr_panel(image: np.ndarray, at: list[float],
                    y_top: float, y_bot: float,
                    y_bot_uncert: float = 0.0) -> list[dict]:
    """Calibrate axes from black frame & ticks, trace blue correlation curve G(r)."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0
    H, W = dark.shape

    col_frac = dark.mean(axis=0)
    col_frac[:30] = 0.0
    c_lines = run_centers(np.where(col_frac > 0.55)[0], gap=5)
    r_lines = run_centers(np.where(dark.mean(axis=1) > 0.55)[0], gap=5)
    c0, c1 = c_lines[0], c_lines[-1]
    r0, r1 = r_lines[0], r_lines[-1]

    # Detect outer y-ticks (3 ticks: top, middle, bottom) and x-ticks (6 ticks: 0..1500 nm)
    yt = run_centers(np.where(dark[max(0, int(r0) - 20):min(H, int(r1) + 20), int(round(c0)) - 15])[0]
                     + max(0, int(r0) - 20), gap=4)
    xt = run_centers(np.where(dark[int(round(r1)) + 15, max(0, int(c0) - 20):min(W, int(c1) + 20)])[0]
                     + max(0, int(c0) - 20), gap=4)
    if len(xt) != 6:
        raise RuntimeError(f"Expected 6 x-ticks (0..1500), got {xt}")
    if len(yt) != 3:
        raise RuntimeError(f"Expected 3 y-ticks, got {yt}")

    xa = D.Axis("linear", list(zip(xt, [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0])))
    # Top and bottom ticks define exact span [y_bot, y_top]; middle tick is midpoint (rounded to 3 dp on print)
    y_mid = 0.5 * (y_top + y_bot)
    ya = D.Axis("linear", [(yt[0], y_top), (yt[1], y_mid), (yt[2], y_bot)])

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    bmask = (B - R > 45.0) & (B - G > 45.0)
    # Restrict to interior of plot box and exclude top-right legend box
    bmask[:int(r0) + 2, :] = False
    bmask[int(r1) - 2:, :] = False
    bmask[:, :int(c0) + 1] = False
    bmask[:, int(c1) - 1:] = False
    bmask[:int(r0 + 0.25 * (r1 - r0)), int(c0 + 0.65 * (c1 - c0)):] = False

    out = []
    for xv in at:
        col_f = xa.to_pixel(float(xv))
        col = int(round(col_f))
        if abs(xv) < 1e-6:
            # At r = 0 nm the curve starts at top-left corner (c0, r0) = (0, 1.000)
            row_c = yt[0]
            half_px = max(2.0, ya.residual_px)
            y_val = y_top
            lo = ya.to_value(row_c + half_px)
            hi = ya.to_value(row_c - half_px)
            out.append({
                "x": float(xv), "y": float(y_val),
                "lo": float(min(lo, hi)), "hi": float(max(lo, hi)),
                "px": [float(col_f), float(row_c)],
                "px_lo": float(row_c - half_px), "px_hi": float(row_c + half_px),
                "note": f"origin (rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px)",
            })
            continue
        rows = np.where(bmask[:, max(0, col - 1):min(W, col + 2)].any(axis=1))[0]
        if len(rows) == 0:
            out.append({"x": float(xv), "y": None, "note": f"no blue curve pixels at col {col}"})
            continue
        row_c = float(rows.mean() + 0.5)
        half_px = max(0.5 * (rows.max() - rows.min() + 1), 2.0) + max(ya.residual_px, 0.5)
        y_val = float(ya.to_value(row_c))
        a = float(ya.to_value(row_c - half_px))
        b = float(ya.to_value(row_c + half_px))
        # Add honest tick-label rounding uncertainty (scales linearly with fractional depth below 1.000)
        frac = max(0.0, min(1.0, (row_c - yt[0]) / max(yt[2] - yt[0], 1.0)))
        extra = y_bot_uncert * frac
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": float(min(a, b) - extra),
            "hi": float(max(a, b) + extra),
            "px": [float(col_f), row_c],
            "px_lo": float(row_c - half_px),
            "px_hi": float(row_c + half_px),
            "note": f"rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px",
        })
    return out
