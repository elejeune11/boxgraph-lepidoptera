"""Calibrated reader for Fig. 2d (p. 4 / 450): dotted black coat-hindered saturating N_peak curve (misreferenced as Fig. 3d)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    NATIVE_QUANT_PX,
    calibrate_fig2d,
    col_runs,
    detect_fig2d_marker_columns,
    value_interval,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2d(image)
    gray = image.astype(float).mean(axis=2)
    marker_cols = detect_fig2d_marker_columns(image, f)

    dark = gray < 165.0
    for mc in marker_cols:
        dark[:, max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30

    # Purely topological branch selection among Fig. 2d's three black curves:
    # - phi < 1.15: dotted curve is the HIGHEST (smallest row) curve.
    # - 1.35 < phi < 2.25: dotted curve is the MIDDLE of the 3 curves (below dashed, above solid).
    # - 2.55 < phi < 3.60: dotted curve is the LOWEST curve (below solid).
    # - phi > 3.60: dotted curve is the ONLY curve in the entire panel.
    def _sample_window(c: int, x: float) -> tuple[float, float, float] | None:
        half_w = 5 if x < 0.5 else 14
        c0 = max(int(f["lc1"]) + 6, c - half_w)
        c1 = min(int(f["rc0"]) - 6, c + half_w + 1)
        proj = dark[r_lo:r_hi, c0:c1].any(axis=1)
        runs = col_runs(proj, r_lo)
        runs = [r for r in runs if (r[3] - r[2] + 1) <= 140]
        if not runs:
            return None
        runs.sort(key=lambda t: t[0])  # top to bottom
        if x > 3.60:
            chosen = runs[0]  # only curve in panel
        elif x > 2.50:
            chosen = runs[-1]  # lowest curve (below solid)
        elif x < 1.15:
            chosen = runs[0]  # highest curve
        else:
            # 1.35 < x < 2.25: lowest run runs[-1] is the continuous solid curve;
            # the dotted curve is the compact dot run immediately above it (runs[-2])!
            if len(runs) >= 3:
                chosen = runs[-2]
            elif len(runs) == 2:
                # Verify runs[0] is the dotted curve (dashed above within +/-45 cols)
                above_dash = dark[r_lo : max(r_lo, int(runs[0][2]) - 50), max(int(f["lc1"]) + 6, c - 45) : min(int(f["rc0"]) - 6, c + 46)]
                chosen = runs[0] if above_dash.any() else None
                if chosen is None:
                    return None
            else:
                return None
        r0, r1 = chosen[2], chosen[3]
        sub = dark[r0 : r1 + 1, c0:c1]
        rr, cc = np.nonzero(sub)
        r_mid = float(rr.mean() + r0 + 0.5)
        c_mid = float(cc.mean() + c0 + 0.5)
        half_h = (r1 - r0 + 1) / 2.0
        return r_mid, c_mid, half_h

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_window(c, float(x))
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "dotted black curve not found at this phi"})
            continue
        rc, cc, half_h = hit
        hl = _sample_window(c - 28, float(x) - 0.11)
        hr = _sample_window(c + 28, float(x) + 0.11)
        if hl is not None and hr is not None and abs(hr[1] - hl[1]) > 5.0:
            signed_slope = (hr[0] - hl[0]) / (hr[1] - hl[1])
            slope = abs(signed_slope)
            rc_at_x = rc + signed_slope * (c + 0.5 - cc)
        else:
            slope = 0.0
            rc_at_x = rc
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + slope * (NATIVE_QUANT_PX + x_axis.residual_px + abs(c + 0.5 - cc) * 0.5)
        y_val, y_lo, y_hi = value_interval(y_axis, rc_at_x, err_px)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": max(0.0, y_lo),
            "hi": y_hi,
            "px": [float(c + 0.5), rc_at_x],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(max(0.0, y_lo))),
            "note": f"black dotted curve N_peak={y_val:.2f}, half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
