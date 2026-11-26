"""Calibrated reader for Fig. 2d (p. 4 / 450): fitted peak bundle size N_peak vs phi (5 black squares + ±44% error bars)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, _group_runs, calibrate_fig2d, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2d(image)
    gray = image.astype(float).mean(axis=2)

    # Detect the 5 grey-filled black squares on 1/4-res via 8x8 square morphological opening
    sub = image[::4, ::4].astype(float).mean(axis=2)
    sint = np.zeros(sub.shape, dtype=bool)
    sint[int(f["tr1"]) // 4 + 3 : int(f["br0"]) // 4 - 3, int(f["lc1"]) // 4 + 3 : int(f["rc0"]) // 4 - 3] = True
    sq8 = np.ones((8, 8), dtype=bool)
    d_sq = ndimage.binary_opening((sub < 185.0) & sint, structure=sq8)
    lbl, _ = ndimage.label(d_sq)

    dark = gray < 165.0
    r_min_int = int(f["tr1"]) + 10
    r_max_int = int(f["br0"]) - 30  # stay above bottom inward ticks

    center_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0
    cap_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0

    readings: list[dict] = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        m = lbl[sl] == i
        cy_q, cx_q = ndimage.center_of_mass(m)
        c0 = int(round((cx_q + sl[1].start) * 4.0 + 2.0))
        r0 = int(round((cy_q + sl[0].start) * 4.0 + 2.0))

        # Refine square center on full-resolution 56x56 patch (square body is ~40x40 px)
        patch = gray[r0 - 28 : r0 + 29, c0 - 28 : c0 + 29] < 185.0
        rw = patch.sum(axis=1) >= 24
        cw = patch.sum(axis=0) >= 24
        body = patch & rw[:, None] & cw[None, :]
        if int(body.sum()) > 200:
            by, bx = ndimage.center_of_mass(body)
            cx_f = float(c0 - 28 + bx + 0.5)
            cy_f = float(r0 - 28 + by + 0.5)
        else:
            cx_f, cy_f = float(c0), float(r0)

        c_int = int(round(cx_f))
        r_int = int(round(cy_f))

        # Trace vertical error bar stem upward and downward from the marker
        r_up = r_int
        while r_up > r_min_int and dark[r_up - 1, c_int - 6 : c_int + 7].any():
            r_up -= 1
        r_dn = r_int
        while r_dn < r_max_int and dark[r_dn + 1, c_int - 6 : c_int + 7].any():
            r_dn += 1

        # Locate outermost horizontal T-caps (symmetric left & right protrusion around stem end)
        def _cap_near_end(r_edge: int, direction: str) -> float:
            if direction == "up":
                rows_scan = range(r_edge, min(r_edge + 24, r_int - 24))
            else:
                rows_scan = range(max(r_edge - 23, r_int + 25), r_edge + 1)
            cap_rows = [
                r for r in rows_scan
                if dark[r, c_int - 15 : c_int - 4].any()
                and dark[r, c_int + 5 : c_int + 16].any()
                and int(dark[r, c_int - 15 : c_int + 16].sum()) >= 18
            ]
            segs = _group_runs(np.array(cap_rows, dtype=int), gap=2)
            if segs:
                chosen = segs[0] if direction == "up" else segs[-1]
                return float((chosen[0] + chosen[1]) / 2.0 + 0.5)
            return float(r_edge + 4.0) if direction == "up" else float(r_edge - 4.0)

        top_cap = _cap_near_end(r_up, "up")
        bot_cap = _cap_near_end(r_dn, "down")

        x_val = float(x_axis.to_value(cx_f))
        y_val, y_lo, y_hi = value_interval(y_axis, cy_f, center_err_px)
        b_hi = float(y_axis.to_value(top_cap))
        b_lo = float(y_axis.to_value(bot_cap))
        b_hi_out = float(y_axis.to_value(top_cap - cap_err_px))
        b_lo_out = float(y_axis.to_value(bot_cap + cap_err_px))

        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": max(0.0, min(b_lo_out, b_hi_out)),
            "bar_hi_outer": max(b_lo_out, b_hi_out),
            "px": [cx_f, cy_f],
            "px_lo": top_cap,
            "px_hi": bot_cap,
            "note": f"black square N_peak={y_val:.2f}, caps [{b_lo:.2f},{b_hi:.2f}] (-{(y_val-b_lo)/y_val*100:.1f}%/+{(b_hi-y_val)/y_val*100:.1f}%); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })

    readings.sort(key=lambda r: r["x"])
    return readings
