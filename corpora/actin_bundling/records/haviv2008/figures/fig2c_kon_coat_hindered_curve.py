"""Calibrated reader for Fig. 2c (p. 4 / 450): red dotted coat-hindered on-rate curve (2*k_on/D propto n_af*k_off)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    NATIVE_QUANT_PX,
    calibrate_fig2c,
    col_runs,
    detect_fig2c_marker_columns,
    find_hline_row_red,
    value_interval,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    hline_r, hline_r0, hline_r1 = find_hline_row_red(image, f)
    marker_cols = detect_fig2c_marker_columns(image, f)

    red = (R - G > 18.0) & (R - B > 15.0)
    red[max(0, hline_r0 - 14) : min(image.shape[0], hline_r1 + 15), :] = False
    for mc in marker_cols:
        red[:, max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30
    x_cross = 1.15

    # Purely geometric branch separation among non-horizontal red curves:
    # - Left of crossing (phi < 1.1): red dotted curve is the HIGHEST red curve above hline.
    # - Between crossing and hline (1.3 < phi < 2.3): red dotted curve is the LOWER red curve above hline.
    # - Right of hline crossing (phi > 2.8): red dotted curve is the ONLY red curve below hline (coalesced with black curve).
    def _sample_window(c: int, x: float) -> tuple[float, float, float, bool] | None:
        # Use a tighter horizontal window on the steep rising flank (phi < 0.5) so we only grab the dot at x
        half_w = 5 if x < 0.5 else 12
        c0 = max(int(f["lc1"]) + 6, c - half_w)
        c1 = min(int(f["rc0"]) - 6, c + half_w + 1)
        if x > 2.7:
            sub = red[hline_r1 + 15 : r_hi, c0:c1]
            rr, cc = np.nonzero(sub)
            if len(rr) == 0:
                return None
            rows = rr + hline_r1 + 15
            cols = cc + c0
            c_mid = int(round(float(cols.mean())))
            r_mid = float(rows.mean()) + 0.5
            comb = gray[max(hline_r1 + 15, int(r_mid) - 30) : min(r_hi, int(r_mid) + 31), c_mid] < 160.0
            c_idx = np.where(comb)[0]
            half_h = (int(c_idx[-1]) - int(c_idx[0]) + 1) / 2.0 if len(c_idx) else (int(rows.max()) - int(rows.min()) + 1) / 2.0
            return r_mid, float(cols.mean()) + 0.5, half_h, True

        proj = red[r_lo : hline_r0 - 14, c0:c1].any(axis=1)
        runs = col_runs(proj, r_lo)
        runs = [r for r in runs if (r[3] - r[2] + 1) <= 100]
        if not runs:
            return None
        if x < x_cross:
            if x < 0.5 and len(runs) >= 2:
                # Steep rising dotted curve above lower red dashed curve: keep topmost dot cluster
                # merging closely stacked dot pixels within 60 px of top run
                top_r0 = runs[0][2]
                cluster = [r for r in runs if r[2] - top_r0 <= 70]
                chosen = (
                    (cluster[0][2] + cluster[-1][3]) / 2.0 + 0.5,
                    (cluster[-1][3] - cluster[0][2] + 1) / 2.0,
                    cluster[0][2],
                    cluster[-1][3],
                )
            else:
                chosen = min(runs, key=lambda t: t[0])
        else:
            if len(runs) >= 2:
                chosen = max(runs, key=lambda t: t[0])
            else:
                rc_cand = runs[0][0]
                above_dash = red[r_lo : max(r_lo, int(rc_cand) - 25), max(int(f["lc1"]) + 6, c - 70) : min(int(f["rc0"]) - 6, c + 71)]
                if above_dash.any():
                    chosen = runs[0]
                else:
                    return None
        r0, r1 = chosen[2], chosen[3]
        sub = red[r0 : r1 + 1, c0:c1]
        rr, cc = np.nonzero(sub)
        r_mid = float(rr.mean() + r0 + 0.5)
        c_mid = float(cc.mean() + c0 + 0.5)
        half_h = (r1 - r0 + 1) / 2.0
        return r_mid, c_mid, half_h, False

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_window(c, float(x))
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "red dotted curve not found at this phi"})
            continue
        rc, cc, half_h, coalesced = hit
        hl = _sample_window(c - 28, float(x) - 0.11)
        hr = _sample_window(c + 28, float(x) + 0.11)
        if hl is not None and hr is not None and abs(hr[1] - hl[1]) > 5.0:
            signed_slope = (hr[0] - hl[0]) / (hr[1] - hl[1])
            slope = abs(signed_slope)
            rc_at_x = rc + signed_slope * (c + 0.5 - cc)
        else:
            slope = 0.0
            rc_at_x = rc
        extra = 4.0 if coalesced else 0.0
        err_px = (
            half_h
            + NATIVE_QUANT_PX
            + max(y_axis.residual_px, 1.0)
            + slope * (NATIVE_QUANT_PX + x_axis.residual_px + abs(c + 0.5 - cc) * 0.5)
            + extra
        )
        y_val, y_lo, y_hi = value_interval(y_axis, rc_at_x, err_px)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": max(0.0, y_lo),
            "hi": y_hi,
            "px": [float(c + 0.5), rc_at_x],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(max(0.0, y_lo))),
            "note": f"red dotted curve (2*kon/D={y_val:.4f}){'; coalesced with black solid curve' if coalesced else ''}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
