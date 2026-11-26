"""Calibrated reader for Fig. 1E (TMR-fascin fluorescence accumulation over time)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=140, row_frac=0.30, col_frac=0.30)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    yt_raw = C.find_ticks(dark, y_spine, (260, int(x_spine[1]) - 120), "y", side=+1, band=(4, 24), frac=0.4, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 4, tol=6.0, min_step=180.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.8, 0.6, 0.4, 0.2])), max_residual_px=2.5)

    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 150, image.shape[1] - 50), "x", side=-1, band=(4, 24), frac=0.4, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 4, tol=6.0, min_step=200.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [20.0, 40.0, 60.0, 80.0])), max_residual_px=2.5)
    return x_axis, y_axis, x_spine, y_spine


def read_series(image: np.ndarray, kind: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if kind == "data":
        red = (R > 120) & (R - G > 55) & (R - B > 55)
        dist2 = ndimage.distance_transform_edt(red[::2, ::2])
        loc_max = (dist2 == ndimage.maximum_filter(dist2, size=23)) & (dist2 >= 12.5)
        lbl, n = ndimage.label(loc_max)
        pts = []
        H, W = red.shape
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
            rad = int(round(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0))
            y0, y1 = max(0, int(r0) - rad), min(H, int(r0) + rad + 1)
            x0, x1 = max(0, int(c0) - rad), min(W, int(c0) + rad + 1)
            yy, xx = np.where(red[y0:y1, x0:x1])
            d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
            keep = d2 <= (rad * 0.90) ** 2
            if keep.sum() > 20:
                pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
            else:
                pts.append((c0 + 0.5, r0 + 0.5))
        pts.sort(key=lambda p: p[0])
        out = []
        for c_mid, r_mid in pts:
            xv = float(x_axis.to_value(c_mid))
            if abs(xv) < 1.0:
                xv = 0.0
            out.append(
                C.make_point_reading(
                    xv,
                    c_mid,
                    r_mid,
                    y_axis,
                    err_px=2.5 + y_axis.residual_px,
                    note=f"Fig. 1E red marker; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                )
            )
        return out

    if kind == "fit":
        blk = (R < 90) & (G < 90) & (B < 90)
        row_min = 120
        row_max = int(x_spine[1]) - 30
        col_min = int(y_spine[2]) + 15
        blk[:row_min, :] = False
        blk[row_max + 1 :, :] = False
        blk[:, :col_min] = False
        return [
            C.sample_solid_curve_y(
                blk, x_axis, y_axis, float(xv), row_min, row_max, search_dx=10, note_prefix="Fig. 1E sigmoidal fit; "
            )
            for xv in at
        ]
    raise ValueError(f"Unknown kind {kind!r}")
