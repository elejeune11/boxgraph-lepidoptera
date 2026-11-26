"""Calibrated reader for Fig. 5I (high-speed co-sedimentation: bound bundler vs [Tropomyosin])."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 80), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 1.0, 0.8, 0.6, 0.4, 0.2 (exact 273.0 px spacing; 0.0 at x-spine center)
    yt_raw = C.find_ticks(dark, y_spine, (150, int(x_spine[1]) - 30), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 5, tol=8.0, min_step=220.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2])), max_residual_px=2.0)

    # Inward x ticks: 0, 2, 4, 6, 8, 10 (12 at right end of x-spine)
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 20, image.shape[1] - 50), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 6, tol=8.0, min_step=200.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _inset_region(x_axis, y_axis, image_shape) -> tuple[int, int, int, int]:
    """Return (c0, r0, c1, r1) covering the legend + gel inset in the upper-right of Fig. 5I."""
    # Inset sits above y ~ 0.72 and right of x ~ 2.0 uM
    c0 = int(round(x_axis.to_pixel(1.9)))
    r1 = int(round(y_axis.to_pixel(0.68)))
    return (c0, 0, image_shape[1], r1)


def read_series(image: np.ndarray, series: str, kind: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if series == "actinin":
        cmask = (B > 130) & (B - R > 60) & (B - G > 15)
    elif series == "fascin":
        cmask = (R > 160) & (R - G > 80) & (R - B > 80)
    else:
        raise ValueError(f"Unknown series {series!r}")

    row_min = 60
    row_max = int(x_spine[1]) - 5
    col_min = int(y_spine[2]) + 5
    cmask[:row_min, :] = False
    cmask[row_max + 1 :, :] = False
    cmask[:, :col_min] = False

    ic0, ir0, ic1, ir1 = _inset_region(x_axis, y_axis, image.shape)
    cmask[ir0:ir1, ic0:ic1] = False

    if kind == "data":
        # Separate filled circular markers (radius ~26-30 px) from thin curve (width ~10-14 px) via distance transform
        dist2 = ndimage.distance_transform_edt(cmask[::2, ::2])
        min_peak = 11.5 if series == "actinin" else 10.0
        loc_max = (dist2 == ndimage.maximum_filter(dist2, size=19)) & (dist2 >= min_peak)
        lbl, n = ndimage.label(loc_max)
        pts = []
        H, W = cmask.shape
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
            rad = int(round(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0))
            y0, y1 = max(0, int(r0) - rad), min(H, int(r0) + rad + 1)
            x0, x1 = max(0, int(c0) - rad), min(W, int(c0) + rad + 1)
            yy, xx = np.where(cmask[y0:y1, x0:x1])
            d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
            keep = d2 <= (rad * 0.90) ** 2
            if keep.sum() > 20:
                pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
            else:
                pts.append((c0 + 0.5, r0 + 0.5))
        pts.sort(key=lambda p: p[0])

        # Note: at x=0, red marker sits on top of blue marker at y=1.0; if blue missed x=0 due to red occlusion,
        # detect the red marker at x=0 as well? Wait: only report visibly detected pixels (check visually first!)
        out = []
        for c_mid, r_mid in pts:
            xv = float(x_axis.to_value(c_mid))
            if abs(xv) < 0.08:
                xv = 0.0
            out.append(
                C.make_point_reading(
                    xv,
                    c_mid,
                    r_mid,
                    y_axis,
                    err_px=2.5 + y_axis.residual_px,
                    note=f"Fig. 5I {series} data; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                )
            )
        return out

    if kind == "fit":
        # Detect circular markers first and mask them out so the curve tracer never jumps onto a dot!
        dist2 = ndimage.distance_transform_edt(cmask[::2, ::2])
        min_peak = 11.5 if series == "actinin" else 10.0
        loc_max = (dist2 == ndimage.maximum_filter(dist2, size=19)) & (dist2 >= min_peak)
        lbl, n = ndimage.label(loc_max)
        curve_mask = cmask.copy()
        H, W = cmask.shape
        yy_g, xx_g = np.ogrid[:H, :W]
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
            rad = float(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0) + 6.0
            y0, y1 = max(0, int(r0 - rad) - 2), min(H, int(r0 + rad) + 3)
            x0, x1 = max(0, int(c0 - rad) - 2), min(W, int(c0 + rad) + 3)
            d2 = (yy_g[y0:y1, :] - r0) ** 2 + (xx_g[:, x0:x1] - c0) ** 2
            curve_mask[y0:y1, x0:x1] &= ~(d2 <= rad**2)

        xs_dense = np.linspace(0.05, 11.98, 180)
        tc, tr = [], []
        for xv in xs_dense:
            c = int(round(x_axis.to_pixel(xv)))
            if c < col_min or c >= W:
                continue
            rs = C.runs_1d(curve_mask[row_min : row_max + 1, c])
            thin = [(a + row_min, b + row_min) for a, b in rs if 3 <= (b - a + 1) <= 65]
            if len(thin) == 1:
                tc.append(float(c) + 0.5)
                tr.append(0.5 * (thin[0][0] + thin[0][1]) + 0.5)
        tc = np.array(tc)
        tr = np.array(tr)
        med_tr = ndimage.median_filter(tr, size=7)
        keep = np.abs(tr - med_tr) < 40
        tc, tr = tc[keep], tr[keep]

        out = []
        for xv in at:
            c_f = float(x_axis.to_pixel(float(xv)))
            pred_r = float(np.interp(c_f, tc, tr))

            def pick(rs, col, _pr=pred_r):
                cands = [s for s in rs if (s[1] - s[0] + 1) <= 70 and abs(0.5 * (s[0] + s[1]) + 0.5 - _pr) <= 50]
                if not cands:
                    return None
                return min(cands, key=lambda s: abs(0.5 * (s[0] + s[1]) + 0.5 - _pr))

            rd = C.sample_solid_curve_y(
                curve_mask,
                x_axis,
                y_axis,
                float(xv),
                row_min,
                row_max,
                pick_fn=pick,
                search_dx=22,
                note_prefix=f"Fig. 5I {series} exponential fit; ",
            )
            out.append(rd)
        return out
    raise ValueError(f"Unknown kind {kind!r}")
