"""Calibrated reader for Fig. S5D (Fascin) and Fig. S5E (alpha-Actinin) single-molecule survival curves."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray, panel: str):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 40), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    if panel == "S5d":
        yt_raw = C.find_ticks(dark, y_spine, (550, int(x_spine[1]) - 50), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
        yt = C.match_ticks_regular(yt_raw, 5, tol=8.0, min_step=400.0)
        y_axis = C.verify_axis("linear", list(zip(yt, [0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=3.5)
        xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 40, image.shape[1] - 200), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
        xt = C.match_ticks_regular(xt_raw, 5, tol=8.0, min_step=420.0)
        x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 10.0, 20.0, 30.0, 40.0])), max_residual_px=3.0)
    elif panel == "S5e":
        yt_raw = C.find_ticks(dark, y_spine, (650, int(x_spine[1]) - 200), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
        yt = C.match_ticks_regular(yt_raw, 4, tol=8.0, min_step=360.0)
        y_axis = C.verify_axis("linear", list(zip(yt, [0.8, 0.6, 0.4, 0.2])), max_residual_px=2.0)
        xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 8, image.shape[1] - 200), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
        xt = C.match_ticks_regular(xt_raw, 5, tol=8.0, min_step=440.0)
        x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 10.0, 20.0, 30.0, 40.0])), max_residual_px=2.5)
    else:
        raise ValueError(f"Unknown panel {panel!r}")
    return x_axis, y_axis, x_spine, y_spine


def _mask_legend(mask: np.ndarray, x_axis, y_axis):
    c0 = int(round(x_axis.to_pixel(11.5)))
    c1 = int(round(x_axis.to_pixel(46.0)))
    r0 = max(0, int(round(y_axis.to_pixel(1.08))))
    r1 = int(round(y_axis.to_pixel(0.24)))
    mask[r0:r1, c0:c1] = False


def read_series(image: np.ndarray, panel: str, series: str, kind: str, at: list[float]) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image, panel)
    # Fast uint8 color masks (avoid full-frame float32 allocations!)
    R = image[:, :, 0].astype(np.int16)
    G = image[:, :, 1].astype(np.int16)
    B = image[:, :, 2].astype(np.int16)

    if series == "control":
        cmask = (R > 150) & (R - G > 70) & (R - B > 70)
    elif series == "domains":
        cmask = (B > 170) & (G > 130) & (R < 130) & (B - R > 70)
    else:
        raise ValueError(f"Unknown series {series!r}")

    row_min = max(10, int(round(y_axis.to_pixel(1.08))))
    row_max = int(x_spine[1]) - 5
    col_min = max(10, int(y_spine[1]) - 60)
    cmask[:row_min, :] = False
    cmask[row_max + 1 :, :] = False
    cmask[:, :col_min] = False
    _mask_legend(cmask, x_axis, y_axis)

    H, W = cmask.shape

    if kind == "data":
        # 4x downsampled distance transform is ~16x faster and still gives 4-px peak seeds refined at full-res!
        sub4 = cmask[::4, ::4]
        dist4 = ndimage.distance_transform_edt(sub4)
        min_peak4 = 4.8 if (panel == "S5e" and series == "control") else 5.2
        loc_max = (dist4 == ndimage.maximum_filter(dist4, size=15)) & (dist4 >= min_peak4)
        lbl, n = ndimage.label(loc_max)
        pts = []
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 4.0, float(cc.mean()) * 4.0
            rad = int(round(dist4[int(round(rr.mean())), int(round(cc.mean()))] * 4.0))
            rad_box = max(rad + 8, 30)
            y0, y1 = max(0, int(r0) - rad_box), min(H, int(r0) + rad_box + 1)
            x0, x1 = max(0, int(c0) - rad_box), min(W, int(c0) + rad_box + 1)
            yy, xx = np.where(cmask[y0:y1, x0:x1])
            d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
            keep = d2 <= (max(rad * 1.1, 24.0)) ** 2
            if keep.sum() > 30:
                pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
        pts.sort(key=lambda p: p[0])
        # Deduplicate any double-hits within 70 px (~1.4 s, since markers are spaced by 2 s = ~100 px)
        dedup = []
        for p in pts:
            if dedup and abs(p[0] - dedup[-1][0]) < 70.0 and abs(p[1] - dedup[-1][1]) < 70.0:
                continue
            dedup.append(p)

        out = []
        for c_mid, r_mid in dedup:
            xv = float(x_axis.to_value(c_mid))
            if abs(xv) < 0.1:
                xv = 0.0
            out.append(
                C.make_point_reading(
                    xv,
                    c_mid,
                    r_mid,
                    y_axis,
                    err_px=3.0 + y_axis.residual_px,
                    note=f"Fig. {panel} {series} data; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
                )
            )
        return out

    if kind == "fit":
        # Fast 4x marker detection -> erase ONLY marker cores (tight radius so underlying curve stays intact!)
        sub4 = cmask[::4, ::4]
        dist4 = ndimage.distance_transform_edt(sub4)
        loc_max = (dist4 == ndimage.maximum_filter(dist4, size=13)) & (dist4 >= 5.2)
        lbl, n = ndimage.label(loc_max)
        curve_mask = cmask.copy()
        yy_g, xx_g = np.ogrid[:H, :W]
        for i in range(1, n + 1):
            rr, cc = np.where(lbl == i)
            r0, c0 = float(rr.mean()) * 4.0, float(cc.mean()) * 4.0
            rad = float(dist4[int(round(rr.mean())), int(round(cc.mean()))] * 4.0) + 4.0
            y0, y1 = max(0, int(r0 - rad) - 2), min(H, int(r0 + rad) + 3)
            x0, x1 = max(0, int(c0 - rad) - 2), min(W, int(c0 + rad) + 3)
            d2 = (yy_g[y0:y1, :] - r0) ** 2 + (xx_g[:, x0:x1] - c0) ** 2
            curve_mask[y0:y1, x0:x1] &= ~(d2 <= rad**2)

        # Notice in Fig. S5E domains (blue), the blue fit curve runs right BELOW the blue triangles at t >= 15 s!
        # So sampling at odd times (t = 1, 3, 5, ..., 47 s) sits cleanly BETWEEN the 2-second markers (t = 0, 2, 4, ...)!
        # For steep early times (t <= 6 s), allow vertical run height up to 160 px; for flat tails, pick the thin curve run.
        out = []
        for tv in at:
            tv = float(tv)
            max_h = 180 if tv <= 6.0 else 45

            def pick(rs, col, _mh=max_h, _tv=tv):
                cands = [s for s in rs if 2 <= (s[1] - s[0] + 1) <= _mh]
                if not cands:
                    return None
                if panel == "S5e" and series == "domains" and _tv >= 14.0:
                    # Blue fit curve runs along the bottom edge (largest row) below the blue triangles
                    return max(cands, key=lambda s: s[1])
                if panel == "S5e" and series == "control" and _tv >= 14.0:
                    # Red fit curve sits on/just above the red circles; pick run closest to curve band
                    return min(cands, key=lambda s: abs(0.5 * (s[0] + s[1]) - y_axis.to_pixel(0.06))) if _tv >= 25.0 else cands[0]
                return max(cands, key=lambda s: s[1] - s[0])

            rd = C.sample_solid_curve_y(
                curve_mask,
                x_axis,
                y_axis,
                tv,
                row_min,
                row_max,
                pick_fn=pick,
                search_dx=18,
                note_prefix=f"Fig. {panel} {series} exponential fit; ",
            )
            out.append(rd)
        return out
    raise ValueError(f"Unknown kind {kind!r}")
