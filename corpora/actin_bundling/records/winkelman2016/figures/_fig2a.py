"""Calibrated reader for Fig. 2A (low-speed sedimentation competition: fascin vs alpha-actinin)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=140, row_frac=0.30, col_frac=0.30)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    yt_raw = C.find_ticks(dark, y_spine, (100, int(x_spine[1]) - 150), "y", side=+1, band=(4, 24), frac=0.4, max_width=40)
    yt = C.match_ticks_regular(yt_raw, 5, tol=6.0, min_step=250.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2])), max_residual_px=3.0)

    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 250, image.shape[1] - 100), "x", side=-1, band=(4, 24), frac=0.4, max_width=40)
    xt = C.match_ticks_regular(xt_raw, 4, tol=8.0, min_step=250.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.25, 0.50, 0.75, 1.00])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _find_legend_box(image: np.ndarray, x_spine, y_spine) -> tuple[int, int, int, int]:
    dark = image.astype(np.float32).mean(axis=2) < 120
    r_lim = int(0.38 * x_spine[0])
    c_lo = int(y_spine[2]) + 200
    sub = dark[:r_lim, c_lo:]
    h_rows = np.where(sub.sum(axis=1) > 400)[0]
    if len(h_rows) < 2:
        return (0, 0, 0, 0)
    r0, r1 = int(h_rows.min()), int(h_rows.max())
    v_cols = np.where(dark[r0 : r1 + 1, c_lo:].sum(axis=0) > 0.7 * (r1 - r0))[0] + c_lo
    if len(v_cols) < 2:
        return (0, 0, 0, 0)
    return (int(v_cols.min()) - 5, r0 - 5, int(v_cols.max()) + 5, r1 + 5)


def read_series(image: np.ndarray, series: str) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if series == "fascin":
        cmask = (R > 120) & (R - G > 50) & (R - B > 50)
    elif series == "actinin":
        cmask = (B > 120) & (B - R > 50) & (B - G > 20)
    else:
        raise ValueError(f"Unknown series {series!r}")

    row_min = 30
    row_max = int(x_spine[1]) - 2
    col_min = int(y_spine[2]) + 2
    cmask[:row_min, :] = False
    cmask[row_max + 1 :, :] = False
    cmask[:, :col_min] = False

    lc0, lr0, lc1, lr1 = _find_legend_box(image, x_spine, y_spine)
    if lc1 > lc0 and lr1 > lr0:
        cmask[lr0:lr1, lc0:lc1] = False

    # Fast 2x-downsampled closing across thin black polyline crossing markers
    blk2 = (R[::2, ::2] < 85) & (G[::2, ::2] < 85) & (B[::2, ::2] < 85)
    cm2 = cmask[::2, ::2]
    near2 = ndimage.binary_dilation(cm2, iterations=3)
    filled2 = ndimage.binary_closing(cm2 | (blk2 & near2), iterations=3)
    dist2 = ndimage.distance_transform_edt(filled2)
    loc_max = (dist2 == ndimage.maximum_filter(dist2, size=21)) & (dist2 >= 9.0)
    lbl, n = ndimage.label(loc_max)
    pts = []
    H, W = cmask.shape
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        r0, c0 = float(rr.mean()) * 2.0, float(cc.mean()) * 2.0
        rad = int(round(dist2[int(round(rr.mean())), int(round(cc.mean()))] * 2.0))
        y0, y1 = max(0, int(r0) - rad), min(H, int(r0) + rad + 1)
        x0, x1 = max(0, int(c0) - rad), min(W, int(c0) + rad + 1)
        sub_c = cmask[y0:y1, x0:x1]
        yy, xx = np.where(sub_c)
        d2 = (yy + y0 - r0) ** 2 + (xx + x0 - c0) ** 2
        keep = d2 <= (rad * 0.95) ** 2
        if keep.sum() > 20:
            pts.append((float((xx[keep] + x0).mean()) + 0.5, float((yy[keep] + y0).mean()) + 0.5))
        else:
            pts.append((c0 + 0.5, r0 + 0.5))
    pts.sort(key=lambda p: p[0])

    rad_mark = 28.0
    yy_grid, xx_grid = np.ogrid[:H, :W]
    out = []
    for idx, (c_mid, r_mid) in enumerate(pts):
        x_val = float(x_axis.to_value(c_mid))
        if abs(x_val) < 0.01:
            x_val = 0.0
        # Mask out OTHER markers so error bar trace never jumps into a neighbor
        stem_mask = cmask.copy()
        for j, (cj, rj) in enumerate(pts):
            if j == idx:
                continue
            y0, y1 = max(0, int(rj) - 36), min(H, int(rj) + 37)
            x0, x1 = max(0, int(cj) - 36), min(W, int(cj) + 37)
            sub_d2 = (yy_grid[y0:y1, :] - rj) ** 2 + (xx_grid[:, x0:x1] - cj) ** 2
            stem_mask[y0:y1, x0:x1] &= ~(sub_d2 <= (rad_mark + 2.0) ** 2)

        u_row, l_row = C.find_vertical_error_bar(
            stem_mask,
            c_mid,
            r_mid,
            row_min=row_min,
            row_max=row_max,
            stem_half_w=8,
            max_gap=10,
            cap_half_w=55,
            min_cap_width=40,
            start_up_offset=22,
            start_dn_offset=22,
        )
        if u_row is not None and (r_mid - u_row) < 30:
            u_row = None
        if l_row is not None and (l_row - r_mid) < 30:
            l_row = None
        if u_row is not None and l_row is None:
            l_row = r_mid + (r_mid - u_row)
        elif l_row is not None and u_row is None:
            u_row = r_mid - (l_row - r_mid)

        note = f"Fig. 2A {series}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
        if u_row is None and l_row is None:
            rd = C.make_point_reading(
                x_val,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=r_mid - rad_mark,
                lower_err_row=r_mid + rad_mark,
                mirror_lower=False,
                err_px=2.0 + y_axis.residual_px,
                note=note + "; SEM bar within marker radius",
            )
        else:
            rd = C.make_point_reading(
                x_val,
                c_mid,
                r_mid,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.0 + y_axis.residual_px,
                note=note,
            )
        out.append(rd)
    return out
