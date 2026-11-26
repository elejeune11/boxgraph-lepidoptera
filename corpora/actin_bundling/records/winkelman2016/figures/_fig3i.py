"""Calibrated reader for Fig. 3I (TMR-fascin fluorescence/um vs [TMR-Fascin])."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

import _calib as C


def _calibrate(image: np.ndarray):
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=100, row_frac=0.35, col_frac=0.35, row_range=(30, 1650), col_range=(300, 2100)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Fig. 3I has no y tick marks on the spine; calibrate from printed y-label digit centers '8', '4', '0'
    # (note: top label '1.2' is a misprint for '12' and lies near top edge of crop)
    c_lbl_hi = int(y_spine[1]) - 15
    c_lbl_lo = max(0, c_lbl_hi - 140)
    txt_rows = (gray[150 : int(x_spine[1]) - 20, c_lbl_lo:c_lbl_hi] < 140).any(axis=1)
    tcs = [c for c in C.clusters(np.where(txt_rows)[0] + 150, gap=30) if (c[2] - c[1] + 1) >= 60]
    if len(tcs) != 3:
        raise RuntimeError(f"Fig. 3I: expected 3 y-label digit clusters (8, 4, 0), got {tcs}")
    y_axis = C.verify_axis("linear", [(tcs[0][0], 8.0), (tcs[1][0], 4.0), (tcs[2][0], 0.0)], max_residual_px=3.0)

    # Inward x ticks: 0, 40, 80, 120 (160 is at right end of x-spine)
    xt_raw = C.find_ticks(dark, x_spine, (int(y_spine[2]) + 40, image.shape[1] - 50), "x", side=-1, band=(4, 24), frac=0.35, max_width=45)
    xt = C.match_ticks_regular(xt_raw, 4, tol=8.0, min_step=300.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 40.0, 80.0, 120.0])), max_residual_px=3.0)
    return x_axis, y_axis, x_spine, y_spine


def _find_legend_box(image: np.ndarray, y_spine) -> tuple[int, int, int, int]:
    dark = image.astype(np.float32).mean(axis=2) < 100
    c_lo = int(y_spine[2]) + 150
    sub = dark[650:1220, c_lo:]
    h_rows = np.where(sub.sum(axis=1) > 600)[0] + 650
    if len(h_rows) < 2:
        return (0, 0, 0, 0)
    r0, r1 = int(h_rows.min()), int(h_rows.max())
    v_cols = np.where(dark[r0 : r1 + 1, c_lo:].sum(axis=0) > 0.7 * (r1 - r0))[0] + c_lo
    if len(v_cols) < 2:
        return (0, 0, 0, 0)
    return (int(v_cols.min()) - 6, r0 - 6, int(v_cols.max()) + 6, r1 + 6)


def read_series(image: np.ndarray, series: str) -> list[dict]:
    x_axis, y_axis, x_spine, y_spine = _calibrate(image)
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)

    if series == "pure_fascin":
        cmask = (B > 150) & (G > 130) & (R < 80)
    elif series == "fascin_domains":
        cmask = (R > 160) & (R - G > 80) & (R - B > 80)
    elif series == "actinin_domains":
        cmask = (B > 150) & (R < 80) & (G < 100)
    else:
        raise ValueError(f"Unknown series {series!r}")

    row_min = 20
    row_max = int(x_spine[1]) - 5
    col_min = int(y_spine[2]) + 5
    cmask[:row_min, :] = False
    cmask[row_max + 1 :, :] = False
    cmask[:, :col_min] = False

    lc0, lr0, lc1, lr1 = _find_legend_box(image, y_spine)
    if lc1 > lc0 and lr1 > lr0:
        cmask[lr0:lr1, lc0:lc1] = False

    lbl, n = ndimage.label(cmask)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        if len(rr) > 250:
            pts.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5))
    pts.sort(key=lambda p: p[0])

    out = []
    for c_mid, r_mid in pts:
        xv = float(x_axis.to_value(c_mid))
        if abs(xv) < 2.0:
            xv = 0.0
        out.append(
            C.make_point_reading(
                xv,
                c_mid,
                r_mid,
                y_axis,
                err_px=4.0 + y_axis.residual_px,
                note=f"Fig. 3I {series}; y calibrated on printed 0/4/8 labels (top '1.2' misprint for 12); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            )
        )
    return out
