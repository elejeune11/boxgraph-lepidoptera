"""Calibrated curve readers for Fig. 1d, Supp Fig. 1, Supp Fig. 7 (TIRF profiles), and Supp Fig. 8."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.ndimage import binary_dilation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C


def _cal_box_subset(
    image: np.ndarray,
    y_vals: list[float],
    x_vals: list[float],
    thresh: float = 160.0,
    col_frac: float = 0.45,
    include_y_ends=(True, True),
    include_x_ends=(True, True),
):
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=thresh, row_frac=0.45, col_frac=col_frac)
    sat = image.max(axis=2) - image.min(axis=2)
    dark = dark_all & (sat < 22)
    yt_l = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_left", band=(4, 20), margin=12, frac=0.4, max_width=24)
    yt_r = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_right", band=(4, 20), margin=12, frac=0.4, max_width=24)
    yc = yt_l + yt_r
    if include_y_ends[0]:
        yc.append(rt[0])
    if include_y_ends[1]:
        yc.append(rb[0])
    ya = C.verify_axis("linear", list(zip(C.match_ticks_regular(yc, len(y_vals), tol=8.0), y_vals)), max_residual_px=3.0)

    xt_t = C.find_inward_ticks(dark, cl, rt, cr, rb, "x_top", band=(4, 20), margin=12, frac=0.4, max_width=24)
    xt_b = C.find_inward_ticks(dark, cl, rt, cr, rb, "x_bot", band=(4, 20), margin=12, frac=0.4, max_width=24)
    xc = xt_t + xt_b
    if include_x_ends[0]:
        xc.append(cl[0])
    if include_x_ends[1]:
        xc.append(cr[0])
    xa = C.verify_axis("linear", list(zip(C.match_ticks_regular(xc, len(x_vals), tol=8.0), x_vals)), max_residual_px=3.0)
    return xa, ya, (cl, rt, cr, rb)


def read_fig1d(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 1d pyrene fluorescence curves (+Actin squares, +VCA+Arp2/3+profilin circles, +CP diamonds)."""
    xa, ya, (cl, rt, cr, rb) = _cal_box_subset(
        image,
        y_vals=[120.0, 80.0, 40.0, 0.0],
        x_vals=[0.0, 4.0, 8.0, 12.0, 16.0],
        thresh=150.0,
        col_frac=0.25,
        include_y_ends=(False, False),
        include_x_ends=(True, True),
    )
    g = (image.sum(axis=2, dtype=np.uint16) // 3).astype(np.int16)
    sat = image.max(axis=2) - image.min(axis=2)

    m_blk = (sat < 15) & (g <= 45)
    m_med = (sat < 15) & (g >= 95) & (g <= 155)
    m_lgt = (sat < 15) & (g >= 182) & (g <= 224)

    # Remove anti-aliased edge halos of darker series from lighter series
    if series == "cp":
        mask = m_blk
    elif series == "vca_arp_prof":
        mask = m_med & (~binary_dilation(m_blk, iterations=5))
    elif series == "actin":
        mask = m_lgt & (~binary_dilation(m_blk | m_med, iterations=5))
    else:
        raise ValueError(series)

    # Mask out legend text/symbols strictly below y = 94 and above y = 52 for x in [3.3, 12.0]
    r_lt, r_lb = int(ya.to_pixel(94.0)), int(ya.to_pixel(52.0))
    c_ll, c_lr = int(xa.to_pixel(3.3)), int(xa.to_pixel(12.0))
    mask[r_lt:r_lb, c_ll:c_lr] = False

    r_min, r_max = int(rt[2]) + 10, int(rb[1]) - 10
    mask[:r_min, :] = False
    mask[r_max:, :] = False
    mask[:, : int(cl[2]) + 2] = False
    mask[:, int(cr[1]) - 16 :] = False
    mask[int(rt[2]) : int(rt[2]) + 22, :] = False
    mask[int(rb[1]) - 22 : int(rb[1]), :] = False

    H, W = mask.shape
    out = []
    for x in at:
        c_f = float(xa.to_pixel(float(x)))
        c = max(int(cl[2]) + 4, min(int(cr[1]) - 18, int(round(c_f))))
        win = mask[r_min:r_max, max(0, c - 14) : min(W, c + 15)]
        rows_any = np.where(win.any(axis=1))[0]
        if len(rows_any) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": f"Fig1d {series}: no marker at t={x}"})
            continue
        cls = C.clusters(rows_any + r_min, gap=95)
        big = [t for t in cls if (t[2] - t[1] + 1) >= 20] or cls
        if series == "actin":
            best = max(big, key=lambda t: t[0])  # lowest marker curve (largest row)
        elif series == "vca_arp_prof":
            best = min(big, key=lambda t: t[0])  # highest marker curve (smallest row)
        else:
            best = max(big, key=lambda t: t[2] - t[1])
        r_mid = 0.5 * (best[1] + best[2]) + 0.5
        half_px = 6.0 + ya.residual_px
        out.append(
            C.make_point_reading(
                float(x),
                c_f,
                r_mid,
                ya,
                half_px=half_px,
                note=f"Fig1d {series} marker center, span={best[2]-best[1]+1}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    return out


def read_suppfig1(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Supp Fig. 1 pyrene fluorescence trace (60 nM CP + 13.5 uM profilin, 0..60 min)."""
    xa, ya, (cl, rt, cr, rb) = _cal_box_subset(
        image,
        y_vals=[5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0],
        x_vals=[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
        thresh=170.0,
    )
    g = (image.sum(axis=2, dtype=np.uint16) // 3).astype(np.int16)
    sat = image.max(axis=2) - image.min(axis=2)
    mask = (sat < 15) & (g >= 115) & (g <= 205)
    r_min, r_max = int(ya.to_pixel(0.8)), int(ya.to_pixel(-0.8))
    mask[:, : int(cl[2]) + 18] = False
    mask[:, int(cr[1]) - 18 :] = False

    out = []
    for x in at:
        sdx = 24 if (x <= 0.5 or x >= 59.5) else 6
        out.append(
            C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=sdx, note_prefix="SuppFig1 60nM, "
            )
        )
    return out


def read_suppfig8(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    """Read Supp Fig. 8 pyrene polymerization kinetics at +20, +40, +60, +80, +180 nM CP."""
    xa, ya, (cl, rt, cr, rb) = _cal_box_subset(
        image,
        y_vals=[120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0, -20.0],
        x_vals=[0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0],
        thresh=165.0,
    )
    g = (image.sum(axis=2, dtype=np.uint16) // 3).astype(np.int16)
    sat = image.max(axis=2) - image.min(axis=2)
    bands = {
        20: (198, 232),
        40: (158, 194),
        60: (110, 148),
        80: (60, 102),
        180: (0, 45),
    }
    lo_g, hi_g = bands[cp_nM]
    mask = (sat < 15) & (g >= lo_g) & (g <= hi_g)
    mask[int(ya.to_pixel(52.0)) :, int(xa.to_pixel(10.2)) :] = False
    r_min, r_max = int(rt[2]) + 6, int(rb[1]) - 6
    mask[:r_min, :] = False
    mask[r_max:, :] = False
    mask[:, : int(cl[2]) + 16] = False
    mask[:, int(cr[1]) - 16 :] = False
    mask[int(rt[2]) : int(rt[2]) + 18, :] = False
    mask[int(rb[1]) - 18 : int(rb[1]), :] = False

    out = []
    for x in at:
        sdx = 22 if (x <= 0.3 or x >= 14.7) else 8
        out.append(
            C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=sdx, note_prefix=f"SuppFig8 {cp_nM}nM, "
            )
        )
    return out


def read_suppfig7(image: np.ndarray, at: list[float], cp_nM: int, time_min: int) -> list[dict]:
    """Read Supp Fig. 7 TIRF domain intensity profile vs domain size (um) at given CP and time point (900 dpi)."""
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=215.0, row_frac=0.45, col_frac=0.45)
    ya = C.verify_axis("linear", [(rt[0], 250.0), (rb[0], 50.0)], max_residual_px=2.0)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = dark_all & (sat < 16)
    xt_t = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_top", band=(3, 18), margin=10, frac=0.35, max_width=20)
    xt_b = C.find_inward_ticks(ndark, cl, rt, cr, rb, "x_bot", band=(3, 18), margin=10, frac=0.35, max_width=20)
    xt_m = C.match_ticks_regular([cl[0]] + xt_t + xt_b, 4, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt_m, [0.0, 5.0, 10.0, 15.0])), max_residual_px=3.0)

    palette = {
        0: (224, 166, 179),
        1: (204, 126, 134),
        3: (178, 88, 92),
        6: (168, 70, 44),
        9: (118, 69, 49),
        20: (82, 36, 32),
    }
    is_red = (image[:, :, 0].astype(np.int16) - image[:, :, 2].astype(np.int16)) > 14
    dark_red = is_red & (image.min(axis=2) < 110)
    keys = list(palette.keys())
    dists = []
    for k in keys:
        t = palette[k]
        r = image[:, :, 0].astype(np.int32) - t[0]
        g = image[:, :, 1].astype(np.int32) - t[1]
        b = image[:, :, 2].astype(np.int32) - t[2]
        dists.append(r * r + g * g + b * b)
    stack = np.stack(dists, axis=0)
    best_i = np.argmin(stack, axis=0)
    best_d = np.min(stack, axis=0)
    idx_t = keys.index(time_min)
    max_d = 42 if time_min == 0 else 55
    mask = is_red & (best_i == idx_t) & (best_d <= max_d * max_d)
    if time_min in (0, 1):
        mask &= ~binary_dilation(dark_red, iterations=3)

    if cp_nM == 0:
        mask[: int(ya.to_pixel(135.0)), int(xa.to_pixel(12.0)) :] = False
    r_min, r_max = int(rt[2]) + 5, int(rb[1])
    mask[:r_min, :] = False
    mask[r_max + 1 :, :] = False
    mask[:, : int(cl[2]) + 5] = False
    mask[:, int(cr[1]) - 5 :] = False

    def pick_fn(rs, col):
        if time_min == 0 and cp_nM != 20:
            return max(rs, key=lambda t: t[0])  # lowest curve
        if time_min == 6 and cp_nM in (0, 20, 40):
            return min(rs, key=lambda t: t[0])  # topmost curve
        return max(rs, key=lambda t: t[1] - t[0])

    out = []
    for x in at:
        out.append(
            C.sample_solid_curve_y(
                mask,
                xa,
                ya,
                float(x),
                r_min,
                r_max,
                pick_fn=pick_fn,
                search_dx=12,
                extra_px=2.5,
                note_prefix=f"SuppFig7 {cp_nM}nM {time_min}min, ",
            )
        )
    return out
