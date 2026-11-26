"""Calibrated bar-chart readers for Supp Fig. 2, Supp Fig. 4 (3 panels), and Supp Fig. 10."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C


def _cal_y_inward(image: np.ndarray, y_vals: list[float], include_top: bool = True, include_bot: bool = True, thresh: float = 160.0):
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=thresh, row_frac=0.45, col_frac=0.45)
    sat = image.max(axis=2) - image.min(axis=2)
    dark = dark_all & (sat < 25)
    yt_l = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=12, frac=0.45, max_width=24)
    yt_r = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_right", band=(4, 18), margin=12, frac=0.45, max_width=24)
    cands = yt_l + yt_r
    if include_top:
        cands.append(rt[0])
    if include_bot:
        cands.append(rb[0])
    ticks = C.match_ticks_regular(cands, len(y_vals), tol=6.0)
    ya = C.verify_axis("linear", list(zip(ticks, y_vals)), max_residual_px=2.5)
    return ya, (cl, rt, cr, rb)


def _voronoi_masks(image: np.ndarray, palette: dict[str, tuple[int, int, int]], sat_min: int = 25, max_dist: int = 120) -> dict[str, np.ndarray]:
    """Assign each colored pixel to its nearest palette entry (eliminating anti-aliasing gaps between stacked bars)."""
    sat = image.max(axis=2) - image.min(axis=2)
    # Include colored pixels or dark non-white pixels within max_dist of a palette entry
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
    valid = (sat >= sat_min) & (best_d <= max_dist * max_dist)
    return {k: valid & (best_i == i) for i, k in enumerate(keys)}


def _segment_reading(
    x_val: float,
    col_center: float,
    row_top: float | None,
    row_bot: float | None,
    ya,
    is_stacked_upper: bool = False,
    note_prefix: str = "",
) -> dict:
    """Build calibrated reading for a bar or stacked bar segment (height = y_top - y_bot)."""
    h_px = 2.0 + ya.residual_px
    dy_px = abs(1.0 / ya.b)
    if row_top is None or row_bot is None:
        y0_row = float(ya.to_pixel(0.0))
        err_val = float(h_px / dy_px)
        return {
            "x": float(x_val),
            "y": 0.0,
            "lo": 0.0,
            "hi": err_val,
            "px": [float(col_center), y0_row],
            "px_lo": y0_row - h_px,
            "px_hi": y0_row + h_px,
            "note": f"{note_prefix}0% (absent segment), cal_res={ya.residual_px:.2f}px",
        }

    y_top = float(ya.to_value(row_top))
    y_bot = float(ya.to_value(row_bot))
    height = max(0.0, y_top - y_bot)
    edge_px = (1.414 * h_px) if is_stacked_upper else h_px
    err_val = float(edge_px / dy_px)
    lo = max(0.0, height - err_val)
    hi = height + err_val
    if is_stacked_upper:
        return {
            "x": float(x_val),
            "y": height,
            "lo": lo,
            "hi": hi,
            "px": [float(col_center), 0.5 * (row_top + row_bot)],
            "px_lo": float(row_top),
            "px_hi": float(row_bot),
            "note": f"{note_prefix}stacked [{y_bot:.1f}..{y_top:.1f}%], cal_res={ya.residual_px:.2f}px",
        }
    return {
        "x": float(x_val),
        "y": height,
        "lo": lo,
        "hi": hi,
        "px": [float(col_center), float(row_top)],
        "px_lo": float(row_top - h_px),
        "px_hi": float(row_top + h_px),
        "note": f"{note_prefix}bar top {y_top:.2f}%, cal_res={ya.residual_px:.2f}px",
    }


def _find_color_span_in_band(mask: np.ndarray, c0: int, c1: int, r_min: int, r_max: int, frac: float = 0.5):
    prof = mask[r_min:r_max, c0:c1].mean(axis=1)
    rs = C.runs_1d(prof >= frac, gap=2)
    if not rs:
        return None, None
    a, b = max(rs, key=lambda t: t[1] - t[0])
    return float(a + r_min), float(b + r_min + 1)


def read_suppfig2(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Supp Fig. 2 vesicle fraction bars at c_CP = 0, 20 nM."""
    ya, (cl, rt, cr, rb) = _cal_y_inward(image, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0], include_top=True, include_bot=True)
    r_min, r_max = int(rt[2]) + 4, int(rb[0])
    w_box = cr[0] - cl[0]
    masks = _voronoi_masks(image, {"org": (242, 146, 0), "blu": (15, 112, 183)}, sat_min=30, max_dist=130)
    # Exclude legend in upper middle (y > 72, x in [0.25..0.70] of box)
    r_leg = int(ya.to_pixel(72.0))
    for m in masks.values():
        m[:r_leg, int(cl[0] + 0.25 * w_box) : int(cl[0] + 0.70 * w_box)] = False
    m_all = masks["org"] | masks["blu"]

    col_prof = m_all[r_min:r_max, int(cl[2]) + 10 : int(cr[1]) - 10].sum(axis=0) > 15
    col_runs = [
        (a + int(cl[2]) + 10, b + int(cl[2]) + 10)
        for a, b in C.runs_1d(col_prof, gap=4)
        if (b - a) > 0.08 * w_box
    ]
    if len(col_runs) != 2:
        raise RuntimeError(f"Expected 2 colored membrane bars in Supp Fig 2, found {col_runs}")

    cp_vals = [0.0, 20.0]
    out = []
    if series == "bulk":
        gray = image.sum(axis=2, dtype=np.uint16)
        sat = image.max(axis=2) - image.min(axis=2)
        ndark = (gray < 420) & (sat < 25)
        bulk_regions = [
            (int(cl[2]) + 15, col_runs[0][0] - 15),
            (int(cl[0] + 0.55 * w_box), col_runs[1][0] - 15),
        ]
        for cp, (x0, x1) in zip(cp_vals, bulk_regions):
            v_prof = ndark[int(rb[1]) - 150 : int(rb[1]) - 10, x0:x1].mean(axis=0)
            walls = C.clusters(np.where(v_prof > 0.7)[0] + x0, gap=3)
            bw0, bw1 = int(walls[0][2]) + 10, int(walls[-1][1]) - 10
            h_prof = ndark[r_min : int(rb[1]) - 10, bw0:bw1].mean(axis=1)
            tops = C.clusters(np.where(h_prof > 0.7)[0] + r_min, gap=3)
            out.append(_segment_reading(cp, 0.5 * (bw0 + bw1), tops[0][0], rb[0], ya, False, f"SuppFig2 bulk {cp:.0f}nM, "))
        return out

    for cp, (bx0, bx1) in zip(cp_vals, col_runs):
        c0, c1 = bx0 + 15, bx1 - 15
        c_mid = 0.5 * (bx0 + bx1)
        if series == "domains":
            rt_s, _ = _find_color_span_in_band(masks["org"], c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig2 domains {cp:.0f}nM, "))
        elif series == "cortex":
            rt_s, rb_s = _find_color_span_in_band(masks["blu"], c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb_s, ya, True, f"SuppFig2 cortex {cp:.0f}nM, "))
        elif series == "membrane_total":
            rt_s, _ = _find_color_span_in_band(m_all, c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig2 membrane {cp:.0f}nM, "))
        else:
            raise ValueError(series)
    return out


def read_suppfig4_domains(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Supp Fig. 4 Domains panel at c_CP = 40, 60, 120, 180 nM."""
    ya, (cl, rt, cr, rb) = _cal_y_inward(image, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0], include_top=False, include_bot=True)
    r_min, r_max = int(rt[2]) + 10, int(rb[0])
    w_box = cr[0] - cl[0]

    masks = _voronoi_masks(image, {"red": (163, 24, 48), "org": (222, 126, 12), "yel": (253, 215, 22)}, sat_min=30, max_dist=130)
    r_leg = int(ya.to_pixel(70.0))
    c_leg = int(cl[0] + 0.55 * w_box)
    for m in masks.values():
        m[:r_leg, c_leg:] = False
    m_all = masks["red"] | masks["org"] | masks["yel"]

    col_prof = m_all[r_min:r_max, int(cl[2]) + 10 : int(cr[1]) - 10].sum(axis=0) > 15
    bar_runs = [
        (a + int(cl[2]) + 10, b + int(cl[2]) + 10)
        for a, b in C.runs_1d(col_prof, gap=5)
        if (b - a) > 0.08 * w_box
    ]
    if len(bar_runs) != 4:
        raise RuntimeError(f"Expected 4 domain bars in Supp Fig 4, found {bar_runs}")

    out = []
    for cp, (bx0, bx1) in zip([40.0, 60.0, 120.0, 180.0], bar_runs):
        c0, c1 = bx0 + 20, bx1 - 20
        c_mid = 0.5 * (bx0 + bx1)
        if series == "protrusion":
            rt_s, _ = _find_color_span_in_band(masks["red"], c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig4 protrusion {cp:.0f}nM, "))
        elif series == "concave_domain":
            rt_s, rb_s = _find_color_span_in_band(masks["org"], c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb_s, ya, True, f"SuppFig4 concave_dom {cp:.0f}nM, "))
        elif series == "flat_domain":
            rt_s, rb_s = _find_color_span_in_band(masks["yel"], c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb_s, ya, True, f"SuppFig4 flat_dom {cp:.0f}nM, "))
        elif series == "domains_total":
            rt_s, _ = _find_color_span_in_band(m_all, c0, c1, r_min, r_max)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig4 domains_tot {cp:.0f}nM, "))
        else:
            raise ValueError(series)
    return out


def read_suppfig4_cortex(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Supp Fig. 4 Closed Cortex panel at c_CP = 40, 60, 120, 180 nM."""
    ya, (cl, rt, cr, rb) = _cal_y_inward(image, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0], include_top=False, include_bot=True)
    r_min, r_max = int(rt[2]) + 10, int(rb[0])
    w_box = cr[0] - cl[0]

    masks = _voronoi_masks(image, {"nav": (38, 54, 126), "cyn": (31, 182, 185)}, sat_min=25, max_dist=130)
    r_leg = int(ya.to_pixel(75.0))
    c_leg = int(cl[0] + 0.50 * w_box)
    for m in masks.values():
        m[:r_leg, c_leg:] = False
    m_all = masks["nav"] | masks["cyn"]

    col_prof = m_all[r_min:r_max, int(cl[2]) + 10 : int(cr[1]) - 10].sum(axis=0) > 40
    tall_runs = [
        (a + int(cl[2]) + 10, b + int(cl[2]) + 10)
        for a, b in C.runs_1d(col_prof, gap=5)
        if (b - a) > 0.08 * w_box
    ]
    if len(tall_runs) != 3:
        raise RuntimeError(f"Expected 3 prominent cortex bars (40,120,180 nM), found {tall_runs}")
    b0, b2, b3 = tall_runs
    bw = int(0.5 * ((b0[1] - b0[0]) + (b2[1] - b2[0])))
    c1_mid = int(0.5 * (0.5 * (b0[0] + b0[1]) + 0.5 * (b2[0] + b2[1])))
    b1 = (c1_mid - bw // 2, c1_mid + bw // 2)
    bar_runs = [b0, b1, b2, b3]

    out = []
    for cp, (bx0, bx1) in zip([40.0, 60.0, 120.0, 180.0], bar_runs):
        c0, c1 = bx0 + 20, bx1 - 20
        c_mid = 0.5 * (bx0 + bx1)
        if series == "flat_cortex":
            rt_s, _ = _find_color_span_in_band(masks["nav"], c0, c1, r_min, r_max, frac=0.4)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig4 flat_cor {cp:.0f}nM, "))
        elif series == "concave_cortex":
            rt_s, rb_s = _find_color_span_in_band(masks["cyn"], c0, c1, r_min, r_max, frac=0.4)
            out.append(_segment_reading(cp, c_mid, rt_s, rb_s, ya, True, f"SuppFig4 concave_cor {cp:.0f}nM, "))
        elif series == "cortex_total":
            rt_s, _ = _find_color_span_in_band(m_all, c0, c1, r_min, r_max, frac=0.4)
            out.append(_segment_reading(cp, c_mid, rt_s, rb[0] if rt_s is not None else None, ya, False, f"SuppFig4 cor_tot {cp:.0f}nM, "))
        else:
            raise ValueError(series)
    return out


def read_suppfig4_bulk(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Supp Fig. 4 Bulk Polymerization panel at c_CP = 40, 60, 120, 180 nM."""
    ya, (cl, rt, cr, rb) = _cal_y_inward(image, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0], include_top=False, include_bot=True)
    r_min, r_max = int(rt[2]) + 10, int(rb[1])
    gray = image.sum(axis=2, dtype=np.uint16)
    ndark = gray < 420
    v_prof = ndark[r_max - 45 : r_max - 6, int(cl[2]) + 20 : int(cr[1]) - 20].mean(axis=0)
    walls = [c[0] + int(cl[2]) + 20 for c in C.clusters(np.where(v_prof > 0.65)[0], gap=4)]
    if len(walls) != 6:
        raise RuntimeError(f"Expected 6 vertical bar walls (3 bulk bars for 60,120,180 nM), found {walls}")
    b60 = (int(walls[0]), int(walls[1]))
    b120 = (int(walls[2]), int(walls[3]))
    b180 = (int(walls[4]), int(walls[5]))
    step = 0.5 * ((b120[0] - b60[0]) + (b180[0] - b120[0]))
    b40 = (int(round(b60[0] - step)), int(round(b60[1] - step)))
    bar_runs = [b40, b60, b120, b180]

    out = []
    for cp, (bx0, bx1) in zip([40.0, 60.0, 120.0, 180.0], bar_runs):
        c0, c1 = bx0 + 20, bx1 - 20
        c_mid = 0.5 * (bx0 + bx1)
        h_prof = ndark[r_min : r_max - 4, c0:c1].mean(axis=1)
        tops = C.clusters(np.where(h_prof > 0.7)[0] + r_min, gap=3)
        if not tops:
            out.append(_segment_reading(cp, c_mid, None, None, ya, False, f"SuppFig4 bulk {cp:.0f}nM, "))
        else:
            out.append(_segment_reading(cp, c_mid, tops[0][0], rb[0], ya, False, f"SuppFig4 bulk {cp:.0f}nM, "))
    return out


def read_suppfig10(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Supp Fig. 10 contractile response stacked bar chart at c_CP = 60, 120 nM (+5 uM NMM II)."""
    ya, (cl, rt, cr, rb) = _cal_y_inward(image, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0], include_top=True, include_bot=True)
    r_min, r_max = int(rt[2]) + 10, int(rb[0])
    w_box = cr[0] - cl[0]

    masks = _voronoi_masks(
        image,
        {
            "invagination": (51, 44, 131),
            "fission": (23, 156, 207),
            "deformation": (165, 190, 105),
            "no_deformation": (237, 228, 30),
        },
        sat_min=25,
        max_dist=130,
    )
    c_leg = int(cl[0] + 0.68 * w_box)
    for m in masks.values():
        m[:, c_leg:] = False
    m_all = masks["invagination"] | masks["fission"] | masks["deformation"] | masks["no_deformation"]

    col_prof = m_all[r_min:r_max, int(cl[2]) + 10 : c_leg].sum(axis=0) > 100
    bar_runs = [
        (a + int(cl[2]) + 10, b + int(cl[2]) + 10)
        for a, b in C.runs_1d(col_prof, gap=6)
        if (b - a) > 0.12 * w_box
    ]
    if len(bar_runs) != 2:
        raise RuntimeError(f"Expected 2 stacked bars in Supp Fig 10, found {bar_runs}")

    is_upper = series != "invagination"
    mask = masks[series]
    out = []
    for cp, (bx0, bx1) in zip([60.0, 120.0], bar_runs):
        c0, c1 = bx0 + 20, bx1 - 20
        c_mid = 0.5 * (bx0 + bx1)
        rt_s, rb_s = _find_color_span_in_band(mask, c0, c1, r_min, r_max)
        if not is_upper:
            rb_s = rb[0] if rt_s is not None else None
        out.append(_segment_reading(cp, c_mid, rt_s, rb_s, ya, is_upper, f"SuppFig10 {series} {cp:.0f}nM, "))
    return out
