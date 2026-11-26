"""Shared calibrated readers for Fig. 5c, 5d, 5e and Supp Fig. 11 (infinite-pool PDE model curves)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C


def _match_regular_subset(endpoints: tuple[float, float], candidates: list[float], vals: list[float], tol: float = 5.0):
    """Match all vals whose expected linear position between endpoints has a detected tick (endpoints always included)."""
    p0, p1 = endpoints
    v0, v1 = vals[0], vals[-1]
    marks = [(p0, v0)]
    for v in vals[1:-1]:
        expected = p0 + (v - v0) / (v1 - v0) * (p1 - p0)
        if candidates:
            near = min(candidates, key=lambda c: abs(c - expected))
            if abs(near - expected) <= tol:
                marks.append((near, v))
    marks.append((p1, v1))
    if len(marks) < min(4, len(vals)):
        raise RuntimeError(f"Too few regular ticks matched ({len(marks)}/{len(vals)}) from {candidates}")
    return marks


def _cal_fig5_axes(image: np.ndarray, y_vals: list[float], x_vals: list[float], thresh: float = 150.0):
    dark_all, gray, cl, rt, cr, rb = C.find_box(image, thresh=thresh, row_frac=0.45, col_frac=0.45)
    sat = image.max(axis=2) - image.min(axis=2)
    dark = dark_all & (sat < 18)
    yt_l = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_left", band=(4, 18), margin=15, frac=0.45, max_width=24)
    yt_r = C.find_inward_ticks(dark, cl, rt, cr, rb, "y_right", band=(4, 18), margin=15, frac=0.45, max_width=24)
    y_marks = _match_regular_subset((rt[0], rb[0]), yt_l + yt_r, y_vals, tol=5.0)
    ya = C.verify_axis("linear", y_marks, max_residual_px=2.5)

    xt_t = C.find_inward_ticks(dark, cl, rt, cr, rb, "x_top", band=(4, 18), margin=15, frac=0.45, max_width=24)
    xt_b = C.find_inward_ticks(dark, cl, rt, cr, rb, "x_bot", band=(4, 18), margin=15, frac=0.45, max_width=24)
    x_marks = _match_regular_subset((cl[0], cr[0]), xt_t + xt_b, x_vals, tol=5.0)
    xa = C.verify_axis("linear", x_marks, max_residual_px=2.5)
    return xa, ya, (cl, rt, cr, rb)


def _gray_curve_mask(image: np.ndarray, target_g: int, tol_g: int = 18, sat_max: int = 14) -> np.ndarray:
    g3 = image.sum(axis=2, dtype=np.int16)
    sat = image.max(axis=2) - image.min(axis=2)
    return (sat <= sat_max) & (np.abs(g3 - 3 * target_g) <= 3 * tol_g)


def _sqdist_rgb(image: np.ndarray, target: tuple[int, int, int]) -> np.ndarray:
    r = image[:, :, 0].astype(np.int32) - target[0]
    g = image[:, :, 1].astype(np.int32) - target[1]
    b = image[:, :, 2].astype(np.int32) - target[2]
    return r * r + g * g + b * b


def read_fig5c(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    xa, ya, (cl, rt, cr, rb) = _cal_fig5_axes(
        image,
        y_vals=[500.0, 400.0, 300.0, 200.0, 100.0, 0.0],
        x_vals=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
    )
    target_g = {40: 196, 80: 128, 120: 77}[cp_nM]
    mask = _gray_curve_mask(image, target_g, tol_g=16)
    c_leg = int(xa.to_pixel(22.0))
    r_leg_top = int(ya.to_pixel(400.0))
    r_leg_bot = int(ya.to_pixel(250.0))
    mask[r_leg_top:r_leg_bot, c_leg:] = False
    r_min = int(rt[2]) + 6
    r_max = int(rb[1]) + 4
    mask[:r_min, :] = False
    mask[:, : int(cl[2]) + 2] = False
    mask[:, int(cr[1]) - 18 :] = False

    out = []
    for x in at:
        sdx = 24 if x >= 39.5 else 6
        out.append(
            C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=sdx, note_prefix=f"Fig5c {cp_nM}nM, "
            )
        )
    return out


def read_suppfig11(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    xa, ya, (cl, rt, cr, rb) = _cal_fig5_axes(
        image,
        y_vals=[500.0, 400.0, 300.0, 200.0, 100.0, 0.0],
        x_vals=[0.0, 20.0, 40.0, 60.0, 80.0, 100.0],
        thresh=155.0,
    )
    target_g = {40: 199, 80: 135, 120: 87}[cp_nM]
    mask = _gray_curve_mask(image, target_g, tol_g=18)
    r_min = int(rt[2]) + 6
    r_max = int(rb[1]) + 4
    mask[:r_min, :] = False
    mask[:, : int(cl[2]) + 2] = False
    mask[:, int(cr[1]) - 18 :] = False
    mask[int(rt[2]) : int(rt[2]) + 22, :] = False
    mask[int(rb[1]) - 20 : int(rb[1]) - 6, :] = False

    out = []
    for x in at:
        out.append(
            C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=24 if x >= 99.0 else 6, note_prefix=f"SuppFig11 {cp_nM}nM, "
            )
        )
    return out


def read_fig5d(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    xa, ya, (cl, rt, cr, rb) = _cal_fig5_axes(
        image,
        y_vals=[100.0, 80.0, 60.0, 40.0, 20.0, 0.0],
        x_vals=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
    )
    d40 = _sqdist_rgb(image, (113, 192, 226))
    d80 = _sqdist_rgb(image, (47, 122, 187))
    d120 = _sqdist_rgb(image, (35, 39, 89))
    is_blue = (image[:, :, 2].astype(np.int16) - image[:, :, 0].astype(np.int16)) > 15
    tol2 = 45 * 45
    if cp_nM == 40:
        mask = is_blue & (d40 < tol2) & (d40 < d80) & (d40 < d120)
    elif cp_nM == 80:
        mask = is_blue & (d80 < tol2) & (d80 < d40) & (d80 < d120)
    else:
        mask = is_blue & (d120 < tol2) & (d120 < d40) & (d120 < d80)

    mask[: int(ya.to_pixel(65.0)), int(xa.to_pixel(22.0)) :] = False
    r_min = int(rt[2]) + 6
    r_max = int(rb[1]) + 4
    mask[:r_min, :] = False
    mask[:, : int(cl[2]) + 2] = False
    mask[:, int(cr[1]) - 18 :] = False

    out = []
    for x in at:
        sdx = 24 if x >= 39.5 else 6
        out.append(
            C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=sdx, note_prefix=f"Fig5d {cp_nM}nM, "
            )
        )
    return out


def read_fig5d_ss(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    """Read horizontal dashed steady-state line Arp2/3*_ss for 80 or 120 nM CP in Fig. 5d."""
    xa, ya, (cl, rt, cr, rb) = _cal_fig5_axes(
        image,
        y_vals=[100.0, 80.0, 60.0, 40.0, 20.0, 0.0],
        x_vals=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
    )
    mask = _gray_curve_mask(image, 47, tol_g=35, sat_max=15)
    c0 = int(xa.to_pixel(17.5))
    c1 = int(xa.to_pixel(22.0))
    if cp_nM == 120:
        r0, r1 = int(ya.to_pixel(23.0)), int(ya.to_pixel(16.0))
    else:
        r0, r1 = int(ya.to_pixel(15.0)), int(ya.to_pixel(8.0))
    prof = mask[r0:r1, c0:c1].mean(axis=1)
    rows = np.where(prof > 0.25)[0]
    if len(rows) == 0:
        raise RuntimeError(f"Could not find dashed Arp2/3*_ss line for {cp_nM} nM CP")
    r_mid = float(rows.mean()) + r0 + 0.5
    thick = float(rows[-1] - rows[0] + 1)
    half_px = max(0.5 * thick, 1.5) + 0.5 + ya.residual_px

    xs = at if at else [30.0]
    out = []
    for x in xs:
        c_f = float(xa.to_pixel(float(x)))
        out.append(
            C.make_point_reading(
                float(x),
                c_f,
                r_mid,
                ya,
                half_px=half_px,
                note=f"Fig5d Arp2/3*_ss {cp_nM}nM dashed line, thick={thick:.0f}px, cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_fig5e(image: np.ndarray, at: list[float], kind: str, cp_nM: int) -> list[dict]:
    """Read Fig. 5e branching (solid green) or capping (dashed orange/red) event rate curves."""
    xa, ya, (cl, rt, cr, rb) = _cal_fig5_axes(
        image,
        y_vals=[160.0, 140.0, 120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0],
        x_vals=[0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
    )
    targets = {
        ("branch", 40): (128, 184, 27),
        ("branch", 80): (11, 139, 52),
        ("branch", 120): (36, 36, 36),
        ("cap", 40): (243, 146, 0),
        ("cap", 80): (225, 53, 36),
        ("cap", 120): (151, 32, 24),
    }
    dists = {k: _sqdist_rgb(image, v) for k, v in targets.items()}
    key = (kind, cp_nM)
    tol = 45 if key != ("branch", 120) else 28
    mask = dists[key] <= tol * tol
    for other_k, other_d in dists.items():
        if other_k != key:
            mask &= dists[key] < other_d

    mask[: int(ya.to_pixel(52.0)), int(xa.to_pixel(20.5)) :] = False
    r_min = int(rt[2]) + 6
    r_max = int(rb[1]) + 4
    mask[:r_min, :] = False
    mask[:, : int(cl[2]) + 2] = False
    mask[:, int(cr[1]) - 18 :] = False
    mask[int(rb[1]) - 20 : int(rb[1]) - 4, :] = False
    mask[int(rt[2]) + 4 : int(rt[2]) + 20, :] = False

    out = []
    for x in at:
        sdx = 24 if x >= 39.5 else 6
        if kind == "branch":
            r = C.sample_solid_curve_y(
                mask, xa, ya, float(x), r_min, r_max, search_dx=sdx, note_prefix=f"Fig5e branch {cp_nM}nM, "
            )
        else:
            r = C.sample_dashed_curve_y(
                mask, xa, ya, float(x), r_min, r_max, max_gap_px=75, note_prefix=f"Fig5e cap {cp_nM}nM, "
            )
        out.append(r)
    return out
