"""Calibrated readers for Fig. 5f (finite-pool cone barbed-end density) and Supp Fig. 14 (total VCA)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C


def _cal_fig5f_colorbar(image: np.ndarray, tick_vals: list[float], x_band: tuple[int, int]):
    """Calibrate a single Fig. 5f colorbar located within column range x_band."""
    H, W = image.shape[:2]
    c0, c1 = x_band
    sub = image[:, c0:c1]
    gray = sub.sum(axis=2, dtype=np.uint16)
    sat = sub.max(axis=2) - sub.min(axis=2)
    ndark = (gray < 360) & (sat < 25)
    vs = C.clusters(np.where(ndark.sum(axis=0) > 0.50 * H)[0] + c0, gap=4)
    if len(vs) < 2:
        raise RuntimeError(f"Fig. 5f colorbar borders not found in {x_band}: {vs}")
    cl, cr = vs[0], vs[-1]
    cr_rows = np.where((image[:, cr[1] : cr[2] + 1].sum(axis=2, dtype=np.uint16) < 360).mean(axis=1) > 0.5)[0]
    rt_row, rb_row = float(cr_rows[0]) + 0.5, float(cr_rows[-1]) + 0.5
    c_mid = int(0.5 * (cl[0] + cr[0]))
    mid_col = image[int(rt_row) + 5 : int(rb_row) - 5, c_mid, :].astype(float)
    band = image[int(rt_row) + 5 : int(rb_row) - 5, cr[1] - 14 : cr[1] - 3, :].astype(float)
    diff = np.linalg.norm(band - mid_col[:, None, :], axis=2).mean(axis=1)
    ts = [c[0] for c in C.clusters(np.where(diff > 24.0)[0] + int(rt_row) + 5, gap=3) if (c[2] - c[1] + 1) <= 20]
    ts.append(rb_row)
    matched = C.match_ticks_regular(ts, len(tick_vals), tol=6.0)
    ya = C.verify_axis("linear", list(zip(matched, tick_vals)), max_residual_px=2.5)
    return ya, cl[0], cr[0], rt_row, rb_row


def read_fig5f_colorbar_peak(image: np.ndarray, at: list[float], cp_nM: int) -> list[dict]:
    """Read the auto-scaled colorbar upper limit (= global B_peak over r, t) in one Fig. 5f panel."""
    H, W = image.shape[:2]
    tick_map = {
        40: [400.0, 300.0, 200.0, 100.0, 0.0],
        80: [250.0, 200.0, 150.0, 100.0, 50.0, 0.0],
        120: [150.0, 100.0, 50.0, 0.0],
    }
    ya, cl, cr, rt_row, rb_row = _cal_fig5f_colorbar(image, tick_map[cp_nM], (int(0.75 * W), W))
    c_mid = 0.5 * (cl + cr)
    half_px = 2.0 + ya.residual_px
    return [
        C.make_point_reading(
            float(cp_nM),
            c_mid,
            rt_row,
            ya,
            half_px=half_px,
            note=f"Fig5f {cp_nM}nM colorbar top (auto-scaled global B_peak), cal_res={ya.residual_px:.2f}px",
        )
    ]


def read_fig5f_all_peaks(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read global B_peak from the 3 colorbars in a wide crop spanning all three Fig. 5f panels."""
    H, W = image.shape[:2]
    # Three panels sit in horizontal thirds [0..W/3], [W/3..2W/3], [2W/3..W]; each colorbar is in the right 25% of its third
    specs = [
        (40, [400.0, 300.0, 200.0, 100.0, 0.0], (int(0.24 * W), int(0.36 * W))),
        (80, [250.0, 200.0, 150.0, 100.0, 50.0, 0.0], (int(0.56 * W), int(0.69 * W))),
        (120, [150.0, 100.0, 50.0, 0.0], (int(0.88 * W), W)),
    ]
    out = []
    for cp_nM, tvals, xb in specs:
        ya, cl, cr, rt_row, rb_row = _cal_fig5f_colorbar(image, tvals, xb)
        out.append(
            C.make_point_reading(
                float(cp_nM),
                0.5 * (cl + cr),
                rt_row,
                ya,
                half_px=2.0 + ya.residual_px,
                note=f"Fig5f {cp_nM}nM colorbar top (global B_peak), cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def _cal_sfig14_colorbar(image: np.ndarray, tick_vals: list[float], x_band: tuple[int, int]):
    """Calibrate a single Supp Fig. 14 colorbar using saturation drop at inward grey ticks."""
    H, W = image.shape[:2]
    c0, c1 = x_band
    sub = image[:, c0:c1]
    gray = sub.sum(axis=2, dtype=np.uint16)
    sat_u8 = sub.max(axis=2) - sub.min(axis=2)
    ndark = (gray < 380) & (sat_u8 < 35)
    vs = C.clusters(np.where(ndark.sum(axis=0) > 0.50 * H)[0] + c0, gap=4)
    if len(vs) < 2:
        raise RuntimeError(f"SuppFig14 colorbar borders not found in {x_band}: {vs}")
    cl, cr = vs[0], vs[-1]
    cr_rows = np.where((image[:, cr[1] : cr[2] + 1].sum(axis=2, dtype=np.uint16) < 380).mean(axis=1) > 0.5)[0]
    rt_row, rb_row = float(cr_rows[0]) + 0.5, float(cr_rows[-1]) + 0.5
    c_mid = int(0.5 * (cl[0] + cr[0]))
    sat_f = image.max(axis=2).astype(float) - image.min(axis=2).astype(float)
    mid_s = sat_f[int(rt_row) + 5 : int(rb_row) - 5, c_mid]
    band_s = sat_f[int(rt_row) + 5 : int(rb_row) - 5, cr[1] - 10 : cr[1] - 1].mean(axis=1)
    sdrop = mid_s - band_s
    pks, _ = find_peaks(sdrop, prominence=14.0, distance=35)
    ts = [float(p + int(rt_row) + 5) + 0.5 for p in pks]
    matched = C.match_ticks_regular(ts, len(tick_vals), tol=7.0)
    ya = C.verify_axis("linear", list(zip(matched, tick_vals)), max_residual_px=3.5)
    return ya, cl[0], cr[0], rt_row, rb_row


def read_sfig14_all_extrema(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read global V_tot max ('max') or min ('min') from all three Supp Fig. 14 colorbars."""
    H, W = image.shape[:2]
    specs = [
        (40, [2000.0, 1800.0, 1600.0, 1400.0, 1200.0, 1000.0], (int(0.24 * W), int(0.35 * W))),
        (80, [1800.0, 1600.0, 1400.0, 1200.0, 1000.0, 800.0], (int(0.56 * W), int(0.68 * W))),
        (120, [1400.0, 1200.0, 1000.0, 800.0], (int(0.88 * W), W)),
    ]
    out = []
    for cp_nM, tvals, xb in specs:
        ya, cl, cr, rt_row, rb_row = _cal_sfig14_colorbar(image, tvals, xb)
        row_pt = rt_row if which == "max" else rb_row
        out.append(
            C.make_point_reading(
                float(cp_nM),
                0.5 * (cl + cr),
                row_pt,
                ya,
                half_px=2.5 + ya.residual_px,
                note=f"SuppFig14 {cp_nM}nM colorbar {which} (total VCA), cal_res={ya.residual_px:.2f}px",
            )
        )
    return out


def read_fig5f_80nM_t40_profile(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the t = 40 s double-peaked radial profile B(r, t=40 s) on the front plane of Fig. 5f (80 nM CP).

    Calibrates the orthographic 3D front plane (t = 40 s) from the left z-axis ticks (0..400),
    the right vertical box edge (r = +5), and the bottom r-axis and t-axis edges.
    """
    H, W = image.shape[:2]
    gray = image.sum(axis=2, dtype=np.uint16)
    sat = image.max(axis=2) - image.min(axis=2)
    ndark = (gray < 300) & (sat < 20)

    # 1. Left vertical z-axis spine and its 5 ticks (0..400)
    v_cnt = ndark[:1300, 650:830].sum(axis=0)
    c_z = int(np.argmax(v_cnt)) + 650
    z_rows = np.where(ndark[:1300, c_z - 2 : c_z + 3].any(axis=1))[0]
    r_ztop, r_zbot = int(z_rows[0]), int(z_rows[-1])
    t_prof = ndark[r_ztop : r_zbot + 5, c_z - 20 : c_z - 5].mean(axis=1)
    yt = [c[0] + r_ztop for c in C.clusters(np.where(t_prof > 0.3)[0], gap=4) if (c[2] - c[1] + 1) < 25]
    yt_m = C.match_ticks_regular(yt + [float(r_zbot) + 0.5], 5, tol=8.0)
    z_ax = C.verify_axis("linear", list(zip(yt_m, [400.0, 300.0, 200.0, 100.0, 0.0])), max_residual_px=3.0)
    dz = 1.0 / z_ax.b

    # 2. Fit the clean bottom-front r-axis line across open columns c = 1650..2130 (no floor occlusion above)
    pts_r = []
    for c in range(1650, 2130, 3):
        rs = C.runs_1d(ndark[1220:1410, c], gap=1)
        for a, b_ in rs:
            y_mid = 0.5 * (a + b_) + 1220.0
            if 4 <= (b_ - a + 1) <= 10 and gray[int(y_mid) - 18, c] > 650:
                pts_r.append((float(c), y_mid))
    pts_r = np.array(pts_r)
    mr, br = np.polyfit(pts_r[:, 0], pts_r[:, 1], 1)
    for _ in range(3):
        inl = np.abs(pts_r[:, 1] - (mr * pts_r[:, 0] + br)) < 2.0
        mr, br = np.polyfit(pts_r[inl, 0], pts_r[inl, 1], 1)

    # 3. Fit bottom-left t-axis line to get the outward tick ray direction (+1, mt)
    pts_t = []
    for c in range(c_z + 30, 1400, 5):
        r_est = int(z_ax.to_pixel(0.0) + 0.56 * (c - c_z))
        ys = np.where(ndark[max(0, r_est - 80) : min(H, r_est + 80), c])[0]
        if len(ys):
            pts_t.append((float(c), float(ys[-1]) + max(0, r_est - 80)))
    pts_t = np.array(pts_t)
    mt, bt = np.polyfit(pts_t[:, 0], pts_t[:, 1], 1)
    for _ in range(3):
        inl_t = np.abs(pts_t[:, 1] - (mt * pts_t[:, 0] + bt)) < 2.5
        mt, bt = np.polyfit(pts_t[inl_t, 0], pts_t[inl_t, 1], 1)

    # 4. Detect all 5 outward ticks along the r-axis (-5, -2.5, 0, +2.5, +5) via the +t ray
    cs = np.arange(1400, 2320)
    s_vals = np.arange(8, 22)
    sc = np.array(
        [ndark[np.round(mr * c + br + s_vals * mt).astype(int), np.round(c + s_vals).astype(int)].mean() for c in cs]
    )
    t_cl = C.clusters(cs[sc > 0.5], gap=5)
    t_cs = [x[0] for x in t_cl if (x[2] - x[1] + 1) <= 20]
    matched = C.match_ticks_regular(t_cs, 5, tol=4.0)
    r_ax = C.verify_axis("linear", list(zip(matched, [-5.0, -2.5, 0.0, 2.5, 5.0])), max_residual_px=2.0)

    # Colored mesh surface (red/orange/yellow: R > B + 40)
    surf = (image[:, :, 0].astype(np.int16) - image[:, :, 2].astype(np.int16)) > 40

    out = []
    for rv in at:
        c_f = r_ax.to_pixel(float(rv))
        r_base = mr * c_f + br
        c_i = int(round(c_f))
        r_top_lim = max(0, int(r_base + 300.0 * dz))
        r_bot_lim = min(H - 1, int(r_base) - 4)
        edges = []
        for dc in range(-4, 5):
            ys = np.where(surf[r_top_lim:r_bot_lim, c_i + dc])[0]
            if len(ys):
                edges.append(float(ys[-1] + r_top_lim) + 0.5)
        if not edges:
            out.append({"x": float(rv), "y": None, "found": False, "note": f"no colored mesh at r={rv}"})
            continue
        r_edge = float(np.median(edges))
        z_val = float((r_edge - r_base) / dz)
        half_px = 6.0 + z_ax.residual_px + r_ax.residual_px
        z_lo = float((r_edge + half_px - r_base) / dz)
        z_hi = float((r_edge - half_px - r_base) / dz)
        out.append(
            {
                "x": float(rv),
                "y": z_val,
                "lo": min(z_lo, z_hi),
                "hi": max(z_lo, z_hi),
                "px": [float(c_f), r_edge],
                "px_lo": r_edge - half_px,
                "px_hi": r_edge + half_px,
                "note": f"Fig5f 80nM t=40s front-plane 3D ortho, z_res={z_ax.residual_px:.2f}px, r_res={r_ax.residual_px:.2f}px",
            }
        )
    return out
