"""Shared calibrated readers for all 12 subpanels (28 curves) of
Figure 4—figure supplement 1 (PDF p. 40, 2400 dpi vector graphics):
  Rows a (kappa=-0.2), b (kappa=-0.5), c (kappa=-0.8) on x in [0, 50];
  Row d (Reduce gamma by 200-fold, chaotic snapshot) on x in [0, 12].
  Col 1: rho/rho_0 (solid black)
  Col 2: sigma_par (solid black) and sigma_perp (dashed black)
  Col 3: active (black solid/dashed) and viscous/hydrodynamic (blue solid/dashed) components.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def _calibrate(image: np.ndarray, row_id: str, col_id: int):
    frame = C.find_box(image, dark=160, max_sat=80, min_frac=0.50)
    xt_c = C.detect_x_ticks(image, frame, side="in", spine="bottom", dark=160, max_sat=80, min_len=25, max_len=150)
    yt_r = C.detect_y_ticks(image, frame, side="in", spine="right", dark=160, max_sat=80, min_len=25, max_len=150)

    # Filter inward ticks on right spine strictly by tick geometry (length ~42 px, width ~7 px)
    yt_inner = sorted(p for p, l, w in yt_r if 32 <= l <= 55 and w <= 11)

    if row_id in ("a", "b", "c"):
        xt = C.select_regular_ticks(xt_c, expected_n=6, tol_px=10.0)
        xa = C.calibrate_axis(xt, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0], max_residual_px=2.5)
        if col_id == 1:
            pts = [frame["r_top"]] + yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [3.5, 3.0, 2.0, 1.0, 0.0], max_residual_px=2.5)
        elif col_id == 2:
            pts = [frame["r_top"]] + yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [2.0, 1.0, 0.0], max_residual_px=2.5)
        elif col_id == 3:
            pts = [frame["r_top"]] + yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [5.0, 4.0, 2.0, -2.0, -3.0], max_residual_px=2.5)
        else:
            raise ValueError(col_id)
    elif row_id == "d":
        xt = C.select_regular_ticks(xt_c, expected_n=7, tol_px=10.0)
        xa = C.calibrate_axis(xt, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0], max_residual_px=2.5)
        if col_id == 1:
            pts = [frame["r_top"]] + yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [2.0, 1.5, 1.0, 0.5], max_residual_px=2.5)
        elif col_id == 2:
            pts = yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [1.2, 1.0, 0.8], max_residual_px=2.5)
        elif col_id == 3:
            pts = [frame["r_top"]] + yt_inner + [frame["r_bot"]]
            ya = C.calibrate_axis(pts, [2.0, 1.0, -1.0], max_residual_px=2.5)
        else:
            raise ValueError(col_id)
    else:
        raise ValueError(row_id)
    return frame, xa, ya


def _interior_mask(raw: np.ndarray, frame: dict, pad_y: int = 6, pad_x: int = 8) -> np.ndarray:
    """Restrict mask to plot interior and erase inward tick stubs along all 4 borders."""
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    m = np.zeros_like(raw)
    m[r0 + pad_y : r1 - pad_y, c0 + pad_x : c1 - pad_x] = raw[r0 + pad_y : r1 - pad_y, c0 + pad_x : c1 - pad_x]
    tick_band = 55
    for band_r0, band_r1 in ((r0 + pad_y, r0 + tick_band), (r1 - tick_band, r1 - pad_y)):
        sub = m[band_r0:band_r1, c0 + pad_x : c1 - pad_x]
        m[band_r0:band_r1, c0 + pad_x : c1 - pad_x] = ndimage.binary_opening(sub, structure=np.ones((1, 13), bool))
    for band_c0, band_c1 in ((c0 + pad_x, c0 + tick_band), (c1 - tick_band, c1 - pad_x)):
        sub = m[r0 + pad_y : r1 - pad_y, band_c0:band_c1]
        m[r0 + pad_y : r1 - pad_y, band_c0:band_c1] = ndimage.binary_opening(sub, structure=np.ones((13, 1), bool))
    return m


def _black_mask(image: np.ndarray, frame: dict) -> np.ndarray:
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    raw = (R < 85) & (G < 85) & (B < 85)
    return _interior_mask(raw, frame, pad_y=6, pad_x=8)


def _blue_mask(image: np.ndarray, frame: dict) -> np.ndarray:
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    raw = (B > 110) & ((B - R) > 45) & ((B - G) > 25)
    return _interior_mask(raw, frame, pad_y=6, pad_x=8)


def _read_two_curve_at(
    mask: np.ndarray,
    xa,
    ya,
    at: list[float],
    frame: dict,
    which_list: list[str],
    solid_is_upper_list: list[bool],
    note_prefix: str = "",
    win: int = 48,
) -> list[dict]:
    """Fast local-window reader for upper/lower branch of a two-curve (solid + dashed) same-color mask.
    Merges vertically overlapping dash fragments on steep slopes so the dashed curve center is
    unbiased, and interpolates across dash gaps on flat troughs/peaks."""
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    r_lo, r_hi = r0 + 4, r1 - 4
    cal_x = max(xa.residual_px, 0.5)
    cal_y = max(ya.residual_px, 0.5)

    # Estimate typical single-stroke half-height from flat/single-run columns across the box
    sample_h = []
    for c in range(c0 + 20, c1 - 20, 25):
        rr = [t for t in C.col_runs(mask, c, r_lo, r_hi) if t[1] >= 3.0]
        if len(rr) == 1:
            sample_h.append(rr[0][1])
    typ_half = float(np.percentile(sample_h, 35)) if sample_h else 13.0

    # Cache column runs only within windows needed for `at`
    cols_runs: dict[int, list[tuple[float, float]]] = {}

    def get_runs(c: int) -> list[tuple[float, float]]:
        if c not in cols_runs:
            if c <= c0 + 4 or c >= c1 - 4:
                cols_runs[c] = []
            else:
                cols_runs[c] = sorted(
                    [t for t in C.col_runs(mask, c, r_lo, r_hi) if t[1] >= 3.0],
                    key=lambda t: t[0],
                )
        return cols_runs[c]

    out = []
    for xv, which, solid_is_upper in zip(at, which_list, solid_is_upper_list):
        pc = xa.to_pixel(xv)
        col = int(round(pc - 0.5))

        sep_cols, sep_upper, sep_lower, sep_hu, sep_hl = [], [], [], [], []
        for c in range(max(c0 + 8, col - win), min(c1 - 8, col + win + 1)):
            rr = get_runs(c)
            if len(rr) >= 2:
                if solid_is_upper:
                    u_c, u_h = rr[0]
                    d_runs = rr[1:]
                    wsum = sum(h for _, h in d_runs)
                    l_c = sum(r * h for r, h in d_runs) / wsum
                    l_h = 0.5 * (max(r + h for r, h in d_runs) - min(r - h for r, h in d_runs))
                else:
                    l_c, l_h = rr[-1]
                    d_runs = rr[:-1]
                    wsum = sum(h for _, h in d_runs)
                    u_c = sum(r * h for r, h in d_runs) / wsum
                    u_h = 0.5 * (max(r + h for r, h in d_runs) - min(r - h for r, h in d_runs))
                sep_cols.append(c + 0.5)
                sep_upper.append(u_c)
                sep_lower.append(l_c)
                sep_hu.append(u_h)
                sep_hl.append(l_h)

        separated = False
        if len(sep_cols) >= 15:
            dx_arr = np.array(sep_cols) - pc
            gaps = np.array(sep_lower) - np.array(sep_upper)
            # Require separated columns to span near/across pc (avoid one-sided flank extrapolation)
            if float(np.median(gaps)) > 1.25 * typ_half and dx_arr.min() <= 16.0 and dx_arr.max() >= -16.0:
                separated = True

        if separated:
            dx = np.array(sep_cols) - pc
            yr = np.array(sep_upper if which == "upper" else sep_lower)
            halfs = np.array(sep_hu if which == "upper" else sep_hl)
            w = np.exp(-0.5 * (dx / 22.0) ** 2)
            # Use degree 2 only when well-bracketed on BOTH sides of pc; otherwise linear
            bracketed = (dx.min() <= -10.0) and (dx.max() >= 10.0)
            deg = 2 if (bracketed and (dx.max() - dx.min()) > 28.0 and len(dx) >= 15) else 1
            coeff = np.polyfit(dx, yr, deg, w=np.sqrt(w))
            center = float(coeff[-1])
            slope = float(coeff[-2]) if deg >= 1 else 0.0
            resid = yr - np.polyval(coeff, dx)
            wrms = float(math.sqrt(np.sum(w * resid * resid) / max(np.sum(w), 1e-12)))
            h_px = max(float(np.median(halfs)), 6.0) + max(wrms, 1.0) + abs(slope) * cal_x + cal_y
            note = f"{note_prefix}separated {which} n={len(dx)} deg={deg} wrms={wrms:.1f}px"
        else:
            rr_here = get_runs(col)
            if not rr_here:
                for dc in (1, -1, 2, -2, 4, -4, 8, -8):
                    rr_here = get_runs(col + dc)
                    if rr_here:
                        col = col + dc
                        break
            if not rr_here:
                out.append({"x": float(xv), "y": None, "found": False, "note": f"{note_prefix}no ink at x={xv:g}"})
                continue
            top_e = min(r - h for r, h in rr_here)
            bot_e = max(r + h for r, h in rr_here)
            center_raw = 0.5 * (top_e + bot_e)
            lc = get_runs(col - 6)
            rc = get_runs(col + 6)
            slope = 0.0
            if lc and rc:
                l_r = min(lc, key=lambda t: abs(t[0] - center_raw))[0]
                r_r = min(rc, key=lambda t: abs(t[0] - center_raw))[0]
                slope = (r_r - l_r) / 12.0

            # On flat regions (troughs/peaks), scan nearby columns to bridge white dash gaps
            if abs(slope) < 0.30:
                win_tops, win_bots, win_fh = [], [], []
                for c in range(max(c0 + 8, col - 70), min(c1 - 8, col + 71), 2):
                    rr = get_runs(c)
                    if rr:
                        te = min(r - h for r, h in rr)
                        be = max(r + h for r, h in rr)
                        win_tops.append(te)
                        win_bots.append(be)
                        win_fh.append(0.5 * (be - te))
                if win_fh:
                    max_fh = float(np.percentile(win_fh, 95))
                    if max_fh > 1.25 * typ_half:
                        sel = [i for i, fh in enumerate(win_fh) if fh >= 0.85 * max_fh]
                        top_e = float(np.median([win_tops[i] for i in sel]))
                        bot_e = float(np.median([win_bots[i] for i in sel]))

            full_h = 0.5 * (bot_e - top_e)
            stroke_h = typ_half * math.hypot(1.0, slope)
            if full_h > 1.22 * stroke_h:
                center = (top_e + stroke_h) if which == "upper" else (bot_e - stroke_h)
                half = max(stroke_h, full_h - stroke_h + 2.0)
                note = f"{note_prefix}touching {which} edge (full_h={full_h:.1f}px)"
            else:
                center = 0.5 * (top_e + bot_e)
                half = full_h
                note = f"{note_prefix}coincident curves (half={half:.1f}px)"
            h_px = half + 1.0 + abs(slope) * (abs(pc - (col + 0.5)) + cal_x) + cal_y

        yv = ya.to_value(center)
        lo, hi = sorted((ya.to_value(center - h_px), ya.to_value(center + h_px)))
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(pc), float(center)],
            "px_lo": float(center - h_px),
            "px_hi": float(center + h_px),
            "note": f"{note}; cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


def read_col1_rho(image: np.ndarray, at: list[float], row_id: str) -> list[dict]:
    """Read density profile rho/rho_0 (single solid black curve in col 1)."""
    frame, xa, ya = _calibrate(image, row_id, 1)
    mask = _black_mask(image, frame)
    return C.read_curve_at(
        mask, xa, ya, at, frame["box"], native_px=1.0, search_cols=4, pick="thickest",
        note_prefix=f"fig4s1_{row_id} rho: ",
    )


def read_col2_tension(image: np.ndarray, at: list[float], row_id: str, comp: str) -> list[dict]:
    """Read total tension sigma_par (comp='par', solid) or sigma_perp (comp='perp', dashed) in col 2."""
    frame, xa, ya = _calibrate(image, row_id, 2)
    mask = _black_mask(image, frame)
    c0_b, r0_b, c1_b, r1_b = [int(v) for v in frame["box"]]

    solid_upper_list: list[bool] = []
    if row_id == "c":
        solid_upper_list = [False] * len(at)
    elif row_id in ("a", "b"):
        y_thresh = 0.60 if row_id == "a" else 0.65
        for xv in at:
            col = int(round(xa.to_pixel(xv) - 0.5))
            rr = C.col_runs(mask, col, r0_b + 4, r1_b - 4)
            y_est = ya.to_value(np.mean([t[0] for t in rr])) if rr else 1.0
            solid_upper_list.append(y_est > y_thresh)
    elif row_id == "d":
        for xv in at:
            col = int(round(xa.to_pixel(xv) - 0.5))
            upper_ys, lower_ys = [], []
            for c in range(max(c0_b + 8, col - 70), min(c1_b - 8, col + 71), 3):
                rr = [t for t in C.col_runs(mask, c, r0_b + 4, r1_b - 4) if t[1] >= 3.0]
                if len(rr) >= 2:
                    upper_ys.append(ya.to_value(min(rr, key=lambda t: t[0])[0]))
                    lower_ys.append(ya.to_value(max(rr, key=lambda t: t[0])[0]))
            if upper_ys:
                u_med, l_med = float(np.median(upper_ys)), float(np.median(lower_ys))
                u_dist = max(0.0, u_med - 1.01) + max(0.0, 0.925 - u_med)
                l_dist = max(0.0, l_med - 1.01) + max(0.0, 0.925 - l_med)
                solid_upper_list.append(u_dist > l_dist)
            else:
                solid_upper_list.append(True)
    else:
        raise ValueError(row_id)

    which_list = [
        ("upper" if su else "lower") if comp == "par" else ("lower" if su else "upper")
        for su in solid_upper_list
    ]
    return _read_two_curve_at(
        mask, xa, ya, at, frame, which_list, solid_upper_list, note_prefix=f"fig4s1_{row_id} {comp}: "
    )


def read_col3_active(image: np.ndarray, at: list[float], row_id: str, comp: str) -> list[dict]:
    """Read active stress component in col 3: black solid sigma_act_par (ALWAYS lower)
    or black dashed sigma_act_perp (ALWAYS upper)."""
    frame, xa, ya = _calibrate(image, row_id, 3)
    mask = _black_mask(image, frame)
    solid_upper_list = [False] * len(at)
    which_list = ["lower" if comp == "par" else "upper"] * len(at)
    return _read_two_curve_at(
        mask, xa, ya, at, frame, which_list, solid_upper_list, note_prefix=f"fig4s1_{row_id} act_{comp}: "
    )


def read_col3_viscous(image: np.ndarray, at: list[float], row_id: str, comp: str) -> list[dict]:
    """Read viscous/hydrodynamic stress component in col 3: blue solid sigma_hydr_par
    or blue dashed sigma_hydr_perp.
    For y < 0 (converging bundles): solid is UPPER, dashed is LOWER.
    For y > 0 (diverging troughs): dashed is UPPER, solid is LOWER.
    """
    frame, xa, ya = _calibrate(image, row_id, 3)
    mask = _blue_mask(image, frame)
    r_zero = ya.to_pixel(0.0)
    c0_b, r0_b, c1_b, r1_b = [int(v) for v in frame["box"]]
    solid_upper_list: list[bool] = []
    for xv in at:
        col = int(round(xa.to_pixel(xv) - 0.5))
        b_rows = []
        for c in range(max(c0_b + 8, col - 50), min(c1_b - 8, col + 51), 3):
            for r_c, h in C.col_runs(mask, c, r0_b + 4, r1_b - 4):
                if h >= 3.0:
                    b_rows.append(r_c)
        if b_rows and float(np.median(b_rows)) > r_zero + 8.0:
            solid_upper_list.append(True)
        elif b_rows and float(np.median(b_rows)) < r_zero - 8.0:
            solid_upper_list.append(False)
        else:
            solid_upper_list.append(True)
    which_list = [
        ("upper" if su else "lower") if comp == "par" else ("lower" if su else "upper")
        for su in solid_upper_list
    ]
    return _read_two_curve_at(
        mask, xa, ya, at, frame, which_list, solid_upper_list, note_prefix=f"fig4s1_{row_id} vis_{comp}: "
    )
