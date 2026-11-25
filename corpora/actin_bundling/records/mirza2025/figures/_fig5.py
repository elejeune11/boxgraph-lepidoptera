"""Shared calibrated readers for Fig. 5(c), (d), (e), (f), and (g) (PDF p. 41)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 5(c): solid horizontal time-average lines for S_0 = 0 and S_0 = 0.3
# --------------------------------------------------------------------------

def read_fig5c(image: np.ndarray, at: list[float], subpanel: str, component: str) -> list[dict]:
    """subpanel: 'S0_0' (top) or 'S0_0p3' (bottom); component: 'par' (blue) or 'perp' (red)."""
    H, W = image.shape[:2]
    mid_r = H // 2 - 200  # split between top box (r~120..1415) and bottom box (r~1552..2839)
    # Use bottom box for x ticks (labeled 2, 3, 4)
    bot_frame = C.find_box(image, dark=150, max_sat=50, min_frac=0.40, row_range=(mid_r, H))
    xt_c = C.detect_x_ticks(image, bot_frame, side="out", spine="bottom", dark=150, max_sat=50, min_len=8, max_len=80)
    xt = C.select_regular_ticks(xt_c, expected_n=3, tol_px=15.0)
    xa = C.calibrate_axis(xt, [2.0, 3.0, 4.0], max_residual_px=4.0)

    if subpanel == "S0_0":
        frame = C.find_box(image, dark=150, max_sat=50, min_frac=0.40, row_range=(0, mid_r))
    else:
        frame = bot_frame

    # Outward left y-ticks: 3 ticks corresponding to 2 (top), 1 (middle), 0 (bottom)
    dn = C.neutral_dark_mask(image, dark=150, max_sat=50)
    r0, r1 = frame["top"][0] - 10, frame["bot"][1] + 10
    cl = frame["left"][0]
    sl = dn[max(0, r0) : min(H, r1), max(0, cl - 40) : cl]
    has = sl[:, -25:].mean(axis=1) > 0.5
    idx = np.where(has)[0]
    yt = [max(0, r0) + (a + b) / 2.0 + 0.5 for a, b in C.groups(idx, gap=3)]
    ya = C.calibrate_axis(yt, [2.0, 1.0, 0.0], max_residual_px=4.0)

    # Locate the solid horizontal time-average line across the full width of the subpanel
    c0, tr0, c1, tr1 = [int(v) for v in frame["box"]]
    patch = image[tr0 + 10 : tr1 - 10, c0 + 20 : c1 - 20]
    R, G, B = patch[:, :, 0].astype(int), patch[:, :, 1].astype(int), patch[:, :, 2].astype(int)
    if component == "par":
        mask = (B > 180) & (R < 90) & (G < 90)
    else:
        mask = (R > 180) & (G < 90) & (B < 90)

    row_frac = mask.mean(axis=1)
    line_rows = np.where(row_frac > 0.70)[0]
    if len(line_rows) == 0:
        raise RuntimeError(f"Fig 5(c) {subpanel} {component}: solid horizontal line not found")
    grps = C.groups(line_rows, gap=2)
    a_r, b_r = max(grps, key=lambda g: g[1] - g[0])
    center_r = (tr0 + 10) + (a_r + b_r) / 2.0 + 0.5
    half_r = (b_r - a_r + 1) / 2.0
    # Native 300 dpi -> 8 px per native pixel at 2400 dpi
    native_px = 8.0
    h_px = max(half_r, 0.5 * native_px) + 0.5 * native_px + max(ya.residual_px, 0.5 * native_px)

    yv = ya.to_value(center_r)
    lo, hi = ya.to_value(center_r - h_px), ya.to_value(center_r + h_px)
    lo, hi = min(lo, hi), max(lo, hi)

    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(pc), float(center_r)],
            "px_lo": float(center_r - h_px),
            "px_hi": float(center_r + h_px),
            "note": f"solid horizontal avg line rows=[{tr0+10+a_r},{tr0+10+b_r}] cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


# --------------------------------------------------------------------------
# Fig. 5(d): normalized mean tension vs normalized density
# --------------------------------------------------------------------------

def _calibrate_fig5d(image: np.ndarray):
    frame = C.find_box(image, dark=170, max_sat=80, min_frac=0.45)
    xt_c = C.detect_x_ticks(image, frame, side="out", spine="bottom", dark=170, max_sat=80, min_len=8, max_len=80)
    yt_c = C.detect_y_ticks(image, frame, side="out", spine="left", dark=170, max_sat=80, min_len=8, max_len=80)
    xt = C.select_regular_ticks(xt_c, expected_n=3, tol_px=12.0)
    yt = C.select_regular_ticks(yt_c, expected_n=5, tol_px=12.0)
    xa = C.calibrate_axis(xt, [1.0, 2.0, 3.0], max_residual_px=3.0)
    ya = C.calibrate_axis(yt, [4.0, 3.0, 2.0, 1.0, 0.0], max_residual_px=3.0)
    return frame, xa, ya


def read_fig5d_series(image: np.ndarray, series: str) -> list[dict]:
    """series: 'S0_0p1' (blue) or 'S0_0p3' (red)."""
    frame, xa, ya = _calibrate_fig5d(image)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    if series == "S0_0p3":
        cmask = (R > 160) & (G < 90) & (B < 90)
    elif series == "S0_0p1":
        cmask = (B > 160) & (R < 90) & (G < 90)
    else:
        raise ValueError(f"Unknown Fig 5(d) series: {series}")

    # Mask out legend region (bottom right)
    leg_c0 = int(c0 + 0.55 * (c1 - c0))
    leg_r0 = int(r0 + 0.45 * (r1 - r0))
    cmask[leg_r0:, leg_c0:] = False

    native_px = 6.0  # 1800 dpi / 300 dpi
    cal_y = max(ya.residual_px, 0.5 * native_px)
    out = []
    for rho_val in (1.0, 2.0, 3.0):
        xc = xa.to_pixel(rho_val)
        # Search vertical window around y = rho_val +- 0.5
        r_search_lo = max(r0 + 5, int(ya.to_pixel(rho_val + 0.5)))
        r_search_hi = min(r1 - 5, int(ya.to_pixel(rho_val - 0.5)))
        caps = C.find_errorbar_caps_at_x(
            cmask, xc, r_search_lo, r_search_hi,
            flank_inner=28, flank_outer=72, min_flank_px=12, max_cap_h=14,
        )
        if len(caps) < 2:
            raise RuntimeError(f"Fig 5(d) {series} rho={rho_val}: expected >=2 caps, found {caps}")
        top_cap, bot_cap = caps[0], caps[-1]
        r_top_cap, h_top = top_cap[0], top_cap[1]
        r_bot_cap, h_bot = bot_cap[0], bot_cap[1]

        # For red (top layer), also measure circular marker core centroid
        r_mid_caps = 0.5 * (r_top_cap + r_bot_cap)
        if series == "S0_0p3":
            sub = cmask[int(r_top_cap + 6) : int(r_bot_cap - 6), int(xc - 36) : int(xc + 36)]
            core = ndimage.uniform_filter(sub.astype(np.float32), size=19) > 0.75
            rr, cc = np.nonzero(core)
            if len(rr) > 20:
                r_dot = int(r_top_cap + 6) + float(rr.mean() + 0.5)
                c_dot = int(xc - 36) + float(cc.mean() + 0.5)
                r_center = r_dot
                c_center = c_dot
                note = f"marker+caps (cap_mid_diff={r_dot - r_mid_caps:+.1f}px)"
            else:
                r_center, c_center = r_mid_caps, xc
                note = "symmetric cap midpoint"
        else:
            r_center, c_center = r_mid_caps, xc
            note = "symmetric cap midpoint (marker partially behind S0=0.3)"

        h_center = 0.5 * native_px + cal_y + 2.0
        h_cap = max(h_top, h_bot, 0.5 * native_px) + 0.5 * native_px + cal_y

        yv = ya.to_value(r_center)
        lo_y, hi_y = sorted((ya.to_value(r_center - h_center), ya.to_value(r_center + h_center)))
        bar_lo, bar_hi = sorted((ya.to_value(r_top_cap), ya.to_value(r_bot_cap)))
        bar_lo_out, bar_hi_out = sorted((ya.to_value(r_top_cap - h_cap), ya.to_value(r_bot_cap + h_cap)))

        out.append({
            "x": float(xa.to_value(c_center)),
            "y": float(yv),
            "lo": float(lo_y),
            "hi": float(hi_y),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bar_lo_out),
            "bar_hi_outer": float(bar_hi_out),
            "px": [float(c_center), float(r_center)],
            "px_lo": float(r_top_cap),
            "px_hi": float(r_bot_cap),
            "note": f"{note}; cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


def read_fig5d_theory_line(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the dashed theoretical line lambda = 1.0 in Fig. 5(d)."""
    frame, xa, ya = _calibrate_fig5d(image)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    blk = C.neutral_dark_mask(image, dark=90, max_sat=40)
    mask = np.zeros_like(blk)
    mask[r0 + 8 : r1 - 8, c0 + 8 : c1 - 8] = blk[r0 + 8 : r1 - 8, c0 + 8 : c1 - 8]
    # Exclude top-left legend ("--- lambda = 1.0") and bottom-right legend
    mask[r0 : int(r0 + 0.28 * (r1 - r0)), c0 : int(c0 + 0.65 * (c1 - c0))] = False
    mask[int(r0 + 0.45 * (r1 - r0)) :, int(c0 + 0.50 * (c1 - c0)) :] = False
    # Exclude marker columns (rho = 1, 2, 3 +- 75 px) so only clean black dashes remain
    for rv in (1.0, 2.0, 3.0):
        xc = int(round(xa.to_pixel(rv)))
        mask[:, max(0, xc - 75) : min(mask.shape[1], xc + 75)] = False

    # Fit straight line to all remaining dash pixels
    rr, cc = np.nonzero(mask)
    if len(rr) < 50:
        raise RuntimeError(f"Fig 5(d) theory line: too few pixels ({len(rr)})")
    # Linear fit r = m * c + b (with one outlier rejection pass)
    coeff = np.polyfit(cc + 0.5, rr + 0.5, 1)
    resid = (rr + 0.5) - np.polyval(coeff, cc + 0.5)
    keep = np.abs(resid) < 25.0
    coeff = np.polyfit(cc[keep] + 0.5, rr[keep] + 0.5, 1)
    resid = (rr[keep] + 0.5) - np.polyval(coeff, cc[keep] + 0.5)
    rms = float(np.sqrt(np.mean(resid ** 2)))
    slope = float(coeff[0])
    native_px = 6.0
    cal_x = max(xa.residual_px, 0.5 * native_px)
    cal_y = max(ya.residual_px, 0.5 * native_px)
    h_px = max(rms, 0.5 * native_px) + 0.5 * native_px + abs(slope) * cal_x + cal_y

    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        pr = float(np.polyval(coeff, pc))
        yv = ya.to_value(pr)
        lo, hi = sorted((ya.to_value(pr - h_px), ya.to_value(pr + h_px)))
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(pc), pr],
            "px_lo": pr - h_px,
            "px_hi": pr + h_px,
            "note": f"linear fit to {int(keep.sum())} dash px, rms={rms:.1f}px",
        })
    return out


# --------------------------------------------------------------------------
# Fig. 5(e): deviatoric tension ratio (sigma_par - sigma_perp)/sigma_bar vs S_0
# --------------------------------------------------------------------------

def _calibrate_fig5e(image: np.ndarray):
    frame = C.find_box(image, dark=170, max_sat=80, min_frac=0.45)
    xt_c = C.detect_x_ticks(image, frame, side="out", spine="bottom", dark=170, max_sat=80, min_len=8, max_len=80)
    yt_c = C.detect_y_ticks(image, frame, side="out", spine="left", dark=170, max_sat=80, min_len=8, max_len=80)
    xt = C.select_regular_ticks(xt_c, expected_n=4, tol_px=12.0)
    yt = C.select_regular_ticks(yt_c, expected_n=4, tol_px=12.0)
    xa = C.calibrate_axis(xt, [0.0, 0.1, 0.2, 0.3], max_residual_px=3.0)
    ya = C.calibrate_axis(yt, [0.2, 0.0, -0.2, -0.4], max_residual_px=3.0)
    return frame, xa, ya


def read_fig5e_series(image: np.ndarray, series: str) -> list[dict]:
    """series: 'rho_1x' (blue), 'rho_2x' (green), 'rho_3x' (red)."""
    frame, xa, ya = _calibrate_fig5e(image)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    if series == "rho_1x":
        cmask = (B > 120) & ((B - R) > 65) & ((B - G) > 15)
    elif series == "rho_2x":
        cmask = (G > 65) & ((G - R) > 30) & ((G - B) > 15)
    elif series == "rho_3x":
        cmask = (R > 140) & ((R - G) > 65) & ((R - B) > 55)
    else:
        raise ValueError(f"Unknown Fig 5(e) series: {series}")

    # Mask out top-right legend (columns > 0.62*W, rows < 0.45*H)
    cmask[: int(r0 + 0.45 * (r1 - r0)), int(c0 + 0.62 * (c1 - c0)) :] = False

    native_px = 6.0
    cal_y = max(ya.residual_px, 0.5 * native_px)
    out = []
    for s0_val in (0.0, 0.1, 0.2, 0.3):
        xc = xa.to_pixel(s0_val)
        caps = C.find_errorbar_caps_at_x(
            cmask, xc, r0 + 5, r1 - 5,
            flank_inner=30, flank_outer=72, min_flank_px=12, max_cap_h=22,
        )
        if len(caps) < 2:
            raise RuntimeError(f"Fig 5(e) {series} S0={s0_val}: expected >=2 caps, found {caps}")
        top_cap, bot_cap = caps[0], caps[-1]
        r_top_cap, h_top = top_cap[0], top_cap[1]
        r_bot_cap, h_bot = bot_cap[0], bot_cap[1]
        r_mid = 0.5 * (r_top_cap + r_bot_cap)

        # Cross-check against unobstructed circular marker core near r_mid
        band_lo = max(r0 + 5, int(r_mid - 45))
        band_hi = min(r1 - 5, int(r_mid + 45))
        sub = cmask[band_lo:band_hi, int(xc - 36) : int(xc + 36)]
        core = ndimage.uniform_filter(sub.astype(np.float32), size=21) > 0.78
        rr, cc = np.nonzero(core)
        if len(rr) > 120:
            r_dot = band_lo + float(rr.mean() + 0.5)
            c_dot = int(xc - 36) + float(cc.mean() + 0.5)
            if abs(r_dot - r_mid) <= 12.0:
                r_center, c_center = r_dot, c_dot
                note = f"marker+caps (cap_mid_diff={r_dot - r_mid:+.1f}px)"
            else:
                r_center, c_center = r_mid, xc
                note = f"symmetric cap midpoint (dot shifted {r_dot - r_mid:+.1f}px by overlap)"
        else:
            r_center, c_center = r_mid, xc
            note = "symmetric cap midpoint (marker partially occluded)"

        h_center = 0.5 * native_px + cal_y + 2.0
        h_cap = max(h_top, h_bot, 0.5 * native_px) + 0.5 * native_px + cal_y
        yv = ya.to_value(r_center)
        lo_y, hi_y = sorted((ya.to_value(r_center - h_center), ya.to_value(r_center + h_center)))
        bar_lo, bar_hi = sorted((ya.to_value(r_top_cap), ya.to_value(r_bot_cap)))
        bar_lo_out, bar_hi_out = sorted((ya.to_value(r_top_cap - h_cap), ya.to_value(r_bot_cap + h_cap)))

        out.append({
            "x": float(xa.to_value(c_center)),
            "y": float(yv),
            "lo": float(lo_y),
            "hi": float(hi_y),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bar_lo_out),
            "bar_hi_outer": float(bar_hi_out),
            "px": [float(c_center), float(r_center)],
            "px_lo": float(r_top_cap),
            "px_hi": float(r_bot_cap),
            "note": f"{note}; cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


def read_fig5e_fit_line(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the black dashed linear regression line (kappa = -1.6) in Fig. 5(e)."""
    frame, xa, ya = _calibrate_fig5e(image)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    blk = C.neutral_dark_mask(image, dark=90, max_sat=40)
    mask = np.zeros_like(blk)
    mask[r0 + 8 : r1 - 8, c0 + 8 : c1 - 8] = blk[r0 + 8 : r1 - 8, c0 + 8 : c1 - 8]
    # Exclude top-right legend and bottom-left legend ("--- kappa = -1.6")
    mask[: int(r0 + 0.45 * (r1 - r0)), int(c0 + 0.62 * (c1 - c0)) :] = False
    mask[int(r1 - 320) :, : int(c0 + 0.78 * (c1 - c0))] = False
    # Exclude narrow bands around the 4 marker x positions
    for sv in (0.0, 0.1, 0.2, 0.3):
        xc = int(round(xa.to_pixel(sv)))
        mask[:, max(0, xc - 60) : min(mask.shape[1], xc + 60)] = False

    rr, cc = np.nonzero(mask)
    if len(rr) < 50:
        raise RuntimeError(f"Fig 5(e) fit line: too few pixels ({len(rr)})")
    coeff = np.polyfit(cc + 0.5, rr + 0.5, 1)
    resid = (rr + 0.5) - np.polyval(coeff, cc + 0.5)
    keep = np.abs(resid) < 25.0
    coeff = np.polyfit(cc[keep] + 0.5, rr[keep] + 0.5, 1)
    resid = (rr[keep] + 0.5) - np.polyval(coeff, cc[keep] + 0.5)
    rms = float(np.sqrt(np.mean(resid ** 2)))
    slope = float(coeff[0])
    native_px = 6.0
    cal_x = max(xa.residual_px, 0.5 * native_px)
    cal_y = max(ya.residual_px, 0.5 * native_px)
    h_px = max(rms, 0.5 * native_px) + 0.5 * native_px + abs(slope) * cal_x + cal_y

    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        pr = float(np.polyval(coeff, pc))
        yv = ya.to_value(pr)
        lo, hi = sorted((ya.to_value(pr - h_px), ya.to_value(pr + h_px)))
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(pc), pr],
            "px_lo": pr - h_px,
            "px_hi": pr + h_px,
            "note": f"linear fit to {int(keep.sum())} dash px, rms={rms:.1f}px",
        })
    return out


# --------------------------------------------------------------------------
# Fig. 5(f): ensemble-averaged S(t) - S_0 vs t / tau_A
# --------------------------------------------------------------------------

def read_fig5f_curve(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """series: 'S0_0p0' (black), 'S0_0p2' (blue), 'S0_0p4' (green)."""
    frame = C.find_box(image, dark=170, max_sat=80, min_frac=0.45)
    xt_c = C.detect_x_ticks(image, frame, side="out", spine="bottom", dark=170, max_sat=80, min_len=8, max_len=80)
    yt_c = C.detect_y_ticks(image, frame, side="out", spine="left", dark=170, max_sat=80, min_len=8, max_len=80)
    xt = C.select_regular_ticks(xt_c, expected_n=4, tol_px=12.0)
    yt = C.select_regular_ticks(yt_c, expected_n=3, tol_px=12.0)
    xa = C.calibrate_axis(xt, [0.0, 2.0, 4.0, 6.0], max_residual_px=3.0)
    ya = C.calibrate_axis(yt, [0.2, 0.1, 0.0], max_residual_px=3.0)

    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    if series == "S0_0p0":
        cmask = (R < 75) & (G < 75) & (B < 75)
    elif series == "S0_0p2":
        cmask = (B > 170) & (R < 75) & (G < 75)
    elif series == "S0_0p4":
        cmask = (G > 75) & (G < 185) & (R < 60) & (B < 60)
    else:
        raise ValueError(f"Unknown Fig 5(f) series: {series}")

    mask = np.zeros_like(cmask)
    mask[r0 + 6 : r1 - 6, c0 + 6 : c1 - 6] = cmask[r0 + 6 : r1 - 6, c0 + 6 : c1 - 6]
    # Exclude top-left legend text
    mask[: int(r0 + 0.48 * (r1 - r0)), : int(c0 + 0.44 * (c1 - c0))] = False

    return C.read_curve_at(
        mask, xa, ya, at, frame["box"], native_px=6.0, search_cols=8, pick="thickest",
        note_prefix=f"fig5f {series}: ",
    )


# --------------------------------------------------------------------------
# Fig. 5(g): turnover-normalized ordering rate S_dot vs S_0
# --------------------------------------------------------------------------

def read_fig5g_series(image: np.ndarray, series: str) -> list[dict]:
    """series: 'kT_pos' (black) or 'kT_zero' (red)."""
    frame = C.find_box(image, dark=170, max_sat=80, min_frac=0.45)
    xt_c = C.detect_x_ticks(image, frame, side="out", spine="bottom", dark=170, max_sat=80, min_len=8, max_len=80)
    yt_c = C.detect_y_ticks(image, frame, side="out", spine="left", dark=170, max_sat=80, min_len=8, max_len=80)
    xt = C.select_regular_ticks(xt_c, expected_n=5, tol_px=12.0)
    yt = C.select_regular_ticks(yt_c, expected_n=5, tol_px=12.0)
    xa = C.calibrate_axis(xt, [0.0, 0.1, 0.2, 0.3, 0.4], max_residual_px=4.0)
    ya = C.calibrate_axis(yt, [0.04, 0.03, 0.02, 0.01, 0.00], max_residual_px=4.0)

    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    if series == "kT_zero":
        cmask = (R > 165) & (G < 85) & (B < 85)
    elif series == "kT_pos":
        cmask = (R < 90) & (G < 90) & (B < 90)
    else:
        raise ValueError(f"Unknown Fig 5(g) series: {series}")

    # Mask out bottom-right legend (columns > 0.52*W, rows > 0.65*H)
    cmask[int(r0 + 0.65 * (r1 - r0)) :, int(c0 + 0.52 * (c1 - c0)) :] = False

    native_px = 6.0
    cal_y = max(ya.residual_px, 0.5 * native_px)
    out = []
    for s0_val in (0.0, 0.1, 0.2, 0.3, 0.4):
        xc = xa.to_pixel(s0_val)
        caps = C.find_errorbar_caps_at_x(
            cmask, xc, r0 + 5, r1 - 5,
            flank_inner=30, flank_outer=72, min_flank_px=12, max_cap_h=22,
        )
        # Note: at S_0 = 0, black bottom cap lies directly underneath the red bottom cap
        if series == "kT_pos" and s0_val == 0.0 and len(caps) == 1:
            r_top_cap, h_top = caps[0][0], caps[0][1]
            r_dot_search = (int(r_top_cap + 60), r1 - 5)
            dark_or_red = (G < 75) & (B < 75)
            caps_any = C.find_errorbar_caps_at_x(
                dark_or_red, xc, r_dot_search[0], r_dot_search[1],
                flank_inner=30, flank_outer=72, min_flank_px=12, max_cap_h=22,
            )
            if caps_any:
                r_bot_cap, h_bot = caps_any[-1][0], caps_any[-1][1]
                cap_note = "black bottom cap co-located under red bottom cap"
            else:
                raise RuntimeError("Fig 5(g) kT_pos S0=0: missing bottom cap")
        else:
            if len(caps) < 2:
                raise RuntimeError(f"Fig 5(g) {series} S0={s0_val}: expected >=2 caps, found {caps}")
            r_top_cap, h_top = caps[0][0], caps[0][1]
            r_bot_cap, h_bot = caps[-1][0], caps[-1][1]
            cap_note = "both caps directly detected"

        r_mid = 0.5 * (r_top_cap + r_bot_cap)
        band_lo = max(r0 + 5, int(r_mid - 45))
        band_hi = min(r1 - 5, int(r_mid + 45))
        sub = cmask[band_lo:band_hi, int(xc - 36) : int(xc + 36)]
        core = ndimage.uniform_filter(sub.astype(np.float32), size=21) > 0.78
        rr, cc = np.nonzero(core)
        if len(rr) > 120:
            r_dot = band_lo + float(rr.mean() + 0.5)
            c_dot = int(xc - 36) + float(cc.mean() + 0.5)
            if abs(r_dot - r_mid) <= 14.0:
                r_center, c_center = r_dot, c_dot
                note = f"marker+caps ({cap_note}; cap_mid_diff={r_dot - r_mid:+.1f}px)"
            else:
                r_center, c_center = r_mid, xc
                note = f"symmetric cap midpoint ({cap_note})"
        else:
            r_center, c_center = r_mid, xc
            note = f"symmetric cap midpoint ({cap_note})"

        h_center = 0.5 * native_px + cal_y + 2.0
        h_cap = max(h_top, h_bot, 0.5 * native_px) + 0.5 * native_px + cal_y
        yv = ya.to_value(r_center)
        lo_y, hi_y = sorted((ya.to_value(r_center - h_center), ya.to_value(r_center + h_center)))
        bar_lo, bar_hi = sorted((ya.to_value(r_top_cap), ya.to_value(r_bot_cap)))
        bar_lo_out, bar_hi_out = sorted((ya.to_value(r_top_cap - h_cap), ya.to_value(r_bot_cap + h_cap)))

        out.append({
            "x": float(xa.to_value(c_center)),
            "y": float(yv),
            "lo": float(lo_y),
            "hi": float(hi_y),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bar_lo_out),
            "bar_hi_outer": float(bar_hi_out),
            "px": [float(c_center), float(r_center)],
            "px_lo": float(r_top_cap),
            "px_hi": float(r_bot_cap),
            "note": f"{note}; cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out
