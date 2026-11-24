"""Shared calibrated axis detection and curve reading helpers for SI Fig. S2 (PDF p. 18 / SI-9)."""
from __future__ import annotations

import math
import numpy as np
from boxgraph import digitize as D


def _clusters(mask: np.ndarray, min_len: int = 1, gap: int = 2) -> list[tuple[int, int, float]]:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    out = []
    for g in np.split(idx, splits):
        if len(g) >= min_len:
            out.append((int(g[0]), int(g[-1]), float(g.mean())))
    return out


def _black_mask(img: np.ndarray, gray_th: float = 120.0, sat_th: float = 30.0) -> np.ndarray:
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    return (gray < gray_th) & (sat < sat_th)


def _blue_mask(img: np.ndarray, shaded: bool = False) -> np.ndarray:
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    if shaded:
        # Shaded background blends RGB slightly toward yellow/pink: blue median ~(75, 140, 180)
        return (r < 118.0) & (g < 172.0) & (b > 142.0) & ((b - g) > 15.0) & ((b - r) > 55.0)
    d = np.sqrt((r - 44.0) ** 2 + (g - 127.0) ** 2 + (b - 184.0) ** 2)
    return (d < 70.0) & ((b - g) > 22.0) & (r < 110.0)


def _teal_mask(img: np.ndarray, shaded: bool = False) -> np.ndarray:
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    if shaded:
        # Over yellow/pink background: teal median ~(143, 200, 192) vs yellow bg ~(248, 249, 228) & pink bg ~(252, 243, 250)
        return (r < 178.0) & (g > 172.0) & (b > 158.0) & ((g - r) > 32.0) & ((b - r) > 22.0) & ((b - g) < 12.0)
    d = np.sqrt((r - 127.0) ** 2 + (g - 205.0) ** 2 + (b - 187.0) ** 2)
    return (d < 68.0) & ((g - r) > 35.0) & (b > 145.0)


def _interval_from_px(ya: D.Axis, r_center: float, half_px: float, slope_px: float = 0.0,
                      x_err_px: float = 1.0):
    total_half = max(math.hypot(half_px, abs(slope_px) * x_err_px), ya.residual_px + 2.0)
    y = ya.to_value(r_center)
    hi = ya.to_value(r_center - total_half)
    lo = ya.to_value(r_center + total_half)
    return y, lo, hi, r_center - total_half, r_center + total_half


def calibrate_panel(img: np.ndarray, panel: str):
    """Calibrate Fig. S2 'left', 'exp', or 'sim' log-log panel from spines and inward ticks."""
    black = _black_mask(img, gray_th=165.0, sat_th=35.0)
    vsp = _clusters(black.mean(axis=0) > 0.55)
    hsp = _clusters(black.mean(axis=1) > 0.55)
    c0, c1 = vsp[0][2], vsp[-1][2]
    r0, r1 = hsp[0][2], hsp[-1][2]
    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1))

    xt_raw = _clusters(black[r1i - 26:r1i - 10, c0i + 15:c1i - 15].mean(axis=0) > 0.65)
    xt = [m + c0i + 15 for _, _, m in xt_raw]
    yt_raw = _clusters(black[r0i + 20:r1i - 15, c0i + 10:c0i + 28].mean(axis=1) > 0.65)
    yt = [m + r0i + 20 for _, _, m in yt_raw]

    if panel == "left":
        if len(xt) != 4:
            raise RuntimeError(f"Fig. S2 left: expected 4 x ticks (100, 1000, 3000, 6000), got {xt}")
        if len(yt) != 4:
            raise RuntimeError(f"Fig. S2 left: expected 4 y ticks (4.5, 3.0, 1.5, 1.0), got {yt}")
        xa = D.Axis("log10", [(xt[0], 100.0), (xt[1], 1000.0), (xt[2], 3000.0), (xt[3], 6000.0)])
        ya = D.Axis("log10", [(yt[0], 4.5), (yt[1], 3.0), (yt[2], 1.5), (yt[3], 1.0)])
    elif panel == "exp":
        if len(xt) != 4:
            raise RuntimeError(f"Fig. S2 exp: expected 4 interior x ticks (100, 1000, 3000, 6000), got {xt}")
        if len(yt) != 4:
            raise RuntimeError(f"Fig. S2 exp: expected 4 y ticks (1.6, 1.4, 1.2, 1.0), got {yt}")
        xa = D.Axis("log10", [(c0, 10.0), (xt[0], 100.0), (xt[1], 1000.0), (xt[2], 3000.0), (xt[3], 6000.0)])
        ya = D.Axis("log10", [(yt[0], 1.6), (yt[1], 1.4), (yt[2], 1.2), (yt[3], 1.0)])
    elif panel == "sim":
        if len(xt) != 4:
            raise RuntimeError(f"Fig. S2 sim: expected 4 interior x ticks (100, 1000, 3000, 6000), got {xt}")
        if len(yt) != 4:
            raise RuntimeError(f"Fig. S2 sim: expected 4 y ticks (4.0, 3.0, 2.0, 1.0), got {yt}")
        xa = D.Axis("log10", [(c0, 10.0), (xt[0], 100.0), (xt[1], 1000.0), (xt[2], 3000.0), (xt[3], 6000.0)])
        ya = D.Axis("log10", [(yt[0], 4.0), (yt[1], 3.0), (yt[2], 2.0), (yt[3], 1.0)])
    else:
        raise ValueError(panel)
    return xa, ya, (c0, r0, c1, r1)


def read_left_alpha_third(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read black dashed alpha = 1/3 reference line in Fig. S2 left panel."""
    xa, ya, (c0, r0, c1, r1) = calibrate_panel(img, "left")
    black = _black_mask(img, gray_th=120.0, sat_th=30.0)
    bw, bh = c1 - c0, r1 - r0
    clean = black.copy()
    clean[:int(r0) + 12, :] = False
    clean[int(r1) - 12:, :] = False
    clean[:, :int(c0) + 12] = False
    clean[:, int(c1) - 12:] = False
    # Exclude legend (top-left)
    clean[:int(r0 + 0.20 * bh), :int(c0 + 0.30 * bw)] = False

    rr, cc = np.nonzero(clean)
    best_score, best_m, best_b = -1.0, -1.0, 3000.0
    for m in np.linspace(-1.15, -0.85, 121):
        proj = np.round(rr - m * cc).astype(int)
        bins = proj // 12
        vals, cnts = np.unique(bins, return_counts=True)
        for idx in np.argsort(cnts)[::-1][:5]:
            bval = vals[idx] * 12 + 6
            sel = np.abs(proj - bval) <= 12
            if sel.sum() < 200:
                continue
            span = cc[sel].max() - cc[sel].min()
            score = float(sel.sum()) * (1.0 if span > 0.75 * bw else 0.1)
            if score > best_score:
                best_score, best_m, best_b = score, m, float(bval)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        pts_c, pts_r = [], []
        for col in range(max(int(c0) + 12, int(round(c - 140))), min(int(c1) - 12, int(round(c + 140))) + 1):
            rp = best_m * col + best_b
            ra, rb = max(int(r0) + 10, int(round(rp - 24))), min(int(r1) - 10, int(round(rp + 24)))
            ys = np.where(clean[ra:rb, col])[0] + ra
            if 3 <= len(ys) <= 28:
                pts_c.append(float(col))
                pts_r.append(float(ys.mean()))
        if len(pts_c) < 10:
            out.append({"x": xv, "y": None, "note": "dashed line not found"})
            continue
        sc, sr = np.array(pts_c), np.array(pts_r)
        sl, ic = np.polyfit(sc, sr, 1)
        resid = np.abs(sr - (sl * sc + ic))
        inliers = resid <= max(4.0, np.percentile(resid, 75) * 1.5)
        sl, ic = np.polyfit(sc[inliers], sr[inliers], 1)
        r_center = float(sl * c + ic)
        rmse = float(np.sqrt(np.mean((sr[inliers] - (sl * sc[inliers] + ic)) ** 2)))
        half_px = max(6.0, 2.0 * rmse + ya.residual_px + 3.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, sl, max(xa.residual_px, 1.0))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.2f}px",
        })
    return out


def read_colored_curve(img: np.ndarray, at: list[float], panel: str, series: str, branch: str = "all") -> list[dict]:
    """Read dark-blue ('low') or teal ('high') curve in Fig. S2 ('left', 'exp', or 'sim').
    For panel='left', branch='sim' picks the upper circles band and branch='exp' picks the lower squares band."""
    xa, ya, (c0, r0, c1, r1) = calibrate_panel(img, panel)
    shaded = panel in ("exp", "sim")
    mask = _blue_mask(img, shaded=shaded) if series == "low" else _teal_mask(img, shaded=shaded)
    bw, bh = c1 - c0, r1 - r0

    clean = mask.copy()
    clean[:int(r0) + 10, :] = False
    clean[int(r1) - 4:, :] = False
    clean[:, :int(c0) + 4] = False
    clean[:, int(c1) - 4:] = False

    # Mask colored text labels ("High" in teal, "Low" in blue) strictly outside the curves!
    if panel == "exp":
        clean[int(r0 + 0.05 * bh):int(r0 + 0.25 * bh), int(c0 + 0.80 * bw):int(c0 + 0.94 * bw)] = False
        clean[int(r0 + 0.42 * bh):int(r0 + 0.59 * bh), int(c0 + 0.83 * bw):int(c0 + 0.96 * bw)] = False
    elif panel == "sim":
        # "High" teal text sits strictly above the High curve (r < 0.13*bh)
        clean[int(r0 + 0.04 * bh):int(r0 + 0.135 * bh), int(c0 + 0.68 * bw):int(c0 + 0.85 * bw)] = False
        # "Low" blue text sits strictly below the Low curve (r in 0.24..0.34*bh, c in 0.78..0.88*bw)
        clean[int(r0 + 0.24 * bh):int(r0 + 0.35 * bh), int(c0 + 0.78 * bw):int(c0 + 0.88 * bw)] = False

    pts_c, pts_r = [], []
    for col in range(int(c0) + 6, int(c1) - 6, 2):
        col_slice = clean[int(r0) + 10:int(r1) - 4, col]
        # Gap=90 keeps upper & lower rims of one open circle/square together (~70-85 px tall),
        # while keeping Sim (upper) and Exp (lower) cleanly separate in panel='left'!
        cls = _clusters(col_slice, min_len=2, gap=95)
        valid = [(a + int(r0) + 10, b + int(r0) + 10, int(col_slice[a:b + 1].sum()))
                 for a, b, _ in cls if col_slice[a:b + 1].sum() >= 6]
        if not valid:
            continue
        if panel == "left":
            if branch == "sim":
                # Uppermost cluster (Sim circles sit well above Exp squares for t >= 100)
                a, b, _ = valid[0]
            else:
                # Bottommost cluster (Exp squares sit along the bottom)
                a, b, _ = valid[-1]
        else:
            # Single curve in exp or sim panel: pick cluster with most colored pixels
            a, b, _ = max(valid, key=lambda g: g[2])
        pts_c.append(float(col))
        pts_r.append(0.5 * (a + b))

    pts_c = np.array(pts_c)
    pts_r = np.array(pts_r)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        win = 42.0 if xv < 120 else 36.0
        sel = np.abs(pts_c - c) <= win
        if sel.sum() < 3:
            # Try slightly wider window for sparse early log-spaced markers (e.g. t=20..80)
            sel = np.abs(pts_c - c) <= 75.0
        if sel.sum() < 2:
            out.append({"x": xv, "y": None, "note": f"no {series} markers near t={xv}"})
            continue
        sc, sr = pts_c[sel], pts_r[sel]
        if len(sc) >= 5:
            sl, ic = np.polyfit(sc, sr, 1)
            resid = np.abs(sr - (sl * sc + ic))
            inliers = resid <= max(14.0, np.percentile(resid, 80) * 1.5)
            sl, ic = np.polyfit(sc[inliers], sr[inliers], 1)
            r_center = float(sl * c + ic)
            rmse = float(np.sqrt(np.mean((sr[inliers] - (sl * sc[inliers] + ic)) ** 2)))
        else:
            sl = 0.0
            r_center = float(np.median(sr))
            rmse = float(np.std(sr))
        half_px = max(12.0, 1.8 * rmse + ya.residual_px + 6.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, sl, max(xa.residual_px, 1.0))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.1f}px",
        })
    return out
