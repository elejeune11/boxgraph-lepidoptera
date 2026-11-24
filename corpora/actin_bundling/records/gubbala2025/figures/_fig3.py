"""Shared calibrated axis detection and curve reading helpers for Fig. 3 (PDF p. 6)."""
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


def _black_mask(img: np.ndarray, gray_th: float = 110.0, sat_th: float = 28.0) -> np.ndarray:
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    return (gray < gray_th) & (sat < sat_th)


def _blue_mask(img: np.ndarray) -> np.ndarray:
    """Dark blue (#2b7fb8, approx RGB 43, 127, 184) for low actin density."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    d = np.sqrt((r - 43.0) ** 2 + (g - 127.0) ** 2 + (b - 184.0) ** 2)
    return (d < 72.0) & ((b - g) > 22.0) & (r < 110.0)


def _teal_mask(img: np.ndarray) -> np.ndarray:
    """Light blue/teal (#7fcdbb, approx RGB 127, 205, 187) for high actin density."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    d = np.sqrt((r - 127.0) ** 2 + (g - 205.0) ** 2 + (b - 187.0) ** 2)
    return (d < 68.0) & ((g - r) > 35.0) & (b > 145.0)


def calibrate_main_panel(img: np.ndarray, panel: str):
    """Calibrate Fig. 3 top ('sim') or bottom ('exp') log-log panel from spines & inward ticks."""
    black = _black_mask(img, gray_th=110.0, sat_th=28.0)
    vsp = _clusters(black.mean(axis=0) > 0.55)
    hsp = _clusters(black.mean(axis=1) > 0.55)
    if len(vsp) < 2 or len(hsp) < 2:
        raise RuntimeError(f"Fig. 3 ({panel}): could not find plot spines")
    c0, c1 = vsp[0][2], vsp[-1][2]
    # Top crop may include top edge of bottom panel near bottom; bottom crop may include bottom edge of top panel near top
    if panel == "sim":
        # Choose the large panel (top ~64, bottom ~2832)
        h_valid = [h for h in hsp if h[2] < img.shape[0] - 15]
        r0, r1 = h_valid[0][2], h_valid[-1][2]
    else:
        # Bottom panel starts below top-panel bottom edge (first two spines are ~32 and ~52)
        h_valid = [h for h in hsp if h[2] > 40]
        r0, r1 = h_valid[0][2], h_valid[-1][2]

    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1))

    # Bottom inward ticks: 1000, 2000, 4000
    xt_raw = _clusters(black[r1i - 28:r1i - 8, c0i + 30:c1i - 30].mean(axis=0) > 0.75)
    xt = [m + c0i + 30 for _, _, m in xt_raw]
    if len(xt) != 3:
        raise RuntimeError(f"Fig. 3 ({panel}): expected 3 interior x ticks (1000, 2000, 4000), got {xt}")
    xa = D.Axis("log10", [(xt[0], 1000.0), (xt[1], 2000.0), (xt[2], 4000.0)])

    # Left inward ticks: top 4 unoccluded ticks 2.1, 1.8, 1.5, 1.2 (y=1.0 is merged with first marker)
    yt_raw = _clusters(black[r0i + 20:r1i - 300, c0i + 8:c0i + 26].mean(axis=1) > 0.75)
    yt = [m + r0i + 20 for _, _, m in yt_raw]
    if len(yt) != 4:
        raise RuntimeError(f"Fig. 3 ({panel}): expected 4 upper y ticks (2.1, 1.8, 1.5, 1.2), got {yt}")
    ya = D.Axis("log10", [(yt[0], 2.1), (yt[1], 1.8), (yt[2], 1.5), (yt[3], 1.2)])

    return xa, ya, (c0, r0, c1, r1)


def _interval_from_px(ya: D.Axis, r_center: float, half_px: float, slope_dy_dx_px: float = 0.0,
                      x_err_px: float = 1.0) -> tuple[float, float, float, float, float]:
    total_half = math.hypot(half_px, abs(slope_dy_dx_px) * x_err_px)
    total_half = max(total_half, ya.residual_px + 2.0)
    y = ya.to_value(r_center)
    # Row grows downward, so hi value is at r_center - total_half, lo value is at r_center + total_half
    hi = ya.to_value(r_center - total_half)
    lo = ya.to_value(r_center + total_half)
    return y, lo, hi, r_center - total_half, r_center + total_half


def read_sim_no_actin(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read black dashed alpha=1/3 reference line in Fig. 3 top panel."""
    xa, ya, (c0, r0, c1, r1) = calibrate_main_panel(img, "sim")
    black = _black_mask(img, gray_th=110.0, sat_th=28.0)

    # First find a rough diagonal guide from bottom-left (720, 1.0) to top-right (7200, 2.154)
    # using ONLY the image: scan clean columns away from insets/text to trace the dashed line!
    # Clean region: rows below top spine + 15, above bottom spine - 15, inside plot box.
    # Mask out the "Low" inset (c in [1380, 2060], r in [120, 940]) and "High" inset (c in [730, 1390], r in [980, 1800])
    # and the "Simulation" text (r > 2500, c > 2500).
    H, W = black.shape
    clean = black.copy()
    clean[:int(r0) + 12, :] = False
    clean[int(r1) - 12:, :] = False
    clean[:, :int(c0) + 12] = False
    clean[:, int(c1) - 12:] = False
    # Exclude insets & labels by relative box coordinates
    bw, bh = c1 - c0, r1 - r0
    clean[int(r0 + 0.03 * bh):int(r0 + 0.35 * bh), int(c0 + 0.24 * bw):int(c0 + 0.48 * bw)] = False  # Low inset
    clean[int(r0 + 0.34 * bh):int(r0 + 0.66 * bh), int(c0 + 0.03 * bw):int(c0 + 0.27 * bw)] = False  # High inset
    clean[int(r0 + 0.88 * bh):, int(c0 + 0.60 * bw):] = False  # Simulation label

    # Locate the dashed line via Hough/line scan across clean pixels:
    # We know it's the long diagonal black dashed line spanning bottom-left to top-right.
    # Find the diagonal band directly from image pixels using a 2D line projection over slopes ~ -0.8 to -1.0.
    rr, cc = np.nonzero(clean)
    best_cnt, best_m, best_b = -1, -0.85, 0.0
    for m in np.linspace(-0.98, -0.75, 116):
        proj = np.round(rr - m * cc).astype(int)
        # bin width 16 px
        bins = proj // 16
        vals, cnts = np.unique(bins, return_counts=True)
        # Prefer lines whose pixels span wide column range (the full diagonal, not short text)
        order = np.argsort(cnts)[::-1][:5]
        for idx in order:
            bval = vals[idx]
            sel = np.abs(proj - (bval * 16 + 8)) <= 18
            span = cc[sel].max() - cc[sel].min() if sel.sum() > 20 else 0
            score = int(sel.sum()) * (1.0 if span > 0.7 * bw else 0.2)
            if score > best_cnt:
                best_cnt, best_m, best_b = score, m, bval * 16 + 8

    # Refine locally around each target x using a +/- 140 px horizontal window within +/- 32 px of the detected line
    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        r_est = best_m * c + best_b
        c_lo = max(int(c0) + 15, int(round(c - 160)))
        c_hi = min(int(c1) - 15, int(round(c + 160)))
        pts_c, pts_r = [], []
        for col in range(c_lo, c_hi + 1):
            r_pred = best_m * col + best_b
            r_a = max(int(r0) + 12, int(round(r_pred - 30)))
            r_b = min(int(r1) - 12, int(round(r_pred + 30)))
            ys = np.where(clean[r_a:r_b, col])[0] + r_a
            if len(ys) >= 4 and len(ys) <= 35:
                pts_c.append(float(col))
                pts_r.append(float(ys.mean()))
        if len(pts_c) < 10:
            out.append({"x": xv, "y": None, "note": "dashed line pixels not found in local window"})
            continue
        pts_c = np.array(pts_c)
        pts_r = np.array(pts_r)
        # Robust linear fit in local window (rejects any stray touching text pixels)
        slope, intercept = np.polyfit(pts_c, pts_r, 1)
        resid = np.abs(pts_r - (slope * pts_c + intercept))
        inliers = resid <= max(4.0, np.percentile(resid, 75) * 1.5)
        slope, intercept = np.polyfit(pts_c[inliers], pts_r[inliers], 1)
        r_center = float(slope * c + intercept)
        rmse = float(np.sqrt(np.mean((pts_r[inliers] - (slope * pts_c[inliers] + intercept)) ** 2)))
        half_px = max(6.0, 2.0 * rmse + ya.residual_px + 3.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, slope, max(xa.residual_px, 1.0))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.2f}px",
        })
    return out


def read_sim_colored(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read dark-blue ('low') or teal ('high') simulation curve+circles in Fig. 3 top panel."""
    xa, ya, (c0, r0, c1, r1) = calibrate_main_panel(img, "sim")
    mask = _blue_mask(img) if series == "low" else _teal_mask(img)
    bw, bh = c1 - c0, r1 - r0
    clean = mask.copy()
    clean[:int(r0) + 10, :] = False
    clean[int(r1) - 6:, :] = False
    clean[:, :int(c0) + 2] = False
    clean[:, int(c1) - 2:] = False
    # Exclude inset boxes
    clean[int(r0 + 0.02 * bh):int(r0 + 0.36 * bh), int(c0 + 0.23 * bw):int(c0 + 0.48 * bw)] = False  # Low inset
    clean[int(r0 + 0.33 * bh):int(r0 + 0.67 * bh), int(c0 + 0.02 * bw):int(c0 + 0.27 * bw)] = False  # High inset

    # All extraneous annotations ('alpha ~ 0.25', inset borders, pointer arrows) lie STRICTLY ABOVE
    # the plotted curves. So at every column, the BOTTOMMOST cluster of colored pixels (with gap=95
    # to keep upper & lower halves of an open circle together, while separating higher text) IS the curve!
    # Also exclude thin pointer arrows (<= 14 px vertical span & separated from circles) by taking the
    # lowest cluster of at least 8 colored pixels.
    pts_c, pts_r = [], []
    for col in range(int(c0) + 6, int(c1) - 6, 3):
        col_slice = clean[int(r0) + 10:int(r1) - 6, col]
        cls = _clusters(col_slice, min_len=3, gap=95)
        if not cls:
            continue
        # Bottommost cluster with enough ink (>= 6 px)
        valid = [(a + int(r0) + 10, b + int(r0) + 10) for a, b, _ in cls if col_slice[a:b + 1].sum() >= 6]
        if not valid:
            continue
        a, b = valid[-1]
        # When pointer arrow touches the top of a circle (around col 1600-1900), b - a can exceed 135 px.
        # In that case, the circle+line occupies the bottom ~110 px of [a, b], so cap top at b - 110!
        if (b - a) > 130:
            a = b - 108
        pts_c.append(float(col))
        pts_r.append(0.5 * (a + b))

    pts_c = np.array(pts_c)
    pts_r = np.array(pts_r)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        win = 70.0 if c < c0 + 100 else 55.0
        sel = np.abs(pts_c - c) <= win
        if sel.sum() < 4:
            out.append({"x": xv, "y": None, "note": f"no {series} curve pixels near t={xv}"})
            continue
        sc, sr = pts_c[sel], pts_r[sel]
        slope, intercept = np.polyfit(sc, sr, 1)
        resid = np.abs(sr - (slope * sc + intercept))
        inliers = resid <= max(8.0, np.percentile(resid, 75) * 1.5)
        slope, intercept = np.polyfit(sc[inliers], sr[inliers], 1)
        r_center = float(slope * c + intercept)
        rmse = float(np.sqrt(np.mean((sr[inliers] - (slope * sc[inliers] + intercept)) ** 2)))
        half_px = max(8.0, 1.5 * rmse + ya.residual_px + 4.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, slope, max(xa.residual_px, 1.0))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.2f}px",
        })
    return out


def read_exp_colored(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read dark-blue ('low') or teal ('high') experimental squares in Fig. 3 bottom panel."""
    xa, ya, (c0, r0, c1, r1) = calibrate_main_panel(img, "exp")
    mask = _blue_mask(img) if series == "low" else _teal_mask(img)
    bw, bh = c1 - c0, r1 - r0
    clean = mask.copy()
    clean[:int(r0) + 10, :] = False
    clean[int(r1) - 4:, :] = False
    clean[:, :int(c0) + 4] = False
    clean[:, int(c1) - 4:] = False
    # Exclude insets (exact relative positions in plot box):
    # High inset (teal border): left side, upper-middle
    clean[int(r0 + 0.34 * bh):int(r0 + 0.65 * bh), int(c0 + 0.03 * bw):int(c0 + 0.29 * bw)] = False
    # Low inset (blue border): lower-right (starts well below the curves at r > 0.56*bh, c > 0.76*bw)
    clean[int(r0 + 0.56 * bh):int(r0 + 0.90 * bh), int(c0 + 0.76 * bw):int(c0 + 0.995 * bw)] = False
    # Exclude thin pointer arrows only in the empty space strictly between inset and curve:
    # Low arrow lies below curve at c in [0.71..0.77]*bw, r in [0.62..0.75]*bh
    clean[int(r0 + 0.62 * bh):int(r0 + 0.80 * bh), int(c0 + 0.71 * bw):int(c0 + 0.77 * bw)] = False
    # High arrow lies above curve at c in [0.27..0.42]*bw, r in [0.45..0.68]*bh
    clean[int(r0 + 0.45 * bh):int(r0 + 0.68 * bh), int(c0 + 0.27 * bw):int(c0 + 0.42 * bw)] = False

    # Collect per-column square band centers across the entire plot
    pts_c, pts_r = [], []
    for col in range(int(c0) + 10, int(c1) - 10, 3):
        col_slice = clean[int(r0) + 10:int(r1) - 4, col]
        cls = _clusters(col_slice, min_len=3, gap=120)
        if not cls:
            continue
        best = max(cls, key=lambda g: int(col_slice[g[0]:g[1] + 1].sum()))
        npix = int(col_slice[best[0]:best[1] + 1].sum())
        if npix < 12:
            continue
        a, b = best[0] + int(r0) + 10, best[1] + int(r0) + 10
        pts_c.append(float(col))
        pts_r.append(0.5 * (a + b))

    pts_c = np.array(pts_c)
    pts_r = np.array(pts_r)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        sel = np.abs(pts_c - c) <= 55
        if sel.sum() < 4:
            out.append({"x": xv, "y": None, "note": f"no {series} exp squares near t={xv}"})
            continue
        sc, sr = pts_c[sel], pts_r[sel]
        slope, intercept = np.polyfit(sc, sr, 1)
        resid = np.abs(sr - (slope * sc + intercept))
        inliers = resid <= max(14.0, np.percentile(resid, 80) * 1.5)
        slope, intercept = np.polyfit(sc[inliers], sr[inliers], 1)
        r_center = float(slope * c + intercept)
        rmse = float(np.sqrt(np.mean((sr[inliers] - (slope * sc[inliers] + intercept)) ** 2)))
        half_px = max(18.0, 2.0 * rmse + ya.residual_px + 10.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, slope, max(xa.residual_px, 1.0))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.2f}px",
        })
    return out


def calibrate_inset(img: np.ndarray):
    """Calibrate Fig. 3 bottom inset ('No actin', log-log: x in [10, 100], y in [1.0, 2.0])."""
    black = _black_mask(img, gray_th=120.0, sat_th=28.0)
    H, W = black.shape
    vsp = _clusters(black[int(0.08 * H):int(0.75 * H), :].mean(axis=0) > 0.75)
    c0, c1 = vsp[0][2], vsp[-1][2]
    c0i, c1i = int(round(c0)), int(round(c1))
    # Restrict horizontal spine search strictly to [c0..c1] so outer panel borders are never picked
    hsp = _clusters(black[:, c0i + 5:c1i - 5].mean(axis=1) > 0.85)
    # Inset top spine and bottom spine are the pair bounding the interior (distance > 0.5*H)
    r0, r1 = hsp[0][2], hsp[-1][2]
    if len(hsp) > 2:
        # Exclude any outer line above the inset box (which spans full image width > 0.9*W)
        inner = [h for h in hsp if black[int(round(h[2])), :].mean() < 0.90]
        if len(inner) >= 2:
            r0, r1 = inner[0][2], inner[-1][2]
    r0i, r1i = int(round(r0)), int(round(r1))

    xt = [m + c0i + 30 for _, _, m in _clusters(black[r1i - 35:r1i - 10, c0i + 30:c1i - 15].mean(axis=0) > 0.75)]
    if len(xt) != 2:
        raise RuntimeError(f"Fig. 3 inset: expected 2 interior x ticks (50, 100), got {xt}")
    xa = D.Axis("log10", [(c0, 10.0), (xt[0], 50.0), (xt[1], 100.0)])

    yt = [m + r0i + 30 for _, _, m in _clusters(black[r0i + 30:r1i - 20, c0i + 12:c0i + 42].mean(axis=1) > 0.75)]
    if len(yt) != 2:
        raise RuntimeError(f"Fig. 3 inset: expected 2 y ticks (2.0, 1.0), got {yt}")
    ya = D.Axis("log10", [(yt[0], 2.0), (yt[1], 1.0)])
    return xa, ya, (c0, r0, c1, r1)


def _find_inset_dashed_line(img: np.ndarray, c0: float, r0: float, c1: float, r1: float) -> tuple[float, float]:
    """Locate the upper black dashed alpha=1/3 line in Fig. 3 inset using white-above unoccluded columns."""
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    dark = (gray < 175.0) & (sat < 35.0)
    bw, bh = c1 - c0, r1 - r0

    # Three clean zones where pure white paper lies strictly above the dashed line:
    # Zone 1: below micrograph (r > r0 + 0.62*bh, c in c0 + [0.03..0.28]*bw)
    # Zone 2: between micrograph and 'alpha = 1/3' (c in c0 + [0.35..0.46]*bw, r > r0 + 0.05*bh)
    # Zone 3: right of 'alpha = 1/3' (c in c0 + [0.80..0.98]*bw, r > r0 + 0.02*bh)
    pts_c, pts_r = [], []
    zones = [
        (int(c0 + 0.03 * bw), int(c0 + 0.28 * bw), int(r0 + 0.62 * bh), int(r1 - 50)),
        (int(c0 + 0.35 * bw), int(c0 + 0.46 * bw), int(r0 + 0.05 * bh), int(r1 - 50)),
        (int(c0 + 0.80 * bw), int(c0 + 0.98 * bw), int(r0 + 0.02 * bh), int(r1 - 50)),
    ]
    for ca, cb, ra, rb in zones:
        for col in range(ca, cb, 2):
            cls = _clusters(dark[ra:rb, col], min_len=2, gap=2)
            if not cls:
                continue
            a, b, m = cls[0]
            if (b - a + 1) <= 12:
                pts_c.append(float(col))
                pts_r.append(float(m + ra))

    pts_c = np.array(pts_c)
    pts_r = np.array(pts_r)
    slope, intercept = np.polyfit(pts_c, pts_r, 1)
    for _ in range(4):
        resid = pts_r - (slope * pts_c + intercept)
        inliers = np.abs(resid) <= max(5.0, np.percentile(np.abs(resid), 60) * 1.5)
        slope, intercept = np.polyfit(pts_c[inliers], pts_r[inliers], 1)
    return float(slope), float(intercept)


def read_inset_ref(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read black dashed alpha=1/3 reference line in Fig. 3 inset."""
    xa, ya, (c0, r0, c1, r1) = calibrate_inset(img)
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    dark = (gray < 175.0) & (sat < 35.0)
    slope, intercept = _find_inset_dashed_line(img, c0, r0, c1, r1)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        pts_c, pts_r = [], []
        for col in range(max(int(c0) + 15, int(round(c - 90))), min(int(c1) - 15, int(round(c + 90))) + 1):
            rp = slope * col + intercept
            ra, rb = max(int(r0) + 6, int(round(rp - 12))), min(int(r1) - 10, int(round(rp + 12)))
            ys = np.where(dark[ra:rb, col])[0] + ra
            if 2 <= len(ys) <= 12:
                pts_c.append(float(col))
                pts_r.append(float(ys.mean()))
        if len(pts_c) < 8:
            out.append({"x": xv, "y": None, "note": "inset dashed line not found"})
            continue
        sl, ic = np.polyfit(pts_c, pts_r, 1)
        r_center = float(sl * c + ic)
        rmse = float(np.sqrt(np.mean((np.array(pts_r) - (sl * np.array(pts_c) + ic)) ** 2)))
        half_px = max(5.0, 2.0 * rmse + 3.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, sl, max(xa.residual_px, 2.5))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px local_rmse={rmse:.2f}px",
        })
    return out


def read_inset_exp(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read black open squares (No actin experiment) in Fig. 3 inset."""
    xa, ya, (c0, r0, c1, r1) = calibrate_inset(img)
    black = _black_mask(img, gray_th=130.0, sat_th=30.0)
    slope, intercept = _find_inset_dashed_line(img, c0, r0, c1, r1)
    bw, bh = c1 - c0, r1 - r0

    clean = black.copy()
    clean[:int(r0) + 15, :] = False
    clean[int(r1) - 12:, :] = False
    clean[:, :int(c0) + 25] = False
    clean[:, int(c1) - 15:] = False
    clean[int(r0 + 0.04 * bh):int(r0 + 0.58 * bh), int(c0 + 0.03 * bw):int(c0 + 0.34 * bw)] = False
    clean[int(r0 + 0.65 * bh):, int(c0 + 0.48 * bw):] = False

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        c_lo = max(int(c0) + 30, int(round(c - 32)))
        c_hi = min(int(c1) - 20, int(round(c + 32)))
        ys_all = []
        for col in range(c_lo, c_hi + 1):
            r_dash = slope * col + intercept
            r_min = int(round(r_dash + 25))
            r_max = min(int(r1) - 12, int(round(r_dash + 220)))
            if r_min >= r_max:
                continue
            ys = np.where(clean[r_min:r_max, col])[0] + r_min
            if len(ys) >= 2:
                ys_all.append(0.5 * (ys.min() + ys.max()))
        if len(ys_all) < 5:
            out.append({"x": xv, "y": None, "note": f"no inset squares near t={xv}"})
            continue
        ys_arr = np.array(ys_all)
        r_center = float(np.median(ys_arr))
        spread = float(np.percentile(ys_arr, 85) - np.percentile(ys_arr, 15)) * 0.5
        half_px = max(14.0, spread + 6.0)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, r_center, half_px, slope, max(xa.residual_px, 2.5))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"x_res={xa.residual_px:.2f}px spread_px={spread:.1f}",
        })
    return out
