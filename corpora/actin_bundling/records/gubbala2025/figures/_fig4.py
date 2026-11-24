"""Shared calibrated axis detection and curve/histogram readers for Fig. 4 (PDF p. 7)."""
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


def _calibrate_hist_panel(img: np.ndarray):
    """Calibrate Fig. 4(A) or 4(B) histogram/KDE panel (linear x in [0, 1], linear y in [0, 3+])."""
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    neutral = (gray < 210.0) & (sat < 30.0)
    H, W = gray.shape

    # Left spine near col 340..430 (ignoring any faint border of left panel near col 15..30)
    lsp = _clusters(neutral[int(0.10 * H):int(0.75 * H), int(0.16 * W):int(0.24 * W)].mean(axis=0) > 0.70)
    c0 = lsp[0][2] + int(0.16 * W)
    # Right spine near col 0.90..0.98*W in upper empty half
    rsp = _clusters(neutral[int(0.10 * H):int(0.38 * H), int(0.90 * W):int(0.98 * W)].mean(axis=0) > 0.70)
    c1 = rsp[-1][2] + int(0.90 * W)
    # Top spine near row 0.02..0.09*H
    tsp = _clusters(neutral[int(0.02 * H):int(0.09 * H), int(c0) + 50:int(c1) - 50].mean(axis=1) > 0.70)
    r0 = tsp[0][2] + int(0.02 * H)
    # Bottom spine: last grey/nonwhite continuous horizontal line near row 0.78..0.88*H
    nonwhite_row = (gray[int(0.78 * H):int(0.88 * H), int(c0) + 30:int(c1) - 30] < 215.0).mean(axis=1) > 0.90
    r1_bot = float(np.where(nonwhite_row)[0][-1] + int(0.78 * H))

    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1_bot))
    # Left inward ticks: 3, 2, 1 (above y=0 bottom spine)
    yt_cls = _clusters(neutral[r0i + 30:r1i - 100, c0i + 8:c0i + 24].mean(axis=1) > 0.55)
    yt = [m + r0i + 30 for _, _, m in yt_cls if (m + r0i + 30) < r1i - 150]
    if len(yt) != 3:
        raise RuntimeError(f"Fig. 4(A/B): expected 3 upper y ticks (3, 2, 1), got {yt}")

    # Linear y axis fitted from the 3 clean ticks (3, 2, 1); verify y=0 sits at bottom spine!
    ya = D.Axis("linear", [(yt[0], 3.0), (yt[1], 2.0), (yt[2], 1.0)])
    # Linear x axis from left spine (0.0) to right spine (1.0)
    xa = D.Axis("linear", [(c0, 0.0), (c1, 1.0)])
    return xa, ya, (c0, r0, c1, r1_bot)


def read_hist_kde(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read solid black KDE curve P(kappa) in Fig. 4(A) or 4(B)."""
    xa, ya, (c0, r0, c1, r1) = _calibrate_hist_panel(img)
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    # Thick pure-black KDE curve (spines and ticks are medium grey > 110)
    kde = (gray < 75.0) & (sat < 40.0)
    r0i, r1i = int(round(r0)) + 8, int(round(r1)) - 2
    c0i, c1i = int(round(c0)) + 4, int(round(c1)) - 4

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        col_lo = max(c0i, int(round(c - 8)))
        col_hi = min(c1i, int(round(c + 8)))
        pts_c, pts_r, hw_list = [], [], []
        for col in range(col_lo, col_hi + 1):
            ys = np.where(kde[r0i:r1i, col])[0] + r0i
            if len(ys) >= 3:
                pts_c.append(float(col))
                pts_r.append(float(np.median(ys)))
                hw_list.append(0.5 * (ys.max() - ys.min()))
        if len(pts_c) < 2:
            out.append({"x": xv, "y": None, "note": f"black KDE line not found at kappa={xv}"})
            continue
        pts_c = np.array(pts_c)
        pts_r = np.array(pts_r)
        if len(pts_c) >= 4:
            sl, ic = np.polyfit(pts_c, pts_r, 1)
            r_center = float(sl * c + ic)
        else:
            sl = 0.0
            r_center = float(np.median(pts_r))
        line_hw = float(np.median(hw_list)) if hw_list else 10.0
        half_px = max(8.0, 0.5 * line_hw + ya.residual_px + 3.0)
        y = max(0.0, ya.to_value(r_center))
        hi = ya.to_value(r_center - half_px)
        lo = max(0.0, ya.to_value(r_center + half_px))
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": r_center - half_px, "px_hi": r_center + half_px,
            "note": f"y_res={ya.residual_px:.2f}px line_hw={line_hw:.1f}px",
        })
    return out


def read_hist_bars(img: np.ndarray, n_bins: int) -> list[dict]:
    """Read colored histogram bars across n_bins equal bins in [0, 1] in Fig. 4(A) or 4(B)."""
    xa, ya, (c0, r0, c1, r1) = _calibrate_hist_panel(img)
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    colored = sat > 45.0
    kde = (gray < 75.0) & (sat < 40.0)
    r0i = int(round(r0)) + 6
    r_zero = ya.to_pixel(0.0)
    r_bot = min(int(round(r1)) - 2, int(round(r_zero)))
    bw = c1 - c0

    out = []
    for i in range(n_bins):
        x_lo = i / n_bins
        x_hi = (i + 1) / n_bins
        x_mid = 0.5 * (x_lo + x_hi)
        ca = int(round(c0 + (x_lo + 0.22 / n_bins) * bw))
        cb = int(round(c0 + (x_hi - 0.22 / n_bins) * bw))
        ca = max(int(round(c0)) + 5, ca)
        cb = min(int(round(c1)) - 5, cb)
        if cb <= ca:
            cb = ca + 1

        # For each column inside the bin center, find topmost colored pixel
        col_tops = []
        kde_blocks = 0
        for col in range(ca, cb + 1):
            ys = np.where(colored[r0i:r_bot, col])[0]
            if len(ys) >= 2:
                top_r = int(ys[0]) + r0i
                # Check if black KDE sits immediately above top_r (partially occluding the bar top)
                if kde[max(r0i, top_r - 12):top_r, col].any():
                    kde_blocks += 1
                col_tops.append(top_r)
            else:
                col_tops.append(r_bot)

        # Taking the 25th percentile (near topmost) picks the unoccluded part of a bar when KDE cuts diagonally across it!
        top_clean = float(np.percentile(col_tops, 25))
        # Check if ALL columns in the bin had KDE directly on top of the bar top
        all_blocked = kde_blocks >= max(1, int(0.8 * len(col_tops)))
        c_mid = xa.to_pixel(x_mid)
        y_val = max(0.0, ya.to_value(top_clean))
        if y_val < 0.015:
            y_val = 0.0
            top_clean = r_zero
        half_px = 20.0 if all_blocked else max(8.0, ya.residual_px + 5.0)
        hi = ya.to_value(top_clean - half_px)
        lo = max(0.0, ya.to_value(top_clean + half_px))
        if y_val < lo:
            lo = y_val
        if y_val > hi:
            hi = y_val
        note = f"bin {i+1}/{n_bins} [{x_lo:.3f},{x_hi:.3f}]" + (" kde_overlap" if all_blocked else "")
        out.append({
            "x": round(x_mid, 5), "y": y_val, "lo": lo, "hi": hi,
            "px": [float(c_mid), float(top_clean)],
            "px_lo": float(top_clean - half_px), "px_hi": float(top_clean + half_px),
            "note": note,
        })
    return out


def _calibrate_4c(img: np.ndarray):
    """Calibrate Fig. 4(C) Experiment or Simulation panel (linear x in [0, 1], linear y in [0, 6])."""
    gray = img.astype(float).mean(axis=2)
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    black = (gray < 140.0) & (sat < 30.0)
    vsp = _clusters(black.mean(axis=0) > 0.55)
    hsp = _clusters(black.mean(axis=1) > 0.55)
    c0, c1 = vsp[0][2], vsp[-1][2]
    r0, r1 = hsp[0][2], hsp[-1][2]
    c0i, c1i, r0i, r1i = int(round(c0)), int(round(c1)), int(round(r0)), int(round(r1))

    # Left inward y ticks: 6, 4, 2 (plus bottom spine r1 = 0)
    yt_cls = _clusters(black[r0i + 30:r1i - 300, c0i + 12:c0i + 32].mean(axis=1) > 0.70)
    yt = [m + r0i + 30 for _, _, m in yt_cls]
    if len(yt) != 3:
        raise RuntimeError(f"Fig. 4(C): expected 3 upper y ticks (6, 4, 2), got {yt}")
    ya = D.Axis("linear", [(yt[0], 6.0), (yt[1], 4.0), (yt[2], 2.0), (r1, 0.0)])

    # X axis: left spine (0.0) to right spine (1.0), cross-checked with 0.75 tick closest to c0 + 0.75*(c1-c0)
    x75_exp = c0 + 0.75 * (c1 - c0)
    xt_cls = _clusters(black[r1i - 32:r1i - 12, int(x75_exp - 40):int(x75_exp + 40)].mean(axis=0) > 0.70)
    marks_x = [(c0, 0.0), (c1, 1.0)]
    if xt_cls:
        best_t = min(xt_cls, key=lambda g: abs(g[2] + int(x75_exp - 40) - x75_exp))
        marks_x.insert(1, (best_t[2] + int(x75_exp - 40), 0.75))
    xa = D.Axis("linear", marks_x)
    return xa, ya, (c0, r0, c1, r1)


def read_4c_curve(img: np.ndarray, at: list[float], panel: str, series: str) -> list[dict]:
    """Read No Actin ('no_actin'), Low ('low'), or High ('high') curve in Fig. 4(C) ('exp' or 'sim')."""
    xa, ya, (c0, r0, c1, r1) = _calibrate_4c(img)
    bw, bh = c1 - c0, r1 - r0
    r_ch = img[:, :, 0].astype(float)
    g_ch = img[:, :, 1].astype(float)
    b_ch = img[:, :, 2].astype(float)
    gray = (r_ch + g_ch + b_ch) / 3.0
    sat = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)

    if series == "no_actin":
        mask = (gray < 85.0) & (sat < 30.0)
    elif series == "low":
        d = np.sqrt((r_ch - 43.0) ** 2 + (g_ch - 126.0) ** 2 + (b_ch - 184.0) ** 2)
        mask = (d < 72.0) & ((b_ch - g_ch) > 20.0) & (r_ch < 115.0)
    elif series == "high":
        d = np.sqrt((r_ch - 129.0) ** 2 + (g_ch - 205.0) ** 2 + (b_ch - 188.0) ** 2)
        mask = (d < 62.0) & ((g_ch - r_ch) > 38.0) & (g_ch > 175.0) & (b_ch > 155.0) & ((b_ch - g_ch) < 6.0)
    else:
        raise ValueError(series)

    clean = mask.copy()
    # Keep plot interior (allow curve right down to y=0 at r1-2 and up near top at r0+6)
    clean[:int(r0) + 6, :] = False
    clean[int(r1) - 1:, :] = False
    clean[:, :int(c0) + 6] = False
    clean[:, int(c1) - 6:] = False

    # Mask out panel title, inset micrographs, inset labels, and pointer arrows!
    # Title ("Experiment" or "Simulation"): top left (r < r0 + 0.10*bh, c < c0 + 0.45*bw)
    clean[:int(r0 + 0.10 * bh), :int(c0 + 0.45 * bw)] = False

    if panel == "exp":
        # High inset + label + arrow: x in [0.02..0.31], y in [0.36..0.79] (r in [0.36..0.79]*bh)
        clean[int(r0 + 0.35 * bh):int(r0 + 0.72 * bh), int(c0 + 0.02 * bw):int(c0 + 0.30 * bw)] = False
        # High arrow below inset: x in [0.23..0.30], r in [0.70..0.78]*bh
        clean[int(r0 + 0.70 * bh):int(r0 + 0.78 * bh), int(c0 + 0.23 * bw):int(c0 + 0.30 * bw)] = False
        # Low inset + label: x in [0.30..0.59], r in [0.30..0.67]*bh
        clean[int(r0 + 0.30 * bh):int(r0 + 0.67 * bh), int(c0 + 0.30 * bw):int(c0 + 0.59 * bw)] = False
        # Low arrow: x in [0.57..0.66], r in [0.60..0.73]*bh
        clean[int(r0 + 0.60 * bh):int(r0 + 0.73 * bh), int(c0 + 0.57 * bw):int(c0 + 0.66 * bw)] = False
        # No Actin inset + label + arrow: x in [0.55..0.88], r in [0.04..0.39]*bh
        clean[int(r0 + 0.04 * bh):int(r0 + 0.39 * bh), int(c0 + 0.55 * bw):int(c0 + 0.88 * bw)] = False
    else:
        # Simulation panel insets:
        # High inset + label: x in [0.02..0.30], r in [0.40..0.74]*bh
        clean[int(r0 + 0.40 * bh):int(r0 + 0.74 * bh), int(c0 + 0.02 * bw):int(c0 + 0.30 * bw)] = False
        # High arrow: x in [0.28..0.34], r in [0.68..0.76]*bh
        clean[int(r0 + 0.68 * bh):int(r0 + 0.76 * bh), int(c0 + 0.28 * bw):int(c0 + 0.34 * bw)] = False
        # Low inset + label: x in [0.34..0.63], r in [0.30..0.66]*bh
        clean[int(r0 + 0.30 * bh):int(r0 + 0.66 * bh), int(c0 + 0.34 * bw):int(c0 + 0.63 * bw)] = False
        # Low arrow: x in [0.61..0.70], r in [0.48..0.59]*bh
        clean[int(r0 + 0.48 * bh):int(r0 + 0.59 * bh), int(c0 + 0.61 * bw):int(c0 + 0.70 * bw)] = False
        # No Actin inset + label + arrow: x in [0.52..0.86], r in [0.02..0.34]*bh
        clean[int(r0 + 0.02 * bh):int(r0 + 0.34 * bh), int(c0 + 0.52 * bw):int(c0 + 0.86 * bw)] = False

    # Also mask out left inward ticks (cols c0..c0+35) for black curve!
    if series == "no_actin":
        clean[:, :int(c0) + 38] = False
        # Mask bottom inward ticks (rows r1-35..r1) ONLY where curve isn't flat at y=0?
        # Actually at kappa < 0.50, black curve IS at y=0 (rows r1-35..r1-2), and bottom ticks are at 0.25, 0.50, 0.75.
        # At 0.75, black curve is high (P ~ 1.5), so excluding bottom ticks at x=0.75 avoids picking the tick!
        for xtick in (0.25, 0.50, 0.75):
            ct = int(round(xa.to_pixel(xtick)))
            if xtick == 0.75:
                clean[int(r1) - 40:, ct - 10:ct + 10] = False

    # Extract per-column curve center!
    # Because each curve is a single-valued function P(kappa) with smooth line + centered markers (~70 px tall),
    # in every column take the largest/main vertical cluster (gap=90 px to keep top & bottom of open marker together).
    pts_c, pts_r = [], []
    for col in range(int(c0) + 8, int(c1) - 8, 2):
        col_slice = clean[int(r0) + 6:int(r1) - 1, col]
        cls = _clusters(col_slice, min_len=2, gap=95)
        if not cls:
            continue
        # Pick cluster with most ink pixels
        best = max(cls, key=lambda g: int(col_slice[g[0]:g[1] + 1].sum()))
        npix = int(col_slice[best[0]:best[1] + 1].sum())
        if npix < 4:
            continue
        a, b = best[0] + int(r0) + 6, best[1] + int(r0) + 6
        pts_c.append(float(col))
        pts_r.append(0.5 * (a + b))

    pts_c = np.array(pts_c)
    pts_r = np.array(pts_r)

    out = []
    for xv in at:
        c = xa.to_pixel(xv)
        win = 26.0 if series == "no_actin" and xv > 0.72 else 34.0
        sel = np.abs(pts_c - c) <= win
        if sel.sum() < 2:
            # Handle y=0 baseline right next to left spine or edge
            if xv <= 0.05 and series == "no_actin":
                r_center = r1 - 12.0
                half_px = 20.0
                y = max(0.0, ya.to_value(r_center))
                out.append({
                    "x": float(xv), "y": y, "lo": 0.0, "hi": ya.to_value(r_center - half_px),
                    "px": [float(c), r_center], "px_lo": r_center - half_px, "px_hi": r1,
                    "note": "baseline near left spine",
                })
                continue
            out.append({"x": xv, "y": None, "note": f"no {series} curve pixels at kappa={xv}"})
            continue
        sc, sr = pts_c[sel], pts_r[sel]
        if len(sc) >= 5:
            sl, ic = np.polyfit(sc, sr, 1)
            resid = np.abs(sr - (sl * sc + ic))
            inliers = resid <= max(10.0, np.percentile(resid, 80) * 1.5)
            sl, ic = np.polyfit(sc[inliers], sr[inliers], 1)
            r_center = float(sl * c + ic)
            rmse = float(np.sqrt(np.mean((sr[inliers] - (sl * sc[inliers] + ic)) ** 2)))
        else:
            sl = 0.0
            r_center = float(np.median(sr))
            rmse = float(np.std(sr))
        # Clamp r_center to bottom spine if slightly below
        r_center = min(r_center, r1)
        half_px = max(10.0, 1.8 * rmse + ya.residual_px + 5.0)
        y = max(0.0, ya.to_value(r_center))
        hi = ya.to_value(r_center - half_px)
        lo = max(0.0, ya.to_value(r_center + half_px))
        if y < lo:
            lo = y
        out.append({
            "x": float(xv), "y": y, "lo": lo, "hi": hi,
            "px": [float(c), r_center], "px_lo": r_center - half_px, "px_hi": min(r1, r_center + half_px),
            "note": f"x_res={xa.residual_px:.2f}px y_res={ya.residual_px:.2f}px local_rmse={rmse:.1f}px",
        })
    return out
