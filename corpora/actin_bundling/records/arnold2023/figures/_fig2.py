"""Shared calibrated reader for Main Fig. 2 (PDF p. 2 / 128402-2): a/a_0 vs. t [s] (log-log).

Resolution-independent via scale factor s = (box_width) / 3463.5 (1.0 at 1200 dpi, 0.5 at 600 dpi).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_fig2(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box
    s = (c1 - c0) / 3463.5
    omin = max(3, int(round(8 * s)))
    omax = max(12, int(round(35 * s)))
    mlen_x = max(4, int(round(8 * s)))
    mlen_y = max(5, int(round(10 * s)))
    margin = 50.0 * s

    # X ticks (log10): 10, 50, 100 matched across bottom and top borders
    xt = [
        p
        for p, _ in match_ticks_both_sides(
            scan_ticks(image, box, "bottom_in", offset_min=omin, offset_max=omax, min_len=mlen_x),
            scan_ticks(image, box, "top_in", offset_min=omin, offset_max=omax, min_len=mlen_x),
            tol_px=max(2.0, 3.0 * s),
        )
        if c0 + margin < p < c1 - margin
    ]
    xa = build_axis("log10", xt, [10.0, 50.0, 100.0], max_residual_px=2.0)

    # Y ticks (log10): 3.0, 2.0, 1.0..0.4 matched across left & right + box edges 5.0 and 0.3
    yt = [
        p
        for p, _ in match_ticks_both_sides(
            scan_ticks(image, box, "left_in", offset_min=omin, offset_max=omax, min_len=mlen_y),
            scan_ticks(image, box, "right_in", offset_min=omin, offset_max=omax, min_len=mlen_y),
            tol_px=max(2.0, 3.0 * s),
        )
        if r0 + margin < p < r1 - margin
    ]
    ya = build_axis("log10", [r0] + yt + [r1], [5.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3], max_residual_px=2.0)
    return box, xa, ya, s


def _green_mask(image: np.ndarray, box, xa, ya, s: float) -> np.ndarray:
    c0, r0, c1, r1 = box
    rgb = image.astype(np.float64)
    g = (rgb[:, :, 1] - rgb[:, :, 0] > 55.0) & (rgb[:, :, 1] - rgb[:, :, 2] > 55.0)
    pad = max(6, int(round(15 * s)))
    valid = np.zeros_like(g)
    valid[int(r0) + pad : int(r1) - pad, int(c0) + pad : int(c1) - pad] = True
    valid[int(ya.to_pixel(0.48)) : int(r1), int(xa.to_pixel(24.0)) : int(c1)] = False
    return g & valid


def _black_interior(image: np.ndarray, box, s: float) -> np.ndarray:
    c0, r0, c1, r1 = box
    rgb = image.astype(np.float64)
    blk = (rgb.mean(axis=2) < 100.0) & ((rgb.max(axis=2) - rgb.min(axis=2)) < 35.0)
    pad = max(10, int(round(25 * s)))
    valid = np.zeros_like(blk)
    valid[int(r0) + pad : int(r1) - pad, int(c0) + pad : int(c1) - pad] = True
    return blk & valid


def _local_slope_center(cols: np.ndarray, rows: np.ndarray, target_col: float, s: float) -> tuple[float, float]:
    """Fit a local line row = m*(col - target_col) + b; return (b, 10-90% residual half-width)."""
    dc = cols - target_col
    min_h = max(4.0, 8.0 * s)
    if len(dc) >= 12 and np.ptp(dc) >= max(6.0, 12.0 * s):
        m, b = np.polyfit(dc, rows, 1)
        resid = rows - (m * dc + b)
        p10, p90 = float(np.percentile(resid, 10)), float(np.percentile(resid, 90))
        return float(b), max(0.5 * (p90 - p10), min_h)
    row_med = float(np.median(rows))
    p10, p90 = float(np.percentile(rows, 10)), float(np.percentile(rows, 90))
    return row_med, max(0.5 * (p90 - p10), min_h)


def _split_vertical_clusters(rr: np.ndarray, cc: np.ndarray, gap_px: float, min_px: int) -> list[tuple[np.ndarray, np.ndarray]]:
    if len(rr) == 0:
        return []
    u = np.sort(np.unique(rr))
    gaps = np.where(np.diff(u) >= gap_px)[0]
    cuts = [0] + [int(g) + 1 for g in gaps] + [len(u)]
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        r_min, r_max = u[a], u[b - 1]
        sel = (rr >= r_min) & (rr <= r_max)
        if sel.sum() >= min_px:
            out.append((rr[sel], cc[sel]))
    return out


def _fit_green_passive_line(green: np.ndarray, box, xa, ya, s: float) -> tuple[float, float]:
    c0, r0, c1, r1 = box
    pad = max(6, int(round(15 * s)))
    hw = max(25, int(round(60 * s)))
    pts_c, pts_r = [], []
    for t in [5.05, 7.55, 10.0, 15.0, 22.5, 33.0, 48.0, 68.0, 92.0, 115.0]:
        col = xa.to_pixel(t)
        ci = int(round(col))
        strip = green[int(r0) + pad : int(r1) - pad, max(int(c0) + pad, ci - hw) : min(int(c1) - pad, ci + hw + 1)]
        rr, cc = np.nonzero(strip)
        clusters = _split_vertical_clusters(rr, cc, gap_px=120.0 * s, min_px=max(5, int(15 * s * s)))
        if clusters:
            rc, _ = clusters[-1]
            pts_c.append(col)
            pts_r.append(float(np.median(rc)) + int(r0) + pad + 0.5)
    m, b = np.polyfit(np.array(pts_c), np.array(pts_r), 1)
    return float(m), float(b)


def _fit_passive_ref_dashed_line(blk: np.ndarray, box, xa, ya, green_mb: tuple[float, float], s: float) -> tuple[float, float, list[tuple[float, float]]]:
    gm, gb = green_mb
    lbl, _ = ndimage.label(blk)
    slices = ndimage.find_objects(lbl)
    cands = []
    s2 = s * s
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if not (14 * s <= h <= 55 * s and 35 * s <= w <= 115 * s):
            continue
        if w / max(h, 1) < 1.25:
            continue
        sub = lbl[sl] == i
        area = int(sub.sum())
        if not (400 * s2 <= area <= 3200 * s2):
            continue
        rr, cc = np.nonzero(sub)
        row = float(rr.mean()) + sl[0].start + 0.5
        col = float(cc.mean()) + sl[1].start + 0.5
        r_green = gm * col + gb
        if not (r_green + 220.0 * s <= row <= r_green + 365.0 * s):
            continue
        if col > xa.to_pixel(25.0) and row > ya.to_pixel(0.50) and row > r_green + 340.0 * s:
            continue
        cands.append((col, row))

    if len(cands) < 6:
        raise RuntimeError(f"Fig. 2 passive ref line: only {len(cands)} dash candidates found")

    pts = np.array(cands)
    best_inliers = np.zeros(len(pts), dtype=bool)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(pts[j, 0] - pts[i, 0]) < 400.0 * s:
                continue
            m = (pts[j, 1] - pts[i, 1]) / (pts[j, 0] - pts[i, 0])
            if not (-0.55 <= m <= -0.30):
                continue
            b = pts[i, 1] - m * pts[i, 0]
            resid = np.abs(pts[:, 1] - (m * pts[:, 0] + b))
            inliers = resid < 10.0 * s
            if inliers.sum() > best_inliers.sum():
                best_inliers = inliers

    in_pts = pts[best_inliers]
    if len(in_pts) < 8:
        raise RuntimeError(f"Fig. 2 passive ref line: RANSAC found only {len(in_pts)} inlier dashes")
    m, b = np.polyfit(in_pts[:, 0], in_pts[:, 1], 1)
    return float(m), float(b), [(float(c), float(r)) for c, r in in_pts]


def _read_green_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya, s = calibrate_fig2(image)
    c0, r0, c1, r1 = box
    pad = max(6, int(round(15 * s)))
    green = _green_mask(image, box, xa, ya, s)
    gm, gb = _fit_green_passive_line(green, box, xa, ya, s)

    out = []
    for t in at:
        col = xa.to_pixel(t)
        ci = int(round(col))
        half_w = int(round((80 if series in ("sim_piv_pe1e2", "sim_piv_pe5e3") else 65) * s))
        c_lo = max(int(c0) + pad, ci - half_w)
        c_hi = min(int(c1) - pad, ci + half_w + 1)
        strip = green[int(r0) + pad : int(r1) - pad, c_lo:c_hi]
        rr, cc = np.nonzero(strip)
        clusters = _split_vertical_clusters(rr, cc, gap_px=110.0 * s, min_px=max(5, int(15 * s * s)))
        if not clusters:
            out.append({"x": float(t), "y": None, "found": False, "note": f"no green pixels near t={t}"})
            continue

        r_pass_guide = gm * col + gb
        c_info = [(float(np.median(rc)) + int(r0) + pad + 0.5, rc, cc_l) for rc, cc_l in clusters]
        picked = None
        if series == "sim_passive":
            near = [c for c in c_info if abs(c[0] - r_pass_guide) <= 200.0 * s]
            if near:
                best = min(near, key=lambda c: abs(c[0] - r_pass_guide))
                picked = (best[1], best[2])
        elif series == "sim_piv_pe1e4":
            near = [c for c in c_info if 220.0 * s <= (r_pass_guide - c[0]) <= 680.0 * s]
            if near:
                best = min(near, key=lambda c: abs((r_pass_guide - c[0]) - 420.0 * s))
                picked = (best[1], best[2])
        elif series == "sim_piv_pe1e2":
            picked = (c_info[0][1], c_info[0][2])
        elif series == "sim_piv_pe5e3":
            high = [c for c in c_info if (r_pass_guide - c[0]) > 750.0 * s]
            if t <= 20.5 and len(high) >= 2:
                picked = (high[1][1], high[1][2])
            elif high:
                picked = (high[0][1], high[0][2])
        else:
            raise ValueError(series)

        if picked is None:
            out.append({"x": float(t), "y": None, "found": False, "note": f"cluster not resolved for {series} at t={t}"})
            continue

        rr_c, cc_c = picked
        rows_abs = rr_c + int(r0) + pad + 0.5
        cols_abs = cc_c + c_lo + 0.5
        row_est, half_px = _local_slope_center(cols_abs, rows_abs, col, s)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_est,
                half_px,
                x_val=float(t),
                note=f"Fig. 2 {series} ({len(rows_abs)} green px, half-width {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out


def _read_active_exp(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya, s = calibrate_fig2(image)
    c0, r0, c1, r1 = box
    pad = max(10, int(round(25 * s)))
    blk = _black_interior(image, box, s)

    lbl, _ = ndimage.label(blk)
    slices = ndimage.find_objects(lbl)
    comp_mask = None
    min_comp_area = int(round(40000 * s * s))
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        sub = lbl[sl] == i
        area = int(sub.sum())
        if area < min_comp_area:
            continue
        t_lo, t_hi = xa.to_value(sl[1].start), xa.to_value(sl[1].stop)
        y_hi = ya.to_value(sl[0].start)
        if y_hi < 1.05:
            continue
        match = False
        if series in ("exp_active_069", "exp_active_074") and t_lo < 7.0 and t_hi > 22.0:
            match = True
        elif series == "exp_active_059" and 6.5 <= t_lo <= 8.5 and t_hi < 18.0:
            match = True
        elif series == "exp_active_066" and 20.0 <= t_lo <= 26.0 and t_hi < 40.0:
            match = True
        elif series == "exp_active_065" and t_lo > 55.0:
            match = True
        if match:
            comp_mask = np.zeros_like(blk)
            comp_mask[sl] = sub
            break

    if comp_mask is None:
        raise RuntimeError(f"Fig. 2: could not isolate connected component for {series}")

    out = []
    half_w = max(12, int(round(30 * s)))
    for t in at:
        col = xa.to_pixel(t)
        ci = int(round(col))
        c_lo = max(int(c0) + pad, ci - half_w)
        c_hi = min(int(c1) - pad, ci + half_w + 1)
        strip = comp_mask[int(r0) + pad : int(r1) - pad, c_lo:c_hi]
        rr, cc = np.nonzero(strip)
        if len(rr) < max(8, int(20 * s * s)):
            out.append({"x": float(t), "y": None, "found": False, "note": f"no pixels for {series} at t={t}"})
            continue
        rows_abs = rr + int(r0) + pad + 0.5
        cols_abs = cc + c_lo + 0.5
        row_est, half_px = _local_slope_center(cols_abs, rows_abs, col, s)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_est,
                half_px,
                x_val=float(t),
                note=f"Fig. 2 {series} ({len(rows_abs)} px, half-width {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out


def _read_passive_ref_line(image: np.ndarray, at: list[float]) -> list[dict]:
    box, xa, ya, s = calibrate_fig2(image)
    green = _green_mask(image, box, xa, ya, s)
    blk = _black_interior(image, box, s)
    green_mb = _fit_green_passive_line(green, box, xa, ya, s)
    m, b, inliers = _fit_passive_ref_dashed_line(blk, box, xa, ya, green_mb, s)
    in_arr = np.array(inliers)
    resid = float(np.max(np.abs(in_arr[:, 1] - (m * in_arr[:, 0] + b))))

    out = []
    for t in at:
        col = xa.to_pixel(t)
        row = m * col + b
        half_px = max(resid + 6.0 * s, 5.0 * s)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row,
                half_px,
                x_val=float(t),
                note=f"Fig. 2 dashed alpha=1/3 ref line ({len(inliers)} dashes, fit max resid {resid:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out


def _read_exp_passive_band(image: np.ndarray, at: list[float]) -> list[dict]:
    box, xa, ya, s = calibrate_fig2(image)
    c0, r0, c1, r1 = box
    pad = max(10, int(round(25 * s)))
    green = _green_mask(image, box, xa, ya, s)
    blk = _black_interior(image, box, s)
    gm, gb = _fit_green_passive_line(green, box, xa, ya, s)
    rm, rb, _ = _fit_passive_ref_dashed_line(blk, box, xa, ya, (gm, gb), s)

    out = []
    half_w = max(14, int(round(35 * s)))
    for t in at:
        col = xa.to_pixel(t)
        ci = int(round(col))
        c_lo = max(int(c0) + pad, ci - half_w)
        c_hi = min(int(c1) - pad, ci + half_w + 1)
        r_top = max(int(r0) + pad, int(round(gm * col + gb - 15.0 * s)))
        r_bot = min(int(r1) - pad, int(round(rm * col + rb - 28.0 * s)))
        strip = blk[r_top:r_bot, c_lo:c_hi]
        rr, cc = np.nonzero(strip)
        if len(rr) < max(8, int(20 * s * s)):
            out.append({"x": float(t), "y": None, "found": False, "note": f"no passive exp pixels at t={t}"})
            continue
        rows_abs = rr + r_top + 0.5
        cols_abs = cc + c_lo + 0.5
        row_est, _ = _local_slope_center(cols_abs, rows_abs, col, s)
        dc = cols_abs - col
        m_loc = 0.5 * (gm + rm)
        resid = rows_abs - (m_loc * dc + row_est)
        p05, p95 = float(np.percentile(resid, 5)), float(np.percentile(resid, 95))
        half_px = max(0.5 * (p95 - p05), 14.0 * s)
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row_est,
                half_px,
                x_val=float(t),
                note=f"Fig. 2 combined passive exp band ({len(rows_abs)} px, band half-width {half_px:.1f} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    if series in ("sim_piv_pe1e2", "sim_piv_pe5e3", "sim_piv_pe1e4", "sim_passive"):
        return _read_green_series(image, at, series)
    if series.startswith("exp_active_"):
        return _read_active_exp(image, at, series)
    if series == "passive_ref_line":
        return _read_passive_ref_line(image, at)
    if series == "exp_passive_band":
        return _read_exp_passive_band(image, at)
    raise ValueError(f"unknown Fig. 2 series {series!r}")
