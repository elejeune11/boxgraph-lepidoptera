"""Shared image calibration and curve/marker/box extraction helpers for bathe2008 figures.

Every coordinate, axis calibration, and reading is computed strictly from the rendered
image pixels and the printed tick labels visible on the figure.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


def _clusters(idx: np.ndarray, gap: int = 8) -> list[tuple[float, float, float]]:
    """Return list of (center, lo, hi) for contiguous index clusters separated by > gap."""
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    out = []
    for g in np.split(idx, splits):
        out.append((float(np.mean(g)), float(g[0]), float(g[-1])))
    return out


def _dark_mask(image: np.ndarray, thr: float = 160.0) -> np.ndarray:
    """Fast uint16 exact equivalent of image.mean(axis=2) < thr without float64 allocation."""
    s = image[:, :, 0].astype(np.uint16) + image[:, :, 1] + image[:, :, 2]
    return s < int(round(thr * 3.0))


def _box_erode(mask: np.ndarray, size: int) -> np.ndarray:
    """Exact separable 2D square box erosion using 1D minimum filters."""
    u = mask.view(np.uint8)
    e = ndimage.minimum_filter1d(u, size=size, axis=0, mode="constant", cval=0)
    e = ndimage.minimum_filter1d(e, size=size, axis=1, mode="constant", cval=0)
    return e.astype(bool)


def _box_dilate(mask: np.ndarray, size: int) -> np.ndarray:
    """Exact separable 2D square box dilation using 1D maximum filters."""
    u = mask.view(np.uint8)
    d = ndimage.maximum_filter1d(u, size=size, axis=0, mode="constant", cval=0)
    d = ndimage.maximum_filter1d(d, size=size, axis=1, mode="constant", cval=0)
    return d.astype(bool)


def find_plot_box(image: np.ndarray, dark_thr: float = 160.0, dark: np.ndarray | None = None) -> tuple[float, float, float, float]:
    """Locate outer plot border (c0, r0, c1, r1) from dark border lines spanning >45% of crop."""
    if dark is None:
        dark = _dark_mask(image, dark_thr)
    rc = _clusters(np.where(dark.mean(axis=1) > 0.45)[0])
    cc = _clusters(np.where(dark.mean(axis=0) > 0.45)[0])
    if len(rc) < 2 or len(cc) < 2:
        raise ValueError(f"Could not detect plot box borders: rc={rc}, cc={cc}")
    return cc[0][0], rc[0][0], cc[-1][0], rc[-1][0]


def _fit_log_axis_from_ticks(
    tick_px: list[float],
    candidates: list[float],
    approx_m: float,
    approx_b: float,
    tol_px: float = 8.0,
) -> tuple[float, float, float]:
    """Match detected tick pixels to expected log10 values and fit linear map log10(v) = m*px + b.

    Returns (m, b, max_resid_px).
    """
    matched_px = []
    matched_lv = []
    for v in candidates:
        lv = math.log10(v)
        px_pred = (lv - approx_b) / approx_m
        best = min(tick_px, key=lambda p: abs(p - px_pred))
        if abs(best - px_pred) <= tol_px:
            matched_px.append(best)
            matched_lv.append(lv)
    if len(matched_px) < 2:
        raise ValueError(f"Matched too few ticks ({len(matched_px)}) for axis calibration")
    m, b = np.polyfit(matched_px, matched_lv, 1)
    resid_px = float(np.max(np.abs((np.array(matched_lv) - b) / m - np.array(matched_px))))
    return float(m), float(b), resid_px


def _both_ticks(
    dark: np.ndarray,
    c0: float,
    r0: float,
    c1: float,
    r1: float,
    off_x: int = 20,
    off_y: int = 20,
    margin: int = 12,
    gap: int = 5,
) -> tuple[list[float], list[float]]:
    """Find x-ticks present on both top & bottom borders, and y-ticks on both left & right borders."""
    xb = [c + int(c0) + margin for c, _, _ in _clusters(np.where(dark[int(r1) - off_x, int(c0) + margin:int(c1) - margin])[0], gap=gap)]
    xt = [c + int(c0) + margin for c, _, _ in _clusters(np.where(dark[int(r0) + off_x, int(c0) + margin:int(c1) - margin])[0], gap=gap)]
    xs = [0.5 * (a + min(xt, key=lambda b: abs(b - a))) for a in xb if xt and min(abs(b - a) for b in xt) < 6]

    yl = [r + int(r0) + margin for r, _, _ in _clusters(np.where(dark[int(r0) + margin:int(r1) - margin, int(c0) + off_y])[0], gap=gap)]
    yr = [r + int(r0) + margin for r, _, _ in _clusters(np.where(dark[int(r0) + margin:int(r1) - margin, int(c1) - off_y])[0], gap=gap)]
    ys = [0.5 * (a + min(yr, key=lambda b: abs(b - a))) for a in yl if yr and min(abs(b - a) for b in yr) < 6]
    return xs, ys


def _reading_from_row_interval(
    x_val: float,
    c_px: float,
    r_mid: float,
    r_lo_px: float,
    r_hi_px: float,
    my: float,
    by: float,
    extra_px: float = 2.0,
    note: str = "",
) -> dict:
    """Convert pixel center and stroke extent [r_lo_px, r_hi_px] (r_lo_px <= r_hi_px) to log10 value & 95% interval.

    Note: rows increase downward, so my < 0; smaller row r_lo_px gives LARGER y (hi).
    """
    r_top_unc = r_lo_px - extra_px
    r_bot_unc = r_hi_px + extra_px
    y = 10.0 ** (my * r_mid + by)
    hi = 10.0 ** (my * r_top_unc + by)
    lo = 10.0 ** (my * r_bot_unc + by)
    if lo > hi:
        lo, hi = hi, lo
    lo = min(lo, y)
    hi = max(hi, y)
    return {
        "x": float(x_val),
        "y": float(y),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(c_px), float(r_mid)],
        "px_lo": float(r_top_unc),
        "px_hi": float(r_bot_unc),
        "note": note,
    }


# ---------------------------------------------------------------------------
# Fig. 3A (p. 4): kappa_B* vs N for alpha in {1e-1, 1e0, 1e1, 1e2, 1e3, 1e4}
# ---------------------------------------------------------------------------

def calibrate_fig3a(image: np.ndarray, dark: np.ndarray | None = None):
    if dark is None:
        dark = _dark_mask(image, 160)
    c0, r0, c1, r1 = find_plot_box(image, 160, dark=dark)
    xs, ys = _both_ticks(dark, c0, r0, c1, r1, off_x=20, off_y=20)
    # Printed axis ticks: N from 1 (left border) to 400 (right border); kappa_B* from 1 (bottom) to 5e5 (top)
    x_cand = [m * 10**e for e in (0, 1) for m in range(2, 10)] + [10, 100, 200, 300]
    x_cand = sorted(set(x_cand))
    mx0 = (math.log10(400.0) - 0.0) / (c1 - c0)
    bx0 = -mx0 * c0
    mx, bx, x_res = _fit_log_axis_from_ticks(xs, x_cand, mx0, bx0, tol_px=12.0)

    y_cand = [m * 10**e for e in (0, 1, 2, 3, 4) for m in range(1, 10)] + [1e5, 2e5, 3e5, 4e5]
    y_cand = sorted(set(v for v in y_cand if v > 1.0))
    my0 = (math.log10(5e5) - 0.0) / (r0 - r1)
    by0 = -my0 * r1
    my, by, y_res = _fit_log_axis_from_ticks(ys, y_cand, my0, by0, tol_px=12.0)
    return (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res)


def read_fig3a_curve(image: np.ndarray, at: list[float], curve_rank_from_bottom: int) -> list[dict]:
    """Read one of the 8 curves in Fig. 3A ordered from bottom (0=decoupled thick) to top (7=coupled thick).

    Rank mapping: 0=decoupled thick, 1=alpha=0.1, 2=alpha=1, 3=alpha=10,
                  4=alpha=100, 5=alpha=1e3, 6=alpha=1e4, 7=fully coupled thick.
    Uses top-down stroke counting at each requested N where all relevant curves are distinct.
    """
    dark = _dark_mask(image, 160)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig3a(image, dark=dark)
    dark_curve = _dark_mask(image, 150)
    out = []
    # Also compute analytical guide for disambiguating merged thick curves at small N
    alpha_map = {1: 0.1, 2: 1.0, 3: 10.0, 4: 100.0, 5: 1e3, 6: 1e4}
    # Margin 50 px clears inward major x-tick marks (~32 px) while staying well outside all curves (>=225 px)
    r_margin = 50
    for x in at:
        c_px = (math.log10(x) - bx) / mx
        ci = int(round(c_px))
        col = dark_curve[int(r0) + r_margin:int(r1) - r_margin, ci]
        segs = [(ctr + int(r0) + r_margin, lo + int(r0) + r_margin, hi + int(r0) + r_margin)
                for ctr, lo, hi in _clusters(np.where(col)[0], gap=4)]
        # Order bottom-to-top (descending row coordinate)
        segs_bt = sorted(segs, key=lambda s: -s[0])
        if len(segs_bt) == 8:
            ctr, lo_r, hi_r = segs_bt[curve_rank_from_bottom]
            note = f"all 8 curves separated; rank {curve_rank_from_bottom} from bottom; cal_res=({x_res:.2f},{y_res:.2f})px"
        elif curve_rank_from_bottom in alpha_map:
            # Match closest thin curve to the expected Eq. 2 region (only used to pick WHICH stroke, never as value)
            a_val = alpha_map[curve_rank_from_bottom]
            k_guide = x * (1.0 + (x - 1.0) / (1.0 + (x + math.sqrt(x)) / a_val))
            r_guide = (math.log10(k_guide) - by) / my
            # Pick stroke whose center is closest to r_guide, provided it is within 40 px
            best = min(segs_bt, key=lambda s: abs(s[0] - r_guide))
            if abs(best[0] - r_guide) > 40.0:
                out.append({"x": float(x), "y": None, "note": "curve merged or not found"})
                continue
            ctr, lo_r, hi_r = best
            # Widen interval if this stroke is thick (merged with a neighbor)
            note = f"{len(segs_bt)} strokes visible; stroke width={hi_r - lo_r + 1:.1f}px; cal_res=({x_res:.2f},{y_res:.2f})px"
        else:
            out.append({"x": float(x), "y": None, "note": f"only {len(segs_bt)} strokes visible"})
            continue
        extra = max(2.0, 2.0 * y_res)
        out.append(_reading_from_row_interval(x, c_px, ctr, lo_r, hi_r, my, by, extra_px=extra, note=note))
    return out


# ---------------------------------------------------------------------------
# Fig. 3B (p. 4): kappa_B* vs alpha for N in {4, 9, 16, 25, 36, 49, 64, 81, 100}
# plus Timoshenko dotted lines for N in {4, 25, 100}
# ---------------------------------------------------------------------------

def calibrate_fig3b(image: np.ndarray, dark: np.ndarray | None = None):
    if dark is None:
        dark = _dark_mask(image, 160)
    c0, r0, c1, r1 = find_plot_box(image, 160, dark=dark)
    xs, ys = _both_ticks(dark, c0, r0, c1, r1, off_x=20, off_y=20)
    x_cand = [m * 10**e for e in (-2, -1, 0, 1, 2) for m in range(1, 10)] + [1e3, 2e3, 3e3]
    x_cand = sorted(set(v for v in x_cand if v > 1e-2))
    mx0 = (math.log10(4e3) - (-2.0)) / (c1 - c0)
    bx0 = -2.0 - mx0 * c0
    mx, bx, x_res = _fit_log_axis_from_ticks(xs, x_cand, mx0, bx0, tol_px=12.0)

    # Bottom is 1.0, top is 2e4 (confirmed by major ticks 10^1, 10^2, 10^3, 10^4)
    y_cand = [m * 10**e for e in (0, 1, 2, 3) for m in range(1, 10)] + [1e4]
    y_cand = sorted(set(v for v in y_cand if v > 1.0))
    my0 = (math.log10(2e4) - 0.0) / (r0 - r1)
    by0 = -my0 * r1
    my, by, y_res = _fit_log_axis_from_ticks(ys, y_cand, my0, by0, tol_px=14.0)
    return (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res)


def _fig3b_inset_rect(c0: float, r0: float, c1: float, r1: float) -> tuple[int, int, int, int]:
    """Pixel bounding box of the inset inside Fig. 3B so we mask it when reading curves."""
    return (int(c0 + 0.02 * (c1 - c0)), int(r0 + 0.01 * (r1 - r0)),
            int(c0 + 0.43 * (c1 - c0)), int(r0 + 0.38 * (r1 - r0)))


def read_fig3b_solid_curve(image: np.ndarray, at: list[float], n_val: int) -> list[dict]:
    """Read one of the 9 solid curves in Fig. 3B (N in {4, 9, 16, 25, 36, 49, 64, 81, 100})."""
    dark_cal = _dark_mask(image, 160)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig3b(image, dark=dark_cal)
    dark = _dark_mask(image, 150)
    ic0, ir0, ic1, ir1 = _fig3b_inset_rect(c0, r0, c1, r1)
    # Build solid-curve mask: remove short isolated dashes of Timoshenko curves via connected component width!
    interior = np.zeros_like(dark)
    interior[int(r0) + 25:int(r1) - 25, int(c0) + 25:int(c1) - 25] = dark[int(r0) + 25:int(r1) - 25, int(c0) + 25:int(c1) - 25]
    interior[ir0:ir1, ic0:ic1] = False
    lbl, nlab = ndimage.label(interior)
    objs = ndimage.find_objects(lbl)
    solid = np.zeros_like(dark)
    for idx, sl in enumerate(objs, 1):
        if sl is None:
            continue
        w_cc = sl[1].stop - sl[1].start
        if w_cc > 300:
            solid[sl] |= (lbl[sl] == idx)

    out = []
    for a in at:
        c_px = (math.log10(a) - bx) / mx
        ci = int(round(c_px))
        col = solid[int(r0) + 25:int(r1) - 25, ci]
        segs = [(ctr + int(r0) + 25, lo + int(r0) + 25, hi + int(r0) + 25)
                for ctr, lo, hi in _clusters(np.where(col)[0], gap=4)]
        # Guide row from Eq. 2 (pinned c=1) purely to identify which of the 9 solid curves belongs to n_val
        k_guide = n_val * (1.0 + (n_val - 1.0) / (1.0 + (n_val + math.sqrt(n_val)) / a))
        r_guide = (math.log10(k_guide) - by) / my
        if not segs:
            out.append({"x": float(a), "y": None, "note": "no solid curve found"})
            continue
        best = min(segs, key=lambda s: abs(s[0] - r_guide))
        if abs(best[0] - r_guide) > 60.0:
            out.append({"x": float(a), "y": None, "note": f"nearest solid stroke {abs(best[0]-r_guide):.1f}px away"})
            continue
        ctr, lo_r, hi_r = best
        extra = max(2.0, 2.0 * y_res)
        out.append(_reading_from_row_interval(
            a, c_px, ctr, lo_r, hi_r, my, by, extra_px=extra,
            note=f"solid curve N={n_val}; cal_res=({x_res:.2f},{y_res:.2f})px"
        ))
    return out


def read_fig3b_timoshenko(image: np.ndarray, at: list[float], n_val: int) -> list[dict]:
    """Read one of the 3 dotted Timoshenko curves in Fig. 3B (N in {4, 25, 100}).

    Collects all short dash components belonging to the dotted curve near the Timoshenko trajectory
    and interpolates across gaps at each requested alpha.
    """
    dark_cal = _dark_mask(image, 160)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig3b(image, dark=dark_cal)
    dark = _dark_mask(image, 150)
    ic0, ir0, ic1, ir1 = _fig3b_inset_rect(c0, r0, c1, r1)
    interior = np.zeros_like(dark)
    interior[int(r0) + 25:int(r1) - 25, int(c0) + 25:int(c1) - 25] = dark[int(r0) + 25:int(r1) - 25, int(c0) + 25:int(c1) - 25]
    interior[ir0:ir1, ic0:ic1] = False

    # Find dash components (small blobs, 5 <= area <= 1500, width < 80 px)
    lbl, _ = ndimage.label(interior)
    objs = ndimage.find_objects(lbl)
    dash_pts = []
    for idx, sl in enumerate(objs, 1):
        if sl is None:
            continue
        w_cc = sl[1].stop - sl[1].start
        h_cc = sl[0].stop - sl[0].start
        if w_cc > 80 or h_cc > 80:
            continue
        comp = lbl[sl] == idx
        if comp.sum() < 10:
            continue
        rr, cc = np.nonzero(comp)
        r_c = float(rr.mean() + sl[0].start)
        c_c = float(cc.mean() + sl[1].start)
        # Check proximity to Timoshenko curve for n_val (for disambiguation among the 3 dotted curves)
        a_c = 10.0 ** (mx * c_c + bx)
        k_t = (n_val ** 2) / (1.0 + n_val / a_c)
        r_t = (math.log10(k_t) - by) / my
        if abs(r_c - r_t) < 45.0:
            dash_pts.append((c_c, r_c, float(rr.min() + sl[0].start), float(rr.max() + sl[0].start)))
    dash_pts.sort(key=lambda p: p[0])

    out = []
    for a in at:
        c_px = (math.log10(a) - bx) / mx
        left = [p for p in dash_pts if p[0] <= c_px]
        right = [p for p in dash_pts if p[0] >= c_px]
        if left and abs(left[-1][0] - c_px) < 12.0:
            p = left[-1]
            r_mid, lo_r, hi_r = p[1], p[2], p[3]
            extra = max(3.0, 2.0 * y_res)
            note = f"direct Timoshenko dash hit N={n_val}"
        elif left and right and (right[0][0] - left[-1][0]) < 150.0:
            p0, p1 = left[-1], right[0]
            frac = (c_px - p0[0]) / max(1e-6, p1[0] - p0[0])
            r_mid = p0[1] + frac * (p1[1] - p0[1])
            half_w = 0.5 * max(p0[3] - p0[2], p1[3] - p1[2])
            lo_r, hi_r = r_mid - half_w, r_mid + half_w
            extra = max(5.0, 2.0 * y_res)
            note = f"interpolated between adjacent dashes (gap={p1[0]-p0[0]:.0f}px) N={n_val}"
        else:
            out.append({"x": float(a), "y": None, "note": "no nearby Timoshenko dashes"})
            continue
        out.append(_reading_from_row_interval(a, c_px, r_mid, lo_r, hi_r, my, by, extra_px=extra, note=note))
    return out


# ---------------------------------------------------------------------------
# Fig. 3B inset (p. 4): crossover alpha** vs N (pinned squares, clamped circles)
# ---------------------------------------------------------------------------

def calibrate_fig3b_inset(image: np.ndarray, dark: np.ndarray | None = None):
    if dark is None:
        dark = _dark_mask(image, 160)
    c0, r0, c1, r1 = find_plot_box(image, 160, dark=dark)
    yl = [r + int(r0) + 12 for r, _, _ in _clusters(np.where(dark[int(r0) + 12:int(r1) - 12, int(c0) + 16])[0], gap=6)]
    yr = [r + int(r0) + 12 for r, _, _ in _clusters(np.where(dark[int(r0) + 12:int(r1) - 12, int(c1) - 16])[0], gap=6)]
    ys = [0.5 * (a + min(yr, key=lambda b: abs(b - a))) for a in yl if yr and min(abs(b - a) for b in yr) < 8]
    if len(ys) != 2:
        raise ValueError(f"Expected 2 major y ticks (10^2, 10^0) in Fig 3B inset, got {ys}")
    my = (2.0 - 0.0) / (ys[0] - ys[1])
    by = 0.0 - my * ys[1]

    xb = [c + int(c0) + 12 for c, _, _ in _clusters(np.where(dark[int(r1) - 16, int(c0) + 12:int(c1) - 12])[0], gap=6)]
    diffs = np.diff(xb)
    i_10 = int(np.argmax(diffs))
    c_10 = xb[i_10]
    c_20 = xb[i_10 + 1]
    mx0 = math.log10(2.0) / (c_20 - c_10)
    bx0 = 1.0 - mx0 * c_10
    x_cand = [4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    mx, bx, x_res = _fit_log_axis_from_ticks(xb, x_cand, mx0, bx0, tol_px=15.0)
    return (c0, r0, c1, r1), (mx, bx, x_res), (my, by, 1.0)


def read_fig3b_inset_series(image: np.ndarray, at: list[float], bc: str, branch: str) -> list[dict]:
    """Read crossover alpha** vs N from Fig. 3B inset.

    bc: 'pinned' (squares) or 'clamped' (circles)
    branch: 'dec' (decoupled-to-intermediate, bottom) or 'fc' (fully coupled-to-intermediate, top)
    """
    dark = _dark_mask(image, 160)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig3b_inset(image, dark=dark)
    c_fac = 1.0 if bc == "pinned" else 4.0

    out = []
    for n_val in at:
        c_px = (math.log10(n_val) - bx) / mx
        if n_val == 4:
            a_guide = c_fac * 3.0
        elif branch == "dec":
            a_guide = c_fac * (n_val + math.sqrt(n_val)) / (n_val - 2.0)
        else:
            a_guide = c_fac * (n_val + math.sqrt(n_val)) * (n_val - 2.0) / n_val
        r_guide = (math.log10(a_guide) - by) / my

        ci = int(round(c_px))
        ri = int(round(r_guide))
        # Narrow window at N=9 where pinned_fc (~8.5) and clamped_dec (~6.5) sit close together
        win = 42 if (n_val == 9 and ((bc == "pinned" and branch == "fc") or (bc == "clamped" and branch == "dec"))) else 80
        r_win0 = max(int(r0) + 25, ri - win)
        r_win1 = min(int(r1) - 25, ri + win)
        c_win0 = max(int(c0) + 20, ci - 22)
        c_win1 = min(int(c1) - 20, ci + 22)
        patch = dark[r_win0:r_win1, c_win0:c_win1]

        row_profile = patch.sum(axis=1)
        active = np.where(row_profile >= 3)[0]
        segs = _clusters(active, gap=3)
        if not segs:
            out.append({"x": float(n_val), "y": None, "note": "no marker found"})
            continue
        if n_val == 9 and len(segs) >= 2 and bc == "clamped" and branch == "dec":
            # Lower marker (larger row) is the clamped-dec circle; upper marker is pinned-fc square
            best = max(segs, key=lambda s: s[0])
        elif n_val == 9 and len(segs) >= 2 and bc == "pinned" and branch == "fc":
            best = min(segs, key=lambda s: s[0])
        else:
            best = min(segs, key=lambda s: abs((s[0] + r_win0) - r_guide))
        lo_r = best[1] + r_win0
        hi_r = best[2] + r_win0
        r_mid = 0.5 * (lo_r + hi_r)
        half_unc = max(6.0, 0.35 * (hi_r - lo_r))
        out.append(_reading_from_row_interval(
            n_val, c_px, r_mid, r_mid - half_unc, r_mid + half_unc, my, by, extra_px=2.0,
            note=f"Fig 3B inset {bc} {branch} marker (span={hi_r-lo_r+1:.0f}px)"
        ))
    return out


# ---------------------------------------------------------------------------
# Fig. 4 (p. 6): experimental symbols + 3 theoretical curves (N=21, 27, 33)
# ---------------------------------------------------------------------------

def calibrate_fig4(image: np.ndarray, dark: np.ndarray | None = None):
    if dark is None:
        dark = _dark_mask(image, 180)
    c0, r0, c1, r1 = find_plot_box(image, 180, dark=dark)
    xs, ys = _both_ticks(dark, c0, r0, c1, r1, off_x=20, off_y=20)
    x_cand = [m * 1e-4 for m in range(3, 10)] + [m * 10**e for e in (-3, -2, -1) for m in range(1, 10)]
    mx0 = (0.0 - math.log10(2e-4)) / (c1 - c0)
    bx0 = 0.0 - mx0 * c1
    mx, bx, x_res = _fit_log_axis_from_ticks(xs, x_cand, mx0, bx0, tol_px=12.0)

    y_cand = [3e-23, 2e-23, 1e-23, 9e-24, 8e-24, 7e-24, 6e-24, 5e-24, 4e-24, 3e-24, 2e-24,
              1e-24, 9e-25, 8e-25, 7e-25, 6e-25]
    my0 = (math.log10(5e-23) - math.log10(5e-25)) / (r0 - r1)
    by0 = math.log10(5e-25) - my0 * r1
    my, by, y_res = _fit_log_axis_from_ticks(ys, y_cand, my0, by0, tol_px=12.0)
    return (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res)


def read_fig4_theory(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read theoretical curves from Fig. 4: which in {'n27' (solid), 'n21' (lower dashed), 'n33' (upper dashed)}."""
    dark = _dark_mask(image, 180)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig4(image, dark=dark)
    interior = np.zeros_like(dark)
    interior[int(r0) + 32:int(r1) - 32, int(c0) + 32:int(c1) - 32] = dark[int(r0) + 32:int(r1) - 32, int(c0) + 32:int(c1) - 32]
    # Fast exact separable equivalent of 22x22 opening followed by 18-px dilation (masks filled symbols)
    opened = _box_dilate(_box_erode(interior, 22), 22)
    markers = _box_dilate(opened, 37)
    curves_only = interior & (~markers)

    lbl, _ = ndimage.label(curves_only)
    objs = ndimage.find_objects(lbl)
    solid_mask = np.zeros_like(dark)
    dashed_mask = np.zeros_like(dark)
    for idx, sl in enumerate(objs, 1):
        if sl is None:
            continue
        w_cc = sl[1].stop - sl[1].start
        comp = lbl[sl] == idx
        if w_cc > 500:
            solid_mask[sl] |= comp
        elif w_cc > 15:
            dashed_mask[sl] |= comp

    solid_cols, solid_rows = [], []
    for ci in range(int(c0) + 40, int(c1) - 40, 4):
        rr = np.where(solid_mask[:, ci])[0]
        if len(rr):
            solid_cols.append(float(ci))
            solid_rows.append(float(np.mean(rr)))
    solid_cols = np.array(solid_cols)
    solid_rows = np.array(solid_rows)

    out = []
    for x in at:
        c_px = (math.log10(x) - bx) / mx
        ci = int(round(c_px))
        if which == "n27":
            rr = np.where(solid_mask[:, ci])[0]
            if len(rr) == 0:
                # Local linear interpolation across occluded marker region
                r_mid = float(np.interp(c_px, solid_cols, solid_rows))
                lo_r, hi_r = r_mid - 6.0, r_mid + 6.0
                out.append(_reading_from_row_interval(x, c_px, r_mid, lo_r, hi_r, my, by, extra_px=4.0,
                                                      note="solid curve interpolated across marker"))
                continue
            r_mid, lo_r, hi_r = float(np.mean(rr)), float(rr.min()), float(rr.max())
            out.append(_reading_from_row_interval(x, c_px, r_mid, lo_r, hi_r, my, by, extra_px=max(2.0, 2.0 * y_res),
                                                  note=f"Fig 4 solid curve N=27; cal_res=({x_res:.2f},{y_res:.2f})px"))
        else:
            pts_c, pts_r, pts_w = [], [], []
            for c_scan in range(max(int(c0) + 35, ci - 90), min(int(c1) - 35, ci + 91), 2):
                rr = np.where(dashed_mask[:, c_scan])[0]
                if len(rr) == 0:
                    continue
                r_s = float(np.interp(c_scan, solid_cols, solid_rows))
                sel = rr[rr < r_s - 25] if which == "n33" else rr[rr > r_s + 25]
                if len(sel):
                    pts_c.append(float(c_scan))
                    pts_r.append(float(np.mean(sel)))
                    pts_w.append(float(sel.max() - sel.min()))
            if len(pts_c) < 2:
                out.append({"x": float(x), "y": None, "note": f"dashed curve {which} not found"})
                continue
            slope, icpt = np.polyfit(pts_c, pts_r, 1)
            r_mid = float(slope * c_px + icpt)
            half_w = 0.5 * float(np.median(pts_w))
            out.append(_reading_from_row_interval(x, c_px, r_mid, r_mid - half_w, r_mid + half_w, my, by,
                                                  extra_px=max(4.0, 2.0 * y_res),
                                                  note=f"Fig 4 dashed curve {which}; cal_res=({x_res:.2f},{y_res:.2f})px"))
    return out


def _detect_fig4_exp_markers(image: np.ndarray):
    """Detect all 20 experimental marker centroids, radii, interior mask, and axis calibration in Fig. 4."""
    dark_cal = _dark_mask(image, 180)
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig4(image, dark=dark_cal)
    dark = _dark_mask(image, 160)
    ri0, ri1 = int(r0) + 35, int(r1) - 35
    ci0, ci1 = int(c0) + 35, int(c1) - 35
    interior = np.zeros_like(dark)
    interior[ri0:ri1, ci0:ci1] = dark[ri0:ri1, ci0:ci1]

    # Eroding by 15x15 removes thin curves (~8 px) and isolates filled marker cores;
    # running EDT only on each padded core bounding box is ~100x faster and exact at full resolution.
    cores = _box_erode(interior[ri0:ri1, ci0:ci1], 15)
    clbl, _ = ndimage.label(cores)
    raw_pts = []
    pad = 50
    for c_idx, sl in enumerate(ndimage.find_objects(clbl), 1):
        if sl is None:
            continue
        br0 = max(ri0, sl[0].start + ri0 - pad)
        br1 = min(ri1, sl[0].stop + ri0 + pad)
        bc0 = max(ci0, sl[1].start + ci0 - pad)
        bc1 = min(ci1, sl[1].stop + ci0 + pad)
        sub_m = interior[br0:br1, bc0:bc1]
        dist = ndimage.distance_transform_edt(sub_m)
        dmax = ndimage.maximum_filter1d(ndimage.maximum_filter1d(dist, size=27, axis=0), size=27, axis=1)
        local_max = (dist == dmax) & (dist >= 12.0)
        lbl, npeaks = ndimage.label(local_max)
        for idx in range(1, npeaks + 1):
            rr, cc = np.nonzero(lbl == idx)
            r_loc = float(rr.mean())
            c_loc = float(cc.mean())
            rad = float(dist[int(round(r_loc)), int(round(c_loc))])
            r_pk = r_loc + br0
            c_pk = c_loc + bc0
            if clbl[int(round(r_pk)) - ri0, int(round(c_pk)) - ci0] != c_idx:
                continue
            r0w, r1w = max(0, int(r_pk - rad * 1.4)), min(dark.shape[0], int(r_pk + rad * 1.4) + 1)
            c0w, c1w = max(0, int(c_pk - rad * 1.4)), min(dark.shape[1], int(c_pk + rad * 1.4) + 1)
            yy, xx = np.ogrid[r0w:r1w, c0w:c1w]
            in_disk = ((yy - r_pk) ** 2 + (xx - c_pk) ** 2) <= (rad * 1.25) ** 2
            m_sub = interior[r0w:r1w, c0w:c1w] & in_disk
            ry, rx = np.nonzero(m_sub)
            r_c = float(ry.mean() + r0w)
            c_c = float(rx.mean() + c0w)
            raw_pts.append((c_c, r_c, rad))

    # Deduplicate twin sub-peaks within 18 px (same diamond/square symbol)
    merged = []
    for c_c, r_c, rad in sorted(raw_pts, key=lambda p: -p[2]):
        if all(math.hypot(c_c - mc, r_c - mr) > 18.0 for mc, mr, _ in merged):
            merged.append((c_c, r_c, rad))

    merged.sort(key=lambda p: p[0])
    return merged, interior, (mx, bx, x_res), (my, by, y_res)


def read_fig4_exp_markers(image: np.ndarray, at: list[float]) -> list[dict]:
    """Extract all 20 experimental data points (filled black symbols) from Fig. 4."""
    merged, _, (mx, bx, x_res), (my, by, y_res) = _detect_fig4_exp_markers(image)
    out = []
    for c_c, r_c, rad in merged:
        x_val = 10.0 ** (mx * c_c + bx)
        unc_px = max(5.0, 0.3 * rad, 2.0 * y_res)
        out.append(_reading_from_row_interval(
            x_val, c_c, r_c, r_c - unc_px, r_c + unc_px, my, by, extra_px=1.0,
            note=f"Fig 4 experimental symbol (rad={rad:.1f}px)"
        ))
    return out


def read_fig4_exp_by_shape(image: np.ndarray, at: list[float], shape_key: str) -> list[dict]:
    """Extract one fascin-spacing symbol series from Fig. 4 purely by geometric shape matching.

    shape_key in {'d40' (squares), 'd56' (circles), 'd68' (diamonds),
                  'd225' (pointed-up triangles), 'd412' (pointed-down triangles)}.
    """
    merged, interior, (mx, bx, x_res), (my, by, y_res) = _detect_fig4_exp_markers(image)
    h_sz = 36
    yy, xx = np.mgrid[-h_sz:h_sz + 1, -h_sz:h_sz + 1]
    templates = {
        "d40":  (np.abs(xx) <= 30) & (np.abs(yy) <= 30),
        "d56":  (xx ** 2 + yy ** 2) <= 31 ** 2,
        "d68":  (np.abs(xx) + np.abs(yy)) <= 35,
        "d225": (yy <= 20) & (yy >= -42 + 1.73 * np.abs(xx)),
        "d412": (yy >= -23) & (yy <= 44 - 1.73 * np.abs(xx)),
    }
    roi = np.zeros_like(templates["d40"])
    for tmpl in templates.values():
        roi |= tmpl
    roi = _box_dilate(roi, 7)

    labels = {
        "d40": "square (delta=40 nm)",
        "d56": "circle (delta=56 nm)",
        "d68": "diamond (delta=68 nm)",
        "d225": "up-triangle (delta=225 nm)",
        "d412": "down-triangle (delta=412 nm)",
    }

    out = []
    for c_c, r_c, rad in merged:
        ri, ci = int(round(r_c)), int(round(c_c))
        patch = interior[ri - h_sz:ri + h_sz + 1, ci - h_sz:ci + h_sz + 1]
        scores = {}
        for k, tmpl in templates.items():
            inter = np.logical_and(patch & roi, tmpl).sum()
            union = np.logical_or(patch & roi, tmpl).sum()
            scores[k] = float(inter) / float(union)
        best_k = max(scores, key=scores.get)
        # Disambiguate touching circles vs squares via 4 diagonal corners at (+/-27, +/-27) px
        # (distance 38.2 px > circle radius 31 px, inside square half-width 30 px)
        if best_k in ("d40", "d56"):
            c4 = sum(bool(patch[h_sz + dr, h_sz + dc]) for dr in (-27, 27) for dc in (-27, 27))
            best_k = "d40" if c4 >= 3 else "d56"
        if best_k != shape_key:
            continue
        x_val = 10.0 ** (mx * c_c + bx)
        unc_px = max(5.0, 0.3 * rad, 2.0 * y_res)
        out.append(_reading_from_row_interval(
            x_val, c_c, r_c, r_c - unc_px, r_c + unc_px, my, by, extra_px=1.0,
            note=f"Fig 4 {labels[shape_key]} (rad={rad:.1f}px)"
        ))
    return out


# ---------------------------------------------------------------------------
# Fig. 5 (p. 6): State diagram (alpha vs N) colored bounding boxes (a)-(g)
# ---------------------------------------------------------------------------

def calibrate_fig5(image: np.ndarray, dark: np.ndarray | None = None):
    if dark is None:
        dark = _dark_mask(image, 210)
    c0, r0, c1, r1 = find_plot_box(image, 210, dark=dark)
    xs, ys = _both_ticks(dark, c0, r0, c1, r1, off_x=22, off_y=22)
    x_cand = [m * 10**e for e in (0, 1, 2, 3) for m in range(1, 10)]
    x_cand = sorted(set(v for v in x_cand if v > 1.0))
    mx0 = 4.0 / (c1 - c0)
    bx0 = -mx0 * c0
    mx, bx, x_res = _fit_log_axis_from_ticks(xs, x_cand, mx0, bx0, tol_px=15.0)

    y_cand = [m * 10**e for e in (-4, -3, -2, -1, 0, 1, 2, 3) for m in range(1, 10)]
    y_cand = sorted(set(v for v in y_cand if v > 1e-4))
    my0 = 8.0 / (r0 - r1)
    by0 = -4.0 - my0 * r1
    my, by, y_res = _fit_log_axis_from_ticks(ys, y_cand, my0, by0, tol_px=15.0)
    return (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res)


def read_fig5_box(image: np.ndarray, at: list[float], box_key: str) -> list[dict]:
    """Read the lower-left (N_min, alpha_min) and upper-right (N_max, alpha_max) bounds of box (a)-(g) in Fig. 5."""
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_fig5(image)
    ir0, ir1 = int(r0) + 15, int(r1) - 15
    ic0, ic1 = int(c0) + 15, int(c1) - 15
    sub = image[ir0:ir1, ic0:ic1]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    masks = {
        "a": (B > 180) & (R < 80) & (G < 80),
        "b": (R > 180) & (G < 60) & (B < 60),
        "c": (R > 200) & (G > 75) & (G < 170) & (B < 60),
        "d": (B > 180) & (G > 140) & (R < 80),
        "e": (G > 180) & (R < 80) & (B < 80),
        "f": (R > 180) & (B > 180) & (G < 80),
        "g": (R > 200) & (G > 200) & (B < 100),
    }
    mask = masks[box_key]
    lbl, nlab = ndimage.label(mask)
    if nlab == 0:
        raise ValueError(f"Fig 5 box {box_key} not found")
    sizes = ndimage.sum(mask, lbl, range(1, nlab + 1))
    keep = [i + 1 for i, sz in enumerate(sizes) if sz >= 80]
    if not keep:
        keep = [int(np.argmax(sizes)) + 1]
    rect = np.isin(lbl, keep)

    rr, cc = np.nonzero(rect)
    c_min_px, c_max_px = int(cc.min()), int(cc.max())
    r_top_px, r_bot_px = int(rr.min()), int(rr.max())
    # Measure exact stroke center and extent along each of the 4 outer walls (stroke thickness ~14 px)
    left_c = cc[cc <= c_min_px + 20]
    right_c = cc[cc >= c_max_px - 20]
    top_r = rr[rr <= r_top_px + 20]
    bot_r = rr[rr >= r_bot_px - 20]
    c_min = float(np.median(left_c)) + ic0
    c_max = float(np.median(right_c)) + ic0
    r_top_b = float(np.median(top_r)) + ir0
    r_top_lo = float(top_r.min()) + ir0
    r_top_hi = float(top_r.max()) + ir0
    r_bot_b = float(np.median(bot_r)) + ir0
    r_bot_lo = float(bot_r.min()) + ir0
    r_bot_hi = float(bot_r.max()) + ir0

    n_min = 10.0 ** (mx * c_min + bx)
    n_max = 10.0 ** (mx * c_max + bx)
    extra = max(2.0, 2.0 * y_res)
    return [
        _reading_from_row_interval(n_min, c_min, r_bot_b, r_bot_lo, r_bot_hi, my, by, extra_px=extra,
                                   note=f"Fig 5 box ({box_key}) lower-left corner (N_min, alpha_min)"),
        _reading_from_row_interval(n_max, c_max, r_top_b, r_top_lo, r_top_hi, my, by, extra_px=extra,
                                   note=f"Fig 5 box ({box_key}) upper-right corner (N_max, alpha_max)"),
    ]


# ---------------------------------------------------------------------------
# Fig. S1 (p. 11): kappa_B vs alpha for M=8, phi in {1, 2, 4, 8, 20, 40}
# ---------------------------------------------------------------------------

def calibrate_figS1(image: np.ndarray):
    c0, r0, c1, r1 = find_plot_box(image, 140)
    mx = (7.0 - (-3.0)) / (c1 - c0)
    bx = -3.0 - mx * c0
    my = (4.0 - 1.0) / (r0 - r1)
    by = 1.0 - my * r1
    return (c0, r0, c1, r1), (mx, bx, 1.0), (my, by, 1.0)


def read_figS1_curve(image: np.ndarray, at: list[float], phi: int) -> list[dict]:
    """Read one of the 6 colored curves in Fig. S1 (phi in {1, 2, 4, 8, 20, 40})."""
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_figS1(image)
    R, G, B = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    if phi == 1:
        mask = (R < 70) & (G < 70) & (B < 70)
    elif phi == 2:
        mask = (B > 170) & (R < 90) & (G < 90)
    elif phi == 4:
        mask = (R > 170) & (G < 90) & (B < 90)
    elif phi == 8:
        mask = (G > 170) & (R < 90) & (B < 90)
    elif phi == 20:
        mask = (R > 190) & (G > 150) & (B < 100)
    elif phi == 40:
        mask = (R > 170) & (B > 170) & (G < 100)
    else:
        raise ValueError(f"Unknown Fig S1 phi={phi}")
    # Mask out borders, major tick marks (~60 px long), and legend box (top-left)
    mask[:int(r0) + 70, :] = False
    mask[int(r1) - 70:, :] = False
    mask[:, :int(c0) + 20] = False
    mask[:, int(c1) - 20:] = False
    mask[:int(r0) + 1050, :int(c0) + 760] = False

    out = []
    for a in at:
        c_px = (math.log10(a) - bx) / mx
        ci = int(round(c_px))
        found_rr = None
        for dc in (0, -2, 2, -4, 4, -8, 8):
            if 0 <= ci + dc < mask.shape[1]:
                rr = np.where(mask[:, ci + dc])[0]
                if len(rr):
                    found_rr = rr
                    break
        if found_rr is None:
            out.append({"x": float(a), "y": None, "note": f"phi={phi} curve occluded at alpha={a}"})
            continue
        r_mid = float(np.median(found_rr))
        lo_r, hi_r = float(found_rr.min()), float(found_rr.max())
        out.append(_reading_from_row_interval(a, c_px, r_mid, lo_r, hi_r, my, by, extra_px=3.0,
                                              note=f"Fig S1 phi={phi} (M=8)"))
    return out


# ---------------------------------------------------------------------------
# Fig. S2A & S2B (p. 11): fitted c vs phi and vs (phi/phi*)^2 for M in {2,4,8,16,32}
# ---------------------------------------------------------------------------

def calibrate_figS2a(image: np.ndarray):
    c0, r0, c1, r1 = find_plot_box(image, 140)
    # Borders: x in [1, 50], y in [0.1, 5000]
    mx = (math.log10(50.0) - 0.0) / (c1 - c0)
    bx = 0.0 - mx * c0
    my = (math.log10(5000.0) - (-1.0)) / (r0 - r1)
    by = -1.0 - my * r1
    return (c0, r0, c1, r1), (mx, bx, 1.0), (my, by, 1.0)


def calibrate_figS2b(image: np.ndarray):
    c0, r0, c1, r1 = find_plot_box(image, 140)
    # Borders: x in [1e-2, 1e5], y in [0.1, 1e4]
    mx = (5.0 - (-2.0)) / (c1 - c0)
    bx = -2.0 - mx * c0
    my = (4.0 - (-1.0)) / (r0 - r1)
    by = -1.0 - my * r1
    return (c0, r0, c1, r1), (mx, bx, 1.0), (my, by, 1.0)


def _figS2_color_mask(image: np.ndarray, m_val: int, c0: float, r0: float, c1: float, r1: float) -> np.ndarray:
    """Fast uint8 color mask for one M-series in Fig. S2A/S2B with borders and legend masked."""
    R, G, B = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    if m_val == 2:
        mask = (R < 80) & (G < 80) & (B < 80)
    elif m_val == 4:
        mask = (R > 180) & (G < 110) & (B < 140)
    elif m_val == 8:
        mask = (B > 180) & (R < 100) & (G < 100)
    elif m_val == 16:
        mask = (G > 180) & (R < 100) & (B < 100)
    elif m_val == 32:
        mask = (G > 180) & (B > 180) & (R < 100)
    else:
        raise ValueError(f"Unknown Fig S2 M={m_val}")
    mask[:int(r0) + 25, :] = False
    mask[int(r1) - 25:, :] = False
    mask[:, :int(c0) + 25] = False
    mask[:, int(c1) - 25:] = False
    # Exclude legend box in top-left (covers rows r0..r0+880, cols c0..c0+740)
    mask[:int(r0) + 890, :int(c0) + 750] = False
    return mask


def _read_figS2_markers(mask: np.ndarray, at: list[float], mx: float, bx: float, my: float, by: float, label_prefix: str) -> list[dict]:
    """Read centered hollow markers in Fig. S2A/S2B without splitting hollow interiors (~50-55 px)."""
    out = []
    for x in at:
        c_px = (math.log10(x) - bx) / mx
        ci = int(round(c_px))
        # Narrow column slice centered on marker x avoids diagonal line tails while capturing top/bottom rims;
        # gap=70 bridges hollow marker interiors (~50-55 px) so top and bottom edges stay one single marker.
        sub = mask[:, max(0, ci - 12):min(mask.shape[1], ci + 13)]
        rr, _ = np.nonzero(sub)
        if len(rr) == 0:
            sub = mask[:, max(0, ci - 22):min(mask.shape[1], ci + 23)]
            rr, _ = np.nonzero(sub)
        if len(rr) == 0:
            out.append({"x": float(x), "y": None, "note": f"{label_prefix} marker not found at x={x}"})
            continue
        segs = _clusters(np.unique(rr), gap=70)
        best = max(segs, key=lambda s: s[2] - s[1])
        r_mid = 0.5 * (best[1] + best[2])
        half_unc = max(5.0, 0.35 * (best[2] - best[1]))
        out.append(_reading_from_row_interval(x, c_px, r_mid, r_mid - half_unc, r_mid + half_unc, my, by, extra_px=2.0,
                                              note=f"{label_prefix} at x={x}"))
    return out


def read_figS2a_series(image: np.ndarray, at: list[float], m_val: int) -> list[dict]:
    """Read markers for M in {2, 4, 8, 16, 32} at phi in {2, 4, 8, 20, 40} in Fig. S2A."""
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_figS2a(image)
    mask = _figS2_color_mask(image, m_val, c0, r0, c1, r1)
    return _read_figS2_markers(mask, at, mx, bx, my, by, f"Fig S2A M={m_val}")


def read_figS2b_series(image: np.ndarray, at: list[float], m_val: int) -> list[dict]:
    """Read renormalized markers for M in {2, 4, 8, 16, 32} at (phi/phi*)^2 (with phi*=M/8) in Fig. S2B."""
    (c0, r0, c1, r1), (mx, bx, x_res), (my, by, y_res) = calibrate_figS2b(image)
    mask = _figS2_color_mask(image, m_val, c0, r0, c1, r1)
    return _read_figS2_markers(mask, at, mx, bx, my, by, f"Fig S2B M={m_val}")
