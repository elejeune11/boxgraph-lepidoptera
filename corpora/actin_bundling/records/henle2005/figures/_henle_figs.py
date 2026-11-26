"""Calibrated figure-reading helpers for Henle & Pincus (2005), Figs. 1 and 2 (p. 3).

All axis calibrations are derived directly from the rendered image's frame spines and
inward tick marks, matched to the printed tick labels on each axis, and checked for
linear fit residual before reading any curve.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import find_objects, label


def _cluster_1d(indices: np.ndarray, gap: int = 8) -> list[tuple[float, int, int]]:
    """Group sorted 1D pixel indices into contiguous clusters -> (center, start, end)."""
    if len(indices) == 0:
        return []
    groups: list[list[int]] = [[int(indices[0])]]
    for idx in indices[1:]:
        v = int(idx)
        if v - groups[-1][-1] > gap:
            groups.append([v])
        else:
            groups[-1].append(v)
    return [(float(np.mean(g)) + 0.5, g[0], g[-1]) for g in groups]


def _find_frame(dark: np.ndarray, left_rng, right_rng, top_rng, bot_rng, min_v_frac=0.5, min_h_frac=0.5):
    """Locate the four plot box borders inside search ranges."""
    H, W = dark.shape
    col_sums = dark.sum(axis=0)
    row_sums = dark.sum(axis=1)

    l_cands = np.where(col_sums[left_rng[0]:left_rng[1]] > H * min_v_frac)[0] + left_rng[0]
    r_cands = np.where(col_sums[right_rng[0]:right_rng[1]] > H * min_v_frac)[0] + right_rng[0]
    t_cands = np.where(row_sums[top_rng[0]:top_rng[1]] > W * min_h_frac)[0] + top_rng[0]
    b_cands = np.where(row_sums[bot_rng[0]:bot_rng[1]] > W * min_h_frac)[0] + bot_rng[0]
    if len(l_cands) == 0 or len(r_cands) == 0 or len(t_cands) == 0 or len(b_cands) == 0:
        raise RuntimeError("Could not detect all four frame borders")

    l_cl = _cluster_1d(l_cands, gap=10)[0]
    r_cl = _cluster_1d(r_cands, gap=10)[-1]
    t_cl = _cluster_1d(t_cands, gap=10)[0]
    b_cl = _cluster_1d(b_cands, gap=10)[-1]
    # Return inner edges (l1, r0, t1, b0) and centers (lc, rc, tc, bc)
    return {
        "l1": l_cl[2] + 1, "r0": r_cl[1] - 1, "t1": t_cl[2] + 1, "b0": b_cl[1] - 1,
        "lc": l_cl[0], "rc": r_cl[0], "tc": t_cl[0], "bc": b_cl[0],
    }


def _match_regular_ticks(cands: list[float], n: int, expected_step_range: tuple[float, float], tol: float = 12.0) -> list[float]:
    """Find an equally spaced sequence of `n` ticks among `cands`."""
    best = None
    best_err = 1e9
    cs = sorted(cands)
    for i in range(len(cs)):
        for j in range(i + 1, len(cs)):
            step = (cs[j] - cs[i]) / (n - 1)
            if not (expected_step_range[0] <= step <= expected_step_range[1]):
                continue
            picked = []
            errs = []
            ok = True
            for k in range(n):
                target = cs[i] + k * step
                nearest = min(cs, key=lambda v: abs(v - target))
                if abs(nearest - target) > tol:
                    ok = False
                    break
                picked.append(nearest)
                errs.append(abs(nearest - target))
            if ok:
                rms = float(np.sqrt(np.mean(np.square(errs))))
                if rms < best_err:
                    best_err = rms
                    best = picked
    if best is None:
        raise RuntimeError(f"Could not match {n} regularly spaced ticks in {cs}")
    return best


def _fit_axis(pxs: list[float], vals: list[float], max_resid_px: float = 5.0):
    """Fit linear pixel -> value calibration and return (slope, intercept, max_resid_val, px_per_unit)."""
    p = np.asarray(pxs, dtype=float)
    v = np.asarray(vals, dtype=float)
    slope, intercept = np.polyfit(p, v, 1)
    resid_val = np.abs(slope * p + intercept - v)
    px_per_unit = abs(1.0 / slope)
    max_resid_px_actual = float(np.max(resid_val) * px_per_unit)
    if max_resid_px_actual > max_resid_px:
        raise RuntimeError(f"Axis calibration residual {max_resid_px_actual:.2f} px exceeds {max_resid_px} px")
    return float(slope), float(intercept), float(np.max(resid_val)), float(px_per_unit)


def _pix_to_val(cal, px: float) -> float:
    return cal[0] * px + cal[1]


def _val_to_pix(cal, val: float) -> float:
    return (val - cal[1]) / cal[0]


def _make_reading(x_val, y_val, half_err_y, col_px, row_px, half_err_px, note=""):
    return {
        "x": float(x_val),
        "y": float(y_val),
        "lo": float(y_val - half_err_y),
        "hi": float(y_val + half_err_y),
        "px": [float(col_px), float(row_px)],
        "px_lo": float(row_px - half_err_px),
        "px_hi": float(row_px + half_err_px),
        "note": note,
    }


# ---------------------------------------------------------------------------
# Panel 1: Fig. 1 main (1200 dpi, crop [1200, 880, 4560, 2680])
# ---------------------------------------------------------------------------

def _calibrate_fig1_main(image: np.ndarray):
    dark = (image.sum(axis=2, dtype=np.uint16) < 420) if image.dtype == np.uint8 else (image.astype(float).mean(axis=2) < 140.0)
    fr = _find_frame(dark, (400, 750), (2900, 3250), (20, 200), (1500, 1750))
    l1, r0, t1, b0 = fr["l1"], fr["r0"], fr["t1"], fr["b0"]

    # Left spine ticks: 13 ticks from 2.25 down to -0.75 every 0.25 (major ticks: 2, 1.5, 1, 0.5, 0, -0.5)
    yt_raw = np.where(dark[t1 + 10 : b0 - 10, l1 + 2 : l1 + 10].mean(axis=1) > 0.5)[0] + t1 + 10
    yt_cl = [c[0] for c in _cluster_1d(yt_raw, gap=8)]
    yt_all = _match_regular_ticks(yt_cl, 13, (95.0, 125.0), tol=10.0)
    y_vals = [2.25 - 0.25 * k for k in range(13)]
    y_cal = _fit_axis(yt_all, y_vals, max_resid_px=4.0)

    # Bottom spine major ticks: 2, 4, 6, 8 (and minor ticks every 0.5)
    xt_raw = np.where(dark[b0 - 9 : b0 - 2, l1 + 10 : r0 - 10].mean(axis=0) > 0.5)[0] + l1 + 10
    xt_cl = [c[0] for c in _cluster_1d(xt_raw, gap=8)]
    # Top spine ticks are cleaner (solid curve touches bottom near x=2)
    xt_top_raw = np.where(dark[t1 + 2 : t1 + 10, l1 + 10 : l1 + 1350].mean(axis=0) > 0.5)[0] + l1 + 10
    xt_top_cl = [c[0] for c in _cluster_1d(xt_top_raw, gap=8)]
    xt_cands = sorted(set(xt_cl + xt_top_cl))
    # Match 4 major ticks 2, 4, 6, 8 (step ~556 px)
    xt_maj = _match_regular_ticks(xt_cands, 4, (530.0, 585.0), tol=10.0)
    x_cal = _fit_axis(xt_maj, [2.0, 4.0, 6.0, 8.0], max_resid_px=4.0)

    # Detect inset frame inside plot to mask it out
    col_sums_top = dark[t1 + 5 : t1 + 780, l1 + 1000 : r0 - 50].sum(axis=0)
    inset_left = l1 + 1000 + int(np.argmax(col_sums_top))
    row_sums_right = dark[t1 + 400 : t1 + 850, inset_left : r0 - 5].sum(axis=1)
    inset_bot = t1 + 400 + int(np.argmax(row_sums_right))

    mask = np.zeros_like(dark, dtype=bool)
    mask[t1 + 14 : b0 - 14, l1 + 14 : r0 - 14] = dark[t1 + 14 : b0 - 14, l1 + 14 : r0 - 14]
    # Mask out inset box and its axis labels (left to inset_left - 340, down to inset_bot + 92)
    mask[: inset_bot + 92, inset_left - 340 :] = False
    return x_cal, y_cal, fr, mask


def _extract_three_curves_components(
    mask: np.ndarray,
    x_cal,
    y_cal,
    solid_min_area: int,
    dash_area_range,
    dash_min_maxdim: int,
    dot_area_range,
    dot_max_dim: int,
):
    """Separate solid (bottom), dashed (middle), and dotted (top) curves via connected components."""
    lbl, _ = label(mask)
    solid_mask = np.zeros_like(mask, dtype=bool)
    dashes = []
    dots = []

    for i, sl in enumerate(find_objects(lbl), start=1):
        if sl is None:
            continue
        r0, c0 = sl[0].start, sl[1].start
        sub = lbl[sl] == i
        rr_loc, cc_loc = np.where(sub)
        area = len(rr_loc)
        h = int(sl[0].stop - r0)
        w = int(sl[1].stop - c0)
        maxdim = max(w, h)
        if area >= solid_min_area:
            solid_mask[sl] |= sub
        elif dash_area_range[0] <= area <= dash_area_range[1] and maxdim >= dash_min_maxdim:
            # Evaluate at the exact midpoint column of the dash to avoid slope centroid bias
            rr = rr_loc + r0
            cc = cc_loc + c0
            c_mid = int(round(0.5 * (cc.min() + cc.max())))
            rs_mid = rr[cc == c_mid]
            if len(rs_mid) == 0:
                c_mid = int(cc[np.argmin(np.abs(cc - 0.5 * (cc.min() + cc.max())))])
                rs_mid = rr[cc == c_mid]
            r_center = float(np.mean(rs_mid)) + 0.5
            r_thick = float(rs_mid.max() - rs_mid.min() + 1)
            dashes.append((float(c_mid) + 0.5, r_center, r_thick, area))
        elif dot_area_range[0] <= area <= dot_area_range[1] and maxdim <= dot_max_dim:
            c_center = float(np.mean(cc_loc + c0)) + 0.5
            r_center = float(np.mean(rr_loc + r0)) + 0.5
            r_thick = float(h)
            dots.append((c_center, r_center, r_thick, area))

    dashes.sort(key=lambda t: t[0])
    dots.sort(key=lambda t: t[0])
    return solid_mask, dashes, dots


def _sample_solid(solid_mask: np.ndarray, at: list[float], x_cal, y_cal, note_prefix: str):
    out = []
    cal_err_y = y_cal[2]
    px_per_y = y_cal[3]
    for xv in at:
        c = int(round(_val_to_pix(x_cal, float(xv))))
        rows = np.where(solid_mask[:, c])[0]
        if len(rows) == 0:
            out.append({"x": float(xv), "y": None, "note": f"{note_prefix}: curve not present at x={xv}"})
            continue
        r_mid = float(np.mean(rows)) + 0.5
        thick = float(rows.max() - rows.min() + 1)
        # Slope estimate across +/- 5 px
        r_l = np.where(solid_mask[:, max(0, c - 5)])[0]
        r_r = np.where(solid_mask[:, min(solid_mask.shape[1] - 1, c + 5)])[0]
        slope_px = abs((float(np.mean(r_r)) - float(np.mean(r_l))) / 10.0) if len(r_l) and len(r_r) else 0.0
        half_px = 0.5 * thick + 1.0 + 0.5 * slope_px + cal_err_y * px_per_y
        half_y = half_px / px_per_y
        yv = _pix_to_val(y_cal, r_mid)
        out.append(_make_reading(xv, yv, half_y, c + 0.5, r_mid, half_px, f"{note_prefix} solid curve"))
    return out


def _format_components(comps, x_cal, y_cal, note_prefix: str, x_min=None, x_max=None):
    out = []
    cal_err_y = y_cal[2]
    px_per_y = y_cal[3]
    for cx, cy, thick, _ in comps:
        xv = _pix_to_val(x_cal, cx)
        if x_min is not None and xv < x_min:
            continue
        if x_max is not None and xv > x_max:
            continue
        yv = _pix_to_val(y_cal, cy)
        half_px = 0.5 * thick + 1.5 + cal_err_y * px_per_y
        half_y = half_px / px_per_y
        out.append(_make_reading(xv, yv, half_y, cx, cy, half_px, note_prefix))
    return out


def read_fig1_main(image: np.ndarray, at: list[float], curve: str):
    x_cal, y_cal, fr, mask = _calibrate_fig1_main(image)
    solid_mask, dashes, dots = _extract_three_curves_components(
        mask, x_cal, y_cal,
        solid_min_area=5000,
        dash_area_range=(550, 1500),
        dash_min_maxdim=44,
        dot_area_range=(420, 640),
        dot_max_dim=34,
    )
    if curve == "0.1uM":
        return _sample_solid(solid_mask, at, x_cal, y_cal, "Fig. 1 n_q=0.1 uM")
    if curve == "10uM":
        return _format_components(dashes, x_cal, y_cal, "Fig. 1 n_q=10 uM dashed", x_min=1.15, x_max=9.85)
    if curve == "1mM":
        return _format_components(dots, x_cal, y_cal, "Fig. 1 n_q=1 mM dotted", x_min=1.15, x_max=9.85)
    raise ValueError(curve)


# ---------------------------------------------------------------------------
# Panel 2: Fig. 1 inset (2400 dpi, crop [5600, 1880, 8960, 3680])
# ---------------------------------------------------------------------------

def _calibrate_fig1_inset(image: np.ndarray):
    dark = (image.sum(axis=2, dtype=np.uint16) < 420) if image.dtype == np.uint8 else (image.astype(float).mean(axis=2) < 140.0)
    fr = _find_frame(dark, (600, 900), (2800, 3150), (10, 150), (1300, 1600), min_v_frac=0.5, min_h_frac=0.5)
    l1, r0, t1, b0 = fr["l1"], fr["r0"], fr["t1"], fr["b0"]

    yt_raw = np.where(dark[t1 + 40 : b0 - 40, l1 + 2 : l1 + 12].mean(axis=1) > 0.5)[0] + t1 + 40
    yt_cl = [c[0] for c in _cluster_1d(yt_raw, gap=10)]
    yt_6 = _match_regular_ticks(yt_cl, 6, (180.0, 215.0), tol=12.0)
    y_cal = _fit_axis(yt_6, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0], max_resid_px=5.0)

    xt_raw = np.where(dark[b0 - 12 : b0 - 2, l1 + 20 : r0 - 20].mean(axis=0) > 0.5)[0] + l1 + 20
    xt_cl = [c[0] for c in _cluster_1d(xt_raw, gap=10)]
    xt_9 = _match_regular_ticks(xt_cl, 9, (210.0, 240.0), tol=12.0)
    x_cal = _fit_axis(xt_9, [10.0 * k for k in range(1, 10)], max_resid_px=5.0)

    r_asymp_max = int(round(_val_to_pix(y_cal, 7.85)))
    return x_cal, y_cal, fr, dark, r_asymp_max


def read_fig1_inset(image: np.ndarray, at: list[float], curve: str):
    x_cal, y_cal, fr, dark, r_asymp_max = _calibrate_fig1_inset(image)
    l1, r0, t1, b0 = fr["l1"], fr["r0"], fr["t1"], fr["b0"]
    if curve == "lam_tot":
        solid_mask = np.zeros_like(dark, dtype=bool)
        solid_mask[r_asymp_max + 5 : b0 - 15, l1 + 15 : r0 - 15] = dark[r_asymp_max + 5 : b0 - 15, l1 + 15 : r0 - 15]
        return _sample_solid(solid_mask, at, x_cal, y_cal, "Fig. 1 inset lambda_tot(R)")
    if curve == "lam_tot_asymp":
        asymp_mask = np.zeros_like(dark, dtype=bool)
        asymp_mask[t1 + 15 : r_asymp_max, l1 + 20 : r0 - 20] = dark[t1 + 15 : r_asymp_max, l1 + 20 : r0 - 20]
        lbl, _ = label(asymp_mask)
        dashes = []
        for i, sl in enumerate(find_objects(lbl), start=1):
            if sl is None:
                continue
            r0, c0 = sl[0].start, sl[1].start
            rr_loc, cc_loc = np.where(lbl[sl] == i)
            if len(rr_loc) >= 500:
                c_mid = float(np.mean(cc_loc + c0)) + 0.5
                r_mid = float(np.mean(rr_loc + r0)) + 0.5
                thick = float(sl[0].stop - r0)
                dashes.append((c_mid, r_mid, thick, len(rr_loc)))
        dashes.sort(key=lambda t: t[0])
        return _format_components(dashes, x_cal, y_cal, "Fig. 1 inset asymptotic lambda_tot dashed line", x_min=5.0, x_max=95.0)
    raise ValueError(curve)


# ---------------------------------------------------------------------------
# Panel 3: Fig. 2 main (1200 dpi, crop [5600, 880, 8960, 2680])
# ---------------------------------------------------------------------------

def _calibrate_fig2_main(image: np.ndarray):
    dark = (image.sum(axis=2, dtype=np.uint16) < 420) if image.dtype == np.uint8 else (image.astype(float).mean(axis=2) < 140.0)
    fr = _find_frame(dark, (400, 750), (2900, 3250), (20, 200), (1500, 1750))
    l1, r0, t1, b0 = fr["l1"], fr["r0"], fr["t1"], fr["b0"]

    yt_raw = np.where(dark[t1 + 10 : b0 - 10, l1 + 2 : l1 + 10].mean(axis=1) > 0.5)[0] + t1 + 10
    yt_cl = [c[0] for c in _cluster_1d(yt_raw, gap=8)]
    yt_9 = _match_regular_ticks(yt_cl, 9, (140.0, 168.0), tol=10.0)
    y_vals = [-5.1 - 0.1 * k for k in range(9)]
    y_cal = _fit_axis(yt_9, y_vals, max_resid_px=4.0)

    xt_raw = np.where(dark[b0 - 9 : b0 - 2, l1 + 10 : r0 - 10].mean(axis=0) > 0.5)[0] + l1 + 10
    xt_cl = [c[0] for c in _cluster_1d(xt_raw, gap=8)]
    xt_5 = _match_regular_ticks(xt_cl, 5, (410.0, 450.0), tol=10.0)
    x_cal = _fit_axis(xt_5, [5.0, 10.0, 15.0, 20.0, 25.0], max_resid_px=4.0)

    col_sums_top = dark[t1 + 5 : t1 + 780, l1 + 1000 : r0 - 50].sum(axis=0)
    inset_left = l1 + 1000 + int(np.argmax(col_sums_top))
    row_sums_right = dark[t1 + 400 : t1 + 850, inset_left : r0 - 5].sum(axis=1)
    inset_bot = t1 + 400 + int(np.argmax(row_sums_right))

    mask = np.zeros_like(dark, dtype=bool)
    mask[t1 + 14 : b0 - 14, l1 + 14 : r0 - 14] = dark[t1 + 14 : b0 - 14, l1 + 14 : r0 - 14]
    mask[: inset_bot + 105, inset_left - 340 :] = False
    return x_cal, y_cal, fr, mask


def read_fig2_main(image: np.ndarray, at: list[float], curve: str):
    x_cal, y_cal, fr, mask = _calibrate_fig2_main(image)
    solid_mask, dashes, dots = _extract_three_curves_components(
        mask, x_cal, y_cal,
        solid_min_area=5000,
        dash_area_range=(620, 1500),
        dash_min_maxdim=44,
        dot_area_range=(420, 640),
        dot_max_dim=34,
    )
    if curve == "0.1uM":
        return _sample_solid(solid_mask, at, x_cal, y_cal, "Fig. 2 n_q=0.1 uM")
    if curve == "10uM":
        return _format_components(dashes, x_cal, y_cal, "Fig. 2 n_q=10 uM dashed", x_min=1.5, x_max=29.2)
    if curve == "1mM":
        return _format_components(dots, x_cal, y_cal, "Fig. 2 n_q=1 mM dotted", x_min=3.0, x_max=29.2)
    raise ValueError(curve)


# ---------------------------------------------------------------------------
# Panel 4: Fig. 2 inset (2400 dpi, crop [14320, 1880, 17760, 3680])
# ---------------------------------------------------------------------------

def _calibrate_fig2_inset(image: np.ndarray):
    dark = (image.sum(axis=2, dtype=np.uint16) < 420) if image.dtype == np.uint8 else (image.astype(float).mean(axis=2) < 140.0)
    fr = _find_frame(dark, (600, 1000), (2800, 3200), (10, 150), (1300, 1600), min_v_frac=0.5, min_h_frac=0.5)
    l1, r0, t1, b0 = fr["l1"], fr["r0"], fr["t1"], fr["b0"]

    yt_raw = np.where(dark[t1 + 20 : b0 - 20, l1 + 2 : l1 + 12].mean(axis=1) > 0.5)[0] + t1 + 20
    yt_cl = [c[0] for c in _cluster_1d(yt_raw, gap=10)]
    yt_10 = _match_regular_ticks(yt_cl, 10, (115.0, 136.0), tol=12.0)
    y_vals = [0.3 - 0.1 * k for k in range(10)]
    y_cal = _fit_axis(yt_10, y_vals, max_resid_px=5.0)

    xt_raw = np.where(dark[t1 + 2 : t1 + 12, l1 + 20 : r0 - 20].mean(axis=0) > 0.5)[0] + l1 + 20
    xt_cl = [c[0] for c in _cluster_1d(xt_raw, gap=10)]
    xt_6 = _match_regular_ticks(xt_cl, 6, (325.0, 365.0), tol=12.0)
    x_cal = _fit_axis(xt_6, [2.0, 4.0, 6.0, 8.0, 10.0, 12.0], max_resid_px=5.0)

    mask = np.zeros_like(dark, dtype=bool)
    mask[t1 + 25 : b0 - 25, l1 + 15 : r0 - 15] = dark[t1 + 25 : b0 - 25, l1 + 15 : r0 - 15]
    return x_cal, y_cal, fr, mask


def read_fig2_inset(image: np.ndarray, at: list[float], curve: str):
    x_cal, y_cal, fr, mask = _calibrate_fig2_inset(image)
    solid_mask, dashes, dots = _extract_three_curves_components(
        mask, x_cal, y_cal,
        solid_min_area=15000,
        dash_area_range=(2500, 14500),
        dash_min_maxdim=90,
        dot_area_range=(400, 2400),
        dot_max_dim=65,
    )
    if curve == "0.1uM":
        return _sample_solid(solid_mask, at, x_cal, y_cal, "Fig. 2 inset n_q=0.1 uM")
    if curve == "10uM":
        # Exclude clipped bottom dash at R/a ~ 3.6
        return _format_components(dashes, x_cal, y_cal, "Fig. 2 inset n_q=10 uM dashed", x_min=4.0, x_max=13.7)
    if curve == "1mM":
        return _format_components(dots, x_cal, y_cal, "Fig. 2 inset n_q=1 mM dotted", x_min=1.1, x_max=13.9)
    raise ValueError(curve)
