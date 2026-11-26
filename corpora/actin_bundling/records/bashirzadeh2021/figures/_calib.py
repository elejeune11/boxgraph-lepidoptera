"""Shared image-based axis calibration and figure reading helpers for bashirzadeh2021.

Every calibration and data reading is derived strictly from the passed RGB `image`
array (no PDF re-rendering inside `read(image, at)` and no hardcoded data values).
All operations use fast vectorized NumPy/SciPy for sub-150ms execution per figure.
"""

from __future__ import annotations

import math
import numpy as np

from boxgraph import digitize as D


# ---------------------------------------------------------------------------
# Fast vectorized luminance, saturation, clustering, and axis calibration
# ---------------------------------------------------------------------------

def _lum_sat(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fast integer-based grayscale (0..255 uint8) and chroma spread (uint8)."""
    r = img[:, :, 0]
    g = img[:, :, 1]
    b = img[:, :, 2]
    gray = ((r.astype(np.uint16) + g.astype(np.uint16) + b.astype(np.uint16)) // 3).astype(np.uint8)
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    sat = (mx - mn).astype(np.uint8)
    return gray, sat


def _cl1d(idx: np.ndarray, gap: int = 4) -> list[tuple[int, int, float]]:
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    parts = np.split(idx, splits)
    return [(int(p[0]), int(p[-1]), 0.5 * (float(p[0]) + float(p[-1]))) for p in parts]


def _fit_linear_axis(ticks_px: list[float], vals: list[float]) -> tuple[D.Axis, float]:
    ax = D.Axis("linear", list(zip(ticks_px, vals)))
    return ax, float(ax.residual_px)


def _select_regular(cands: list[float], n_expected: int, tol: float = 3.5) -> list[float]:
    """Select the most evenly spaced subset of exact length `n_expected` from candidates."""
    c = sorted(cands)
    if len(c) == n_expected:
        return c
    if len(c) < n_expected:
        raise RuntimeError(f"Expected {n_expected} ticks, only found {len(c)}: {c}")
    best: list[float] | None = None
    best_err = float("inf")
    best_span = -1.0
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < 25.0:
                continue
            pts = []
            max_e = 0.0
            ok = True
            for k in range(n_expected):
                target = c[i] + k * step
                nearest = min(c, key=lambda v: abs(v - target))
                e = abs(nearest - target)
                if e > max(tol, 0.03 * step):
                    ok = False
                    break
                pts.append(nearest)
                if e > max_e:
                    max_e = e
            if ok and len(set(pts)) == n_expected:
                span = pts[-1] - pts[0]
                if (max_e <= tol and (best is None or best_err > tol or span > best_span + 5.0 or (abs(span - best_span) <= 5.0 and max_e < best_err))) or (
                    best is None or (best_err > tol and max_e < best_err)
                ):
                    best = pts
                    best_err = max_e
                    best_span = span
    if best is None:
        raise RuntimeError(f"Could not find {n_expected} evenly spaced ticks among {c}")
    return best


def _find_l_axes(
    gray: np.ndarray,
    sat: np.ndarray,
    *,
    x_r_lo_frac: float = 0.55,
    x_r_hi_frac: float = 0.95,
    y_c_lo_frac: float = 0.10,
    y_c_hi_frac: float = 0.48,
    dark_thr: int = 140,
    sat_thr: int = 40,
) -> tuple[float, float, int, int, int, int]:
    """Find L-shaped bottom x-axis `r_x` and LEFTMOST y-axis spine `c_y`, plus their spans.

    Returns (r_x, c_y, x0, x1, y0, y1).
    """
    dark = (gray < dark_thr) & (sat < sat_thr)
    H, W = dark.shape

    # 1. Bottom x-axis: strongest dark horizontal line in [x_r_lo_frac .. x_r_hi_frac]
    r0, r1 = int(H * x_r_lo_frac), int(H * x_r_hi_frac)
    row_frac = (gray[r0:r1, :] < dark_thr).mean(axis=1)
    r_best = int(np.argmax(row_frac)) + r0
    r_win0, r_win1 = max(0, r_best - 8), min(H, r_best + 9)
    row_run = _cl1d(np.where((gray[r_win0:r_win1, :] < dark_thr).mean(axis=1) > 0.5 * row_frac[r_best - r0])[0], gap=2)
    r_x = float(r_win0 + max(row_run, key=lambda t: t[1] - t[0])[2]) if row_run else float(r_best)
    ri = int(round(r_x))
    x_cols = np.where((gray[max(0, ri - 3) : min(H, ri + 4), :] < dark_thr).any(axis=0))[0]
    x_span = max(_cl1d(x_cols, gap=14), key=lambda t: t[1] - t[0])

    # 2. Left y-axis: LEFTMOST tall vertical dark spine in [y_c_lo_frac .. y_c_hi_frac] above r_x
    c0, c1 = max(int(W * y_c_lo_frac), x_span[0] - 30), min(int(W * y_c_hi_frac), x_span[0] + 150)
    if c0 >= c1:
        c0, c1 = int(W * y_c_lo_frac), int(W * y_c_hi_frac)
    r_lo_y = int(H * 0.04)
    col_frac = dark[r_lo_y : ri + 4, c0:c1].mean(axis=0)
    mx_cf = float(col_frac.max())
    # Take the LEFTMOST cluster with >= 65% of max vertical coverage (avoids picking bar side walls!)
    v_cls = _cl1d(np.where(col_frac >= 0.65 * mx_cf)[0], gap=4)
    c_y = float(c0 + v_cls[0][2])
    ci = int(round(c_y))
    y_rows = np.where(dark[r_lo_y : ri + 6, max(0, ci - 3) : min(W, ci + 4)].any(axis=1))[0] + r_lo_y
    y_span = max(_cl1d(y_rows, gap=15), key=lambda t: t[1] - t[0])
    return float(r_x), float(c_y), int(x_span[0]), int(x_span[1]), int(y_span[0]), int(y_span[1])


def _outward_yticks(
    dark: np.ndarray,
    c_y: float,
    r_lo: int,
    r_hi: int,
    n_expected: int,
    *,
    d0: int = 8,
    d1: int = 22,
    thr: float = 0.45,
) -> list[float]:
    ci = int(round(c_y))
    ra, rb = max(0, r_lo), min(dark.shape[0], r_hi)
    prof = dark[ra:rb, max(0, ci - d1) : max(1, ci - d0)].mean(axis=1)
    cands = [ra + m for _, _, m in _cl1d(np.where(prof > thr)[0], gap=5)]
    return _select_regular(cands, n_expected)


def _inward_yticks(
    dark: np.ndarray,
    c_y: float,
    r_lo: int,
    r_hi: int,
    n_expected: int,
    *,
    d0: int = 6,
    d1: int = 16,
    thr: float = 0.45,
) -> list[float]:
    ci = int(round(c_y))
    ra, rb = max(0, r_lo), min(dark.shape[0], r_hi)
    prof = dark[ra:rb, min(dark.shape[1] - 1, ci + d0) : min(dark.shape[1], ci + d1)].mean(axis=1)
    cands = [ra + m for _, _, m in _cl1d(np.where(prof > thr)[0], gap=5)]
    return _select_regular(cands, n_expected)


def _outward_xticks(
    dark: np.ndarray,
    r_x: float,
    c_lo: int,
    c_hi: int,
    n_expected: int,
    *,
    d0: int = 8,
    d1: int = 22,
    thr: float = 0.45,
) -> list[float]:
    ri = int(round(r_x))
    ca, cb = max(0, c_lo), min(dark.shape[1], c_hi)
    prof = dark[min(dark.shape[0] - 1, ri + d0) : min(dark.shape[0], ri + d1), ca:cb].mean(axis=0)
    cands = [ca + m for _, _, m in _cl1d(np.where(prof > thr)[0], gap=5)]
    return _select_regular(cands, n_expected)


def _inward_xticks(
    dark: np.ndarray,
    r_x: float,
    c_lo: int,
    c_hi: int,
    n_expected: int,
    *,
    d0: int = 5,
    d1: int = 15,
    thr: float = 0.45,
) -> list[float]:
    ri = int(round(r_x))
    ca, cb = max(0, c_lo), min(dark.shape[1], c_hi)
    prof = dark[max(0, ri - d1) : max(1, ri - d0), ca:cb].mean(axis=0)
    cands = [ca + m for _, _, m in _cl1d(np.where(prof > thr)[0], gap=5)]
    return _select_regular(cands, n_expected)


def _make_row(
    x: float,
    row_px: float,
    col_px: float,
    y_ax: D.Axis,
    half_px: float,
    *,
    note: str = "",
    bar_lo_px: float | None = None,
    bar_hi_px: float | None = None,
    bar_half_px: float = 2.0,
) -> dict:
    """Build a BoxGraph code-mode reading dict with pixel uncertainty and optional error bars."""
    y = float(y_ax.to_value(row_px))
    v1 = float(y_ax.to_value(row_px - half_px))
    v2 = float(y_ax.to_value(row_px + half_px))
    lo, hi = (min(v1, v2), max(v1, v2))
    out: dict = {
        "x": float(x),
        "y": y,
        "lo": lo,
        "hi": hi,
        "px": [float(col_px), float(row_px)],
        "px_lo": float(row_px - half_px),
        "px_hi": float(row_px + half_px),
        "note": note,
    }
    if bar_lo_px is not None and bar_hi_px is not None:
        y_hi = float(y_ax.to_value(bar_hi_px))
        y_lo = float(y_ax.to_value(bar_lo_px))
        b_lo, b_hi = (min(y_lo, y_hi), max(y_lo, y_hi))
        dy_bar = abs(float(y_ax.to_value(bar_hi_px - bar_half_px)) - y_hi)
        out["bar_lo"] = b_lo
        out["bar_hi"] = b_hi
        out["bar_lo_outer"] = b_lo - dy_bar
        out["bar_hi_outer"] = b_hi + dy_bar
    return out


# ---------------------------------------------------------------------------
# Fast integer color masks tailored to Bashirzadeh et al. (2021) plots
# ---------------------------------------------------------------------------

def _color_mask(img: np.ndarray, kind: str) -> np.ndarray:
    R = img[:, :, 0].astype(np.int16)
    G = img[:, :, 1].astype(np.int16)
    B = img[:, :, 2].astype(np.int16)
    if kind == "blue":
        return (B > 155) & ((B - R) > 70) & ((B - G) > 70)
    if kind == "magenta":
        return (R > 160) & (B > 160) & (G < 110) & ((R - G) > 70)
    if kind == "olive":
        return (R > 85) & (G > 85) & (B < 65) & (np.abs(R - G) < 35) & (R < 200)
    if kind == "teal":
        return (G > 115) & (B > 115) & (R < 70) & ((G - R) > 65)
    if kind == "maroon":
        return (R > 110) & (R < 225) & (G < 60) & (B < 60) & ((R - G) > 75)
    if kind == "black":
        mx = np.maximum(np.maximum(R, G), B)
        mn = np.minimum(np.minimum(R, G), B)
        return ((R + G + B) < 225) & ((mx - mn) < 35)
    if kind == "red":
        return (R > 140) & ((R - np.maximum(G, B)) > 55)
    if kind == "sfig3_blue":
        return (B > 120) & ((B - R) > 45) & (R < 130)
    if kind == "sfig3_red":
        return (R > 140) & ((R - B) > 70) & (B < 120)
    if kind == "sfig8_green":
        return (G > R + 12) & (G > B + 8) & (G < 165) & (R < 125)
    if kind == "sfig8_magenta":
        return (R > G + 40) & (B > G + 10) & (R > 130)
    if kind == "sfig8e_green":
        return (G > R + 5) & (G > B + 4) & ((R + G + B) < 435)
    if kind == "sfig8e_magenta":
        return (R > G + 22) & (R > B + 5) & (R > 115)
    if kind == "sfig8e_blue":
        return (B > R + 12) & (B > G + 8) & (B > 85) & (R < 115)
    raise ValueError(f"Unknown color kind {kind!r}")


# ---------------------------------------------------------------------------
# Generic curve reader on clean color mask
# ---------------------------------------------------------------------------

def _read_curve_on_mask(
    mask: np.ndarray,
    x_ax: D.Axis,
    y_ax: D.Axis,
    at: list[float],
    r_min: int,
    r_max: int,
    *,
    win_cols: int = 3,
    calib_resid: float = 0.5,
    note_prefix: str = "",
) -> list[dict]:
    H, W = mask.shape
    rows_out: list[dict] = []
    for xv in at:
        col = float(x_ax.to_pixel(xv))
        ci = int(round(col))
        best_r = None
        best_half = None
        for hw in (win_cols, win_cols + 4, win_cols + 9, win_cols + 15):
            c0 = max(0, ci - hw)
            c1 = min(W, ci + hw + 1)
            sub = mask[r_min:r_max, c0:c1]
            rr, _ = np.where(sub)
            if len(rr) > 0:
                rr = rr + r_min
                r_med = float(np.median(rr))
                cls = _cl1d(np.unique(rr), gap=8)
                chosen = min(cls, key=lambda t: abs(t[2] - r_med))
                r_sel = rr[(rr >= chosen[0]) & (rr <= chosen[1])]
                best_r = float(np.mean(r_sel))
                line_half = 0.5 * (float(r_sel.max()) - float(r_sel.min()))
                best_half = math.hypot(max(line_half, 1.5), calib_resid + (0.0 if hw == win_cols else 0.5 * (hw - win_cols)))
                break
        if best_r is None:
            rows_out.append({"x": float(xv), "y": None, "found": False, "note": "no curve pixels found"})
        else:
            rows_out.append(
                _make_row(xv, best_r, col, y_ax, best_half, note=f"{note_prefix}calib_resid={calib_resid:.2f}px")
            )
    return rows_out


# ---------------------------------------------------------------------------
# 1. Stacked grayscale bar charts: Fig. 1b, Fig. 1g, Fig. 2b
# ---------------------------------------------------------------------------

def _read_stacked_gray_bars(
    img: np.ndarray,
    *,
    y_vals: list[float],
    x_vals: list[float],
    boundary: str,
    n_yticks: int,
) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.65, x_r_hi_frac=0.95)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, n_yticks)
    y_ax, y_res = _fit_linear_axis(yticks, y_vals)
    xticks = _outward_xticks(dark, r_x, int(round(c_y)) + 15, x1 + 10, len(x_vals))

    r_top = int(round(yticks[0]))
    r_bot = int(round(yticks[-1]))
    results: list[dict] = []

    for xv, xc in zip(x_vals, xticks):
        ci = int(round(xc))
        strip = gray[r_top - 4 : r_bot - 4, ci - 12 : ci + 13].mean(axis=1)
        rows_abs = np.arange(r_top - 4, r_bot - 4)
        light_mask = strip > 175.0
        med_mask = (strip >= 75.0) & (strip <= 175.0)

        if boundary in ("ring", "periph"):
            light_rows = rows_abs[light_mask]
            if len(light_rows) == 0:
                row_val = float(r_bot)
            else:
                top_of_light = _cl1d(light_rows, gap=4)[-1][0]
                rel = top_of_light - (r_top - 4)
                b_win = strip[max(0, rel - 12) : rel + 1]
                b_hits = np.where(b_win < 75.0)[0]
                if len(b_hits) > 0:
                    row_val = float(max(0, rel - 12) + np.mean(b_hits) + (r_top - 4))
                else:
                    row_val = float(top_of_light)
        elif boundary == "ring_plus_ringnet":
            # Require contiguous medium-grey block (>= 6 px tall) strictly below top spine border (r_top + 4)
            # to ignore 1-2 px anti-aliasing along the top black bar edge.
            valid_med = med_mask & (rows_abs >= r_top + 4)
            med_rows = rows_abs[valid_med]
            med_cls = [c for c in _cl1d(med_rows, gap=3) if (c[1] - c[0] + 1) >= 6]
            if not med_cls:
                # Tiny or absent medium-grey segment: check if any medium-grey exists directly above light segment
                light_rows = rows_abs[light_mask]
                row_val = float(_cl1d(light_rows, gap=4)[-1][0]) if len(light_rows) else float(r_bot)
            else:
                top_of_med = med_cls[0][0]
                row_val = float(top_of_med) - 2.0
        else:
            raise ValueError(boundary)

        half_px = math.hypot(2.5, y_res)
        results.append(_make_row(xv, row_val, xc, y_ax, half_px, note=f"{boundary}; y_resid={y_res:.2f}px"))
    return results


def read_fig1b(img: np.ndarray, boundary: str) -> list[dict]:
    return _read_stacked_gray_bars(
        img, y_vals=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0], x_vals=[1.0, 2.0, 3.0], boundary=boundary, n_yticks=6
    )


def read_fig1g(img: np.ndarray) -> list[dict]:
    return _read_stacked_gray_bars(
        img, y_vals=[100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0], x_vals=[0.1, 0.2, 0.3], boundary="periph", n_yticks=11
    )


def read_fig2b(img: np.ndarray, boundary: str) -> list[dict]:
    return _read_stacked_gray_bars(
        img, y_vals=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0], x_vals=[1.0, 2.0, 3.0], boundary=boundary, n_yticks=6
    )


# ---------------------------------------------------------------------------
# 2. Hollow colored bar charts with error bars + replicate diamonds:
#    Fig. 1c, Fig. 1f, Fig. 2d, Fig. 2f
# ---------------------------------------------------------------------------

def _measure_hollow_bar_and_err(
    mask: np.ndarray,
    xc_est: float,
    bar_w: float,
    r_100: float,
    r_0: float,
    *,
    one_sided: bool = False,
) -> tuple[float, float, float | None, float | None, str]:
    H, W = mask.shape
    r_top_lim = max(0, int(round(r_100)) - 45)
    r_base = int(round(r_0)) - 4
    hw = int(round(0.62 * bar_w))
    c_s = max(0, int(round(xc_est)) - hw)
    c_e = min(W, int(round(xc_est)) + hw + 1)

    near_base = mask[max(0, r_base - 14) : r_base, c_s:c_e]
    base_cols = np.where(near_base.mean(axis=0) > 0.35)[0]
    if len(base_cols) == 0:
        return float(r_0), float(xc_est), None, None, "zero-height bar at baseline"

    wall_cls = _cl1d(base_cols, gap=4)
    cL_rel, cR_rel = wall_cls[0][0], wall_cls[-1][1]
    if (cR_rel - cL_rel + 1) < 0.65 * bar_w:
        return float(r_0), float(xc_est), None, None, "near-zero bar (no full-width walls at baseline)"

    cL = c_s + cL_rel
    cR = c_s + cR_rel
    xc = 0.5 * (cL + cR)

    left_w = mask[r_top_lim:r_base, cL : min(W, cL + 5)].sum(axis=1) >= 3
    right_w = mask[r_top_lim:r_base, max(0, cR - 4) : cR + 1].sum(axis=1) >= 3
    both_w = left_w & right_w
    idx = len(both_w) - 1
    while idx > 0 and (both_w[idx] or (idx >= 1 and both_w[idx - 1]) or (idx >= 2 and both_w[idx - 2])):
        idx -= 1
    r_wall_top = r_top_lim + idx + 1

    ra = max(r_top_lim, r_wall_top - 4)
    rb = min(r_base, r_wall_top + 18)
    horiz_frac = mask[ra:rb, cL + 2 : max(cL + 3, cR - 1)].mean(axis=1)
    h_hits = np.where(horiz_frac > 0.85)[0]
    r_bar = float(ra + _cl1d(h_hits, gap=2)[0][2]) if len(h_hits) > 0 else float(r_wall_top)

    ci_chk = int(round(xc))
    r0_i = int(round(r_0))
    # Reject a diamond centered at y = 0% (symmetrical colored ink below the black x-axis with no bar extending above the diamond)
    below_axis_ink = mask[min(H - 1, r0_i + 5) : min(H, r0_i + 15), max(0, ci_chk - 6) : min(W, ci_chk + 7)].mean()
    above_diam_ink = mask[max(0, r0_i - 26) : max(1, r0_i - 18), max(0, ci_chk - 6) : min(W, ci_chk + 7)].mean()
    if r_0 - r_bar < 4.0 or (r_0 - r_bar < 12.0 and below_axis_ink > 0.45 and above_diam_ink < 0.15):
        return float(r_0), float(xc), None, None, "zero-height bar at baseline"

    ci = int(round(xc))
    if not one_sided:
        r_s0 = max(r_top_lim, int(round(r_100)) - 25)
        r_s1 = r_base - 2
        cap_band = mask[r_s0:r_s1, max(0, ci - 10) : min(W, ci + 11)].mean(axis=1)
        cand_rows = np.where(cap_band > 0.85)[0] + r_s0
        strict_above, strict_below = [], []
        loose_above, loose_below = [], []
        for r in cand_rows:
            if abs(r - r_bar) < 8 or r - 7 < 0 or r + 7 >= H:
                continue
            air_ul = mask[r - 7, max(0, ci - 10) : max(1, ci - 6)].mean() < 0.25
            air_dl = mask[r + 7, max(0, ci - 10) : max(1, ci - 6)].mean() < 0.25
            air_ur = mask[r - 7, min(W - 1, ci + 7) : min(W, ci + 11)].mean() < 0.25
            air_dr = mask[r + 7, min(W - 1, ci + 7) : min(W, ci + 11)].mean() < 0.25
            if air_ul and air_dl and air_ur and air_dr:
                if r < r_bar - 7:
                    strict_above.append(int(r))
                elif r > r_bar + 7:
                    strict_below.append(int(r))
            elif r < r_bar - 7 and air_dl and air_dr:
                w_cap = mask[r, max(0, ci - 22) : min(W, ci + 23)].sum()
                if w_cap <= 34:
                    loose_above.append(int(r))
            elif r > r_bar + 7 and air_ul and air_ur:
                w_cap = mask[r, max(0, ci - 15) : min(W, ci + 16)].sum()
                if w_cap <= 31:
                    loose_below.append(int(r))
        above = [m for _, _, m in _cl1d(np.array(strict_above, dtype=int), gap=3)]
        below = [m for _, _, m in _cl1d(np.array(strict_below, dtype=int), gap=3)]
        if not above and not below:
            above = [m for _, _, m in _cl1d(np.array(loose_above, dtype=int), gap=3)]
            below = [m for _, _, m in _cl1d(np.array(loose_below, dtype=int), gap=3)]
        eb_hi = above[-1] if above else None
        eb_lo = below[0] if below else None
        note = ""
        if eb_hi is not None and eb_lo is None:
            eb_lo = min(float(r_0), r_bar + (r_bar - eb_hi))
            note = "lower SEM cap occluded by diamond; mirrored from upper cap"
        elif eb_lo is not None and eb_hi is None:
            eb_hi = r_bar - (eb_lo - r_bar)
            note = "upper SEM cap occluded by diamond; mirrored from lower cap"
        elif eb_hi is not None and eb_lo is not None:
            d_up, d_dn = r_bar - eb_hi, eb_lo - r_bar
            if abs(d_up - d_dn) > max(8.0, 0.25 * max(d_up, d_dn)):
                stem_up = mask[int(round(eb_hi)) : int(round(r_bar)), ci - 3 : ci + 4].any(axis=1).mean()
                stem_dn = mask[int(round(r_bar)) : int(round(eb_lo)), ci - 3 : ci + 4].any(axis=1).mean()
                if stem_up >= stem_dn:
                    eb_lo = min(float(r_0), r_bar + d_up)
                else:
                    eb_hi = r_bar - d_dn
                note = "reconciled asymmetric cap against central stem"
        return r_bar, float(xc), eb_hi, eb_lo, note
    else:
        stem = mask[:, max(0, ci - 3) : min(W, ci + 4)].any(axis=1)
        r_cur = int(round(r_bar)) - 4
        while r_cur > r_top_lim and stem[r_cur]:
            r_cur -= 1
        if r_bar - r_cur < 6:
            return r_bar, float(xc), None, None, "upper error whisker hidden inside bar/diamond"
        blob_w = mask[r_cur : min(H, r_cur + 40), max(0, ci - 24) : min(W, ci + 25)].sum(axis=1)
        if blob_w.max() > 26:
            eb_hi = float(r_cur + int(np.argmax(blob_w)))
            note = "one-sided upward whisker terminates at top diamond (~mean+SD; lower mirrored)"
        else:
            eb_hi = float(r_cur + 1)
            note = "one-sided upward whisker (lower mirrored)"
        eb_lo = min(float(r_0), r_bar + (r_bar - eb_hi))
        return r_bar, float(xc), eb_hi, eb_lo, note


def _read_grouped_bar_panel(
    img: np.ndarray,
    *,
    ratio: str,
    n_groups: int,
    n_yticks: int,
    one_sided: bool = False,
) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.65, x_r_hi_frac=0.92)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, n_yticks)
    y_vals = [100.0 - k * (100.0 / (n_yticks - 1)) for k in range(n_yticks)]
    y_ax, y_res = _fit_linear_axis(yticks, y_vals)
    r_100, r_0 = yticks[0], yticks[-1]

    m_b = _color_mask(img, "blue")
    m_m = _color_mask(img, "magenta")
    m_o = _color_mask(img, "olive")
    r_lo = max(0, int(round(r_100)) - 35)
    r_hi = int(round(r_0)) - 3
    c_lo = int(round(c_y)) + 15
    c_hi = img.shape[1] - 15

    union = (m_b | m_m | m_o)[r_lo:r_hi, c_lo:c_hi]
    col_has = np.where(union.sum(axis=0) > 8)[0] + c_lo
    sub_cls = _cl1d(col_has, gap=12)
    grp_cls = _cl1d(col_has, gap=42)
    if len(grp_cls) != n_groups:
        g0, g1 = float(col_has[0]), float(col_has[-1])
        gw = (g1 - g0) / n_groups
        grp_centers = [g0 + (k + 0.5) * gw for k in range(n_groups)]
    else:
        grp_centers = [m for _, _, m in grp_cls]
    widths = [b - a + 1 for a, b, _ in sub_cls if 35 <= (b - a + 1) <= 85]
    bar_w = float(np.median(widths)) if widths else 52.0

    sub_offsets = []
    for gc in grp_centers:
        win = int(round(2.1 * bar_w))
        for mm in (m_b, m_o):
            sub = mm[int(round(r_100)) : r_hi, max(0, int(gc) - win) : min(img.shape[1], int(gc) + win)]
            cols = np.where(sub.sum(axis=0) > 15)[0]
            if len(cols) > 0:
                sub_offsets.append(abs(max(0, int(gc) - win) + 0.5 * (cols[0] + cols[-1]) - gc))
    sub_step = float(np.median(sub_offsets)) if sub_offsets else 1.25 * bar_w

    masks_map = {"0.1": (m_b, -1), "0.2": (m_m, 0), "0.3": (m_o, +1)}
    mask, slot = masks_map[ratio]

    out: list[dict] = []
    for gi, gc in enumerate(grp_centers):
        xc_est = gc + slot * sub_step
        r_bar, xc, eb_hi, eb_lo, note = _measure_hollow_bar_and_err(
            mask, xc_est, bar_w, r_100, r_0, one_sided=one_sided
        )
        half_px = math.hypot(2.5, y_res)
        row = _make_row(
            float(gi + 1),
            r_bar,
            xc,
            y_ax,
            half_px,
            note=f"ratio={ratio}; {note}; y_resid={y_res:.2f}px".strip("; "),
            bar_lo_px=eb_lo,
            bar_hi_px=eb_hi,
            bar_half_px=3.0,
        )
        if eb_hi is None:
            row["bar_lo"] = row["y"]
            row["bar_hi"] = row["y"]
            row["bar_lo_outer"] = row["lo"]
            row["bar_hi_outer"] = row["hi"]
        out.append(row)
    return out


def read_fig1c(img: np.ndarray, ratio: str) -> list[dict]:
    return _read_grouped_bar_panel(img, ratio=ratio, n_groups=3, n_yticks=11, one_sided=False)


def read_fig1f(img: np.ndarray, ratio: str) -> list[dict]:
    return _read_grouped_bar_panel(img, ratio=ratio, n_groups=2, n_yticks=6, one_sided=True)


def read_fig2d(img: np.ndarray, ratio: str) -> list[dict]:
    return _read_grouped_bar_panel(img, ratio=ratio, n_groups=3, n_yticks=11, one_sided=False)


def read_fig2f(img: np.ndarray, ratio: str) -> list[dict]:
    return _read_grouped_bar_panel(img, ratio=ratio, n_groups=2, n_yticks=6, one_sided=True)


# ---------------------------------------------------------------------------
# 3. Stacked hollow colored bar charts (Network bottom + Ring/Network top):
#    Fig. 1d, Fig. 2e
# ---------------------------------------------------------------------------

def _read_stacked_hollow_panel(img: np.ndarray, *, ratio: str, layer: str) -> list[dict]:
    """Read Fig. 1d or Fig. 2e stacked hollow bars using exact outer-wall pairs."""
    gray, sat = _lum_sat(img)
    dark = (gray < 140) & (sat < 40)
    H, W = dark.shape
    c_y_col = int(np.argmax(dark[int(H * 0.44) : int(H * 0.88), int(W * 0.15) : int(W * 0.42)].mean(axis=0))) + int(W * 0.15)
    c_y = float(c_y_col)
    y_rows = np.where(dark[int(H * 0.40) : int(H * 0.92), c_y_col - 3 : c_y_col + 4].any(axis=1))[0] + int(H * 0.40)
    y_span = max(_cl1d(y_rows, gap=12), key=lambda t: t[1] - t[0])
    yticks = _outward_yticks(dark, c_y, y_span[0] - 10, y_span[1] + 15, 7)
    y_ax, y_res = _fit_linear_axis(yticks, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0])
    r_120, r_0 = int(round(yticks[0])), int(round(yticks[-1]))

    c_lo = c_y_col + 15
    c_hi = W - 15

    masks_bot = {"0.1": _color_mask(img, "teal"), "0.2": _color_mask(img, "maroon"), "0.3": _color_mask(img, "black")}
    masks_top = {"0.1": _color_mask(img, "blue"), "0.2": _color_mask(img, "magenta"), "0.3": _color_mask(img, "olive")}
    union = np.zeros((H, W), dtype=bool)
    for m in list(masks_bot.values()) + list(masks_top.values()):
        union |= m
    col_has = np.where(union[r_120 : r_0 - 3, c_lo:c_hi].sum(axis=0) > 10)[0] + c_lo
    grp_cls = _cl1d(col_has, gap=38)
    for g in (35, 30, 45, 25):
        if len(grp_cls) == 3:
            break
        grp_cls = _cl1d(col_has, gap=g)
    if len(grp_cls) != 3:
        raise RuntimeError(f"Expected 3 groups in stacked hollow panel, found {len(grp_cls)}")

    g2_c0, g2_c1, _ = grp_cls[2]
    walls_ref = {}
    for rk in ("0.1", "0.2", "0.3"):
        mm = (masks_bot[rk] | masks_top[rk])[r_120 : r_0 - 3, g2_c0 - 5 : g2_c1 + 6]
        tall = np.where(mm.sum(axis=0) > 160)[0]
        w_cls = _cl1d(tall, gap=6)
        cL, cR = g2_c0 - 5 + w_cls[0][2], g2_c0 - 5 + w_cls[-1][2]
        walls_ref[rk] = (cL, cR, cR - cL)
    bar_w = float(np.median([v[2] for v in walls_ref.values()]))
    step_slot = 0.5 * (0.5 * (walls_ref["0.3"][0] + walls_ref["0.3"][1]) - 0.5 * (walls_ref["0.1"][0] + walls_ref["0.1"][1]))

    grp_mid_cols = []
    for g0, g1, gm in grp_cls:
        mm02 = (masks_bot["0.2"] | masks_top["0.2"])[r_120 : r_0 - 3, g0 - 5 : g1 + 6]
        cc = np.where(mm02.sum(axis=0) > 15)[0]
        grp_mid_cols.append(g0 - 5 + 0.5 * (cc[0] + cc[-1]) if len(cc) > 0 else gm)
    g_step = 0.5 * (grp_mid_cols[2] - grp_mid_cols[0])
    grp_mid_cols = [grp_mid_cols[2] - 2 * g_step, grp_mid_cols[2] - g_step, grp_mid_cols[2]]

    slot_offset = {"0.1": -1, "0.2": 0, "0.3": +1}[ratio]

    def _find_wall_top_in_mask(mask: np.ndarray, cL: int, cR: int, r_start: int, r_limit: int, *, allow_detached: bool) -> float:
        """Find top edge of a hollow bar with walls at `cL` and `cR` inside [r_limit..r_start]."""
        if r_start - r_limit < 8:
            return float(r_start)
        lw = mask[r_limit : r_start + 1, max(0, cL - 4) : min(W, cL + 5)].any(axis=1)
        rw = mask[r_limit : r_start + 1, max(0, cR - 4) : min(W, cR + 5)].any(axis=1)
        bw = lw & rw
        if not allow_detached:
            n = len(bw)
            if not any(bw[max(0, n - 5) : n]):
                return float(r_start)
            i = n - 1
            while i > 0 and (bw[i] or (i >= 1 and bw[i - 1]) or (i >= 2 and bw[i - 2])):
                i -= 1
            r_wt = r_limit + i + 1
            ra = max(r_limit, r_wt - 3)
            rb = min(r_start, r_wt + 4)
            hf = mask[ra:rb, cL + 3 : max(cL + 4, cR - 2)].mean(axis=1)
            hh = np.where(hf > 0.72)[0]
            r_top = float(ra + _cl1d(hh, gap=2)[0][2]) if len(hh) > 0 else float(r_wt)
            return float(r_start) if (r_start - r_top < 4.0) else r_top
        else:
            bw_idx = np.where(bw)[0]
            if len(bw_idx) == 0:
                return float(r_start)
            cls = _cl1d(bw_idx, gap=5)
            valid_cls = [
                c
                for c in cls
                if (r_limit + c[1]) >= r_start - 45
                and ((c[1] - c[0] + 1) >= 12 or ((r_limit + c[1]) >= r_start - 6 and (c[1] - c[0] + 1) >= 6))
            ]
            if not valid_cls:
                return float(r_start)
            top_cl = valid_cls[-1]
            r_wt = r_limit + top_cl[0]
            ra = max(r_limit, r_wt - 2)
            rb = min(r_start, r_wt + 8)
            hf = mask[ra:rb, cL + 2 : max(cL + 3, cR - 1)].mean(axis=1)
            hh = np.where(hf > 0.78)[0]
            if len(hh) > 0:
                return float(ra + _cl1d(hh, gap=2)[0][2])
            return float(r_wt)

    out: list[dict] = []
    for gi, gmid in enumerate(grp_mid_cols):
        xc = gmid + slot_offset * step_slot
        cL = int(round(xc - 0.5 * bar_w))
        cR = int(round(xc + 0.5 * bar_w))
        r_net = _find_wall_top_in_mask(masks_bot[ratio], cL, cR, r_0 - 3, r_120, allow_detached=False)
        if r_net >= r_0 - 6:
            r_net = float(r_0)
        if layer == "network":
            r_use = r_net
            note = f"ratio={ratio} Network bottom bar"
        elif layer == "total":
            r_102 = int(round(y_ax.to_pixel(103.0)))
            r_tot = _find_wall_top_in_mask(masks_top[ratio], cL, cR, int(round(r_net)) - 1, r_102, allow_detached=True)
            ci_m = int(round(xc))
            below_ax = masks_top[ratio][min(H - 1, r_0 + 5) : min(H, r_0 + 15), max(0, ci_m - 6) : min(W, ci_m + 7)].mean()
            if r_net >= r_0 - 1 and (r_0 - r_tot < 13.0) and below_ax > 0.45:
                r_use = float(r_0)
            else:
                r_use = r_tot if (r_net - r_tot >= 4.0) else r_net
            note = f"ratio={ratio} Network+Ring/Network stacked top"
        else:
            raise ValueError(layer)
        half_px = math.hypot(2.5, y_res)
        out.append(_make_row(float(gi + 1), r_use, xc, y_ax, half_px, note=f"{note}; y_resid={y_res:.2f}px"))
    return out


def read_fig1d(img: np.ndarray, ratio: str, layer: str) -> list[dict]:
    return _read_stacked_hollow_panel(img, ratio=ratio, layer=layer)


def read_fig2e(img: np.ndarray, ratio: str, layer: str) -> list[dict]:
    return _read_stacked_hollow_panel(img, ratio=ratio, layer=layer)


# ---------------------------------------------------------------------------
# 4. Grouped white/grey bar charts: Fig. 3c, Fig. 3d
# ---------------------------------------------------------------------------

def _read_fig3cd(img: np.ndarray, series: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 140) & (sat < 40)
    # Pure black border ink (< 75) cleanly excludes medium-grey Periphery bar fill (gray ~ 125)!
    black_border = (gray < 75) & (sat < 35)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.55, x_r_hi_frac=0.80)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, 6)
    y_ax, y_res = _fit_linear_axis(yticks, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])
    xticks = _outward_xticks(dark, r_x, int(round(c_y)) + 20, x1, 2)
    ri_0 = int(round(r_x))

    H, W = dark.shape
    leg_rows = np.where(dark[int(H * 0.15) : int(round(yticks[2])), int(round(c_y)) + 20 : W - 20].mean(axis=1) > 0.55)[0] + int(H * 0.15)
    r_lim = int(leg_rows[-1]) + 18 if len(leg_rows) > 0 else int(round(yticks[2]))

    out: list[dict] = []
    for gi, xt in enumerate(xticks):
        ci = int(round(xt))
        c0, c1 = (ci - 85, ci - 20) if series == "center" else (ci + 20, ci + 85)
        hf = black_border[r_lim : ri_0 - 4, c0:c1].mean(axis=1)
        hits = np.where(hf > 0.80)[0]
        r_bar = float(r_lim + _cl1d(hits, gap=3)[0][2])
        xc = 0.5 * (c0 + c1)
        half_px = math.hypot(2.0, y_res)
        out.append(_make_row(float(gi + 1), r_bar, xc, y_ax, half_px, note=f"{series}; y_resid={y_res:.2f}px"))
    return out


def read_fig3c(img: np.ndarray, series: str) -> list[dict]:
    return _read_fig3cd(img, series)


def read_fig3d(img: np.ndarray, series: str) -> list[dict]:
    return _read_fig3cd(img, series)


# ---------------------------------------------------------------------------
# 5. Fig. 3b & Supp. Fig. 9b: black & red line scans with outward ticks
# ---------------------------------------------------------------------------

def read_fig3b(img: np.ndarray, at: list[float], ratio: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.60, x_r_hi_frac=0.82)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, 6)
    y_ax, y_res = _fit_linear_axis(yticks, [200.0, 160.0, 120.0, 80.0, 40.0, 0.0])
    # 6 major outward x-ticks: 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
    xticks = _outward_xticks(dark, r_x, int(round(c_y)) - 15, x1 + 10, 6)
    x_ax, x_res = _fit_linear_axis(xticks, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    mask = _color_mask(img, "black" if ratio == "0.3" else "red")
    r_min = int(round(yticks[0])) + 5
    r_max = int(round(r_x)) - 5
    mask[:, : int(round(c_y)) + 8] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, r_min, r_max, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"ratio={ratio}; "
    )


def read_sfig9b(img: np.ndarray, at: list[float], ratio: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.65, x_r_hi_frac=0.85)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, 6)
    y_ax, y_res = _fit_linear_axis(yticks, [200.0, 160.0, 120.0, 80.0, 40.0, 0.0])
    xticks = _outward_xticks(dark, r_x, int(round(c_y)) - 15, x1 + 15, 9, d0=10, d1=22)
    x_ax, x_res = _fit_linear_axis(xticks, [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0])

    mask = _color_mask(img, "black" if ratio == "0.3" else "red")
    r_min = int(round(yticks[0])) + 6
    r_max = int(round(r_x)) - 5
    mask[:, : int(round(c_y)) + 10] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, r_min, r_max, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"fascin={ratio}; "
    )


# ---------------------------------------------------------------------------
# 6. Fig. 4c & Fig. 4d: MATLAB-style inward-tick profiles (4d has dual y-axes!)
# ---------------------------------------------------------------------------

def read_fig4c(img: np.ndarray, at: list[float], species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.70, x_r_hi_frac=0.88)
    yticks = _inward_yticks(dark, c_y, y0 - 10, int(round(r_x)) - 15, 4, d0=6, d1=17)
    y_ax, y_res = _fit_linear_axis(yticks, [2.0, 1.5, 1.0, 0.5])
    ri = int(round(r_x))
    mid_est = int(round(0.5 * (c_y + x1)))
    xt_prof = dark[ri - 17 : ri - 6, mid_est - 60 : mid_est + 60].mean(axis=0)
    xt_mid_cls = _cl1d(np.where(xt_prof > 0.45)[0], gap=4)
    if xt_mid_cls:
        c_25 = mid_est - 60 + xt_mid_cls[0][2]
        x_ax, x_res = _fit_linear_axis([c_y, c_25, float(x1)], [0.0, 25.0, 50.0])
    else:
        x_ax, x_res = _fit_linear_axis([c_y, float(x1)], [0.0, 50.0])

    mask = _color_mask(img, "magenta" if species == "fascin" else "blue")
    r_min = int(round(yticks[0])) - 10
    r_max = int(round(r_x)) - 4
    mask[: int(round(yticks[1])) - 40, int(round(x_ax.to_pixel(32.0))) :] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, r_min, r_max, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"{species}; "
    )


def read_fig4d(img: np.ndarray, at: list[float], species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.70, x_r_hi_frac=0.88)
    xticks = _inward_xticks(dark, r_x, int(round(c_y)) - 10, x1 + 10, 3, d0=6, d1=17)
    x_ax, x_res = _fit_linear_axis(xticks, [0.0, 200.0, 400.0])
    yticks = _inward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 10, 6, d0=6, d1=17)
    if species == "fascin":
        y_vals = [1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0]
        mask = _color_mask(img, "magenta")
    else:
        y_vals = [1000.0, 900.0, 800.0, 700.0, 600.0, 500.0]
        mask = _color_mask(img, "blue")
    y_ax, y_res = _fit_linear_axis(yticks, y_vals)

    r_min = int(round(yticks[0])) + 5
    r_max = int(round(r_x)) - 5
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, r_min, r_max, win_cols=4, calib_resid=max(x_res, y_res), note_prefix=f"{species}; "
    )


# ---------------------------------------------------------------------------
# 7. Fig. 4e: Bundle persistence length L_p bars + symmetric SEM error bars
# ---------------------------------------------------------------------------

def read_fig4e(img: np.ndarray, condition: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.72, x_r_hi_frac=0.88)
    ci = int(round(c_y))
    sub_g = gray[y0 - 15 : int(round(r_x)) + 15, ci - 24 : ci - 8]
    sub_d = dark[y0 - 15 : int(round(r_x)) + 15, ci - 24 : ci - 8]
    yt_prof = np.where(sub_d, 255.0 - sub_g.astype(np.float32), 0.0).mean(axis=1)
    yticks = []
    for a, b, _ in _cl1d(np.where(yt_prof > 60.0)[0], gap=6):
        w = yt_prof[a : b + 1]
        yticks.append(float(y0 - 15 + np.sum(np.arange(a, b + 1) * w) / np.sum(w)))
    yticks = _select_regular(yticks, 6)
    y_ax, y_res = _fit_linear_axis(yticks, [250.0, 200.0, 150.0, 100.0, 50.0, 0.0])

    xticks = _outward_xticks(dark, r_x, ci + 20, x1 + 10, 5)
    blue = _color_mask(img, "blue")
    mag = _color_mask(img, "magenta")
    r_lim = int(round(yticks[0]))
    r_base = int(round(r_x)) - 5

    items = [(0, 0.1, blue), (1, 0.3, blue)] if condition == "no_fascin" else [(2, 0.1, mag), (3, 0.2, mag), (4, 0.3, mag)]

    out: list[dict] = []
    for bi, ratio_val, mask in items:
        xc = xticks[bi]
        ci_b = int(round(xc))
        bar_frac = mask[r_lim:r_base, ci_b - 75 : ci_b + 75].mean(axis=1)
        top_cl = _cl1d(np.where(bar_frac > 0.88)[0], gap=4)[0]
        r_bar = float(r_lim + top_cl[2])

        cap_frac = mask[r_lim:r_base, ci_b - 40 : ci_b + 40].mean(axis=1)
        cap_cls = [r_lim + m for _, _, m in _cl1d(np.where(cap_frac > 0.85)[0], gap=4) if abs((r_lim + m) - r_bar) > 10]
        above = [r for r in cap_cls if r < r_bar]
        below = [r for r in cap_cls if r > r_bar]
        eb_hi = above[-1] if above else None
        if eb_hi is not None and below:
            d_up = r_bar - eb_hi
            eb_lo = min(below, key=lambda r: abs((r - r_bar) - d_up))
        elif below:
            eb_lo = below[0]
            eb_hi = r_bar - (eb_lo - r_bar)
        else:
            eb_lo = None

        half_px = math.hypot(2.0, y_res)
        out.append(
            _make_row(
                ratio_val,
                r_bar,
                xc,
                y_ax,
                half_px,
                note=f"Fig. 4e {condition} ratio={ratio_val}; y_resid={y_res:.2f}px",
                bar_lo_px=eb_lo,
                bar_hi_px=eb_hi,
                bar_half_px=2.5,
            )
        )
    return out


# ---------------------------------------------------------------------------
# 8. Fig. 4h (top: Single component, bottom: Combined) Binder cumulant vs time
# ---------------------------------------------------------------------------

def read_fig4h(img: np.ndarray, species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 140) & (sat < 40)
    H, W = dark.shape

    hcls = _cl1d(np.where(dark.mean(axis=1) > 0.50)[0], gap=10)
    rT, rB = hcls[0][2], hcls[-1][2]
    # Vertical box spines have > 85% coverage over [rT..rB] (cleanly rejects rotated y-label text!)
    vcls = _cl1d(np.where(dark[int(rT) : int(rB), :].mean(axis=0) > 0.85)[0], gap=10)
    cL, cR = vcls[0][2], vcls[-1][2]

    # 6 interior inward y-ticks (0.995..0.970) + top spine (1.000)
    yt_in = _inward_yticks(dark, cL, int(rT) + 40, int(rB) - 18, 6, d0=10, d1=22)
    yticks = [float(rT)] + yt_in
    y_ax, y_res = _fit_linear_axis(yticks, [1.000, 0.995, 0.990, 0.985, 0.980, 0.975, 0.970])

    # 4 clean interior bottom x-ticks: 20, 40, 60, 80 s
    xt_in = _inward_xticks(dark, rB, int(cL) + 40, int(cR) - 40, 4, d0=10, d1=22)
    x_ax, x_res = _fit_linear_axis(xt_in, [20.0, 40.0, 60.0, 80.0])

    red = _color_mask(img, "red")
    blk = _color_mask(img, "black")

    interior = np.zeros((H, W), dtype=bool)
    interior[int(rT) + 32 : int(rB) - 32, int(cL) + 32 : int(cR) - 32] = True
    legend_mask = np.zeros((H, W), dtype=bool)
    legend_mask[int(rT + 0.55 * (rB - rT)) :, int(cL + 0.52 * (cR - cL)) :] = True
    valid = interior & (~legend_mask)

    red_p = red & valid
    blk_p = blk & valid

    red_cols = np.where(red_p.any(axis=0))[0]
    blk_cols = np.where(blk_p.any(axis=0))[0]
    c_first = 0.5 * (float(red_cols[0]) + float(blk_cols[0])) + 5.0
    c_last = 0.5 * (float(red_cols[-1]) + float(blk_cols[-1])) - 5.0

    est_cols = np.linspace(c_first, c_last, 10)
    station_cols = []
    for ce in est_cols:
        ci = int(round(ce))
        prof_r = red_p[:, max(0, ci - 14) : min(W, ci + 15)].sum(axis=0)
        if prof_r.max() > 35:
            top_c = np.where(prof_r >= 0.85 * prof_r.max())[0]
            station_cols.append(float(max(0, ci - 14) + np.mean(top_c)))
        else:
            station_cols.append(float(ce))
    k_idx = np.arange(10, dtype=float)
    st_slope, st_int = np.polyfit(k_idx, station_cols, 1)
    station_cols = [float(st_int + st_slope * k) for k in k_idx]

    mask = red_p if species == "alpha_actinin" else blk_p
    out: list[dict] = []
    for k, cs in enumerate(station_cols):
        ci = int(round(cs))
        t_val = float(x_ax.to_value(cs))
        col_rows = np.where(mask[:, max(0, ci - 2) : min(W, ci + 3)].any(axis=1))[0]
        if len(col_rows) == 0:
            out.append({"x": t_val, "y": None, "found": False, "note": "no pixels at station"})
            continue
        cls = _cl1d(col_rows, gap=12)
        main_cl = max(cls, key=lambda t: t[1] - t[0])
        r_eb_top, r_eb_bot = float(main_cl[0]), float(main_cl[1])
        r_center = 0.5 * (r_eb_top + r_eb_bot)

        if (r_eb_bot - r_eb_top) < 18.0:
            eb_half = max(0.0, 0.5 * (r_eb_bot - r_eb_top) - 4.0)
            eb_hi_px = r_center - max(eb_half, 2.0)
            eb_lo_px = r_center + max(eb_half, 2.0)
            note = f"Fig. 4h {species}; error bar near line thickness; resid=({x_res:.2f},{y_res:.2f})px"
        else:
            eb_hi_px = r_eb_top + 1.5
            eb_lo_px = r_eb_bot - 1.5
            note = f"Fig. 4h {species}; symmetric SD error bar; resid=({x_res:.2f},{y_res:.2f})px"

        half_px = math.hypot(3.0, y_res)
        out.append(
            _make_row(
                round(t_val, 2),
                r_center,
                cs,
                y_ax,
                half_px,
                note=note,
                bar_lo_px=eb_lo_px,
                bar_hi_px=eb_hi_px,
                bar_half_px=3.0,
            )
        )
    return out


# ---------------------------------------------------------------------------
# 9. Supplementary Fig. 3c: Ring radial deviation d/d_av vs s/p
# ---------------------------------------------------------------------------

def read_sfig3c(img: np.ndarray, at: list[float], size: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.68, x_r_hi_frac=0.86)
    xticks = _inward_xticks(dark, r_x, int(round(c_y)) - 15, x1 + 15, 3, d0=5, d1=14)
    x_ax, x_res = _fit_linear_axis(xticks, [0.0, 0.5, 1.0])
    ci = int(round(c_y))
    mid_r_est = int(round(0.5 * (y0 + r_x)))
    yt_prof = dark[mid_r_est - 60 : mid_r_est + 60, ci + 5 : ci + 15].mean(axis=1)
    yt_mid = mid_r_est - 60 + _cl1d(np.where(yt_prof > 0.45)[0], gap=4)[0][2]
    y_ax, y_res = _fit_linear_axis([float(y0) + 3.0, float(yt_mid), float(r_x)], [1.5, 1.0, 0.5])

    mask = _color_mask(img, "sfig3_blue" if size == "small" else "sfig3_red")
    r_leg = int(round(y_ax.to_pixel(1.22)))
    mask[:r_leg, :] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, r_leg, int(round(r_x)) - 10, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"{size}; "
    )


# ---------------------------------------------------------------------------
# 10. Supplementary Fig. 7b & 7d
# ---------------------------------------------------------------------------

def read_sfig7b(img: np.ndarray, at: list[float]) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.70, x_r_hi_frac=0.88)
    yticks = _inward_yticks(dark, c_y, y0 - 15, int(round(r_x)) - 5, 3, d0=6, d1=17)
    y_ax, y_res = _fit_linear_axis(yticks, [200.0, 150.0, 100.0])
    ri = int(round(r_x))
    mid_est = int(round(0.5 * (c_y + x1)))
    xt_prof = dark[ri - 17 : ri - 6, mid_est - 60 : mid_est + 60].mean(axis=0)
    xt_mid_cls = _cl1d(np.where(xt_prof > 0.45)[0], gap=4)
    if xt_mid_cls:
        c_5 = mid_est - 60 + xt_mid_cls[0][2]
        x_ax, x_res = _fit_linear_axis([c_y, c_5, float(x1)], [0.0, 5.0, 10.0])
    else:
        x_ax, x_res = _fit_linear_axis([c_y, float(x1)], [0.0, 10.0])

    mask = _color_mask(img, "black")
    mask[:, : int(round(c_y)) + 8] = False
    mask[int(round(r_x)) - 6 :, :] = False
    c_d0, c_d1 = int(round(x_ax.to_pixel(4.15))), int(round(x_ax.to_pixel(5.05)))
    r_d0, r_d1 = int(round(y_ax.to_pixel(125.0))), int(round(y_ax.to_pixel(103.0)))
    mask[r_d0:r_d1, c_d0:c_d1] = False

    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, int(round(yticks[0])) - 10, int(round(r_x)) - 6, win_cols=3, calib_resid=max(x_res, y_res), note_prefix="SFig7b; "
    )


def read_sfig7d(img: np.ndarray, series: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 140) & (sat < 40)
    blk_ink = img[:, :, 0] < 85
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.75, x_r_hi_frac=0.88)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, 7)
    y_ax, y_res = _fit_linear_axis(yticks, [3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0])
    xticks = _outward_xticks(dark, r_x, int(round(c_y)) + 40, x1 + 10, 4)

    ri_0 = int(round(r_x))
    ri_top = int(round(yticks[0]))

    def _read_bar_at(xc: float, ratio_val: float, label: str) -> dict:
        ci = int(round(xc))
        hf_bar = blk_ink[ri_top : ri_0 - 10, ci - 48 : ci + 49].mean(axis=1)
        r_bar = float(ri_top + _cl1d(np.where(hf_bar > 0.85)[0], gap=4)[0][2])
        hf_cap = blk_ink[ri_top : ri_0 - 10, ci - 22 : ci + 23].mean(axis=1)
        cap_cls = [ri_top + m for _, _, m in _cl1d(np.where(hf_cap > 0.85)[0], gap=4) if abs((ri_top + m) - r_bar) > 15]
        above = [r for r in cap_cls if r < r_bar]
        below = [r for r in cap_cls if r > r_bar]
        eb_hi = above[-1] if above else None
        eb_lo = below[0] if below else None
        if eb_hi is not None and eb_lo is None:
            eb_lo = r_bar + (r_bar - eb_hi)
        half_px = math.hypot(2.5, y_res)
        return _make_row(
            ratio_val,
            r_bar,
            xc,
            y_ax,
            half_px,
            note=f"SFig7d {label}; y_resid={y_res:.2f}px",
            bar_lo_px=eb_lo,
            bar_hi_px=eb_hi,
            bar_half_px=3.0,
        )

    if series == "bulk":
        return [
            _read_bar_at(xticks[2], 0.05, "bulk 0.05"),
            _read_bar_at(xticks[1], 0.1, "bulk 0.1"),
            _read_bar_at(xticks[0], 0.5, "bulk 0.5"),
        ]
    elif series == "aster":
        return [_read_bar_at(xticks[3], 0.1, "central aster 0.1")]
    raise ValueError(series)


# ---------------------------------------------------------------------------
# 11. Supplementary Fig. 8b, 8c, 8e
# ---------------------------------------------------------------------------

def read_sfig8b(img: np.ndarray, at: list[float], species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.72, x_r_hi_frac=0.88)
    xticks = _inward_xticks(dark, r_x, int(round(c_y)) - 15, x1, 4, d0=5, d1=16)
    x_ax, x_res = _fit_linear_axis(xticks, [0.0, 100.0, 200.0, 300.0])
    yticks = _inward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 10, 4, d0=6, d1=17)
    y_ax, y_res = _fit_linear_axis(yticks, [2.0, 1.5, 1.0, 0.5])

    mask = _color_mask(img, "sfig8_green" if species == "actin" else "sfig8_magenta")
    mask[: int(round(yticks[1])) - 40, int(round(x_ax.to_pixel(275.0))) :] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, int(round(yticks[0])) - 15, int(round(r_x)) - 4, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"SFig8b {species}; "
    )


def read_sfig8c(img: np.ndarray, at: list[float], species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.72, x_r_hi_frac=0.88)
    xticks = _inward_xticks(dark, r_x, int(round(c_y)) - 15, x1 + 15, 4, d0=5, d1=14)
    x_ax, x_res = _fit_linear_axis(xticks, [0.0, 5.0, 10.0, 15.0])

    ci = int(round(c_y))
    prof = dark[y0 : int(round(r_x)) + 8, ci + 5 : ci + 15].mean(axis=1)
    cands = [y0 + m for _, _, m in _cl1d(np.where(prof > 0.45)[0], gap=5)]
    if not any(abs(c - r_x) < 6 for c in cands):
        cands.append(float(r_x))
    yticks = _select_regular(cands, 5, tol=4.5)
    y_ax, y_res = _fit_linear_axis(yticks, [2.0, 1.5, 1.0, 0.5, 0.0])

    mask = _color_mask(img, "sfig8_green" if species == "actin" else "sfig8_magenta")
    mask[: int(round(y_ax.to_pixel(1.3))), int(round(x_ax.to_pixel(13.0))) :] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, max(10, int(round(y_ax.to_pixel(2.45)))), int(round(r_x)) - 4, win_cols=3, calib_resid=max(x_res, y_res), note_prefix=f"SFig8c {species}; "
    )


def read_sfig8e(img: np.ndarray, at: list[float], species: str) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 150) & (sat < 45)
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.74, x_r_hi_frac=0.88)
    yticks = _inward_yticks(dark, c_y, y0 + 20, int(round(r_x)) - 15, 4, d0=6, d1=17)
    y_ax, y_res = _fit_linear_axis(yticks, [1.4, 1.2, 1.0, 0.8])
    ri = int(round(r_x))
    mid_est = int(round(0.5 * (c_y + x1)))
    xt_prof = dark[ri - 17 : ri - 6, mid_est - 60 : mid_est + 60].mean(axis=0)
    xt_mid_cls = _cl1d(np.where(xt_prof > 0.45)[0], gap=4)
    if xt_mid_cls:
        c_2 = mid_est - 60 + xt_mid_cls[0][2]
        x_ax, x_res = _fit_linear_axis([c_y, c_2, float(x1)], [0.0, 2.0, 4.0])
    else:
        x_ax, x_res = _fit_linear_axis([c_y, float(x1)], [0.0, 4.0])

    cmap = {"actin": "sfig8e_green", "alpha_actinin": "sfig8e_magenta", "fascin": "sfig8e_blue"}
    mask = _color_mask(img, cmap[species])
    mask[: int(round(y_ax.to_pixel(1.26))), int(round(x_ax.to_pixel(2.75))) :] = False
    mask[:, : int(round(c_y)) + 8] = False
    return _read_curve_on_mask(
        mask, x_ax, y_ax, at, max(10, int(round(y_ax.to_pixel(1.56)))), int(round(r_x)) - 4, win_cols=4, calib_resid=max(x_res, y_res), note_prefix=f"SFig8e {species}; "
    )


# ---------------------------------------------------------------------------
# 12. Supplementary Fig. 9c: Bundle persistence length L_p bars + SEM error bars
# ---------------------------------------------------------------------------

def read_sfig9c(img: np.ndarray) -> list[dict]:
    gray, sat = _lum_sat(img)
    dark = (gray < 140) & (sat < 40)
    blk_ink = gray < 95
    r_x, c_y, x0, x1, y0, y1 = _find_l_axes(gray, sat, x_r_lo_frac=0.65, x_r_hi_frac=0.82)
    yticks = _outward_yticks(dark, c_y, y0 - 15, int(round(r_x)) + 15, 7)
    y_ax, y_res = _fit_linear_axis(yticks, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0])

    # Detect exact bar center columns directly from the 4 vertical bar side walls (avoids the shifted x-tick under 0.3!)
    ri_0 = int(round(r_x))
    ri_top = int(round(yticks[0])) - 20
    wall_cols = np.where(blk_ink[int(round(yticks[3])) : ri_0 - 10, int(round(c_y)) + 40 : x1].mean(axis=0) > 0.80)[0] + int(round(c_y)) + 40
    wall_cls = _cl1d(wall_cols, gap=10)
    if len(wall_cls) != 4:
        raise RuntimeError(f"Expected 4 vertical bar walls in SFig9c, found {wall_cls}")
    bar_centers = [0.5 * (wall_cls[0][2] + wall_cls[1][2]), 0.5 * (wall_cls[2][2] + wall_cls[3][2])]

    out: list[dict] = []
    for ratio_val, xc in zip([0.1, 0.3], bar_centers):
        ci = int(round(xc))
        hf_bar = blk_ink[ri_top : ri_0 - 10, ci - 115 : ci + 116].mean(axis=1)
        r_bar = float(ri_top + _cl1d(np.where(hf_bar > 0.85)[0], gap=4)[0][2])
        hf_cap = blk_ink[ri_top : ri_0 - 10, ci - 20 : ci + 21].mean(axis=1)
        cap_cls = [ri_top + m for _, _, m in _cl1d(np.where(hf_cap > 0.85)[0], gap=4) if abs((ri_top + m) - r_bar) > 20]
        above = [r for r in cap_cls if r < r_bar]
        below = [r for r in cap_cls if r > r_bar]
        eb_hi = above[-1] if above else None
        eb_lo = below[0] if below else None
        half_px = math.hypot(2.5, y_res)
        out.append(
            _make_row(
                ratio_val,
                r_bar,
                xc,
                y_ax,
                half_px,
                note=f"SFig9c fascin={ratio_val}; y_resid={y_res:.2f}px",
                bar_lo_px=eb_lo,
                bar_hi_px=eb_hi,
                bar_half_px=3.0,
            )
        )
    return out
