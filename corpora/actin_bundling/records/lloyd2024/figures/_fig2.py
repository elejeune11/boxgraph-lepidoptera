"""Shared readers for Fig. 2D/F linear regression + 95% CI curves and Fig. 2J-L notched boxplots & adult reference lines (PDF p. 4, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D
from _common import find_runs, run_centers, y_interval


def read_lm_and_ci(image: np.ndarray, xa: D.Axis, ya: D.Axis,
                   c0: float, r1: float, r_top: int, c_right: int,
                   series: str, at: list[float],
                   x_fit_range: tuple[float, float] | None = None) -> list[dict]:
    """Robustly extract a ggplot2 geom_smooth(method='lm') straight regression line and its 95% CI ribbon."""
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0
    H, W = gray.shape

    if series == "blue":
        line_m = (B - R > 90.0) & (R < 115.0) & (G < 185.0)
        rib_m = (B - R > 12.0) & (B > 170.0)
    else:
        line_m = (gray < 105.0) & (np.abs(B - R) < 15.0)
        rib_m = (gray < 230.0) & (np.abs(B - R) < 15.0)

    c_min_box = int(round(c0)) + 5
    r_max_box = int(round(r1)) - 5
    line_m[:r_top, :] = False
    line_m[r_max_box:, :] = False
    line_m[:, :c_min_box] = False
    line_m[:, c_right:] = False
    rib_m[:r_top, :] = False
    rib_m[r_max_box:, :] = False
    rib_m[:, :c_min_box] = False
    rib_m[:, c_right:] = False

    # Collect clean thin stroke runs across the line's horizontal extent
    x_lo, x_hi = x_fit_range if x_fit_range is not None else (min(at), max(at))
    c_start = max(c_min_box, int(round(xa.to_pixel(x_lo))))
    c_end = min(c_right - 1, int(round(xa.to_pixel(x_hi))))
    cols_ok, rows_ok = [], []
    for c in range(c_start, c_end + 1, 2):
        rs = np.where(line_m[:, c])[0]
        if len(rs) == 0:
            continue
        for s, e, rc in find_runs(rs, gap=3):
            if 2 <= (e - s + 1) <= 14:
                cols_ok.append(float(c))
                rows_ok.append(rc)
    cols_ok = np.array(cols_ok, dtype=float)
    rows_ok = np.array(rows_ok, dtype=float)
    if len(cols_ok) < 20:
        raise RuntimeError(f"Too few regression line samples for {series}: {len(cols_ok)}")

    # Deterministic RANSAC to isolate the exact straight regression line ignoring marker overlaps
    rng = np.random.default_rng(0)
    best_cnt, best_mb = -1, (0.0, 0.0)
    n = len(cols_ok)
    for _ in range(500):
        i, j = rng.integers(0, n, size=2)
        if abs(cols_ok[i] - cols_ok[j]) < 50.0:
            continue
        m_try = (rows_ok[j] - rows_ok[i]) / (cols_ok[j] - cols_ok[i])
        b_try = rows_ok[i] - m_try * cols_ok[i]
        cnt = int((np.abs(rows_ok - (m_try * cols_ok + b_try)) < 3.5).sum())
        if cnt > best_cnt:
            best_cnt = cnt
            best_mb = (m_try, b_try)
    m, b = best_mb
    inl = np.abs(rows_ok - (m * cols_ok + b)) < 4.0
    m, b = np.polyfit(cols_ok[inl], rows_ok[inl], 1)
    line_res = float(np.median(np.abs(rows_ok[inl] - (m * cols_ok[inl] + b))))

    # Measure symmetric 95% CI ribbon half-width across dense columns:
    # At columns where the other series' ribbon overlaps one side, or a marker sticks out on one side,
    # nonwhite extent on the unoccluded side gives the clean ribbon boundary.
    nonwhite = (R < 246.0) | (G < 246.0) | (B < 246.0)
    nonwhite[:r_top, :] = False
    nonwhite[r_max_box:, :] = False
    nonwhite[:, :c_min_box] = False
    nonwhite[:, c_right:] = False

    # Opposite series mask so we know when one side of the CI ribbon is occluded by the other ribbon
    if series == "blue":
        other_m = (gray < 232.0) & (np.abs(B - R) <= 15.0)
    else:
        other_m = (B - R > 15.0) & (B > 150.0)

    c_grid = np.linspace(c_start + 8.0, c_end - 8.0, 140)
    hw_up = np.zeros(len(c_grid), dtype=float)
    hw_dn = np.zeros(len(c_grid), dtype=float)
    occ_up = np.zeros(len(c_grid), dtype=bool)
    occ_dn = np.zeros(len(c_grid), dtype=bool)
    for idx_c, cf in enumerate(c_grid):
        c = int(round(cf))
        r_l = int(round(m * cf + b))
        up = 0
        while r_l - up > r_top and (rib_m[r_l - up, c] or (r_l - up - 2 > r_top and rib_m[r_l - up - 2, c])):
            up += 1
        if r_l - up - 2 > r_top and other_m[max(r_top, r_l - up - 3):r_l - up + 1, c].any():
            occ_up[idx_c] = True
        dn = 0
        while r_l + dn < r_max_box and (rib_m[r_l + dn, c] or (r_l + dn + 2 < r_max_box and rib_m[r_l + dn + 2, c])):
            dn += 1
        if r_l + dn + 2 < r_max_box and other_m[r_l + dn:min(r_max_box, r_l + dn + 4), c].any():
            occ_dn[idx_c] = True
        hw_up[idx_c] = float(up)
        hw_dn[idx_c] = float(dn)

    # Use unoccluded side whenever the opposite side hits the other series' ribbon; otherwise take min(up, dn)
    # so stray scatter markers extending outside one side do not inflate the ribbon width.
    hws = np.minimum(hw_up, hw_dn)
    hws = np.where(occ_up & ~occ_dn, hw_dn, hws)
    hws = np.where(occ_dn & ~occ_up, hw_up, hws)
    clipped = (hws < 18.0) & (np.maximum(hw_up, hw_dn) > 25.0)
    hws = np.where(clipped, np.maximum(hw_up, hw_dn), hws)

    # Fit positive quadratic to hw^2 (the exact mathematical shape of an OLS 95% confidence band: hw^2 = a0 + a1*u + a2*u^2!)
    c_mid = 0.5 * (c_start + c_end)
    c_scale = max(0.5 * (c_end - c_start), 1.0)
    u = (c_grid - c_mid) / c_scale
    both_occ = occ_up & occ_dn
    valid = (~both_occ) & (hws >= 15.0) & (hws <= 320.0)
    for _ in range(5):
        if valid.sum() < 12:
            break
        p2 = np.polyfit(u[valid], hws[valid] ** 2, 2)
        pred_hw = np.sqrt(np.maximum(np.polyval(p2, u), 100.0))
        err = np.abs(hws - pred_hw)
        thr = max(5.0, float(np.percentile(err[valid], 70)) * 1.4)
        valid = (~both_occ) & (hws >= 15.0) & (hws <= 320.0) & (err < thr)
    p2 = np.polyfit(u[valid], hws[valid] ** 2, 2)
    # Guarantee positive definite (convex upward) CI band
    if p2[0] < 0.0:
        p1 = np.polyfit(u[valid], hws[valid] ** 2, 1)
        p2 = np.array([0.0, p1[0], p1[1]])

    out = []
    for xv in at:
        cf = float(xa.to_pixel(float(xv)))
        rf = float(m * cf + b)
        uf = (cf - c_mid) / c_scale
        hw_px = float(np.sqrt(max(np.polyval(p2, uf), 100.0)))
        r_top_ci = rf - hw_px
        r_bot_ci = rf + hw_px

        y_val, lo, hi, _, _ = y_interval(ya, rf, max(2.5, line_res + 1.5))
        bar_hi = float(ya.to_value(r_top_ci))
        bar_lo = float(ya.to_value(r_bot_ci))
        h_ci = 5.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_top_ci - h_ci))
        bar_lo_outer = float(ya.to_value(r_bot_ci + h_ci))
        out.append({
            "x": float(xv),
            "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(bar_lo, bar_hi)), "bar_hi": float(max(bar_lo, bar_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [cf, rf],
            "px_lo": float(r_top_ci), "px_hi": float(r_bot_ci),
            "note": f"{series} lm +- 95% CI (line_res={line_res:.2f} px, rx={xa.residual_px:.2f} px, ry={ya.residual_px:.2f} px)",
        })
    return out


def calibrate_fig2_scatter(image: np.ndarray, x_labels: list[float],
                           y_labels: list[float]) -> tuple[D.Axis, D.Axis, float, float]:
    gray = image.astype(float).mean(axis=2)
    c0 = 250.0 + float(np.argmin(gray[450:1850, 250:500].mean(axis=0)))
    r1 = 1850.0 + float(np.argmin(gray[1850:2080, 450:2000].mean(axis=1)))
    yt = run_centers(np.where(gray[420:int(r1) - 10, int(round(c0)) - 7] < 180.0)[0] + 420, gap=4)
    xt = run_centers(np.where(gray[int(round(r1)) + 7, int(round(c0)) + 10:2150] < 180.0)[0] + int(round(c0)) + 10, gap=4)
    if len(xt) != len(x_labels) or len(yt) != len(y_labels):
        raise RuntimeError(f"Fig. 2 scatter ticks mismatch: xt={xt} vs {x_labels}, yt={yt} vs {y_labels}")
    xa = D.Axis("linear", list(zip(xt, x_labels)))
    ya = D.Axis("linear", list(zip(yt, y_labels)))
    return xa, ya, c0, r1


def _calibrate_fig2_jkl(image: np.ndarray, y_labels: list[float]) -> tuple[D.Axis, float, float, float, float, list[float]]:
    gray = image.astype(float).mean(axis=2)
    H, W = gray.shape
    c0 = 250.0 + float(np.argmin(gray[200:2000, 250:600].mean(axis=0)))
    r1 = 2000.0 + float(np.argmin(gray[2000:2250, 600:2300].mean(axis=1)))
    # Top frame border r0 (in rows 50..250)
    r0 = 50.0 + float(np.argmin(gray[50:250, int(c0) + 100:W - 100].mean(axis=1)))
    c1 = W - 200.0 + float(np.argmin(gray[int(r0) + 100:int(r1) - 100, W - 200:W - 5].mean(axis=0)))
    yt = run_centers(np.where(gray[int(r0) + 5:int(r1) - 10, int(round(c0)) - 15] < 150.0)[0] + int(r0) + 5, gap=4)
    xt = run_centers(np.where(gray[int(round(r1)) + 15, int(round(c0)) + 10:int(c1) - 10] < 150.0)[0] + int(round(c0)) + 10, gap=4)
    if len(yt) != len(y_labels) or len(xt) != 2:
        raise RuntimeError(f"Fig. 2J-L ticks mismatch: yt={yt}, xt={xt}")
    ya = D.Axis("linear", list(zip(yt, y_labels)))
    return ya, c0, r0, c1, r1, xt


def read_fig2_boxplot(image: np.ndarray, y_labels: list[float], what: str,
                      has_dashed_lines: bool = False) -> list[dict]:
    """Read notched boxplot medians and full whisker/dot ranges for black (x=1) and blue (x=2) in Fig. 2J, K, or L."""
    ya, c0, r0, c1, r1, xt = _calibrate_fig2_jkl(image, y_labels)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    # Find significance bracket strictly below r0 + 15 (horizontal bar spanning between xt[0] and xt[1])
    span_mean = gray[int(r0) + 15:int(r0) + 320, int(xt[0]) + 100:int(xt[1]) - 100].mean(axis=1)
    if float(span_mean.min()) < 200.0:
        r_br = int(np.argmin(span_mean)) + int(r0) + 15
        # Walk down past the vertical bracket legs at both xt[0] and xt[1]
        r_top_start = r_br
        while r_top_start < r_br + 180 and (
            (gray[r_top_start, int(xt[0]) - 8:int(xt[0]) + 9].min() < 180.0)
            or (gray[r_top_start, int(xt[1]) - 8:int(xt[1]) + 9].min() < 180.0)
        ):
            r_top_start += 1
        r_top_start += 10
    else:
        r_top_start = int(r0) + 40

    # If Fig. 2L has horizontal dashed lines, locate their row bands so they are excluded from median/dot detection
    dash_mask = np.zeros(image.shape[:2], dtype=bool)
    if has_dashed_lines:
        bl_dash = (B - R > 60.0) & (G > 140.0) & (R < 130.0)
        bk_dash = (R < 50.0) & (G < 50.0) & (B < 50.0)
        for dm in (bl_dash, bk_dash):
            dr = np.where(dm[r_top_start:int(r1) - 20, int(c0) + 20:int(c0) + 200].any(axis=1))[0] + r_top_start
            if len(dr):
                dash_mask[max(0, dr.min() - 4):min(image.shape[0], dr.max() + 5), :] = True

    out = []
    for gi, (xc, gname) in enumerate(zip(xt, ["black", "blue"]), start=1):
        # Median line: thickest dark horizontal bar (R < 105 works even where translucent blue dot overlaps!)
        # spanning the notch waist [xc-115..xc+115]
        sub_r = R[r_top_start:int(r1) - 20, int(round(xc)) - 115:int(round(xc)) + 115] < 105.0
        if has_dashed_lines:
            sub_r = sub_r & ~dash_mask[r_top_start:int(r1) - 20, int(round(xc)) - 115:int(round(xc)) + 115]
        med_cands = np.where(sub_r.mean(axis=1) > 0.72)[0] + r_top_start
        runs_m = find_runs(med_cands, gap=3)
        if not runs_m:
            raise RuntimeError(f"Fig. 2 {what} ({gname}): could not locate median bar")
        # Thickest horizontal bar is the median (typically 10-16 px thick vs 3-4 px box border)
        best_m = max(runs_m, key=lambda t: t[1] - t[0] + 1)
        r_med = best_m[2]
        med_thick = 0.5 * (best_m[1] - best_m[0] + 1)

        # Whisker stem strictly below bracket along xc +- 5 px
        stem = (gray[r_top_start:int(r1) - 20, int(round(xc)) - 5:int(round(xc)) + 6] < 110.0)
        if has_dashed_lines:
            stem = stem & ~dash_mask[r_top_start:int(r1) - 20, int(round(xc)) - 5:int(round(xc)) + 6]
        st_rows = np.where(stem.any(axis=1))[0] + r_top_start
        r_w_top = float(st_rows.min() + 0.5)
        r_w_bot = float(st_rows.max() + 0.5)

        # Full range including jittered data points (non-white ink within xc +- 125 px, excluding dashed lines)
        nw = (gray[r_top_start:int(r1) - 20, int(round(xc)) - 125:int(round(xc)) + 125] < 235.0)
        if has_dashed_lines:
            nw = nw & ~dash_mask[r_top_start:int(r1) - 20, int(round(xc)) - 125:int(round(xc)) + 125]
        nw_rows = np.where(nw.any(axis=1))[0] + r_top_start
        # Extreme dot circles have radius ~24 px; adjust dot edges inward by 24 px if dots extend past whiskers
        r_d_top = float(nw_rows.min() + 24.0) if nw_rows.min() + 24.0 < r_w_top else r_w_top
        r_d_bot = float(nw_rows.max() - 24.0) if nw_rows.max() - 24.0 > r_w_bot else r_w_bot

        y_med, lo, hi, _, _ = y_interval(ya, r_med, max(med_thick, 2.0))
        d_hi = float(ya.to_value(r_d_top))
        d_lo = float(ya.to_value(r_d_bot))
        w_hi = float(ya.to_value(r_w_top))
        w_lo = float(ya.to_value(r_w_bot))
        h_px = 4.0 + max(ya.residual_px, 0.5)
        bar_hi_outer = float(ya.to_value(r_d_top - h_px))
        bar_lo_outer = float(ya.to_value(r_d_bot + h_px))

        out.append({
            "x": float(gi),
            "y": float(y_med), "lo": float(lo), "hi": float(hi),
            "bar_lo": float(min(d_lo, d_hi)), "bar_hi": float(max(d_lo, d_hi)),
            "bar_lo_outer": float(min(bar_lo_outer, bar_hi_outer)),
            "bar_hi_outer": float(max(bar_lo_outer, bar_hi_outer)),
            "px": [float(xc), float(r_med)],
            "px_lo": float(r_d_top), "px_hi": float(r_d_bot),
            "note": f"{gname} {what}: med={y_med:.4g}, whisk=[{w_lo:.4g},{w_hi:.4g}], range=[{d_lo:.4g},{d_hi:.4g}], ry={ya.residual_px:.2f} px",
        })
    return out


def read_fig2l_adult_lines(image: np.ndarray) -> list[dict]:
    """Read the two horizontal dashed adult cover ridge-number reference lines in Fig. 2L (1=black, 2=blue)."""
    ya, c0, r0, c1, r1, xt = _calibrate_fig2_jkl(image, y_labels=[50.0, 40.0, 30.0])
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    # Scan clean left margin between c0+20 and c0+200 (completely clear of boxplots and points)
    ml, mh = int(round(c0)) + 20, int(round(c0)) + 200
    bk_dash = (R < 50.0) & (G < 50.0) & (B < 50.0)
    bl_dash = (B - R > 60.0) & (G > 140.0) & (R < 130.0)

    out = []
    for idx, (mask, name, xc_mark) in enumerate([(bk_dash, "black adult cover ridge number", xt[0]),
                                                  (bl_dash, "blue adult cover ridge number", xt[1])], start=1):
        rows = np.where(mask[int(r0) + 100:int(r1) - 30, ml:mh].any(axis=1))[0] + int(r0) + 100
        if len(rows) == 0:
            raise RuntimeError(f"Could not find dashed line for {name}")
        rc = float(rows.mean() + 0.5)
        half_px = max(0.5 * (rows.max() - rows.min() + 1), 2.0)
        y_val, lo, hi, px_lo, px_hi = y_interval(ya, rc, half_px)
        out.append({
            "x": float(idx), "y": float(y_val), "lo": float(lo), "hi": float(hi),
            "px": [float(xc_mark), rc], "px_lo": px_lo, "px_hi": px_hi,
            "note": f"{name} dashed line (ry={ya.residual_px:.2f} px)",
        })
    return out
