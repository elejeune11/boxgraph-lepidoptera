"""Shared fast readers for Fig. 1(c), Fig. 1(d), Fig. 2(b), Fig. 2(b) inset, Fig. 2(c), and Fig. 2(d)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _helpers as H


def _find_all_inner_frames(img: np.ndarray, box: dict, min_w_frac: float = 0.10, min_h_frac: float = 0.08) -> list[tuple[int, int, int, int]]:
    """Detect all dark-bordered rectangular boxes (legends / insets) inside the plot box."""
    r0i, r1i = box["r0_e"] + 15, box["r1_s"] - 15
    c0i, c1i = box["c0_e"] + 15, box["c1_s"] - 15
    step = 4
    sub = img[r0i:r1i:step, c0i:c1i:step]
    dark = sub.max(axis=2) < 110
    H_s, W_s = dark.shape
    min_w = max(8, int(min_w_frac * W_s))
    min_h = max(8, int(min_h_frac * H_s))
    frame = H.open_rect(dark, 1, min_w) | H.open_rect(dark, min_h, 1)
    boxes = []
    for _, rsl, csl, _ in H.labeled_components(frame):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if w >= min_w and h >= min_h:
            boxes.append((
                c0i + csl.start * step - 10,
                r0i + rsl.start * step - 10,
                c0i + csl.stop * step + 10,
                r0i + rsl.stop * step + 10,
            ))
    return boxes


def _read_hollow_marker_at_col(
    mask: np.ndarray,
    c_target: float,
    r_min: int,
    r_max: int,
    xa: H.AxisCal,
    ya: H.AxisCal,
    half_win: int = 85,
    min_body_span: int = 48,
    min_body_h: int = 24,
    is_triangle: bool = False,
    is_square: bool = False,
    note: str = "",
    force_x: float | None = None,
    max_bar_gap: int = 65,
) -> dict | None:
    """Locate one hollow marker (square/circle/triangle) + its vertical error bar near column c_target.

    Works reliably even when markers of different colors overlap at the same x position.
    """
    H_im, W_im = mask.shape
    c0 = max(0, int(round(c_target)) - half_win)
    c1 = min(W_im, int(round(c_target)) + half_win + 1)
    strip = mask[r_min:r_max, c0:c1]
    if strip.sum() < 120:
        return None

    spans = np.zeros(r_max - r_min, dtype=int)
    c_mids = np.full(r_max - r_min, np.nan)
    for idx in range(r_max - r_min):
        cs = np.nonzero(strip[idx])[0]
        if len(cs) >= 2:
            spans[idx] = int(cs[-1] - cs[0] + 1)
            c_mids[idx] = c0 + 0.5 * (cs[0] + cs[-1]) + 0.5
        elif len(cs) == 1:
            spans[idx] = 1
            c_mids[idx] = c0 + float(cs[0]) + 0.5

    wide_rows = np.where(spans >= min_body_span)[0]
    clusters = H.clusters_1d(wide_rows, max_gap=18)
    body_cands = [cl for cl in clusters if (cl[1] - cl[0] + 1) >= min_body_h]
    if not body_cands:
        wide_rows = np.where(spans >= max(30, min_body_span - 14))[0]
        clusters = H.clusters_1d(wide_rows, max_gap=18)
        body_cands = [cl for cl in clusters if (cl[1] - cl[0] + 1) >= min_body_h]
    if not body_cands:
        return None
    bs, be, _ = max(body_cands, key=lambda cl: (cl[1] - cl[0] + 1) * float(spans[cl[0]:cl[1] + 1].max()))
    r_top = r_min + bs
    r_bot = r_min + be
    max_w = int(spans[bs:be + 1].max())
    if is_triangle:
        max_apex_rise = int(round(0.45 * max_w))
        for k in range(1, max_apex_rise + 1):
            ri = bs - k
            if ri < 0:
                break
            if spans[ri] >= 12 and spans[ri] <= spans[min(bs, ri + 4)] + 8:
                r_top = r_min + ri
            else:
                break
    elif is_square and (r_bot - r_top + 1) < int(round(0.82 * max_w)):
        # Bottom half of square is occluded by an overlaid marker (e.g. red triangle/black circle):
        # since a square has height == width, extend r_bot to r_top + max_w - 1
        r_bot = r_top + max_w - 1

    c_vals = c_mids[bs:be + 1]
    c_vals = c_vals[~np.isnan(c_vals)]
    c_m = float(np.median(c_vals)) if len(c_vals) else float(c_target)
    r_m = 0.5 * (r_top + r_bot) + 0.5

    # Trace vertical error bar upward from r_top and downward from r_bot in narrow strip around c_m
    r_up, r_dn = H.trace_vertical_ebar(
        mask, c_m, r_top - 1, r_bot + 1,
        r_min_bound=r_min, r_max_bound=r_max - 1,
        col_half_win=22, max_gap=max_bar_gap,
    )
    if r_up is not None and (r_top - r_up) < 5.0:
        r_up = None
    if r_dn is not None and (r_dn - r_bot) < 5.0:
        r_dn = None
    # Note: if an overlapping marker of another color hides one or both short ends inside this marker
    has_bar = True
    full_note = note
    if r_up is None and r_dn is None:
        r_up, r_dn = float(r_top), float(r_bot)
        full_note += " (error bar hidden inside marker)"
    elif r_up is None:
        r_up = float(r_top)
        full_note += " (upper bar end hidden inside marker)"
    elif r_dn is None:
        r_dn = float(r_bot)
        full_note += " (lower bar end hidden inside marker)"

    return H.make_marker_reading(
        xa, ya, c_m, r_m, r_up, r_dn,
        marker_half_h=0.5 * (r_bot - r_top + 1),
        note=full_note, force_x=force_x, include_ebar=has_bar,
    )


# ---------------------------------------------------------------------------
# Fig. 1(c): Normalized Intensity vs Distance (um) — Actin (red) & FUS (black)
# ---------------------------------------------------------------------------
def _cal_fig1c(img: np.ndarray):
    box = H.find_box(img)
    xt = H.find_ticks(img, box, "top", min_len=12)
    if len(xt) != 3:
        raise ValueError(f"fig1c: expected 3 top x ticks (2,4,6), got {len(xt)}: {xt}")
    xpairs = [(box["c0"], 0.0), (xt[0][0], 2.0), (xt[1][0], 4.0), (xt[2][0], 6.0)]
    lt = H.find_ticks(img, box, "left", min_len=12)
    lt_clean = [t for t in lt if t[0] < box["r1_s"] - 120]
    if len(lt_clean) != 4:
        raise ValueError(f"fig1c: expected 4 left y ticks, got {len(lt_clean)}: {lt_clean}")
    ypairs = [(lt_clean[0][0], 1.0), (lt_clean[1][0], 0.75), (lt_clean[2][0], 0.5), (lt_clean[3][0], 0.25), (box["r1"], 0.0)]
    return box, H.AxisCal(xpairs, max_res_px=3.0), H.AxisCal(ypairs, max_res_px=3.0)


def read_fig1c(img: np.ndarray, series: str) -> list[dict]:
    """Read Fig. 1(c) midplane fluorescence intensity profile ('actin' red triangles or 'fus' black circles)."""
    box, xa, ya = _cal_fig1c(img)
    valid = np.zeros(img.shape[:2], dtype=bool)
    # Valid data sits below y ~ 1.04 and stops well inside the right spine ticks (data ends at x ~ 8.3 um)
    r_y_top = max(box["r0_e"] + 10, int(round(ya.to_pixel(1.04))))
    c_x_max = min(box["c1_s"] - 50, int(round(xa.to_pixel(8.5))))
    valid[r_y_top:box["r1_s"] - 2, box["c0_e"] + 10:c_x_max] = True
    # Exclude left-spine inward ticks (44 px long) at y = 0.25, 0.5, 0.75, 1.0
    valid[r_y_top:box["r1_s"] - 65, box["c0_e"]:box["c0_e"] + 50] = False
    # Exclude ALL inner boxes (both upper-left legend and upper-right inset micrograph)
    for c0, r0, c1, r1 in _find_all_inner_frames(img, box, min_w_frac=0.10, min_h_frac=0.08):
        valid[max(0, r0):r1 + 1, max(0, c0):c1 + 1] = False

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    if series == "actin":
        mask = valid & (r_ch > 150) & (r_ch - g_ch > 80) & (r_ch - b_ch > 80)
    elif series == "fus":
        mask = valid & (mx < 85) & ((mx - mn) < 30)
    else:
        raise ValueError(f"Unknown fig1c series {series!r}")

    opened = H.open_rect(mask, 8, 8)
    out = []
    for _, rsl, csl, loc in H.labeled_components(opened):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        area = int(loc.sum())
        if area < 120:
            continue
        if w <= 56 and h <= 56:
            rr, cc = np.nonzero(loc)
            c_m = csl.start + float(cc.mean()) + 0.5
            r_m = rsl.start + float(rr.mean()) + 0.5
            out.append(H.make_marker_reading(xa, ya, c_m, r_m, None, None, 0.5 * h, f"{series} marker", include_ebar=False))
        else:
            k = max(2, int(round(w / 34.0)))
            edges = np.linspace(0, w, k + 1).astype(int)
            for j in range(k):
                sub = loc[:, edges[j]:edges[j + 1]]
                if sub.sum() < 60:
                    continue
                rr, cc = np.nonzero(sub)
                c_m = csl.start + edges[j] + float(cc.mean()) + 0.5
                r_m = rsl.start + float(rr.mean()) + 0.5
                out.append(H.make_marker_reading(xa, ya, c_m, r_m, None, None, 14.0, f"{series} marker (split {j+1}/{k})", include_ebar=False))
    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 1(d): I_inside / I_outside vs [Actin] (uM)
# ---------------------------------------------------------------------------
def _cal_fig1d(img: np.ndarray):
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    major_x = [t[0] for t in tt if t[1] >= 34]
    if len(major_x) != 3:
        raise ValueError(f"fig1d: expected 3 major top x ticks (0.01, 0.1, 1), got {len(major_x)}: {major_x}")
    xa = H.AxisCal([(major_x[0], 0.01), (major_x[1], 0.1), (major_x[2], 1.0)], scale="log10", max_res_px=2.5)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 8:
        raise ValueError(f"fig1d: expected 8 right y ticks, got {len(rt)}: {rt}")
    ypairs = [(rt[i][0], 80.0 - 10.0 * i) for i in range(8)] + [(box["r1"], 0.0)]
    ya = H.AxisCal(ypairs, max_res_px=2.5)
    return box, xa, ya, major_x


def read_fig1d(img: np.ndarray, series: str) -> list[dict]:
    """Read Fig. 1(d): '90nm' (blue triangles), '135nm' (red circles), '270nm' (black squares)."""
    box, xa, ya, major_x = _cal_fig1d(img)
    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    r0_i, r1_i = box["r0_e"] + 46, box["r1_s"] - 46
    c0_i, c1_i = box["c0_e"] + 46, box["c1_s"] - 46
    valid = np.zeros(img.shape[:2], dtype=bool)
    valid[r0_i:r1_i, c0_i:c1_i] = True
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        valid[lr0:lr1 + 1, lc0:lc1 + 1] = False

    if series == "90nm":
        mask = valid & (b_ch > 105) & (b_ch - r_ch > 45) & (b_ch - g_ch > 15)
        min_span, is_tri, bar_gap = 56, True, 65
    elif series == "135nm":
        mask = valid & (r_ch > 150) & (r_ch - g_ch > 80) & (r_ch - b_ch > 80)
        min_span, is_tri, bar_gap = 62, False, 65
    elif series == "270nm":
        mask = valid & (mx < 95) & ((mx - mn) < 30)
        min_span, is_tri, bar_gap = 58, False, 20
    else:
        raise ValueError(f"Unknown fig1d series {series!r}")

    out = []
    for c_tick in major_x:
        rd = _read_hollow_marker_at_col(
            mask, c_tick, r0_i, r1_i, xa, ya,
            half_win=95, min_body_span=min_span, min_body_h=28,
            is_triangle=is_tri, note=f"L~{series} ±1 SD", max_bar_gap=bar_gap,
        )
        if rd is not None:
            out.append(rd)
    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 2(b) main: Normalized Area vs Rescaled Time — single exponential e^(-t/tau)
# ---------------------------------------------------------------------------
def _cal_fig2b(img: np.ndarray):
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    if len(tt) != 4:
        raise ValueError(f"fig2b: expected 4 top x ticks (1,2,3,4), got {len(tt)}: {tt}")
    xpairs = [(box["c0"], 0.0), (tt[0][0], 1.0), (tt[1][0], 2.0), (tt[2][0], 3.0), (tt[3][0], 4.0)]
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 5:
        raise ValueError(f"fig2b: expected 5 right y ticks, got {len(rt)}: {rt}")
    ypairs = [(rt[0][0], 1.0), (rt[1][0], 0.75), (rt[2][0], 0.5), (rt[3][0], 0.25), (rt[4][0], 0.0)]
    return box, H.AxisCal(xpairs, max_res_px=3.0), H.AxisCal(ypairs, max_res_px=3.0)


def read_fig2b_exp_curve(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read the black dashed single-exponential master curve e^(-t/tau) in Fig. 2(b)."""
    box, xa, ya = _cal_fig2b(img)
    mx = img.max(axis=2).astype(np.int16)
    mn = img.min(axis=2).astype(np.int16)
    black = (mx < 75) & ((mx - mn) < 25)
    leg = H.find_inner_frame_box(img, box, min_w_frac=0.25, min_h_frac=0.25)
    if leg is not None:
        c0, r0, c1, r1 = leg
        black[max(0, r0 - 60):r1 + 190, max(0, c0 - 260):c1 + 60] = False
    black[:box["r0_e"] + 8, :] = False
    black[box["r1_s"] - 8:, :] = False
    black[:, :box["c0_e"] + 8] = False
    black[:, box["c1_s"] - 8:] = False
    # Exclude small compact scatter squares (~18x18 px) and any stray dots far from the dashed line
    filtered = np.zeros_like(black)
    for _, rsl, csl, loc in H.labeled_components(black):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        area = int(loc.sum())
        fill = area / max(1, h * w)
        if area < 40:
            continue
        if max(h, w) <= 24 and fill > 0.68:
            continue
        # At t/tau > 2.5 the dashed curve is nearly horizontal (h <= 22, w >= 18)
        c_mid = 0.5 * (csl.start + csl.stop)
        if xa.to_value(c_mid) > 2.5 and h > 25:
            continue
        filtered[rsl, csl] |= loc
    return H.sample_curve_at(filtered, xa, ya, at, box, fit_half_win=85, max_gap_px=180)


# ---------------------------------------------------------------------------
# Fig. 2(b) inset: Area (um^2) vs Time (s) — green triangles & black dashed fit
# ---------------------------------------------------------------------------
def _cal_fig2b_inset(img: np.ndarray):
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=6, max_width=40)
    if len(tt) != 7:
        raise ValueError(f"fig2b_inset: expected 7 top x ticks (0..60), got {len(tt)}: {tt}")
    xpairs = [(tt[i][0], 10.0 * i) for i in range(7)]
    rt = H.find_ticks(img, box, "right", min_len=6, max_width=40)
    if len(rt) != 9:
        raise ValueError(f"fig2b_inset: expected 9 right y ticks (24..8), got {len(rt)}: {rt}")
    ypairs = [(rt[i][0], 24.0 - 2.0 * i) for i in range(9)]
    # At 2400 dpi (2x normal 1200 dpi), 5.0 px corresponds to ~0.2 s (< 0.3% of axis)
    return box, H.AxisCal(xpairs, max_res_px=5.0), H.AxisCal(ypairs, max_res_px=5.0)


def read_fig2b_inset(img: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read Fig. 2(b) inset for t >= 0: which='data' (green triangles) or 'fit' (black dashed curve)."""
    box, xa, ya = _cal_fig2b_inset(img)
    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    c_zero = int(round(xa.to_pixel(0.0)))
    valid = np.zeros(img.shape[:2], dtype=bool)
    valid[box["r0_e"] + 15:box["r1_s"] - 15, max(box["c0_e"] + 10, c_zero - 8):box["c1_s"] - 15] = True
    if which == "data":
        mask = valid & (g_ch > 105) & (g_ch - r_ch > 45) & (g_ch - b_ch > 30)
    elif which == "fit":
        mask = valid & (mx < 70) & ((mx - mn) < 25)
    else:
        raise ValueError(f"Unknown fig2b_inset mode {which!r}")
    return H.sample_curve_at(mask, xa, ya, at, box, fit_half_win=70, max_gap_px=160)


# ---------------------------------------------------------------------------
# Fig. 2(c): tau (s) vs sqrt(A_i) - sqrt(A_f) (um)
# ---------------------------------------------------------------------------
def _cal_fig2c(img: np.ndarray):
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    if len(tt) != 5:
        raise ValueError(f"fig2c: expected 5 top x ticks (0.2..0.6), got {len(tt)}: {tt}")
    xpairs = [(tt[i][0], 0.2 + 0.1 * i) for i in range(5)]
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 4:
        raise ValueError(f"fig2c: expected 4 right y ticks (20..5), got {len(rt)}: {rt}")
    ypairs = [(rt[i][0], 20.0 - 5.0 * i) for i in range(4)] + [(box["r1"], 0.0)]
    return box, H.AxisCal(xpairs, max_res_px=2.5), H.AxisCal(ypairs, max_res_px=2.5)


def _fig2c_squares_mask(img: np.ndarray, box: dict):
    dark = img.max(axis=2) < 95
    m = np.zeros_like(dark)
    m[box["r0_e"] + 50:box["r1_s"] - 50, box["c0_e"] + 50:box["c1_s"] - 50] = dark[box["r0_e"] + 50:box["r1_s"] - 50, box["c0_e"] + 50:box["c1_s"] - 50]
    filled = ndimage.binary_fill_holes(m)
    sq = H.open_rect(filled, 45, 45)
    return m, sq


def read_fig2c_data(img: np.ndarray) -> list[dict]:
    """Read the 9 open black squares (and their small vertical error bars) in Fig. 2(c)."""
    box, xa, ya = _cal_fig2c(img)
    m, sq = _fig2c_squares_mask(img, box)
    out = []
    for _, rsl, csl, loc in H.labeled_components(sq):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if not (55 <= w <= 110 and 55 <= h <= 110):
            continue
        rr, cc = np.nonzero(loc)
        c_m = csl.start + float(cc.mean()) + 0.5
        r_m = rsl.start + float(rr.mean()) + 0.5
        r_up, r_dn = H.trace_vertical_ebar(
            m, c_m, rsl.start - 2, rsl.stop + 1,
            r_min_bound=max(box["r0_e"] + 50, rsl.start - 60),
            r_max_bound=min(box["r1_s"] - 50, rsl.stop + 60),
            col_half_win=6, max_gap=8,
        )
        if r_up is None or r_dn is None:
            ri0, ri1 = rsl.start + 12, rsl.stop - 12
            ci0, ci1 = int(round(c_m)) - 12, int(round(c_m)) + 13
            inner = m[ri0:ri1, ci0:ci1]
            row_hits = np.where(inner.sum(axis=1) >= 10)[0]
            if len(row_hits) > 0:
                if r_up is None:
                    r_up = ri0 + float(row_hits[0]) + 0.5
                if r_dn is None:
                    r_dn = ri0 + float(row_hits[-1]) + 0.5
        has_bar = (r_up is not None) and (r_dn is not None)
        note = "coalescence event square with small unstated vertical error bar" if has_bar else "coalescence event square"
        out.append(H.make_marker_reading(xa, ya, c_m, r_m, r_up, r_dn, 0.5 * h, note, include_ebar=has_bar))
    return sorted(out, key=lambda d: d["x"])


def read_fig2c_fit(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read the black dashed linear fit line in Fig. 2(c)."""
    box, xa, ya = _cal_fig2c(img)
    m, sq = _fig2c_squares_mask(img, box)
    line_m = m.copy()
    for _, rsl, csl, _ in H.labeled_components(sq):
        line_m[max(0, rsl.start - 35):rsl.stop + 35, max(0, csl.start - 35):csl.stop + 35] = False
    return H.sample_curve_at(line_m, xa, ya, at, box, fit_half_win=120, max_gap_px=260)


# ---------------------------------------------------------------------------
# Fig. 2(d): eta/gamma (s/um) vs [Actin] (uM) — broken x-axis (0 // 0.01..1)
# ---------------------------------------------------------------------------
def read_fig2d(img: np.ndarray) -> list[dict]:
    """Read the 4 filled black squares and SEM error bars in Fig. 2(d)."""
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    major_x = [t[0] for t in tt if t[1] >= 36]
    if len(major_x) != 4:
        raise ValueError(f"fig2d: expected 4 major top x ticks (0, 0.01, 0.1, 1), got {len(major_x)}: {major_x}")
    c_zero = major_x[0]
    xa_log = H.AxisCal([(major_x[1], 0.01), (major_x[2], 0.1), (major_x[3], 1.0)], scale="log10", max_res_px=3.0)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 7:
        raise ValueError(f"fig2d: expected 7 right y ticks (30..0), got {len(rt)}: {rt}")
    ya = H.AxisCal([(rt[i][0], 30.0 - 5.0 * i) for i in range(7)], max_res_px=3.5)

    mx = img.max(axis=2)
    m = np.zeros_like(mx < 95)
    r0_i, r1_i = box["r0_e"] + 8, box["r1_s"] - 8
    c0_i, c1_i = box["c0_e"] + 8, box["c1_s"] - 8
    m[r0_i:r1_i, c0_i:c1_i] = (mx[r0_i:r1_i, c0_i:c1_i] < 95)

    sq = H.open_rect(m, 28, 28)
    out = []
    for _, rsl, csl, loc in H.labeled_components(sq):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if not (30 <= w <= 85 and 30 <= h <= 85):
            continue
        rr, cc = np.nonzero(loc)
        c_m = csl.start + float(cc.mean()) + 0.5
        r_m = rsl.start + float(rr.mean()) + 0.5
        r_up, r_dn = H.trace_vertical_ebar(
            m, c_m, rsl.start - 2, rsl.stop + 1,
            r_min_bound=r0_i, r_max_bound=r1_i,
            col_half_win=6, max_gap=10,
        )
        if abs(c_m - c_zero) < 40.0:
            x_val = 0.0
            note = "[Actin]=0 uM (left of axis break), SEM from linear fit"
        else:
            x_val = xa_log.to_value(c_m)
            note = "log10-calibrated [Actin], SEM from linear fit"
        out.append(H.make_marker_reading(
            xa_log, ya, c_m, r_m, r_up, r_dn, 0.5 * h, note,
            force_x=x_val, include_ebar=(r_up is not None or r_dn is not None),
        ))
    return sorted(out, key=lambda d: d["x"])
