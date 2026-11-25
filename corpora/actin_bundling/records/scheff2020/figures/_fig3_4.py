"""Shared fast readers for Figs. 3(c), 3(d), 3(f), 3(g), 4(b), and 4(c)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _helpers as H
from _fig1_2 import _read_hollow_marker_at_col


# ---------------------------------------------------------------------------
# Shared dense black-square + vertical-stem extractor for Figs. 3(c) & 3(f)
# ---------------------------------------------------------------------------
def _extract_stemmed_and_isolated_squares(
    img: np.ndarray,
    box: dict,
    xa: H.AxisCal,
    ya: H.AxisCal,
    y_text_cut: float,
    sq_k: int,
    sq_w: float,
    stem_h: int,
) -> list[dict]:
    """Extract dense binned black squares (separated by vertical error-bar stems) plus isolated squares."""
    r0_i, r1_i = box["r0_e"] + 12, box["r1_s"] - 10
    c0_i, c1_i = box["c0_e"] + 12, box["c1_s"] - 10
    black = np.zeros(img.shape[:2], dtype=bool)
    black[r0_i:r1_i, c0_i:c1_i] = img[r0_i:r1_i, c0_i:c1_i].max(axis=2) < 95

    sq_blk = H.open_rect(black, sq_k, sq_k)
    r_text = int(round(ya.to_pixel(y_text_cut)))
    sq_blk[:max(0, r_text), :] = False
    # Exclude any inward ticks on spines that might survive small opening
    sq_blk[:box["r0_e"] + 46, :] = False
    sq_blk[:, :box["c0_e"] + 46] = False
    sq_blk[:, box["c1_s"] - 46:] = False

    # Isolate vertical stems outside the square bodies
    stems_v = H.open_rect(black & ~ndimage.binary_dilation(sq_blk, iterations=2), stem_h, 1)
    clean_stems = np.zeros_like(stems_v)
    sq_col_has = sq_blk.any(axis=0)
    for _, rsl, csl, loc in H.labeled_components(stems_v):
        if (csl.stop - csl.start) <= 12 and sq_col_has[max(0, csl.start - 4):csl.stop + 4].any():
            sq_rows = np.nonzero(sq_blk[:, max(0, csl.start - 4):csl.stop + 4].any(axis=1))[0]
            if len(sq_rows) and (rsl.stop >= sq_rows.min() - 16 and rsl.start <= sq_rows.max() + 16):
                clean_stems[rsl, csl] |= loc

    col_hits = np.where(clean_stems.sum(axis=0) >= stem_h)[0]
    stem_cols = H.clusters_1d(col_hits, max_gap=2)

    out = []
    used_cols = []
    # 1. Read every stemmed bin centered at each vertical stem column
    for sc_s, sc_e, c_m in stem_cols:
        c_idx = int(round(c_m))
        # Find the solid square body at this column in sq_blk (search +-3 px around c_idx)
        rs_sq = np.nonzero(sq_blk[:, max(0, c_idx - 3):c_idx + 4].any(axis=1))[0]
        if len(rs_sq) == 0:
            continue
        # If multiple vertical runs exist (rare), pick the run closest to the stem midpoint
        runs = H.clusters_1d(rs_sq, max_gap=3)
        st_rows = np.nonzero(clean_stems[:, max(0, sc_s - 1):sc_e + 2].any(axis=1))[0]
        st_mid = float(st_rows.mean()) if len(st_rows) else runs[0][2]
        r_top, r_bot, r_m = min(runs, key=lambda rn: abs(rn[2] - st_mid))
        # Clamp marker height if two stacked squares barely merged
        half_h = min(0.5 * sq_w, 0.5 * (r_bot - r_top + 1))

        # Trace vertical error bar strictly along +-3 px of c_m
        r_up, r_dn = H.trace_vertical_ebar(
            black, c_m, r_top - 2, r_bot + 2,
            r_min_bound=r0_i, r_max_bound=r1_i,
            col_half_win=3, max_gap=14,
        )
        if r_up is not None and (r_top - r_up) < 8.0:
            r_up = None
        if r_dn is not None and (r_dn - r_bot) < 8.0:
            r_dn = None
        has_bar = (r_up is not None) or (r_dn is not None)
        note = "binned mean ±1 SD" if has_bar else "binned mean"
        if r_up is not None and r_dn is None and (r1_i - r_bot) < (r_top - r_up):
            r_dn = float(r1_i)
            note += " (lower end clipped at bottom axis)"
        out.append(H.make_marker_reading(xa, ya, c_m, r_m, r_up, r_dn, half_h, note, include_ebar=has_bar))
        used_cols.append(c_m)

    # 2. Read isolated squares (typically single-droplet bins at larger area with no error bars)
    for _, rsl, csl, loc in H.labeled_components(sq_blk):
        # Skip any blob that already contains a detected stem column
        if any((csl.start - 6) <= uc <= (csl.stop + 6) for uc in used_cols):
            continue
        w = csl.stop - csl.start
        h = rsl.stop - rsl.start
        if w < sq_k or h < sq_k:
            continue
        if w <= 1.55 * sq_w and h <= 1.55 * sq_w:
            rr, cc = np.nonzero(loc)
            c_m = csl.start + float(cc.mean()) + 0.5
            r_m = rsl.start + float(rr.mean()) + 0.5
            out.append(H.make_marker_reading(xa, ya, c_m, r_m, None, None, 0.5 * h, "single-droplet bin (no error bar)", include_ebar=False))
        else:
            # Split horizontally touching isolated squares (e.g. 36 and 37 um^2 in Fig. 3f)
            k = max(2, int(round(w / sq_w)))
            edges = np.linspace(0, w, k + 1).astype(int)
            for j in range(k):
                sub = loc[:, edges[j]:edges[j + 1]]
                if sub.sum() < 0.3 * sq_k * sq_k:
                    continue
                rr, cc = np.nonzero(sub)
                c_m = csl.start + edges[j] + float(cc.mean()) + 0.5
                r_m = rsl.start + float(rr.mean()) + 0.5
                out.append(H.make_marker_reading(xa, ya, c_m, r_m, None, None, 0.5 * sq_w, "single-droplet bin (split)", include_ebar=False))

    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 3(c): Aspect Ratio vs Area (1 uM actin, L ~ 270 nm, n=1180)
# ---------------------------------------------------------------------------
def read_fig3c(img: np.ndarray) -> list[dict]:
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    if len(tt) != 9:
        raise ValueError(f"fig3c: expected 9 top x ticks (0..120), got {len(tt)}: {tt}")
    xa = H.AxisCal([(tt[i][0], 15.0 * i) for i in range(9)], max_res_px=4.0)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 6:
        raise ValueError(f"fig3c: expected 6 right y ticks (2.25..1.0), got {len(rt)}: {rt}")
    ya = H.AxisCal([(box["r0"], 2.5)] + [(rt[i][0], 2.25 - 0.25 * i) for i in range(6)], max_res_px=4.0)
    return _extract_stemmed_and_isolated_squares(img, box, xa, ya, y_text_cut=2.10, sq_k=17, sq_w=28.0, stem_h=18)


# ---------------------------------------------------------------------------
# Fig. 3(f): Aspect Ratio vs Area (1 uM actin, L ~ 135 nm, n=1308)
# ---------------------------------------------------------------------------
def read_fig3f(img: np.ndarray) -> list[dict]:
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    if len(tt) != 5:
        raise ValueError(f"fig3f: expected 5 top x ticks (0..40), got {len(tt)}: {tt}")
    xa = H.AxisCal([(tt[i][0], 10.0 * i) for i in range(5)], max_res_px=4.0)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 5:
        raise ValueError(f"fig3f: expected 5 right y ticks (2.0..1.0), got {len(rt)}: {rt}")
    ya = H.AxisCal([(rt[i][0], 2.0 - 0.25 * i) for i in range(5)], max_res_px=4.0)
    return _extract_stemmed_and_isolated_squares(img, box, xa, ya, y_text_cut=1.70, sq_k=24, sq_w=36.0, stem_h=10)


def _symmetrize_sd_if_occluded(rd: dict, ya: H.AxisCal, max_asym_ratio: float = 0.35) -> dict:
    """For explicitly symmetric ±1 SD error bars in Figs. 3(d)/(g), if one arm is truncated
    by an overlapping marker/bar of another series (asymmetry > 35% and longer arm extends clearly
    beyond the marker), mirror the longer unobstructed arm."""
    if "bar_lo" not in rd or "bar_hi" not in rd:
        return rd
    if "error bar hidden inside marker" in rd.get("note", ""):
        return rd
    y = rd["y"]
    r_m = rd["px"][1]
    r_top, r_bot = sorted((rd["px_lo"], rd["px_hi"]))
    up_px = r_m - r_top
    dn_px = r_bot - r_m
    longer_px = max(up_px, dn_px)
    shorter_px = min(up_px, dn_px)
    if longer_px > 55.0 and (longer_px - shorter_px) / longer_px > max_asym_ratio:
        half_px = longer_px
        pad = 2.0 + ya.residual_px
        b1, b2 = ya.to_value(r_m - half_px), ya.to_value(r_m + half_px)
        bo1, bo2 = ya.to_value(r_m - half_px - pad), ya.to_value(r_m + half_px + pad)
        rd["bar_lo"] = min(b1, b2, y)
        rd["bar_hi"] = max(b1, b2, y)
        rd["bar_lo_outer"] = min(bo1, bo2, rd["lo"])
        rd["bar_hi_outer"] = max(bo1, bo2, rd["hi"])
        rd["px_lo"] = float(r_m - half_px)
        rd["px_hi"] = float(r_m + half_px)
        rd["note"] = rd.get("note", "") + " (occluded arm mirrored symmetrically)"
    return rd


# ---------------------------------------------------------------------------
# Fig. 3(d): Aspect Ratio vs [Actin] (uM), L ~ 270 nm (3 droplet size bins)
# ---------------------------------------------------------------------------
def read_fig3d(img: np.ndarray, series: str) -> list[dict]:
    """series: 'small' (4-6 um^2 blue squares), 'medium' (14-16 um^2 red triangles), 'large' (>20 um^2 black circles)."""
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    major_x = [t[0] for t in tt if t[1] >= 30]
    if len(major_x) != 3:
        raise ValueError(f"fig3d: expected 3 major top x ticks (0.01, 0.1, 1), got {len(major_x)}: {major_x}")
    xa = H.AxisCal([(major_x[0], 0.01), (major_x[1], 0.1), (major_x[2], 1.0)], scale="log10", max_res_px=3.0)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 4:
        raise ValueError(f"fig3d: expected 4 right y ticks (1.6..1.0), got {len(rt)}: {rt}")
    ya = H.AxisCal([(box["r0"], 1.8)] + [(rt[i][0], 1.6 - 0.2 * i) for i in range(4)], max_res_px=3.5)

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    r0_i, r1_i = box["r0_e"] + 40, box["r1_s"] - 40
    c0_i, c1_i = box["c0_e"] + 40, box["c1_s"] - 40
    valid = np.zeros(img.shape[:2], dtype=bool)
    valid[r0_i:r1_i, c0_i:c1_i] = True
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        valid[lr0:lr1 + 1, lc0:lc1 + 1] = False

    if series == "small":
        mask = valid & (b_ch > 105) & (b_ch - r_ch > 45) & (b_ch - g_ch > 15)
        min_span, is_tri, is_sq, bar_gap, label = 54, False, True, 80, "4-6 um^2 ±1 SD"
    elif series == "medium":
        mask = valid & (r_ch > 150) & (r_ch - g_ch > 80) & (r_ch - b_ch > 80)
        min_span, is_tri, is_sq, bar_gap, label = 50, True, False, 80, "14-16 um^2 ±1 SD"
    elif series == "large":
        mask = valid & (mx < 95) & ((mx - mn) < 30)
        min_span, is_tri, is_sq, bar_gap, label = 54, False, False, 80, ">20 um^2 ±1 SD"
    else:
        raise ValueError(f"Unknown fig3d series {series!r}")

    out = []
    for c_tick in major_x:
        rd = _read_hollow_marker_at_col(
            mask, c_tick, r0_i, r1_i, xa, ya,
            half_win=90, min_body_span=min_span, min_body_h=24,
            is_triangle=is_tri, is_square=is_sq, note=label, max_bar_gap=bar_gap,
        )
        if rd is not None:
            out.append(_symmetrize_sd_if_occluded(rd, ya))
    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 3(g): Aspect Ratio vs Approx. Filament Length (nm), 1 uM actin (broken axis)
# ---------------------------------------------------------------------------
def read_fig3g(img: np.ndarray, series: str) -> list[dict]:
    """series: 'small' (4-6 um^2), 'medium' (14-16 um^2), 'large' (>20 um^2). Lat plotted at x=0 nm."""
    box = H.find_box(img)
    tt = H.find_ticks(img, box, "top", min_len=10)
    if len(tt) != 7:
        raise ValueError(f"fig3g: expected 7 top x ticks (Lat, //, 100..250), got {len(tt)}: {tt}")
    c_lat = tt[0][0]
    xa = H.AxisCal([(tt[3][0], 100.0), (tt[4][0], 150.0), (tt[5][0], 200.0), (tt[6][0], 250.0), (box["c1"], 300.0)], max_res_px=3.5)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 4:
        raise ValueError(f"fig3g: expected 4 right y ticks (2.5..1.0), got {len(rt)}: {rt}")
    ya = H.AxisCal([(box["r0"], 3.0), (rt[0][0], 2.5), (rt[1][0], 2.0), (rt[2][0], 1.5), (rt[3][0], 1.0)], max_res_px=3.5)

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    r0_i, r1_i = box["r0_e"] + 42, box["r1_s"] - 42
    c0_i, c1_i = box["c0_e"] + 42, box["c1_s"] - 42
    valid = np.zeros(img.shape[:2], dtype=bool)
    valid[r0_i:r1_i, c0_i:c1_i] = True
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        valid[lr0:lr1 + 1, lc0:lc1 + 1] = False

    if series == "small":
        mask = valid & (b_ch > 105) & (b_ch - r_ch > 45) & (b_ch - g_ch > 15)
        min_span, is_tri, is_sq, bar_gap, label = 54, False, True, 85, "4-6 um^2 ±1 SD"
    elif series == "medium":
        mask = valid & (r_ch > 150) & (r_ch - g_ch > 80) & (r_ch - b_ch > 80)
        min_span, is_tri, is_sq, bar_gap, label = 50, True, False, 85, "14-16 um^2 ±1 SD"
    elif series == "large":
        mask = valid & (mx < 95) & ((mx - mn) < 30)
        min_span, is_tri, is_sq, bar_gap, label = 54, False, False, 85, ">20 um^2 ±1 SD"
    else:
        raise ValueError(f"Unknown fig3g series {series!r}")

    col_targets = [
        (c_lat, 0.0, label + " (Lat = monomeric actin before axis break)"),
        (xa.to_pixel(90.0), None, label),
        (xa.to_pixel(135.0), None, label),
        (xa.to_pixel(270.0), None, label),
    ]
    out = []
    for c_t, fx, nt in col_targets:
        rd = _read_hollow_marker_at_col(
            mask, c_t, r0_i, r1_i, xa, ya,
            half_win=85, min_body_span=min_span, min_body_h=22,
            is_triangle=is_tri, is_square=is_sq, note=nt, force_x=fx, max_bar_gap=bar_gap,
        )
        if rd is not None:
            out.append(_symmetrize_sd_if_occluded(rd, ya))
    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 4(b): Fraction Tactoid (%) vs Area (um^2) — variable-width histogram bars
# ---------------------------------------------------------------------------
def read_fig4b(img: np.ndarray) -> list[dict]:
    box = H.find_box(img)
    bt = H.find_ticks(img, box, "bottom", outward=True, min_len=8, max_len=60)
    major_x = [t[0] for t in bt if t[1] >= 20]
    if len(major_x) != 4:
        raise ValueError(f"fig4b: expected 4 major bottom outward x ticks (5,10,15,20), got {len(major_x)}: {major_x}")
    xa = H.AxisCal([(major_x[0], 5.0), (major_x[1], 10.0), (major_x[2], 15.0), (major_x[3], 20.0), (box["c1"], 25.0)], max_res_px=4.5)
    lt = H.find_ticks(img, box, "left", outward=True, min_len=8, max_len=60)
    if len(lt) != 7:
        raise ValueError(f"fig4b: expected 7 left outward y ticks (70..10), got {len(lt)}: {lt}")
    ypairs = [(lt[i][0], 70.0 - 10.0 * i) for i in range(7)] + [(box["r1"], 0.0)]
    ya = H.AxisCal(ypairs, max_res_px=4.5)

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    r0_i, r1_i = box["r0_e"] + 10, box["r1_s"] - 4
    c0_i, c1_i = box["c0_e"] + 4, box["c1_s"] - 4
    grey = np.zeros(img.shape[:2], dtype=bool)
    grey[r0_i:r1_i, c0_i:c1_i] = (
        (mn[r0_i:r1_i, c0_i:c1_i] > 110) &
        (mx[r0_i:r1_i, c0_i:c1_i] < 205) &
        ((mx - mn)[r0_i:r1_i, c0_i:c1_i] < 20)
    )
    opened = H.open_rect(grey, 15, 15)
    dark = img.max(axis=2) < 100

    out = []
    for _, rsl, csl, _ in H.labeled_components(opened):
        c_mid = 0.5 * (csl.start + csl.stop)
        col_mid = int(round(c_mid))
        r_edge = rsl.start
        r_outer = r_edge
        for r in range(r_edge - 1, max(r0_i, r_edge - 25) - 1, -1):
            if dark[r, max(0, col_mid - 4):col_mid + 5].mean() > 0.5:
                r_outer = r
            else:
                break
        r_top = 0.5 * (r_edge + r_outer)
        x_lo = xa.to_value(csl.start - 6)
        x_hi = xa.to_value(csl.stop + 6)
        x_mid = xa.to_value(c_mid)
        y_val = ya.to_value(r_top)
        dh_px = 0.5 * (r_edge - r_outer) + 1.5 + ya.residual_px
        ya1, ya2 = ya.to_value(r_top - dh_px), ya.to_value(r_top + dh_px)
        out.append({
            "x": x_mid,
            "y": y_val,
            "lo": min(ya1, ya2),
            "hi": max(ya1, ya2),
            "px": [float(c_mid), float(r_top)],
            "px_lo": float(r_top - dh_px),
            "px_hi": float(r_top + dh_px),
            "note": f"histogram bin [{x_lo:.1f}, {x_hi:.1f}] um^2",
        })
    return sorted(out, key=lambda d: d["x"])


# ---------------------------------------------------------------------------
# Fig. 4(c): Droplet shape phase diagram — Area vs Approx. Filament Length
# ---------------------------------------------------------------------------
def read_fig4c(img: np.ndarray, shape: str) -> list[dict]:
    """Read Fig. 4(c) classified points: shape in ('tactoid', 'ellipsoid', 'sphere')."""
    box = H.find_box(img)
    bt = H.find_ticks(img, box, "bottom", min_len=10)
    if len(bt) != 8:
        raise ValueError(f"fig4c: expected 8 bottom x ticks, got {len(bt)}: {bt}")
    c_lat = bt[0][0]
    xa = H.AxisCal([(bt[3][0], 100.0), (bt[4][0], 150.0), (bt[5][0], 200.0), (bt[6][0], 250.0), (bt[7][0], 300.0)], max_res_px=3.5)
    rt = H.find_ticks(img, box, "right", min_len=10)
    if len(rt) != 9:
        raise ValueError(f"fig4c: expected 9 right y ticks (45..5), got {len(rt)}: {rt}")
    ya = H.AxisCal([(rt[i][0], 45.0 - 5.0 * i) for i in range(9)] + [(box["r1"], 0.0)], max_res_px=3.5)

    r_ch = img[:, :, 0].astype(np.int16)
    g_ch = img[:, :, 1].astype(np.int16)
    b_ch = img[:, :, 2].astype(np.int16)
    mx = np.maximum(np.maximum(r_ch, g_ch), b_ch)
    mn = np.minimum(np.minimum(r_ch, g_ch), b_ch)

    r0_i, r1_i = box["r0_e"] + 50, box["r1_s"] - 50
    c0_i, c1_i = box["c0_e"] + 50, box["c1_s"] - 50
    valid = np.zeros(img.shape[:2], dtype=bool)
    valid[r0_i:r1_i, c0_i:c1_i] = True
    leg = H.find_inner_frame_box(img, box)
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        valid[lr0:lr1 + 1, lc0:lc1 + 1] = False

    col_mask = np.zeros(img.shape[1], dtype=bool)
    col_bands = [
        (c_lat, 0.0),
        (xa.to_pixel(90.0), None),
        (xa.to_pixel(135.0), None),
        (xa.to_pixel(270.0), None),
    ]
    for cc, _ in col_bands:
        ci = int(round(cc))
        col_mask[max(0, ci - 55):min(img.shape[1], ci + 56)] = True
    valid &= col_mask[None, :]

    if shape == "tactoid":
        mask = valid & (r_ch > 180) & (r_ch - g_ch > 110) & (r_ch - b_ch > 110)
        m_h = 40.0
    elif shape == "ellipsoid":
        mask = valid & (mx < 75) & ((mx - mn) < 25)
        m_h = 38.0
    elif shape == "sphere":
        mask = valid & (b_ch > 165) & (g_ch > 110) & (r_ch < 70)
        m_h = 38.0
    else:
        raise ValueError(f"Unknown fig4c shape {shape!r}")

    opened = H.open_rect(mask, 8, 8)
    out = []
    for _, rsl, csl, loc in H.labeled_components(opened):
        h = rsl.stop - rsl.start
        w = csl.stop - csl.start
        if int(loc.sum()) < 220 or w < 18 or w > 95:
            continue
        if h <= 62:
            chunks = [(rsl.start, rsl.stop, loc)]
        else:
            k = max(2, int(round(h / 38.0)))
            edges = np.linspace(0, h, k + 1).astype(int)
            chunks = [(rsl.start + edges[j], rsl.start + edges[j + 1], loc[edges[j]:edges[j + 1], :]) for j in range(k)]
        for rs0, rs1, sub in chunks:
            if sub.sum() < 120:
                continue
            rr, cc = np.nonzero(sub)
            c_m = csl.start + float(cc.mean()) + 0.5
            r_m = rs0 + float(rr.mean()) + 0.5
            if abs(c_m - c_lat) < 55.0:
                fx = 0.0
                nt = f"majority {shape} (Lat = monomeric actin before axis break)"
            else:
                fx = None
                nt = f"majority {shape}"
            out.append(H.make_marker_reading(xa, ya, c_m, r_m, None, None, 0.5 * m_h, nt, force_x=fx, include_ebar=False))
    return sorted(out, key=lambda d: (round(d["x"], 0), d["y"]))
