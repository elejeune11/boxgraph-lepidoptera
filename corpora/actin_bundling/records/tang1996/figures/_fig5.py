"""Figure 5a & 5b (p. 5 / 8560) readers for tang1996.

Fig. 5a: Scattering intensity vs [Co(NH3)6^3+] (mM) at 30, 50, 100, and 150 mM KCl.
Fig. 5b: Bundling onset [Mn2+] (4 filled circles + fit line) and [Co(NH3)6^3+]
         (12 hollow symbols across 3 independent preparations: circles, up-triangles,
         down-triangles + fit line) vs [KCl] (mM).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def _fig5b_ex(cal: C.PanelCal) -> list[tuple[int, int, int, int]]:
    return [
        cal.data_box_to_px(35.0, 3.3, 125.0, 5.7),   # [Mn2+]=3.11 x [KCl]^0.20
        cal.data_box_to_px(35.0, 0.50, 155.0, 0.90), # [Co(NH3)6^3+]=0.068 x [KCl]^0.74
    ]


def read_fig5b_mn_points(image: np.ndarray) -> list[dict]:
    cal = P.calibrate("fig5b", image)
    pts = [p for p in C.detect_filled_markers(image, cal, exclude_boxes_px=_fig5b_ex(cal)) if cal.px_to_data(p[0], p[1])[1] >= 5.0]
    h_px = max(3.5, 2.0 + cal.resid_y_px)
    return [C.make_reading(cal, c, r, h_px, f"Fig. 5b Mn2+ filled circle (cal resid {cal.resid_y_px:.2f} px)") for c, r, _ in pts]


def read_fig5b_fit_line(image: np.ndarray, at: list[float], species: str) -> list[dict]:
    cal = P.calibrate("fig5b", image)
    ex = _fig5b_ex(cal)
    if species == "mn":
        return C.read_line_at(image, cal, at, y_band=(5.5, 9.5), exclude_boxes_px=ex, note_prefix="Fig. 5b Mn2+ power-law fit line")
    if species == "cohex":
        return C.read_line_at(image, cal, at, y_band=(0.8, 3.0), exclude_boxes_px=ex, note_prefix="Fig. 5b Co(NH3)6^3+ power-law fit line")
    raise ValueError(f"Unknown Fig. 5b fit species: {species}")


def read_fig5b_cohex_prep(image: np.ndarray, shape: str) -> list[dict]:
    """Read one of the 3 independent Co(NH3)6^3+ preparations (circle, uptri, downtri)
    at the 4 automatically detected KCl clusters in the lower half of Fig. 5b."""
    cal = P.calibrate("fig5b", image)
    ex = _fig5b_ex(cal)
    co_blobs = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex, max_bbox_2x=95, max_hole_area_2x=1500, fuse_px_2x=12)
    co_blobs = [p for p in co_blobs if cal.px_to_data(p[0], p[1])[1] < 4.5]
    clusters: list[list] = []
    for p in sorted(co_blobs, key=lambda t: t[0]):
        if not clusters or p[0] - clusters[-1][-1][0] > 60:
            clusters.append([p])
        else:
            clusters[-1].append(p)
    dark = C.mask_boxes(C._to_dark(image), ex)
    pos, neg = C.make_geom_template(shape, radius_px=29.0, stroke_px=7.5)
    neg_w = 0.35 if shape == "circle" else 0.20
    out = []
    for cl in clusters:
        c_mid = float(np.mean([p[0] for p in cl]))
        r_mid = float(np.mean([p[1] for p in cl]))
        sc, bc, br = C.find_template_in_band(dark, (int(c_mid) - 50, int(c_mid) + 50), (int(r_mid) - 120, int(r_mid) + 120), pos, neg, neg_w=neg_w)
        h_px = max(4.5, 2.5 + cal.resid_y_px)
        out.append(C.make_reading(cal, bc, br, h_px, f"Fig. 5b Co(NH3)6^3+ {shape} (score {sc:.2f}, cal resid {cal.resid_y_px:.2f} px)"))
    return out


def read_fig5a(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read one of the 4 KCl curves in Fig. 5a (30, 50, 100, 150 mM KCl) at experimental x positions `at`."""
    cal = P.calibrate("fig5a", image)
    ex = [cal.data_box_to_px(3.15, 0.075, 5.4, 0.45)]
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    cross = C.row_clusters(dark, cal, 0.25, x_band=(0.4, 3.2), half_r=4, max_gap=6)
    rank_map = {"kcl_30mM": 0, "kcl_50mM": 1, "kcl_100mM": 2, "kcl_150mM": 3}
    rk = rank_map[series]
    c_cross = cross[rk][0]
    r_cross = cal.row_for_y(0.25, c_cross)
    x_cross, _ = cal.px_to_data(c_cross, r_cross)

    shape_map = {"kcl_30mM": "circle", "kcl_50mM": "square", "kcl_100mM": "uptri", "kcl_150mM": "diamond"}
    pos_t, neg_t = C.make_geom_template(shape_map[series], radius_px=23.0, stroke_px=6.5)

    def follow_rise(r_start, c_start, step_r, y_lim):
        c, r = c_start, r_start
        while True:
            r_next = r + step_r
            _, yv = cal.px_to_data(c, r_next)
            if (step_r < 0 and yv >= y_lim) or (step_r > 0 and yv <= y_lim):
                break
            if r_next < cal.inner_box[1] + 10 or r_next > cal.inner_box[3] - 10:
                break
            ri = int(round(r_next))
            c_lo, c_hi = max(cal.inner_box[0] + 4, int(round(c)) - 45), min(cal.inner_box[2] - 4, int(round(c)) + 46)
            hits = np.nonzero(dark[ri - 2:ri + 3, c_lo:c_hi].any(axis=0))[0]
            if len(hits) == 0:
                break
            groups = [c_lo + 0.5 * (g[0] + g[-1]) for g in np.split(hits, np.where(np.diff(hits) > 4)[0] + 1)]
            c_best = min(groups, key=lambda cm: abs(cm - c))
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    if 16 <= groups[j] - groups[i] <= 55 and abs(0.5 * (groups[i] + groups[j]) - c) <= 22:
                        c_best = 0.5 * (groups[i] + groups[j])
            c, r = c_best, r_next
        return c, r

    c_top, r_top = follow_rise(r_cross, c_cross, -18.0, 0.58)
    c_bot, r_bot = follow_rise(r_cross, c_cross, +18.0, 0.085)
    x_bot, _ = cal.px_to_data(c_bot, r_bot)
    x_top, _ = cal.px_to_data(c_top, r_top)

    xs = sorted(float(x) for x in at)
    hi_map: dict[int, tuple[float, float, float]] = {}
    c_cur_hi, r_cur_hi = c_top, r_top
    for i, xv in enumerate(xs):
        if xv < x_top - 0.05:
            continue
        c_tgt = cal.col_for_x(xv, r_cur_hi)
        if abs(xv - x_top) <= 0.25:
            near = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 0.12 and cal.px_to_data(p[0], p[1])[1] >= 0.45]
            if near:
                best = min(near, key=lambda p: abs(p[1] - r_cur_hi))
                c_ref, r_ref, err = best[0], best[1], max(3.5, 2.0 + cal.resid_y_px)
            else:
                _, r_tr = C.trace_step_along_curve(dark, cal, c_cur_hi, r_cur_hi, c_tgt, search_r_px=65, step_c=14.0)
                c_ref, r_ref, err = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=30)
        else:
            _, r_tr = C.trace_step_along_curve(dark, cal, c_cur_hi, r_cur_hi, c_tgt, search_r_px=60, step_c=14.0)
            near = [p for p in hpts if abs(p[0] - c_tgt) <= 32 and abs(p[1] - r_tr) <= 45]
            if near:
                best = min(near, key=lambda p: abs(p[1] - r_tr))
                c_ref, r_ref, err = best[0], best[1], max(3.5, 2.0 + cal.resid_y_px)
            else:
                sc, c_ref, r_ref = C.find_template_in_band(dark, (int(c_tgt) - 32, int(c_tgt) + 32), (int(r_tr) - 55, int(r_tr) + 55), pos_t, neg_t, neg_w=0.25)
                if sc < 0.40:
                    c_ref, r_ref, err = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=30)
                else:
                    err = max(4.0, 2.5 + cal.resid_y_px)
        hi_map[i] = (c_ref, r_ref, err)
        c_cur_hi, r_cur_hi = c_ref, r_ref

    mid_map: dict[int, tuple[float, float, float]] = {}
    for i, xv in enumerate(xs):
        if i in hi_map or xv < x_bot - 0.08:
            continue
        # Search for detected hollow marker on the rising segment (y in [0.06, 0.65])
        near = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 0.06 and 0.065 <= cal.px_to_data(p[0], p[1])[1] <= 0.65]
        frac = max(0.0, min(1.0, (xv - x_bot) / max(1e-6, (x_top - x_bot))))
        r_est = r_bot + frac * (r_top - r_bot)
        if near:
            best = min(near, key=lambda p: abs(p[1] - r_est))
            mid_map[i] = (best[0], best[1], max(3.5, 2.0 + cal.resid_y_px))
        else:
            c_est = cal.col_for_x(xv, r_est)
            sc, bc, br = C.find_template_in_band(dark, (int(c_est) - 35, int(c_est) + 35), (int(r_est) - 160, int(r_est) + 160), pos_t, neg_t, neg_w=0.25)
            if sc > 0.35:
                mid_map[i] = (bc, br, max(4.0, 2.5 + cal.resid_y_px))
            else:
                c_ref, r_ref, err = C.refine_hollow_or_ink_at(dark, cal, c_est, r_est, win=34)
                mid_map[i] = (c_ref, r_ref, err)

    lo_map: dict[int, tuple[float, float, float]] = {}
    # Seed baseline right at the last baseline marker (x_bot)
    c_cur_lo, r_cur_lo = c_bot, r_bot
    for i in range(len(xs) - 1, -1, -1):
        if i in hi_map or i in mid_map:
            continue
        xv = xs[i]
        c_tgt = cal.col_for_x(xv, r_cur_lo)
        near = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 0.10 and cal.px_to_data(p[0], p[1])[1] <= 0.095]
        if len(near) == 1:
            c_ref, r_ref, err = near[0][0], near[0][1], max(3.5, 2.0 + cal.resid_y_px)
        else:
            cls = C.column_clusters(dark, cal, xv, y_band=(0.055, 0.095), half_c=3, max_gap=18)
            r_tr = cls[0][0] if cls else r_cur_lo
            if near:
                best = min(near, key=lambda p: abs(p[1] - r_tr))
                c_ref, r_ref, err = best[0], best[1], max(3.5, 2.0 + cal.resid_y_px)
            else:
                c_ref, r_ref, err = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=24)
        lo_map[i] = (c_ref, r_ref, err)
        c_cur_lo, r_cur_lo = c_ref, r_ref

    out = []
    for i in range(len(xs)):
        c_f, r_f, err = hi_map.get(i) or mid_map.get(i) or lo_map[i]
        out.append(C.make_reading(cal, c_f, r_f, err, f"Fig. 5a {series} (cal resid {cal.resid_y_px:.2f} px)"))
    return out
