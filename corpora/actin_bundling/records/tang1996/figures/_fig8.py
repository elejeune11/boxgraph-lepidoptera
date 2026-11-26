"""Figure 8a & 8b (p. 6 / 8561) readers for tang1996.

Fig. 8a: Scattering intensity vs [Calponin] (uM) for 2.5 uM F-actin at 30, 50, and 75 mM KCl.
Fig. 8b: Scattering intensity vs [Adenosine Phosphate] (mM) disassembling 2.5 uM F-actin +
         3 uM calponin bundles (50 mM KCl) by AMP (circles), ADP (diamonds), and ATP (downtris).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def _connected_components(dark_dil: np.ndarray, pts: list[tuple[float, float, int]], max_dist: float = 950.0):
    n = len(pts)
    adj: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dx = pts[j][0] - pts[i][0]
            dy = pts[j][1] - pts[i][1]
            dist = (dx * dx + dy * dy) ** 0.5
            if dx > 12 and dist < max_dist:
                if C.fast_segment_ink(dark_dil, pts[i][:2], pts[j][:2], margin_px=34.0) >= 0.85:
                    adj[i].append(j)
                    adj[j].append(i)
    vis = [False] * n
    comps = []
    for i in range(n):
        if not vis[i]:
            st = [i]
            vis[i] = True
            comp = []
            while st:
                u = st.pop()
                comp.append(u)
                for v in adj[u]:
                    if not vis[v]:
                        vis[v] = True
                        st.append(v)
            comps.append(sorted(comp, key=lambda idx: pts[idx][0]))
    return comps


def read_fig8a(image: np.ndarray, series: str) -> list[dict]:
    cal = P.calibrate("fig8a", image)
    ex = [cal.data_box_to_px(3.35, 0.042, 5.95, 0.165)]
    dark = C.mask_boxes(C._to_dark(image), ex)
    dark_dil = C.dilate_for_trace(dark, radius=4)
    all_pts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
    # Only build components on markers past x > 1.1 uM where all 3 curves are strictly disjoint!
    sep_pts = [p for p in all_pts if cal.px_to_data(p[0], p[1])[0] > 1.1]
    comps = [c for c in _connected_components(dark_dil, sep_pts) if len(c) >= 4]
    def sort_key(comp):
        hi = [sep_pts[i][0] for i in comp if cal.px_to_data(sep_pts[i][0], sep_pts[i][1])[1] > 0.2]
        return float(np.mean(hi)) if hi else sep_pts[comp[0]][0]
    comps.sort(key=sort_key)
    idx_map = {"kcl_30mM": 0, "kcl_50mM": 1, "kcl_75mM": 2}
    comp = comps[idx_map[series]]
    chain = [(sep_pts[i][0], sep_pts[i][1]) for i in comp]

    # Recover points at x <= 1.1 uM (0.92, 0.46, 0.0 uM) by tracing leftward from first marker
    x0_det, _ = cal.px_to_data(*chain[0])
    x1_det, _ = cal.px_to_data(*chain[1])
    dx_step = x1_det - x0_det
    c_cur, r_cur = chain[0]
    x_cur = x0_det
    prepended = []
    while x_cur - dx_step > -0.15:
        x_next = max(0.0, x_cur - dx_step)
        c_tgt = cal.col_for_x(x_next, r_cur)
        # All 3 curves lie on the low baseline (y <= 0.11) for x <= 1.0 uM!
        near = [p for p in all_pts if abs(cal.px_to_data(p[0], p[1])[0] - x_next) <= 0.15 and cal.px_to_data(p[0], p[1])[1] <= 0.11]
        if near:
            # 30 mM is top of baseline at 0.92 (~0.088), 50 mM middle (~0.072), 75 mM lowest (~0.068)
            near.sort(key=lambda p: p[1])  # ascending r (descending y)
            pick_i = min(idx_map[series], len(near) - 1)
            c_ref, r_ref = near[pick_i][0], near[pick_i][1]
        else:
            cls = C.column_clusters(dark, cal, x_next, y_band=(0.050, 0.105), half_c=4, max_gap=18)
            r_tr = cls[0][0] if cls else cal.row_for_y(0.06, c_tgt)
            c_ref, r_ref, _ = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=24)
        prepended.append((c_ref, r_ref))
        c_cur, r_cur, x_cur = c_ref, r_ref, x_next
        if x_next <= 0.05:
            break
    chain = list(reversed(prepended)) + chain
    h_px = max(3.5, 2.0 + cal.resid_y_px)
    return [C.make_reading(cal, c, r, h_px, f"Fig. 8a {series} (cal resid {cal.resid_y_px:.2f} px)") for c, r in chain]


def read_fig8b(image: np.ndarray, series: str) -> list[dict]:
    cal = P.calibrate("fig8b", image)
    ex = [
        cal.data_box_to_px(1.1, 1.12, 22.0, 1.85),   # title text
        cal.data_box_to_px(5.8, 0.34, 22.0, 0.95),   # legend box (+AMP/+ADP/+ATP)
    ]
    dark = C.mask_boxes(C._to_dark(image), ex)
    dark_dil = C.dilate_for_trace(dark, radius=4)
    pts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
    comps = [c for c in _connected_components(dark_dil, pts) if len(c) >= 5]
    # Order top-to-bottom by average y at x in [2, 6] mM: AMP (top), ADP (middle), ATP (bottom)
    def mid_y(comp):
        ys = [cal.px_to_data(pts[i][0], pts[i][1])[1] for i in comp if 1.8 <= cal.px_to_data(pts[i][0], pts[i][1])[0] <= 6.0]
        return float(np.mean(ys)) if ys else 0.0
    comps.sort(key=mid_y, reverse=True)
    idx_map = {"amp": 0, "adp": 1, "atp": 2}
    comp = comps[idx_map[series]]
    chain = [(pts[i][0], pts[i][1]) for i in comp]

    # Recover leftmost point at x = 0.50 mM by tracing leftward from first detected marker
    c_first, r_first = chain[0]
    x_first, _ = cal.px_to_data(c_first, r_first)
    if x_first > 0.7:
        c_tgt = cal.col_for_x(0.50, r_first)
        _, r_tr = C.trace_step_along_curve(dark, cal, c_first, r_first, c_tgt, search_r_px=60, step_c=16.0)
        c_ref, r_ref, _ = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=30)
        chain.insert(0, (c_ref, r_ref))

    # For ADP, locate rightmost diamond (~15.5 mM) along its flat plateau (y ~ 0.14) without drifting up to +AMP circle (~0.16)
    if series == "adp":
        c_last, r_last = chain[-1]
        x_last, _ = cal.px_to_data(c_last, r_last)
        if x_last < 13.0:
            c_tgt = cal.col_for_x(15.5, r_last)
            pos_d, neg_d = C.make_geom_template("diamond", radius_px=25.0, stroke_px=6.5)
            _, bc, br = C.find_template_in_band(dark, (int(c_tgt) - 45, int(c_tgt) + 45), (int(r_last) - 24, int(r_last) + 24), pos_d, neg_d, neg_w=0.45)
            chain.append((bc, br))

    h_px = max(3.5, 2.0 + cal.resid_y_px)
    return [C.make_reading(cal, c, r, h_px, f"Fig. 8b +{series.upper()} (cal resid {cal.resid_y_px:.2f} px)") for c, r in chain]
