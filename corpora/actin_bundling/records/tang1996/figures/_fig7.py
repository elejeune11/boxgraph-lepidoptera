"""Figure 7 (p. 5 / 8560) reader for tang1996.

Scattering intensity vs peptide concentration (mM, log-log) for Lys42 (circles),
MARCKS peptide (diamonds), and Lys18 (squares).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def _connected_components(dark_dil: np.ndarray, pts: list[tuple[float, float, int]], max_dist: float = 900.0):
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


def read_fig7(image: np.ndarray, series: str) -> list[dict]:
    cal = P.calibrate("fig7", image)
    ex = [
        cal.data_box_to_px(0.00013, 0.08, 0.00038, 0.145),  # "Lys42"
        cal.data_box_to_px(0.00018, 0.032, 0.0036, 0.070),  # "MARCKS peptide"
        cal.data_box_to_px(0.0085, 0.022, 0.030, 0.068),    # "Lys18"
    ]
    dark = C.mask_boxes(C._to_dark(image), ex)
    dark_dil = C.dilate_for_trace(dark, radius=4)
    pts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
    comps = [c for c in _connected_components(dark_dil, pts) if len(c) >= 4]
    # Order the 3 main curves left-to-right by their starting x!
    comps.sort(key=lambda comp: pts[comp[0]][0])
    h_px = max(3.5, 2.0 + cal.resid_y_px)

    if series == "lys18":
        # Rightmost curve (8 clean open squares)
        comp = comps[2]
        return [C.make_reading(cal, pts[i][0], pts[i][1], h_px, f"Fig. 7 Lys18 square (cal resid {cal.resid_y_px:.2f} px)") for i in comp]

    if series == "lys42":
        # Leftmost curve (open circles): also recover the leftmost point touching the left axis border
        # (x = 0.0001 mM) by tracing leftward along the connecting line from pts[comp[0]], and check
        # the dense plateau markers around 0.008 and 0.010 mM where Lys42 and MARCKS converge.
        comp = comps[0]
        chain = [(pts[i][0], pts[i][1]) for i in comp]
        # Leftmost point on left frame border (x = 0.0001 mM): trace leftward from first detected circle
        c_first, r_first = chain[0]
        c_left_edge = cal.col_for_x(0.0001, r_first)
        _, r_left = C.trace_step_along_curve(dark, cal, c_first, r_first, c_left_edge + 12.0, search_r_px=65, step_c=14.0)
        dx_l = (c_left_edge + 12.0) - c_first
        slope_l = (r_left - r_first) / dx_l if abs(dx_l) > 1.0 else 0.0
        r_border = r_left + slope_l * (-12.0)
        chain.insert(0, (c_left_edge, r_border))
        # Lys42 plateau runs horizontally from 0.005 mM (~1.41) directly to the circle at 0.010 mM (~1.44),
        # passing below the MARCKS diamonds at 0.008 mM (~1.50) and 0.010 mM (~1.62).
        chain = [p for p in chain if cal.px_to_data(p[0], p[1])[0] < 0.007]
        pos_c, neg_c = C.make_geom_template("circle", radius_px=24.0, stroke_px=6.5)
        c_e = cal.col_for_x(0.010, chain[-1][1])
        sc, bc, br = C.find_template_in_band(
            dark,
            (int(c_e) - 45, int(c_e) + 45),
            (int(chain[-1][1]) - 30, int(chain[-1][1]) + 35),
            pos_c, neg_c, neg_w=0.45,
        )
        if sc > 0.35:
            chain.append((bc, br))
        return [C.make_reading(cal, c, r, h_px, f"Fig. 7 Lys42 circle (cal resid {cal.resid_y_px:.2f} px)") for c, r in chain]

    if series == "marcks":
        # Middle curve (open diamonds): extends from comp[1] (0.0007..0.006 mM) to diamonds at ~0.008, 0.010, 0.012 mM
        comp = comps[1]
        chain = [(pts[i][0], pts[i][1]) for i in comp]
        pos_d, neg_d = C.make_geom_template("diamond", radius_px=24.0, stroke_px=6.5)
        # Rightmost detected point at ~0.012 mM in pts is a clean diamond; check 0.008, 0.010, 0.012
        for xv in (0.008, 0.010, 0.012):
            c_e = cal.col_for_x(xv, chain[-1][1])
            # Prefer isolated detected hollow marker at 0.012 if present
            near = [p for p in pts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 0.0008 and cal.px_to_data(p[0], p[1])[1] >= 1.2]
            if near and xv >= 0.011:
                chain.append((near[0][0], near[0][1]))
            else:
                sc, bc, br = C.find_template_in_band(dark, (int(c_e) - 45, int(c_e) + 45), (int(chain[-1][1]) - 85, int(chain[-1][1]) + 35), pos_d, neg_d, neg_w=0.30)
                if sc > 0.45:
                    chain.append((bc, br))
        return [C.make_reading(cal, c, r, h_px, f"Fig. 7 MARCKS peptide diamond (cal resid {cal.resid_y_px:.2f} px)") for c, r in chain]

    raise ValueError(f"Unknown Fig. 7 series: {series}")
