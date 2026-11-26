"""Figure 4 (p. 4 / 8559) reader for tang1996.

Scattering intensity vs [Lys18] (uM) for 0.2 mg/ml F-actin at six gelsolin:actin ratios:
no gel (uptri), 1:500 (circle), 1:100 (diamond), 1:50 (downtri), 1:30 (square), 1:10 (dotcircle).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def read_fig4(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    cal = P.calibrate("fig4", image)
    ex = [
        cal.data_box_to_px(-15.0, 0.56, 31.0, 1.18),  # main legend block
        cal.data_box_to_px(31.0, 0.89, 57.0, 1.07),   # "(gel:actin)" overhang
    ]
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    # At x = 44.4 uM, all 6 curves are cleanly separated vertically from bottom to top:
    # 0=1:10, 1=1:30, 2=1:50, 3=1:100, 4=1:500, 5=no_gel
    rank_map = {
        "gel_1_10": 0,
        "gel_1_30": 1,
        "gel_1_50": 2,
        "gel_1_100": 3,
        "gel_1_500": 4,
        "no_gel": 5,
    }
    rk = rank_map[series]
    pts44 = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - 44.4) <= 2.5]
    pts44.sort(key=lambda p: -p[1])  # ascending y (bottom to top)
    assert len(pts44) == 6, f"Expected 6 separated markers at x=44.4 uM, got {len(pts44)}"
    c_seed, r_seed, _ = pts44[rk]

    # Templates to distinguish no_gel (uptri) vs 1:500 (circle) where they meet/cross at x=55.5..83.3
    shape_map = {
        "gel_1_10": "dotcircle",
        "gel_1_30": "square",
        "gel_1_50": "downtri",
        "gel_1_100": "diamond",
        "gel_1_500": "circle",
        "no_gel": "uptri",
    }
    pos_t, neg_t = C.make_geom_template(shape_map[series], radius_px=26.0, stroke_px=6.5)

    xs = sorted(float(x) for x in at)
    idx_s = min(range(len(xs)), key=lambda i: abs(xs[i] - 44.4))
    res_px: dict[int, tuple[float, float, float]] = {}
    res_px[idx_s] = (c_seed, r_seed, max(3.5, 2.0 + cal.resid_y_px))

    def snap_next(cc: float, rr: float, x_next: float) -> tuple[float, float, float]:
        c_tgt = cal.col_for_x(x_next, rr)
        # If tracing leftward from x=33.3 -> 22.2 on the steep upper curves (rk >= 3),
        # give a downward slope hint (smaller x -> lower y -> larger r) and wider vertical search window!
        x_cur, _ = cal.px_to_data(cc, rr)
        if x_next < x_cur and x_next <= 24.0 and rk >= 3:
            # Find ink strictly below y=0.20 at x=22.2 (all 6 curves lie below y=0.18 at x<=22.2!)
            cls = C.column_clusters(dark, cal, x_next, y_band=(0.03, 0.19), half_c=4, max_gap=18)
            r_tr = cls[0][0] if (cls and rk >= 4) else (cls[-1][0] if cls else cal.row_for_y(0.08, c_tgt))
            if x_next <= 13.0 and cls:
                r_tr = float(np.median([c_[0] for c_ in cls]))
            return C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=26)
        # At x = 111.1 uM where 5 markers are cleanly detected (no_gel top, 1:500 2nd, 1:100/1:50 3rd, 1:30 4th, 1:10 5th)
        if x_next >= 108.0:
            pts111 = sorted([p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - 111.2) <= 3.0], key=lambda p: -p[1])
            if len(pts111) == 5:
                m_idx = {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4}[rk]
                return pts111[m_idx][0], pts111[m_idx][1], max(4.0, 2.5 + cal.resid_y_px)
        _, r_tr = C.trace_step_along_curve(dark, cal, cc, rr, c_tgt, search_r_px=95, step_c=12.0)
        # Prefer isolated detected hollow marker near traced position
        near = [p for p in hpts if abs(p[0] - c_tgt) <= 35 and abs(p[1] - r_tr) <= 45]
        if len(near) == 1:
            return near[0][0], near[0][1], max(3.5, 2.0 + cal.resid_y_px)
        # At x = 66.7 uM for top two curves (no_gel uptri slightly above 1:500 circle)
        if abs(x_next - 66.7) <= 2.5 and rk in (4, 5):
            sc, bc, br = C.find_template_in_band(dark, (int(c_tgt) - 35, int(c_tgt) + 35), (int(r_tr) - 70, int(r_tr) + 70), pos_t, neg_t, neg_w=0.30)
            if sc > 0.45:
                return bc, br, max(4.5, 3.0 + cal.resid_y_px)
        c_ref, r_ref, err = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=30)
        return c_ref, r_ref, err

    # Trace rightward from x = 44.4
    cc, rr = c_seed, r_seed
    for i in range(idx_s + 1, len(xs)):
        res_px[i] = snap_next(cc, rr, xs[i])
        cc, rr = res_px[i][0], res_px[i][1]

    # Trace leftward from x = 44.4
    cc, rr = c_seed, r_seed
    for i in range(idx_s - 1, -1, -1):
        res_px[i] = snap_next(cc, rr, xs[i])
        cc, rr = res_px[i][0], res_px[i][1]

    return [C.make_reading(cal, *res_px[i], f"Fig. 4 {series} (cal resid {cal.resid_y_px:.2f} px)") for i in range(len(xs))]
