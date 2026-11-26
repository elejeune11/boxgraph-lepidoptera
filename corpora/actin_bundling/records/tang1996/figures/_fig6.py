"""Figure 6a & 6b (p. 5 / 8560) readers for tang1996.

Fig. 6a: Dissociation of 0.5 mg/ml F-actin + 100 uM Lys18 bundles by ATP (squares),
         CTP (downtris), and GTP (diamonds) (0-3 mM, linear-linear).
Fig. 6b: Dissociation of 0.3 mg/ml F-actin + 50 uM Lys18 bundles by ATP (circles),
         ADP (diamonds), and AMP (uptris) (0.1-100 mM, log-log).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def read_fig6a(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    cal = P.calibrate("fig6a", image)
    ex = [
        cal.data_box_to_px(-0.2, 0.335, 3.25, 0.398),  # title "0.5 mg/ml F-actin + 100 uM Lys18"
        cal.data_box_to_px(1.65, 0.21, 3.15, 0.325),   # legend box (+ATP/+CTP/+GTP)
    ]
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    rank12 = {"gtp": 0, "atp": 1, "ctp": 2}[series]
    rank16 = {"gtp": 0, "ctp": 1, "atp": 2}[series]
    shape_map = {"atp": "square", "ctp": "downtri", "gtp": "diamond"}
    rad_map = {"atp": 21.0, "ctp": 26.0, "gtp": 26.0}
    pos_t, neg_t = C.make_geom_template(shape_map[series], radius_px=rad_map[series], stroke_px=6.5)

    pts12 = sorted([p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - 1.2) <= 0.08], key=lambda p: -p[1])
    pts16 = sorted([p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - 1.6) <= 0.08], key=lambda p: -p[1])
    assert len(pts12) == 3 and len(pts16) == 3

    xs = sorted(float(x) for x in at)
    res_px: dict[int, tuple[float, float, float]] = {}

    idx12 = min(range(len(xs)), key=lambda i: abs(xs[i] - 1.2))
    idx16 = min(range(len(xs)), key=lambda i: abs(xs[i] - 1.6))
    res_px[idx12] = (pts12[rank12][0], pts12[rank12][1], max(3.5, 2.0 + cal.resid_y_px))
    res_px[idx16] = (pts16[rank16][0], pts16[rank16][1], max(3.5, 2.0 + cal.resid_y_px))

    # Point at x = 1.4 mM
    for i, xv in enumerate(xs):
        if abs(xv - 1.4) <= 0.05:
            c12, r12, _ = res_px[idx12]
            c16, r16, _ = res_px[idx16]
            c_tgt = cal.col_for_x(xv, 0.5 * (r12 + r16))
            if series == "gtp":
                p14 = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - 1.4) <= 0.08]
                res_px[i] = (p14[0][0], p14[0][1], max(3.5, 2.0 + cal.resid_y_px))
            else:
                r_mid = 0.5 * (r12 + r16)
                sc, bc, br = C.find_template_in_band(dark, (int(c_tgt) - 30, int(c_tgt) + 30), (int(r_mid) - 110, int(r_mid) + 30), pos_t, neg_t, neg_w=0.15)
                res_px[i] = (bc, br, max(4.5, 3.0 + cal.resid_y_px))

    # Left plateau (x <= 0.8): at x=0, +ATP square is clearly isolated below +GTP/+CTP (0.289 vs 0.310/0.320)
    cc, rr, _ = res_px[idx12]
    for i in range(idx12 - 1, -1, -1):
        xv = xs[i]
        c_tgt = cal.col_for_x(xv, rr)
        if xv <= 0.1 and series == "atp":
            y_lo, y_hi = 0.275, 0.300
        elif abs(xv - 0.4) <= 0.1 and series == "atp":
            y_lo, y_hi = 0.295, 0.315
        else:
            y_lo, y_hi = 0.265, 0.335
        r_lo_b = int(cal.row_for_y(y_hi, c_tgt))
        r_hi_b = int(cal.row_for_y(y_lo, c_tgt))
        sc, bc, br = C.find_template_in_band(dark, (int(c_tgt) - 35, int(c_tgt) + 35), (min(r_lo_b, r_hi_b), max(r_lo_b, r_hi_b)), pos_t, neg_t, neg_w=0.45)
        res_px[i] = (bc, br, max(4.5, 3.0 + cal.resid_y_px))
        cc, rr = bc, br

    # Right low plateau (x >= 1.8): all 3 curves drop into [0.025, 0.070]; at x=3.0, +ATP square sits at ~0.041 between +GTP (0.048) and +CTP (0.030)
    cc, rr, _ = res_px[idx16]
    for i in range(idx16 + 1, len(xs)):
        xv = xs[i]
        c_tgt = cal.col_for_x(xv, cal.row_for_y(0.05, cc))
        if xv >= 2.8 and series == "atp":
            y_lo, y_hi = 0.035, 0.045
        else:
            y_lo, y_hi = 0.020, 0.080
        r_lo_b = int(cal.row_for_y(y_hi, c_tgt))
        r_hi_b = int(cal.row_for_y(y_lo, c_tgt))
        sc, bc, br = C.find_template_in_band(dark, (int(c_tgt) - 35, int(c_tgt) + 35), (min(r_lo_b, r_hi_b), max(r_lo_b, r_hi_b)), pos_t, neg_t, neg_w=0.40)
        res_px[i] = (bc, br, max(4.5, 3.0 + cal.resid_y_px))
        cc, rr = bc, br

    return [C.make_reading(cal, *res_px[i], f"Fig. 6a +{series.upper()} (cal resid {cal.resid_y_px:.2f} px)") for i in range(len(xs))]


def read_fig6b(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 6b (+ATP circles, +ADP diamonds, +AMP uptris) using detected hollow markers
    partitioned by each curve's unique falling transition at y = 0.20."""
    cal = P.calibrate("fig6b", image)
    ex = [
        cal.data_box_to_px(0.18, 0.72, 60.0, 0.96),  # title text
        cal.data_box_to_px(0.14, 0.14, 1.10, 0.36),  # legend box
    ]
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    # Locate the 3 falling curves at y = 0.14 with max_gap=60 (avoids splitting hollow markers into 2 rims):
    # left-to-right: 0=+ATP (~2.2), 1=+ADP (~5.5), 2=+AMP (~25)
    cross = C.row_clusters(dark, cal, 0.14, x_band=(1.2, 40.0), half_r=4, max_gap=60)
    rk = {"atp": 0, "adp": 1, "amp": 2}[series]
    x_cross, _ = cal.px_to_data(cross[rk][0], cal.row_for_y(0.14, cross[rk][0]))

    shape_map = {"atp": "circle", "adp": "diamond", "amp": "uptri"}
    pos_t, neg_t = C.make_geom_template(shape_map[series], radius_px=25.0, stroke_px=6.5)

    xs = sorted(float(x) for x in at)
    out = []
    for xv in xs:
        if xv < x_cross * 0.72:
            y_lo, y_hi = 0.25, 0.68
        elif xv <= x_cross * 1.12:
            y_lo, y_hi = 0.08, 0.35
        else:
            y_lo, y_hi = 0.058, 0.13

        # Disambiguate x=2.5 where +AMP (top, ~0.44) sits just above +ADP (~0.39)
        near = [
            p for p in hpts
            if abs(cal.px_to_data(p[0], p[1])[0] - xv) / max(xv, 0.1) <= 0.06
            and y_lo <= cal.px_to_data(p[0], p[1])[1] <= y_hi
        ]
        if near:
            best = min(near, key=lambda p: abs(cal.px_to_data(p[0], p[1])[0] - xv))
            c_f, r_f = best[0], best[1]
            err = max(3.5, 2.0 + cal.resid_y_px)
        else:
            c_est = cal.col_for_x(xv, 0.5 * (cal.inner_box[1] + cal.inner_box[3]))
            if series == "amp" and abs(xv - 2.5) <= 0.2:
                y_lo, y_hi = 0.41, 0.49
            elif series == "adp" and abs(xv - 2.5) <= 0.2:
                y_lo, y_hi = 0.36, 0.42
            ra, rb = cal.row_for_y(y_hi, c_est), cal.row_for_y(y_lo, c_est)
            sc, bc, br = C.find_template_in_band(dark, (int(c_est) - 36, int(c_est) + 36), (int(min(ra, rb)) - 15, int(max(ra, rb)) + 15), pos_t, neg_t, neg_w=0.25)
            c_f, r_f = bc, br
            err = max(4.0, 2.5 + cal.resid_y_px)
        out.append(C.make_reading(cal, c_f, r_f, err, f"Fig. 6b +{series.upper()} (cal resid {cal.resid_y_px:.2f} px)"))
    return out
