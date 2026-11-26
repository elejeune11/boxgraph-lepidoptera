"""Figure 1 (p. 2 / 8557) reader for tang1996.

Scattering intensity vs concentration (0.1-100 mM, log-log) for 0.5 mg/ml F-actin
induced by 8 cations: Arg6 (filled circle), His6 (filled uptri), His5 (open downtri),
[Co(NH3)6]3+ (open diamond), Co2+ (filled downtri), Mn2+ (open uptri),
Ca2+ (open square), and Mg2+ (open circle).
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def _fig1_ex(cal: C.PanelCal) -> list[tuple[int, int, int, int]]:
    return [
        cal.data_box_to_px(0.22, 1.35, 0.85, 2.30),   # "Arg6"
        cal.data_box_to_px(0.85, 1.12, 1.75, 1.95),   # "His6"
        cal.data_box_to_px(2.00, 1.38, 18.0, 2.40),   # "[Co(NH3)6]3+"
        cal.data_box_to_px(1.55, 0.16, 3.00, 0.35),   # "His5"
        cal.data_box_to_px(2.00, 0.03, 95.0, 0.056),  # "Co2+ Mn2+ Ca2+ Mg2+" labels
    ]


def read_fig1(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read one of the 8 cation curves in Fig. 1 at its experimental marker concentrations `at`."""
    cal = P.calibrate("fig1", image)
    ex = _fig1_ex(cal)
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
    fpts = C.detect_filled_markers(image, cal, exclude_boxes_px=ex, erode_2x=8, min_edt_px=8.5, min_core_size_2x=14)
    # Remove the 4 arrow tips (all lie strictly inside y in [0.11, 0.28] for x > 4.0)
    fpts = [p for p in fpts if not (cal.px_to_data(p[0], p[1])[0] > 4.0 and 0.11 <= cal.px_to_data(p[0], p[1])[1] <= 0.28)]

    # Locate the 8 rising edges unambiguously left-to-right at y = 0.35
    cross = C.row_clusters(dark, cal, 0.35, x_band=(0.15, 60.0), half_r=3, max_gap=6)
    rank_map = {
        "arg6": 0, "his6": 1, "his5": 2, "cohex": 3,
        "co2": 4, "mn2": 5, "ca2": 6, "mg2": 7,
    }
    rk = rank_map[series]
    c_cross = cross[rk][0]
    x_cross, _ = cal.px_to_data(c_cross, cal.row_for_y(0.35, c_cross))

    is_filled = series in ("arg6", "his6", "co2")
    shape_map = {
        "arg6": "filled_circle", "his6": "filled_uptri", "co2": "filled_downtri",
        "his5": "downtri", "cohex": "diamond", "mn2": "uptri",
        "ca2": "square", "mg2": "circle",
    }
    rad = 22.0 if is_filled else 24.0
    pos_t, neg_t = C.make_geom_template(shape_map[series], radius_px=rad, stroke_px=6.5)

    xs = sorted(float(x) for x in at)
    pts_pool = fpts if is_filled else hpts
    out = []

    for xv in xs:
        c_est = cal.col_for_x(xv, 0.5 * (cal.inner_box[1] + cal.inner_box[3]))
        if xv < x_cross * 0.96 or (series == "cohex" and xv <= 3.05):
            if series == "ca2" and xv <= 6.0:
                # Open square at x=5 mM sits directly below the solid Co2+ downtri (~0.082)
                y_lo, y_hi = 0.060, 0.076
            elif series == "mg2" and xv <= 12.0:
                # Open circle at x=10 mM sits slightly above the open Ca2+ square (~0.066)
                y_lo, y_hi = 0.070, 0.082
            elif series == "his5" and xv <= 0.7:
                # Open downtri at x=0.6 mM overlaps solid His6 uptri (~0.060)
                y_lo, y_hi = 0.056, 0.070
            else:
                y_lo, y_hi = 0.055, 0.115
        else:
            if series == "his5":
                y_lo, y_hi = 0.30, 0.65
            elif series == "ca2" and xv >= 24.0:
                # Ca2+ upper plateau squares sit strictly above Mg2+ circles (1.15..1.45 vs 0.80..1.12)
                y_lo, y_hi = 1.15, 1.45
            elif series == "mg2":
                # Mg2+ upper plateau circles sit strictly below Ca2+ squares (0.72..1.16)
                y_lo, y_hi = 0.72, 1.16
            else:
                y_lo, y_hi = 0.52, 1.45

        force_tmpl = (
            (series == "ca2" and xv <= 6.0)
            or (series == "mg2" and xv <= 12.0)
            or (series == "his5" and xv <= 0.7)
        )
        near = [] if force_tmpl else [
            p for p in pts_pool
            if abs(cal.px_to_data(p[0], p[1])[0] - xv) / max(xv, 1e-6) <= 0.07
            and y_lo <= cal.px_to_data(p[0], p[1])[1] <= y_hi
        ]
        if len(near) == 1:
            c_f, r_f = near[0][0], near[0][1]
            err = max(3.5, 2.0 + cal.resid_y_px)
        else:
            ra, rb = cal.row_for_y(y_hi, c_est), cal.row_for_y(y_lo, c_est)
            pad_r = 8 if force_tmpl else 20
            r0_b, r1_b = int(round(min(ra, rb))) - pad_r, int(round(max(ra, rb))) + pad_r
            neg_w = 0.65 if (not is_filled and not (series == "his5" and xv <= 0.7)) else 0.35
            sc, bc, br = C.find_template_in_band(dark, (int(c_est) - 34, int(c_est) + 34), (r0_b, r1_b), pos_t, neg_t, neg_w=neg_w)
            if sc > 0.30:
                c_f, r_f = bc, br
                err = max(4.0, 2.5 + cal.resid_y_px)
            elif near:
                best = min(near, key=lambda p: abs(p[0] - c_est))
                c_f, r_f = best[0], best[1]
                err = max(4.5, 3.0 + cal.resid_y_px)
            else:
                cls = C.column_clusters(dark, cal, xv, y_band=(y_lo, y_hi), half_c=3, max_gap=4)
                r_mid = cls[0][0] if cls else 0.5 * (r0_b + r1_b)
                c_f, r_f, err = C.refine_hollow_or_ink_at(dark, cal, c_est, r_mid, win=30)
        out.append(C.make_reading(cal, c_f, r_f, err, f"Fig. 1 {series} (cal resid {cal.resid_y_px:.2f} px)"))
    return out
