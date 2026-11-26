"""Figure 3a & 3b (p. 3 / 8558) readers for tang1996.

Fig. 3a: Scattering intensity vs [Lys18] (uM) for 6.0 mg/ml (squares), 2.0 mg/ml (diamonds),
         and 0.1 mg/ml (three overlapping replicate open-circle traces) F-actin in 150 mM KCl.
Fig. 3b: Bundling onset [Lys18] (uM) vs [Actin] (uM) — 6 measured open-circle points and
         the plotted linear fit line [Lys18] = 19.4 + 0.26 x [Actin].
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def read_fig3b_points(image: np.ndarray) -> list[dict]:
    cal = P.calibrate("fig3b", image)
    ex = [cal.data_box_to_px(5.0, 48.0, 115.0, 59.0)]  # equation text box
    pts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
    h_px = max(3.5, 2.0 + cal.resid_y_px)
    return [C.make_reading(cal, c, r, h_px, f"Fig. 3b onset open circle (cal resid {cal.resid_y_px:.2f} px)") for c, r, _ in pts]


def read_fig3b_fit_line(image: np.ndarray, at: list[float]) -> list[dict]:
    cal = P.calibrate("fig3b", image)
    ex = [cal.data_box_to_px(5.0, 48.0, 115.0, 59.0)]
    return C.read_line_at(image, cal, at, y_band=(15.0, 58.0), exclude_boxes_px=ex, note_prefix="Fig. 3b linear fit line")


def read_fig3a(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    cal = P.calibrate("fig3a", image)
    ex = [cal.data_box_to_px(-15.0, 1.5, 48.0, 9.5)]  # legend box
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    if series == "actin_6p0":
        # Open squares (6.0 mg/ml): cavity area ~ 2000-2700 px (vs ~1200-1350 for diamonds);
        # start from top marker at x=0 and trace across `at`, snapping to large-area square markers
        # or square template at touching crossing near x=39, 67, 83.
        pos_sq, neg_sq = C.make_geom_template("square", radius_px=29.0, stroke_px=7.0)
        sq_pts = [p for p in hpts if p[2] >= 1850 and cal.px_to_data(p[0], p[1])[1] >= 0.65]
        # Match each x in `at` to closest detected square or template search near traced curve
        out = []
        last_r = cal.row_for_y(0.80, cal.col_for_x(0.0, 0.5 * (cal.inner_box[1] + cal.inner_box[3])))
        for xv in sorted(float(x) for x in at):
            c_x = cal.col_for_x(xv, last_r)
            near = [p for p in sq_pts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 2.2]
            if near:
                best = min(near, key=lambda p: abs(p[1] - last_r))
                c_f, r_f = best[0], best[1]
                err = max(3.5, 2.0 + cal.resid_y_px)
            else:
                # Template search near traced vertical position
                sc, c_f, r_f = C.find_template_in_band(dark, (int(c_x) - 35, int(c_x) + 35), (int(last_r) - 180, int(last_r) + 120), pos_sq, neg_sq, neg_w=0.35)
                err = max(4.5, 3.0 + cal.resid_y_px)
            last_r = r_f
            out.append(C.make_reading(cal, c_f, r_f, err, f"Fig. 3a 6.0 mg/ml square (cal resid {cal.resid_y_px:.2f} px)"))
        return out

    if series == "actin_2p0":
        # Open diamonds (2.0 mg/ml): cavity area ~ 1150-1450 px, starting around y=0.29 at x=0
        # and rising above 6.0 mg/ml for x in [40..62], touching at x=66.7 and 83.3.
        pos_d, neg_d = C.make_geom_template("diamond", radius_px=27.0, stroke_px=7.0)
        dia_pts = [p for p in hpts if 1050 <= p[2] <= 1500 and cal.px_to_data(p[0], p[1])[1] >= 0.22]
        out = []
        last_r = cal.row_for_y(0.30, cal.col_for_x(0.0, 0.5 * (cal.inner_box[1] + cal.inner_box[3])))
        for xv in sorted(float(x) for x in at):
            c_x = cal.col_for_x(xv, last_r)
            near = [p for p in dia_pts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 2.2]
            if xv < 25.0:
                near = [p for p in near if cal.px_to_data(p[0], p[1])[1] >= 0.25]
            elif xv >= 45.0:
                near = [p for p in near if cal.px_to_data(p[0], p[1])[1] >= 1.5]
            if near:
                best = min(near, key=lambda p: abs(p[1] - last_r))
                c_f, r_f = best[0], best[1]
                err = max(3.5, 2.0 + cal.resid_y_px)
            else:
                # Search near last_r (monotonic increasing curve -> r_f <= last_r + 35)
                up_win = 100 if xv >= 60.0 else 240
                sc, c_f, r_f = C.find_template_in_band(dark, (int(c_x) - 35, int(c_x) + 35), (int(last_r) - up_win, int(last_r) + 35), pos_d, neg_d, neg_w=0.25)
                err = max(4.5, 3.0 + cal.resid_y_px)
            last_r = r_f
            out.append(C.make_reading(cal, c_f, r_f, err, f"Fig. 3a 2.0 mg/ml diamond (cal resid {cal.resid_y_px:.2f} px)"))
        return out

    if series == "actin_0p1":
        # Three closely overlapping replicate open-circle curves at 0.1 mg/ml (bottom group,
        # rising from ~0.025 at x=0 to ~0.37 plateau for x >= 33 uM).
        out = []
        for xv in sorted(float(x) for x in at):
            y_hi_band = 0.22 if xv <= 25.0 else (0.31 if xv <= 29.0 else 0.45)
            y_lo_band = 0.015 if xv <= 25.0 else (0.20 if xv <= 29.0 else 0.28)
            c_mid = cal.col_for_x(xv, 0.5 * (cal.inner_box[1] + cal.inner_box[3]))
            ra = cal.row_for_y(y_hi_band, c_mid)
            rb = cal.row_for_y(y_lo_band, c_mid)
            r0_b, r1_b = max(cal.inner_box[1] + 8, int(round(min(ra, rb)))), min(cal.inner_box[3] - 8, int(round(max(ra, rb))))
            ci = int(round(c_mid))
            near = [p for p in hpts if abs(cal.px_to_data(p[0], p[1])[0] - xv) <= 2.2 and y_lo_band <= cal.px_to_data(p[0], p[1])[1] <= y_hi_band]
            if len(near) >= 2:
                cs = [p[0] for p in near]
                rs = [p[1] for p in near]
                c_m = float(np.mean(cs))
                r_m = 0.5 * (min(rs) + max(rs))
                half_px = max(4.0, 0.5 * (max(rs) - min(rs)) + 2.5 + cal.resid_y_px)
            else:
                strip = dark[r0_b:r1_b, max(cal.inner_box[0] + 4, ci - 16):min(cal.inner_box[2] - 4, ci + 17)]
                row_cnt = strip.sum(axis=1)
                hits = np.nonzero(row_cnt >= 9)[0]
                if len(hits) == 0:
                    hits = np.nonzero(row_cnt >= 1)[0]
                r_top = r0_b + float(hits[0]) + 22.0
                r_bot = r0_b + float(hits[-1]) - 22.0
                r_m = 0.5 * (r_top + r_bot)
                half_px = max(4.0, (0.5 * (r_bot - r_top) if r_bot > r_top else 0.0) + 3.0 + cal.resid_y_px)
                c_m = cal.col_for_x(xv, r_m)
            out.append(C.make_reading(cal, c_m, r_m, half_px, f"Fig. 3a 0.1 mg/ml 3-replicate envelope (cal resid {cal.resid_y_px:.2f} px)"))
        return out

    raise ValueError(f"Unknown Fig. 3a series: {series}")
