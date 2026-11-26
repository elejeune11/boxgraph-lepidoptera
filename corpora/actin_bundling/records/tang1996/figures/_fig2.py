"""Figure 2a & 2b (p. 3 / 8558) readers for tang1996.

Fig. 2a: Scattering intensity vs [Cu2+] (mM) with (+20mM MgCl2, filled circles)
         and without (no MgCl2, hollow circles) 20 mM MgCl2.
Fig. 2b: Scattering intensity vs [Mg2+] (mM) preincubated with 0 (downtri),
         10 uM (square), 50 uM (circle), and 100 uM (diamond) CuCl2.
"""

from __future__ import annotations

import numpy as np
import _common as C
import _panels as P


def read_fig2a(image: np.ndarray, series: str) -> list[dict]:
    cal = P.calibrate("fig2a", image)
    ex = [
        cal.data_box_to_px(0.50, 1.80, 0.96, 3.00),  # "+20mM MgCl2"
        cal.data_box_to_px(0.55, 0.50, 0.92, 0.95),  # "no MgCl2"
    ]
    if series == "cu_with_20mM_mg":
        pts = C.detect_filled_markers(image, cal, exclude_boxes_px=ex, erode_2x=6, min_edt_px=6.0)
        h_px = max(3.5, 2.0 + cal.resid_y_px)
        return [C.make_reading(cal, c, r, h_px, f"Fig. 2a +20mM MgCl2 solid circle (cal resid {cal.resid_y_px:.2f} px)") for c, r, _ in pts]
    if series == "cu_no_mg":
        pts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)
        h_px = max(3.5, 2.0 + cal.resid_y_px)
        return [C.make_reading(cal, c, r, h_px, f"Fig. 2a no MgCl2 open circle (cal resid {cal.resid_y_px:.2f} px)") for c, r, _ in pts]
    raise ValueError(f"Unknown Fig. 2a series: {series}")


def _fig2b_pick_curve(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 2b curve at its experimental [Mg2+] positions `at` using marker cavity detection,
    curve-order disambiguation on separated columns, and bidirectional thin-line tracing."""
    cal = P.calibrate("fig2b", image)
    ex = [cal.data_box_to_px(-1.5, 1.05, 22.0, 4.0)]  # legend box
    dark = C.mask_boxes(C._to_dark(image), ex)
    hpts = C.detect_hollow_markers(image, cal, exclude_boxes_px=ex)

    # Helper: snap to a detected hollow marker within dx_data <= 1.2 and within y_band, else trace/column
    def marker_in(x_val: float, y_lo: float, y_hi: float):
        cands = []
        for c, r, area in hpts:
            xv, yv = cal.px_to_data(c, r)
            if abs(xv - x_val) <= 1.0 and y_lo <= yv <= y_hi:
                cands.append((c, r, xv, yv, area))
        return sorted(cands, key=lambda t: t[3])  # ascending y

    # Each series starts from an unambiguous separated column (x=40 for no_cu/10uM; x=15 for 50uM/100uM)
    # and traces cleanly across every x in `at`!
    xs = sorted(float(x) for x in at)

    def trace_from_seed(x_seed: float, r_seed: float, c_seed: float) -> list[dict]:
        # Trace both left and right from seed across xs
        idx_s = min(range(len(xs)), key=lambda i: abs(xs[i] - x_seed))
        res_px: dict[int, tuple[float, float, float]] = {}
        res_px[idx_s] = (c_seed, r_seed, max(3.5, 2.0 + cal.resid_y_px))

        # Trace right
        cc, rr = c_seed, r_seed
        for i in range(idx_s + 1, len(xs)):
            c_tgt = cal.col_for_x(xs[i], rr)
            _, r_tr = C.trace_step_along_curve(dark, cal, cc, rr, c_tgt, search_r_px=65, step_c=16.0)
            c_ref, r_ref, err_px = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=36)
            res_px[i] = (c_ref, r_ref, err_px)
            cc, rr = c_ref, r_ref

        # Trace left
        cc, rr = c_seed, r_seed
        for i in range(idx_s - 1, -1, -1):
            c_tgt = cal.col_for_x(xs[i], rr)
            _, r_tr = C.trace_step_along_curve(dark, cal, cc, rr, c_tgt, search_r_px=65, step_c=16.0)
            c_ref, r_ref, err_px = C.refine_hollow_or_ink_at(dark, cal, c_tgt, r_tr, win=36)
            res_px[i] = (c_ref, r_ref, err_px)
            cc, rr = c_ref, r_ref

        return [C.make_reading(cal, *res_px[i], f"Fig. 2b {series} (cal resid {cal.resid_y_px:.2f} px)") for i in range(len(xs))]

    if series == "mg_no_cu":
        # Top curve at x = 40 mM (y > 1.6) AND second from bottom at x = 15..20 mM
        # Note the steep jump between 25 and 30 mM; trace plateau [30..50] from x=40 top marker
        # and baseline [0..25] from x=15 second-lowest marker!
        top40 = marker_in(40.0, 1.6, 3.0)
        mid15 = marker_in(15.0, 0.09, 0.14)
        out_hi = [r for r in trace_from_seed(40.0, top40[-1][1], top40[-1][0]) if r["x"] >= 28.0]
        out_lo = [r for r in trace_from_seed(15.0, mid15[0][1], mid15[0][0]) if r["x"] <= 27.0]
        return sorted(out_lo + out_hi, key=lambda d: d["x"])

    if series == "mg_10uM_cu":
        # 2nd from top at x = 40 mM (1.0 < y < 1.6) and bottom curve at x = 15 mM (0.06 < y < 0.09)
        mid40 = marker_in(40.0, 1.0, 1.6)
        bot15 = marker_in(15.0, 0.06, 0.09)
        sq30 = marker_in(30.0, 0.85, 1.45)
        out_hi = [r for r in trace_from_seed(40.0, mid40[0][1], mid40[0][0]) if r["x"] >= 32.0]
        if sq30:
            out_hi.append(C.make_reading(cal, sq30[0][0], sq30[0][1], max(3.5, 2.0 + cal.resid_y_px), f"Fig. 2b mg_10uM_cu square (cal resid {cal.resid_y_px:.2f} px)"))
        out_lo = [r for r in trace_from_seed(15.0, bot15[0][1], bot15[0][0]) if r["x"] <= 27.0]
        return sorted(out_lo + out_hi, key=lambda d: d["x"])

    if series == "mg_50uM_cu":
        # Unambiguous rising point at x = 15 mM (0.2 < y < 0.5)
        p15 = marker_in(15.0, 0.2, 0.5)[0]
        return trace_from_seed(15.0, p15[1], p15[0])

    if series == "mg_100uM_cu":
        # Unambiguous top point at x = 15 mM (0.5 < y < 0.9); baseline diamonds at x=0, 5 lie in [0.06, 0.13]
        p15 = marker_in(15.0, 0.5, 0.9)[0]
        out = trace_from_seed(15.0, p15[1], p15[0])
        pos_d, neg_d = C.make_geom_template("diamond", radius_px=31.0, stroke_px=7.5)
        for idx, rd in enumerate(out):
            xv = xs[idx]
            if xv <= 6.0:
                c_e = cal.col_for_x(xv, cal.row_for_y(0.085, cal.inner_box[0] + 100))
                r_top = int(cal.row_for_y(0.13, c_e))
                r_bot = int(cal.row_for_y(0.06, c_e))
                sc, bc, br = C.find_template_in_band(dark, (int(c_e) - 45, int(c_e) + 45), (r_top, r_bot), pos_d, neg_d, neg_w=0.30)
                out[idx] = C.make_reading(cal, bc, br, max(3.5, 2.0 + cal.resid_y_px), f"Fig. 2b mg_100uM_cu baseline diamond (cal resid {cal.resid_y_px:.2f} px)")
            elif xv >= 28.0:
                c_e, r_e = rd["px"]
                sc, bc, br = C.find_template_in_band(dark, (int(c_e) - 45, int(c_e) + 45), (int(r_e) - 80, int(r_e) + 60), pos_d, neg_d, neg_w=0.30)
                if sc > 0.50:
                    out[idx] = C.make_reading(cal, bc, br, max(3.5, 2.0 + cal.resid_y_px), f"Fig. 2b mg_100uM_cu diamond (cal resid {cal.resid_y_px:.2f} px)")
        return out

    raise ValueError(f"Unknown Fig. 2b series: {series}")


def read_fig2b(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    return _fig2b_pick_curve(image, at, series)
