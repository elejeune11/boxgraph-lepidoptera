"""Calibrated readers for SI Figs. S4(g), S4(h), and S5(c)."""

from __future__ import annotations

import numpy as np
import _calib as C


def _read_static_and_triangles(
    image: np.ndarray,
    at: list[float],
    x_ticks_vals: list[float],
    y_ticks_vals: list[float],
    static_max_x: float,
    tri_ref_xs: list[float],
    label: str,
) -> list[dict]:
    """Shared reader for Mathematica plots with static blue squares at y=0 and traveling/oscillatory red triangles."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, len(x_ticks_vals)), x_ticks_vals)))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, len(y_ticks_vals)), y_ticks_vals)))

    bmask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    r_min, r_max = ht[2] + 4, hb[1] - 4
    y0_px = float(y_axis.to_pixel(0.0))

    # Calibrate blue square baseline offset from static squares at x = 0.4, 0.8, 1.0
    sq_mids = []
    for xs in (0.4, 0.8, 1.0):
        cs = int(round(x_axis.to_pixel(xs)))
        rs = np.nonzero(bmask[r_min:r_max, cs])[0] + r_min
        if len(rs):
            sq_mids.append(0.5 * (rs.min() + rs.max()) + 0.5)
    sq_offset = float(np.median(sq_mids) - y0_px) if sq_mids else 0.0

    # Measure exact uncontaminated triangle height H_tri as the lower-quartile/mode of triangle runs
    # (taking runs in [45, 72] px so curve-merged runs > 72 px never bias H_tri upward).
    hs = []
    for xv in tri_ref_xs:
        cx = int(round(x_axis.to_pixel(xv)))
        for cc in range(cx - 5, cx + 1):
            if 0 <= cc < rmask.shape[1]:
                rs = C.runs_1d(rmask[r_min:r_max, cc], gap=2)
                for a, b in rs:
                    h = b - a + 1
                    if 48 <= h <= 72:
                        hs.append(h)
    H_tri = float(np.median(hs)) if hs else 65.0

    out = []
    for xv in at:
        cf = float(x_axis.to_pixel(xv))
        ci = int(round(cf))
        if xv <= static_max_x:
            rs = np.nonzero(bmask[r_min:r_max, ci])[0] + r_min
            if len(rs) == 0:
                out.append({"x": float(xv), "y": None, "note": f"no blue square at x={xv:g}"})
                continue
            r_mid = 0.5 * (rs.min() + rs.max()) + 0.5 - sq_offset
            out.append(C.make_marker_reading(
                xv, r_mid, cf, x_axis, y_axis, half_px=4.0,
                note=f"{label} static blue square (sq_off={sq_offset:+.1f}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
            ))
        else:
            # Find highest apex (minimum r_top) among runs of height >= 0.7 * H_tri in [ci-5, ci+1].
            # Special case: if a fitted curve runs separately ABOVE the triangle (like Fig. S5c at x=1.6),
            # the curve run is thin (h < 35 px) and automatically rejected by h >= 0.7 * H_tri!
            best_top = None
            best_c = ci
            for cc in range(ci - 5, ci + 2):
                if 0 <= cc < rmask.shape[1]:
                    rs = C.runs_1d(rmask[r_min:r_max, cc], gap=2)
                    for a, b in rs:
                        if (b - a + 1) >= 0.7 * H_tri:
                            rt = r_min + a
                            if best_top is None or rt < best_top:
                                best_top = rt
                                best_c = cc
            if best_top is None:
                out.append({"x": float(xv), "y": None, "note": f"no red triangle at x={xv:g}"})
                continue
            r_anc = best_top + H_tri / 3.0 + 0.5
            out.append(C.make_marker_reading(
                xv, r_anc, cf, x_axis, y_axis, half_px=6.5,
                note=f"{label} red triangle (apex={best_top}, H_tri={H_tri:.0f}px, c={best_c}, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
            ))
    return out


def read_figS4g_sim(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_static_and_triangles(
        image, at, [0.0, 1.0, 2.0, 3.0], [0.02, 0.01, 0.0],
        static_max_x=1.41, tri_ref_xs=[1.5, 2.0, 2.4, 2.6, 2.8, 3.0], label="Fig. S4(g)",
    )


def read_figS4g_theory(image: np.ndarray, at: list[float]) -> list[dict]:
    """Sample solid red 1D theoretical curve v = sqrt(delta^2 - delta_c^2) between markers in Fig. S4(g)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 4), [0.0, 1.0, 2.0, 3.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.02, 0.01, 0.0])))
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    r_min, r_max = ht[2] + 4, hb[1] - 4
    return [
        C.sample_curve_y(rmask, x_axis, y_axis, xv, r_min, r_max, note_prefix="Fig. S4(g) red solid curve: ")
        for xv in at
    ]


def read_figS4h_sim(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_static_and_triangles(
        image, at, [0.0, 1.0, 2.0, 3.0], [0.003, 0.002, 0.001, 0.0],
        static_max_x=1.45, tri_ref_xs=[1.7, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0], label="Fig. S4(h)",
    )


def read_figS5c_sim(image: np.ndarray, at: list[float]) -> list[dict]:
    return _read_static_and_triangles(
        image, at, [0.0, 2.0, 4.0], [0.02, 0.01, 0.0],
        static_max_x=1.45, tri_ref_xs=[1.6, 1.8, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4], label="Fig. S5(c)",
    )


def read_figS5c_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    """Sample solid red fitted curve v = sqrt(delta^2 - (1.5*kappa)^2) between triangles in Fig. S5(c)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 3), [0.0, 2.0, 4.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.02, 0.01, 0.0])))
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    r_min, r_max = ht[2] + 4, hb[1] - 4
    out = []
    for xv in at:
        # Between triangles (at odd tenths like 1.7, 2.1, 2.5, 2.9, 3.3, 3.7), the solid red line sits cleanly
        # below the triangle corners (bottommost thin run of height <= 32 px).
        rec = C.sample_curve_y(
            rmask, x_axis, y_axis, xv, r_min, r_max,
            pick_fn=lambda rs: max(rs, key=lambda ab: ab[1]),
            search_dc=3,
            note_prefix="Fig. S5(c) red fitted curve: ",
        )
        out.append(rec)
    return out
