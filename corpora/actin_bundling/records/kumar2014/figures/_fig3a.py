"""Shared calibrated reader for Fig. 3(a): 2D ARD phase diagram in the (beta, P) plane."""

from __future__ import annotations

import numpy as np

from _calib import (
    clusters,
    find_frame,
    find_inward_ticks,
    verify_axis,
)


def _calibrate_fig3a(image: np.ndarray):
    dark, _, _, v_lines = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=100.0)
    v_left, v_right = v_lines[0], v_lines[-1]
    # Detect top/bottom plot-box borders strictly inside [v_left, v_right] and aligned with v_left/v_right endpoints
    v_rows = np.where(dark[:, v_left[1] : v_left[2] + 1].any(axis=1))[0]
    r_vtop, r_vbot = int(v_rows[0]), int(v_rows[-1])
    h_box = clusters(np.where(dark[r_vtop - 5 : r_vbot + 6, v_left[1] : v_right[2] + 1].mean(axis=1) > 0.85)[0] + (r_vtop - 5))
    h_top, h_bot = h_box[0], h_box[-1]
    _, xt_bot, yt_left, _ = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(6, 20), margin=25, max_width=18, frac=0.60
    )
    if len(xt_bot) != 3:
        raise RuntimeError(f"Fig. 3(a): expected 3 bottom x-ticks (1.5, 2.0, 2.5), found {xt_bot}")
    if len(yt_left) != 2:
        raise RuntimeError(f"Fig. 3(a): expected 2 left y-ticks (7, 6), found {yt_left}")
    x_axis = verify_axis("linear", [(xt_bot[0], 1.5), (xt_bot[1], 2.0), (xt_bot[2], 2.5)], max_residual_px=1.0)
    y_axis = verify_axis("linear", [(yt_left[0], 7.0), (yt_left[1], 6.0)], max_residual_px=1.0)
    return h_top, h_bot, v_left, v_right, x_axis, y_axis


def _find_legend_box(image: np.ndarray, h_top, h_bot, v_left, v_right):
    """Locate the inner rectangular legend frame in Fig. 3(a) so legend markers are ignored."""
    dark = image.astype(float).mean(axis=2) < 80.0
    sub = dark[h_top[2] + 20 : h_bot[1] - 10, v_left[2] + 20 : v_right[1] - 20]
    # Legend horizontal borders span ~55% of the plot width
    h_leg = clusters(np.where(sub.mean(axis=1) > 0.45)[0])
    if len(h_leg) >= 2:
        r0 = h_leg[0][1] + h_top[2] + 20
        r1 = h_leg[-1][2] + h_top[2] + 20
        # Find left/right vertical borders of the legend box within [r0, r1]
        v_leg = clusters(np.where(dark[r0 : r1 + 1, v_left[2] + 20 : v_right[1] - 20].mean(axis=0) > 0.70)[0])
        if len(v_leg) >= 2:
            c0 = v_leg[0][1] + v_left[2] + 20
            c1 = v_leg[-1][2] + v_left[2] + 20
            return (c0 - 5, r0 - 5, c1 + 5, r1 + 5)
    return None


def read_stability_boundary(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Homogeneous-to-Patterned onset boundary P_c(beta) from Fig. 3(a) simulation grid."""
    h_top, h_bot, v_left, v_right, x_axis, y_axis = _calibrate_fig3a(image)
    leg = _find_legend_box(image, h_top, h_bot, v_left, v_right)

    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    # Patterned (colored) markers: blue hexagons, orange squares, yellow triangles, pink diamonds
    colored = (
        ((B > 180) & (B - R > 60))
        | ((R > 200) & (G >= 55) & (G <= 185) & (B < 40))
        | ((R > 210) & (G >= 155) & (G <= 210) & (B > 180))
    )
    # Dark ink for detecting open white circle rims (Homogeneous)
    dark = (R < 80) & (G < 80) & (B < 80)
    white = (R > 235) & (G > 235) & (B > 235)

    r_top_lim = h_top[2] + 15
    r_bot_lim = h_bot[1] - 15
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        colored[lr0:lr1, lc0:lc1] = False
        dark[lr0:lr1, lc0:lc1] = False

    out = []
    for bv in at:
        c_f = float(x_axis.to_pixel(float(bv)))
        c = int(round(c_f))

        # Colored marker blobs along this beta column
        col_rows = np.where(colored[r_top_lim:r_bot_lim, c - 15 : c + 16].any(axis=1))[0] + r_top_lim
        pat_blobs = [cl for cl in clusters(col_rows, gap=10) if (cl[2] - cl[1] + 1) >= 24]
        if not pat_blobs:
            out.append({"x": float(bv), "y": None, "found": False, "note": "no patterned marker found"})
            continue
        # Lowest-P patterned marker has largest row index
        r_pat_bot = pat_blobs[-1][0]

        # Open white circles (Homogeneous) strictly below r_pat_bot
        # Use strip wider than circle diameter (~66 px) so left/right arcs connect top & bottom rims
        search_r0 = int(round(r_pat_bot)) + 45
        d_rows = np.where(dark[search_r0:r_bot_lim, c - 38 : c + 39].any(axis=1))[0] + search_r0
        hom_blobs = []
        for cl in clusters(d_rows, gap=6):
            h_span = cl[2] - cl[1] + 1
            if 45 <= h_span <= 85:
                rm = int(round(cl[0]))
                # Confirm bright white interior and absence of colored fill
                if white[rm - 8 : rm + 9, c - 8 : c + 9].mean() > 0.85 and not colored[rm - 12 : rm + 13, c - 12 : c + 13].any():
                    hom_blobs.append(cl)
        if not hom_blobs:
            out.append({"x": float(bv), "y": None, "found": False, "note": "no homogeneous marker below patterned"})
            continue
        # Highest-P homogeneous marker has smallest row index
        r_hom_top = hom_blobs[0][0]

        P_pat = float(y_axis.to_value(r_pat_bot))
        P_hom = float(y_axis.to_value(r_hom_top))
        P_mid = 0.5 * (P_pat + P_hom)
        r_mid = 0.5 * (r_pat_bot + r_hom_top)
        lo, hi = min(P_hom, P_pat), max(P_hom, P_pat)
        out.append(
            {
                "x": float(bv),
                "y": P_mid,
                "lo": lo,
                "hi": hi,
                "px": [c_f, float(r_mid)],
                "px_lo": float(r_pat_bot),
                "px_hi": float(r_hom_top),
                "note": f"Fig. 3(a) onset bracket [{lo:.3f}, {hi:.3f}], cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            }
        )
    return out


def read_drifting_to_honeycomb(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Drifting-square (pink diamond) to Oscillating-honeycomb (yellow triangle) boundary beta(P)."""
    h_top, h_bot, v_left, v_right, x_axis, y_axis = _calibrate_fig3a(image)
    leg = _find_legend_box(image, h_top, h_bot, v_left, v_right)

    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    pink = (R > 210) & (G >= 155) & (G <= 210) & (B > 180)
    yellow = (R > 200) & (G >= 120) & (G <= 185) & (B < 40)

    c_left_lim = v_left[2] + 15
    c_right_lim = v_right[1] - 15
    if leg is not None:
        lc0, lr0, lc1, lr1 = leg
        pink[lr0:lr1, lc0:lc1] = False
        yellow[lr0:lr1, lc0:lc1] = False

    out = []
    for pv in at:
        r_f = float(y_axis.to_pixel(float(pv)))
        r = int(round(r_f))

        p_cols = np.where(pink[r - 18 : r + 19, c_left_lim:c_right_lim].any(axis=0))[0] + c_left_lim
        y_cols = np.where(yellow[r - 18 : r + 19, c_left_lim:c_right_lim].any(axis=0))[0] + c_left_lim

        p_blobs = [cl for cl in clusters(p_cols, gap=10) if (cl[2] - cl[1] + 1) >= 20]
        y_blobs = [cl for cl in clusters(y_cols, gap=10) if (cl[2] - cl[1] + 1) >= 20]

        if not p_blobs or not y_blobs:
            out.append({"x": float(pv), "y": None, "found": False, "note": f"missing pink/yellow markers at P={pv}"})
            continue

        # Rightmost pink diamond and leftmost yellow triangle
        c_pink_max = p_blobs[-1][0]
        c_yel_min = y_blobs[0][0]
        b_pink = float(x_axis.to_value(c_pink_max))
        b_yel = float(x_axis.to_value(c_yel_min))
        b_mid = 0.5 * (b_pink + b_yel)
        c_mid = 0.5 * (c_pink_max + c_yel_min)
        lo, hi = min(b_pink, b_yel), max(b_pink, b_yel)
        out.append(
            {
                "x": float(pv),
                "y": b_mid,
                "lo": lo,
                "hi": hi,
                "px": [float(c_mid), r_f],
                "note": f"Fig. 3(a) drift->honeycomb bracket [{lo:.3f}, {hi:.3f}], cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            }
        )
    return out
