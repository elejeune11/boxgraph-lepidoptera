"""Shared calibrated reader for Fig. 1(d): 1D phase diagram in the (beta, P) plane."""

from __future__ import annotations

import numpy as np

from _calib import (
    clusters,
    find_frame,
    find_inward_ticks,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_fig1d(image: np.ndarray):
    dark, _, h_lines, v_lines = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]
    _, xt_bot, yt_left, _ = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(4, 11), margin=15, max_width=14, frac=0.50
    )
    if len(xt_bot) != 5:
        raise RuntimeError(f"Fig. 1(d): expected 5 bottom x-ticks (-1..3), found {xt_bot}")
    if len(yt_left) != 2:
        raise RuntimeError(f"Fig. 1(d): expected 2 interior left y-ticks (2, 1), found {yt_left}")
    x_axis = verify_axis(
        "linear",
        [(xt_bot[0], -1.0), (xt_bot[1], 0.0), (xt_bot[2], 1.0), (xt_bot[3], 2.0), (xt_bot[4], 3.0)],
        max_residual_px=1.0,
    )
    y_axis = verify_axis(
        "linear",
        [(yt_left[0], 2.0), (yt_left[1], 1.0), (h_bot[0], 0.0)],
        max_residual_px=1.0,
    )
    return h_top, h_bot, v_left, v_right, x_axis, y_axis


def read_linear_boundary(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read solid black linear stability boundary P_c(beta) at unobstructed beta positions."""
    h_top, h_bot, v_left, v_right, x_axis, y_axis = _calibrate_fig1d(image)
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    kmask = (R < 40) & (G < 40) & (B < 40)
    # Exclude borders
    kmask[: h_top[2] + 10, :] = False
    kmask[h_bot[1] - 10 :, :] = False
    kmask[:, : v_left[2] + 5] = False
    kmask[:, v_right[1] - 5 :] = False

    rmin = h_top[2] + 10
    rmax = h_bot[1] - 10

    def pick_thick_lowest(runs, col):
        # Solid linear stability curve is >= 20 px thick vertically (marker rims & box borders are <= 10 px)
        thick = [rn for rn in runs if (rn[1] - rn[0] + 1) >= 18]
        return thick[-1] if thick else None

    return [
        sample_solid_curve_y(
            kmask,
            x_axis,
            y_axis,
            float(bv),
            rmin,
            rmax,
            pick_fn=pick_thick_lowest,
            note_prefix="Fig. 1(d) linear boundary, ",
        )
        for bv in at
    ]


def read_osc_to_stat_boundary(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read nonlinear simulation oscillatory-to-stationary transition P(beta) from marker grid."""
    h_top, h_bot, _, _, x_axis, y_axis = _calibrate_fig1d(image)
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    omask = (R > 200) & (G > 50) & (G < 130) & (B < 40)
    bmask = (B > 180) & (B - R > 60)

    r_top_lim = h_top[2] + 10
    r_osc_floor = int(round(y_axis.to_pixel(1.13)))

    out = []
    for bv in at:
        c_f = float(x_axis.to_pixel(float(bv)))
        c = int(round(c_f))
        o_rows = np.where(omask[r_top_lim:r_osc_floor, c - 12 : c + 13].any(axis=1))[0] + r_top_lim
        b_rows = np.where(bmask[r_top_lim:r_osc_floor, c - 12 : c + 13].any(axis=1))[0] + r_top_lim

        o_blobs = [cl for cl in clusters(o_rows, gap=5) if (cl[2] - cl[1] + 1) >= 14]
        b_blobs = [cl for cl in clusters(b_rows, gap=5) if (cl[2] - cl[1] + 1) >= 14]

        if not o_blobs:
            out.append({"x": float(bv), "y": None, "found": False, "note": "no oscillating marker found"})
            continue
        # Highest P oscillating marker has smallest row coordinate
        r_org_top = o_blobs[0][0]
        # Lowest P stationary marker strictly above the highest oscillating marker has largest row < r_org_top
        b_above = [cl for cl in b_blobs if cl[0] < r_org_top - 15.0]
        if not b_above:
            out.append({"x": float(bv), "y": None, "found": False, "note": "no stationary marker above oscillating"})
            continue
        r_blu_bot = b_above[-1][0]

        P_org = float(y_axis.to_value(r_org_top))
        P_blu = float(y_axis.to_value(r_blu_bot))
        P_mid = 0.5 * (P_org + P_blu)
        r_mid = 0.5 * (r_org_top + r_blu_bot)
        lo, hi = min(P_org, P_blu), max(P_org, P_blu)
        out.append(
            {
                "x": float(bv),
                "y": P_mid,
                "lo": lo,
                "hi": hi,
                "px": [c_f, float(r_mid)],
                "px_lo": float(r_blu_bot),
                "px_hi": float(r_org_top),
                "note": f"Fig. 1(d) osc-to-stat bracket [{lo:.3f}, {hi:.3f}], cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            }
        )
    return out
