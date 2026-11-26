"""Shared calibrated reader for Fig. 1(e): oscillation period T/tau vs beta at P = 1.2."""

from __future__ import annotations

import numpy as np

from _calib import (
    clusters,
    find_frame,
    find_inward_ticks,
    runs_1d,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_fig1e(image: np.ndarray):
    dark, _, h_lines, v_lines = find_frame(image, row_frac=0.50, col_frac=0.50, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]
    _, xt_bot, _, _ = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(8, 26), margin=25, max_width=24, frac=0.60
    )
    # y-ticks include corner ticks at top (40) and bottom (0), so search full [h_top[1], h_bot[2]]
    yt_left = [
        c[0]
        for c in clusters(
            np.where(dark[h_top[1] : h_bot[2] + 1, v_left[2] + 8 : v_left[2] + 26].mean(axis=1) > 0.60)[0]
            + h_top[1]
        )
        if (c[2] - c[1] + 1) <= 24
    ]
    if len(xt_bot) != 3:
        raise RuntimeError(f"Fig. 1(e): expected 3 bottom x-ticks (1, 2, 3), found {xt_bot}")
    if len(yt_left) != 3:
        raise RuntimeError(f"Fig. 1(e): expected 3 left y-ticks (40, 20, 0), found {yt_left}")
    x_axis = verify_axis("linear", [(xt_bot[0], 1.0), (xt_bot[1], 2.0), (xt_bot[2], 3.0)], max_residual_px=1.0)
    y_axis = verify_axis("linear", [(yt_left[0], 40.0), (yt_left[1], 20.0), (yt_left[2], 0.0)], max_residual_px=1.0)
    return h_top, h_bot, v_left, v_right, x_axis, y_axis


def read_pde_period(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the 19 orange numerical simulation period markers T/tau vs beta in Fig. 1(e)."""
    h_top, h_bot, v_left, v_right, x_axis, y_axis = _calibrate_fig1e(image)
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)
    omask = (R > 200) & (G > 45) & (G < 135) & (B < 45)
    omask[: h_top[2] + 5, :] = False
    omask[h_bot[1] - 5 :, :] = False
    omask[:, : v_left[2] + 5] = False
    omask[:, v_right[1] - 5 :] = False

    # Measure reference unobstructed orange circle vertical diameter directly from image
    spans = []
    for bv in at:
        c = int(round(x_axis.to_pixel(float(bv))))
        rs = np.where(omask[:, c - 2 : c + 3].any(axis=1))[0]
        if len(rs) > 20:
            spans.append(int(rs[-1] - rs[0] + 1))
    if not spans:
        raise RuntimeError("Fig. 1(e): no orange circle markers detected")
    diam_ref = float(np.percentile(spans, 80))
    rad_ref = 0.5 * diam_ref

    out = []
    for bv in at:
        c_f = float(x_axis.to_pixel(float(bv)))
        c = int(round(c_f))
        rs = np.where(omask[:, c - 2 : c + 3].any(axis=1))[0]
        if len(rs) < 15:
            out.append({"x": float(bv), "y": None, "found": False, "note": f"no orange marker at col {c}"})
            continue
        r_top, r_bot = int(rs[0]), int(rs[-1])
        span = float(r_bot - r_top + 1)
        if span >= 0.90 * diam_ref:
            # Both top and bottom orange edges unobstructed (even if black line crosses middle)
            r_mid = 0.5 * (r_top + r_bot) + 0.5
            half_px = 3.0 + y_axis.residual_px
            mode_note = f"two-edge span={span:.0f}px"
        else:
            # Bottom clipped by the solid black linear curve; top edge is unobstructed
            r_mid = float(r_top) + 0.5 + rad_ref
            half_px = 5.5 + y_axis.residual_px
            mode_note = f"top-edge+R (span={span:.0f}px < ref {diam_ref:.0f}px)"

        y_val = float(y_axis.to_value(r_mid))
        y_a = float(y_axis.to_value(r_mid - half_px))
        y_b = float(y_axis.to_value(r_mid + half_px))
        lo, hi = min(y_a, y_b), max(y_a, y_b)
        out.append(
            {
                "x": float(bv),
                "y": y_val,
                "lo": lo,
                "hi": hi,
                "px": [c_f, float(r_mid)],
                "px_lo": float(r_mid - half_px),
                "px_hi": float(r_mid + half_px),
                "note": f"Fig. 1(e) PDE symbol, {mode_note}, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
            }
        )
    return out


def read_linear_period(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read the solid black linear-stability curve 2*pi/Im[tau*lambda_+(k_1)] in Fig. 1(e)."""
    h_top, h_bot, v_left, v_right, x_axis, y_axis = _calibrate_fig1e(image)
    R = image[:, :, 0].astype(int)
    G = image[:, :, 1].astype(int)
    B = image[:, :, 2].astype(int)

    kmask = (R < 40) & (G < 40) & (B < 40)
    omask = (R > 200) & (G > 45) & (G < 135) & (B < 45)

    rmin = h_top[2] + 10
    rmax = h_bot[1] - 10
    kmask[:rmin, :] = False
    kmask[rmax + 1 :, :] = False
    kmask[:, : v_left[2] + 10] = False
    kmask[:, v_right[1] - 10 :] = False
    # Exclude "(e)" panel label in upper-right corner
    kmask[: int(round(y_axis.to_pixel(28.0))), int(round(x_axis.to_pixel(2.3))) :] = False

    def pick_linear_curve(runs, col):
        # 1. If a thick black run (>= 24 px) lies strictly below all orange marker pixels in this column
        #    (beta <= 1.55 or beta > 3.05), that unobstructed run is unambiguously the linear curve.
        org_rows = np.where(omask[rmin : rmax + 1, col])[0] + rmin
        r_org_bot = int(org_rows[-1]) if len(org_rows) else rmin
        below = [rn for rn in runs if (rn[1] - rn[0] + 1) >= 24 and rn[1] > r_org_bot + 5]
        if below:
            return below[-1]
        # 2. For beta >= 2.25 where the solid black curve cuts cleanly through the orange markers,
        #    pick the thick black run (>= 22 px) flanked immediately above and below by orange interior.
        flanked = []
        for r0, r1 in runs:
            thick = r1 - r0 + 1
            if thick >= 22:
                has_o_above = r0 >= 10 and bool(omask[max(rmin, r0 - 12) : r0, col].any())
                has_o_below = r1 + 10 < kmask.shape[0] and bool(omask[r1 + 1 : min(rmax, r1 + 13), col].any())
                if has_o_above and has_o_below:
                    flanked.append((r0, r1))
        if flanked:
            return max(flanked, key=lambda t: t[1] - t[0])
        thick_runs = [rn for rn in runs if (rn[1] - rn[0] + 1) >= 24]
        return thick_runs[-1] if thick_runs else None

    return [
        sample_solid_curve_y(
            kmask,
            x_axis,
            y_axis,
            float(bv),
            rmin,
            rmax,
            pick_fn=pick_linear_curve,
            note_prefix="Fig. 1(e) linear curve, ",
        )
        for bv in at
    ]
