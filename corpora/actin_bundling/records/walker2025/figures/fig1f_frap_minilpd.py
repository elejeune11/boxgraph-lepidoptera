"""Calibrated reader for Fig. 1F: mini-Lpd FRAP normalized intensity vs time (mean +/- SD, n=6)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def _read_frap_sd_at_col(
    img: np.ndarray,
    col_f: float,
    curve_mask: np.ndarray,
    black_mask: np.ndarray,
    r_top_lim: int,
    r_bot_lim: int,
    y_axis,
    x_val: float,
    x_axis,
    note_prefix: str = "",
) -> dict:
    """Read a colored FRAP curve point and its vertical black SD error bar at column col_f."""
    c_int = int(round(col_f))
    # Search narrow column window around col_f for curve pixels (green/magenta)
    c0 = max(0, c_int - 2)
    c1 = min(img.shape[1] - 1, c_int + 2)
    sub_c = curve_mask[r_top_lim : r_bot_lim + 1, c0 : c1 + 1]
    ys_curve = np.nonzero(sub_c.any(axis=1))[0] + r_top_lim
    if len(ys_curve) == 0:
        return {"x": float(x_val), "y": None, "note": f"no curve pixels near col {col_f:.1f}"}
    # Group vertical clusters and pick the primary curve stroke
    cs = C.clusters_1d(ys_curve, gap=2)
    # Use the run at the exact column if present
    exact_ys = np.nonzero(curve_mask[r_top_lim : r_bot_lim + 1, c_int])[0] + r_top_lim
    if len(exact_ys) > 0:
        exact_cs = C.clusters_1d(exact_ys, gap=2)
        best_c = max(exact_cs, key=lambda t: t[2] - t[1])
    else:
        best_c = max(cs, key=lambda t: t[2] - t[1])
    y_center = best_c[0]
    half_px = 0.5 * (best_c[2] - best_c[1] + 1)

    # Find vertical black error bar stem + horizontal caps centered near col_f.
    # Each 3-s timepoint in Fig. 1F / 2G has a crisp black stem +/- 4 px around col_f.
    # Find the exact stem column within [-6, +6] of col_f by counting black pixels near y_center.
    win_l = max(0, c_int - 6)
    win_r = min(img.shape[1] - 1, c_int + 6)
    # Restrict vertical search to +/- 140 px around y_center (SD < 0.15 -> < 180 px)
    sr0 = max(r_top_lim, int(round(y_center)) - 200)
    sr1 = min(r_bot_lim, int(round(y_center)) + 200)

    # Walk upward and downward along a narrow 5-px column strip centered at col_f
    # bridging across the curve itself (up to 18 px) and tiny 2-px gaps.
    stem_strip = black_mask[sr0 : sr1 + 1, max(0, c_int - 3) : min(img.shape[1], c_int + 4)]
    stem_on = stem_strip.any(axis=1)

    # Walk UP from top of curve stroke (best_c[1])
    r_start_up = max(sr0, best_c[1] - 1) - sr0
    up_last = None
    miss = 0
    for idx in range(r_start_up, -1, -1):
        if stem_on[idx]:
            up_last = sr0 + idx
            miss = 0
        else:
            miss += 1
            if up_last is None and miss > 14:
                break
            if up_last is not None and miss > 3:
                break

    # Walk DOWN from bottom of curve stroke (best_c[2])
    r_start_dn = min(sr1, best_c[2] + 1) - sr0
    dn_last = None
    miss = 0
    for idx in range(r_start_dn, len(stem_on)):
        if stem_on[idx]:
            dn_last = sr0 + idx
            miss = 0
        else:
            miss += 1
            if dn_last is None and miss > 14:
                break
            if dn_last is not None and miss > 3:
                break

    # Refine cap row center (cap line is ~4-6 px thick): take center of top/bottom cap run
    bar_hi_row = None
    bar_lo_row = None
    note_extra = ""
    if up_last is not None:
        # Cap spans up_last .. up_last + cap_thick
        cap_w = [int(black_mask[r, max(0, c_int - 14) : min(img.shape[1], c_int + 15)].sum())
                 for r in range(up_last, min(sr1, up_last + 10))]
        cap_rows = [up_last + i for i, w in enumerate(cap_w) if w >= 10]
        if cap_rows:
            bar_hi_row = 0.5 * (min(cap_rows) + max(cap_rows)) + 0.5
        else:
            bar_hi_row = float(up_last) + 2.0
    if dn_last is not None:
        cap_w = [int(black_mask[r, max(0, c_int - 14) : min(img.shape[1], c_int + 15)].sum())
                 for r in range(max(sr0, dn_last - 9), dn_last + 1)]
        cap_rows = [max(sr0, dn_last - 9) + i for i, w in enumerate(cap_w) if w >= 10]
        if cap_rows:
            bar_lo_row = 0.5 * (min(cap_rows) + max(cap_rows)) + 0.5
        else:
            bar_lo_row = float(dn_last) - 2.0

    # If error bar is hidden inside curve thickness (e.g., pre-bleach at t=-6), bound by curve stroke
    if bar_hi_row is None and bar_lo_row is None:
        bar_hi_row = float(best_c[1])
        bar_lo_row = float(best_c[2]) + 1.0
        note_extra = "; SD bar hidden within curve thickness"
    elif bar_hi_row is None:
        bar_hi_row = 2.0 * y_center - bar_lo_row
        note_extra = "; upper SD cap occluded, mirrored from lower"
    elif bar_lo_row is None:
        bar_lo_row = 2.0 * y_center - bar_hi_row
        note_extra = "; lower SD cap occluded, mirrored from upper"

    # Slope uncertainty across x calibration residual
    slope = 0.0
    if c_int - 3 >= 0 and c_int + 3 < img.shape[1]:
        yl = np.nonzero(curve_mask[r_top_lim : r_bot_lim + 1, c_int - 3])[0]
        yr = np.nonzero(curve_mask[r_top_lim : r_bot_lim + 1, c_int + 3])[0]
        if len(yl) > 0 and len(yr) > 0:
            slope = ((yr.mean() + r_top_lim) - (yl.mean() + r_top_lim)) / 6.0

    h_px = half_px + 0.5 + abs(slope) * max(x_axis.residual_px, 0.5)
    note = f"{note_prefix}res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px{note_extra}"
    return C.reading_point(
        x_val=x_val,
        col=col_f,
        row=y_center,
        y_axis=y_axis,
        half_px=h_px,
        note=note,
        bar_hi_row=bar_hi_row,
        bar_lo_row=bar_lo_row,
        bar_err_px=2.0,
        clamp_min=0.0,
    )


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    ax = C.find_l_axes(image, dark_thresh=110.0, min_h_frac=0.5, min_v_frac=0.5, max_tick_len=40)
    xt = C.match_regular_ticks(ax["x_ticks"], 5, tol=4.0)
    yt = C.match_regular_ticks(ax["y_ticks"], 7, tol=4.0)
    xa = C.make_axis("linear", list(zip(xt, [-10.0, 40.0, 90.0, 140.0, 190.0])), max_res_px=3.5)
    ya = C.make_axis("linear", list(zip(yt, [1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=3.5)

    rgb = image.astype(np.float64)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    # Green curve in Fig. 1F: G strongly exceeds R and B
    green = (G - R > 40.0) & (G - B > 40.0) & (G > 70.0)
    # Black error bars: low brightness and low saturation (excluding axes)
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    black = (gray < 100.0) & (sat < 25.0)

    r_top = int(round(ya.to_pixel(1.15)))
    r_bot = ax["bot"]["lo"] - 2

    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        out.append(
            _read_frap_sd_at_col(
                image, pc, green, black, r_top, r_bot, ya, float(xv), xa, note_prefix="Fig. 1F mini-Lpd FRAP; "
            )
        )
    return out
