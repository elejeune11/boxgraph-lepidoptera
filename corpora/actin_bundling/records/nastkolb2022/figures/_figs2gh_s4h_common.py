"""Shared reader for black-square-with-error-bars VASP titration panels: Fig. S2G, Fig. S2H, and Fig. S4H (1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import math
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    LinearCal,
    detect_ticks_along_spine,
    find_black_spines,
    select_major_ticks,
)


def read_vasp_marker_errorbars(
    image: np.ndarray,
    x_ticks_expected: list[float],
    y_ticks_expected: list[float],
    n_points: int,
    x_coord_min_offset: int = -10,
    y_coord_max_offset: int = 10,
) -> list[dict]:
    """Read black square markers and their vertical error bars (with T-caps) on a clean OriginLab panel."""
    im = image.astype(float)
    H, W, _ = im.shape
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=48, min_len=12, coord_min=vl[0] + x_coord_min_offset, coord_max=W - 10)
    xt = select_major_ticks(xt_raw, expected=len(x_ticks_expected), frac=0.72)
    x_cal = LinearCal(xt, x_ticks_expected, max_res_px=2.5)

    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=12, coord_min=50, coord_max=hr[1] + y_coord_max_offset)
    yt = select_major_ticks(yt_raw, expected=len(y_ticks_expected), frac=0.72)
    y_cal = LinearCal(yt, y_ticks_expected, max_res_px=2.5)

    # Interior plot region strictly excluding spines
    r_top = 40
    r_bot = hr[0] - 4
    c_left = vl[1] + 10
    c_right = W - 15
    sub = np.zeros((H, W), dtype=bool)
    sub[r_top:r_bot, c_left:c_right] = dark[r_top:r_bot, c_left:c_right]

    # Black square markers are ~38-48 px wide; opening with 20x20 erases all 6-10 px stems and caps cleanly!
    squares = ndimage.binary_opening(sub, structure=np.ones((20, 20), dtype=bool))
    lab, n = ndimage.label(squares)

    pts = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 400:
            continue
        r0, r1 = int(rr.min()), int(rr.max())
        c0, c1 = int(cc.min()), int(cc.max())
        pts.append((0.5 * (c0 + c1) + 0.5, 0.5 * (r0 + r1) + 0.5, r0, r1, c0, c1, len(rr)))
    pts.sort(key=lambda p: p[0])
    if len(pts) != n_points:
        raise RuntimeError(f"Expected {n_points} square markers, found {len(pts)}: {pts}")

    out = []
    for c_mid, r_mid, r0, r1, c0, c1, area in pts:
        ci = int(round(c_mid))
        # Walk up from r0 - 1 and down from r1 + 1 along central stem (ci ± 5)
        stem = sub[:, max(c_left, ci - 5) : min(c_right, ci + 6)].any(axis=1)
        cap_half = max(12, int(0.45 * (c1 - c0)))
        cap_frac = sub[:, max(c_left, ci - cap_half) : min(c_right, ci + cap_half + 1)].mean(axis=1)

        r_up = r0 - 1
        last_up = r0
        miss = 0
        while r_up >= r_top and miss <= 4:
            if stem[r_up]:
                last_up = r_up
                miss = 0
            else:
                miss += 1
            r_up -= 1

        r_dn = r1 + 1
        last_dn = r1
        miss = 0
        while r_dn < r_bot and miss <= 4:
            if stem[r_dn]:
                last_dn = r_dn
                miss = 0
            else:
                miss += 1
            r_dn += 1

        hidden_note = []
        if last_up >= r0 - 3:
            r_bar_hi = float(r0) + 0.5
            hidden_note.append("upper")
        else:
            cu = np.where(cap_frac[max(r_top, last_up - 3) : min(r_bot, last_up + 12)] > 0.60)[0]
            r_bar_hi = float(cu.mean()) + 0.5 + max(r_top, last_up - 3) if len(cu) > 0 else float(last_up) + 0.5

        if last_dn <= r1 + 3:
            r_bar_lo = float(r1) + 0.5
            hidden_note.append("lower")
        else:
            cd = np.where(cap_frac[max(r_top, last_dn - 12) : min(r_bot, last_dn + 4)] > 0.60)[0]
            r_bar_lo = float(cd.mean()) + 0.5 + max(r_top, last_dn - 12) if len(cd) > 0 else float(last_dn) + 0.5

        xv = x_cal.to_val(c_mid)
        yv = y_cal.to_val(r_mid)
        bar_hi = max(yv, y_cal.to_val(r_bar_hi))
        bar_lo = min(yv, y_cal.to_val(r_bar_lo))
        err_y = abs(y_cal.b) * (2.5 + y_cal.res_px)
        bar_pad = abs(y_cal.b) * (3.0 + y_cal.res_px)

        # Snap x to exact incubated VASP concentrations (0.05, 0.125, 0.25, 0.5, 1.0, 2.0 uM) if within 0.01 uM
        known = [0.05, 0.125, 0.25, 0.5, 1.0, 2.0]
        nearest = min(known, key=lambda v: abs(v - xv))
        if abs(nearest - xv) < 0.015:
            x_rep = nearest
        else:
            x_rep = round(xv, 3)

        note = f"marker={area}px, x_raw={xv:.4f}, bar=[{bar_lo:.4g}, {bar_hi:.4g}], rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px"
        if hidden_note:
            note += f"; error bar ({' and '.join(hidden_note)}) within marker: taken as at most marker extent"
        out.append({
            "x": x_rep,
            "y": yv,
            "lo": yv - err_y,
            "hi": yv + err_y,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo - bar_pad,
            "bar_hi_outer": bar_hi + bar_pad,
            "px": [c_mid, r_mid],
            "px_lo": r_bar_hi,
            "px_hi": r_bar_lo,
            "note": note,
        })
    return out
