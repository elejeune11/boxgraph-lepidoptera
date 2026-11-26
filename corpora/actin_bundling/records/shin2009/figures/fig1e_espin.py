"""Calibrated reader for Fig. 1(e) (p. 2): measured espin-actin bundle overtwist vs R (7 black squares)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import black_mask, calibrate_fig1e, runs_1d, value_interval


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB) = calibrate_fig1e(image)
    blk = black_mask(image, gray_max=80.0, sat_max=35.0)
    r_ch = image[:, :, 0].astype(float)
    g_ch = image[:, :, 1].astype(float)
    b_ch = image[:, :, 2].astype(float)
    red = (r_ch > 140.0) & (r_ch - g_ch > 60.0) & (r_ch - b_ch > 60.0)

    # Interior mask excluding legend region (R > 0.23, overtwist < 0.65 deg)
    interior = np.zeros_like(blk)
    interior[rT + 25 : rB - 35, cL + 25 : cR - 25] = blk[rT + 25 : rB - 35, cL + 25 : cR - 25]
    c_leg = int(round(x_axis.to_pixel(0.23)))
    r_leg = int(round(y_axis.to_pixel(0.65)))
    interior[r_leg:, c_leg:] = False

    # Strip thin error-bar stems (~20 px) while preserving filled squares (~95-106 x 74-82 px)
    opened = ndimage.binary_opening(interior, structure=np.ones((30, 30), bool))
    lab, n_comp = ndimage.label(opened)

    # First pass: collect isolated single squares to measure nominal marker dimensions
    isolated = []
    merged = []
    for i in range(1, n_comp + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 1000:
            continue
        w = int(cc.max() - cc.min() + 1)
        h = int(rr.max() - rr.min() + 1)
        if w <= 125 and h <= 100:
            isolated.append((int(cc.min()), int(cc.max()), int(rr.min()), int(rr.max()), len(rr)))
        else:
            merged.append((int(cc.min()), int(cc.max()), int(rr.min()), int(rr.max()), i))

    if not isolated:
        raise RuntimeError("Fig. 1(e) espin: no isolated black square markers found")
    nom_w = float(np.median([c1 - c0 + 1 for c0, c1, _, _, _ in isolated]))
    nom_h = float(np.median([r1 - r0 + 1 for _, _, r0, r1, _ in isolated]))

    # Decompose any touching-square cluster using the measured single-square width/height
    markers: list[tuple[float, float, int, int, str]] = []
    for c0, c1, r0, r1, area in isolated:
        markers.append((0.5 * (c0 + c1) + 0.5, 0.5 * (r0 + r1) + 0.5, r0, r1, f"isolated square {c1-c0+1}x{r1-r0+1}px"))

    for c0, c1, r0, r1, cid in merged:
        sub = (lab == cid)
        # Check for an upper offset square protruding above the bottom row of height nom_h
        r_split = r1 - int(round(nom_h)) - 3
        if r0 < r_split - 20:
            up_rr, up_cc = np.nonzero(sub[:r_split, :])
            if len(up_cc) > 0:
                uc0, uc1 = int(up_cc.min()), int(up_cc.max())
                ur0 = int(up_rr.min())
                ur1 = ur0 + int(round(nom_h)) - 1
                markers.append((0.5 * (uc0 + uc1) + 0.5, 0.5 * (ur0 + ur1) + 0.5, ur0, ur1,
                                f"upper square of cluster ({uc1-uc0+1}x{ur1-ur0+1}px)"))
        # Check the lower horizontal row near r1 - 20 (below the upper square's error-bar cap)
        probe_r = r1 - 20
        lo_cols = np.where(sub[probe_r, :])[0]
        if len(lo_cols) > 0:
            lc0, lc1 = int(lo_cols.min()), int(lo_cols.max())
            span_w = lc1 - lc0 + 1
            lr0 = r1 - int(round(nom_h)) + 1
            if span_w > 1.35 * nom_w:
                # Two touching horizontal squares: anchor left and right centers by half nominal width
                markers.append((lc0 + 0.5 * nom_w, 0.5 * (lr0 + r1) + 0.5, lr0, r1,
                                f"left square of touching pair (span {span_w}px)"))
                markers.append((lc1 + 1.0 - 0.5 * nom_w, 0.5 * (lr0 + r1) + 0.5, lr0, r1,
                                f"right square of touching pair (span {span_w}px)"))
            else:
                markers.append((0.5 * (lc0 + lc1) + 0.5, 0.5 * (lr0 + r1) + 0.5, lr0, r1,
                                f"lower square of cluster ({span_w}px)"))

    cal_y = max(y_axis.residual_px, 1.0)
    center_err_px = 2.5 + cal_y
    cap_err_px = 2.0 + cal_y

    readings = []
    for c_c, r_c, r_sq_top, r_sq_bot, m_note in markers:
        c_int = int(round(c_c))
        # Check if a red fascin marker sits directly above (e.g. R = 0.33); if so, any upward
        # vertical stem belongs to the red fascin marker's lower error bar, not to this espin square.
        red_above = red[max(rT + 10, r_sq_top - 200) : max(rT + 10, r_sq_top - 25),
                        max(cL, c_int - 50) : min(cR, c_int + 51)].any()

        # Trace thin vertical error-bar stem above and below the square
        has_up = (not red_above) and blk[r_sq_top - 12 : r_sq_top - 4, c_int - 10 : c_int + 11].mean() > 0.35
        has_dn = blk[r_sq_bot + 4 : r_sq_bot + 12, c_int - 10 : c_int + 11].mean() > 0.35

        def find_cap(r_start: int, direction: int) -> float | None:
            r_cur = r_start
            while rT + 15 < r_cur < rB - 15 and blk[r_cur, c_int - 10 : c_int + 11].any():
                r_cur += direction
            r_end = r_cur - direction
            # Look for horizontal cap (width 38..92 px) within 25 px of stem end
            lo_r = min(r_end, r_end - direction * 25)
            hi_r = max(r_end, r_end - direction * 25)
            cap_rs = [r for r in range(max(rT + 10, lo_r), min(rB - 10, hi_r + 1))
                      if 38 <= blk[r, c_int - 55 : c_int + 56].sum() <= 92]
            if cap_rs:
                return float(np.mean(cap_rs) + 0.5)
            return float(r_end + 0.5)

        bar_notes = []
        if has_up:
            top_px = find_cap(r_sq_top - 5, -1)
        else:
            top_px = float(r_sq_top + 0.5)
            bar_notes.append("upper")

        if has_dn:
            bot_px = find_cap(r_sq_bot + 5, +1)
        else:
            bot_px = float(r_sq_bot + 0.5)
            bar_notes.append("lower")

        x_val = float(x_axis.to_value(c_c))
        y_val, y_lo, y_hi = value_interval(y_axis, r_c, center_err_px)
        bar_hi = float(y_axis.to_value(top_px))
        bar_lo = float(y_axis.to_value(bot_px))
        bar_hi_outer = float(y_axis.to_value(top_px - cap_err_px))
        bar_lo_outer = float(y_axis.to_value(bot_px + cap_err_px))

        b_note = (f"error bar ({' and '.join(bar_notes)} end) not protruding outside marker: bounded by marker vertical extent"
                  if bar_notes else f"error bar caps ({top_px:.1f}, {bot_px:.1f}) px")
        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "bar_lo": min(bar_lo, bar_hi),
            "bar_hi": max(bar_lo, bar_hi),
            "bar_lo_outer": min(bar_lo_outer, bar_hi_outer),
            "bar_hi_outer": max(bar_lo_outer, bar_hi_outer),
            "px": [c_c, r_c],
            "px_lo": top_px,
            "px_hi": bot_px,
            "note": f"espin {m_note}; {b_note}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })

    readings.sort(key=lambda d: d["x"])
    return readings
