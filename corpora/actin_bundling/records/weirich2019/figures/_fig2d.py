"""Shared calibrated reader for Main Fig. 2D (initial myosin centering speed vs distance, with error bars)."""

from __future__ import annotations

import numpy as np
from _calib import (
    detect_template_mp,
    find_blobs,
    find_box_spines,
    find_ticks,
    make_point_reading,
    match_ticks_regular,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # 6 interior x-ticks at 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 um; left = 0.0, right = 3.5.
    xt = match_ticks_regular(find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5), 6)
    x_marks = [(left[0], 0.0)] + [(xt[i], 0.5 * (i + 1)) for i in range(6)] + [(right[0], 3.5)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # IMPORTANT: y=-2 tick sits ~21 px ABOVE bottom spine (bottom spine is at y ~ -2.20);
    # top spine is exact y=10, and the 6 interior left/right ticks are 8, 6, 4, 2, 0, -2.
    yt = match_ticks_regular(find_ticks(dark, left, (top[2] + 10, bot[1] - 5), "y", side=+1, band=(5, 25), frac=0.5), 6)
    y_marks = [(top[0], 10.0)] + [(yt[i], 8.0 - 2.0 * i) for i in range(6)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def _trace_vertical_caps(
    mask: np.ndarray,
    col_c: float,
    r_center: float,
    r_clear_up: int,
    r_clear_dn: int,
    r_min_b: int,
    r_max_b: int,
    y_axis,
    cap_half_w: int = 18,
    gap: int = 6,
):
    """Trace vertical error-bar stem from center upward and downward, refining horizontal cap centers."""
    H, W = mask.shape
    ci = int(round(col_c))
    c0, c1 = max(0, ci - 4), min(W, ci + 5)

    # Walk upward from r_clear_up
    r_up = r_clear_up
    missed = 0
    for r in range(r_clear_up, max(r_min_b - 1, -1), -1):
        if mask[r, c0:c1].any():
            r_up = r
            missed = 0
        else:
            missed += 1
            if missed > gap:
                break

    # Walk downward from r_clear_dn
    r_dn = r_clear_dn
    missed = 0
    for r in range(r_clear_dn, min(r_max_b + 1, H)):
        if mask[r, c0:c1].any():
            r_dn = r
            missed = 0
        else:
            missed += 1
            if missed > gap:
                break

    def refine_cap(r_end: int, direction: int) -> float:
        cw0, cw1 = max(0, ci - cap_half_w), min(W, ci + cap_half_w + 1)
        widths = []
        for dr in range(0, 12):
            rr = r_end + direction * dr
            if 0 <= rr < H:
                widths.append((rr, int(mask[rr, cw0:cw1].sum())))
        if not widths:
            return float(r_end) + 0.5
        best_w = max(w for _, w in widths)
        if best_w >= 14:
            cap_rs = [rr for rr, w in widths if w >= 0.75 * best_w]
            return float(np.mean(cap_rs)) + 0.5
        return float(r_end) + 0.5

    row_hi = refine_cap(r_up, +1)
    row_lo = refine_cap(r_dn, -1)

    end_unc = 2.5 + y_axis.residual_px
    y_hi = float(y_axis.to_value(row_hi))
    y_lo = float(y_axis.to_value(row_lo))
    y_hi_out = float(y_axis.to_value(row_hi - end_unc))
    y_lo_out = float(y_axis.to_value(row_lo + end_unc))
    return (min(y_lo, y_hi), max(y_lo, y_hi), min(y_lo_out, y_hi_out), max(y_lo_out, y_hi_out), row_hi, row_lo)


def read_series(image: np.ndarray, series: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 8 : bot[1] - 8, left[2] + 8 : right[1] - 8] = True
    # Exclude legend box (top-left)
    roi[: top[2] + 490, : left[2] + 770] = False

    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
    out = []

    if series == "LD3p6":
        # 2 black 'x' markers with vertical error bars (exclude axis ticks near spines)
        inner = roi.copy()
        inner[: top[2] + 35, :] = False
        inner[bot[1] - 35 :, :] = False
        inner[:, : left[2] + 35] = False
        inner[:, right[1] - 35 :] = False
        black = (R < 60) & (G < 60) & (B < 60) & inner
        blobs = sorted([b for b in find_blobs(black, min_area=1500) if 55 <= b["w"] <= 90], key=lambda b: b["c"])
        if len(blobs) != 2:
            raise RuntimeError(f"Fig. 2D LD=3.6: expected 2 black x's, found {len(blobs)}")
        for b in blobs:
            # Find the x marker center from the wide diagonal-stroke rows (width > 38 px, wider than 30-px cap)
            sub = black[b["r0"] : b["r1"] + 1, b["c0"] : b["c1"] + 1]
            # Span of nonzero pixels in each row
            x_rows = []
            for rr in range(sub.shape[0]):
                cs = np.where(sub[rr])[0]
                if len(cs) >= 2 and (cs[-1] - cs[0] + 1) >= 42:
                    x_rows.append(b["r0"] + rr)
            r_mid = float(np.mean(x_rows)) + 0.5
            c_mid = 0.5 * (b["c0"] + b["c1"]) + 0.5
            # Caps sit cleanly at top and bottom of the blob
            eb = _trace_vertical_caps(black, c_mid, r_mid, int(r_mid - 20), int(r_mid + 20), top[2] + 5, bot[1] - 5, y_axis)
            out.append(
                make_point_reading(
                    c_mid,
                    r_mid,
                    x_axis,
                    y_axis,
                    half_px_y=3.5,
                    note=f"black x, caps=({eb[4]:.1f},{eb[5]:.1f}), {cal_note}",
                    error_bars=eb[:4],
                )
            )
        return out

    elif series == "LD5p1":
        # 4 green filled diamonds (83x83 px rotated square |dr|+|dc| <= 41) with vertical error bars
        green = (G > 70) & ((G - R) > 35) & ((G - B) > 25) & roi
        rad = 41
        yy, xx = np.ogrid[-rad : rad + 1, -rad : rad + 1]
        tpl = ((np.abs(xx) + np.abs(yy)) <= rad).astype(float)
        pts = detect_template_mp(green, tpl, min_full_frac=0.75, min_res_frac=0.60, min_dist=25, refine_r=6, max_hits=10)
        pts.sort(key=lambda p: p[1])
        if len(pts) != 4:
            raise RuntimeError(f"Fig. 2D LD=5.1: expected 4 green diamonds, found {len(pts)}")
        for r_c, c_c in pts:
            # Trace stem starting just outside diamond tips (r_c +- 38); gap=22 bridges blue-circle crossing
            eb = _trace_vertical_caps(green, c_c, r_c, int(r_c - 38), int(r_c + 38), top[2] + 5, bot[1] - 5, y_axis, gap=22)
            out.append(
                make_point_reading(
                    c_c,
                    r_c,
                    x_axis,
                    y_axis,
                    half_px_y=3.5,
                    note=f"green diamond, caps=({eb[4]:.1f},{eb[5]:.1f}), {cal_note}",
                    error_bars=eb[:4],
                )
            )
        return out

    elif series == "LD4p4":
        # 8 blue open circles (outer R=28.5, inner R=17.5) with vertical error bars
        blue = (B > 110) & ((B - R) > 40) & ((B - G) > 35) & roi
        rad_i = 29
        yy, xx = np.ogrid[-rad_i : rad_i + 1, -rad_i : rad_i + 1]
        dist = np.sqrt(xx * xx + yy * yy)
        tpl = ((dist >= 17.0) & (dist <= 28.5)).astype(float)
        # Note: one circle (d ~ 1.53) is ~30% overlapped by a green diamond, so min_full_frac=0.55 captures it cleanly
        pts = detect_template_mp(blue, tpl, min_full_frac=0.55, min_res_frac=0.45, min_dist=25, refine_r=6, max_hits=15)
        pts.sort(key=lambda p: p[1])
        if len(pts) != 8:
            raise RuntimeError(f"Fig. 2D LD=4.4: expected 8 blue circles, found {len(pts)}")
        for r_c, c_c in pts:
            ci = int(round(c_c))
            ri = int(round(r_c))
            # Check if a genuine vertical stem extends OUTSIDE the ring (width >= 4 px centered at ci)
            has_ext_up = bool((blue[max(0, ri - 37) : max(0, ri - 30), max(0, ci - 2) : min(W, ci + 3)].sum(axis=1) >= 4).any())
            has_ext_dn = bool((blue[min(H, ri + 31) : min(H, ri + 38), max(0, ci - 2) : min(W, ci + 3)].sum(axis=1) >= 4).any())
            if has_ext_up or has_ext_dn:
                # Trace from inside ring rim outward
                eb = _trace_vertical_caps(blue, c_c, r_c, ri - 27, ri + 27, top[2] + 5, bot[1] - 5, y_axis)
                row_hi, row_lo = eb[4], eb[5]
                up_len = r_c - row_hi
                dn_len = row_lo - r_c
                if abs(up_len - dn_len) > 5.0 and not has_ext_up:
                    row_hi = r_c - dn_len
                    end_unc = 3.0 + y_axis.residual_px
                    y_hi, y_lo = float(y_axis.to_value(row_hi)), float(y_axis.to_value(row_lo))
                    y_hi_o, y_lo_o = float(y_axis.to_value(row_hi - end_unc)), float(y_axis.to_value(row_lo + end_unc))
                    eb = (min(y_lo, y_hi), max(y_lo, y_hi), min(y_lo_o, y_hi_o), max(y_lo_o, y_hi_o), row_hi, row_lo)
                elif abs(up_len - dn_len) > 5.0 and not has_ext_dn:
                    row_lo = r_c + up_len
                    end_unc = 3.0 + y_axis.residual_px
                    y_hi, y_lo = float(y_axis.to_value(row_hi)), float(y_axis.to_value(row_lo))
                    y_hi_o, y_lo_o = float(y_axis.to_value(row_hi - end_unc)), float(y_axis.to_value(row_lo + end_unc))
                    eb = (min(y_lo, y_hi), max(y_lo, y_hi), min(y_lo_o, y_hi_o), max(y_lo_o, y_hi_o), row_hi, row_lo)
            else:
                # Error bar caps lie inside the open ring or flush with its 57-px wall.
                # Inside the white hole (|dr| <= 16), a horizontal cap fills dc=-10..+10 continuously.
                # Find the inner-facing edge of any continuous horizontal cap (strictly inside ring inner rim |dr|<=16).
                cap_dn = [ri + dr for dr in range(4, 17) if blue[ri + dr, ci - 10 : ci + 11].all()]
                cap_up = [ri - dr for dr in range(4, 17) if blue[ri - dr, ci - 10 : ci + 11].all()]
                # Distinguish a detached inner cap (has a white gap below/above it before the ring wall)
                # from a cap flush with the inner ring wall (first solid row starts at |dr| >= 13 and stays solid).
                dn_isolated = [rr for rr in cap_dn if rr <= ri + 13]
                up_isolated = [rr for rr in cap_up if rr >= ri - 13]
                if dn_isolated and not blue[min(H - 1, max(dn_isolated) + 2), ci - 6 : ci + 7].all():
                    # Clean isolated lower cap inside hole (e.g. Blue #5, dr=+6..+12)
                    row_lo = float(np.mean(dn_isolated))
                    row_hi = 2.0 * r_c - row_lo
                    end_unc = 3.5 + y_axis.residual_px
                elif up_isolated and not blue[max(0, min(up_isolated) - 2), ci - 6 : ci + 7].all():
                    row_hi = float(np.mean(up_isolated))
                    row_lo = 2.0 * r_c - row_hi
                    end_unc = 3.5 + y_axis.residual_px
                elif cap_dn:
                    # Cap merges into lower ring wall: inner edge is cap_dn[0], 7-px cap center is cap_dn[0] + 3.0
                    row_lo = float(cap_dn[0]) + 3.0
                    row_hi = 2.0 * r_c - row_lo
                    end_unc = 4.5 + y_axis.residual_px
                elif cap_up:
                    row_hi = float(cap_up[0]) - 3.0
                    row_lo = 2.0 * r_c - row_hi
                    end_unc = 4.5 + y_axis.residual_px
                else:
                    # Cap hidden entirely inside ring wall (center ~+-23 px, outer rim +-28.5 px).
                    row_hi, row_lo = r_c - 23.0, r_c + 23.0
                    end_unc = 5.5 + y_axis.residual_px
                y_hi, y_lo = float(y_axis.to_value(row_hi)), float(y_axis.to_value(row_lo))
                y_hi_o, y_lo_o = float(y_axis.to_value(row_hi - end_unc)), float(y_axis.to_value(row_lo + end_unc))
                eb = (min(y_lo, y_hi), max(y_lo, y_hi), min(y_lo_o, y_hi_o), max(y_lo_o, y_hi_o), row_hi, row_lo)

            out.append(
                make_point_reading(
                    c_c,
                    r_c,
                    x_axis,
                    y_axis,
                    half_px_y=3.5,
                    note=f"blue circle, caps=({eb[4]:.1f},{eb[5]:.1f}), {cal_note}",
                    error_bars=eb[:4],
                )
            )
        return out
    else:
        raise ValueError(f"Unknown series {series!r}")
