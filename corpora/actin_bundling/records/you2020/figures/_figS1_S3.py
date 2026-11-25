"""Calibrated readers for SI Figs. S1(a)-(b), S2, and S3(a)-(c)."""

from __future__ import annotations

import numpy as np
import _calib as C


def read_figS1(image: np.ndarray, at: list[float], panel: str, species: str) -> list[dict]:
    """Read Fig. S1(a) or S1(b) Fourier mode amplitudes |phi_hat_mu^j| for species A (red) or B (blue)."""
    y_vals = [0.14, 0.07, 0.0] if panel == "a" else [0.10, 0.05, 0.0]
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 3), [0.0, 5.0, 10.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), y_vals)))

    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    bmask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)

    # Calibrate Mathematica circular marker radius and baseline font-glyph vertical offset
    # directly from the unoccluded y=0 blue dots at j = 6, 7, 8, 9 (below legend rows).
    y0_px = float(y_axis.to_pixel(0.0))
    r_leg_cut = int(round(y_axis.to_pixel(0.025)))
    base_tops, base_bots = [], []
    for jb in (6, 7, 8, 9):
        cb = int(round(x_axis.to_pixel(jb)))
        rs = np.nonzero(bmask[r_leg_cut : hb[1] - 2, cb])[0] + r_leg_cut
        if len(rs):
            base_tops.append(int(rs.min()))
            base_bots.append(int(rs.max()))
    if not base_tops:
        raise RuntimeError("Could not locate baseline reference dots at j=6..9 in Fig. S1")
    r_top_base = float(np.median(base_tops))
    r_bot_base = float(np.median(base_bots))
    top_to_y0 = y0_px - r_top_base
    bot_to_y0 = r_bot_base - y0_px

    mask = rmask if species == "A" else bmask
    out = []
    for xv in at:
        j = int(round(xv))
        cj = int(round(x_axis.to_pixel(xv)))
        # Exclude legend region (around x in [4, 7], y in [0.04, 0.12])
        r_start_search = ht[2] + 4 if j <= 2 else r_leg_cut
        rs = np.nonzero(mask[r_start_search : hb[1] - 2, cj])[0] + r_start_search
        if len(rs) == 0:
            out.append({"x": float(xv), "y": None, "note": f"no {species} pixels at j={j}"})
            continue

        if species == "B":
            if j == 0:
                # Black inward x-tick sits at exact cj; sample 6 px right inside the 22-px-radius circle
                rs_off = np.nonzero(mask[r_leg_cut : hb[1] - 2, cj + 6])[0] + r_leg_cut
                # Correct circular chord vertical sag sqrt(R^2 - dx^2) - R ~ 0.8 px (< 1 px)
                r_data = float(rs_off.max()) - bot_to_y0 + 0.8 if len(rs_off) else float(rs.max()) - bot_to_y0
                h_px = 2.5 + max(y_axis.residual_px, 0.5)
                note = f"Fig. S1({panel}) B j=0 bot-edge off-tick zero-ref (cal_res={y_axis.residual_px:.2f}px)"
            elif j == 2:
                # Rising line on left side of j=2 leaves bottom edge completely clean (and no x-tick at j=2)
                r_data = float(rs.max()) - bot_to_y0
                h_px = 2.5 + max(y_axis.residual_px, 0.5)
                note = f"Fig. S1({panel}) B j=2 bot-edge zero-ref (cal_res={y_axis.residual_px:.2f}px)"
            else:
                # For j=1, 3, and j>=4 (including x-tick positions j=5, 10), top edge of dot is completely unoccluded!
                r_data = float(rs.min()) + top_to_y0
                h_px = 2.5 + max(y_axis.residual_px, 0.5)
                note = f"Fig. S1({panel}) B j={j} top-edge zero-ref (cal_res={y_axis.residual_px:.2f}px)"
        else:
            # Species A (red): drawn underneath blue
            if j == 1:
                # Peak j=1 is far below blue j=1 peak, so top edge of red dot is completely unoccluded
                r_data = float(rs.min()) + top_to_y0
                h_px = 2.0 + max(y_axis.residual_px, 0.5)
                note = f"Fig. S1({panel}) A j=1 top-edge zero-ref (cal_res={y_axis.residual_px:.2f}px)"
            else:
                # At j=3, red dot sits slightly below blue dot so its bottom crescent is unoccluded
                r_data = float(rs.max()) - bot_to_y0
                h_px = 3.0 + max(y_axis.residual_px, 0.5)
                note = f"Fig. S1({panel}) A j={j} bot-crescent zero-ref (cal_res={y_axis.residual_px:.2f}px)"

        y_val = float(y_axis.to_value(r_data))
        v1 = float(y_axis.to_value(r_data - h_px))
        v2 = float(y_axis.to_value(r_data + h_px))
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [float(cj), float(r_data)],
            "px_lo": float(r_data - h_px),
            "px_hi": float(r_data + h_px),
            "note": note,
        })
    return out


def read_figS2(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    """Read Fig. S2 Re[lambda_+^h], Re[lambda_-^h] (blue), or Im[lambda_+^h] (red)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 3), [0.0, 2.0, 4.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.02, 0.0, -0.02])))

    bmask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)

    # Legend occupies x in [0.4, 2.2], y in [-0.019, -0.008]; mask it out for all branches
    r_leg_top = int(round(y_axis.to_pixel(-0.008)))
    c_leg_r = int(round(x_axis.to_pixel(2.25)))
    bmask[r_leg_top:, :c_leg_r] = False
    rmask[r_leg_top:, :c_leg_r] = False

    r_min, r_max = ht[2] + 5, hb[1] - 5
    r_zero = int(round(y_axis.to_pixel(0.0)))

    out = []
    for xv in at:
        if branch == "re_plus":
            # Upper blue branch (topmost blue run)
            rec = C.sample_curve_y(
                bmask, x_axis, y_axis, xv, r_min, r_max,
                pick_fn=lambda rs: min(rs, key=lambda ab: ab[0]),
                note_prefix="Fig. S2 Re[lambda_+^h]: ",
            )
        elif branch == "re_minus":
            # Lower blue branch (bottommost blue run; coincides with re_plus when merged)
            rec = C.sample_curve_y(
                bmask, x_axis, y_axis, xv, r_min, r_max,
                pick_fn=lambda rs: max(rs, key=lambda ab: ab[1]),
                note_prefix="Fig. S2 Re[lambda_-^h]: ",
            )
        elif branch == "im_plus":
            # Upper red branch (y >= 0 -> row <= r_zero + 15)
            rec = C.sample_curve_y(
                rmask, x_axis, y_axis, xv, r_min, r_zero + 15,
                pick_fn=lambda rs: min(rs, key=lambda ab: ab[0]),
                note_prefix="Fig. S2 Im[lambda_+^h]: ",
            )
        else:
            raise ValueError(f"Unknown Fig. S2 branch: {branch}")
        out.append(rec)
    return out


def read_figS3(image: np.ndarray, at: list[float], panel: str, branch: str) -> list[dict]:
    """Read Fig. S3(a), (b), or (c) blue ('plus') or red ('minus') branch vs delta/delta_c."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 3), [0.0, 1.0, 2.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.04, 0.0, -0.04])))

    color = (0, 0, 255) if branch == "plus" else (255, 0, 0)
    mask = C.mask_interior(C.rgb_mask(image, color, tol=90, min_chroma=60), ht, hb, vl, vr)

    # Mask out legend box strictly (never overlaps curves):
    # S3a, S3b: legend is in upper-left (x < 0.9, y > 0.010)
    # S3c: legend is in upper-right (x > 1.1, y > 0.010)
    r_leg_bot = int(round(y_axis.to_pixel(0.010)))
    if panel in ("a", "b"):
        c_leg_r = int(round(x_axis.to_pixel(0.95)))
        mask[:r_leg_bot, :c_leg_r] = False
    else:
        c_leg_l = int(round(x_axis.to_pixel(1.05)))
        mask[:r_leg_bot, c_leg_l:] = False

    r_min, r_max = ht[2] + 5, hb[1] - 5
    out = []
    for xv in at:
        rec = C.sample_curve_y(
            mask, x_axis, y_axis, xv, r_min, r_max,
            search_dc=5,
            note_prefix=f"Fig. S3({panel}) {branch}: ",
        )
        out.append(rec)
    return out
