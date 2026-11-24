"""Calibrated readers for Fig. 9(a)-(b) (p. 9, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 9(a): autocatalytic coarsening rate lambda/k_0 vs V_bar/w^3
# --------------------------------------------------------------------------


def _calib_9a(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "bottom", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 2)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", [(xt_maj[0], 1e4), (xt_maj[1], 1e5)]),
        xt_all,
        [1e3, 1e4, 1e5],
    )

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.match_ticks_regular(C.filter_ticks_by_len(yl_all, 32, 48), 2)
    y_marks = [(yl_maj[0], 1e-3), (yl_maj[1], 1e-4), (box["r_bot"], 1e-5)]
    y_axis = C.verify_log_minor_ticks(C.fit_axis("log10", y_marks), yl_all, [1e-5, 1e-4, 1e-3])
    return box, x_axis, y_axis


def read_fig9a(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 9(a) series:
      - 'eq24_small': lower solid blue line (phi_0^B = phi_+^B, J = 0)
      - 'eq24_large': upper solid blue line (phi_0^B = 0, J = 4*pi*D_B*gamma*beta)
      - 'pde_small': orange circular discs (V_s/w^3 = 1.5e5)
      - 'pde_large': orange squares (V_s/w^3 = 8.8e5)
    """
    box, x_axis, y_axis = _calib_9a(img)
    if series in ("eq24_small", "eq24_large"):
        # Dark blue lines only (excludes pale-blue shaded band where R > 200)
        r = img[:, :, 0].astype(float)
        b = img[:, :, 2].astype(float)
        dark_blue = C.interior_mask((r < 75) & (b > 135) & (b - r > 70), box, inset=5)
        # Exclude legend swatch in upper-right (V > 7e3, lambda > 3.5e-4)
        r_leg = int(math.ceil(y_axis.to_pixel(3.5e-4)))
        c_leg = int(math.floor(x_axis.to_pixel(7e3)))
        dark_blue[:r_leg, c_leg:] = False
        # Strictly separate upper vs lower Eq. (24) line across the geometric midpoint
        # y_mid(V) = 1.25e-4 * (1e4 / V) so occlusion by orange squares can never alias lines.
        H, W = dark_blue.shape
        for col in range(box["left_in"], box["right_in"]):
            xv_col = x_axis.to_value(col + 0.5)
            r_mid = int(round(y_axis.to_pixel(1.25e-4 * (1e4 / xv_col))))
            if series == "eq24_large":
                dark_blue[max(0, r_mid) :, col] = False
            else:
                dark_blue[: min(H, r_mid), col] = False
        return C.read_curve_at(dark_blue, x_axis, y_axis, at, search=4, min_run=2, max_run=40)

    # Orange PDE markers (exclude legend orange circle at V > 1e4, lambda > 7e-4)
    om = C.interior_mask(C.mask_orange(img), box, inset=5)
    r_leg = int(math.ceil(y_axis.to_pixel(7e-4)))
    c_leg = int(math.floor(x_axis.to_pixel(1e4)))
    om[:r_leg, c_leg:] = False

    lab, n = ndimage.label(om)
    pts = []
    my = max(y_axis.residual_px, 0.5)
    err_px = 2.0 + my
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        area = len(rr)
        if area < 400 or area > 2500:
            continue
        r0, r1, c0, c1 = int(rr.min()), int(rr.max()), int(cc.min()), int(cc.max())
        # Distinguish squares (upper series) vs discs (lower series) using BOTH bounding-box
        # corner occupancy AND area/bbox fill ratio (squares have corners_orange >= 3 or fill > 0.90)
        corners = 0
        for cr, cc_ in [(r0 + 3, c0 + 3), (r0 + 3, c1 - 3), (r1 - 3, c0 + 3), (r1 - 3, c1 - 3)]:
            if 0 <= cr < om.shape[0] and 0 <= cc_ < om.shape[1] and om[cr, cc_]:
                corners += 1
        fill = area / float((r1 - r0 + 1) * (c1 - c0 + 1))
        # Also note: one square at (791, 481) is partially clipped on its right edge by the adjacent
        # square (826, 503), so check whether its centroid sits above the geometric midpoint line!
        cx = (c0 + c1) / 2.0 + 0.5
        cy = (r0 + r1) / 2.0 + 0.5
        xv = x_axis.to_value(cx)
        yv = y_axis.to_value(cy)
        # Geometric midpoint between lower Eq. 24 (~7e-5 at 1e4) and upper Eq. 24 (~2.1e-4 at 1e4):
        # squares lie above 1.25e-4 * (1e4 / V); discs lie cleanly below!
        y_mid = 1.25e-4 * (1e4 / xv)
        is_square = (corners >= 3) or (fill >= 0.90) or (yv > y_mid)
        if series == "pde_large" and not is_square:
            continue
        if series == "pde_small" and is_square:
            continue
        # If partially clipped (w < 36 or h < 36), widen reading uncertainty honestly
        extra = max(0.0, (39.0 - min(r1 - r0 + 1, c1 - c0 + 1)) / 2.0)
        ep = err_px + extra
        a, b = y_axis.to_value(cy - ep), y_axis.to_value(cy + ep)
        pts.append(
            {
                "x": float(xv),
                "y": float(yv),
                "lo": float(min(a, b)),
                "hi": float(max(a, b)),
                "px": [cx, cy],
                "px_lo": float(cy - ep),
                "px_hi": float(cy + ep),
                "note": f"{'square' if is_square else 'disc'} area={area} bbox={c1 - c0 + 1}x{r1 - r0 + 1}",
            }
        )
    pts.sort(key=lambda d: d["x"])
    return pts


# --------------------------------------------------------------------------
# Fig. 9(b): stability diagram Q/(w^3 k_0) vs phibar
# --------------------------------------------------------------------------


def _calib_9b(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = [p for p, _ in C.find_inward_ticks(img, box, "top", 32, 48)]
    xt_reg = C.match_ticks_regular([box["c_left"]] + xt_all, 3)
    x_axis = C.fit_axis("linear", [(xt_reg[0], 0.0), (xt_reg[1], 0.01), (xt_reg[2], 0.02)])

    yl_all = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "left", 32, 48), 32, 48))
    if len(yl_all) != 3:
        raise RuntimeError(f"Fig. 9b expected 3 left major y-ticks (300, 200, 100), got {yl_all}")
    y_axis = C.fit_axis("linear", [(yl_all[0], 300.0), (yl_all[1], 200.0), (yl_all[2], 100.0), (box["r_bot"], 0.0)])
    return box, x_axis, y_axis


def read_fig9b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 9(b) curves:
      - 'Qstab_line': red dashed horizontal line (Eq. 25, phi_0^B = 0)
      - 'multi_boundary': black curve separating single droplet (dark blue) from multiple droplets (light blue)
      - 'droplet_boundary': black curve separating no droplets (white) from single/multiple droplets (blue)
    """
    box, x_axis, y_axis = _calib_9b(img)
    r_top, r_bot = box["top_in"] + 4, box["bot_in"] - 4

    if series == "Qstab_line":
        red = C.interior_mask(C.mask_red(img), box, inset=5)
        rr, cc = np.nonzero(red)
        if len(rr) < 500:
            raise RuntimeError("Could not detect red dashed line in Fig. 9b")
        # Measure global row center & half-thickness across all red dashes
        # (and locally near each requested `at` if a dash sits there)
        global_r = float(np.median(rr) + 0.5)
        half_global = float((np.percentile(rr, 95) - np.percentile(rr, 5)) / 2.0) + 1.0
        out = []
        for xv in at:
            pc = x_axis.to_pixel(xv)
            col = int(round(pc))
            local = np.nonzero(red[:, max(box["left_in"], col - 50) : min(box["right_in"], col + 51)])[0]
            r_use = float(np.median(local) + 0.5) if len(local) >= 50 else global_r
            out.append(
                C.make_reading(xv, pc, r_use, half_global, 0.0, x_axis, y_axis, col_used=None, extra_px=0.5)
            )
        return out

    r = img[:, :, 0].astype(float)
    b = img[:, :, 2].astype(float)
    blk = C.interior_mask(C.mask_black(img, thresh=75), box, inset=3)
    dark_blue = (b - r > 60) & (r < 80) & (b > 120)
    any_blue = (b - r > 20) & (b > 120)

    if series == "multi_boundary":
        out = []
        for xv in at:
            pc = x_axis.to_pixel(xv)
            col = int(round(pc))
            db_rows = np.nonzero(dark_blue[r_top:r_bot, col])[0] + r_top
            if len(db_rows) < 20:
                out.append({"x": float(xv), "y": None, "note": f"no single-droplet region at phibar={xv:g}"})
                continue
            top_db = int(db_rows[0])
            # Black boundary curve sits directly above top_db (within 18 px above)
            band = blk[max(r_top, top_db - 20) : top_db + 3, col]
            b_idx = np.nonzero(band)[0]
            if len(b_idx) == 0:
                out.append({"x": float(xv), "y": None, "note": f"no black boundary above dark blue at {xv:g}"})
                continue
            s, e = C._runs_1d(b_idx)[-1]
            center = max(r_top, top_db - 20) + (s + e) / 2.0 + 0.5
            half = (e - s + 1) / 2.0
            # Estimate local slope from col +/- 10
            slopes = []
            for sgn in (-10, 10):
                c2 = col + sgn
                if box["left_in"] + 4 <= c2 < box["right_in"] - 4:
                    db2 = np.nonzero(dark_blue[r_top:r_bot, c2])[0] + r_top
                    if len(db2) >= 20:
                        slopes.append((c2, float(db2[0])))
            slope = (slopes[1][1] - slopes[0][1]) / (slopes[1][0] - slopes[0][0]) if len(slopes) == 2 else 0.0
            out.append(C.make_reading(xv, pc, center, half, slope, x_axis, y_axis, col_used=col, extra_px=0.5))
        return out

    elif series == "droplet_boundary":
        # Left black curve separating white ("no droplets") from blue ("single" or "multiple droplets").
        # Close any_blue vertically by 18 px to bridge the red dashed line cleanly, then take the single
        # largest connected blue region. Its bottom edge in column `col` sits directly above the black boundary.
        blue_full = np.zeros_like(any_blue, dtype=bool)
        blue_full[r_top:r_bot, box["left_in"] + 2 : box["right_in"] - 2] = any_blue[
            r_top:r_bot, box["left_in"] + 2 : box["right_in"] - 2
        ]
        blue_closed = ndimage.binary_closing(blue_full, structure=np.ones((21, 1), dtype=bool))
        lab, n = ndimage.label(blue_closed)
        best = max(range(1, n + 1), key=lambda i: int((lab == i).sum()))
        region = lab == best
        out = []
        for xv in at:
            pc = x_axis.to_pixel(xv)
            col = int(round(pc))
            b_rows = np.nonzero(region[r_top:r_bot, col])[0] + r_top
            if len(b_rows) < 5:
                out.append({"x": float(xv), "y": None, "note": f"no blue droplet region at phibar={xv:g}"})
                continue
            bot_b = int(b_rows[-1])
            # Black boundary curve sits directly below bot_b
            band = blk[max(r_top, bot_b - 4) : min(r_bot, bot_b + 85), col]
            b_idx = np.nonzero(band)[0]
            if len(b_idx) == 0:
                out.append({"x": float(xv), "y": None, "note": f"no black droplet boundary at {xv:g}"})
                continue
            s, e = C._runs_1d(b_idx)[0]
            center = max(r_top, bot_b - 4) + (s + e) / 2.0 + 0.5
            half = (e - s + 1) / 2.0
            slopes = []
            for sgn in (-6, 6):
                c2 = col + sgn
                if box["left_in"] + 4 <= c2 < box["right_in"] - 4:
                    b2 = np.nonzero(region[r_top:r_bot, c2])[0] + r_top
                    if len(b2) >= 5:
                        slopes.append((c2, float(b2[-1])))
            slope = (slopes[1][1] - slopes[0][1]) / (slopes[1][0] - slopes[0][0]) if len(slopes) == 2 else 0.0
            out.append(C.make_reading(xv, pc, center, half, slope, x_axis, y_axis, col_used=col, extra_px=1.0))
        return out
    raise ValueError(f"Unknown Fig. 9b series {series!r}")
