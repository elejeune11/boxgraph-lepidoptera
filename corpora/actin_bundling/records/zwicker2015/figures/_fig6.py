"""Calibrated readers for Fig. 6(a)-(b) (p. 7, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 6(a): N_max vs k_f/k_0, chi = 1 line, and k_f = k_f^c vertical line
# --------------------------------------------------------------------------


def _calib_6a(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "top", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 3)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", list(zip(xt_maj, [1e-5, 1e-4, 1e-3]))),
        xt_all,
        [1e-6, 1e-5, 1e-4, 1e-3],
    )

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.match_ticks_regular(C.filter_ticks_by_len(yl_all, 32, 48), 2)
    y_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", [(yl_maj[0], 10.0), (yl_maj[1], 1.0)]),
        yl_all,
        [0.1, 1.0, 10.0],
    )
    return box, x_axis, y_axis


def read_fig6a(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_6a(img)
    if series == "Nmax":
        # Dark blue step outline (excludes light blue fill under the curve where R > 150)
        r = img[:, :, 0].astype(float)
        g = img[:, :, 1].astype(float)
        b = img[:, :, 2].astype(float)
        dark_blue = (r < 70) & (b > 135) & (b - r > 70)
        mask = C.interior_mask(dark_blue, box, inset=5)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=3, min_run=2, max_run=60)

    elif series == "chi1":
        # Diagonal black line chi = 1 running from bottom-left (~4e-6, 0.7) to top (~2.3e-4, 48).
        # Extract the single largest black diagonal connected component (excludes text "chi=1" and slope triangle).
        blk = C.interior_mask(C.mask_black(img, thresh=65), box, inset=6)
        # Remove the vertical k_f^c line on the right (x > 5e-4)
        c_cut = int(math.floor(x_axis.to_pixel(5e-4)))
        blk[:, c_cut:] = False
        # Remove the slope triangle region (it is disconnected from the diagonal line!)
        lab, n = ndimage.label(blk)
        # The diagonal line chi=1 is the component with the largest vertical span (nearly full box height!)
        best_i = None
        best_span = 0
        for i in range(1, n + 1):
            rr = np.nonzero((lab == i).any(axis=1))[0]
            if len(rr) > 0:
                span = int(rr[-1] - rr[0])
                if span > best_span:
                    best_span = span
                    best_i = i
        if best_i is None or best_span < 400:
            raise RuntimeError("Could not isolate diagonal chi=1 black line in Fig. 6a")
        line_mask = lab == best_i
        return C.read_curve_at(line_mask, x_axis, y_axis, at, search=4, min_run=2, max_run=60)

    elif series == "kfc":
        # Vertical black line k_f = k_f^c on the right side (x ~ 1.2e-3).
        # Here `at` holds N_max heights (y-axis values) at which we read k_f^c / k_0 (returned as y!).
        blk = C.interior_mask(C.mask_black(img, thresh=65), box, inset=6)
        c0 = int(math.floor(x_axis.to_pixel(8e-4)))
        c1 = int(math.ceil(x_axis.to_pixel(2e-3)))
        out = []
        mx = max(x_axis.residual_px, 0.5)
        for nv in at:
            pr = y_axis.to_pixel(nv)
            row = int(round(pr))
            cols = np.nonzero(blk[row, c0:c1])[0] + c0
            if len(cols) == 0:
                out.append({"x": float(nv), "y": None, "note": f"no vertical kfc line at N={nv:g}"})
                continue
            runs = C._runs_1d(cols)
            # Pick widest run (or closest to median across whole vertical line)
            s, e = max(runs, key=lambda t: t[1] - t[0])
            cc = (s + e) / 2.0 + 0.5
            half = (e - s + 1) / 2.0
            h = half + 0.5 + mx
            kf_val = x_axis.to_value(cc)
            a, b = x_axis.to_value(cc - h), x_axis.to_value(cc + h)
            out.append(
                {
                    "x": float(nv),
                    "y": float(kf_val),
                    "lo": float(min(a, b)),
                    "hi": float(max(a, b)),
                    "px": [float(cc), float(pr)],
                    "note": f"vertical line width {e - s + 1} px",
                }
            )
        return out
    raise ValueError(f"Unknown Fig. 6a series {series!r}")


# --------------------------------------------------------------------------
# Fig. 6(b): Phase diagram k_b/k_0 vs k_f/k_0
# --------------------------------------------------------------------------


def _calib_6b(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "bottom", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 3)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", list(zip(xt_maj, [1e-3, 1e-2, 1e-1]))),
        xt_all,
        [1e-4, 1e-3, 1e-2, 1e-1],
    )

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.match_ticks_regular(C.filter_ticks_by_len(yl_all, 32, 48), 3)
    y_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", list(zip(yl_maj, [1e-1, 1e-2, 1e-3]))),
        yl_all,
        [1e-4, 1e-3, 1e-2, 1e-1],
    )
    return box, x_axis, y_axis


def read_fig6b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read boundaries of Fig. 6(b):
      - 'Bdroplet_upper': solid black upper boundary of B-droplet region (k_b = 16 k_0 (phi_0^B)^3)
      - 'Bdroplet_lower': solid black lower boundary of B-droplet region (k_b = 3 k_f)
      - 'spinodal_upper': dashed black upper boundary of homogeneous-unstable region (Eq. F2)
      - 'spinodal_lower': dashed black lower boundary of homogeneous-unstable region (Eq. F2)
    """
    box, x_axis, y_axis = _calib_6b(img)
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    blk = C.interior_mask(C.mask_black(img, thresh=85), box, inset=4)
    blue_fill = C.interior_mask((b - r > 15) & (b > 110), box, inset=4)
    orange_hatch = C.interior_mask((r - b > 30) & (r > 180), box, inset=4)

    r_top, r_bot = box["top_in"] + 5, box["bot_in"] - 5

    if series in ("Bdroplet_upper", "Bdroplet_lower"):
        # For each requested k_f, locate the blue B-droplet fill above the diagonal k_b = k_f,
        # then snap directly onto the solid black border line immediately adjacent to the blue fill!
        out = []
        for xv in at:
            pc = x_axis.to_pixel(xv)
            col = int(round(pc))
            r_diag = int(round(y_axis.to_pixel(xv)))
            # Collect blue pixels above diagonal in columns col-2..col+2
            brows = np.nonzero(blue_fill[r_top:r_diag, col - 2 : col + 3].any(axis=1))[0] + r_top
            if len(brows) < 5:
                out.append({"x": float(xv), "y": None, "note": f"no B-droplet blue region at kf={xv:g}"})
                continue
            if series == "Bdroplet_upper":
                edge_r = int(brows[0])
                # Find the black boundary run sitting right above/at edge_r (within 14 px above)
                band = blk[max(r_top, edge_r - 16) : edge_r + 3, col]
                b_idx = np.nonzero(band)[0]
                if len(b_idx) == 0:
                    center, half = float(edge_r), 3.0
                else:
                    s, e = C._runs_1d(b_idx)[-1]
                    center = max(r_top, edge_r - 16) + (s + e) / 2.0 + 0.5
                    half = (e - s + 1) / 2.0
            else:
                edge_r = int(brows[-1])
                band = blk[edge_r - 2 : min(r_bot, edge_r + 16), col]
                b_idx = np.nonzero(band)[0]
                if len(b_idx) == 0:
                    center, half = float(edge_r), 3.0
                else:
                    s, e = C._runs_1d(b_idx)[0]
                    center = (edge_r - 2) + (s + e) / 2.0 + 0.5
                    half = (e - s + 1) / 2.0
            # Estimate local slope from adjacent columns ±6 px
            slope = -1.0  # log-log slope ~ 1 (row decreases as col increases -> slope ~ -1)
            out.append(
                C.make_reading(
                    xv, pc, center, half, slope, x_axis, y_axis, col_used=col, extra_px=1.0
                )
            )
        return out

    elif series in ("spinodal_upper", "spinodal_lower"):
        # The dashed black spinodal boundary encloses the orange hatched region AND (for upper branch
        # at small k_f < 7e-3) extends slightly inside the blue B-droplet wedge above k_b = 3 k_f.
        # Method:
        #   - For k_f >= 1.2e-2 (upper) or ALL k_f >= 1.8e-3 (lower): trace the exact edge of the
        #     orange-hatched region (averaged across a small +/-6 px horizontal window to bridge hatch
        #     stripes), then snap to the nearby dashed black segment within +/-12 px if present!
        #   - For k_f < 1.2e-2 on the upper branch (inside blue wedge above 3*k_f): detect the dashed
        #     black segments lying strictly inside the blue B-droplet region between 3.1*k_f and 4.5*k_f
        #     (well below the "B-droplets" text!).
        out = []
        for xv in at:
            pc = x_axis.to_pixel(xv)
            col = int(round(pc))
            if series == "spinodal_upper" and xv < 8e-3:
                # Inside blue B-droplet region: find dashed black ink in corridor [3.15*kf, 4.3*kf]
                # (avoiding solid black 3*kf line below and "B-droplets" text above!)
                r_hi = int(math.floor(y_axis.to_pixel(4.3 * xv)))
                r_lo = int(math.ceil(y_axis.to_pixel(3.15 * xv)))
                # Search horizontally +/-20 px along the dashed line (~45 deg slope -> shift row with col!)
                pts_r = []
                for dc in range(-22, 23):
                    c_ = col + dc
                    if not (box["left_in"] + 4 <= c_ < box["right_in"] - 4):
                        continue
                    # Shift search window by dc along -1 slope (row moves by -dc)
                    ra = max(r_top, r_hi - dc)
                    rb = min(r_bot, r_lo - dc)
                    if rb <= ra:
                        continue
                    rows = np.nonzero(blk[ra:rb, c_])[0] + ra
                    for s, e in C._runs_1d(rows):
                        if 2 <= (e - s + 1) <= 14:
                            # Project back to column `col` along slope ~ -1
                            pts_r.append((c_, (s + e) / 2.0 + 0.5))
                if len(pts_r) < 3:
                    out.append({"x": float(xv), "y": None, "note": f"no upper spinodal dash near kf={xv:g}"})
                    continue
                cs = np.array([p[0] for p in pts_r])
                rs = np.array([p[1] for p in pts_r])
                # Linear fit through local dash pixels to interpolate smoothly across dash gaps!
                m, b0 = np.polyfit(cs, rs, 1)
                center = float(b0 + m * pc)
                resid = float(np.median(np.abs(rs - (b0 + m * cs))))
                out.append(
                    C.make_reading(
                        xv,
                        pc,
                        center,
                        max(3.5, resid + 2.0),
                        float(m),
                        x_axis,
                        y_axis,
                        col_used=None,
                        extra_px=1.5,
                        note=f"dashed spinodal inside blue wedge ({len(pts_r)} px)",
                    )
                )
            else:
                # Edge of orange-hatched region + snap to dashed black border!
                # First compute orange envelope across a +/-14 px window with local linear fit
                env_pts = []
                for dc in range(-18, 19):
                    c_ = col + dc
                    if not (box["left_in"] + 4 <= c_ < box["right_in"] - 4):
                        continue
                    orows = np.nonzero(orange_hatch[r_top:r_bot, c_])[0] + r_top
                    if len(orows) >= 3:
                        env_pts.append((c_, float(orows[0] if series == "spinodal_upper" else orows[-1])))
                if len(env_pts) < 5:
                    out.append({"x": float(xv), "y": None, "note": f"no orange spinodal hatch near kf={xv:g}"})
                    continue
                cs = np.array([p[0] for p in env_pts])
                rs = np.array([p[1] for p in env_pts])
                # Use outer envelope (min row for upper boundary, max row for lower boundary) via robust fit:
                # Because diagonal orange stripes hit the boundary periodically every ~12 px, take the
                # 15th percentile (upper) / 85th percentile (lower) relative to a preliminary slope fit.
                m0, b0 = np.polyfit(cs, rs, 1)
                res0 = rs - (b0 + m0 * cs)
                if series == "spinodal_upper":
                    keep = res0 <= np.percentile(res0, 35)
                else:
                    keep = res0 >= np.percentile(res0, 65)
                m_env, b_env = np.polyfit(cs[keep], rs[keep], 1)
                row_env = float(b_env + m_env * pc)

                # Now look for dashed black border pixels within +/-12 px of the orange envelope line!
                dash_pts = []
                for dc in range(-32, 33):
                    c_ = col + dc
                    if not (box["left_in"] + 4 <= c_ < box["right_in"] - 4):
                        continue
                    pred_r = int(round(b_env + m_env * c_))
                    # Dashed line lies just outside (above for upper, below for lower) the orange hatch
                    if series == "spinodal_upper":
                        ra, rb = max(r_top, pred_r - 15), min(r_bot, pred_r + 5)
                    else:
                        ra, rb = max(r_top, pred_r - 5), min(r_bot, pred_r + 15)
                    rows = np.nonzero(blk[ra:rb, c_])[0] + ra
                    for s, e in C._runs_1d(rows):
                        if 2 <= (e - s + 1) <= 14:
                            dash_pts.append((c_, (s + e) / 2.0 + 0.5))
                if len(dash_pts) >= 4:
                    dcs = np.array([p[0] for p in dash_pts])
                    drs = np.array([p[1] for p in dash_pts])
                    md, bd = np.polyfit(dcs, drs, 1)
                    center = float(bd + md * pc)
                    slope = float(md)
                    half = 3.5
                    note = f"dashed spinodal border ({len(dash_pts)} px)"
                else:
                    center = row_env
                    slope = float(m_env)
                    half = 5.0
                    note = "orange hatch boundary envelope"
                out.append(
                    C.make_reading(
                        xv,
                        pc,
                        center,
                        half,
                        slope,
                        x_axis,
                        y_axis,
                        col_used=None,
                        extra_px=1.5,
                        note=note,
                    )
                )
        return out
    raise ValueError(f"Unknown Fig. 6b series {series!r}")
