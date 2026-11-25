"""Calibrated readers for main-text Figs. 1A, 2A, 2B, 3B, and 3C."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


# ---------------------------------------------------------------------------
# Fig. 1A (p. 2): State diagram in (delta/kappa, |chi_A|)
# ---------------------------------------------------------------------------

def _calibrate_fig1a(image: np.ndarray):
    """Calibrate Fig. 1A frame and ticks (semi-transparent shaded regions blend borders)."""
    _, gray, ht, hb, vl, vr = C.find_box(image, thresh=170.0, chroma_max=60.0)
    # Use relative local darkening in a thin 10-px strip just inside the left and bottom borders
    # to detect ticks cleanly across cyan, pink, and gray shaded regions.
    l_strip = gray[ht[2] + 6 : hb[1] - 6, vl[2] + 1 : vl[2] + 11]
    l_bg = ndimage.median_filter(l_strip, size=(41, 1))
    yt_cands = [
        c[0] for c in C.clusters(np.where((l_bg - l_strip > 25.0).mean(axis=1) > 0.5)[0] + ht[2] + 6)
        if (c[2] - c[1] + 1) <= 18
    ]
    yt = C.match_ticks_regular(yt_cands, 4, tol=4.0)

    b_strip = gray[hb[1] - 11 : hb[1] - 1, vl[2] + 6 : vr[1] - 6]
    b_bg = ndimage.median_filter(b_strip, size=(1, 41))
    xt_cands = [
        c[0] for c in C.clusters(np.where((b_bg - b_strip > 25.0).mean(axis=0) > 0.5)[0] + vl[2] + 6)
        if (c[2] - c[1] + 1) <= 18
    ]
    # Note: ticks at 0, 1, 2 are cleanly inside the box; tick 3 sits right near the right border.
    # Match the 4 regular ticks (0, 1, 2, 3) with slightly wider tolerance for border clipping at 3.
    xt = C.match_ticks_regular(xt_cands, 4, tol=6.0)
    x_axis = C.verify_axis("linear", list(zip(xt, [0.0, 1.0, 2.0, 3.0])), max_residual_px=2.5)
    y_axis = C.verify_axis("linear", list(zip(yt, [0.06, 0.05, 0.04, 0.03])), max_residual_px=1.5)
    return ht, hb, vl, vr, x_axis, y_axis


def read_fig1a(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read Fig. 1A boundaries: 'diff_solid', 'diff_dashed', 'osc', or 'dc_bdy'."""
    ht, hb, vl, vr, x_axis, y_axis = _calibrate_fig1a(image)
    arr = image.astype(float)
    r_min, r_max = ht[2] + 6, hb[1] - 6

    if curve == "diff_solid":
        # Dark navy blue solid boundary RGB ~ (0, 19, 183) (distinct from cyan bg and bright blue markers B~236)
        navy = (arr[:, :, 0] < 70) & (arr[:, :, 1] < 70) & (arr[:, :, 2] > 100) & (arr[:, :, 2] < 218)
        navy = C.mask_interior(navy, ht, hb, vl, vr, inset=6)
        out = []
        for xv in at:
            # Restrict to lower half (|chi_A| in [0.033, 0.047]) to avoid any dark pixels near top
            r_lo = int(round(y_axis.to_pixel(0.047)))
            r_hi = int(round(y_axis.to_pixel(0.032)))
            rec = C.sample_curve_y(
                navy, x_axis, y_axis, xv, r_lo, r_hi,
                pick_fn=lambda rs: max(rs, key=lambda ab: ab[1]),
                search_dc=8,
                extra_px=2.5,
                note_prefix="Fig. 1A blue solid diffusive boundary: ",
            )
            out.append(rec)
        return out

    if curve == "diff_dashed":
        # Dark navy dashed curve inside pink region (delta/kappa > sqrt(2), |chi_A| > 0.045)
        navy = (arr[:, :, 0] < 90) & (arr[:, :, 1] < 70) & (arr[:, :, 2] > 100)
        navy = C.mask_interior(navy, ht, hb, vl, vr, inset=6)
        r_lo = ht[2] + 6
        r_hi = int(round(y_axis.to_pixel(0.044)))
        c_start = int(round(x_axis.to_pixel(1.43)))
        # Collect all navy dash column centers in the pink region for clean interpolation across dash gaps
        xs_px, ys_px, hs_px = [], [], []
        for cc in range(c_start, vr[1] - 6):
            rs = C.runs_1d(navy[r_lo : r_hi + 1, cc], gap=3)
            if rs:
                best = max(rs, key=lambda ab: ab[1] - ab[0])
                xs_px.append(cc + 0.5)
                ys_px.append(r_lo + 0.5 * (best[0] + best[1]) + 0.5)
                hs_px.append(0.5 * (best[1] - best[0] + 1))
        xs_px = np.array(xs_px, dtype=float)
        ys_px = np.array(ys_px, dtype=float)
        hs_px = np.array(hs_px, dtype=float)
        out = []
        for xv in at:
            cf = float(x_axis.to_pixel(xv))
            near = np.abs(xs_px - cf) <= 45.0
            if near.sum() < 4:
                out.append({"x": float(xv), "y": None, "note": f"too few dashed pixels near x={xv:g}"})
                continue
            m, b = np.polyfit(xs_px[near], ys_px[near], 1)
            r_eval = float(m * cf + b)
            half_s = float(np.median(hs_px[near]))
            h_px = half_s + 1.0 + 2.5 + abs(m) * max(x_axis.residual_px, 0.5) + max(y_axis.residual_px, 0.5)
            yv = float(y_axis.to_value(r_eval))
            v1, v2 = float(y_axis.to_value(r_eval - h_px)), float(y_axis.to_value(r_eval + h_px))
            out.append({
                "x": float(xv),
                "y": yv,
                "lo": min(v1, v2),
                "hi": max(v1, v2),
                "px": [cf, r_eval],
                "px_lo": r_eval - h_px,
                "px_hi": r_eval + h_px,
                "note": f"Fig. 1A blue dashed boundary (local fit over {int(near.sum())} cols, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
            })
        return out

    if curve == "osc":
        # Pure red horizontal oscillatory boundary at |chi_A| = 0.045 (narrow band excludes pink-region triangles at >=0.046)
        red = (arr[:, :, 0] > 200) & (arr[:, :, 1] < 60) & (arr[:, :, 2] < 60)
        red = C.mask_interior(red, ht, hb, vl, vr, inset=6)
        r_lo = int(round(y_axis.to_pixel(0.0456)))
        r_hi = int(round(y_axis.to_pixel(0.0442)))
        out = []
        for xv in at:
            rec = C.sample_curve_y(
                red, x_axis, y_axis, xv, r_lo, r_hi,
                search_dc=4,
                note_prefix="Fig. 1A red oscillatory boundary: ",
            )
            out.append(rec)
        return out

    if curve == "dc_bdy":
        # Vertical black boundary delta_c/kappa as a function of |chi_A| (at gives |chi_A| values).
        # Note blue simulation squares at delta/kappa=1.4 slightly tint B near the line, while R,G stay <45 on the black line.
        blk = (arr[:, :, 0] < 45) & (arr[:, :, 1] < 45) & (arr[:, :, 2] < 225)
        blk = C.mask_interior(blk, ht, hb, vl, vr, inset=6)
        c_lo = int(round(x_axis.to_pixel(1.36)))
        c_hi = int(round(x_axis.to_pixel(1.48)))
        out = []
        for chi_val in at:
            rf = float(y_axis.to_pixel(chi_val))
            ri = int(round(rf))
            rs = C.runs_1d(blk[ri, c_lo : c_hi + 1], gap=2)
            if not rs:
                out.append({"x": float(chi_val), "y": None, "note": f"no black line at |chi_A|={chi_val:g}"})
                continue
            # Rightmost run in [1.36, 1.48] is the vertical black boundary (any blue square sits to its left)
            best = max(rs, key=lambda ab: ab[1])
            c_ctr = c_lo + 0.5 * (best[0] + best[1]) + 0.5
            half_w = 0.5 * (best[1] - best[0] + 1)
            h_px = half_w + 0.5 + max(x_axis.residual_px, 0.5)
            dc_val = float(x_axis.to_value(c_ctr))
            v1, v2 = float(x_axis.to_value(c_ctr - h_px)), float(x_axis.to_value(c_ctr + h_px))
            out.append({
                "x": float(chi_val),
                "y": dc_val,
                "lo": min(v1, v2),
                "hi": max(v1, v2),
                "px": [float(c_ctr), rf],
                "px_lo": rf - 8.0,
                "px_hi": rf + 8.0,
                "note": f"Fig. 1A black vertical boundary delta_c/kappa (half_w={half_w:.1f}px, cal_res={x_axis.residual_px:.2f}px)",
            })
        return out

    raise ValueError(f"Unknown Fig. 1A curve: {curve}")


# ---------------------------------------------------------------------------
# Fig. 2A and Fig. 2B (p. 4)
# ---------------------------------------------------------------------------

def _measure_tri_height(rmask: np.ndarray, x_axis, test_xs: list[float], r_min: int, r_max: int) -> float:
    """Measure the uncontaminated height H of Mathematica filled/empty upward triangle markers."""
    hs = []
    for xv in test_xs:
        cx = int(round(x_axis.to_pixel(xv)))
        for cc in range(cx - 5, cx + 1):
            if 0 <= cc < rmask.shape[1]:
                rs = C.runs_1d(rmask[r_min:r_max, cc], gap=2)
                if rs:
                    best = max(rs, key=lambda ab: ab[1] - ab[0])
                    h = best[1] - best[0] + 1
                    if 40 <= h <= 90:
                        hs.append(h)
    if not hs:
        return 66.0
    return float(np.median(hs))


def read_fig2a(image: np.ndarray, at: list[float], branch: str) -> list[dict]:
    """Read Fig. 2A first Fourier mode amplitudes: 'rhoA_static', 'rhoB_static', 'rhoA_travel', 'rhoB_travel'."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 4), [0.0, 1.0, 2.0, 3.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.2, 0.1, 0.0])))

    bmask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    r_min, r_max = ht[2] + 4, hb[1] - 4
    r_split_red = int(round(y_axis.to_pixel(0.085)))  # cleanly separates rho_B^t (>0.095) from rho_A^t (~0.071)

    if branch in ("rhoA_static", "rhoB_static"):
        # Blue text labels rho_B^s and rho_A^s sit below the curve / inside the loop;
        # use the solid/empty blue squares and curves:
        # At each x in [0, 1.4], rho_A^s is the LOWER curve (~0.075..0.124) and rho_B^s is the UPPER curve (~0.124..0.201).
        out = []
        for xv in at:
            cf = float(x_axis.to_pixel(xv))
            ci = int(round(cf))
            if branch == "rhoA_static":
                # Filled blue square centered on rho_A^s (y in [0.065, 0.130])
                r_lo = int(round(y_axis.to_pixel(0.132)))
                r_hi = int(round(y_axis.to_pixel(0.065)))
                # Exclude label "rho_A^s" which sits around x in [0.65, 1.05], y < 0.095 (curve is at y >= 0.098 there!)
                if 0.60 <= xv <= 1.05:
                    r_hi = int(round(y_axis.to_pixel(0.094)))
                rs = C.runs_1d(bmask[r_lo:r_hi, ci], gap=3)
                if not rs:
                    out.append({"x": float(xv), "y": None, "note": f"no blue rho_A^s pixels at x={xv:g}"})
                    continue
                # Thickest run is the filled blue square (H ~ 54 px)
                best = max(rs, key=lambda ab: ab[1] - ab[0])
                r_top = r_lo + best[0]
                r_bot = r_lo + best[1]
                r_ctr = 0.5 * (r_top + r_bot) + 0.5
                h_px = 4.0 + max(y_axis.residual_px, 0.5)
                out.append(C.make_marker_reading(
                    xv, r_ctr, cf, x_axis, y_axis, half_px=h_px,
                    note=f"Fig. 2A filled blue square rho_A^s (H={r_bot-r_top+1}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
                ))
            else:
                # Empty blue square + dashed curve on rho_B^s (y in [0.135, 0.212])
                r_lo = int(round(y_axis.to_pixel(0.212)))
                r_hi = int(round(y_axis.to_pixel(0.135)))
                # Collect all blue pixels in the upper empty square around ci (top rim to bottom rim!)
                sub = bmask[r_lo:r_hi, ci]
                rs = C.runs_1d(sub, gap=4)
                if not rs:
                    out.append({"x": float(xv), "y": None, "note": f"no blue rho_B^s pixels at x={xv:g}"})
                    continue
                # Exclude label "rho_B^s" below the upper curve (around x in [0.55, 0.95], y < 0.165)
                if 0.50 <= xv <= 1.05:
                    r_cut = int(round(y_axis.to_pixel(0.168))) - r_lo
                    rs = [ab for ab in rs if ab[0] < r_cut]
                if not rs:
                    out.append({"x": float(xv), "y": None, "note": f"no blue rho_B^s rim at x={xv:g}"})
                    continue
                # Outer top rim and outer bottom rim of the empty square (span ~ 50 px)
                r_top = r_lo + rs[0][0]
                r_bot = r_lo + rs[-1][1]
                if r_bot - r_top > 65:
                    # Just keep runs within 58 px of topmost rim
                    kept = [ab for ab in rs if ab[1] - rs[0][0] <= 60]
                    r_bot = r_lo + kept[-1][1]
                r_ctr = 0.5 * (r_top + r_bot) + 0.5
                h_px = 5.0 + max(y_axis.residual_px, 0.5)
                out.append(C.make_marker_reading(
                    xv, r_ctr, cf, x_axis, y_axis, half_px=h_px,
                    note=f"Fig. 2A empty blue square rho_B^s (span={r_bot-r_top+1}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
                ))
        return out

    if branch == "rhoA_travel":
        # Filled red triangles + horizontal solid red line (rows > r_split_red, above label rho_A^t)
        r_hi = int(round(y_axis.to_pixel(0.055)))
        H_tri = _measure_tri_height(rmask, x_axis, [1.6, 2.0, 2.5, 3.0], r_split_red, r_hi)
        out = []
        for xv in at:
            cf = float(x_axis.to_pixel(xv))
            ci = int(round(cf))
            # Find highest red apex in [ci-5, ci+1]
            best_top = None
            best_c = ci
            for cc in range(ci - 5, ci + 2):
                if 0 <= cc < rmask.shape[1]:
                    rs = C.runs_1d(rmask[r_split_red:r_hi, cc], gap=2)
                    for a, b in rs:
                        if (b - a + 1) >= 0.7 * H_tri:
                            rt = r_split_red + a
                            if best_top is None or rt < best_top:
                                best_top = rt
                                best_c = cc
            if best_top is None:
                out.append({"x": float(xv), "y": None, "note": f"no red rho_A^t triangle at x={xv:g}"})
                continue
            # Mathematica upward triangle anchor sits at r_apex + H_tri/3
            r_anc = best_top + H_tri / 3.0 + 0.5
            h_px = 5.0 + max(y_axis.residual_px, 0.5)
            out.append(C.make_marker_reading(
                xv, r_anc, cf, x_axis, y_axis, half_px=h_px,
                note=f"Fig. 2A filled red triangle rho_A^t (apex={best_top}, H_tri={H_tri:.0f}px, c={best_c})",
            ))
        return out

    if branch == "rhoB_travel":
        # Empty red triangles + dashed red curve (rows < r_split_red)
        # Exclude label "rho_B^t" sitting above curve around x in [1.9, 2.3], y > 0.135
        H_tri = 64.0
        out = []
        for xv in at:
            cf = float(x_axis.to_pixel(xv))
            ci = int(round(cf))
            r_lo = int(round(y_axis.to_pixel(0.185)))
            if 1.85 <= xv <= 2.35:
                r_lo = int(round(y_axis.to_pixel(0.135)))
            # Search [ci-5, ci+1] for empty triangle span (top apex to bottom base)
            best_top = None
            best_bot = None
            best_c = ci
            for cc in range(ci - 5, ci + 2):
                if 0 <= cc < rmask.shape[1]:
                    rs = C.runs_1d(rmask[r_lo:r_split_red, cc], gap=3)
                    if rs:
                        rt = r_lo + rs[0][0]
                        rb = r_lo + rs[-1][1]
                        span = rb - rt + 1
                        if 40 <= span <= 85:
                            if best_top is None or rt < best_top:
                                best_top, best_bot, best_c = rt, rb, cc
            if best_top is None:
                out.append({"x": float(xv), "y": None, "note": f"no empty red triangle rho_B^t at x={xv:g}"})
                continue
            span = float(best_bot - best_top + 1)
            r_anc = best_top + span / 3.0 + 0.5
            h_px = 6.0 + max(y_axis.residual_px, 0.5)
            out.append(C.make_marker_reading(
                xv, r_anc, cf, x_axis, y_axis, half_px=h_px,
                note=f"Fig. 2A empty red triangle rho_B^t (apex={best_top}, base={best_bot}, span={span:.0f}px, c={best_c})",
            ))
        return out

    raise ValueError(f"Unknown Fig. 2A branch: {branch}")


def read_fig2b(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read Fig. 2B traveling velocity: 'sim' (blue static squares + red traveling triangles)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 4), [0.0, 1.0, 2.0, 3.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.02, 0.01, 0.0])))

    bmask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)
    rmask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
    r_min, r_max = ht[2] + 4, hb[1] - 4
    y0_px = float(y_axis.to_pixel(0.0))

    # Calibrate blue square baseline center offset from static squares at x = 0.4, 0.8, 1.0 (offset +8px from ticks)
    sq_mids = []
    for xs in (0.4, 0.8, 1.0):
        cs = int(round(x_axis.to_pixel(xs))) + 8
        rs = np.nonzero(bmask[r_min:r_max, cs])[0] + r_min
        if len(rs):
            sq_mids.append(0.5 * (rs.min() + rs.max()) + 0.5)
    sq_offset = float(np.median(sq_mids) - y0_px) if sq_mids else 0.0

    # Measure triangle height H_tri directly from traveling triangles
    H_tri = _measure_tri_height(rmask, x_axis, [1.6, 2.0, 2.4, 2.8, 3.0], r_min, r_max)

    out = []
    for xv in at:
        cf = float(x_axis.to_pixel(xv))
        ci = int(round(cf))
        if xv <= 1.41:
            # Sample 8 px right of exact center (still well inside the ~50px square) to avoid inward black x-ticks at x=0, 1
            cs = min(vr[1] - 4, ci + 8)
            rs = np.nonzero(bmask[r_min:r_max, cs])[0] + r_min
            if len(rs) == 0:
                rs = np.nonzero(bmask[r_min:r_max, ci])[0] + r_min
            if len(rs) == 0:
                out.append({"x": float(xv), "y": None, "note": f"no blue square at x={xv:g}"})
                continue
            r_mid = 0.5 * (rs.min() + rs.max()) + 0.5 - sq_offset
            out.append(C.make_marker_reading(
                xv, r_mid, cf, x_axis, y_axis, half_px=4.0,
                note=f"Fig. 2B static blue square (offset={sq_offset:+.1f}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
            ))
        else:
            best_top = None
            best_c = ci
            for cc in range(ci - 5, ci + 2):
                if 0 <= cc < rmask.shape[1]:
                    rs = C.runs_1d(rmask[r_min:r_max, cc], gap=2)
                    for a, b in rs:
                        if (b - a + 1) >= 0.7 * H_tri:
                            rt = r_min + a
                            if best_top is None or rt < best_top:
                                best_top = rt
                                best_c = cc
            if best_top is None:
                out.append({"x": float(xv), "y": None, "note": f"no red triangle at x={xv:g}"})
                continue
            r_anc = best_top + H_tri / 3.0 + 0.5
            out.append(C.make_marker_reading(
                xv, r_anc, cf, x_axis, y_axis, half_px=6.0,
                note=f"Fig. 2B traveling red triangle (apex={best_top}, H_tri={H_tri:.0f}px, c={best_c})",
            ))
    return out


# ---------------------------------------------------------------------------
# Fig. 3B and Fig. 3C (p. 4)
# ---------------------------------------------------------------------------

def read_fig3b(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Fig. 3B analytical velocity ratio v_A^s / v_B^s (blue/green/red parabola)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    # 5 bottom ticks include delta_c/kappa between 1 and 2; select 4 evenly spaced [0, 1, 2, 3]
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 4), [0.0, 1.0, 2.0, 3.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 6), [8.0, 6.0, 4.0, 2.0, 0.0, -2.0])))

    bmask = C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60)
    gmask = C.rgb_mask(image, (0, 255, 0), tol=110, min_chroma=60)
    rmask = C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60)
    mask = C.mask_interior(bmask | gmask | rmask, ht, hb, vl, vr, inset=5)

    # Mask out the schematic arrows in upper-left (x < 2.25, y > 2.2) — curve stays below y=4.06 for x<2.25,
    # specifically arrows sit at y in [2.6, 6.6] for x < 2.20!
    r_arr_bot = int(round(y_axis.to_pixel(2.3)))
    c_arr_r = int(round(x_axis.to_pixel(2.15)))
    # Note at x = 2.0, true curve is at y = 3.0; at x = 2.15, true curve is at y = 3.62, while arrows end by x = 2.2 at y = 6!
    # Mask precisely each of the 3 arrow bands: blue/green/red arrows all have y > 2.5 and x < 1.85, plus red arrow y > 5.3 for x < 2.25.
    mask[: int(round(y_axis.to_pixel(2.5))), : int(round(x_axis.to_pixel(1.85)))] = False
    mask[: int(round(y_axis.to_pixel(5.2))), : int(round(x_axis.to_pixel(2.30)))] = False

    r_min, r_max = ht[2] + 5, hb[1] - 5
    out = []
    for xv in at:
        rec = C.sample_curve_y(
            mask, x_axis, y_axis, xv, r_min, r_max,
            pick_fn=lambda rs: max(rs, key=lambda ab: ab[1]),
            search_dc=5,
            note_prefix="Fig. 3B v_A^s/v_B^s: ",
        )
        out.append(rec)
    return out


def read_fig3c(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read Fig. 3C 'vratio' (red), 'ampratio' (blue), or 'psit_line' (vertical dashed psi^t)."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    # 4 bottom ticks: [0, 0.3, psi^t, 0.6]; select 3 evenly spaced [0, 0.3, 0.6]
    xt_reg = C.match_ticks_regular(xt, 3)
    x_axis = C.verify_axis("linear", list(zip(xt_reg, [0.0, 0.3, 0.6])))
    # 3 left ticks: [1.4, 1.0, 0.4] (non-uniform spacing!)
    yt_sorted = sorted(yt)
    if len(yt_sorted) != 3:
        raise RuntimeError(f"Expected 3 y ticks in Fig. 3C, got {yt_sorted}")
    y_axis = C.verify_axis("linear", list(zip(yt_sorted, [1.4, 1.0, 0.4])))

    if curve == "psit_line":
        # The 4th bottom inward tick (not in xt_reg) marks psi^t directly!
        extra = [c for c in xt if min(abs(c - r) for r in xt_reg) > 20.0]
        if not extra:
            raise RuntimeError(f"Could not locate psi^t tick among {xt}")
        c_psit = extra[0]
        psi_val = float(x_axis.to_value(c_psit))
        h_px = 2.0 + max(x_axis.residual_px, 0.5)
        v1, v2 = float(x_axis.to_value(c_psit - h_px)), float(x_axis.to_value(c_psit + h_px))
        r_mid = float(y_axis.to_pixel(0.7))
        return [
            {
                "x": float(xv),
                "y": psi_val,
                "lo": min(v1, v2),
                "hi": max(v1, v2),
                "px": [float(c_psit), r_mid],
                "px_lo": r_mid - 10.0,
                "px_hi": r_mid + 10.0,
                "note": f"Fig. 3C vertical dashed line psi^t tick (c={c_psit:.1f}, cal_res={x_axis.residual_px:.2f}px)",
            }
            for xv in at
        ]

    if curve == "vratio":
        mask = C.mask_interior(C.rgb_mask(image, (255, 0, 0), tol=90, min_chroma=60), ht, hb, vl, vr)
        # Mask out legend in upper-right (x > 0.35, y > 1.25)
        mask[: int(round(y_axis.to_pixel(1.26))), int(round(x_axis.to_pixel(0.35))) :] = False
        r_min, r_max = int(round(y_axis.to_pixel(1.30))), int(round(y_axis.to_pixel(0.95)))
        prefix = "Fig. 3C red v_A/v_B: "
    elif curve == "ampratio":
        mask = C.mask_interior(C.rgb_mask(image, (0, 0, 255), tol=90, min_chroma=60), ht, hb, vl, vr)
        # Mask out legend in upper-right (y > 0.60)
        r_min, r_max = int(round(y_axis.to_pixel(0.55))), hb[1] - 4
        prefix = "Fig. 3C blue rho_A/rho_B: "
    else:
        raise ValueError(f"Unknown Fig. 3C curve: {curve}")

    out = []
    for xv in at:
        rec = C.sample_curve_y(mask, x_axis, y_axis, xv, r_min, r_max, search_dc=4, note_prefix=prefix)
        out.append(rec)
    return out
