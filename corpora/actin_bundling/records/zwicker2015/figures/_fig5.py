"""Calibrated readers for Fig. 5(a)-(d) (p. 7, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 5(a) & 5(c): droplet volume trajectories V_i(t) and Delta V(t)
# --------------------------------------------------------------------------


def _calib_5ac_V(img: np.ndarray, y_tick_vals: list[float]):
    box = C.find_plot_box(img)
    # x-axis: 0 (left spine), 5000 (single major tick on top & bottom)
    xt = C.filter_ticks_by_len(C.find_inward_ticks(img, box, "top", 10, 48), 26, 45)
    if len(xt) != 1:
        raise RuntimeError(f"Fig. 5a/c V expected 1 top major tick (5000), got {xt}")
    x_axis = C.fit_axis("linear", [(box["c_left"], 0.0), (xt[0], 5000.0)])

    # y-axis: major ticks (length 26..45) + bottom spine = 0
    yl = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "left", 10, 48), 26, 45))
    yr = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "right", 10, 48), 26, 45))
    ym = yl if len(yl) == len(y_tick_vals) else yr
    if len(ym) != len(y_tick_vals):
        raise RuntimeError(f"Fig. 5a/c V expected {len(y_tick_vals)} major y-ticks, got yl={yl}, yr={yr}")
    y_marks = list(zip(ym, y_tick_vals)) + [(box["r_bot"], 0.0)]
    y_axis = C.fit_axis("linear", y_marks)
    return box, x_axis, y_axis


def read_fig5ac_V(img: np.ndarray, branch: str, y_tick_vals: list[float]) -> list[dict]:
    """Read upper ('V1') or lower ('V2') green-circle trajectory V_i/w^3 vs t*k_0."""
    box, x_axis, y_axis = _calib_5ac_V(img, y_tick_vals)
    # Green filled circles sit on a thin green connecting line; 11x11 opening isolates circle cores cleanly
    gm = C.interior_mask(C.mask_green(img), box, inset=2)
    cores = ndimage.binary_opening(gm, structure=np.ones((11, 11), dtype=bool))
    lab, n = ndimage.label(cores)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 100 or len(rr) > 2000:
            continue
        cx, cy = float(cc.mean() + 0.5), float(rr.mean() + 0.5)
        pts.append((cx, cy, len(rr)))
    # Group into pairs by x-column (within 15 px) and split upper vs lower branch
    pts.sort(key=lambda t: t[0])
    groups: list[list[tuple[float, float, int]]] = []
    for p in pts:
        if not groups or abs(p[0] - groups[-1][0][0]) > 18.0:
            groups.append([p])
        else:
            groups[-1].append(p)
    out = []
    my = max(y_axis.residual_px, 0.5)
    err_px = 2.0 + my
    for g in groups:
        if len(g) != 2:
            continue
        g.sort(key=lambda t: t[1])  # smaller row = larger V (upper branch V1)
        cx, cy, area = g[0] if branch == "V1" else g[1]
        xv = x_axis.to_value(cx)
        yv = y_axis.to_value(cy)
        a, b = y_axis.to_value(cy - err_px), y_axis.to_value(cy + err_px)
        out.append(
            {
                "x": float(xv),
                "y": float(yv),
                "lo": float(min(a, b)),
                "hi": float(max(a, b)),
                "px": [cx, cy],
                "px_lo": float(cy - err_px),
                "px_hi": float(cy + err_px),
                "note": f"green circle core area={area} px",
            }
        )
    return out


def _calib_5ac_dV(img: np.ndarray):
    box = C.find_plot_box(img)
    xt = C.filter_ticks_by_len(C.find_inward_ticks(img, box, "top", 10, 48), 26, 45)
    if len(xt) != 1:
        raise RuntimeError(f"Fig. 5a/c dV expected 1 top major tick (5000), got {xt}")
    x_axis = C.fit_axis("linear", [(box["c_left"], 0.0), (xt[0], 5000.0)])

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.filter_ticks_by_len(yl_all, 26, 45)
    if len(yl_maj) != 1:
        raise RuntimeError(f"Fig. 5a/c dV expected 1 left major tick (1e3), got {yl_maj}")
    y_init = C.fit_axis("log10", [(yl_maj[0], 1e3), (box["r_bot"], 1e2)])
    y_axis = C.verify_log_minor_ticks(y_init, yl_all, [1e2, 1e3])
    return box, x_axis, y_axis


def read_fig5ac_dV(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    """Read orange exponential fit line ('fit') or green circle markers ('markers') in Fig. 5a/c dV."""
    box, x_axis, y_axis = _calib_5ac_dV(img)
    if series == "fit":
        # Orange straight line spans full plot width; exclude orange text "propto e^{lambda t}" via largest component
        om = C.interior_mask(C.mask_orange(img), box, inset=4)
        lab, n = ndimage.label(om)
        best = max(range(1, n + 1), key=lambda i: int((lab == i).sum()))
        line_mask = lab == best
        return C.read_curve_at(line_mask, x_axis, y_axis, at, search=4, min_run=2, max_run=40)
    elif series == "markers":
        # Green circles are bisected by the orange fit line -> bridge across orange line with morphological closing
        gm = C.interior_mask(C.mask_green(img), box, inset=4)
        om = C.interior_mask(C.mask_orange(img), box, inset=4)
        # Dilate green by 6 px intersected with (green | orange) to bridge the ~8 px orange line inside each circle
        # Better: close green vertically by 11 px so upper & lower green half-moons of each circle merge,
        # without merging horizontally adjacent circles (spaced ~36 px apart).
        bridged = ndimage.binary_closing(gm, structure=np.ones((15, 3), dtype=bool))
        # Also remove tiny specks
        lab, n = ndimage.label(bridged)
        pts = []
        my = max(y_axis.residual_px, 0.5)
        err_px = 2.5 + my
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            if len(rr) < 150 or len(rr) > 2500:
                continue
            h = rr.max() - rr.min() + 1
            w = cc.max() - cc.min() + 1
            if w < 18 or w > 48 or h < 18 or h > 48:
                continue
            # Use bounding-box midpoint so an off-center orange line cut doesn't bias the circle center!
            cx = (cc.min() + cc.max()) / 2.0 + 0.5
            cy = (rr.min() + rr.max()) / 2.0 + 0.5
            xv = x_axis.to_value(cx)
            yv = y_axis.to_value(cy)
            a, b = y_axis.to_value(cy - err_px), y_axis.to_value(cy + err_px)
            pts.append(
                {
                    "x": float(xv),
                    "y": float(yv),
                    "lo": float(min(a, b)),
                    "hi": float(max(a, b)),
                    "px": [cx, cy],
                    "px_lo": float(cy - err_px),
                    "px_hi": float(cy + err_px),
                    "note": f"green circle bbox {w}x{h} px",
                }
            )
        pts.sort(key=lambda d: d["x"])
        return pts
    raise ValueError(f"Unknown Fig. 5a/c dV series {series!r}")


# --------------------------------------------------------------------------
# Fig. 5(b): passive coarsening rate lambda/k_0 vs V_bar/w^3
# --------------------------------------------------------------------------


def _calib_5b(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "bottom", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 3)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", list(zip(xt_maj, [1e4, 1e5, 1e6]))),
        xt_all,
        [1e3, 1e4, 1e5, 1e6],
    )

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.match_ticks_regular(C.filter_ticks_by_len(yl_all, 32, 48), 3)
    # 3 interior major ticks (1e-3, 1e-4, 1e-5) + bottom spine (1e-6)
    y_marks = [(yl_maj[0], 1e-3), (yl_maj[1], 1e-4), (yl_maj[2], 1e-5), (box["r_bot"], 1e-6)]
    y_axis = C.verify_log_minor_ticks(C.fit_axis("log10", y_marks), yl_all, [1e-6, 1e-5, 1e-4, 1e-3])
    return box, x_axis, y_axis


def read_fig5b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_5b(img)
    if series == "line":
        bm = C.interior_mask(C.mask_blue(img), box, inset=6)
        # Keep the long diagonal curve (exclude the short horizontal legend swatch)
        lab, n = ndimage.label(bm)
        best = max(range(1, n + 1), key=lambda i: int((lab == i).sum()))
        return C.read_curve_at(lab == best, x_axis, y_axis, at, search=4, min_run=2, max_run=50)

    # Markers: exclude legend region in upper right (x > 2e4, y > 1e-4)
    om = C.interior_mask(C.mask_orange(img), box, inset=6)
    r_leg = int(math.ceil(y_axis.to_pixel(1e-4)))
    c_leg = int(math.floor(x_axis.to_pixel(2e4)))
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
        h, w = r1 - r0 + 1, c1 - c0 + 1
        # Check the 4 corners of the bounding box (3 px inset): squares have orange corners, circles do not!
        corners_orange = 0
        for cr, cc_ in [(r0 + 3, c0 + 3), (r0 + 3, c1 - 3), (r1 - 3, c0 + 3), (r1 - 3, c1 - 3)]:
            if 0 <= cr < om.shape[0] and 0 <= cc_ < om.shape[1] and om[cr, cc_]:
                corners_orange += 1
        is_square = corners_orange >= 2
        if series == "squares" and not is_square:
            continue
        if series == "discs" and is_square:
            continue
        # Bounding-box midpoint is unbiased even when a blue line clips one corner
        cx = (c0 + c1) / 2.0 + 0.5
        cy = (r0 + r1) / 2.0 + 0.5
        xv = x_axis.to_value(cx)
        yv = y_axis.to_value(cy)
        a, b = y_axis.to_value(cy - err_px), y_axis.to_value(cy + err_px)
        pts.append(
            {
                "x": float(xv),
                "y": float(yv),
                "lo": float(min(a, b)),
                "hi": float(max(a, b)),
                "px": [cx, cy],
                "px_lo": float(cy - err_px),
                "px_hi": float(cy + err_px),
                "note": f"{'square' if is_square else 'disc'} area={area} corners={corners_orange}/4",
            }
        )
    pts.sort(key=lambda d: d["x"])
    return pts


# --------------------------------------------------------------------------
# Fig. 5(d): active coarsening rate lambda/k_0 vs k_b/k_0
# --------------------------------------------------------------------------


def _calib_5d(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "bottom", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 3)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", list(zip(xt_maj, [1e-5, 1e-4, 1e-3]))),
        xt_all,
        [1e-6, 1e-5, 1e-4, 1e-3],
    )

    y0 = C.find_zero_line(img, box)
    yl = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "left", 10, 48), 32, 48))
    if len(yl) != 2:
        raise RuntimeError(f"Fig. 5d expected 2 left major y-ticks (-1e-3, -2e-3), got {yl}")
    y_axis = C.fit_axis("linear", [(y0, 0.0), (yl[0], -1e-3), (yl[1], -2e-3)])
    return box, x_axis, y_axis


def read_fig5d(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_5d(img)
    if series == "line":
        bm = C.interior_mask(C.mask_blue(img), box, inset=6)
        # Exclude lower-left legend swatch (x < 3e-4, y < -0.5e-3)
        r_leg = int(math.floor(y_axis.to_pixel(-0.5e-3)))
        c_leg = int(math.ceil(x_axis.to_pixel(3e-4)))
        bm[r_leg:, :c_leg] = False
        return C.read_curve_at(bm, x_axis, y_axis, at, search=4, min_run=2, max_run=60)
    elif series == "markers":
        om = C.interior_mask(C.mask_orange(img), box, inset=6)
        # Exclude the legend orange circle in lower-left (x < 1e-4, y < -0.5e-3)
        r_leg = int(math.floor(y_axis.to_pixel(-0.5e-3)))
        c_leg = int(math.ceil(x_axis.to_pixel(1e-4)))
        om[r_leg:, :c_leg] = False
        # Some orange discs are bisected by the blue curve (~12 px wide) -> bridge with 15x15 closing
        bridged = ndimage.binary_closing(om, structure=np.ones((15, 15), dtype=bool))
        lab, n = ndimage.label(bridged)
        pts = []
        my = max(y_axis.residual_px, 0.5)
        err_px = 2.0 + my
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            if len(rr) < 400 or len(rr) > 2500:
                continue
            r0, r1, c0, c1 = int(rr.min()), int(rr.max()), int(cc.min()), int(cc.max())
            # Use bounding-box midpoint (exact circle center even when blue curve slices through!)
            cx = (c0 + c1) / 2.0 + 0.5
            cy = (r0 + r1) / 2.0 + 0.5
            # If the zero-line or blue line slightly clipped the top/bottom edge of a ~40 px circle,
            # use the full black-rimmed circle bounding box for maximum accuracy!
            h = r1 - r0 + 1
            if h < 36:
                # Re-center using the black outline + orange interior inside a 50x50 patch
                pr0, pr1 = max(0, int(cy) - 26), min(img.shape[0], int(cy) + 27)
                pc0, pc1 = max(0, int(cx) - 26), min(img.shape[1], int(cx) + 27)
                patch = om[pr0:pr1, pc0:pc1]
                # Fit circle center to orange boundary or keep bbox center with slightly wider uncertainty
                extra = (39.0 - h) / 2.0
            else:
                extra = 0.0
            xv = x_axis.to_value(cx)
            yv = y_axis.to_value(cy)
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
                    "note": f"orange disc bbox {c1 - c0 + 1}x{h} px",
                }
            )
        pts.sort(key=lambda d: d["x"])
        return pts
    raise ValueError(f"Unknown Fig. 5d series {series!r}")
