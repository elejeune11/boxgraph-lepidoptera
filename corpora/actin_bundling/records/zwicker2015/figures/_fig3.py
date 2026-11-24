"""Calibrated readers for Fig. 3(a)-(d) (p. 5, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def _calib_3ac(img: np.ndarray):
    box = C.find_plot_box(img)
    # x-axis: 0 (left spine), 10, 20, 30, 40 (top/bottom major ticks)
    xt = C.match_ticks_regular(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "top"), 30, 48), 4)
    x_marks = [(box["c_left"], 0.0), (xt[0], 10.0), (xt[1], 20.0), (xt[2], 30.0), (xt[3], 40.0)]
    x_axis = C.fit_axis("linear", x_marks)

    # y-axis: right ticks at 15e-3, 10e-3, 5e-3, -5e-3 + horizontal zero line y=0
    yr = C.filter_ticks_by_len(C.find_inward_ticks(img, box, "right"), 30, 48)
    y0 = C.find_zero_line(img, box)
    # Sort right ticks top-to-bottom: 15, 10, 5 above y0; -5 below y0
    above = sorted(r for r in yr if r < y0 - 40)
    below = sorted(r for r in yr if r > y0 + 40)
    if len(above) != 3 or len(below) != 1:
        raise RuntimeError(f"Fig. 3a/c unexpected right y-ticks: above={above}, below={below}, y0={y0}")
    y_marks = [
        (above[0], 15e-3),
        (above[1], 10e-3),
        (above[2], 5e-3),
        (y0, 0.0),
        (below[0], -5e-3),
    ]
    y_axis = C.fit_axis("linear", y_marks)
    return box, x_axis, y_axis


def _largest_component(mask: np.ndarray, min_width_frac: float = 0.25) -> np.ndarray:
    """Keep connected components spanning at least min_width_frac of image width (removes text)."""
    lab, n = ndimage.label(mask)
    W = mask.shape[1]
    out = np.zeros_like(mask, dtype=bool)
    for i in range(1, n + 1):
        comp = lab == i
        cols = np.nonzero(comp.any(axis=0))[0]
        if len(cols) > 0 and (cols[-1] - cols[0]) >= min_width_frac * W:
            out |= comp
    return out


def read_fig3ac(img: np.ndarray, at: list[float], curve: str) -> list[dict]:
    box, x_axis, y_axis = _calib_3ac(img)
    if curve == "passive":
        raw = C.interior_mask(C.mask_orange(img), box, inset=6)
        # Keep the single continuous orange curve (excludes orange text "no chemical reactions")
        mask = _largest_component(raw, min_width_frac=0.35)
        # In the steep left flank (R/w < 4), blue overlaps orange partially; use slightly wider search
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=90)
    elif curve == "active":
        raw = C.interior_mask(C.mask_blue(img), box, inset=6)
        mask = _largest_component(raw, min_width_frac=0.45)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=90)
    raise ValueError(f"Unknown Fig. 3a/c curve {curve!r}")


def _calib_3bd(img: np.ndarray):
    box = C.find_plot_box(img)
    # x-axis: log10, 4 major ticks 1e-6, 1e-5, 1e-4, 1e-3 on top & bottom
    xt_all = C.find_inward_ticks(img, box, "bottom", min_len=10, max_len=48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 4)
    x_marks = list(zip(xt_maj, [1e-6, 1e-5, 1e-4, 1e-3]))
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", x_marks), xt_all, [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    )

    # y-axis: log10, bottom spine = 1, major ticks at 10 and 100
    yl_all = C.find_inward_ticks(img, box, "left", min_len=10, max_len=48)
    yl_maj = sorted(C.filter_ticks_by_len(yl_all, 32, 48))
    if len(yl_maj) != 2:
        raise RuntimeError(f"Fig. 3b/d expected 2 left major ticks (100, 10), got {yl_maj}")
    y_marks = [(yl_maj[0], 100.0), (yl_maj[1], 10.0), (box["r_bot"], 1.0)]
    y_axis = C.verify_log_minor_ticks(C.fit_axis("log10", y_marks), yl_all, [1.0, 10.0, 100.0])
    return box, x_axis, y_axis


def read_fig3bd(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_3bd(img)
    blue = C.interior_mask(C.mask_blue(img), box, inset=6)
    if series == "Rbar":
        # Solid blue curve is the single huge connected component (> 50% width)
        mask = _largest_component(blue, min_width_frac=0.50)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=80)
    elif series == "Rcrit":
        # Dotted blue curve lies at R/w ~ 1.83 for k_f < 2e-4 and rises to ~2.45 near k_f = 1e-3.
        # Mask out the "R_crit" text label sitting above R/w = 2.4 on the left half.
        r_hi_left = int(math.floor(y_axis.to_pixel(2.35)))
        r_hi_right = int(math.floor(y_axis.to_pixel(2.65)))
        c_split = int(math.floor(x_axis.to_pixel(2e-4)))
        r_max_cut = int(math.ceil(y_axis.to_pixel(1.3)))
        c_max_cut = int(math.floor(x_axis.to_pixel(1.12e-3)))
        dot_mask = blue.copy()
        dot_mask[:r_hi_left, :c_split] = False
        dot_mask[:r_hi_right, c_split:] = False
        dot_mask[r_max_cut:, :] = False
        dot_mask[:, c_max_cut:] = False
        return C.read_dotted_curve_at(
            dot_mask, x_axis, y_axis, at, min_area=15, max_area=400, max_aspect=3.0, extra_px=1.5
        )
    elif series == "bifurcation":
        # Blue filled circle surrounded by black outline near k_f/k_0 ~ 1.3e-3, R/w ~ 2.75
        r0 = int(math.floor(y_axis.to_pixel(5.0)))
        r1 = int(math.ceil(y_axis.to_pixel(1.5)))
        c0 = int(math.floor(x_axis.to_pixel(6e-4)))
        c1 = int(math.ceil(x_axis.to_pixel(3e-3)))
        sub = np.zeros_like(blue, dtype=bool)
        sub[r0:r1, c0:c1] = blue[r0:r1, c0:c1]
        # The circle interior is cleanly separated from the solid curve by its black outline!
        pts = C.read_filled_markers(sub, x_axis, y_axis, min_area=600, max_area=6000, shape="circle")
        if len(pts) != 1:
            raise RuntimeError(f"Expected 1 bifurcation dot in Fig. 3b/d, found {len(pts)}")
        return pts
    raise ValueError(f"Unknown Fig. 3b/d series {series!r}")
