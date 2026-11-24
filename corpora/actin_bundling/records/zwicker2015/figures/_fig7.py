"""Calibrated readers for Fig. 7(a)-(b) (p. 8, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 7(a): single autocatalytic droplet growth rate dR/dt/(w*k_0) vs R/w
# --------------------------------------------------------------------------


def _calib_7a(img: np.ndarray):
    box = C.find_plot_box(img)
    # x-axis: 0 (left spine), 50, 100 (top/bottom major ticks)
    xt = C.match_ticks_regular(
        [box["c_left"]] + [p for p, _ in C.find_inward_ticks(img, box, "top", 32, 48)], 3
    )
    x_axis = C.fit_axis("linear", [(xt[0], 0.0), (xt[1], 50.0), (xt[2], 100.0)])

    # y-axis: 2e-3, 1e-3 (left/right major ticks) + black horizontal zero line y=0
    y0 = C.find_zero_line(img, box)
    yl = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "left", 32, 48), 32, 48))
    if len(yl) != 2:
        raise RuntimeError(f"Fig. 7a expected 2 left major y-ticks (2e-3, 1e-3), got {yl}")
    y_axis = C.fit_axis("linear", [(yl[0], 2e-3), (yl[1], 1e-3), (y0, 0.0)])
    return box, x_axis, y_axis


def _largest_wide_comp(mask: np.ndarray, min_width_frac: float = 0.25) -> np.ndarray:
    lab, n = ndimage.label(mask)
    W = mask.shape[1]
    out = np.zeros_like(mask, dtype=bool)
    for i in range(1, n + 1):
        comp = lab == i
        cols = np.nonzero(comp.any(axis=0))[0]
        if len(cols) > 0 and (cols[-1] - cols[0]) >= min_width_frac * W:
            out |= comp
    return out


def read_fig7a(img: np.ndarray, at: list[float], curve: str) -> list[dict]:
    box, x_axis, y_axis = _calib_7a(img)
    if curve == "inf":
        raw = C.interior_mask(C.mask_orange(img), box, inset=5)
        mask = _largest_wide_comp(raw, min_width_frac=0.18)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=90)
    elif curve == "fin":
        raw = C.interior_mask(C.mask_blue(img), box, inset=5)
        mask = _largest_wide_comp(raw, min_width_frac=0.18)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=100)
    elif curve == "core":
        raw = C.interior_mask(C.mask_green(img), box, inset=5)
        mask = _largest_wide_comp(raw, min_width_frac=0.18)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=100)
    raise ValueError(f"Unknown Fig. 7a curve {curve!r}")


# --------------------------------------------------------------------------
# Fig. 7(b): bifurcation diagram R/w vs phibar
# --------------------------------------------------------------------------


def _calib_7b(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = C.find_inward_ticks(img, box, "bottom", 10, 48)
    xt_maj = C.match_ticks_regular(C.filter_ticks_by_len(xt_all, 32, 48), 2)
    x_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", [(xt_maj[0], 1e-2), (xt_maj[1], 1e-1)]),
        xt_all,
        [1e-3, 1e-2, 1e-1],
    )

    yl_all = C.find_inward_ticks(img, box, "left", 10, 48)
    yl_maj = C.match_ticks_regular(C.filter_ticks_by_len(yl_all, 32, 48), 2)
    y_axis = C.verify_log_minor_ticks(
        C.fit_axis("log10", [(yl_maj[0], 1e2), (yl_maj[1], 1e1)]),
        yl_all,
        [1.0, 10.0, 100.0],
    )
    return box, x_axis, y_axis


def read_fig7b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_7b(img)
    if series == "core_Rbar":
        raw = C.interior_mask(C.mask_green(img), box, inset=5)
        mask = _largest_wide_comp(raw, min_width_frac=0.35)
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=100)
    elif series == "fin_Rbar":
        # Solid blue upper branch merges behind the green curve for phibar > 0.035;
        # read cleanly where visible (phibar in [0.0145, 0.035]), or take blue|green above 0.035!
        # Wait: let's keep purely blue curve pixels (0.0145..0.035) AND note green merges above!
        # Actually, check at phibar <= 0.030 where blue is completely distinct from green.
        raw = C.interior_mask(C.mask_blue(img), box, inset=5)
        # Keep only the upper solid curve (R/w > 50, excluding blue text "finite system" below R/w=75 at phibar>0.018)
        # Even cleaner: the solid blue curve is the single largest blue component (> 1200 px)!
        lab, n = ndimage.label(raw)
        best = max(range(1, n + 1), key=lambda i: int((lab == i).sum()))
        mask = lab == best
        return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=100)
    elif series == "fin_Rcrit":
        # Dotted blue lower branch: phibar in [0.014, 0.068], R/w < 55 (below text "finite system" at R/w ~ 35-80?)
        # Wait: "finite system" text sits at phibar ~ 0.019..0.038, R/w ~ 32..80; dotted blue curve is at R/w < 20 for phibar >= 0.018!
        raw = C.interior_mask(C.mask_blue(img), box, inset=5)
        # Exclude the text region: phibar > 0.0175 and R/w > 28 (keep everything below the text and left of it)
        dot_mask = raw.copy()
        r_txt = int(math.floor(y_axis.to_pixel(28.0)))
        c_txt = int(math.floor(x_axis.to_pixel(0.0175)))
        dot_mask[:r_txt, c_txt:] = False
        # Also exclude the solid upper branch (R/w > 52)
        r_up = int(math.floor(y_axis.to_pixel(52.0)))
        dot_mask[:r_up, :] = False
        return C.read_dotted_curve_at(
            dot_mask, x_axis, y_axis, at, min_area=12, max_area=400, max_aspect=3.0, extra_px=1.5
        )
    elif series == "inf_Rcrit":
        # Dotted orange curve: phibar in [0.0106, 0.0155], R/w in [30, 350]
        # Exclude orange text "infinite system" sitting to the left (phibar < 0.0085)!
        raw = C.interior_mask(C.mask_orange(img), box, inset=5)
        c_cut = int(math.floor(x_axis.to_pixel(0.0095)))
        raw[:, :c_cut] = False
        return C.read_dotted_curve_at(
            raw, x_axis, y_axis, at, min_area=12, max_area=400, max_aspect=3.0, extra_px=2.0
        )
    raise ValueError(f"Unknown Fig. 7b series {series!r}")
