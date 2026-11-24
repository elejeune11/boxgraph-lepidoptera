"""Calibrated readers for Fig. 10(a)-(b) (p. 10, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


# --------------------------------------------------------------------------
# Fig. 10(a): unequal-core two-droplet phase portrait R_2/w vs R_1/w
# --------------------------------------------------------------------------


def _calib_10a(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = [p for p, _ in C.find_inward_ticks(img, box, "top", 32, 48)]
    xt_reg = C.match_ticks_regular([box["c_left"]] + xt_all, 3)
    x_axis = C.fit_axis("linear", [(xt_reg[0], 0.0), (xt_reg[1], 50.0), (xt_reg[2], 100.0)])

    yr_all = [p for p, _ in C.find_inward_ticks(img, box, "right", 32, 48)]
    yr_reg = C.match_ticks_regular(yr_all + [box["r_bot"]], 3)
    y_axis = C.fit_axis("linear", [(yr_reg[0], 100.0), (yr_reg[1], 50.0), (yr_reg[2], 0.0)])
    return box, x_axis, y_axis


def read_fig10a_nullcline(img: np.ndarray, at: list[float], which: str) -> list[dict]:
    box, x_axis, y_axis = _calib_10a(img)
    if which == "R1":
        mask = C.interior_mask(C.mask_blue(img), box, inset=6)
    elif which == "R2":
        mask = C.interior_mask(C.mask_orange(img), box, inset=6)
    else:
        raise ValueError(f"Unknown Fig. 10a nullcline {which!r}")
    return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=3, max_run=140, extra_px=0.5)


def read_fig10a_fp(img: np.ndarray) -> list[dict]:
    """Read the unique stable fixed-point black disk (R_1/w, R_2/w) in Fig. 10(a)."""
    box, x_axis, y_axis = _calib_10a(img)
    c0, c1 = int(math.floor(x_axis.to_pixel(95.0))), int(math.ceil(x_axis.to_pixel(135.0)))
    r0, r1 = int(math.floor(y_axis.to_pixel(75.0))), int(math.ceil(y_axis.to_pixel(35.0)))
    black = C.mask_black(img, thresh=60)
    sub = np.zeros_like(black, dtype=bool)
    sub[r0:r1, c0:c1] = black[r0:r1, c0:c1]
    disk = ndimage.binary_opening(sub, structure=np.ones((27, 27), dtype=bool))
    pts = C.read_filled_markers(disk, x_axis, y_axis, min_area=1000, max_area=12000, shape="circle", extra_px=1.0)
    if len(pts) != 1:
        raise RuntimeError(f"Expected 1 filled fixed-point disk in Fig. 10a, got {len(pts)}")
    return pts


# --------------------------------------------------------------------------
# Fig. 10(b): steady-state droplet volumes V_1, V_2, V_1 + V_2 vs Q_1/(w^3 k_0)
# --------------------------------------------------------------------------


def _calib_10b(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = [p for p, _ in C.find_inward_ticks(img, box, "top", 30, 48)]
    xt_reg = C.match_ticks_regular([box["c_left"]] + xt_all, 5)
    x_axis = C.fit_axis(
        "linear",
        [(xt_reg[0], 0.0), (xt_reg[1], 200.0), (xt_reg[2], 400.0), (xt_reg[3], 600.0), (xt_reg[4], 800.0)],
    )

    yl_all = sorted(C.filter_ticks_by_len(C.find_inward_ticks(img, box, "left", 30, 48), 30, 48))
    if len(yl_all) != 3:
        raise RuntimeError(f"Fig. 10b expected 3 left major y-ticks (6e6, 4e6, 2e6), got {yl_all}")
    y_axis = C.fit_axis("linear", [(yl_all[0], 6e6), (yl_all[1], 4e6), (yl_all[2], 2e6), (box["r_bot"], 0.0)])
    return box, x_axis, y_axis


def read_fig10b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, x_axis, y_axis = _calib_10b(img)
    # Exclude legend rectangle (Q_1 in [500, 820], V in [2.0e6, 5.0e6])
    c_l0, c_l1 = int(math.floor(x_axis.to_pixel(500.0))), int(math.ceil(x_axis.to_pixel(820.0)))
    r_l0, r_l1 = int(math.floor(y_axis.to_pixel(5.0e6))), int(math.ceil(y_axis.to_pixel(2.0e6)))

    if series == "V1":
        # For Q_1 <= 209, the blue curve V_1 = 0 is clipped along the bottom spine (only its upper half
        # protrudes above the black spine at bot_lo); extend any run touching bot_lo symmetrically to r_bot!
        bm = C.mask_blue(img)
        mask = np.zeros_like(bm, dtype=bool)
        mask[box["top_in"] + 4 : box["bot_hi"] + 2, box["left_in"] + 4 : box["right_in"] - 4] = bm[
            box["top_in"] + 4 : box["bot_hi"] + 2, box["left_in"] + 4 : box["right_in"] - 4
        ]
        mask[r_l0:r_l1, c_l0:c_l1] = False
        rs = C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=2, max_run=90)
        for r in rs:
            if r.get("y") is not None and r["px"][1] >= box["bot_lo"] - 10:
                # Clipped against bottom axis spine y = 0
                r_top_edge = r["px_lo"]
                r["px"][1] = float(box["r_bot"])
                r["y"] = 0.0
                r["lo"] = 0.0
                r["hi"] = float(y_axis.to_value(r_top_edge))
                r["px_lo"] = float(r_top_edge)
                r["px_hi"] = float(box["r_bot"])
                r["note"] = "blue curve clipped on bottom axis spine (V_1 = 0)"
        return rs
    elif series == "V2":
        om = C.interior_mask(C.mask_orange(img), box, inset=4)
        om[r_l0:r_l1, c_l0:c_l1] = False
        return C.read_curve_at(om, x_axis, y_axis, at, search=4, min_run=2, max_run=90)
    elif series == "Vtot":
        # Dotted black curve V_1 + V_2 lies strictly above V = 6.5e6 (well above legend and orange/blue curves)
        blk = C.interior_mask(C.mask_black(img, thresh=75), box, inset=6)
        r_lo = int(math.ceil(y_axis.to_pixel(6.5e6)))
        r_hi = int(math.floor(y_axis.to_pixel(7.6e6)))
        c_lo = int(math.floor(x_axis.to_pixel(210.0)))
        dot_mask = np.zeros_like(blk, dtype=bool)
        dot_mask[r_hi:r_lo, c_lo : box["right_in"] - 6] = blk[r_hi:r_lo, c_lo : box["right_in"] - 6]
        return C.read_dotted_curve_at(
            dot_mask, x_axis, y_axis, at, min_area=12, max_area=300, max_aspect=3.0, extra_px=1.5
        )
    raise ValueError(f"Unknown Fig. 10b series {series!r}")
