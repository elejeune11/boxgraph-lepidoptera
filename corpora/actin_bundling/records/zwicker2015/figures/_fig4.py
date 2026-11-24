"""Calibrated readers for Fig. 4(a)-(b) (p. 6, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def _calib_fig4(img: np.ndarray):
    box = C.find_plot_box(img)
    # x-axis: 0 (left spine), 5, 10 (top spine major ticks)
    xt_all = [p for p, _ in C.find_inward_ticks(img, box, "top", min_len=35, max_len=48)]
    xt_reg = C.match_ticks_regular([box["c_left"]] + xt_all, 3)
    x_marks = [(xt_reg[0], 0.0), (xt_reg[1], 5.0), (xt_reg[2], 10.0)]
    x_axis = C.fit_axis("linear", x_marks)

    # y-axis: 0 (bottom spine), 5, 10 (right spine major ticks)
    yr_all = [p for p, _ in C.find_inward_ticks(img, box, "right", min_len=35, max_len=48)]
    yr_reg = C.match_ticks_regular(yr_all + [box["r_bot"]], 3)
    y_marks = [(yr_reg[0], 10.0), (yr_reg[1], 5.0), (yr_reg[2], 0.0)]
    y_axis = C.fit_axis("linear", y_marks)
    return box, x_axis, y_axis


def read_fig4_nullcline_R1(img: np.ndarray, at: list[float]) -> list[dict]:
    """Read the medium-blue nullcline (dR_1/dt = 0) as R_2/w vs R_1/w."""
    box, x_axis, y_axis = _calib_fig4(img)
    blue = C.interior_mask(C.mask_blue(img), box, inset=8)
    return C.read_curve_at(blue, x_axis, y_axis, at, search=4, min_run=3, max_run=120, extra_px=0.5)


def read_fig4_sym_fp(img: np.ndarray, filled: bool) -> list[dict]:
    """Read the symmetric large two-droplet fixed point (R_sym/w, R_sym/w) on the diagonal."""
    box, x_axis, y_axis = _calib_fig4(img)
    # Search in the diagonal neighborhood R1, R2 in [5.5, 9.0]
    c0, c1 = int(math.floor(x_axis.to_pixel(5.5))), int(math.ceil(x_axis.to_pixel(9.0)))
    r0, r1 = int(math.floor(y_axis.to_pixel(9.0))), int(math.ceil(y_axis.to_pixel(5.5)))

    if filled:
        # Solid black disk (Fig. 4b): erode thin black trajectory lines away to isolate the large disk
        black = C.mask_black(img, thresh=60)
        sub = np.zeros_like(black, dtype=bool)
        sub[r0:r1, c0:c1] = black[r0:r1, c0:c1]
        # Trajectory lines are ~6-8 px wide; disk diameter is ~65 px -> opening by 21x21 isolates disk cleanly!
        disk = ndimage.binary_opening(sub, structure=np.ones((21, 21), dtype=bool))
        pts = C.read_filled_markers(disk, x_axis, y_axis, min_area=500, max_area=12000, extra_px=1.0)
        if len(pts) != 1:
            raise RuntimeError(f"Expected 1 filled symmetric fixed-point disk in Fig. 4b, got {len(pts)}")
        return pts
    else:
        # Open black circle with white interior (Fig. 4a): find white region enclosed by black ring
        black = C.mask_black(img, thresh=80)
        # Fill holes of black mask inside [r0:r1, c0:c1] -> the circle's hollow interior is the unique circular hole!
        sub_b = black[r0:r1, c0:c1]
        # Dilate black slightly by 2 px so ring is guaranteed closed
        closed = ndimage.binary_dilation(sub_b, iterations=2)
        filled_b = ndimage.binary_fill_holes(closed)
        holes = filled_b & ~closed
        lab, n = ndimage.label(holes)
        cands = []
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            area = len(rr)
            if 400 <= area <= 8000:
                h = rr.max() - rr.min() + 1
                w = cc.max() - cc.min() + 1
                if max(h, w) / max(1, min(h, w)) < 1.35:
                    cands.append((c0 + float(cc.mean() + 0.5), r0 + float(rr.mean() + 0.5), area))
        if len(cands) != 1:
            raise RuntimeError(f"Expected 1 open symmetric fixed-point circle in Fig. 4a, got {cands}")
        cx, cy, area = cands[0]
        err_px = 2.0 + max(y_axis.residual_px, 0.5)
        yv = y_axis.to_value(cy)
        a, b = y_axis.to_value(cy - err_px), y_axis.to_value(cy + err_px)
        return [
            {
                "x": float(x_axis.to_value(cx)),
                "y": float(yv),
                "lo": float(min(a, b)),
                "hi": float(max(a, b)),
                "px": [cx, cy],
                "px_lo": float(cy - err_px),
                "px_hi": float(cy + err_px),
                "note": f"open circle interior area={area} px",
            }
        ]
