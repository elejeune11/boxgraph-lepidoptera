"""Calibrated readers for Fig. 8(a)-(b) (p. 8, Zwicker et al. 2015)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def _calib_fig8(img: np.ndarray):
    box = C.find_plot_box(img)
    xt_all = [p for p, _ in C.find_inward_ticks(img, box, "top", 32, 48)]
    xt_reg = C.match_ticks_regular([box["c_left"]] + xt_all, 3)
    x_axis = C.fit_axis("linear", [(xt_reg[0], 0.0), (xt_reg[1], 50.0), (xt_reg[2], 100.0)])

    yr_all = [p for p, _ in C.find_inward_ticks(img, box, "right", 32, 48)]
    yr_reg = C.match_ticks_regular(yr_all + [box["r_bot"]], 3)
    y_axis = C.fit_axis("linear", [(yr_reg[0], 100.0), (yr_reg[1], 50.0), (yr_reg[2], 0.0)])
    return box, x_axis, y_axis


def read_fig8_nullcline(img: np.ndarray, at: list[float], which: str) -> list[dict]:
    """Read 'R1' (blue nullcline dR_1/dt = 0) or 'R2' (orange nullcline dR_2/dt = 0)."""
    box, x_axis, y_axis = _calib_fig8(img)
    if which == "R1":
        mask = C.interior_mask(C.mask_blue(img), box, inset=8)
    elif which == "R2":
        mask = C.interior_mask(C.mask_orange(img), box, inset=8)
    else:
        raise ValueError(f"Unknown Fig. 8 nullcline {which!r}")
    return C.read_curve_at(mask, x_axis, y_axis, at, search=4, min_run=3, max_run=140, extra_px=0.5)


def read_fig8_sym_fp(img: np.ndarray, filled: bool) -> list[dict]:
    """Read the symmetric large two-droplet fixed point (R_sym/w, R_sym/w) on the diagonal."""
    box, x_axis, y_axis = _calib_fig8(img)
    c0, c1 = int(math.floor(x_axis.to_pixel(70.0))), int(math.ceil(x_axis.to_pixel(115.0)))
    r0, r1 = int(math.floor(y_axis.to_pixel(115.0))), int(math.ceil(y_axis.to_pixel(70.0)))

    if filled:
        black = C.mask_black(img, thresh=60)
        sub = np.zeros_like(black, dtype=bool)
        sub[r0:r1, c0:c1] = black[r0:r1, c0:c1]
        disk = ndimage.binary_opening(sub, structure=np.ones((21, 21), dtype=bool))
        pts = C.read_filled_markers(disk, x_axis, y_axis, min_area=500, max_area=12000, extra_px=1.0)
        if len(pts) != 1:
            raise RuntimeError(f"Expected 1 filled fixed-point disk in Fig. 8b, got {len(pts)}")
        return pts
    else:
        black = C.mask_black(img, thresh=80)
        sub_b = black[r0:r1, c0:c1]
        closed = ndimage.binary_dilation(sub_b, iterations=2)
        holes = ndimage.binary_fill_holes(closed) & ~closed
        lab, n = ndimage.label(holes)
        cands = []
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            area = len(rr)
            if 400 <= area <= 8000:
                h, w = rr.max() - rr.min() + 1, cc.max() - cc.min() + 1
                if max(h, w) / max(1, min(h, w)) < 1.35:
                    cands.append((c0 + float(cc.mean() + 0.5), r0 + float(rr.mean() + 0.5), area))
        if len(cands) != 1:
            raise RuntimeError(f"Expected 1 open symmetric fixed-point circle in Fig. 8a, got {cands}")
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
