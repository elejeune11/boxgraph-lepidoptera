"""Calibrated reader for Fig. 6D (p. 9): K_obs vs [CP] without VASP (filled black circles)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6d(image)
    cL, rT, cR, rB = fr["box"]
    blk = C.black_ink_mask(image, gray_max=105.0, sat_max=35.0)
    # Exclude both legend regions (top-right and bottom-middle K_+ annotations)
    blk[: int(ya.to_pixel(0.036)), int(xa.to_pixel(7.5)) :] = False
    blk[int(ya.to_pixel(0.013)) :, int(xa.to_pixel(7.5)) :] = False
    sub = np.zeros_like(blk)
    sub[int(rT) + 4 : int(rB) - 4, int(cL) + 2 : int(cR) - 4] = blk[int(rT) + 4 : int(rB) - 4, int(cL) + 2 : int(cR) - 4]
    lab, n = ndimage.label(ndimage.binary_opening(sub, structure=np.ones((9, 9), bool)))
    cal_err = max(ya.residual_px, 0.5)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 220:
            cm, rm = float(cc.mean()) + 0.5, float(rr.mean()) + 0.5
            xv = float(xa.to_value(cm))
            y, lo, hi = C.val_interval(ya, rm, 2.0 + cal_err)
            pts.append({
                "x": round(xv, 3), "y": y, "lo": lo, "hi": hi,
                "px": [cm, rm], "px_lo": rm - 2.0 - cal_err, "px_hi": rm + 2.0 + cal_err,
                "note": f"filled circle {len(rr)}px",
            })
    pts.sort(key=lambda p: p["x"])
    return pts
