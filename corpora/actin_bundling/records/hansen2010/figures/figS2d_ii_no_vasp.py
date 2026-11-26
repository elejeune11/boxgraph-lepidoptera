"""Calibrated reader for Fig. S2D(ii) (p. 16): K_obs vs [CP] at 100 mM KCl without VASP (filled circles)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_ii(image)
    cL, rT, cR, rB = fr["box"]
    blk = C.black_ink_mask(image, gray_max=110.0, sat_max=35.0)
    # Exclude top-right legend AND bottom-middle K_+ legend symbols (x >= 3.2, y <= 0.012)
    blk[: int(ya.to_pixel(0.038)), int(xa.to_pixel(4.5)) :] = False
    blk[int(ya.to_pixel(0.012)) :, int(xa.to_pixel(3.2)) :] = False
    sub = np.zeros_like(blk)
    sub[int(rT) + 4 : int(rB) - 4, int(cL) + 4 : int(cR) - 4] = blk[int(rT) + 4 : int(rB) - 4, int(cL) + 4 : int(cR) - 4]
    lab, n = ndimage.label(ndimage.binary_opening(sub, structure=np.ones((8, 8), bool)))
    cal_err = max(ya.residual_px, 0.5)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 35:
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
