"""Calibrated reader for Fig. S2D(i) (p. 16): subunits added before capping at 100 mM KCl without VASP (filled circles)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_i(image)
    g = C.gray(image)
    cL, rT, cR, rB = [int(round(v)) for v in fr["box"]]
    blk = g < 90.0
    blk[: int(ya.to_pixel(3900.0)), int(xa.to_pixel(5.0)) :] = False
    sub = np.zeros_like(blk)
    sub[rT + 4 : rB - 4, cL + 4 : cR - 4] = blk[rT + 4 : rB - 4, cL + 4 : cR - 4]
    lab, n = ndimage.label(ndimage.binary_opening(sub, structure=np.ones((9, 9), bool)))
    cal_err = max(ya.residual_px, 0.5)
    raw = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 280:
            raw.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5, len(rr)))
    raw.sort(key=lambda t: t[0])
    nominal = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0]
    if len(raw) != len(nominal):
        raise ValueError(f"Fig. S2D(i) -VASP: expected 6 filled circles, found {len(raw)}")
    out = []
    for cp, (cm, rm, area) in zip(nominal, raw):
        y, lo, hi = C.val_interval(ya, rm, 2.0 + cal_err)
        out.append({
            "x": cp, "y": y, "lo": lo, "hi": hi,
            "px": [cm, rm], "px_lo": rm - 2.0 - cal_err, "px_hi": rm + 2.0 + cal_err,
            "note": f"filled circle {area}px",
        })
    return out
