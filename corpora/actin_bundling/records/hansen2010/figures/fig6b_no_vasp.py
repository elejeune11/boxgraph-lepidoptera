"""Calibrated reader for Fig. 6B (p. 9): subunits added before capping vs [CP] without VASP (filled circles)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6b(image)
    g = C.gray(image)
    cL, rT, cR, rB = [int(round(v)) for v in fr["box"]]
    blk = g < 90.0
    # Exclude upper-right legend and top-left steep open-circle overlap region
    blk[: int(ya.to_pixel(4800.0)), int(xa.to_pixel(6.0)) :] = False
    blk[: int(ya.to_pixel(3200.0)), :] = False
    sub = np.zeros_like(blk)
    sub[rT + 4 : rB - 4, cL + 4 : cR - 4] = blk[rT + 4 : rB - 4, cL + 4 : cR - 4]
    lab, n = ndimage.label(ndimage.binary_opening(sub, structure=np.ones((10, 10), bool)))
    cal_err = max(ya.residual_px, 0.5)
    raw = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 200:
            cm, rm = float(cc.mean()) + 0.5, float(rr.mean()) + 0.5
            # Merged [0, 0.125] pair at top-left has area > 950; read 0 nM from its top half & 0.125 nM from bottom half
            if len(rr) > 950 and ya.to_value(rm) > 2500.0:
                r_mid_split = float(np.median(rr))
                top_m = rr <= r_mid_split
                bot_m = rr > r_mid_split
                raw.append((float(cc[top_m].mean()) + 0.5, float(rr[top_m].mean()) + 0.5, "merged 0 nM top"))
                raw.append((float(cc[bot_m].mean()) + 0.5, float(rr[bot_m].mean()) + 0.5, "merged 0.125 nM bot"))
            else:
                raw.append((cm, rm, f"filled circle {len(rr)}px"))
    raw.sort(key=lambda t: t[0])
    # Match the 10 nominal [CP] concentrations in Fig. 6B
    nominal = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0]
    if len(raw) != len(nominal):
        raise ValueError(f"Fig. 6B -VASP: expected {len(nominal)} circles, found {len(raw)}")
    out = []
    for cp, (cm, rm, note) in zip(nominal, raw):
        y, lo, hi = C.val_interval(ya, rm, 2.5 + cal_err)
        out.append({
            "x": cp, "y": y, "lo": lo, "hi": hi,
            "px": [cm, rm], "px_lo": rm - 2.5 - cal_err, "px_hi": rm + 2.5 + cal_err,
            "note": note,
        })
    return out
