"""Calibrated reader for Fig. 4I (p. 7): actin monomers delivered per VASP tetramer (filled circles, right axis)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, _, ya_r, fr = C.calibrate_fig4i(image)
    cL, rT, cR, rB = fr["box"]
    g = C.gray(image)
    dk = np.zeros_like(g, dtype=bool)
    dk[int(rT) + 8 : int(rB) - 8, int(cL) + 8 : int(cR) - 8] = g[int(rT) + 8 : int(rB) - 8, int(cL) + 8 : int(cR) - 8] < 140
    # Filled circles survive 9x9 binary opening; thin error bars and open-circle outlines do not
    blobs = ndimage.binary_opening(dk, structure=np.ones((9, 9), bool))
    lab, n = ndimage.label(blobs)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 100:
            pts.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5, int(len(rr))))
    pts.sort(key=lambda p: p[0])
    if len(pts) != 5:
        raise ValueError(f"Fig. 4I monomers: expected 5 filled markers, found {len(pts)}")
    cal_err = max(ya_r.residual_px, 0.5)
    out = []
    for c_m, r_m, area in pts:
        x_val = float(xa.to_value(c_m))
        y, lo, hi = C.val_interval(ya_r, r_m, 1.5 + cal_err)
        eb = C.measure_error_bar_vertical(
            image, c_m, r_m, int(rT) + 4, int(rB) - 4, ya_r,
            dark_thresh=135.0, stem_half_w=3, gap_tol=3, mirror_lower_if_hidden=False,
        )
        out.append({
            "x": round(x_val, 3), "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_m, r_m], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"filled circle {area}px; {eb['note']}",
        })
    return out
