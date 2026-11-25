"""Calibrated reader for SI Fig. S3B (p. 15): log-scale reflected intensity I = A^2 vs imaging depth z (42% dev).

Note: On the printed panel, imaging depth z runs vertically downward (0..105+ um) and log10 intensity
runs horizontally (10^9..10^12 arb. units). Because I(z) is a single-valued function of depth z (whereas
z(I) is multi-valued across the 4 reflection peaks), this reading indexes by x = axial_scan_depth (um)
and reads y = reflected_intensity (dimensionless, log10 scale), reporting each point's exact [col, row]
pixel location on the blue curve for visual verification.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_box, find_ticks_mask, interval_from_px, neutral_dark_mask


def read(image, at):
    c0, r0, c1, r1 = find_box(image, dark_thr=120.0, neutral_only=True)
    dark = neutral_dark_mask(image, dark_thr=120.0)

    # Bottom inward major x-ticks (longer than minor log ticks): 10^9, 10^10, 10^11 + right border 10^12
    xt_maj = find_ticks_mask(dark, r1, c0 + 10, c1 - 10, "x", -1, offset_lo=20, offset_hi=34)
    if len(xt_maj) != 3:
        raise RuntimeError(f"Fig. S3B expected 3 interior major log10 ticks (10^9..10^11), found {len(xt_maj)}: {xt_maj}")
    i_axis = Axis("log10", [(xt_maj[0], 1e9), (xt_maj[1], 1e10), (xt_maj[2], 1e11), (c1, 1e12)])

    # Left inward y-ticks: 0, 20, 40, 60, 80, 100 um (depth increases downward!)
    yt = find_ticks_mask(dark, c0, r0 + 15, r1 - 15, "y", +1, offset_lo=10, offset_hi=22)
    if len(yt) != 6:
        raise RuntimeError(f"Fig. S3B expected 6 depth y-ticks (0..100 um), found {len(yt)}: {yt}")
    z_axis = Axis("linear", list(zip(yt, [0.0, 20.0, 40.0, 60.0, 80.0, 100.0])))

    # Blue curve mask
    r_ch = image[:, :, 0].astype(float)
    g_ch = image[:, :, 1].astype(float)
    b_ch = image[:, :, 2].astype(float)
    blue = (b_ch - r_ch > 50.0) & (g_ch > r_ch + 20.0) & (b_ch > 150.0)

    zs = [float(v) for v in at] if at else [
        1.0, 3.0, 5.0, 7.0, 8.64, 10.0, 12.0, 15.0, 18.66, 21.0, 24.0, 28.29,
        33.0, 36.0, 39.0, 41.29, 44.0, 48.0, 52.0, 56.0, 60.0, 64.0, 67.54,
        72.0, 76.0, 80.0, 84.0, 88.0, 91.98, 96.0, 100.0, 104.0, 108.0,
    ]
    c_lo, c_hi = int(math.ceil(c0)) + 8, int(math.floor(c1)) - 8
    cal_i = max(i_axis.residual_px, 0.5)
    cal_z = max(z_axis.residual_px, 0.5)

    # First pass: locate curve columns at each depth row
    raw = []
    for zv in zs:
        rf = float(z_axis.to_pixel(zv))
        row = int(round(rf))
        if row <= int(r0) + 5 or row >= int(r1) - 5:
            raw.append((zv, rf, None, None))
            continue
        cols = np.nonzero(blue[row, c_lo:c_hi])[0] + c_lo
        if cols.size == 0:
            raw.append((zv, rf, None, None))
            continue
        cmid = 0.5 * (int(cols.min()) + int(cols.max())) + 0.5
        half_stroke = 0.5 * (int(cols.max()) - int(cols.min()) + 1)
        raw.append((zv, rf, cmid, half_stroke))

    out = []
    for i, (zv, rf, cmid, half_stroke) in enumerate(raw):
        if cmid is None:
            out.append({"x": float(zv), "y": None, "note": "blue curve not found at this depth"})
            continue
        prev_pt = next(((raw[j][1], raw[j][2]) for j in range(i - 1, -1, -1) if raw[j][2] is not None), None)
        next_pt = next(((raw[j][1], raw[j][2]) for j in range(i + 1, len(raw)) if raw[j][2] is not None), None)
        if prev_pt is not None and next_pt is not None and abs(next_pt[0] - prev_pt[0]) > 1e-6:
            slope = abs((next_pt[1] - prev_pt[1]) / (next_pt[0] - prev_pt[0]))
        else:
            slope = 0.0
        half_px = max(half_stroke, 2.0) + 0.5 + cal_i + min(slope * (0.5 + cal_z), 15.0)
        iv, lo, hi = interval_from_px(i_axis, cmid, half_px)
        out.append({
            "x": float(zv),
            "y": iv,
            "lo": lo,
            "hi": hi,
            "px": [cmid, rf],
            # Omit px_lo/px_hi because intensity uncertainty is horizontal, not vertical
            "note": f"horizontal stroke {2 * half_stroke:.1f} px; cal res z={z_axis.residual_px:.2f} logI={i_axis.residual_px:.2f} px",
        })
    return out
