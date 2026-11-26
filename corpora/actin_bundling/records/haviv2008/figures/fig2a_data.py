"""Calibrated reader for Fig. 2a (p. 4 / 450): measured distribution P(N_b) at phi = 3/1 (29 blue circles)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, calibrate_fig2a, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2a(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    blue = (B - R > 30.0) & (B - G > 30.0)
    # Restrict to plot region (allowing markers centered on bottom and right spines)
    mask = np.zeros_like(blue)
    mask[int(f["tr1"]) + 4 : int(f["br1"]) + 45, int(f["lc1"]) + 4 : int(f["rc1"]) + 45] = True
    blue &= mask

    lbl, _ = ndimage.label(blue)
    centers: list[tuple[float, float, int]] = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        sub = lbl[sl] == i
        area = int(sub.sum())
        if area < 400:
            continue
        w = sl[1].stop - sl[1].start
        if w > 105:
            # Two adjacent equal-bin zero-count markers touching along the x-axis: split horizontally
            mid = w // 2
            for part_slice, c_off in ((slice(0, mid), 0), (slice(mid, w), mid)):
                part = sub[:, part_slice]
                if int(part.sum()) < 400:
                    continue
                cy, cx = ndimage.center_of_mass(part)
                centers.append((float(cx + sl[1].start + c_off + 0.5), float(cy + sl[0].start + 0.5), int(part.sum())))
        else:
            cy, cx = ndimage.center_of_mass(sub)
            centers.append((float(cx + sl[1].start + 0.5), float(cy + sl[0].start + 0.5), area))

    err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0
    readings: list[dict] = []
    for cx, cy, area in sorted(centers, key=lambda t: t[0]):
        x_val = float(x_axis.to_value(cx))
        y_val, y_lo, y_hi = value_interval(y_axis, cy, err_px)
        # Markers sitting on the bottom x-axis spine represent zero-count bins (P >= 0)
        on_axis = abs(cy - float(f["B"])) <= 8.0
        if on_axis or y_val < 0.0:
            y_val = max(0.0, y_val)
            y_lo = 0.0
            y_hi = max(y_hi, y_val + float(abs(y_axis.b) * err_px))
        else:
            y_lo = max(0.0, y_lo)
        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "px": [cx, cy],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(y_lo)),
            "note": f"blue circle area={area} px{' (zero bin on x-axis)' if on_axis else ''}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
