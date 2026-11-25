"""1D profile best-fit slope flattening rotation and peak-height mean/SD (SI p. 11, Fig. S10B-C).

Rotates (x_nm, h_nm) in physical nm coordinates by -arctan(m) where m is the OLS slope,
locates peaks via scipy.signal.find_peaks (analogous to MATLAB findpeaks), identifies the
minimum trough on each side between adjacent peaks, and computes mean and sample SD peak height.
"""

import math
import numpy as np
from scipy.signal import find_peaks


def run(inputs, settings):
    curve = inputs["h_raw"]
    x_um = np.asarray(curve["x"], dtype=float)
    h_nm = np.asarray(curve["y"], dtype=float)
    if x_um.size < 5:
        return {
            "mean_peak_height": float("nan"),
            "sd_peak_height": float("nan"),
            "slope_angle": float("nan"),
        }

    order = np.argsort(x_um)
    x_um, h_nm = x_um[order], h_nm[order]
    x_nm = x_um * 1000.0

    # OLS linear fit in physical nm coordinates: h = m * x + b
    m, _b = np.polyfit(x_nm, h_nm, 1)
    theta = math.atan(float(m))
    c, s = math.cos(-theta), math.sin(-theta)
    x_rot = c * x_nm - s * h_nm
    h_rot = s * x_nm + c * h_nm
    h_rot = h_rot - float(np.min(h_rot))

    # Resample onto uniform grid for distance-based peak detection
    n_grid = max(256, x_rot.size)
    xg = np.linspace(x_rot[0], x_rot[-1], n_grid)
    hg = np.interp(xg, x_rot, h_rot)
    dx_um = float((xg[1] - xg[0]) / 1000.0)

    min_dist_pts = max(1, int(round(float(settings.get("min_peak_distance_um", 0.8)) / dx_um)))
    min_prom = float(settings.get("min_peak_prominence_nm", 15.0))
    peaks, _ = find_peaks(hg, prominence=min_prom, distance=min_dist_pts)

    if peaks.size == 0:
        return {
            "mean_peak_height": float("nan"),
            "sd_peak_height": float("nan"),
            "slope_angle": theta,
        }

    heights = []
    for i, pk in enumerate(peaks):
        left_bound = 0 if i == 0 else int(peaks[i - 1])
        right_bound = (hg.size - 1) if i == (peaks.size - 1) else int(peaks[i + 1])
        tl = left_bound + int(np.argmin(hg[left_bound : pk + 1]))
        tr = pk + int(np.argmin(hg[pk : right_bound + 1]))
        if xg[tr] == xg[tl]:
            continue
        h_base = hg[tl] + ((xg[pk] - xg[tl]) / (xg[tr] - xg[tl])) * (hg[tr] - hg[tl])
        heights.append(float(hg[pk] - h_base))

    if not heights:
        return {
            "mean_peak_height": float("nan"),
            "sd_peak_height": float("nan"),
            "slope_angle": theta,
        }

    arr = np.asarray(heights, dtype=float)
    mean_h = float(np.mean(arr))
    sd_h = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
    return {
        "mean_peak_height": mean_h,
        "sd_peak_height": sd_h,
        "slope_angle": theta,
    }
