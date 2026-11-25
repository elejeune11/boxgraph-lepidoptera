"""AFM ridge height extraction from transverse surface height profile z(x) (SI p. 7, Fig. S12)."""

import numpy as np
from scipy.signal import find_peaks


def run(inputs, settings):
    prof = inputs["profile"]
    x = np.asarray(prof["x"], dtype=float)
    z = np.asarray(prof["y"], dtype=float)
    order = np.argsort(x)
    x = x[order]
    z = z[order]

    z_min = float(np.min(z))
    h_max = float(np.max(z) - z_min)

    dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
    min_dist_um = float(settings.get("min_distance_um", 0.4))
    min_prom_um = float(settings.get("min_prominence_um", 0.2))
    dist_pts = max(1, int(round(min_dist_um / max(dx, 1e-9))))

    peaks, _ = find_peaks(z, prominence=min_prom_um, distance=dist_pts)
    if len(peaks) == 0:
        h_mean = h_max
    else:
        h_mean = float(np.mean(z[peaks] - z_min))

    return {"H_mean": h_mean, "H_max": h_max}
