"""First non-zero local maximum of 1D correlation function G_z(y) giving ridge spacing (SI Fig. 1E-F)."""

import numpy as np
from scipy.signal import find_peaks


def run(inputs, settings):
    curve = inputs["corr_curve"]
    r = np.asarray(curve["x"], dtype=float)
    g = np.asarray(curve["y"], dtype=float)
    order = np.argsort(r)
    r = r[order]
    g = g[order]

    min_lag = float(settings.get("min_lag_nm", 200.0))
    min_prom = float(settings.get("min_prominence", 0.0002))

    mask = r >= min_lag
    r_sub, g_sub = r[mask], g[mask]
    peaks, _ = find_peaks(g_sub, prominence=min_prom)
    if len(peaks) == 0:
        idx = int(np.argmax(g_sub))
    else:
        idx = int(peaks[0])

    # Quadratic sub-sample peak refinement when interior point
    if 0 < idx < len(r_sub) - 1:
        x0, x1, x2 = r_sub[idx - 1], r_sub[idx], r_sub[idx + 1]
        y0, y1, y2 = g_sub[idx - 1], g_sub[idx], g_sub[idx + 1]
        denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
        if abs(denom) > 1e-18:
            a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
            b = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / denom
            if a < 0:
                r_peak = -b / (2.0 * a)
                if x0 <= r_peak <= x2:
                    return {"s_ridge": float(r_peak)}

    return {"s_ridge": float(r_sub[idx])}
