import math
import numpy as np
from scipy.interpolate import splprep, splev


def run(inputs, settings):
    R = float(inputs["R"])
    ar = float(inputs["aspect_ratio"])
    n_edge = int(round(settings.get("n_edge_pixels", 40)))
    oversample = int(round(settings.get("oversample_factor", 5)))
    spline_s = float(settings.get("spline_s", 0.0))

    # Area-preserving ellipse (area = pi * R^2; circle at aspect_ratio = 1)
    a = R * math.sqrt(ar)
    b = R / math.sqrt(ar)
    theta = np.linspace(0.0, 2.0 * np.pi, n_edge, endpoint=False)
    x_e = a * np.cos(theta)
    y_e = b * np.sin(theta)

    # Periodic smoothing spline fit to closed contour
    x_closed = np.append(x_e, x_e[0])
    y_closed = np.append(y_e, y_e[0])
    tck, _ = splprep([x_closed, y_closed], per=True, s=spline_s, k=3)

    m = max(16, oversample * n_edge)
    u = np.linspace(0.0, 1.0, m, endpoint=False)
    du = 1.0 / m
    x, y = splev(u, tck)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Central finite differences for tangent, unit normal n, and curvature div(n)
    dx = (np.roll(x, -1) - np.roll(x, 1)) / (2.0 * du)
    dy = (np.roll(y, -1) - np.roll(y, 1)) / (2.0 * du)
    speed = np.hypot(dx, dy)
    tx = dx / speed
    ty = dy / speed
    nx = ty
    ny = -tx

    dnx = (np.roll(nx, -1) - np.roll(nx, 1)) / (2.0 * du)
    dny = (np.roll(ny, -1) - np.roll(ny, 1)) / (2.0 * du)
    kappa = (dnx * tx + dny * ty) / speed
    abs_k = np.abs(kappa)

    kappa_mean = float(np.mean(abs_k))
    kappa_max = float(np.max(abs_k))
    kappa_norm_min = float(np.min(abs_k) / kappa_max) if kappa_max > 0.0 else 0.0
    return {
        "kappa_mean": kappa_mean,
        "kappa_max": kappa_max,
        "kappa_norm_min": kappa_norm_min,
    }
