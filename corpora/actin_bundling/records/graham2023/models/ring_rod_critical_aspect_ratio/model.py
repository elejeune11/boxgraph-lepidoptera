import math
from scipy.optimize import brentq


def _ellipse_perimeter(a, b, n_quad):
    dtheta = (0.5 * math.pi) / n_quad
    total = 0.0
    for k in range(n_quad):
        theta = (k + 0.5) * dtheta
        s = math.sin(theta)
        c = math.cos(theta)
        total += math.sqrt((a * s) ** 2 + (b * c) ** 2)
    return 4.0 * total * dtheta


def run(inputs, settings):
    R = inputs["R"]
    L = inputs["L"]
    n_quad = max(10, int(round(settings.get("n_quad", 400))))
    root_tol = float(settings.get("root_tol", 1e-10))

    target = 2.0 * L

    def residual(ar):
        s = math.sqrt(ar)
        return _ellipse_perimeter(R * s, R / s, n_quad) - target

    if residual(1.0) >= 0.0:
        return {"ar_crit": 1.0}
    hi = 2.0
    while residual(hi) < 0.0:
        hi *= 2.0
    ar_star = brentq(residual, 1.0, hi, xtol=root_tol, rtol=root_tol)
    return {"ar_crit": float(ar_star)}
