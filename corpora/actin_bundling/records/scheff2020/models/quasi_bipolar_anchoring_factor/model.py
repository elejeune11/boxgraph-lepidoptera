import math
import numpy as np


def run(inputs, settings):
    x = float(inputs["x"])
    y = float(inputs["y"])
    n_xi = int(settings.get("n_xi", 48))
    eta0 = 2.0 * math.atan(1.0 / x)
    s0 = math.sin(eta0)
    c0 = math.cos(eta0)
    t = 1.0 / (y * y)
    u, wu = np.polynomial.legendre.leggauss(n_xi)
    s = 0.5 * (u + 1.0)
    ws = 0.5 * wu
    # Quadratic clustering near xi = 0 on [0, pi/2] (symmetric about pi/2)
    xi = 0.5 * math.pi * (s * s)
    wxi = math.pi * s * ws
    sxi = np.sin(xi)
    cxi = np.cos(xi)
    Z0 = 1.0 + sxi * c0
    t2N0 = (t * sxi * c0 + 0.5 * (1.0 - t) * Z0) ** 2 + t * (sxi * sxi) * (s0 * s0)
    h_surf = (sxi * s0) / (Z0 * Z0)
    integrand = h_surf * (((1.0 - t) ** 2) * (cxi * cxi) * (s0 * s0)) / (4.0 * t2N0)
    # 2*pi from phi, 2 from xi symmetry [0, pi/2] -> [0, pi]
    return {"f_an": float(4.0 * math.pi * np.sum(wxi * integrand))}
