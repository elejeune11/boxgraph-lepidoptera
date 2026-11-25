import math
import numpy as np


def run(inputs, settings):
    x = float(inputs["x"])
    n_xi = int(settings.get("n_xi", 32))
    eta0 = 2.0 * math.atan(1.0 / x)
    s0 = math.sin(eta0)
    c0 = math.cos(eta0)
    u, wu = np.polynomial.legendre.leggauss(n_xi)
    # By symmetry xi <-> pi - xi, integrate on [0, pi/2] and multiply by 2
    xi = 0.25 * math.pi * (u + 1.0)
    wxi = 0.25 * math.pi * wu
    sxi = np.sin(xi)
    Z0 = 1.0 + sxi * c0
    integrand = (sxi * s0) / (Z0 * Z0)
    # 2*pi from phi, 2 from symmetry on xi
    f_s = float(4.0 * math.pi * np.sum(wxi * integrand))
    return {"f_s": f_s}
