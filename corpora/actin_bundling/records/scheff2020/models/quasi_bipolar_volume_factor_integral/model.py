import math
import numpy as np


def run(inputs, settings):
    x = float(inputs["x"])
    n_xi = int(settings.get("n_xi", 32))
    n_eta = int(settings.get("n_eta", 24))
    eta0 = 2.0 * math.atan(1.0 / x)
    u, wu = np.polynomial.legendre.leggauss(n_xi)
    v_n, wv = np.polynomial.legendre.leggauss(n_eta)
    xi = 0.25 * math.pi * (u + 1.0)
    wxi = 0.25 * math.pi * wu
    eta = 0.5 * eta0 * (v_n + 1.0)
    weta = 0.5 * eta0 * wv
    XI = xi[:, None]
    ETA = eta[None, :]
    W = wxi[:, None] * weta[None, :]
    sxi = np.sin(XI)
    seta = np.sin(ETA)
    Z = 1.0 + sxi * np.cos(ETA)
    integrand = (sxi * sxi * seta) / (Z * Z * Z)
    vol = float(4.0 * math.pi * np.sum(W * integrand))
    return {"v": vol}
